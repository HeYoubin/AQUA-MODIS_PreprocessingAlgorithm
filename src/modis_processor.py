import os
import sys
import re
import warnings
import numpy as np
from datetime import date
from pyhdf.SD import SD, SDC
from osgeo import gdal, osr

# ==========================================
# 0. 环境配置与警告屏蔽 (Environment Setup)
# ==========================================
# 忽略一些计算中预期的警告（如除以零产生的 NaN）
warnings.filterwarnings("ignore", category=RuntimeWarning)

# 1. 解决 GDAL 升级带来的 FutureWarning
gdal.UseExceptions()

# 2. 解决 PROJ 库版本冲突（强制使用当前 conda 环境的 proj.db）
# 应对 ERROR 1: PROJ: proj_create_from_database 报错
conda_proj_path = os.path.join(sys.prefix, 'Library', 'share', 'proj')
os.environ['PROJ_LIB'] = conda_proj_path
os.environ['PROJ_DATA'] = conda_proj_path

# ==========================================
# 模块 1：文件与 IO 操作 (IO Utilities)
# 负责 HDF 数据的读取和 GeoTIFF 的生成
# ==========================================
def search_modis_files(data_dir, start_date=None, end_date=None):
    """
    搜索指定日期范围内的 MOD03 定位文件。
    如果不指定日期，默认返回当月至今的数据。
    """
    if start_date is None or end_date is None:
        end_date = date.today()
        start_date = end_date.replace(day=1)

    pattern = re.compile(r".*_(\d{4})_(\d{2})_(\d{2})_.*\.MOD03\.hdf$")
    matched = []

    for fname in os.listdir(data_dir):
        if not fname.endswith('.MOD03.hdf'):
            continue
        m = pattern.match(fname)
        if m:
            year, month, day = map(int, m.groups())
            try:
                file_date = date(year, month, day)
                if start_date <= file_date <= end_date:
                    matched.append((file_date, os.path.join(data_dir, fname)))
            except ValueError:
                continue

    matched.sort(key=lambda x: x[0])
    return [path for _, path in matched]

def parse_filename_info(filename):
    """
    解析文件名，提取平台、轨道时间等关键信息。
    """
    base = os.path.basename(filename)
    tokens = base.split('.')[0].split('_')
    if len(tokens) >= 7:
        platform = tokens[0] # 如 AQUA
        datetime_str = f"{tokens[2]}_{tokens[3]}_{tokens[4]}_{tokens[5]}_{tokens[6]}" # YYYY_MM_DD_HH_MM
        orbit_time = f"{tokens[2]}{tokens[3]}{tokens[4]}_{tokens[5]}{tokens[6]}"
        ad_type = tokens[7] if len(tokens) > 7 else "UNKNOWN"
        return platform, orbit_time, datetime_str, ad_type
    return 'UNKNOWN', 'UNKNOWN_TIME', 'UNKNOWN_Date', 'UNKNOWN'

def match_science_data(orbit_data_path, datetime_str):
    """
    根据时间戳匹配对应的 L1B 科学数据 (MOD021KM/MYD021KM)。
    """
    for fname in os.listdir(orbit_data_path):
        if datetime_str in fname and ("MOD021KM" in fname or "MYD021KM" in fname) and fname.endswith('.hdf'):
            return os.path.join(orbit_data_path, fname)
    return None


def read_hdf_dataset(sd_obj, dataset_name):
    """
    通用函数：读取 HDF 子数据集及其定标属性。
    返回: (数据数组, radiance_scales, radiance_offsets, reflectance_scales, reflectance_offsets)
    """
    sds = sd_obj.select(dataset_name)
    arr = sds.get().astype(np.float32)
    attrs = sds.attributes()

    # 尝试获取辐射和反射率定标系数，如果没有则返回 None
    rad_scales = np.array(attrs.get('radiance_scales')) if 'radiance_scales' in attrs else None
    rad_offsets = np.array(attrs.get('radiance_offsets')) if 'radiance_offsets' in attrs else None
    ref_scales = np.array(attrs.get('reflectance_scales')) if 'reflectance_scales' in attrs else None
    ref_offsets = np.array(attrs.get('reflectance_offsets')) if 'reflectance_offsets' in attrs else None

    return arr, rad_scales, rad_offsets, ref_scales, ref_offsets


def write_geotiff(out_path, data_arr, geo_transform, epsg_code=4326):
    """
    通用函数：将 Numpy 数组保存为带有坐标系的 GeoTIFF 文件。
    """
    # 自动识别数据维度 (单波段还是多波段)
    if len(data_arr.shape) == 3:
        nb, lines, samples = data_arr.shape
    else:
        lines, samples = data_arr.shape
        nb = 1
        data_arr = data_arr[np.newaxis, :, :]  # 统一转为三维方便处理

    # 自动匹配 GDAL 数据类型
    dtype = data_arr.dtype
    if dtype == np.int16:
        gdal_type = gdal.GDT_Int16
    elif dtype == np.uint8:
        gdal_type = gdal.GDT_Byte
    else:
        gdal_type = gdal.GDT_Float32

    # 创建驱动并写入
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(out_path, samples, lines, nb, gdal_type)
    ds.SetGeoTransform(geo_transform)

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg_code)
    ds.SetProjection(srs.ExportToWkt())

    for i in range(nb):
        ds.GetRasterBand(i + 1).WriteArray(data_arr[i])

    ds.FlushCache()
    ds = None
    print(f"[写入成功] -> {out_path}")

# ==========================================
# 模块 2：物理定标核心 (Physics Core)
# 数学、物理公式，从DN值到物理量的转换
# ==========================================

def calibrate_radiance(dn_arr, scales, offsets):
    """
    将原始 DN 值转换为光谱辐亮度 (Radiance)。
    公式: L = scale * (DN - offset)
    """
    if scales is None or offsets is None:
        return None

    # 保护机制：过滤无效值 (MODIS DN 有效范围通常为 0-32767)
    invalid_mask = (dn_arr < 0) | (dn_arr > 32767)

    # 定标计算 (利用 Numpy 广播机制处理 3D 数组)
    rad = (dn_arr - offsets[:, None, None]) * scales[:, None, None]

    # 将无效值设为 NaN，防止污染后续计算
    rad[invalid_mask] = np.nan
    return rad


def calibrate_reflectance(dn_arr, scales, offsets):
    """
    将原始 DN 值转换为大气顶层表观反射率 (TOA Reflectance)。
    公式: rho = scale * (DN - offset)
    """
    if scales is None or offsets is None:
        return None

    # 与辐射定标相同的逻辑，但使用的是反射率系数
    invalid_mask = (dn_arr < 0) | (dn_arr > 32767)
    ref = (dn_arr - offsets[:, None, None]) * scales[:, None, None]
    ref[invalid_mask] = np.nan
    return ref


def calculate_brightness_temperature(rad_arr, band_wl_um):
    """
    利用普朗克定律逆函数计算亮度温度 (Brightness Temperature)。
    输入: rad_arr (光谱辐亮度)
    返回: 开尔文 (K) 温度数组
    """
    # 普朗克定律相关辐射常数
    C1 = 1.4387685e4  # um·K
    C2 = 1.19104356e8  # um^5·W·m-2·sr-1

    bt = np.empty_like(rad_arr)
    for i, wl in enumerate(band_wl_um):
        L = rad_arr[i]

        # 使用 errstate 忽略除以 0 或对数计算中的负值警告
        # 产生这种极值时，Numpy 会自动将其置为 NaN
        with np.errstate(invalid='ignore', divide='ignore'):
            bt[i] = (C1 / wl) / np.log(1 + C2 / (wl ** 5) / L)

    return bt


# ==========================================
# 模块 3：几何校正引擎 (Geometry Engine)
# 支持前向十字填充与 GDAL Warp 后向投影
# ==========================================

def project_forward_cross_fill(data_arr, lat_arr, lon_arr, pixel_size=0.01):
    """
    方法 A：前向投影与十字填充（快速、保留原始数值）
    """
    lat_min, lat_max = np.nanmin(lat_arr), np.nanmax(lat_arr)
    lon_min, lon_max = np.nanmin(lon_arr), np.nanmax(lon_arr)

    lines = int(np.ceil((lat_max - lat_min) / pixel_size))
    samples = int(np.ceil((lon_max - lon_min) / pixel_size))

    # 初始化画布
    proj_data = np.zeros((lines, samples), dtype=data_arr.dtype)

    flat_lat = lat_arr.ravel()
    flat_lon = lon_arr.ravel()
    vals = data_arr.ravel()

    # 过滤无效经纬度或 NaN 数据
    valid = (~np.isnan(flat_lat)) & (~np.isnan(flat_lon)) & (~np.isnan(vals))
    lat_v = flat_lat[valid]
    lon_v = flat_lon[valid]
    val_v = vals[valid]

    r_idx = np.floor((lat_v - lat_min) / pixel_size).astype(int)
    c_idx = np.floor((lon_v - lon_min) / pixel_size).astype(int)
    mask = (r_idx >= 1) & (r_idx < lines - 1) & (c_idx >= 1) & (c_idx < samples - 1)

    r, c, v = r_idx[mask], c_idx[mask], val_v[mask]

    # 两步十字填充
    for t in (1, 2):
        if t == 1:
            proj_data[r - 1, c] = v
            proj_data[r + 1, c] = v
            proj_data[r, c - 1] = v
            proj_data[r, c + 1] = v
        else:
            proj_data[r, c] = v

    # 简单补洞
    for _ in range(3):
        mask0 = (proj_data == 0)
        proj_data[mask0] = np.roll(np.roll(proj_data, 1, axis=0), 1, axis=1)[mask0]

    # 纬度翻转 (纠正上北下南)
    proj_data = np.flip(proj_data, axis=0)  # 注意：对于 2D 数组，行轴是 axis=0

    # 生成 GeoTransform
    geo_transform = (lon_min, pixel_size, 0, lat_max, 0, -pixel_size)
    return proj_data, geo_transform


def project_backward_gdal_warp(data_arr, lat_arr, lon_arr, pixel_size=0.01, resample_alg=gdal.GRA_NearestNeighbour):
    """
    方法 B：使用 GDAL Warp 进行后向投影（内存/跨平台安全版）
    """
    rows, cols = data_arr.shape

    # 1. 自动匹配数据类型
    if data_arr.dtype == np.float32:
        gdal_type = gdal.GDT_Float32
    elif data_arr.dtype == np.int16:
        gdal_type = gdal.GDT_Int16
    else:
        gdal_type = gdal.GDT_Float32

    src_ds = gdal.GetDriverByName('MEM').Create('', cols, rows, 1, gdal_type)
    src_ds.GetRasterBand(1).WriteArray(data_arr)

    # 2. 将经纬度数组写入 GDAL 虚拟内存 (/vsimem/)
    lat_vsi = f'/vsimem/lat_{id(data_arr)}.tif'
    lon_vsi = f'/vsimem/lon_{id(data_arr)}.tif'

    driver_tiff = gdal.GetDriverByName('GTiff')

    ds_lat = driver_tiff.Create(lat_vsi, cols, rows, 1, gdal.GDT_Float32)
    ds_lat.GetRasterBand(1).WriteArray(lat_arr)
    ds_lat.GetRasterBand(1).SetNoDataValue(np.nan)  # 【修复1】明确告知 GDAL NaN 是无效数据
    ds_lat.FlushCache()

    ds_lon = driver_tiff.Create(lon_vsi, cols, rows, 1, gdal.GDT_Float32)
    ds_lon.GetRasterBand(1).WriteArray(lon_arr)
    ds_lon.GetRasterBand(1).SetNoDataValue(np.nan)  # 【修复1】明确告知 GDAL NaN 是无效数据
    ds_lon.FlushCache()

    # 3. 组装 GeoLoc 元数据
    # 生成 GDAL 期望的完整 WKT 投影字符串
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    wkt_srs = srs.ExportToWkt()

    # 3. 组装 GeoLoc 元数据
    geoloc_metadata = {
        'SRS': wkt_srs,  # 【完美修复】不再使用简写，使用完整的 WKT 字符串！
        'X_DATASET': lon_vsi,
        'X_BAND': '1',
        'Y_DATASET': lat_vsi,
        'Y_BAND': '1',
        'PIXEL_OFFSET': '0', 'LINE_OFFSET': '0', 'PIXEL_STEP': '1', 'LINE_STEP': '1'
    }
    src_ds.SetMetadata(geoloc_metadata, 'GEOLOCATION')

    # 【核心修复2】：手动计算经纬度边界（转换为纯 Python float 以防报错）
    lon_min, lon_max = float(np.nanmin(lon_arr)), float(np.nanmax(lon_arr))
    lat_min, lat_max = float(np.nanmin(lat_arr)), float(np.nanmax(lat_arr))

    # 配置 Warp 选项
    warp_options = gdal.WarpOptions(
        format='MEM', dstSRS='EPSG:4326',
        xRes=pixel_size, yRes=pixel_size,
        outputBounds=(lon_min, lat_min, lon_max, lat_max),
        resampleAlg=resample_alg, geoloc=True,
        srcNodata=0, dstNodata=0
    )

    # 4. 执行投影
    dest_ds = gdal.Warp('', src_ds, options=warp_options)
    proj_arr = dest_ds.GetRasterBand(1).ReadAsArray()
    geo_transform = dest_ds.GetGeoTransform()

    # 5. 释放并清理虚拟内存
    src_ds, dest_ds, ds_lat, ds_lon = None, None, None, None
    gdal.Unlink(lat_vsi)
    gdal.Unlink(lon_vsi)

    return proj_arr, geo_transform

def reproject_modis(data_arr, geo_info, method='forward', pixel_size=0.01):
    """
    几何校正主控函数：自由切换投影算法。
    """
    if len(data_arr.shape) == 3:
        nb = data_arr.shape[0]
        out_bands = []
        out_gt = None
        for i in range(nb):
            arr_2d, gt = reproject_modis(data_arr[i], geo_info, method, pixel_size)
            out_bands.append(arr_2d)
            if out_gt is None: out_gt = gt
        return np.stack(out_bands, axis=0), out_gt

    # 【统一接口提取】：现在两种方法都需要 Numpy 数组格式的经纬度
    lat_arr = geo_info['lat']
    lon_arr = geo_info['lon']

    if method == 'forward':
        return project_forward_cross_fill(data_arr, lat_arr, lon_arr, pixel_size)
    elif method == 'gdal_warp':
        # 【修改点】：不再传文件路径，直接传经纬度数组
        return project_backward_gdal_warp(data_arr, lat_arr, lon_arr, pixel_size)
    else:
        raise ValueError(f"未知的方法: {method}")


# ==========================================
# 模块 4：可视化合成 (Visualization)
# 负责生成 MNR 假彩色与 RGB 真彩色
# ==========================================

def generate_mnr_false_color(proj_data, band_indices=(23, 1, 0), max_val=3400):
    """
    生成 MNR 假彩色合成数组，专门用于凸显高温热源（火点）。
    - proj_data: 投影后的全波段三维数组 (Bands, Rows, Cols)
    - band_indices: (红通道波段索引, 绿通道波段索引, 蓝通道波段索引)
      注意：这里的索引需要根据传入 proj_data 的实际波段顺序来定！
    """
    r_idx, g_idx, b_idx = band_indices

    # 提取并转为浮点型以防计算溢出
    R = proj_data[r_idx].astype(np.float32)
    G = proj_data[g_idx].astype(np.float32)
    B = proj_data[b_idx].astype(np.float32)

    # 1. 红通道 (中红外) 特殊拉伸：99% 分位数截断
    low_q = 0.99
    valid_R = R[R > 0]  # 排除背景 0 值

    if valid_R.size > 0:
        a = float(np.quantile(valid_R, low_q))
        b = float(np.max(valid_R))  # 也可以用 1.0 的分位数
    else:
        a, b = 0, max_val

    # 经验阈值限制 (基于 MODIS 数据特性)
    if a < 3000: a = 3000.0
    if b > 3400: b = 3400.0
    if a >= b: b = a + 1.0

    # 截断与归一化
    R_clipped = np.clip(R, a, b)
    R_scaled = (R_clipped - a) / (b - a) * max_val
    R_scaled = np.where(R > 0, R_scaled, 0.0)  # 保持背景为 0
    R_out = np.clip(np.round(R_scaled), 0, max_val).astype(np.int16)

    # 2. 绿/蓝通道线性增强
    boost = 1.2
    G_out = np.clip(np.round(G * boost), 0, max_val).astype(np.int16)
    B_out = np.clip(np.round(B * boost), 0, max_val).astype(np.int16)

    # 将三个通道堆叠为 (3, Rows, Cols) 的数组
    return np.stack([R_out, G_out, B_out], axis=0)


def generate_true_color(proj_data, band_indices=(0, 3, 2)):
    """
    生成常规 RGB 真彩色数组。
    这里直接提取波段，保留原始定标值，便于后续在 GIS 软件中调节。
    """
    r_idx, g_idx, b_idx = band_indices
    R = proj_data[r_idx]
    G = proj_data[g_idx]
    B = proj_data[b_idx]

    return np.stack([R, G, B], axis=0)


# ==========================================
# 主控流水线 (Main Pipeline Execution)
# ==========================================
def process_modis_data(orbit_data_path, img_save_path, method='forward'):
    """
    重构后的 MODIS 预处理主函数
    - method: 'forward' (十字填充) 或 'gdal_warp' (标准后向投影)
    """
    os.makedirs(img_save_path, exist_ok=True)

    # 【模块 1：数据读取】找文件
    geo_files = search_modis_files(orbit_data_path)
    print(f"找到 {len(geo_files)} 个 MOD03 定位文件。")

    for geo_path in geo_files:
        platform, orbit_time, date_str, ad_type = parse_filename_info(geo_path)
        data_path = match_science_data(orbit_data_path, date_str)

        if not data_path:
            print(f"[{orbit_time}] 未找到对应的科学数据 MOD021KM，跳过。")
            continue

        print(f"\n开始处理: {platform} - {orbit_time} (使用 {method} 投影)")

        # 【模块 1：数据读取】读取原始数据
        sd_data = SD(data_path, SDC.READ)
        arr250, s250, o250, ref_s250, ref_o250 = read_hdf_dataset(sd_data, 'EV_250_Aggr1km_RefSB')
        arr500, s500, o500, ref_s500, ref_o500 = read_hdf_dataset(sd_data, 'EV_500_Aggr1km_RefSB')
        arr1kmR, s1kmR, o1kmR, ref_s1kmR, ref_o1kmR = read_hdf_dataset(sd_data, 'EV_1KM_RefSB')
        arr1kmE, s1kmE, o1kmE, _, _ = read_hdf_dataset(sd_data, 'EV_1KM_Emissive')
        sd_data.end()

        # 【模块 2：定标预处理】物理定标
        # 反射率校正 (可见光/近红外) -> 注意传入的是原始 arr
        ref250 = calibrate_reflectance(arr250, ref_s250, ref_o250)
        ref500 = calibrate_reflectance(arr500, ref_s500, ref_o500)
        ref1kmR = calibrate_reflectance(arr1kmR, ref_s1kmR, ref_o1kmR)

        # 辐射定标与亮温反演 (热红外)
        rad1kmE = calibrate_radiance(arr1kmE, s1kmE, o1kmE)
        wl_list = [3.750, 3.959, 3.959, 4.050, 4.4655, 4.5155, 6.715, 7.325, 8.550, 9.730, 11.03, 12.02, 13.335, 13.635,
                   13.935, 14.235]
        bt1kmE = calculate_brightness_temperature(rad1kmE, wl_list)

        # 准备合并数据进行投影 (缩放为 Int16)
        ref250_i = np.round(np.nan_to_num(ref250) * 10000).astype(np.int16)
        ref500_i = np.round(np.nan_to_num(ref500) * 10000).astype(np.int16)
        ref1kmR_i = np.round(np.nan_to_num(ref1kmR) * 10000).astype(np.int16)
        bt1kmE_i = np.round(np.nan_to_num(bt1kmE) * 10).astype(np.int16)

        full_swath_data = np.concatenate([ref250_i, ref500_i, ref1kmR_i, bt1kmE_i], axis=0)

        # 【模块 3：几何投影】投影重采样
        sd_geo = SD(geo_path, SDC.READ)

        # 1. 先把数据取出来
        lat_arr = sd_geo.select('Latitude').get().astype(np.float32)
        lon_arr = sd_geo.select('Longitude').get().astype(np.float32)

        # 2. 【新增修复】严格过滤无效的经纬度（去除填充值造成的极端异常）
        lat_arr[(lat_arr < -90) | (lat_arr > 90)] = np.nan
        lon_arr[(lon_arr < -180) | (lon_arr > 180)] = np.nan

        # 3. 存入字典
        geo_info = {
            'lat': lat_arr,
            'lon': lon_arr,
            'mod03_path': geo_path  # 用于 gdal_warp
        }
        sd_geo.end()

        # 核心路由：根据 method 切换底层算法！
        proj_all, geo_transform = reproject_modis(full_swath_data, geo_info, method=method, pixel_size=0.01)

        # 【模块 1 & 4：包装与入库】写出 GeoTIFF
        # 1. 存全波段
        out_all = os.path.join(img_save_path, f'{platform}_MODIS_1000M_{orbit_time}_{ad_type}_L2_ALL.tif')
        write_geotiff(out_all, proj_all, geo_transform)

        # 2. 存假彩色 MNR (假设新合成数组中波段23在中红外对应索引，这里需按你 concatenate 的顺序仔细核对索引！)
        # 需要确认 concatenate 后的实际索引，因为我需要中波段23所以给23赋红
        mnr_arr = generate_mnr_false_color(proj_all, band_indices=(23, 1, 0))
        out_mnr = os.path.join(img_save_path, f'{platform}_MODIS_1000M_{orbit_time}_{ad_type}_L2_MNR.tif')
        write_geotiff(out_mnr, mnr_arr, geo_transform)

        # 3. 存真彩色RGB ,band_indices=(0, 3, 2) 默认
        rgb_arr = generate_true_color(proj_all)
        out_rgb = os.path.join(img_save_path, f'{platform}_MODIS_1000M_{orbit_time}_{ad_type}_L2_RGB.tif')
        write_geotiff(out_rgb, rgb_arr, geo_transform)

        print(f"[{orbit_time}] 处理完成！")

if __name__ == '__main__':
    # HDF文件命名格式： AQUA_X_2026_03_02_14_15_A_G.MOD03/AQUA_X_2026_03_02_14_15_A_G.MOD021KM
    orbit_dir = r"Your input folder"
    result_dir = r"Your output folder"
    process_modis_data(orbit_dir, result_dir, method='gdal_warp')  # 或者换成 'gdal_warp'（传统商业软件结果）