# AQUA/MODIS 卫星影像高级预处理系统：算法理论基础文档 (ATBD)

**文档版本:** 1.0  
**适用卫星/传感器:** NASA AQUA (EOS PM-1) / MODIS  
**处理层级:** L1B $\rightarrow$ L2 (Geocoded, Radiometrically Calibrated)  

---

## 1. 算法概述 (Algorithm Overview)

本项目提供了一套面向 AQUA/MODIS L1B 级观测数据的自动化预处理算法。核心目标是将原始数字量化值（Digital Number, DN）精确转换为具有物理意义的参量（表观反射率、亮度温度），并解决宽扫描带传感器特有的空间几何畸变（（蝴蝶效应）蝴蝶结效应）。

本系统独创性地集成了**双引擎几何校正架构**，并针对**高温异常点（火情）**定制了极值截断假彩色合成策略，适用于科学研究与业务化灾害监测。

---

## 2. 输入数据规范 (Input Data Specifications)

算法依托基于时间戳的精准匹配机制，摄取以下两种标准的 HDF4 格式数据：

1. **科学数据集 (L1B Science Data, `MOD021KM` / `MYD021KM`)**
   * 提取 `EV_250_Aggr1km_RefSB`、`EV_500_Aggr1km_RefSB`、`EV_1KM_RefSB`（反射波段）以及 `EV_1KM_Emissive`（发射/热红外波段）。
2. **地理定位数据集 (Geolocation Data, `MOD03` / `MYD03`)**
   * 提取像元级经纬度坐标 (`Latitude`, `Longitude`)。

---

## 3. 物理定标与参数反演 (Physical Calibration & Inversion)

### 3.1 辐射定标 (Radiance Calibration)
对于热红外（发射）波段，需将 DN 值转换为光谱辐亮度。
* **物理模型：** 采用线性定标模型。
  $$L = \text{scale} \times (\text{DN} - \text{offset})$$
  *(式中 $L$ 为辐亮度，$W \cdot m^{-2} \cdot sr^{-1} \cdot \mu m^{-1}$)*
* **质量控制：** 严格屏蔽有效物理范围（0~32767）之外的异常背景值。

### 3.2 大气顶层反射率校正 (TOA Reflectance)
对于可见光与近红外波段，直接由原始 DN 值计算大气顶层（TOA）表观反射率。
* **物理模型：**
  $$\rho = \text{scale}_{ref} \times (\text{DN} - \text{offset}_{ref})$$
  *(注：HDF 属性中的定标系数已内含太阳天顶角余弦项的预校正)*

### 3.3 亮度温度反演 (Brightness Temperature)
基于普朗克定律（Planck's Law）的逆函数，由光谱辐亮度反演地表亮温。
* **物理模型：**
  $$T = \frac{C_1}{\lambda \cdot \ln\left(\frac{C_2}{\lambda^5 \cdot L} + 1\right)}$$
  *(式中 $T$ 为亮度温度，$\lambda$ 为中心波长，$C_1, C_2$ 为普朗克辐射常数)*
* **计算鲁棒性：** 引入 NumPy 的 `errstate` 上下文管理器，安全捕获并规避因深空背景导致 $L \le 0$ 时触发的对数域与除零崩溃。

---

## 4. 双引擎几何校正机制 (Dual-Engine Geometric Correction)

MODIS 因其扫描几何特性，在扫描带边缘会产生严重的“蝴蝶结效应（Bowtie Effect）”。本系统提供两种可自由切换的投影方案：

### 4.1 方法 A：前向启发式十字填充 (Forward Cross-Filling)
* **核心逻辑：** 建立目标网格（$0.01^\circ$ 分辨率），计算像元行列索引 $r_{idx}$ 与 $c_{idx}$。
* **两步填充策略：**
  1. **$T=1$ (膨胀):** 将观测值向上下左右邻域扩散，修复分辨率差异导致的采样黑洞。
  2. **$T=2$ (修正):** 重新覆盖中心点，确保最高空间几何精度的数据置于顶层。
* **方向纠正：** 由于 AQUA 卫星升轨采样的矩阵映射呈“南上北下”倒置，系统自动执行行轴翻转（`np.flip`），恢复标准地图投影（North-Up）。

### 4.2 方法 B：GDAL Warp 后向重采样 (Backward Resampling)
* **工业级标准：** 调用 GDAL 底层 C++ 引擎进行重采样，获取平滑的制图级结果。
* **无盘化虚拟内存技术：** 将经纬度矩阵写入 GDAL 虚拟文件系统（`/vsimem/`），以 `GEOLOCATION` 数组形式传参。**该技术解决了因 GDAL 缺失 HDF4 驱动或 Windows 路径解析冲突导致的数据读取失败问题。**
