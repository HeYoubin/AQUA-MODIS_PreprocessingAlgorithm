# 🌍 AQUA/MODIS 卫星影像高级预处理算法

![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue.svg)
![GDAL](https://img.shields.io/badge/GDAL-3.0%2B-green.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)

## 📖 项目简介
本项目是一个基于 Python 和 GDAL 开发的 AQUA/MODIS L1B 卫星数据自动化预处理系统。它不仅实现了标准的辐射定标、反射率校正与亮温反演，还针对 MODIS 宽扫描带数据特有的“蝴蝶效应（蝴蝶结效应）”提供了**双引擎几何校正方案**。

本算法特别针对**高温异常点（如森林火灾）监测**优化了假彩色合成算法，是遥感业务化运行与科研算法探索的理想框架。

## ✨ 核心特性

- **🚀 双引擎投影架构 (Dual-Engine Projection):**
  - `forward`: 极速前向投影与启发式十字填充，有效消除边缘蝴蝶结空洞，保留原始科学数值（适用于快速浏览与定制算法）。
  - `gdal_warp`: 基于 GDAL 底层 C++ 引擎的后向投影重采样，利用虚拟内存（`/vsimem/`）实现无盘化安全投影，规避传统 HDF4 驱动缺失问题。
- **🛡️ 内存与极值安全防线:** 彻底修复了全范围网格生成时的内存溢出（OOM）Bug，以及 GDAL 乘法溢出（Overflow）问题，保障系统工业级稳定性。
- **🔥 业务级可视化合成:**
  - **MNR 假彩色:** 利用 99% 分位数截断策略拉伸中红外波段，精准凸显高温火点（红色），同时区分云层（白色）与植被（绿色）。
  - **RGB 真彩色:** 自动合成标准化自然色彩底图。
  - 用户可自行补充部分，针对彩色合成进行自定义拉伸，使得预处理影像清晰，可视化效果增强。

## 📂 项目结构

* 📁 **AQUA-MODIS_PreprocessingAlgorithm/** (项目根目录)
  * 📁 **src/** (源代码目录)
    * 📄 `modis_processor.py` —— 核心处理脚本（包含 IO、定标、几何引擎与可视化四大解耦模块）
  * 📁 **docs/** (文档目录)
    * 📄 `ATBD.md` —— 算法理论基础文档
    * 📄 `learning_notes.md` —— 编写算法的记录博客
  * 📁 **sample_output/** (示例输出)
    * 🖼️ `*_L2_MNR.png` —— 假彩色火点增强示例图
    * 🖼️ `*_L2_RGB.png` —— 真彩色示例图
  * 📄 `README.md` —— 项目说明文档


## 🛠️ 快速开始

### 1. 环境依赖
推荐使用 Conda 部署环境，以解决空间库冲突问题：
\`\`\`bash
conda create -n modis_env python=3.9 gdal pyhdf numpy scikit-image
conda activate modis_env
\`\`\`

### 2. 运行代码
修改 `modis_processor.py` 底部的路径配置，即可一键运行：
\`\`\`python
orbit_dir = r"你的原始MODIS_HDF数据目录"
result_dir = r"输出目录"

# 自由切换算法：'forward' 或 'gdal_warp'
process_modis_data(orbit_dir, result_dir, method='gdal_warp') 
\`\`\`

## 📚 算法理论文档 (ATBD)
如果您对底层的辐射定标公式、普朗克定律反演或十字填充算法感兴趣，请查阅本项目配套的：
👉 **[算法理论基础文档 (ATBD)](./docs/ATBD.md)** <br>
如果您对编写算法的过程中我遇到的坑和我的一些思考感兴趣，请查阅本项目配套的：
👉 **[AQUA/MODIS卫星数据预处理博客 (Blog)](./docs/learning_notes.md)**

---
*Developed as an open-source initiative for remote sensing algorithm exploration.*