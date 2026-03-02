# LiverCEST_MRI_PostProcessing
人体肝脏CEST磁共振图像后处理MATLAB代码包，支持fasting/meal两类CEST数据的洛伦兹拟合、ROI分析、多层批量处理，适配0.7uTneg/0.7uTpos/2uT三种数据类型。
- [功能特性](#功能特性)
- [环境要求](#环境要求)
- [文件夹结构](#文件夹结构)
- [快速开始](#快速开始)
- [参数配置](#参数配置)
- [运行流程](#运行流程)
- [结果输出](#结果输出)
- [注意事项](#注意事项)
- [作者信息](#作者信息)
## ✨ 功能特性
- 支持空腹(fasting)/餐后(meal)两类数据处理
- 适配0.7uTneg/0.7uTpos/2uT三种频率维度数据
- 支持单切片手动ROI勾画 + 多层批量自动处理
- 提供洛伦兹拟合（LD）的约束/无约束两种模式
- 输出定量图谱（MTRrex）、统计结果、掩膜文件等多类结果
## 🛠️ 环境要求
- MATLAB R2020b 及以上版本
- 需安装MATLAB工具箱：
  - 图像处理工具箱 (Image Processing Toolbox)
  - 曲线拟合工具箱 (Curve Fitting Toolbox)
  - 统计和机器学习工具箱 (Statistics and Machine Learning Toolbox)
## 📂 文件夹结构
代码包需严格遵循以下目录结构（仅根路径可自定义）：
```
LiverCEST_MRI_PostProcessing/
├── 数据/                # 原始扫描数据
│   ├── fasting/         # 空腹数据
│   └── meal/            # 餐后数据
├── 结果/                # 分析结果输出
│   ├── fasting/         # 空腹结果
│   └── meal/            # 餐后结果
└── Process_LiverMainFunction0301.m  # 主执行代码
```
## 🚀 快速开始
### 1. 路径配置
将代码包整个文件夹添加到MATLAB路径：
```matlab
% MATLAB命令行执行
addpath(genpath('你的代码包根路径'));
```
### 2. 基础参数修改
在主代码中修改核心路径和数据类型：
```matlab
% 数据类型选择：0=0.7utpos/1=0.7utneg/2=2uT
NegOrPos = 1;
% 数据根路径（根据实际修改）
datapath = 'C:\LiverCEST\数据\fasting';
% 结果输出路径（根据实际修改）
resultpath = 'C:\LiverCEST\结果\fasting';
```
### 3. 首次运行（单切片）
```matlab
% 单切片配置
N_slice = 21;          % 分析中间切片（示例：21层）
use_exist_ROI = 0;     % 首次运行重新勾画ROI
multi_slice = 0;       % 单切片模式
```
运行代码后，按提示完成：
1. 绘制三角形去噪 → 2. 勾画肝脏轮廓 → 3. 勾画ROI区域
## ⚙️ 参数配置
### 核心必配参数
| 参数名 | 取值范围 | 说明 |
|--------|----------|------|
| `NegOrPos` | 0/1/2 | 0=0.7utpos，1=0.7utneg，2=2uT |
| `mode` | 0/1 | 0=LD无约束（推荐测试），1=LD加约束 |
| `datapath` | 字符串 | 原始数据根路径 |
| `resultpath` | 字符串 | 结果输出根路径 |
| `DataFormat` | 0/1 | 0=bruker格式，1=dicom格式（默认） |
| `nROI` | 1~4 | 需分析的ROI数量（最多4个） |

### 切片控制参数
| 参数名 | 取值 | 说明 |
|--------|------|------|
| `N_slice` | 数值/0 | 单切片填具体层数，多层批量填0 |
| `multi_slice` | 0/1 | 0=单切片，1=多层批量 |
| `use_exist_ROI` | 0/1 | 0=重新勾画ROI，1=使用已有ROI（多层必选） |
| `Thmask_datapath` | 字符串 | 多层模式下，指向单切片生成的mask路径 |

### 进阶参数（默认无需修改）
| 参数名 | 默认值 | 说明 |
|--------|--------|------|
| `register` | 0 | 配准模式：0=不配准，1=刚体配准，2=RPC |
| `correct` | 2 | 校正模式：1=B0校正，2=自校正，3=2-Pools LD校正 |
| `snr` | 15 | 信噪比阈值 |

## 📝 运行流程
### 步骤1：单切片分析（首次运行）
1. 配置单切片参数（`multi_slice=0`、`use_exist_ROI=0`）
2. 运行代码，完成交互式ROI勾画
3. 确认`resultpath`下生成mask文件（如`Thmask.mat`）

### 步骤2：多层批量分析
1. 复制单切片mask路径（示例）：
   ```matlab
   Thmask_datapath = 'C:\LiverCEST\结果\fasting\no_lb_neg_interp_0227\slice21';
   ```
2. 修改批量参数：
   ```matlab
   N_slice = 0;          % 多层模式
   multi_slice = 1;      % 启用多层批量
   use_exist_ROI = 1;    % 使用已有ROI
   total_slice_num = 41; % 数据总层数
   ```
3. 运行代码，自动处理所有切片

## 📊 结果输出
运行完成后，`resultpath`对应目录下生成以下核心文件：
| 文件类型 | 说明 |
|----------|------|
| `liverMask.mat/fig` | 肝脏掩膜文件/可视化图 |
| `Thmask.mat` | ROI掩膜文件 |
| `Lorentzian fits.mat` | 洛伦兹拟合结果 |
| `MTRrex_map.mat` | MTRrex定量图谱 |
| `sliceXX_mean_neg3p5.mat` | 切片统计结果 |
| `sliceXX LD NOE-3.5.fig` | NOE可视化图 |

## ⚠️ 注意事项
1. ROI勾画不可过小，否则会触发程序报错
2. 2uT数据需重新计算背景频率索引（默认适配0.7ut）：
   ```matlab
   % 2uT数据需修改此处索引值
   background_frequency = 101; % 示例值，需按实际计算
   ```
3. 多层分析前必须完成单切片分析，确保`Thmask_datapath`路径正确
4. 数据/结果文件夹的子目录结构必须完全一致

## ✍️ 作者信息
- 作者：Li Kaixiang
- 版本更新时间：20260302
```
