# TCGA_SKCM_Basophil_Analysis
Transcriptome and  clinical analysis of basophil in SKCM
# TCGA-SKCM 嗜碱性粒细胞预后与转录组分析
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# 项目简介
本研究基于TCGA-SKCM队列，通过特征基因集计算嗜碱性粒细胞浸润度，构建LASSO-Cox预后模型，并联合全转录组相关性分析与GO富集，揭示嗜碱性粒细胞在黑色素瘤中的保护作用及免疫机制。

# 数据获取
**注意：由于TCGA数据使用协议限制，本仓库不包含原始数据文件。**
请自行从GDC Portal下载以下文件，并放置于脚本同目录：
| **文件名** | **下载地址** | **说明** |
|:---|:---|:---|
| `TCGA-SKCM.star_counts.tsv.gz` | [GDC Portal](https://portal.gdc.cancer.gov/) | STAR计数数据 |
| `TCGA-SKCM.star_tpm.tsv.gz` | [GDC Portal](https://portal.gdc.cancer.gov/) | STAR-TPM表达数据 |
| `TCGA-SKCM.clinical.tsv.gz` | [GDC Portal](https://portal.gdc.cancer.gov/) | 临床随访数据 |

**下载步骤**：
1. 访问 [GDC Data Portal](https://portal.gdc.cancer.gov/)
2. 选择 **Repository** → 筛选 **Cases** → **TCGA-SKCM**
3. 选择 **Files** → 筛选 **Data Category**: `Transcriptome Profiling` → **Experimental Strategy**: `RNA-Seq`
4. 下载 `star_counts` 和 `star_tpm` 文件
5. 选择 **Clinical** → 下载临床数据

# 使用方法

1. 下载TCGA-SKCM数据文件至脚本同目录：
   - `TCGA-SKCM.star_counts.tsv.gz`
   - `TCGA-SKCM.star_tpm.tsv.gz`
   - `TCGA-SKCM.clinical.tsv.gz`
2. 运行 `TCGA_SKCM_Basophil_Analysis.R`

# 依赖环境
- 本研究基于R4.5.3版本
- 包："tidyverse", "glmnet", "survival", "survminer", "ggplot2", "ggrepel", "rms", "Hmisc", "survivalROC", "clusterProfiler", "org.Hs.eg.db"

# 主要发现
- **临床**：低风险组生存率显著更高（P < 0.0001）
- **转录组**：发现2079个正相关基因，97个负相关基因
- **功能**：富集于白细胞趋化、迁移、免疫激活通路（P达10⁻²⁹）
- **分子**：CPA3、CCL2、CSF2RB等关键基因验证

# 作者

GitHub: [@qiao-yuhui](https://github.com/qiao-yuhui)
