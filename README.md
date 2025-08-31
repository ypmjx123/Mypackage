# Mypackage: 生物信息数据分析工具包

## 项目简介
Mypackage是一个全面的生物信息数据分析R包，专注于癌症基因组学、表达谱分析和临床数据整合分析。该包提供了一系列高效的工具函数，帮助研究人员快速处理、分析和可视化各种生物信息学数据。

## 主要功能和特点

### 1. 突变数据分析
- MAF格式数据处理和分析
- 突变与临床数据相关性分析
- 突变网络构建与可视化
- 突变与表达数据整合分析

### 2. 表达谱分析
- 基因表达差异分析
- GSEA富集分析与可视化
- 免疫相关基因表达分析与可视化

### 3. 临床数据与生存分析
- 生存曲线绘制
- COX森林图分析
- 临床特征与分子特征关联分析

### 4. 基因组不稳定性分析
- HRD(同源重组缺陷)评分计算
- CNV(拷贝数变异)和SV(结构变异)获取

### 5. 数据可视化工具
- 多种热图绘制功能
- 通路与突变整合可视化
- 免疫浸润可视化

## 安装指南

### 从本地安装
```R
# 安装devtools包（如果尚未安装）
install.packages("devtools")

# 从本地安装mypackage
devtools::install_local("d:/R/libs/mypackage")
或install.packages("./Mypackage_1.3.1.tar.gz", repos = NULL, type = "source")

#在线安装mypackage
remotes::install_git("https://gitee.com/YPM2022/mypackage")

# 加载包
library(Mypackage)
```

### 依赖包
mypackage依赖于以下R包，安装时会自动安装：
- ggplot2: 数据可视化
- maftools: MAF文件处理与分析
- survival & survminer: 生存分析
- limma: 表达差异分析
- dplyr: 数据处理
- tidyr: 数据整理
- stringr: 字符串处理
- biomaRt: 基因注释
- pheatmap: 热图绘制

## 内置数据集

mypackage包含多个内置数据集，便于用户学习和测试功能：

| 数据集名称 | 描述 |
|----------|------|
| `mutation_CRC` | 结直肠癌突变数据 |
| `clin_TCGA` | TCGA临床数据 |
| `exp_raw` | 基因表达原始数据 |
| `Gene_group_CRC1` | CRC相关基因集 |
| `pathway_data` | 通路相关数据 |
| `hotsgenes` | 热点基因列表 |
| `tumor_ploidy` | 肿瘤倍性数据 |

## 主要功能模块

### 1. 突变数据分析模块

#### maf_cor函数
分析突变数据与临床特征之间的相关性。
```R
# 示例用法
result <- maf_cor(mutation_data = mutation_CRC, clin = clin_TCGA, top=20,corrplot_method = "circle")
```

#### mut_network函数
构建突变基因网络。
```R
# 示例用法
network <- mut_network(SNV=mutation_CRC,top=20,pValue=0.01,customdata=NULL)
```

#### mutToMAF函数
将自定义突变格式转换为MAF格式。
```R
# 示例用法
maf_obj <- mutToMAF(root_dir=root_dir,clin=clindata,tumor_t=10,site_depth=100,hotspot_vaf=0.009,
                 non_hotspot_vaf=0.045,hotspotloss_vaf=0.095,non_hotspotloss_vaf=0.195)
```

### 2. 表达谱分析模块

#### limma.dif.visual函数
差异表达分析与可视化。
```R
# 示例用法
dif_results <- limma.dif.visual(exprdata=exp_CRC[,-1],
                           pdata=Gene_group_CRC1,datatype="TPM",
                           Type=c("Wild", "Mut"),diff_method="limma",
                           contrastfml= "Wild - Mut",
                           tumor="CRC TCGA",
                           P.Value=0.05,
                           logFC=0.5,tidyHeatmap=TRUE,
                           color= NULL,
                           ann_colors = list(regulate = c(Down = "#1B9E77", Up = "#D95F02"),
                           PREX2 = c(Wild = "#757575", Mut = "#FF4040")),
                           Regulate=c("Up","Down"),GO=TRUE,GO.plot="dotplot",split=TRUE,
                           KEGG=TRUE,KEGG.plot="dotplot",rel_heights= c(1.5, 0.5, 1))
```

#### GSEAplot2函数
GSEA富集分析与高级可视化。
```R
# 示例用法
gsea_plot <- GSEAplot2(x = gseaResult)
```

### 3. 生存分析模块

#### ggsurvplots函数
生存曲线绘制与比较。
```R
# 示例用法
surv_plot <- ggsurvplotsggsurvplots(data = clin_TCGA, conf.int = FALSE,time_col = "PFS_MONTHS",
                 status_col = "PFS_STATUS", group_col = "Status", pvalue_table = TRUE,
                 palette = ggsci::pal_ucscgb()(4), risk.table = FALSE, title = NULL,
                 legend.labs = c("no KRAS or TP53", "TP53", "KRAS", "KRAS&TP53"),
                 xlab = "PFS_MONTHS", ylab = "Survival probability",
                 surv.median.line="hv",surv.scale="default",legend=FALSE)
```

#### cox_forest函数
COX回归森林图绘制。
```R
# 示例用法
cox = cox_forest(data=aa,
                          time_col = "OS_MONTHS",
                          status_col = "OS_STATUS",
                         Univariate=T,
                         univar_predictors=colnames(aa)[c(6:7,18:21,30,24,25,33)],
                         Multivariate=T,
                         multivar_predictors = colnames(aa)[c(6:7,18:21,30,24,25,33)],
                        show_plots = T,xticks1=NULL,#c(0,0.25,0.5,0.75,1.00,1.25,1.5,6.5,11),
                        xticks2=NULL,#c(0,0.25,0.5,0.75,1.00,2,2.5,6,15),
                        title_univar = "OS Univariate",
                        title_multivar = "OS Multivariate",
                        use_baseline_table = TRUE,all =T,forestplot=F,ci_pch=16,ci_col="red",ci_line="blue",zero_col="#e22e2a",
                        log2=T,footnote=NULL)
```

### 4. 基因组不稳定性分析模块

#### HRDscore函数
计算HRD(同源重组缺陷)评分。
```R
# 示例用法
hrd_scores <- HRDscore(dirpath = "d:/R/libs/mypackage/inst/HRD/")
```

#### CNV_SV_get函数
获取和分析CNV与SV数据。
```R
# 示例用法
cnv_sv_data <- CNV_SV_get(data=data,Type=c("重排/融合","扩增"))
```

### 5. 可视化工具模块

#### sig_Heatmap函数
特征热图绘制。
```R
# 示例用法
heatmap <- sig_Heatmap(input = input, features = feas,ID ="SAMPLE_ID",show_plot=F,
            condiction=condiction,id_condiction=colnames(condiction)[[1]],col_condiction=colnames(condiction)[[2]],
             cols_group=c("#757575","#FF4040"),row_group=c("red","green"),
             legend_show=TRUE,column_title_size=10,row_title_size=8,
             heatmap_col=NULL,
             #heatmap_col=c("#0505FA", "#FFFFFF", "#FA050D"),
             group = "PREX2",row_title="Regulate", scale = TRUE,name="Expression")
```

#### immu_visual函数
免疫浸润可视化。
```R
# 示例用法
immu_plot <- immu_visual(im=NULL,exp=exp_CRC[,-1],
                    method = 'epic',
                    sample_group=Gene_group_CRC1,
                    tumor="CRC TCGA",heatmap=TRUE,
                    Type=c("Wild", "Mut"),
                    color=c("#757575", "#FF4040"),
                    geom_text=TRUE,
                    test = "wilcox.test")
```

## 完整使用示例

以下是一个完整的分析流程示例，展示如何使用mypackage进行从数据处理到可视化的全流程分析：

```R
# 加载包和数据
library(Mypackage)
data(mutation_CRC)
data(clin_TCGA)
data(exp_raw)
# 1. 突变数据分析
mut_cor_result <- maf_cor(mutation_data = mutation_CRC, clin = clin_TCGA, gene = "TP53")
print(mut_cor_result)

# 2. 表达差异分析
dif_genes <- limma.dif.visual(exp_matrix = exp_raw, group_info = sample_group)

# 3. GSEA分析
gsea_results <- GSEAplot2(gene_list = dif_genes$logFC, pathway_db = "KEGG")

# 4. 生存分析
survival_plot <- ggsurvplots(surv_data = clin_TCGA, gene_expr = exp_raw, gene = "TP53")

# 5. 综合可视化
heatmap_plot <- sig_Heatmap(expr_matrix = exp_raw, group_info = sample_group, top_genes = 50)
```

## 项目结构

mypackage的项目结构组织如下：

```
mypackage/
├── R/                # R函数源码文件
├── data/             # 内置数据集
├── inst/             # 额外的安装文件（HRD数据、示例数据等）
├── man/              # 帮助文档
├── vignettes/        # 长篇文档和教程
├── DESCRIPTION       # 包元数据和依赖信息
├── NAMESPACE         # 命名空间定义
└── README.md         # 项目说明文档（当前文件）
```

## 函数列表

mypackage包含以下主要函数：

| 函数名 | 描述 |
|-------|------|
| `maf_cor` | 突变与临床数据相关性分析 |
| `mut_network` | 突变网络构建 |
| `mutToMAF` | 转换为MAF格式 |
| `limma.dif.visual` | 差异表达分析与可视化 |
| `GSEAplot2` | GSEA富集分析 |
| `ggsurvplots` | 生存曲线绘制 |
| `cox_forest` | COX森林图分析 |
| `HRDscore` | HRD评分计算 |
| `CNV_SV_get` | CNV和SV数据获取(非分析) |
| `sig_Heatmap` | 特征热图绘制 |
| `immu_visual` | 免疫浸润可视化 |
| `gene.mut_exp` | 基因突变与表达整合分析 |
| `path_mut_visual` | 通路与突变整合可视化 |
| `get_annotation` | 基因注释获取 |
| `exp_geneIDtoSYMBOL` | 表达数据ID转换 |

## 文档和帮助

### 查看帮助文档
```R
# 查看特定函数的帮助
doc <- help("maf_cor")

# 查看包的总体帮助
doc <- help(package = "Mypackage")
```

### 浏览 vignettes
```R
# 查看所有可用的vignettes
browseVignettes("Mypackage")

# 打开特定的vignette
vignette("Mypackge")
```

## 依赖包

Mypackage依赖于以下R包：
- R (>= 4.4.3)
- ggplot2
- maftools
- survival
- survminer
- limma
- dplyr
- tidyr
- stringr
- biomaRt
- pheatmap
- gridExtra
- RColorBrewer

## 贡献指南

我们欢迎社区贡献，如果您想为mypackage做出贡献，请按照以下步骤：

1. Fork 本仓库
2. 创建您的特性分支 (`git checkout -b feature/AmazingFeature`)
3. 提交您的更改 (`git commit -m 'Add some AmazingFeature'`)
4. 推送到分支 (`git push origin feature/AmazingFeature`)
5. 开启一个Pull Request

## 参考文献
1.  Yao J, Sun Q, Wu H, et al. Decoding the molecular landscape: HER2 and PD-L1 in advanced gastric cancer. Front Immunol. 2025;16:1567308. <doi:10.3389/fimmu.2025.1567308>
2.  Qiu Q, Tan D, Chen Q, et al. Clinical implications of PD-L1 expression and pathway-related molecular subtypes in advanced Asian colorectal cancer patients. Am J Cancer Res. 2024;14(2):796-808. <doi:10.62347/FSSF9938>
3.  Ding W, Yang P, Zhao X, et al. Unraveling EGFR-TKI resistance in lung cancer with high PD-L1 or TMB in EGFR-sensitive mutations. Respir Res. 2024;25(1):40. <doi:10.1186/s12931-023-02656-3>
4.  Peng H, Ying J, Zang J, et al. Specific Mutations in APC, with Prognostic Implications in Metastatic Colorectal Cancer. Cancer Res Treat. 2023;55(4):1270-1280. <doi:10.4143/crt.2023.415>
5.  Jiang Y, Mai G, Zhao X, et al. Molecular characterization and prognostic implications of KRAS mutations in pancreatic cancer patients: insights from multi-cohort analysis. NPJ Precis Oncol. 2025;9(1):299. Published 2025 Aug 22.      <doi:10.1038/s41698-025-01087-1>
5.  Cancer Genome Atlas Network. Comprehensive molecular characterization of human colon and rectal cancer. Nature. 2012;487(7407):330-337. Published 2012 Jul 18. <doi:10.1038/nature11252>

## 许可证

本项目采用MIT许可证 - 详见LICENSE文件

## 联系信息

如有问题或建议，请联系：
- 邮箱：[ypm_yy@outlook.com]
- Gitee：[https://github.com/ypmjx123/Mypackage.git]
