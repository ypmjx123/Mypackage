# Mypackage: A Bioinformatic Data Analysis Toolkit

## Project Overview
Mypackage is a comprehensive bioinformatic data analysis R package focusing on cancer genomics, expression profile analysis, and integrated analysis of clinical data.  This package provides a series of efficient tool functions to help researchers quickly process, analyze, and visualize various bioinformatic data.

## Main Functions and Features

### 1. Mutation Data Analysis
- MAF format data processing and analysis
- Correlation analysis between mutation and clinical data
- Mutation network construction and visualization
- Integrated analysis of mutation and expression data

### 2. Expression Profile Analysis
- Differential gene expression analysis
- GSEA enrichment analysis and visualization
- Visualization of immune-related gene expression analysis

### 3. Clinical Data and Survival Analysis
- Survival curve plotting
- COX forest plot analysis
- Correlation analysis between clinical and molecular features

### 4. Genomic Instability Analysis
- HRD (Homologous Recombination Deficiency) score calculation
- CNV (Copy Number Variation) and SV (Structural Variation) acquisition

### 5. Data Visualization Tools
- Multiple heatmap plotting functions
- Integrated visualization of pathways and mutations
- Immune infiltration visualization

## Installation Guide

### Local Installation
```R
# Install the devtools package (if not already installed)
install.packages("devtools")

# Install mypackage from local directory
devtools::install_local("d:/R/libs/mypackage")
or install.packages("./Mypackage_1.3.1.tar.gz", repos = NULL, type = "source")

# Install mypackage online
remotes::install_github(ypmjx123/Mypackage")

# Load the package
library(Mypackage)
```

### Dependent Packages
Mypackage depends on the following R packages, which will be automatically installed during installation:
-ggplot2: Data visualization
-maftools: MAF file processing and analysis
-survival & survminer: Survival analysis
-limma: Differential expression analysis
-dplyr: Data manipulation
-tidyr: Data tidying
-stringr: String manipulation
-biomaRt: Gene annotation
-pheatmap: Heatmap plotting

## Built-in Datasets

Mypackage includes several built-in datasets for users to learn and test functions:：

| Dataset Name | Description |
|----------|------|
| `mutation_CRC` | Colorectal cancer mutation data |
| `clin_TCGA` | TCGA clinical data |
| `exp_raw` | Raw gene expression data |
| `Gene_group_CRC1` | CRC-related gene set |
| `pathway_data` | Pathway-related data |
| `hotsgenes` | Hotspot gene list |
| `tumor_ploidy` | Tumor ploidy data |

## Main Function Modules

### 1. Mutation Data Analysis Module

#### maf_cor
Analyzes the correlation between mutation data and clinical characteristics.
```R
# Example usage
result <- maf_cor(mutation_data = mutation_CRC, clin = clin_TCGA, top=20,corrplot_method = "circle")
```

#### mut_network
Constructs a mutation gene network.
```R
# Example usage
network <- mut_network(SNV=mutation_CRC,top=20,pValue=0.01,customdata=NULL)
```

#### mutToMAF
Converts custom mutation format to MAF format.
```R
# Example usage
maf_obj <- mutToMAF(root_dir=root_dir,clin=clindata,tumor_t=10,site_depth=100,hotspot_vaf=0.009,
                 non_hotspot_vaf=0.045,hotspotloss_vaf=0.095,non_hotspotloss_vaf=0.195)
```

### 2. Expression Profile Analysis Module

#### limma.dif.visual
Differential expression analysis and visualization.
```R
# Example usage
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

#### GSEAplot2
GSEA enrichment analysis and advanced visualization.
```R
# Example usage
gsea_plot <- GSEAplot2(x = gseaResult)
```

### 3. Survival Analysis Module

#### ggsurvplots
Survival curve plotting and comparison.
```R
# Example usage
surv_plot <- ggsurvplotsggsurvplots(data = clin_TCGA, conf.int = FALSE,time_col = "PFS_MONTHS",
                 status_col = "PFS_STATUS", group_col = "Status", pvalue_table = TRUE,
                 palette = ggsci::pal_ucscgb()(4), risk.table = FALSE, title = NULL,
                 legend.labs = c("no KRAS or TP53", "TP53", "KRAS", "KRAS&TP53"),
                 xlab = "PFS_MONTHS", ylab = "Survival probability",
                 surv.median.line="hv",surv.scale="default",legend=FALSE)
```

#### cox_forest
COX regression forest plot plotting.
```R
# Example usage
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

### 4. Genome instability analysis module

#### HRDscore
Calculate the HRD (homologous recombination defect) score.
```R
# Example usage
hrd_scores <- HRDscore(dirpath = "d:/R/libs/mypackage/inst/HRD/")
```

#### CNV_SV_get
Acquire and analyze CNV and SV data.
```R
# Example usage
cnv_sv_data <- CNV_SV_get(data=data,Type=c("rearrangement/fusion","amplification"))
```

### 5. Visual tools module

#### sig_Heatmap
Feature heat map drawing.
```R
# Example usage
heatmap <- sig_Heatmap(input = input, features = feas,ID ="SAMPLE_ID",show_plot=F,
            condiction=condiction,id_condiction=colnames(condiction)[[1]],col_condiction=colnames(condiction)[[2]],
             cols_group=c("#757575","#FF4040"),row_group=c("red","green"),
             legend_show=TRUE,column_title_size=10,row_title_size=8,
             heatmap_col=NULL,
             #heatmap_col=c("#0505FA", "#FFFFFF", "#FA050D"),
             group = "PREX2",row_title="Regulate", scale = TRUE,name="Expression")
```

#### immu_visual
Immune infiltration visualization.
```R
# Example usage
immu_plot <- immu_visual(im=NULL,exp=exp_CRC[,-1],
                    method = 'epic',
                    sample_group=Gene_group_CRC1,
                    tumor="CRC TCGA",heatmap=TRUE,
                    Type=c("Wild", "Mut"),
                    color=c("#757575", "#FF4040"),
                    geom_text=TRUE,
                    test = "wilcox.test")
```

## Complete usage example

Here is an example of a complete analysis process that shows how to use mypackage for the full range of analysis from data processing to visualization:

```R
# Load package and data
library(Mypackage)
data(mutation_CRC)
data(clin_TCGA)
data(exp_raw)
# 1. Mutagenesis data analysis
mut_cor_result <- maf_cor(mutation_data = mutation_CRC, clin = clin_TCGA, gene = "TP53")
print(mut_cor_result)

# 2. Expression difference analysi
dif_genes <- limma.dif.visual(exp_matrix = exp_raw, group_info = sample_group)

# 3. GSEA analyse
gsea_results <- GSEAplot2(gene_list = dif_genes$logFC, pathway_db = "KEGG")

# 4.survival analysis
survival_plot <- ggsurvplots(surv_data = clin_TCGA, gene_expr = exp_raw, gene = "TP53")

# 5. Comprehensive visualization
heatmap_plot <- sig_Heatmap(expr_matrix = exp_raw, group_info = sample_group, top_genes = 50)
```

## Project structure

The project structure of Mypackage is organized as follows:

```
mypackage/
├── R/                # R Function source file
├── data/             # Built-in data set
├── inst/             # Additional installation files (HRD data, sample data, etc.)
├── man/              # Help document
├── vignettes/        # Long documents and tutorials
├── DESCRIPTION       # Metadata and dependency information
├── NAMESPACE         # Definition of namespaces
└── README.md         # Project Description Document
```

## Function list

Mypackage contains the following main functions:

| function name | characterization |
|-------|------|
| `maf_cor` | Correlation analysis between mutation and clinical data |
| `mut_network` | Construction of mutation networks |
| `mutToMAF` | Convert to MAF format |
| `limma.dif.visual` | Differential expression analysis and visualization |
| `GSEAplot2` | GSEA enrichment analysis |
| `ggsurvplots` | Survial curve mapping |
| `cox_forest` | COX forest diagram analysis |
| `HRDscore` | Calculation of the HRD score |
| `CNV_SV_get` | CNV and SV data acquisition (non-analysis) |
| `sig_Heatmap` | Feature heat map drawing |
| `immu_visual` | Immune infiltration visualization |
| `gene.mut_exp` | Analysis of gene mutation and expression integration |
| `path_mut_visual` | Visualization of pathway and mutation integration |
| `get_annotation` | Gene annotation acquisition |
| `exp_geneIDtoSYMBOL` | Express data ID conversion |

## Documentation and help

### Check out the help documentation
```R
# View help for specific functions
doc <- help("maf_cor")

# View the overall help for packages
doc <- help(package = "Mypackage")
```

### Browse vignettes
```R
# View all available vignettes
browseVignettes("Mypackage")

# Open a specific vignette
vignette("Mypackge")
```

## Dependent package

Mypackage depends on the following R packages:
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

## Contribution guide

We welcome contributions from the community. If you want to contribute to my package, follow these steps:

1. Fork warehouse
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Submit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Pushed to branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## References
1.  Yao J, Sun Q, Wu H, et al. Decoding the molecular landscape: HER2 and PD-L1 in advanced gastric cancer. Front Immunol. 2025;16:1567308. <doi:10.3389/fimmu.2025.1567308>
2.  Qiu Q, Tan D, Chen Q, et al. Clinical implications of PD-L1 expression and pathway-related molecular subtypes in advanced Asian colorectal cancer patients. Am J Cancer Res. 2024;14(2):796-808. <doi:10.62347/FSSF9938>
3.  Ding W, Yang P, Zhao X, et al. Unraveling EGFR-TKI resistance in lung cancer with high PD-L1 or TMB in EGFR-sensitive mutations. Respir Res. 2024;25(1):40. <doi:10.1186/s12931-023-02656-3>
4.  Peng H, Ying J, Zang J, et al. Specific Mutations in APC, with Prognostic Implications in Metastatic Colorectal Cancer. Cancer Res Treat. 2023;55(4):1270-1280. <doi:10.4143/crt.2023.415>
5.  Jiang Y, Mai G, Zhao X, et al. Molecular characterization and prognostic implications of KRAS mutations in pancreatic cancer patients: insights from multi-cohort analysis. NPJ Precis Oncol. 2025;9(1):299. Published 2025 Aug 22.      <doi:10.1038/s41698-025-01087-1>
5.  Cancer Genome Atlas Network. Comprehensive molecular characterization of human colon and rectal cancer. Nature. 2012;487(7407):330-337. Published 2012 Jul 18. <doi:10.1038/nature11252>

## Licence

This project uses the MIT license-see LICENSE for details

## Contact information

For questions or suggestions, please contact:
- Mailbox：[ypm_yy@outlook.com]
- Github：[https://github.com/ypmjx123/Mypackage.git]
