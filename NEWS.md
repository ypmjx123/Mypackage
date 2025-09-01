# Mypackage 1.3.1

* 修改cox_forest函数

# Mypackage 1.3.0

* 增加cox_forest函数,Cox Proportional Hazards Univariate and Multivariate Forest Plot Generator

# Mypackage 1.2.11

* 增强immu_visual和path_mut_visual可视化阅读性

# Mypackage 1.2.10

* 增加使用说明文件

# Mypackage 1.2.9

* 修改HRDscore函数，根据样本ploidy数据进行HRD评分矫正

# Mypackage 1.2.8

* 增加HRDscore函数(引自scarHRD包),用于Determining genomic scar score (telomeric allelic imbalance, loss-off heterozigosity, large-scle transitions), signs of homologous recombination deficiency

# Mypackage 1.2.7

* 增加mut_network函数,用于Gene Mutation Co_Occurence/Mutually_Exclusive analysis and visualization

# Mypackage 1.2.6

* 修改limma.dif.visual函数,增加c("limma","DESeq2","edgeR")基因表达差异分析方法

# Mypackage 1.2.5

* 修改sig_Heatmap函数,annotation_tile替换原来的add_tile函数

# Mypackage 1.2.4

* 修改ggsurvplots函数,使中位值颜色和组别颜色一致，更易观察

# Mypackage 1.2.3

* 修改limma.dif.visual函数bug

# Mypackage 1.2.2

* 修改immu_visual函数bug

# Mypackage 1.2.1

* 增加maf_cor函数，用于Analyze Mutation and Clinical Data Correlation

# Mypackage 1.2.0

* 增加ScRNA函数，用于单细胞数据分析

# Mypackage 1.1.15

* 增加immu_visual和gene.mut_exp函数的数据清洗步骤

# Mypackage 1.1.13

* 增加mutToMAF函数输出表格数据的列名label，便于理解数据

# Mypackage 1.1.13

* 更新了mutToMAF函数的位点筛选规则

# Mypackage 1.1.12

* 修改了path_mut_visual函数

# Mypackage 1.1.11

* 修改了ggsurvplots函数

# Mypackage 1.1.10

* 修改了gene.mut_exp函数

# Mypackage 1.1.9

* 修改了immu_visual函数

# Mypackage 1.1.8

* 修改了ggsurvplots函数

# Mypackage 1.1.7

* 加入了ggsurvplots函数

# Mypackage 1.1.6

* 优化了limma.dif.visual函数,可以根据需要生成不同形式的富集图c("dotplot","barplot","heatplot","cnetplot","gseaplot")
* 优化了limma.dif.visual函数,让火山图更好看

# Mypackage 1.1.5

* 加入了GSEAplot2函数,优化了enrichplot::gseaplot2展现形式
* 优化了limma.dif.visual函数,可以根据需要生成不同形式的富集图c("dotplot","heatplot","cnetplot","gseaplot")

# Mypackage 1.1.4

* 优化了mutToMAF函数的基因突变筛选规则
* 删掉了tidyHeatmap::heatmap自动生成PDF
* 修改了tidyHeatmap::heatmap根据需要重新构建的生成方式，避免影响运行环境
* 修改了merge函数的by.x和by.y，更加具有普遍性

# Mypackage 1.1.3

* path_mut_visual,immu_visual和limma.dif.visual函数加入了tidyHeatmap::heatmap
* 修改immu_visual函数代码

# Mypackage 1.1.2

* 修改immu_visual函数代码，可以运行八种免疫评分
* 加入了sig_Heatmap函数，可以生成tidyHeatmap::heatmap热图

# Mypackage 1.1.1

* 加入了immu_visual函数
* 加入了gene.mut_exp函数
* 加入了exp_geneIDtoSYMBOL函数

# Mypackage 1.1.0

* 加入了limma.dif.visual函数

# Mypackage 0.1.1

* path_mut_visual通路突变差异可视化等基本功能
* Initial release
