#' Analyze Mutation and Clinical Data Correlation with MAF Format Data
#'
#' @description
#' This function analyzes mutation and clinical data correlation with MAF format data. It calculates the VAF (Variant Allele Frequency) for each gene and selects the top genes based on the specified criteria. Then, it calculates the correlation between the selected genes and a specified clinical variable. The function also generates a correlation plot to visualize the relationships between the genes and the clinical variable.
#' @param mutation_data 明确要求MAF格式数据框，且至少包含Tumor_Sample_Barcode和Hugo_Symbol两列
#' @param clin 临床数据框，至少包含样本ID
#' @param gene 基因列表，默认为NULL，当top不为NULL时，会自动选择前多少个基因，当top为NULL时，会使用gene列表中的基因
#' @param top 前多少个基因
#' @param cin_col 临床变量 clin数据框中的列名，需要第一列为样本ID
#' @param corrplot_method 相关系数图 corrplot方法 circle, square, ellipse, number, shade, color, pie
#'
#' @return  
#' A list containing the following elements: 
#' \item{cor_data}{A data frame containing the correlation between the selected genes and the clinical variable.}
#' @export  
#' @author Pengmin Yang
maf_cor<-function(mutation_data,clin,gene,top=10,cin_col,corrplot_method=c("circle", "square", "ellipse", "number",
                                                                           "shade", "color", "pie")){
  if (!is.null(top)){
    genes<-unique(mutation_data$Hugo_Symbol)
    gene_mut<-data.frame(Gene=genes)
    gene_mut$mut<-NA
    start_time <- Sys.time()
    for(Gene in genes){
      vaf<-length(unique(mutation_data$Tumor_Sample_Barcode[mutation_data$Hugo_Symbol==Gene]))/length(unique(mutation_data$Tumor_Sample_Barcode))
      gene_mut[gene_mut$Gene==Gene,"mut"]<- if (is.na(vaf)) 0 else as.numeric(vaf)
    }
    end_time <- Sys.time()
    total_time <- end_time - start_time
    print(total_time)
    gene_mut<-gene_mut[order(gene_mut$mut, decreasing = TRUE),]
    rownames(gene_mut)<-NULL
    gene<-gene_mut[1:top,"Gene"]
    gene_mut<-gene_mut[gene_mut$Gene%in%gene,]
  }else{
    gene_mut<-data.frame(Gene=gene)
    gene_mut$mut<-NA
    start_time <- Sys.time()
    for(Gene in gene){
      vaf<-length(unique(mutation_data$Tumor_Sample_Barcode[mutation_data$Hugo_Symbol==Gene]))/length(unique(mutation_data$Tumor_Sample_Barcode))
      gene_mut[gene_mut$Gene==Gene,"mut"]<- if (is.na(vaf)) 0 else as.numeric(vaf)
    }
    end_time <- Sys.time()
    total_time <- end_time - start_time
    print(total_time)
  }
  cor_data <- data.frame(Samples = unique(mutation_data$Tumor_Sample_Barcode), stringsAsFactors = FALSE)

  # 确保为基因动态增加列
  for (Gene in gene) {
    cor_data[[Gene]] <- NA  # 初始化为 NA
  }

  for (sample in unique(mutation_data$Tumor_Sample_Barcode)) {
    for (Gene in gene) {
      # 计算特定样本和基因的突变数量
      count <- length(mutation_data$Hugo_Symbol[mutation_data$Hugo_Symbol == Gene & mutation_data$Tumor_Sample_Barcode == sample])
      # 使用 which 寻找当前样本的行索引
      row_index <- which(cor_data$Samples == sample)
      cor_data[row_index, Gene] <- count
    }
  }
  if(is.null(clin)){
    cor_data<-cor_data
  }else if(!is.null(clin)){
    cor_data<-merge(cor_data,clin[,cin_col],by.x="Samples",by.y=colnames(clin[,cin_col])[[1]],all.x=T)
  }
  cor_data <- na.omit(cor_data)
  rownames(cor_data)<-cor_data[,1]
  cor_data<-cor_data[,-1]
  corr <- round(cor(cor_data), 1)
  res <- cor.mtest(cor_data, conf.level = .95)
  p <- res$p
  library(corrplot)
  library(ggplot2)
  #上三角图添加显著性水平星号：
  corrplot(corr, method = corrplot_method,
           type = c('upper'),
           col = NULL,
           outline = 'grey',
           order = c('AOE'),
           diag = TRUE,
           tl.cex = 1,
           tl.col = 'black',
           tl.pos = 'd',
           p.mat = p,
           sig.level = c(.001, .01, .05),
           insig = "label_sig", #显著性标注样式："pch", "p-value", "blank", "n", "label_sig"
           pch.cex = 1.2, #显著性标记大小
           pch.col = 'black' #显著性标记颜色
  )
  corrplot(corr, add = TRUE,
           method = c('number'),
           type = c('lower'),
           col = NULL,
           order = c('AOE'),
           diag = FALSE,
           number.cex = 0.9,
           tl.pos = 'n',
           cl.pos = 'n',
           p.mat = p,
           insig = "pch"
  )
  library(ggplot2)
  library(ggcor)
  set_scale()
  p<-quickcor(cor_data, cor.test = TRUE) +
    geom_square(data = get_data(type = "lower", show.diag = FALSE)) +
    geom_mark(data = get_data(type = "upper", show.diag = FALSE), size = 2.5) +
    geom_abline(slope = -1, intercept = length(gene)+length(cin_col))
  print(p)
  return(cor_data)
}
