#' Gene Mutation and Expression Data Visualization
#'
#' This function visualizes the  Gene Expression Difference with gene mutated status.
#'
#' @param mutation_data    A data frame containing SNV data.即MAF格式的数据
#' @param exp              is the data of gene expression for samples
#' @param colnum           exp 列名为 sample开始的列数
#' @param exp_cleaning     是否对数据进行清洗，默认FALSE
#' @param gene             gene name,gene=c("PREX2","TP53","APC")
#' @param top              突变率前x的基因, top=10
#' @param visual           可视化否，默认TURE
#' @param test_type        ggstatsplot::ggbetweenstats,test_type= "parametric"
#' @param title            ggstatsplot::ggbetweenstats,title = NULL
#' @param test             参考ggpubr::compare_means，根据需要选择
#' @param gene_vaf         gene_vaf defalut is FALSE
#' @param color_0.05/color_0.01/color_0.001/color_0.0001 axis.text.x颜色，根据不同显著性显示
#' @param only_red         axis.text.x颜色，根据显著性显示
#' @param bar              defalut is FALSE
#'
#' @return A list of Gene Expression Difference with gene mutated status and visualizing
#' @export
##' @author Pengmin Yang
#' @examples
#' data("exp_raw")
#' data("mutation_CRC")
#' exp_CRC<-exp_geneIDtoSYMBOL(exp=exp_raw,genecoltype="ENTREZID")
#' exp_CRC<-exp_CRC$data
#' result<-gene.mut_exp(mutation_data=mutation_CRC,exp=exp_CRC[,-1],colnum=2,gene=NULL,top=10,visual=TRUE,
#'             test_type= "parametric",title= "CRC TCGA",test="wilcox.test",only_red=TRUE,gene_vaf=F,
#'              color_0.05= "#FF34B3",color_0.01="#9400D3",color_0.001="#CD3700",color_0.0001="red",bar=T)
gene.mut_exp <- function(mutation_data, exp, colnum = 3, exp_cleaning=FALSE,gene,top=10, visual = TRUE,
                         test_type = "parametric", title = "CRC TCGA", test = "wilcox.test",only_red=TRUE,gene_vaf=F,
                         color_0.05= "#FF34B3",color_0.01="#9400D3",color_0.001="#CD3700",color_0.0001="red",bar=F) {
  library(ggplot2)
  if(exp_cleaning){
    n <- ncol(exp)
    pseudo_count <- min(c(exp[,colnum:n][exp[,colnum:n] > 0])) / 2
    exp[,colnum:n] <- log10(exp[,colnum:n] + pseudo_count)
  }
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
  results<-list()
  for(i in 1:length(gene)){
    samples <- colnames(exp[, c(colnum:length(colnames(exp)))])
    samples <- data.frame(SAMPLE_ID = samples)
    samples$Gene <- "Wild"
    gene_SNV <- mutation_data[mutation_data$Hugo_Symbol == gene[[i]], c("Tumor_Sample_Barcode", "Hugo_Symbol")]
    samples$Gene <- ifelse(samples$SAMPLE_ID %in% gene_SNV$Tumor_Sample_Barcode, "Mut", samples$Gene)
    genecol <- c("SYMBOL", "Hugo_Symbol")
    if (colnames(exp)[[1]] %in% genecol) {
      gene_exp <- exp[exp[[colnames(exp)[[1]]]] == gene[[i]], ]
      rownames(gene_exp) <- NULL
      gene_exp <- gene_exp[, c(colnum:length(colnames(exp)))]
      gene_exp <- t(gene_exp)
      gene_exp <-as.data.frame(gene_exp)
      gene_exp$SAMPLE_ID<-rownames(gene_exp)
      gene_exp <- gene_exp[order(gene_exp$SAMPLE_ID), ]
      gene_exp <- gene_exp[,c(2,1)]
      samples$gene_exp <- gene_exp$V1
      if (visual) {
        if(length(gene)==1){
          p <- ggstatsplot::ggbetweenstats(
            data = samples,
            x = Gene,
            y = gene_exp,
            ylab = paste(gene[[i]], "Expression"),
            xlab = paste(gene[[i]], "Mutation Status"),
            type = test_type,
            title = title
          )
        }
        colnames(samples) <- c("SAMPLE_ID", gene[[i]], paste0(gene[[i]], "_exp"))
      } else {
        colnames(samples) <- c("SAMPLE_ID", gene[[i]], paste0(gene[[i]], "_exp"))
      }
    } else {
      stop("\nERROR: \"ExprData\" loss a gene colum. Gene column name must be SYMBOL or Hugo_Symbol")
      results <- NULL
      p<-NULL
    }
    results<-c(results,samples)
    results<-as.data.frame(results)
  }
  if (visual) {
    if(length(gene)>1){
      library(dplyr)
      library(tidyr)
      results<- results %>%
        select(SAMPLE_ID, everything()) %>%
        select(-matches("^SAMPLE_ID\\..*"))
      colnames(results) <- gsub("\\.", "-", colnames(results))
      data_long <- pivot_longer(results,
                                cols = gene,
                                names_to = "Gene",
                                values_to = "Mutation_Status")
      data_long_exp <- pivot_longer(results,
                                    cols = paste0(gene, "_exp"),
                                    names_to = "Gene_Exp",
                                    values_to = "Expression_Value")
      final_data <- cbind(data_long[, c("SAMPLE_ID", "Gene", "Mutation_Status")],
                          data_long_exp[, "Expression_Value"])
      final_data$Gene <- gsub("_exp", "", final_data$Gene)
      gene_mut$Gene1<-paste(gene_mut$Gene, "(", paste0(round(gene_mut$mut, 4) * 100, "%", ")"), sep = "")
      final_data<-merge(final_data,gene_mut[,c("Gene","Gene1","mut")],by.x="Gene",by.y="Gene",all.x=T)
      if(gene_vaf){
        p<-ggplot(final_data, aes(x = Gene1,
                                  y = Expression_Value,
                                  fill = Mutation_Status))
      }else{
        p<-ggplot(final_data, aes(x = Gene,
                                  y = Expression_Value,
                                  fill = Mutation_Status))
      }
      p<-p +
        geom_boxplot(stat = "boxplot",outliers = F,
                     position = "dodge2")+
        labs(x = " ", y = "Expression", fill = "Mutation status",
             title = title)+theme(legend.position ="top",
                                  axis.text.y = element_text(angle = 0,color = "black"),
                                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                  axis.line = element_line(linewidth = 1),
                                  panel.background = element_rect(fill = "transparent"))
      filtered_annotation <- function(annotations, x_min, x_max){
        valid_indices <- which(!annotations %in% c("NS.", "ns"))
        list(
          annotations = annotations[valid_indices],
          x_min = x_min[valid_indices],
          x_max = x_max[valid_indices]
        )
      }
      get_annotation <- function(data, types) {
        annotations <- c()
        data_subset_filtered <- data[data$Mutation_Status %in% types, ]#
        test_result <- ggpubr::compare_means(Expression_Value ~ Mutation_Status, data_subset_filtered, #
                                             method = test)# "kruskal.test"
        test_result <- as.data.frame(test_result)
        if (is.data.frame(test_result) && nrow(test_result) == 0) {
          annotations <- c(annotations, "NS.")
        } else if (is.data.frame(test_result) && nrow(test_result) != 0
                   && test_result$p.signif == "?") {
          annotations <- c(annotations, "NS.")
        } else if (is.data.frame(test_result) && nrow(test_result) != 0
                   && test_result$p.signif != "?") {
          annotations <- c(annotations, test_result$p.signif)
        }
        return(annotations)
      }
      A<-as.numeric(length(gene))-1
      i=0
      Type=unique(final_data$Mutation_Status)
      gene<-sort(gene)
      Gene<-sort(unique(final_data$Gene1))
      p_0.05<-c()
      color_vector<-rep("black", length(gene))
      while(i<=A ){
        result1<-final_data[final_data$Gene1==Gene[i+1],]
        if(length(unique(result1$Mutation_Status))>=2){
          annotation1 <- get_annotation(result1, c(Type[[1]],Type[[2]]))
          all_annotation <- c(annotation1)
          all_x_min <- c(0.80)+i
          all_x_max <- c(1.20)+i
          y_positions <- c()
          myQuart1 =quantile(as.numeric(result1$Expression_Value[result1$Mutation_Status == Type[[1]]]), na.rm = TRUE)
          IQR1 =  myQuart1[4] - myQuart1[2]
          y_position1<-myQuart1[4] + 1.5 * IQR1
          myQuart2 =quantile(as.numeric(result1$Expression_Value[result1$Mutation_Status == Type[[2]]]), na.rm = TRUE)
          IQR2 =  myQuart2[4] - myQuart2[2]
          y_position2<-myQuart2[4] + 1.5 * IQR2
          y_position<-pmax(y_position1,y_position2)
          y_positions <- c(y_positions, y_position)
          all_y_position <- c(y_positions)
        }else{all_annotation <- c("NS.")
        all_x_min <- i
        all_x_max <- i
        y_positions <- c()
        myQuart1 =quantile(as.numeric(result1$Expression_Value[result1$Mutation_Status == "Wild"]), na.rm = TRUE)
        IQR1 =  myQuart1[4] - myQuart1[2]
        y_position1<-myQuart1[4] + 1.5 * IQR1
        y_positions <- c(y_positions, y_position1)
        all_y_position <- c(y_positions)}
        filtereds <- filtered_annotation(all_annotation, all_x_min,all_x_max)
        if(length(filtereds[["annotations"]])>0){
          all_y_position<-all_y_position[1:length(filtereds[["annotations"]])]
          p<-p + ggpubr::geom_signif(
            size = 0.25,
            tip_length =c(0.01,0.01),
            vjust = 1.0,
            y_position = all_y_position,
            xmin = filtereds[["x_min"]],
            xmax = filtereds[["x_max"]],
            annotation = filtereds[["annotations"]]
          )
          if(filtereds[["annotations"]]=="*"){
            color_vector[i+1] <- color_0.05
            cat("\033[31m", paste(i+1,":","Oh yeah!",gene[i+1],"有组间显著性差异(p<0.05)"), "\033[0m\n")
          }else if(filtereds[["annotations"]]=="**"){
            color_vector[i+1] <- color_0.01
            cat("\033[31m", paste(i+1,":","Oh yeah!",gene[i+1],"有组间显著性差异(p<0.01)"), "\033[0m\n")
          }else if(filtereds[["annotations"]]=="***"){
            color_vector[i+1] <- color_0.001
            cat("\033[31m", paste(i+1,":","Oh yeah!",gene[i+1],"有组间显著性差异(p<0.001)"), "\033[0m\n")
          }else if(filtereds[["annotations"]]=="****"){
            color_vector[i+1] <- color_0.0001
            cat("\033[31m", paste(i+1,":","Oh yeah!",gene[i+1],"有组间显著性差异(p<0.0001)"), "\033[0m\n")
          }
          p_0.05<-c(p_0.05,i+1)
        }else{cat("\033[32m", paste(i+1,":","Oopps!",gene[i+1],"没有组间显著性差异"), "\033[0m\n")}
        i=i+1
      }
    }
    if(only_red){
      color_vector[p_0.05] <- "red"
    }
    if(length(gene)<=5){
      p<-p+theme(axis.text.x = element_text(angle = 0,color = color_vector))
    }else{
      p<-p+theme(axis.text.x = element_text(angle = 45,hjust = 1,color = color_vector))
    }
    if(bar){
      gene_mut <- gene_mut %>% arrange(Gene)
      gene_mut$VAF<-paste0(round(gene_mut$mut, 3) * 100, "%")
      p<-p+geom_bar(data = gene_mut,aes(x = Gene,y = mut*-10), stat = "identity", fill = "gray", alpha = 0.5)
      for( i in 1:nrow(gene_mut)){
        p<-p+ggpubr::geom_signif(xmin = i,xmax=i,y_position = gene_mut[i,"mut"]*-10,
                                 annotation = gene_mut[i,"VAF"],size = 0.25,
                                 tip_length =c(0.00,0.00),vjust = 1.2)
      }
    }
  }else{
    results-results
    gene_mut<-gene_mut
    p_0.05<-NULL
    p<-NULL
  }
  result_all<-list(
    geneMutexp=results,
    genemut= gene_mut,
    genewithsig=gene[p_0.05],
    plot=p
  )
  return(result_all)
}


