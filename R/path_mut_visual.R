#' Pathway Mutation Visualization
#'
#' This function visualizes the mutation rate of tumor pathways based on SNV and gene data.
#'
#' @param result_data       A data frame containing mutation rate of tumor pathways
#' @param SNV               A data frame containing SNV data.即MAF格式的数据
#' @param gene              A character vector representing the gene of interest.例如BRAF基因或分组名称Group
#' @param Gene_group        A data frame that contains groups of genes.例如data.frame(SMPLE_ID,Group)
#' @param Type              A character vector indicating the type.Type can length from 2 to 5,eg Type=c("Wild", "Mut");(default = Type[[1]] is the control)
#' @param pathway_gene_data A data frame that contains pathway gene information.例如data.frame(Genes,Pathway)
#' @param tumor             A character vector indicating tumor types.例如CRC
#' @param heatmap          defalut heatmap=FALSE
#' @param heatmap_col      defalut heatmap_col=NULL
#' @param color             A character vector for defining colors.例如color=c("#757575", "#FF4040")
#' @param test              参考ggpubr::compare_means，根据需要选择
#' @param ns                没有显著性差异的结果展示与否，默认FALSE
#' @param p_0.05              筛选p值小于0.05，默认FALSE
#' @param p_0.01              筛选p值小于0.01，默认FALSE
#' @param p_0.001             筛选p值小于0.001，默认FALSE
#' @param p_0.0001            筛选p值小于0.0001，默认FALSE
#' @return A list of path mutation data,path mutation frequency in each type and ggplot object visualizing the mutation rates.
#' @export
##' @author Pengmin Yang
#' @examples
#' data("mutation_CRC")
#' data("gene_group_data")
#' data("pathway_data")
#' gene_of_interest<-colnames(gene_group_data)[[2]]
#' tumor_type="CRC TCGA"
#' color_vector=c("#757575", "#FF4040")
#' result <- path_mut_visual(result_data=NULL,
#'                            SNV = mutation_CRC,
#'                            gene = gene_of_interest,
#'                            Gene_group = gene_group_data,
#'                            Type = c("Wild", "Mut"),
#'                            pathway_gene_data = pathway_data,
#'                            tumor = tumor_type,
#'                            heatmap=TRUE,
#'                            heatmap_col=NULL,
#'                            color = color_vector,
#'                            test = "wilcox.test",ns=FALSE,
#'                            p_0.05=FALSE,p_0.01=FALSE,p_0.001=FALSE,p_0.0001=FALSE)
#' print(result)
path_mut_visual<-function(result_data=NULL,SNV,gene,Gene_group,Type,pathway_gene_data,
                          tumor,heatmap=FALSE,heatmap_col=NULL,color,test="wilcox.test",ns=FALSE,
                          p_0.05=FALSE,p_0.01=FALSE,p_0.001=FALSE,p_0.0001=FALSE){
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(ggpubr)
  if(is.null(result_data)){
    sample_gene_data <- SNV[c("Tumor_Sample_Barcode", "Hugo_Symbol")]
    unique_samples <- unique(sample_gene_data$Tumor_Sample_Barcode)
    pathway_gene_data<-pathway_gene_data[pathway_gene_data$Genes!=gene,]
    pathways<-unique(pathway_gene_data$Pathway)
    result_data <- data.frame(Tumor_Sample_Barcode = unique_samples)
    for (pathway in pathways) {
      mutation_rates <- numeric(length(unique_samples))
      for (i in 1:length(unique_samples)) {
        sample <- unique_samples[i]
        sample_genes <- unique(sample_gene_data$Hugo_Symbol[sample_gene_data$Tumor_Sample_Barcode == sample])
        pathway_counts <- sum(sample_genes %in% pathway_gene_data$Genes[pathway_gene_data$Pathway== pathway])
        pathway_mutation_rate <- (pathway_counts / length(pathway_gene_data$Genes[pathway_gene_data$Pathway== pathway]))*100
        mutation_rates[i] <- pathway_mutation_rate
      }
      result_data[[pathway]] <- mutation_rates
    }
    result_data<-merge(Gene_group,result_data,by.x=colnames(Gene_group)[[1]],by.y=colnames(result_data)[[1]],all.y=T)
    result_data[,2][is.na(result_data[,2])]<-Type[[1]]
    result_data[[colnames(result_data)[2]]]<-as.factor(result_data[[colnames(result_data)[2]]])
  }else{
    result_data<-result_data
  }
  if(heatmap){
    heatmap<-sig_Heatmap(input = result_data,
                         features = colnames(result_data[,3:length(colnames(result_data))]),
                         ID =colnames(result_data)[[1]],show_plot=F,#heatmap_col=c("#0505FA", "#FFFFFF", "#FA050D"),
                         cols_group= color ,heatmap_col=NULL,
                         group = gene, scale = T,column_title=NULL,column_title_size=8,
                         row_title=NULL,row_title_size=8,name="Pathway Frequency")
  }else{
    heatmap<-NULL
  }


  columns_to_delete <- c(1:2)
  cols_to_convert <- setdiff(names(result_data), names(result_data)[columns_to_delete])
  result_data[cols_to_convert] <- lapply(result_data[cols_to_convert], as.numeric)
  result <- data.frame(
    Types = character(),
    Pathway = character(),
    Mutation_rate = numeric(),
    SE = numeric(),
    stringsAsFactors = FALSE
  )
  Pathways<- names(result_data)[-columns_to_delete]
  unique_types <- unique(result_data[,2])
  for (type in unique_types) {
    for (pathway in Pathways){
      mutation_values <- result_data[result_data[,2] == type, pathway]
      mutation_values  <- as.numeric(mutation_values)
      se <- sd(mutation_values, na.rm = TRUE)
      #计算突变率
      mutation_rate<-round(100 * sum(mutation_values != 0) / sum((result_data[,2] == type)), 2)
      mutation_rate<- as.numeric(mutation_rate)
      # 添加到结果数据框中
      result <- rbind(result, data.frame(
        Types = type,
        Pathway = pathway,
        Mutation_rate = mutation_rate,
        SE = se))
    }
  }
  result$Types <- factor(result$Types, levels = Type)
  Pathways1<-unique(pathway_gene_data$Pathway)
  Pathways2<-unique(result$Pathway)
  data_long1 <- result_data %>%
    pivot_longer(cols = -columns_to_delete,
                 names_to = "pathway", values_to = "mut")
  data_long1[[gene]]<- factor(data_long1[[gene]], levels = Type)#
  data_long1_path<-unique(data_long1$pathway)
  replacement_map1 <- data.frame(
    original = Pathways2,
    replacement = Pathways1
  )
  replacement_map2 <- data.frame(
    original = data_long1_path,
    replacement = Pathways1
  )
  result$Pathway <- replacement_map1$replacement[match(result$Pathway, replacement_map1$original)]
  data_long1$pathway<- replacement_map2$replacement[match(data_long1$pathway, replacement_map2$original)]
  p<-ggplot(data = result, mapping = aes(x = Pathway, y = Mutation_rate, fill = Types)) +
    geom_bar(position = "dodge", width = 0.8, stat = "identity") +
    theme(legend.position = "bottom",#c(0.9,0.9), #legend.title = element_blank(),
          #legend.title = element_text(),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5,color = "black"),
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1,color = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(linewidth = 1),
          panel.background = element_rect(fill = "transparent"), # 设定底色为透明
          plot.title = element_text(hjust = 0.5, margin = margin(b = 20))) +
    labs(x = "Tumor signaling pathways", y = "Mutation rate(%)",fill =paste(gene,"Group"),
         title = if(is.null(tumor))"Frequency of tumor signaling  pathway alterations" else paste("Frequency of tumor signaling  pathway alterations in",tumor)) +
    coord_flip() +
    geom_text(aes(label = paste0(Mutation_rate, "%")),
              position = position_dodge2(width = 0.8,preserve = "total",reverse = F),
              size = 2.5, vjust = 0.5, hjust = 1)+
    scale_fill_manual(values =color)#c("#757575", "#FF4040")
  if(ns){
    filtered_annotation <- function(annotations, x_min, x_max){
      valid_indices <- which(annotations %in% c("NS.", "ns","*","**","***","****"))
      list(
        annotations = annotations[valid_indices],
        x_min = x_min[valid_indices],
        x_max = x_max[valid_indices]
      )
    }
  }else if (p_0.05){
    filtered_annotation <- function(annotations, x_min, x_max){
      valid_indices <- which(!annotations %in% c("NS.", "ns"))
      list(
        annotations = annotations[valid_indices],
        x_min = x_min[valid_indices],
        x_max = x_max[valid_indices]
      )
    }
  }else if (p_0.01){
    filtered_annotation <- function(annotations, x_min, x_max){
      valid_indices <- which(!annotations %in% c("NS.", "ns","*"))
      list(
        annotations = annotations[valid_indices],
        x_min = x_min[valid_indices],
        x_max = x_max[valid_indices]
      )
    }
  }else if (p_0.001){
    filtered_annotation <- function(annotations, x_min, x_max){
      valid_indices <- which(!annotations %in% c("NS.", "ns","*","**"))
      list(
        annotations = annotations[valid_indices],
        x_min = x_min[valid_indices],
        x_max = x_max[valid_indices]
      )
    }
  }else if (p_0.0001){
    filtered_annotation <- function(annotations, x_min, x_max){
      valid_indices <- which(!annotations %in% c("NS.", "ns","*","**","***"))
      list(
        annotations = annotations[valid_indices],
        x_min = x_min[valid_indices],
        x_max = x_max[valid_indices]
      )
    }
  }
  get_annotation <- function(data, types) {
    annotations <- c()
    data_subset_filtered <- data[data[[gene]] %in% types, ]#
    test_result <- ggpubr::compare_means(as.formula(paste0(colnames(data)[[4]] , " ~ ", gene)), data_subset_filtered, #
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
  pathways <- unique(data_long1$pathway)
  A<-as.numeric(length(pathways))-1
  pathwithsig<-c()
  i=0
  p_0.05<-c()
  color_vector<-rep("black", length(pathways))
  while(i<=A ){
    result1<-data_long1[data_long1$pathway==pathways[i+1],]
    if(length(Type)==4){
      annotation1 <- get_annotation(result1, c(Type[[1]],Type[[2]]))
      annotation2 <- get_annotation(result1, c(Type[[1]],Type[[3]]))
      annotation3 <- get_annotation(result1, c(Type[[1]],Type[[4]]))
      annotation4 <- get_annotation(result1, c(Type[[2]],Type[[3]]))
      annotation5 <- get_annotation(result1, c(Type[[2]],Type[[4]]))
      annotation6 <- get_annotation(result1, c(Type[[3]],Type[[4]]))
      all_annotation <- c(annotation1,annotation2,annotation3,annotation4,annotation5,annotation6)
      all_x_min <- c(0.70,0.70,0.70,0.90,0.90,1.10)+i
      all_x_max <- c(0.90,1.10,1.30,1.10,1.30,1.30)+i
      filtereds <- filtered_annotation(all_annotation, all_x_min,all_x_max)
      y_positions <- c()
      y_position1<-result$Mutation_rate[result$Types==Type[[1]]&result$Pathway==pathways[i+1]]
      y_position2<-result$Mutation_rate[result$Types==Type[[2]]&result$Pathway==pathways[i+1]]
      y_position3<-result$Mutation_rate[result$Types==Type[[3]]&result$Pathway==pathways[i+1]]
      y_position4<-result$Mutation_rate[result$Types==Type[[4]]&result$Pathway==pathways[i+1]]
      y_position<-pmax(y_position1,y_position2,y_position3,y_position4)
      y_positions <- c(y_positions, y_position)
      all_y_position <- c(y_positions+2)
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
        cat("\033[31m", paste(i+1,":","Oh yeah!",pathways[i+1],"有组间显著性差异"), "\033[0m\n")
        pathwithsig<-c(pathwithsig,pathways[i+1])
        p_0.05<-c(p_0.05,i+1)
      }else{cat("\033[32m", paste(i+1,":","Oopps!",pathways[i+1],"没有组间显著性差异"), "\033[0m\n")}
      i=i+1
    } else if(length(Type)==3){
      annotation1 <- get_annotation(result1, c(Type[[1]],Type[[2]]))
      annotation2 <- get_annotation(result1, c(Type[[1]],Type[[3]]))
      annotation3 <- get_annotation(result1, c(Type[[2]],Type[[3]]))
      all_annotation <- c(annotation1,annotation2,annotation3)
      all_x_min <- c(0.80,0.80,1.00)+i
      all_x_max <- c(1.00,1.20,1.20)+i
      filtereds <- filtered_annotation(all_annotation, all_x_min,all_x_max)
      y_positions <- c()
      y_position1<-result$Mutation_rate[result$Types==Type[[1]]&result$Pathway==pathways[i+1]]
      y_position2<-result$Mutation_rate[result$Types==Type[[2]]&result$Pathway==pathways[i+1]]
      y_position3<-result$Mutation_rate[result$Types==Type[[3]]&result$Pathway==pathways[i+1]]
      y_position<-pmax(y_position1,y_position2,y_position3)
      y_positions <- c(y_positions, y_position)
      all_y_position <- c(y_positions+2)
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
        cat("\033[31m", paste(i+1,":","Oh yeah!",pathways[i+1],"有组间显著性差异"), "\033[0m\n")
        pathwithsig<-c(pathwithsig,pathways[i+1])
        p_0.05<-c(p_0.05,i+1)
      }else{cat("\033[32m", paste(i+1,":","Oopps!",pathways[i+1],"没有组间显著性差异"), "\033[0m\n")}
      i=i+1
    }else if(length(Type)==2){
      annotation1 <- get_annotation(result1, c(Type[[1]],Type[[2]]))
      all_annotation <- c(annotation1)
      all_x_min <- c(0.80)+i
      all_x_max <- c(1.20)+i
      filtereds <- filtered_annotation(all_annotation, all_x_min,all_x_max)
      y_positions <- c()
      y_position1<-result$Mutation_rate[result$Types==Type[[1]]&result$Pathway==pathways[i+1]]
      y_position2<-result$Mutation_rate[result$Types==Type[[2]]&result$Pathway==pathways[i+1]]
      y_position<-pmax(y_position1,y_position2)
      y_positions <- c(y_positions, y_position)
      all_y_position <- c(y_positions+2)
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
          cat("\033[31m", paste(i+1,":","Oh yeah!",pathways[i+1],"有组间显著性差异(p<0.05)"), "\033[0m\n")
        }else if(filtereds[["annotations"]]=="**"){
          cat("\033[31m", paste(i+1,":","Oh yeah!",pathways[i+1],"有组间显著性差异(p<0.01)"), "\033[0m\n")
        }else if(filtereds[["annotations"]]=="***"){
          cat("\033[31m", paste(i+1,":","Oh yeah!",pathways[i+1],"有组间显著性差异(p<0.001)"), "\033[0m\n")
        }else if(filtereds[["annotations"]]=="****"){
          cat("\033[31m", paste(i+1,":","Oh yeah!",pathways[i+1],"有组间显著性差异(p<0.0001)"), "\033[0m\n")
        }
        p_0.05<-c(p_0.05,i+1)
        pathwithsig<-c(pathwithsig,pathways[i+1])
      }else{cat("\033[32m", paste(i+1,":","Oopps!",pathways[i+1],"没有组间显著性差异"), "\033[0m\n")}
      i=i+1
    }else if(length(Type)==5){
      annotation01 <- get_annotation(result1, c(Type[[1]],Type[[2]]))
      annotation02 <- get_annotation(result1, c(Type[[1]],Type[[3]]))
      annotation03 <- get_annotation(result1, c(Type[[1]],Type[[4]]))
      annotation04 <- get_annotation(result1, c(Type[[1]],Type[[5]]))
      annotation05 <- get_annotation(result1, c(Type[[2]],Type[[3]]))
      annotation06 <- get_annotation(result1, c(Type[[2]],Type[[4]]))
      annotation07 <- get_annotation(result1, c(Type[[2]],Type[[5]]))
      annotation08 <- get_annotation(result1, c(Type[[3]],Type[[4]]))
      annotation09 <- get_annotation(result1, c(Type[[3]],Type[[5]]))
      annotation10 <- get_annotation(result1, c(Type[[4]],Type[[5]]))
      all_annotation <- c(annotation01,annotation02,annotation03,annotation04,annotation05,annotation06,
                          annotation07,annotation08,annotation09,annotation10)
      all_x_min <- c(0.68,0.68,0.68,0.68,0.84,0.84,0.84,1.00,1.00,1.16)+i
      all_x_max <- c(0.84,1.00,1.16,1.32,1.00,1.16,1.32,1.16,1.32,1.32)+i
      filtereds <- filtered_annotation(all_annotation, all_x_min,all_x_max)
      y_positions <- c()
      y_position1<-result$Mutation_rate[result$Types==Type[[1]]&result$Pathway==pathways[i+1]]
      y_position2<-result$Mutation_rate[result$Types==Type[[2]]&result$Pathway==pathways[i+1]]
      y_position3<-result$Mutation_rate[result$Types==Type[[3]]&result$Pathway==pathways[i+1]]
      y_position4<-result$Mutation_rate[result$Types==Type[[4]]&result$Pathway==pathways[i+1]]
      y_position5<-result$Mutation_rate[result$Types==Type[[5]]&result$Pathway==pathways[i+1]]
      y_position<-pmax(y_position1,y_position2,y_position3,y_position4,y_position5)
      y_positions <- c(y_positions, y_position)
      all_y_position <- c(y_positions+2)
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
        cat("\033[31m", paste(i+1,":","Oh yeah!",pathways[i+1],"有组间显著性差异"), "\033[0m\n")
        pathwithsig<-c(pathwithsig,pathways[i+1])
        p_0.05<-c(p_0.05,i+1)
      }else{cat("\033[32m", paste(i+1,":","Oopps!",pathways[i+1],"没有组间显著性差异"), "\033[0m\n")}
      i=i+1
    }
  }
  color_vector[p_0.05] <- "red"
  p<-p+theme(axis.text.y = element_text(angle = 0,color = color_vector))
  results <- list(
    path_mut = result_data,
    heatmap =heatmap,
    path_mut_Frequency = result,
    path_mut_plot = p,
    path_sig=pathwithsig
  )
  return(results)
}


