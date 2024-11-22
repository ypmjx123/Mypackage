#' Immune score Visualization
#'
#' This function visualizes the mutation rate of immune score based on provided immune score or gene expression data.
#'
#' @param im           is the data of immune score for samples;if exp=NULL,im must be provided
#' @param exp          is the data of gene expression for samples;please ensure gene coloum with unique values;if im=NULL,exp must be provided
#' @param method       IOBR::deconvo_tme 参数或者根据其它方法得到的样本免疫数据
#' @param method       im=NULL method =c('mcpcounter', 'epic',  'cibersort', 'estimate')
#' @param sample_group is a dataframe for sample grouping
#' @param Type         can length from 2 to 5,eg Type=c("Wild", "Mut");(default = Type[[1]] is the control)
#' @param color        for plot of Type, eg color=c("#757575", "#FF4040")
#' @param geom_text    text annotation for ggplot
#' @return A list of immune score for data samples,method and ggplot object visualizing the immune score.
#' @export
#' @examples
#' #im=NULL
#' data("exp_CRC")
#' data("Gene_group_CRC1")
#' result<-immu_visual(im=NULL,exp=exp_CRC[,-2],
#'                     method = 'epic',
#'                     sample_group=Gene_group_CRC1,
#'                     tumor="CRC TCGA",
#'                     Type=c("Wild", "Mut"),
#'                     color=c("#757575", "#FF4040"),
#'                     geom_text=T)
#' print(result)
#'
#' #exp=NULL
#' data("im")
#' data("sample_group")
#' result<-immu_visual(im=im,exp=NULL,
#'                     method = ' ',
#'                     sample_group=sample_group,
#'                     tumor="GC TCGA",
#'                     Type=c("Wild", "Mut"),
#'                     color=c("#757575", "#FF4040"),
#'                     geom_text=F)
#' print(result)

immu_visual<-function(im,exp,method,sample_group,tumor,Type,color,geom_text=T){
  # im is the data of immune score for samples
  # exp is the data of gene expression for samples;please ensure gene coloum with unique values
  # method deconvo_tme 参数或者根据其它方法得到的样本免疫数据
  # if im=NULL method =c('mcpcounter', 'epic',  'cibersort', 'estimate')
  # sample_group is a dataframe for sample grouping
  # Type can length from 2 to 5,eg Type=c("Wild", "Mut");(default = Type[[1]] is the control)
  # color for plot of Type, eg color=c("#757575", "#FF4040")
  # if im=NULL,exp must be provided;if exp=NULL,im must be provided
  start_time <- Sys.time()
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(ggpubr)
  library(tools)
  if(is.null(im)&&!is.null(exp)){
    library(IOBR)
    #method=c('mcpcounter', 'epic',  'cibersort', 'estimate')
    exp<-exp[!is.na(exp[,2]),]
    if (any(duplicated(exp[[1]]))) {
      exp <- aggregate(as.formula(paste0(".", " ~ ", colnames(exp)[[1]])), data = exp, FUN = mean)
    }
    rownames(exp)<-exp[[1]]
    exp<-exp[,-1]
    im <- deconvo_tme(eset = exp,
                      method = method,
                      arrays = F,#plot = T,
                      perm = 1000
    )
    im<-merge(sample_group,im,by.x=colnames(sample_group)[[1]],by.y=colnames(im)[[1]],all.y=T)
    im[,2][is.na(im[, 2])] <- Type[[1]]
    if(method=='cibersort'){
      im<-im[,-c(length(colnames(im))-2,length(colnames(im))-1,length(colnames(im)))]
      data_long <- im%>%
        pivot_longer(cols = c(3:length(colnames(im))), names_to = method, values_to = paste0(method,"_score"))
      data_long[[method]] <- gsub("_CIBERSORT$", "", data_long[[method]])
      data_long[[method]] <- gsub("_", " ", data_long[[method]])
    }else if(method=='mcpcounter'){
      data_long <- im%>%
        pivot_longer(cols = c(3:length(colnames(im))), names_to = method, values_to = paste0(method,"_score"))
      data_long[[method]] <- gsub("_MCPcounter$", "", data_long[[method]])
      data_long[[method]] <- gsub("_", " ", data_long[[method]])
    }else if(method=='epic'){
      data_long <- im%>%
        pivot_longer(cols = c(3:length(colnames(im))), names_to = method, values_to = paste0(method,"_score"))
      data_long[[method]] <- gsub("_EPIC$", "", data_long[[method]])
      data_long[[method]] <- gsub("_", " ", data_long[[method]])
    }else{
      data_long <- im%>%
        pivot_longer(cols = c(3:length(colnames(im))), names_to = method, values_to = paste0(method,"_score"))
      data_long[[method]] <- gsub(paste0("_",method,"$"), "", data_long[[method]])
      data_long[[method]] <- gsub("_", " ", data_long[[method]])
    }
  }else{
    im<-merge(sample_group,im,by.x=colnames(sample_group)[[1]],by.y=colnames(im)[[1]],all.y=T)
    im[,2][is.na(im[, 2])] <- Type[[1]]
    data_long <- im%>%
      pivot_longer(cols = c(3:length(colnames(im))), names_to = method, values_to = paste0(method,"_score"))
    data_long[[method]] <- gsub(paste0("_",method,"$"), "", data_long[[method]])
    data_long[[method]] <- gsub("_", " ", data_long[[method]])
  }
  data_long[[colnames(sample_group)[[2]]]] <- factor(data_long[[colnames(sample_group)[[2]]]], levels = Type)
  data_long[[method]] <- factor(data_long[[method]] , levels = unique(data_long[[method]]))
  colnames(data_long)[colnames(data_long) == method] <-"method"
  colnames(data_long)[colnames(data_long) == paste0(method,"_score")] <-"method_score"
  colnames(data_long)[colnames(data_long) == colnames(data_long)[[2]]] <-"Group"
  max<-max(data_long$method_score)
  p1<-ggplot(data_long, aes(x=method, y=method_score, fill=Group)) +
    geom_boxplot(stat = "boxplot",outliers = F,
                 position = "dodge2")+
    theme(legend.position ="top",
          #axis.text.x = element_text(angle = 45,hjust = 1,color = "black"),
          axis.text.y = element_text(angle = 0,color = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(size = 1),
          panel.background = element_rect(fill = "transparent"))+
    scale_fill_manual(values=color)
  if(length(colnames(im[,-c(1,2)]))<=5){
    p1<-p1+theme(axis.text.x = element_text(angle = 0,color = "black"))
  }else{
    p1<-p1+theme(axis.text.x = element_text(angle = 45,hjust = 1,color = "black"))
  }
  if(method=='mcpcounter'){
    p1<-p1+labs(x = "", y ="MCPcounter Score",fill = paste(colnames(sample_group)[[2]],"Group"))
  }else{
    p1<-p1+labs(x = "", y = toTitleCase(paste(method,"score")),fill = paste(colnames(sample_group)[[2]],"Group"))
  }
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
    data_subset_filtered <- data[data$Group %in% types, ]#
    test_result <- ggpubr::compare_means(method_score ~ Group, data_subset_filtered, #
                                         method = "wilcox.test")# "kruskal.test"
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
  ims <- unique(data_long$method)
  A<-as.numeric(length(ims))-1
  i=0
  if(length(Type)==2){
    while(i<=A ){
      result1<-data_long[data_long$method==ims[i+1],]
      annotation1 <- get_annotation(result1, c(Type[[1]],Type[[2]]))
      all_annotation <- c(annotation1)
      all_x_min <- c(0.80)+i
      all_x_max <- c(1.20)+i
      filtereds <- filtered_annotation(all_annotation, all_x_min,all_x_max)
      y_positions <- c()
      y_positions <- c()
      myQuart1 =quantile(as.numeric(result1$method_score[result1$Group == Type[[1]]]), na.rm = TRUE)
      IQR1 =  myQuart1[4] - myQuart1[2]
      y_position1<-myQuart1[4] + 1.5 * IQR1
      myQuart2 =quantile(as.numeric(result1$method_score[result1$Group == Type[[2]]]), na.rm = TRUE)
      IQR2 =  myQuart2[4] - myQuart2[2]
      y_position2<-myQuart2[4] + 1.5 * IQR2
      y_position<-pmax(y_position1,y_position2)
      y_positions <- c(y_positions, y_position)
      all_y_position <- c(y_positions)
      if(length(filtereds[["annotations"]])>0){
        all_y_position<-all_y_position[1:length(filtereds[["annotations"]])]
        p1<-p1 + ggpubr::geom_signif(
          size = 0.25,
          tip_length =c(0.01,0.01),
          vjust = 1.0,
          y_position = all_y_position,
          xmin = filtereds[["x_min"]],
          xmax = filtereds[["x_max"]],
          annotation = filtereds[["annotations"]]
        )
        if(filtereds[["annotations"]]=="*"){
          cat("\033[31m", paste(i+1,":","Oh yeah!",ims[i+1],"有组间显著性差异(p<0.05)"), "\033[0m\n")
        }else if(filtereds[["annotations"]]=="**"){
          cat("\033[31m", paste(i+1,":","Oh yeah!",ims[i+1],"有组间显著性差异(p<0.01)"), "\033[0m\n")
        }else if(filtereds[["annotations"]]=="***"){
          cat("\033[31m", paste(i+1,":","Oh yeah!",ims[i+1],"有组间显著性差异(p<0.001)"), "\033[0m\n")
        }else if(filtereds[["annotations"]]=="****"){
          cat("\033[31m", paste(i+1,":","Oh yeah!",ims[i+1],"有组间显著性差异(p<0.0001)"), "\033[0m\n")
        }
      }else{cat("\033[32m", paste(i+1,":","Oopps!",ims[i+1],"没有组间显著性差异"), "\033[0m\n")}
      i=i+1
    }
  }else if(length(Type)==3){
    result1<-data_long[data_long$method==ims[i+1],]
    annotation1 <- get_annotation(result1, c(Type[[1]],Type[[2]]))
    annotation2 <- get_annotation(result1, c(Type[[1]],Type[[3]]))
    annotation3 <- get_annotation(result1, c(Type[[2]],Type[[3]]))
    all_annotation <- c(annotation1,annotation2,annotation3)
    all_x_min <- c(0.80,0.80,1.00)+i
    all_x_max <- c(1.00,1.20,1.20)+i
    filtereds <- filtered_annotation(all_annotation, all_x_min,all_x_max)
    y_positions <- c()
    myQuart1 =quantile(as.numeric(result1$method_score[result1$Group == Type[[1]]]), na.rm = TRUE)
    IQR1 =  myQuart1[4] - myQuart1[2]
    y_position1<-myQuart1[4] + 1.5 * IQR1
    myQuart2 =quantile(as.numeric(result1$method_score[result1$Group == Type[[2]]]), na.rm = TRUE)
    IQR2 =  myQuart2[4] - myQuart2[2]
    y_position2<-myQuart2[4] + 1.5 * IQR2
    myQuart3 =quantile(as.numeric(result1$method_score[result1$Group == Type[[3]]]), na.rm = TRUE)
    IQR3 =  myQuart3[4] - myQuart3[2]
    y_position3<-myQuart3[4] + 1.5 * IQR3
    y_position<-pmax(y_position1,y_position2,y_position3)
    y_positions <- c(y_positions, y_position)
    all_y_position <- c(y_positions)
    if(length(filtereds[["annotations"]])>0){
      all_y_position<-all_y_position[1:length(filtereds[["annotations"]])]
      p1<-p1 + ggpubr::geom_signif(
        size = 0.25,
        tip_length =c(0.01,0.01),
        vjust = 1.0,
        y_position = all_y_position,
        xmin = filtereds[["x_min"]],
        xmax = filtereds[["x_max"]],
        annotation = filtereds[["annotations"]]
      )
      cat("\033[31m", paste(i+1,":","Oh yeah!",ims[i+1],"有组间显著性差异"), "\033[0m\n")
    }else{cat("\033[32m", paste(i+1,":","Oopps!",ims[i+1],"没有组间显著性差异"), "\033[0m\n")}
    i=i+1
  }else if(length(Type)==4){
    result1<-data_long[data_long$method==ims[i+1],]
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
    myQuart1 =quantile(as.numeric(result1$method_score[result1$Group == Type[[1]]]), na.rm = TRUE)
    IQR1 =  myQuart1[4] - myQuart1[2]
    y_position1<-myQuart1[4] + 1.5 * IQR1
    myQuart2 =quantile(as.numeric(result1$method_score[result1$Group == Type[[2]]]), na.rm = TRUE)
    IQR2 =  myQuart2[4] - myQuart2[2]
    y_position2<-myQuart2[4] + 1.5 * IQR2
    myQuart3 =quantile(as.numeric(result1$method_score[result1$Group == Type[[3]]]), na.rm = TRUE)
    IQR3 =  myQuart3[4] - myQuart3[2]
    y_position3<-myQuart3[4] + 1.5 * IQR3
    myQuart4 =quantile(as.numeric(result1$method_score[result1$Group == Type[[4]]]), na.rm = TRUE)
    IQR4 =  myQuart3[4] - myQuart3[2]
    y_position4<-myQuart4[4] + 1.5 * IQR4
    y_position<-pmax(y_position1,y_position2,y_position3,y_position4)
    y_positions <- c(y_positions, y_position)
    all_y_position <- c(y_positions)
    if(length(filtereds[["annotations"]])>0){
      all_y_position<-all_y_position[1:length(filtereds[["annotations"]])]
      p1<-p1 + ggpubr::geom_signif(
        size = 0.25,
        tip_length =c(0.01,0.01),
        vjust = 1.0,
        y_position = all_y_position,
        xmin = filtereds[["x_min"]],
        xmax = filtereds[["x_max"]],
        annotation = filtereds[["annotations"]]
      )
      cat("\033[31m", paste(i+1,":","Oh yeah!",ims[i+1],"有组间显著性差异"), "\033[0m\n")
    }else{cat("\033[32m", paste(i+1,":","Oopps!",ims[i+1],"没有组间显著性差异"), "\033[0m\n")}
    i=i+1
  }else if(length(Type)==5){
    result1<-data_long[data_long$method==ims[i+1],]
    annotation1 <- get_annotation(result1, c(Type[[1]],Type[[2]]))
    annotation2 <- get_annotation(result1, c(Type[[1]],Type[[3]]))
    annotation3 <- get_annotation(result1, c(Type[[1]],Type[[4]]))
    annotation4 <- get_annotation(result1, c(Type[[1]],Type[[5]]))
    annotation5 <- get_annotation(result1, c(Type[[2]],Type[[3]]))
    annotation6 <- get_annotation(result1, c(Type[[2]],Type[[4]]))
    annotation7 <- get_annotation(result1, c(Type[[2]],Type[[5]]))
    annotation8 <- get_annotation(result1, c(Type[[3]],Type[[4]]))
    annotation9 <- get_annotation(result1, c(Type[[3]],Type[[5]]))
    annotation10 <- get_annotation(result1, c(Type[[4]],Type[[5]]))
    all_annotation <- c(annotation1,annotation2,annotation3,annotation4,annotation5,annotation6,
                        annotation7,annotation8,annotation9,annotation10)
    all_x_min <- c(0.68,0.68,0.68,0.68,0.84,0.84,0.84,1.00,1.00,1.16)+i
    all_x_max <- c(0.84,1.00,1.16,1.32,1.00,1.16,1.32,1.16,1.32,1.32)+i
    filtereds <- filtered_annotation(all_annotation, all_x_min,all_x_max)
    y_positions <- c()
    myQuart1 =quantile(as.numeric(result1$method_score[result1$Group == Type[[1]]]), na.rm = TRUE)
    IQR1 =  myQuart1[4] - myQuart1[2]
    y_position1<-myQuart1[4] + 1.5 * IQR1
    myQuart2 =quantile(as.numeric(result1$method_score[result1$Group == Type[[2]]]), na.rm = TRUE)
    IQR2 =  myQuart2[4] - myQuart2[2]
    y_position2<-myQuart2[4] + 1.5 * IQR2
    myQuart3 =quantile(as.numeric(result1$method_score[result1$Group == Type[[3]]]), na.rm = TRUE)
    IQR3 =  myQuart3[4] - myQuart3[2]
    y_position3<-myQuart3[4] + 1.5 * IQR3
    myQuart4 =quantile(as.numeric(result1$method_score[result1$Group == Type[[4]]]), na.rm = TRUE)
    IQR4 =  myQuart3[4] - myQuart3[2]
    y_position4<-myQuart4[4] + 1.5 * IQR4
    myQuart5 =quantile(as.numeric(result1$method_score[result1$Group == Type[[5]]]), na.rm = TRUE)
    IQR5 =  myQuart5[4] - myQuart5[2]
    y_position5<-myQuart5[4] + 1.5 * IQR5
    y_position<-pmax(y_position1,y_position2,y_position3,y_position4,y_position5)
    y_positions <- c(y_positions, y_position)
    all_y_position <- c(y_positions)
    if(length(filtereds[["annotations"]])>0){
      all_y_position<-all_y_position[1:length(filtereds[["annotations"]])]
      p1<-p1 + ggpubr::geom_signif(
        size = 0.25,
        tip_length =c(0.01,0.01),
        vjust = 1.0,
        y_position = all_y_position,
        xmin = filtereds[["x_min"]],
        xmax = filtereds[["x_max"]],
        annotation = filtereds[["annotations"]]
      )
      cat("\033[31m", paste(i+1,":","Oh yeah!",ims[i+1],"有组间显著性差异"), "\033[0m\n")
    }else{cat("\033[32m", paste(i+1,":","Oopps!",ims[i+1],"没有组间显著性差异"), "\033[0m\n")}
    i=i+1
  }
  if(geom_text==T){
    if(length(colnames(im[,-c(1,2)]))<=10){
    p1<-p1  +geom_text(data = data.frame(x = 0.8, y = max, label = tumor),
                       mapping = aes(x = x, y = y, label = label), inherit.aes = FALSE)
  }else if(length(colnames(im[,-c(1,2)]))>10&&length(colnames(im[,-c(1,2)]))<=30){
    p1<-p1  +geom_text(data = data.frame(x = 1.8, y = max, label = tumor),
                       mapping = aes(x = x, y = y, label = label), inherit.aes = FALSE)
  }else if(length(colnames(im[,-c(1,2)]))>30&&length(colnames(im[,-c(1,2)]))<=60){
    p1<-p1  +geom_text(data = data.frame(x = 2.8, y = max, label = tumor),
                       mapping = aes(x = x, y = y, label = label), inherit.aes = FALSE)
  }else if(length(colnames(im[,-c(1,2)]))>60){
    p1<-p1  +geom_text(data = data.frame(x = 3.8, y = max, label = tumor),
                       mapping = aes(x = x, y = y, label = label), inherit.aes = FALSE)
  }
  }
  end_time <- Sys.time()
  total_time <- end_time - start_time
  print(total_time)
  results <- list(
    imm_data = im,
    method = method,
    imm_plot = p1
  )
  return(results)
}
