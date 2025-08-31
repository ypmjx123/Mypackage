#' Immune score Visualization
#'
#' This function visualizes the mutation rate of immune score based on provided immune score or gene expression data.
#'
#' @param im           is the data of immune score for samples;if exp=NULL,im must be provided
#' @param exp          is the data of gene expression for samples;please ensure gene coloum with unique values;if im=NULL,exp must be provided
#' @param method       IOBR::deconvo_tme 参数或者根据其它方法得到的样本免疫浸润评分数据; method =c('mcpcounter', 'epic',  'cibersort', 'estimate','ips','quantiseq','xCell','timer')
#' @param exp_cleaning/im_cleaning 是否对数据进行清洗，默认FALSE
#' @param scale_mrna   logical. If FALSE, disable correction for mRNA content of different cell types. This is supported by methods that compute an absolute score (EPIC and quanTIseq)
#' @param tumor.TCGA    method='timer',group_list参数，这里导入的group_list必须是TIMER能够识别的，比如TCGA中33肿瘤类型;Accepted indications are 'kich', 'blca', 'brca', 'cesc', 'gbm', 'hnsc', 'kirp', 'lgg','lihc', 'luad', 'lusc', 'prad', 'sarc', 'pcpg', 'paad', 'tgct','ucec', 'ov', 'skcm', 'dlbc', 'kirc', 'acc', 'meso', 'thca','uvm', 'ucs', 'thym', 'esca', 'stad', 'read', 'coad', 'chol'
#' @param sample_group is a dataframe for sample grouping
#' @param heatmap         defalut heatmap=FALSE
#' @param heatmap_col    defalut heatmap_col=NULL
#' @param tumor        data source, "CRC TCGA"
#' @param Type         can length from 2 to 5,eg Type=c("Wild", "Mut");(default = Type[[1]] is the control)
#' @param color        for plot of Type, eg color=c("#757575", "#FF4040")
#' @param geom_text    text annotation for ggplot
#' @param test         参考ggpubr::compare_means，根据需要选择
#' @return A list of immune score for data samples,method and heatmap of immune score for samples/ggplot object visualizing the immune score.
#' @export
##' @author Pengmin Yang
#' @examples
#' #im=NULL
#' data("exp_raw")
#' data("Gene_group_CRC1")
#' exp_CRC<-exp_geneIDtoSYMBOL(exp=exp_raw,genecoltype="ENTREZID")
#' exp_CRC<-exp_CRC$data
#' result<-immu_visual(im=NULL,exp=exp_CRC[,-1],
#'                     method = 'epic',
#'                     sample_group=Gene_group_CRC1,
#'                     tumor="CRC TCGA",heatmap=TRUE,
#'                     Type=c("Wild", "Mut"),
#'                     color=c("#757575", "#FF4040"),
#'                     geom_text=TRUE,
#'                     test = "wilcox.test")
#' print(result)
#'
#' #exp=NULL
#' data("im")
#' data("sample_group")
#' result<-immu_visual(im=im,exp=NULL,
#'                     method = ' ',
#'                     sample_group=sample_group,
#'                     tumor="GC TCGA",heatmap=TRUE,
#'                     Type=c("Wild", "Mut"),
#'                     color=c("#757575", "#FF4040"),
#'                     geom_text=FALSE,
#'                     test = "wilcox.test")
#' print(result)
immu_visual<-function(im,exp,method,exp_cleaning=FALSE,im_cleaning=FALSE,scale_mrna,tumor.TCGA,
                      sample_group,heatmap=FALSE,heatmap_col=NULL, tumor,Type,color,geom_text,test){
  start_time <- Sys.time()
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(ggpubr)
  library(tools)
  if(is.null(im)&&!is.null(exp)){
    library(IOBR)
    exp<-exp[!is.na(exp[,2]),]
    if(!is.na(as.numeric(rownames(exp)[[length(rownames(exp))]]))){
      if (any(duplicated(exp[[1]]))) {
        exp <- aggregate(as.formula(paste0(".", " ~ ", colnames(exp)[[1]])), data = exp, FUN = mean)
      }
      rownames(exp)<-exp[[1]]
      exp<-exp[,-1]
    }
    if(exp_cleaning){
      n <- ncol(exp)  # 获取总列数
      pseudo_count <- min(c(exp[,2:n][exp[,2:n] > 0]), na.rm = TRUE) / 2  # 添加na.rm=TRUE处理NA值
      exp[,2:n] <- log10(exp[,2:n] + pseudo_count)
    }
    if(method =='quantiseq'){
      im <- deconvo_tme(eset = exp,
                        method = 'quantiseq',
                        arrays = F,
                        perm = 1000,scale_mrna=scale_mrna
      )
    }else if(method =='xCell'){
      library(xCell)
      im <- xCellAnalysis(exp)
      im<-as.data.frame(t(im))
      im$ID<-rownames(im)
      rownames(im)<-NULL
      im<-im[,c(68,1:67)]
    }else if(method =='timer'){
      im <- deconvo_tme(eset = exp,
                        method = 'timer',
                        arrays = F,
                        perm = 1000,group_list = rep(tumor.TCGA,dim(exp)[2])
      )
    }else{
      im <- deconvo_tme(eset = exp,
                        method = method,
                        arrays = F,
                        perm = 1000
      )
    }
    im<-merge(sample_group,im,by.x=colnames(sample_group)[[1]],by.y=colnames(im)[[1]],all.y=T)
    im[,2][is.na(im[, 2])] <- Type[[1]]
    if(im_cleaning){
      n <- ncol(im)
      pseudo_count <- min(c(im[,3:n][im[,3:n] > 0])) / 2
      im[,3:n] <- log10(im[,3:n] + pseudo_count)
    }
    if(method=='cibersort'){
      im1<-im[,-c(length(colnames(im))-2,length(colnames(im))-1,length(colnames(im)))]
      data_long <- im1%>%
        pivot_longer(cols = c(3:length(colnames(im1))), names_to = method, values_to = paste0(method,"_score"))
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
    }else if(method=='ips'){
      data_long <- im%>%
        pivot_longer(cols = c(3:length(colnames(im))), names_to = method, values_to = paste0(method,"_score"))
      data_long[[method]] <- gsub("_IPS$", "", data_long[[method]])
      data_long[[method]] <- gsub("_", " ", data_long[[method]])
    }else if(method=='quantiseq'){
      data_long <- im%>%
        pivot_longer(cols = c(3:length(colnames(im))), names_to = method, values_to = paste0(method,"_score"))
      data_long[[method]] <- gsub("_quantiseq$", "", data_long[[method]])
      data_long[[method]] <- gsub("_", " ", data_long[[method]])
    }else if(method=='xCell'){
      im1<-im[,-c(length(colnames(im))-2,length(colnames(im))-1,length(colnames(im)))]
      data_long <- im1%>%
        pivot_longer(cols = c(3:length(colnames(im1))), names_to = method, values_to = paste0(method,"_score"))
      data_long[[method]] <- gsub("_xCell$", "", data_long[[method]])
      data_long[[method]] <- gsub("_", " ", data_long[[method]])
    }else if(method=='timer'){
      data_long <- im%>%
        pivot_longer(cols = c(3:length(colnames(im))), names_to = method, values_to = paste0(method,"_score"))
      data_long[[method]] <- gsub("_TIMER$", "", data_long[[method]])
      data_long[[method]] <- gsub("_", " ", data_long[[method]])
    }else if(method=='estimate'){
      im1<-im[,-c(length(colnames(im)))]
      data_long <- im1 %>%
        pivot_longer(cols = c(3:length(colnames(im1))), names_to = method, values_to = paste0(method,"_score"))
      data_long[[method]] <- gsub("_estimate$", "", data_long[[method]])
      data_long[[method]] <- gsub("_", " ", data_long[[method]])
    }else{
      data_long <- im%>%
        pivot_longer(cols = c(3:length(colnames(im))), names_to = method, values_to = paste0(method,"_score"))
      data_long[[method]] <- gsub(paste0("_",method,"$"), "", data_long[[method]])
      data_long[[method]] <- gsub("_", " ", data_long[[method]])
    }
  }else if (!is.null(im)&&is.null(exp)){
    im<-merge(sample_group,im,by.x=colnames(sample_group)[[1]],by.y=colnames(im)[[1]],all.y=T)
    im[,2][is.na(im[, 2])] <- Type[[1]]
    if(im_cleaning){
      n <- ncol(im)
      pseudo_count <- min(c(im[,3:n][im[,3:n] > 0])) / 2
      im[,3:n] <- log10(im[,3:n] + pseudo_count)
    }
    if(method=='cibersort'){
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
    }else if(method=='ips'){
      data_long <- im%>%
        pivot_longer(cols = c(3:length(colnames(im))), names_to = method, values_to = paste0(method,"_score"))
      data_long[[method]] <- gsub("_IPS$", "", data_long[[method]])
      data_long[[method]] <- gsub("_", " ", data_long[[method]])
    }else if(method=='quantiseq'){
      data_long <- im%>%
        pivot_longer(cols = c(3:length(colnames(im))), names_to = method, values_to = paste0(method,"_score"))
      data_long[[method]] <- gsub("_quantiseq$", "", data_long[[method]])
      data_long[[method]] <- gsub("_", " ", data_long[[method]])
    }else if(method=='xCell'){
      data_long <- im%>%
        pivot_longer(cols = c(3:length(colnames(im))), names_to = method, values_to = paste0(method,"_score"))
      data_long[[method]] <- gsub("_xCell$", "", data_long[[method]])
      data_long[[method]] <- gsub("_", " ", data_long[[method]])
    }else if(method=='timer'){
      data_long <- im%>%
        pivot_longer(cols = c(3:length(colnames(im))), names_to = method, values_to = paste0(method,"_score"))
      data_long[[method]] <- gsub("_TIMER$", "", data_long[[method]])
      data_long[[method]] <- gsub("_", " ", data_long[[method]])
    }else if(method=='estimate'){
      data_long <- im%>%
        pivot_longer(cols = c(3:length(colnames(im))), names_to = method, values_to = paste0(method,"_score"))
      data_long[[method]] <- gsub("_estimate$", "", data_long[[method]])
      data_long[[method]] <- gsub("_", " ", data_long[[method]])
    }else{
      data_long <- im%>%
        pivot_longer(cols = c(3:length(colnames(im))), names_to = method, values_to = paste0(method,"_score"))
      data_long[[method]] <- gsub(paste0("_",method,"$"), "", data_long[[method]])
      data_long[[method]] <- gsub("_", " ", data_long[[method]])
    }
  }
  data_long[[colnames(data_long)[2]]]<-as.factor(data_long[[colnames(data_long)[2]]])
  # 去除空格
  data_long[[paste0(method,"_score")]] <- trimws(data_long[[paste0(method,"_score")]])
  # 转为数值
  data_long[[paste0(method,"_score")]] <- as.numeric(data_long[[paste0(method,"_score")]])
  data_long <- data_long[!is.na(data_long[[paste0(method,"_score")]]),]
  if(heatmap){
    input <- data_long %>%
      pivot_wider(names_from = !!sym(method), values_from = paste0(method, "_score"))
    feas = colnames(input)[3:length(colnames(input))]
    if(method=='cibersort'){Method="CIBERSORT"}
    else if (method=='mcpcounter'){Method="MCPcounter"}
    else if (method=='epic'){Method="EPIC"}
    else if (method=='ips'){Method="IPS"}
    else if (method=='quantiseq'){Method="quanTIseq"}
    else if (method=='xCell'){Method="xCell"}
    else if (method=='timer'){Method="TIMER"}
    else if (method=='estimate'){Method="ESTIMATE"}
    else {Method=method}
    p<- Mypackage::sig_Heatmap(input = input, features = feas,ID =colnames(input)[[1]],show_plot=F,heatmap_col=NULL,
                               cols_group= color,
                               group = colnames(sample_group)[[2]], scale = TRUE,column_title=NULL,column_title_size=8,
                               row_title=Method,
                               row_title_size=8,name="Score")
  }else(p<-NULL)
  data_long[[colnames(sample_group)[[2]]]] <- factor(data_long[[colnames(sample_group)[[2]]]], levels = Type)
  data_long[[method]] <- factor(data_long[[method]] , levels = unique(data_long[[method]]))
  colnames(data_long)[colnames(data_long) == method] <-"method"
  colnames(data_long)[colnames(data_long) == paste0(method,"_score")] <-"method_score"
  colnames(data_long)[colnames(data_long) == colnames(data_long)[[2]]] <-"Group"
  result0 <- data_long %>%
    group_by(method, Group) %>%
    mutate(
      Q1 = quantile(method_score, 0.25),  # 第一四分位数
      Q3 = quantile(method_score, 0.75),  # 第三四分位数
      IQR = Q3 - Q1,                      # 四分位距
      LowerLimit = Q1 - 1.5 * IQR,        # 下限
      UpperLimit = Q3 + 1.5 * IQR         # 上限
    ) %>%
    filter(method_score >= LowerLimit & method_score <= UpperLimit)
  max<-max(result0$method_score)
  print(result0)
  p1<-ggplot(data_long, aes(x=method, y=method_score, fill=Group)) +
    geom_boxplot(stat = "boxplot",outliers = F,
                 position = "dodge2")+
    theme(legend.position ="top",
          axis.text.y = element_text(angle = 0,color = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(linewidth = 1),
          panel.background = element_rect(fill = "transparent"))+
    scale_fill_manual(values=color)
  if(length(colnames(im[,-c(1,2)]))<=5){
    p1<-p1+theme(axis.text.x = element_text(angle = 0,color = "black"))
  }else{
    p1<-p1+theme(axis.text.x = element_text(angle = 45,hjust = 1,color = "black"))
  }
  if(method=='mcpcounter'){
    p1<-p1+labs(x = "", y ="MCPcounter Score",fill = paste(colnames(sample_group)[[2]],"Group"))
  }else if(method=='xCell'){
    p1<-p1+labs(x = "", y ="xCell Score",fill = paste(colnames(sample_group)[[2]],"Group"))
  }else if(method %in% c('ips','estimate','timer','cibersort','epic')){
    p1<-p1+labs(x = "", y =paste(toupper(method),toTitleCase("score")),fill = paste(colnames(sample_group)[[2]],"Group"))
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
  ims <- unique(data_long$method)
  A<-as.numeric(length(ims))-1
  i=0
  p_0.05<-c()
  color_vector<-rep("black", length(ims))
  if(length(Type)==2){
    while(i<=A ){
      result1<-data_long[data_long$method==ims[i+1],]
      if(length(unique(result1$Group))>=2){
        annotation1 <- get_annotation(result1, c(Type[[1]],Type[[2]]))
        all_annotation <- c(annotation1)
        all_x_min <- c(0.80)+i
        all_x_max <- c(1.20)+i
        filtereds <- filtered_annotation(all_annotation, all_x_min,all_x_max)
        y_positions <- c()
        y_positions <- c()
        result <- result1 %>%
          group_by(method, Group) %>%
          mutate(
            Q1 = quantile(method_score, 0.25),  # 第一四分位数
            Q3 = quantile(method_score, 0.75),  # 第三四分位数
            IQR = Q3 - Q1,                      # 四分位距
            LowerLimit = Q1 - 1.5 * IQR,        # 下限
            UpperLimit = Q3 + 1.5 * IQR         # 上限
          ) %>%
          filter(method_score >= LowerLimit & method_score <= UpperLimit)
        y_position<-max(result$method_score)+0.040 * (max(result0$method_score, na.rm = TRUE) - min(result0$method_score, na.rm = TRUE))
        y_positions <- c(y_positions, y_position)
        all_y_position <- c(y_positions)
      }else{
        all_annotation <- c("NS.")
        all_x_min <- i
        all_x_max <- i
        filtereds <- filtered_annotation(all_annotation, all_x_min,all_x_max)
        y_positions <- c()
        y_positions <- c()
        result <- result1 %>%
          group_by(method, Group) %>%
          mutate(
            Q1 = quantile(method_score, 0.25),  # 第一四分位数
            Q3 = quantile(method_score, 0.75),  # 第三四分位数
            IQR = Q3 - Q1,                      # 四分位距
            LowerLimit = Q1 - 1.5 * IQR,        # 下限
            UpperLimit = Q3 + 1.5 * IQR         # 上限
          ) %>%
          filter(method_score >= LowerLimit & method_score <= UpperLimit)
        y_position1<-max(result$method_score)+0.040 * (max(result0$method_score, na.rm = TRUE) - min(result0$method_score, na.rm = TRUE))
        y_positions <- c(y_positions, y_position1)
        all_y_position <- c(y_positions)
      }
      if(length(filtereds[["annotations"]])>0){
        all_y_position<-all_y_position[1:length(filtereds[["annotations"]])]
        p1<-p1 + ggpubr::geom_signif(
          size = 0.25,
          tip_length =c(0.01,0.01),
          vjust = 0.6,
          y_position = all_y_position,
          xmin = filtereds[["x_min"]],
          xmax = filtereds[["x_max"]],
          annotation = filtereds[["annotations"]]
        )
        if(filtereds[["annotations"]]=="*"){
          cat("\033[31m", paste(i+1,":","Oh yeah!",ims[i+1],"有组间显著性差异(p<0.05)"), "\033[0m\n")
          p_0.05<-c(p_0.05,i+1)
        }else if(filtereds[["annotations"]]=="**"){
          cat("\033[31m", paste(i+1,":","Oh yeah!",ims[i+1],"有组间显著性差异(p<0.01)"), "\033[0m\n")
          p_0.05<-c(p_0.05,i+1)
        }else if(filtereds[["annotations"]]=="***"){
          cat("\033[31m", paste(i+1,":","Oh yeah!",ims[i+1],"有组间显著性差异(p<0.001)"), "\033[0m\n")
          p_0.05<-c(p_0.05,i+1)
        }else if(filtereds[["annotations"]]=="****"){
          cat("\033[31m", paste(i+1,":","Oh yeah!",ims[i+1],"有组间显著性差异(p<0.0001)"), "\033[0m\n")
          p_0.05<-c(p_0.05,i+1)
        }
      }else{cat("\033[32m", paste(i+1,":","Oopps!",ims[i+1],"没有组间显著性差异"), "\033[0m\n")}

      i=i+1
    }
  }else if(length(Type)==3){
    while(i<=A ){
      result1<-data_long[data_long$method==ims[i+1],]
      annotation1 <- get_annotation(result1, c(Type[[1]],Type[[2]]))
      annotation2 <- get_annotation(result1, c(Type[[1]],Type[[3]]))
      annotation3 <- get_annotation(result1, c(Type[[2]],Type[[3]]))
      all_annotation <- c(annotation1,annotation2,annotation3)
      all_x_min <- c(0.80,0.80,1.00)+i
      all_x_max <- c(1.00,1.20,1.20)+i
      filtereds <- filtered_annotation(all_annotation, all_x_min,all_x_max)
      y_positions <- c()
      y_positions <- c()
      result <- result1 %>%
        group_by(method, Group) %>%
        mutate(
          Q1 = quantile(method_score, 0.25),  # 第一四分位数
          Q3 = quantile(method_score, 0.75),  # 第三四分位数
          IQR = Q3 - Q1,                      # 四分位距
          LowerLimit = Q1 - 1.5 * IQR,        # 下限
          UpperLimit = Q3 + 1.5 * IQR         # 上限
        ) %>%
        filter(method_score >= LowerLimit & method_score <= UpperLimit)
      y_position<-max(result$method_score)+0.050 * (max(result0$method_score, na.rm = TRUE) - min(result0$method_score, na.rm = TRUE))
      y_positions <- c(y_positions, y_position)
      all_y_position <- c(y_positions)
      if(length(filtereds[["annotations"]])>0){
        all_y_position<-all_y_position[1:length(filtereds[["annotations"]])]
        p1<-p1 + ggpubr::geom_signif(
          size = 0.25,
          tip_length =c(0.01,0.01),
          vjust = 0.6,
          y_position = all_y_position,
          xmin = filtereds[["x_min"]],
          xmax = filtereds[["x_max"]],
          annotation = filtereds[["annotations"]]
        )
        cat("\033[31m", paste(i+1,":","Oh yeah!",ims[i+1],"有组间显著性差异"), "\033[0m\n")
        p_0.05<-c(p_0.05,i+1)
      }else{cat("\033[32m", paste(i+1,":","Oopps!",ims[i+1],"没有组间显著性差异"), "\033[0m\n")}
      i=i+1
    }
  }else if(length(Type)==4){
    while(i<=A ){
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
      result <- result1 %>%
        group_by(method, Group) %>%
        mutate(
          Q1 = quantile(method_score, 0.25),  # 第一四分位数
          Q3 = quantile(method_score, 0.75),  # 第三四分位数
          IQR = Q3 - Q1,                      # 四分位距
          LowerLimit = Q1 - 1.5 * IQR,        # 下限
          UpperLimit = Q3 + 1.5 * IQR         # 上限
        ) %>%
        filter(method_score >= LowerLimit & method_score <= UpperLimit)
      y_position<-max(result$method_score)+0.050 * (max(result0$method_score, na.rm = TRUE) - min(result0$method_score, na.rm = TRUE))
      y_positions <- c(y_positions, y_position)
      all_y_position <- c(y_positions)
      if(length(filtereds[["annotations"]])>0){
        all_y_position<-all_y_position[1:length(filtereds[["annotations"]])]
        p1<-p1 + ggpubr::geom_signif(
          size = 0.25,
          tip_length =c(0.01,0.01),
          vjust = 0.6,
          y_position = all_y_position,
          xmin = filtereds[["x_min"]],
          xmax = filtereds[["x_max"]],
          annotation = filtereds[["annotations"]]
        )
        cat("\033[31m", paste(i+1,":","Oh yeah!",ims[i+1],"有组间显著性差异"), "\033[0m\n")
        p_0.05<-c(p_0.05,i+1)
      }else{cat("\033[32m", paste(i+1,":","Oopps!",ims[i+1],"没有组间显著性差异"), "\033[0m\n")}
      i=i+1
    }
  }else if(length(Type)==5){
    while(i<=A ){
      result1<-data_long[data_long$method==ims[i+1],]
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
      y_positions <- c()
      result <- result1 %>%
        group_by(method, Group) %>%
        mutate(
          Q1 = quantile(method_score, 0.25),  # 第一四分位数
          Q3 = quantile(method_score, 0.75),  # 第三四分位数
          IQR = Q3 - Q1,                      # 四分位距
          LowerLimit = Q1 - 1.5 * IQR,        # 下限
          UpperLimit = Q3 + 1.5 * IQR         # 上限
        ) %>%
        filter(method_score >= LowerLimit & method_score <= UpperLimit)
      y_position<-max(result$method_score)+0.050 * (max(result0$method_score, na.rm = TRUE) - min(result0$method_score, na.rm = TRUE))
      y_positions <- c(y_positions, y_position)
      all_y_position <- c(y_positions)
      if(length(filtereds[["annotations"]])>0){
        all_y_position<-all_y_position[1:length(filtereds[["annotations"]])]
        p1<-p1 + ggpubr::geom_signif(
          size = 0.25,
          tip_length =c(0.01,0.01),
          vjust = 0.6,
          y_position = all_y_position,
          xmin = filtereds[["x_min"]],
          xmax = filtereds[["x_max"]],
          annotation = filtereds[["annotations"]]
        )
        cat("\033[31m", paste(i+1,":","Oh yeah!",ims[i+1],"有组间显著性差异"), "\033[0m\n")
        p_0.05<-c(p_0.05,i+1)
      }else{cat("\033[32m", paste(i+1,":","Oopps!",ims[i+1],"没有组间显著性差异"), "\033[0m\n")}
      i=i+1
    }
  }
  color_vector[p_0.05] <- "red"
  if(length(colnames(im[,-c(1,2)]))<=5){
    p1<-p1+theme(axis.text.x = element_text(angle = 0,color = color_vector))
  }else{
    p1<-p1+theme(axis.text.x = element_text(angle = 45,hjust = 1,color = color_vector))
  }
  if(geom_text){
    if(length(colnames(im[,-c(1,2)]))<=10){
      p1<-p1  +geom_text(data = data.frame(x = 0.8, y = max, label = tumor),
                         mapping = aes(x = x, y = y, label = label), inherit.aes = FALSE)
    }else if(length(colnames(im[,-c(1,2)]))>10&&length(colnames(im[,-c(1,2)]))<=30){
      p1<-p1  +geom_text(data = data.frame(x = 2.2, y = max, label = tumor),
                         mapping = aes(x = x, y = y, label = label), inherit.aes = FALSE)
    }else if(length(colnames(im[,-c(1,2)]))>30&&length(colnames(im[,-c(1,2)]))<=60){
      p1<-p1  +geom_text(data = data.frame(x = 3.2, y = max, label = tumor),
                         mapping = aes(x = x, y = y, label = label), inherit.aes = FALSE)
    }else if(length(colnames(im[,-c(1,2)]))>60){
      p1<-p1  +geom_text(data = data.frame(x = 4.2, y = max, label = tumor),
                         mapping = aes(x = x, y = y, label = label), inherit.aes = FALSE)
    }
  }
  end_time <- Sys.time()
  total_time <- end_time - start_time
  print(total_time)
  results <- list(
    imm_data = im,
    heatmap=p,
    method = method,
    imm_plot = p1
  )
  return(results)
}

