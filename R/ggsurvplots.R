#' Drawing Survival Curves Using ggplot2
#'
#' @param data a dataset used to fit survival curves. If not supplied then data will be extracted from 'fit' object.
#' @param time_col the colname of time
#' @param status_col the colname of status
#' @param group_col the colname of sample group
#' @param pvalue_table table of Log-rank p test of samples
#' @param palette the color palette to be used. Allowed values include "hue" for the default hue color scale;
#'                "grey" for grey color palettes; brewer palettes e.g. "RdBu", "Blues", ...;
#'                or custom color palette e.g. c("blue", "red");
#'                and scientific journal palettes from ggsci R package, e.g.: "npg", "aaas", "lancet", "jco", "ucscgb", "uchicago", "simpsons" and "rickandmorty".
#'                See details section for more information. Can be also a numeric vector of length(groups);
#'                in this case a basic color palette is created using the function palette.
#' @param risk.table Allowed values include:
#'                   TRUE or FALSE specifying whether to show or not the risk table. Default is FALSE.
#'                  "absolute" or "percentage". Shows the absolute number and the percentage of subjects at risk by time, respectively.
#'                  "abs_pct" to show both absolute number and percentage.
#'                  "nrisk_cumcensor" and "nrisk_cumevents". Show the number at risk and, the cumulative number of censoring and events, respectively.
#' @param title     main title.
#' @param legend.labs character vector specifying legend labels. currently max length is 20.
#'                    Used to replace the names of the strata from the fit. Should be given in the same order as those strata.
#' @param xlab,ylab  xlab, ylab: x and y axis labels, respectively.
#' @param surv.median.line	 character vector for drawing a horizontal/vertical line at median survival. Allowed values include one of c("none", "hv", "h", "v"). v: vertical, h:horizontal.
#' @param surv.scale surv.scale: scale transformation of survival curves. Allowed values are "default" or "percent".
#' @param legend      lenged show or not, Default is TRUE.
#'
#' @return A ggplot object.
#' @export
##' @author Pengmin Yang
#' @examples
#' data("clin_TCGA")
#' p <- ggsurvplots(data = clin_TCGA, conf.int = FALSE,time_col = "PFS_MONTHS",
#'                  status_col = "PFS_STATUS", group_col = "Status", pvalue_table = TRUE,
#'                  palette = ggsci::pal_ucscgb()(4), risk.table = FALSE, title = NULL,
#'                  legend.labs = c("no KRAS or TP53", "TP53", "KRAS", "KRAS&TP53"),
#'                  xlab = "PFS_MONTHS", ylab = "Survival probability",
#'                  surv.median.line="hv",surv.scale="default",legend=FALSE)
#' print(p)
#'
ggsurvplots<-function(data,conf.int = FALSE, time_col="OS_MONTHS",
                      status_col="OS_SATUTS",group_col="Risk",pvalue_table=F,
                      palette=NULL,risk.table=T,title="MSk cohort",
                      legend.labs=c("High-risk","Low-risk"),xlab="OS_MONTHS",
                      ylab="Survial probability",surv.median.line="none",
                      surv.scale="default",legend=TRUE){
  library(survival)
  library(survminer)
  data<-data[!is.na(data[[time_col]]),]
  data[[time_col]]<-as.numeric(data[[time_col]])
  data[[status_col]]<-as.numeric(data[[status_col]])
  data[[group_col]] <- factor(data[[group_col]], levels = legend.labs)
  data <- data[order(data[[group_col]]), ]
  colnames(data)[colnames(data)==time_col]<-"Time"
  colnames(data)[colnames(data)==status_col]<-"Time_Status"
  colnames(data)[colnames(data)==group_col]<-"Group"
  fit <- survfit(Surv(time = Time, event = Time_Status) ~ Group, data =data)
  length<-length(unique(data$Group))
  p<-ggsurvplot(fit,
                data = data,
                pval = F, pval.method = F,
                conf.int = conf.int,
                risk.table = risk.table,
                tables.y.text= if (risk.table) T else F,
                tables.col="Group",
                risk.table.height = 0.25,
                palette = palette,
                legend.title = group_col,
                legend.labs = legend.labs,
                xlab=xlab,#if (risk.table) xlab=NULL else xlab,
                ylab=ylab ,
                title=title,
                surv.median.line = surv.median.line,
                surv.scale=surv.scale,
                ggtheme = theme_classic2()
  )
  if (risk.table){
    p[["plot"]][["labels"]][["x"]]<-" "
    p[["table"]][["labels"]][["y"]]<-NULL
  }
  median<-surv_median(fit)
  median<-median[order(median$median), ]
  data.survdiff <- survdiff(Surv(time = Time, event = Time_Status) ~ Group, data =data)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[2]/data.survdiff$exp[2])/(data.survdiff$obs[1]/data.survdiff$exp[1])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[2]+1/data.survdiff$exp[1]))
  ci <- paste0(sprintf("%.3f",HR)," [",sprintf("%.3f",low95),", ",sprintf("%.3f",up95),"]")
  if(surv.median.line%in%c("none", "hv", "h", "v")){
    if(length==2){
      if(!is.na(surv_median(fit)[["median"]][[1]])&&!is.na(surv_median(fit)[["median"]][[2]])){
        p$plot <- p$plot +
          annotate("text", x = max(data$Time)/2, y = 0.75,
                   label = if (p.val < 0.0001) paste0(" Log-rank\n p < 0.0001","\n HR (95% CI) = ",ci,"\n","n = ",length(data$Time_Status)) else paste0(" Log-rank\n","p = ", round(p.val,4),"\n HR (95% CI) = ",ci,"\n","n = ",length(data$Time_Status)),   ###添加P和HR 95%CI
                   size = 5, color = "black", hjust = 0,vjust=0)+
          theme(text = element_text(size = 15))+annotate("text", x = surv_median(fit)[["median"]][[1]], y = 0,
                                                         label = sprintf("%.3f",surv_median(fit)[["median"]][[1]]),   ###添加P和HR 95%CI
                                                         size = 5, color = "black", hjust =if(surv_median(fit)[["median"]][[1]]<surv_median(fit)[["median"]][[2]]) 1 else 0,vjust=0)+
          theme(text = element_text(size = 15))+annotate("text", x=surv_median(fit)[["median"]][[2]], y = 0,
                                                         label = sprintf("%.3f",surv_median(fit)[["median"]][[2]]),   ###添加P和HR 95%CI
                                                         size = 5, color = "black", hjust = if(surv_median(fit)[["median"]][[1]]<surv_median(fit)[["median"]][[2]]) 0 else 1,vjust=0)+
          theme(text = element_text(size = 15))
      }else if(is.na(surv_median(fit)[["median"]][[1]])&&!is.na(surv_median(fit)[["median"]][[2]])){
        p$plot <- p$plot +
          annotate("text", x = max(data$Time)/2, y = 0.75,
                   label = if (p.val < 0.0001) paste0(" Log-rank\n p < 0.0001","\n HR (95% CI) = ",ci,"\n","n = ",length(data$Time_Status)) else paste0(" Log-rank\n","p = ", round(p.val,4),"\n HR (95% CI) = ",ci,"\n","n = ",length(data$Time_Status)),   ###添加P和HR 95%CI
                   size = 5, color = "black", hjust = 0,vjust=0)+
          theme(text = element_text(size = 15))+annotate("text", x=surv_median(fit)[["median"]][[2]], y = 0,
                                                         label = sprintf("%.3f",surv_median(fit)[["median"]][[2]]),   ###添加P和HR 95%CI
                                                         size = 5, color = "black", hjust = 1,vjust=0)+
          theme(text = element_text(size = 15))
      }else if(!is.na(surv_median(fit)[["median"]][[1]])&&is.na(surv_median(fit)[["median"]][[2]])){
        p$plot <- p$plot +
          annotate("text", x = max(data$Time)/2, y = 0.75,
                   label = if (p.val < 0.0001) paste0(" Log-rank\n p < 0.0001","\n HR (95% CI) = ",ci) else paste0(" Log-rank\n","p = ", round(p.val,4),"\n HR (95% CI) = ",ci),   ###添加P和HR 95%CI
                   size = 5, color = "black", hjust = 0,vjust=0)+
          theme(text = element_text(size = 15))+annotate("text", x=surv_median(fit)[["median"]][[1]], y = 0,
                                                         label = sprintf("%.3f",surv_median(fit)[["median"]][[1]]),   ###添加P和HR 95%CI
                                                         size = 5, color = "black", hjust = 1,vjust=0)+
          theme(text = element_text(size = 15))
      }else if(is.na(surv_median(fit)[["median"]][[1]])&&is.na(surv_median(fit)[["median"]][[2]])){
        p$plot <- p$plot +
          annotate("text", x = max(data$Time)/2, y = 0.75,
                   label = if (p.val < 0.0001) paste0(" Log-rank\n p < 0.0001","\n HR (95% CI) = ",ci,"\n","n = ",length(data$Time_Status)) else paste0(" Log-rank\n","p = ", round(p.val,4),"\n HR (95% CI) = ",ci,"\n","n = ",length(data$Time_Status)),   ###添加P和HR 95%CI
                   size = 5, color = "black", hjust = 0,vjust=0)+
          theme(text = element_text(size = 15))
      }
    }else{
      calculate_survdiff_metrics <- function(data, group1, group2) {
        data_subset <- data[data$Group %in% c(group1, group2), ]
        data_survdiff <- survdiff(Surv(time = Time, event = Time_Status) ~ Group, data = data_subset)
        p_val <- 1 - pchisq(data_survdiff$chisq, length(data_survdiff$n) - 1)
        HR <- (data_survdiff$obs[2] / data_survdiff$exp[2]) / (data_survdiff$obs[1] / data_survdiff$exp[1])
        up95 <- exp(log(HR) + qnorm(0.975) * sqrt(1 / data_survdiff$exp[2] + 1 / data_survdiff$exp[1]))
        low95 <- exp(log(HR) - qnorm(0.975) * sqrt(1 / data_survdiff$exp[2] + 1 / data_survdiff$exp[1]))
        ci <- paste0(sprintf("%.3f", HR), " [", sprintf("%.3f", low95), ", ", sprintf("%.3f", up95), "]")
        return(c(sprintf("%.5f", p_val), ci))
      }
      # Prepare the groups for comparison
      group_combinations <- combn(legend.labs, 2, simplify = FALSE)
      # Initialize a results dataframe
      results <- data.frame(group = character(), pvalue = numeric(), HR_95_CI = character(), stringsAsFactors = FALSE)
      # Loop through the combinations to calculate metrics for each pair
      i=1
      for (groups in group_combinations) {
        metrics <- calculate_survdiff_metrics(data,  groups[1], groups[2])
        results <- rbind(results, data.frame(group = paste(i,groups[1], "vs.", groups[2]),
                                             pvalue = metrics[1],
                                             HR_95_CI = metrics[2],
                                             stringsAsFactors = FALSE))
        i=i+1
      }
      colnames(results)[colnames(results) == "HR_95_CI"] <-"HR(95% CI)"
      colnames(results)[colnames(results) == "pvalue"] <-"Log-rank pvalue"
      rownames(results)<-results$group
      results$group<-NULL
      rownames(results) <- gsub("^\\d+\\s+", "", rownames(results))
      library(gridExtra)
      library(gridExtra)
      library(grid)
      #results$type<-rownames(results)
      #rownames(results)<-NULL
      #results<-results[,c(3,1,2)]
      #colnames(results)[[1]]<-" "
      if(length<=5){
        results$type<-rownames(results)
        rownames(results)<-NULL
        results<-results[,c(3,1,2)]
        colnames(results)[[1]]<-" "
        fontface_vector <- rep("plain", nrow(results))  # 默认设置为 "plain"
        fontface_vector[which(results$`Log-rank pvalue` < 0.05)] <- "bold.italic"
        col_vector<-rep(NA, nrow(results))
        col_vector[which(results$`Log-rank pvalue` < 0.05)] <- "#6BAED6"
        hjust_vector<-rep(0.5, nrow(results))
        hjust_vector[which(results$`Log-rank pvalue` < 0.05)] <- 1
        x_vector<-rep(0.5, nrow(results))
        x_vector[which(results$`Log-rank pvalue` < 0.05)] <- 1
        tp <- gridExtra::tableGrob(results,
                                   theme =ttheme_minimal(
                                     core = list(fg_params = list(fontface = fontface_vector,cex = 0.7,hjust=hjust_vector, x=x_vector),
                                                 bg_params = list(fill = col_vector, col = NA) # 使用自定义字体面设置
                                     ),
                                     colhead = list(fg_params = list(col = "black", fontface = "bold",cex = 0.8)),
                                     rowhead = list(fg_params = list(col = "transparent", fontface = "plain",cex = 0.7),
                                                    bg_params = list(fill = "transparent", col = NA) )
                                   )
        )
      }else if(length>5){
        results<- results[results$`Log-rank pvalue` < 0.05,]
        tp <- gridExtra::tableGrob(results,
                                   theme =ttheme_minimal(
                                     core = list(fg_params = list(fontface = "plain",cex = 0.7,hjust=0.5),
                                                 bg_params = list(fill = NA, col = NA) # 使用自定义字体面设置
                                     ),
                                     colhead = list(fg_params = list(col = "black", fontface = "bold",cex = 0.8)),
                                     rowhead = list(fg_params = list(col = "black",fontface = "bold",cex = 0.8)
                                     )
                                   ))
      }
      median$strata<-gsub("Group=","",median$strata)
      add_annotations <- function(p, p_val, median, length) {
        p <- p + annotate("text", x = 0, y = 0.3,
                          label = if (p_val < 0.0001) paste0("p < 0.0001") else paste0("p = ", round(p_val,4)), #sprintf("%.4f", p_val)
                          size = 5, color =if (p_val < 0.05)"red" else "black", hjust = 0, vjust = 0) +
          theme(text = element_text(size = 15))

        x=c()
        y=c()
        label=c()
        hjust=c()
        palettes=c()
        for (i in seq_len(length)) {
          y_position <- if (i <= 4) {
            if (i %% 2 == 0) 0.1 else 0.05
          } else if (i >= 5 && i <= 8) {
            if (i %% 2 == 0) 0.20 else 0.15
          } else if (i >= 9 && i <= 12) {
            if (i %% 2 == 0) 0.30 else 0.25
          } else if (i >= 13 && i <= 16) {
            if (i %% 2 == 0) 0.40 else 0.35
          } else if (i >= 17 && i <= 20) {
            if (i %% 2 == 0) 0.10 else 0.05
          }
          x_i = median[["median"]][[i]]
          hjust_i <- if ((i - 1) %/% 2 %% 2 == 0) 1 else 0
          label_i=sprintf("%.3f", median[["median"]][[i]])
          color=data.frame(palette = palette,
                           legend.labs = legend.labs)
          palettes_i=color[color$legend.labs==median[["strata"]][[i]],"palette"]
          x=c(x,x_i)
          y=c(y,y_position)
          label=c(label,label_i)
          hjust=c(hjust,hjust_i)
          palettes=c(palettes,palettes_i)
        }
        print(palettes)
        if(surv.median.line%in%c("hv", "v")){
          p <- p + annotate("text", x =x, y = y,
                            label = label,
                            size = 5, color = palettes, hjust = hjust, vjust = 0) +
            theme(text = element_text(size = 15))

          return(p)
        }
      }

      # 使用简化后的函数添加注释
      if (length >= 3) {
        p$plot <- add_annotations(p$plot, p.val, median, length)
      }
      if(pvalue_table){
        if(surv.median.line%in%c("hv", "v")){
          p$plot<-p$plot+ annotation_custom(tp,
                                            xmin=max(data$Time)*0.50,
                                            xmax=max(data$Time)*0.90,
                                            ymin=0.5,ymax=1.0)
        }
      }
    }
    if(risk.table){
      library(patchwork)
      theme <- theme(plot.title = element_text(hjust = 0.5,size=20),
                     axis.text.x = element_text(hjust = 0.5,size=16,colour = "black"),
                     axis.text.y = element_text(size = 16,colour = "black"),
                     axis.title.x = element_text(size = 18),
                     axis.title.y = element_text(size = 18),
                     axis.line = element_line(colour = "black",size = 1),
                     legend.text = element_text(colour = "black",size = 12),
                     #legend.title = element_blank(),
                     legend.position = if(pvalue_table) "top" else c(0.8, 0.8))
      p$plot <- p$plot + theme+ theme(axis.ticks.x = element_blank(), ## 删去所有刻度线
                                      axis.text.x = element_blank(),
                                      axis.title.x = element_blank(),panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank())
      p$table <- p$table + theme(plot.title = element_blank(),
                                 legend.position = "none")
      plot.down <- p$table+ theme(axis.ticks.x = element_blank(), ## 删去所有刻度线
                                  axis.text.x = element_blank(),
                                  axis.title.x = element_blank(),panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank())
      p$plot+plot.down + plot_layout(ncol = 1,heights = c(5,2.5),widths = c(2, 1.3))##调整中间空白
    }else {
      p$plot <- p$plot+theme(legend.position = if(pvalue_table) "top" else c(0.8, 0.8))
    }
    if(legend){
      p<-p
    }else{
      x=c()
      y=c()
      label=c()
      for(id in legend.labs){
        x_id<-p[["data.survplot"]][p[["data.survplot"]]$Group==id,1][length(p[["data.survplot"]][p[["data.survplot"]]$Group==id,1])]
        y_id<-p[["data.survplot"]][p[["data.survplot"]]$Group==id,5][length(p[["data.survplot"]][p[["data.survplot"]]$Group==id,5])]
        label_id<-id
        x=c(x,x_id)
        y=c(y,y_id)
        label=c(label,label_id)
      }
      data = data.frame(x=x,y=y,label=label)
      p$plot<-p$plot+ggrepel::geom_label_repel(data =data,max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                                               aes(x, y, label = label),colour = palette,
                                               arrow = arrow(length = unit(0.02, "npc"),
                                                             ends = "first",
                                                             type = "closed"),
                                               box.padding = 1,
                                               size = 5, fill="white")+theme(
                                                 legend.position = "none")
    }
    p$plot <- p$plot+theme(axis.text.x = element_text(hjust = 0.5,size=16,colour = "black"),
                           axis.text.y = element_text(size = 16,colour = "black"),
                           axis.title.x = element_text(size = 18),
                           axis.title.y = element_text(size = 18),
                           axis.line = element_line(colour = "black",linewidth = 1),
                           legend.text = element_text(colour = "black",size = 12))
    return(p)
  }else{print("surv.median.line should be in one of c(\"none\", \"hv\", \"h\", \"v\")")}
}
