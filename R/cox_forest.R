#' Title
#' Cox Proportional Hazards Univariate and Multivariate Forest Plot Generator
#'
#' @param data A data.frame containing survival data and predictors.
#' @param time_col Character. Column name that stores the survival (time-to-event) data. Default = "OS_MONTHS".
#' @param status_col Character. Column name that stores the event indicator (0 = censored, 1 = event). Default = "OS_STATUS".
#' @param Univariate Logical. Whether to perform univariate Cox regression for each predictor. Default = TRUE.
#' @param univar_predictors Character vector. Variables to include in univariate analysis. Non-NULL columns must exist in data. Default = selected columns from an object called clin.
#' @param Multivariate Logical. Whether to perform multivariate Cox regression using all specified predictors. Default = TRUE.
#' @param multivar_predictors Character vector. Variables to include in multivariate analysis. Non-NULL columns must exist in data. Default = selected columns from an object called clin.
#' @param show_plots Logical. If TRUE, forest plots are produced for both univariate and multivariate results. Default = TRUE.
#' @param xticks1 Numeric vector. Custom x-axis tick positions for the univariate forest plot. If NULL, ticks are auto-generated. Default = NULL.
#' @param xticks2 Numeric vector. Custom x-axis tick positions for the multivariate forest plot. If NULL, ticks are auto-generated. Default = NULL.
#' @param title_univar Character. Title for the univariate forest plot. Default = "OS Univariate".
#' @param title_multivar Character. Title for the multivariate forest plot. Default = "OS Multivariate".
#' @param use_baseline_table Logical. Whether to generate baseline characteristic tables for the selected predictors. Default = TRUE.
#' @param all Logical.show missing data or not. Default = TRUE.
#' @param forestplot plot with R package 'forestplot' or 'forestploter'. Default = TRUE.
#' @param ci_pch the shape of ci frame
#' @param ci_col the color of ci frame
#' @param ci_line the color of ci line
#' @param zero_col the color of reference line
#' @param log2    Default = FALSE
#' @param footnote Default = NULL
#'
#' @returns A named list with four top-level elements:
#' \describe{
#'   \item{univariate}{A list containing:
#'     \itemize{
#'       \item \code{univariate}: data.frame of univariate Cox results (HR, 95%CI, p-value, N(%)).
#'       \item \code{uniforest}: \code{forestplot} object (or NULL if \code{show_plots = FALSE}).
#'       \item \code{baseline_table}: data.frame summarising baseline characteristics used in univariate analysis (or NULL if \code{use_baseline_table = FALSE}).
#'     }
#'   }
#'   \item{multivariate}{A list containing:
#'     \itemize{
#'       \item \code{multivariate}: data.frame of multivariate Cox results (HR, 95%CI, p-value, N(%)).
#'       \item \code{multiforest}: \code{forestplot} object (or NULL if \code{show_plots = FALSE}).
#'       \item \code{baseline_table}: data.frame summarising baseline characteristics used in multivariate analysis (or NULL if \code{use_baseline_table = FALSE}).
#'     }
#'   }
#'   \item{data_used}{The cleaned data.frame that was actually used in the analyses (factors converted, NAs handled).}
#' }
#' @export
##' @author Pengmin Yang
#' @examples
#' data("aa")
#' cox = cox_forest(data=aa,
#'                           time_col = "OS_MONTHS",
#'                           status_col = "OS_STATUS",
#'                          Univariate=T,
#'                          univar_predictors=colnames(aa)[c(6:7,18:21,30,24,25,33)],
#'                          Multivariate=T,
#'                          multivar_predictors = colnames(aa)[c(6:7,18:21,30,24,25,33)],
#'                         show_plots = T,xticks1=NULL,#c(0,0.25,0.5,0.75,1.00,1.25,1.5,6.5,11),
#'                         xticks2=NULL,#c(0,0.25,0.5,0.75,1.00,2,2.5,6,15),
#'                         title_univar = "OS Univariate",
#'                         title_multivar = "OS Multivariate",
#'                         use_baseline_table = TRUE,all =T,forestplot=F,ci_pch=16,ci_col="red",ci_line="blue",zero_col="#e22e2a",
#'                         log2=T,footnote=paste("\nHRD1: without HRD value adjusted", "HRD5: with HRD value adjusted (median)", "HRD adjusted = LST-15.5*ploidy+LOH+TAI ",sep = "\n"))
#'
#' print(cox[["univariate"]][["uniforest"]])
#' grid::grid.newpage()
#' print(cox[["multivariate"]][["multiforest"]])
#'
cox_forest <- function(data,
                        time_col = "OS_MONTHS",
                        status_col = "OS_STATUS",
                        Univariate=T,
                        univar_predictors=colnames(clin)[c(3:6,18:22)],
                        Multivariate=T,
                        multivar_predictors = colnames(clin)[c(3:4,6,18:22)],
                        show_plots = TRUE,xticks1,xticks2,
                        title_univar = "OS Univariate",
                        title_multivar = "OS Multivariate",
                        use_baseline_table = TRUE,
                        all=T,forestplot=T,
                        ci_pch=16,ci_col="red",ci_line="blue",zero_col = "#e22e2a",
                        log2=FALSE,footnote=NULL) {
  require(survival)
  require(forestplot)
  require(forestploter)
  require(tableone)
  require(plyr)
  require(stringr)
  require(tibble)

  # 1) 数据预处理
  aa <- data
  #统一将本是数值型变量却是character的变量变成numeric
  is_all_numeric_chr <- function(x) {
    !is.factor(x) && is.character(x) && all(!is.na(suppressWarnings(as.numeric(x[complete.cases(x)]))))
  }
  cols_to_fix <- sapply(aa, is_all_numeric_chr)
  aa[cols_to_fix] <- lapply(aa[cols_to_fix], function(x) as.numeric(x))
  #aa[[time_col]]<-as.numeric(aa[[time_col]])
  #aa[[status_col]]<-as.numeric(aa[[status_col]])
  for (col in unique(c(univar_predictors, multivar_predictors))) {
    if (!is.numeric(aa[[col]])) {
      if (anyNA(aa[[col]])) {
        aa[[col]]<-as.factor(aa[[col]])
        levels(aa[[col]]) <- c(levels(aa[[col]]), "missing")
        # 统一把 NA 替换为 "missing"
        aa[[col]][is.na(aa[[col]])] <- "missing"
      }
    }else{
      aa[[col]]<-as.numeric(aa[[col]])
    }
  }

  # Surv 对象
  if (!(time_col %in% colnames(aa)) || !(status_col %in% colnames(aa))) {
    stop("请确保数据中存在 time_col 和 status_col 指定的列。")
  }

  if(Univariate){
    y <- Surv(time = aa[[time_col]], event = aa[[status_col]])
    # 2) 单因素分析：Uni_cox_model(x)
    Uni_cox_model<- function(x){
      FML <- as.formula(paste0 ("y~",sprintf("`%s`", x)))
      cox<- coxph(FML,data=aa)
      cox1<-summary(cox)
      HR <- round(cox1$coefficients[,2],2)
      PValue <- round(cox1$coefficients[,5],3)
      CI5 <-round(cox1$conf.int[,3],2)
      CI95 <-round(cox1$conf.int[,4],2)
      Uni_cox_model<- data.frame('Characteristics' = x,
                                 'HR' = HR,
                                 'CI5' = CI5,
                                 'CI95' = CI95,
                                 'P' = PValue)
      return(Uni_cox_model)}
    variable.names <- univar_predictors
    # 过滤掉不存在的变量名
    variable.names <- variable.names[variable.names %in% colnames(aa)]
    Uni_cox<- lapply(variable.names, Uni_cox_model)
    Uni_cox<- ldply(Uni_cox,data.frame);Uni_cox
    if(log2){
      Uni_cox$HR<-round(log2(Uni_cox$HR),3)
      Uni_cox$CI5<-round(log2(Uni_cox$CI5),3)
      Uni_cox$CI95<-round(log2(Uni_cox$CI95),3)
    }
    if (!is.null(Uni_cox) && nrow(Uni_cox) > 0) {
      if(log2){
        Uni_cox$HR.CI95 <- paste0(Uni_cox$HR, " (", Uni_cox$CI5, " ~ ", Uni_cox$CI95, ")")
      }else {
        Uni_cox$HR.CI95 <- paste0(Uni_cox$HR, " (", Uni_cox$CI5, " - ", Uni_cox$CI95, ")")
      }
      Uni_cox$P[Uni_cox$P == 0] <- "<0.001"
    }
    Uni_cox<-rbind(c("Characteristics",NA,NA,NA,"P value","HR (95%CI)"),
                   Uni_cox[1:nrow(Uni_cox),])
    factor_col<-c()
    for(col in univar_predictors){
      if (!is.numeric(aa[[col]])) {
        factor_col<-c(factor_col,col)
      }
    }
    for (col in factor_col) {
      #pos <- which(col %in% Uni_cox$Characteristics)[1]
      pos <- which(grepl(col, Uni_cox$Characteristics))[1]
      Uni_cox <- Uni_cox %>%
        add_row(Characteristics = col , .before = pos)
      Uni_cox[pos,"HR.CI95"]<-"Reference"
      #pos1 <- which(col %in% Uni_cox$Characteristics)[1]
      pos1 <- which(grepl(col, Uni_cox$Characteristics))[1]
      Uni_cox <- Uni_cox %>%
        add_row(Characteristics = col, .before = pos1)
    }

    uBaselineTable <- NULL
    if (use_baseline_table) {
      # 以 predictor_vars 作为变量做基线表
      myVars <- variable.names[variable.names %in% colnames(aa)]
      catVars<-c()
      for (col in myVars) {
        if (!is.numeric(aa[[col]])) {
          catVars<-c(catVars,col)
        }
      }
      if (length(myVars) > 0) {
        tbl <- print(CreateTableOne(vars=myVars,
                                    data = aa,
                                    factorVars = catVars),
                     showAllLevels=TRUE)
        uBaselineTable <-as.data.frame(tbl)
      }
    }
    #极线表处理
    for (col in factor_col) {
      pos <- which(grepl(col, rownames(uBaselineTable)))[1]
      uBaselineTable <- uBaselineTable %>%
        add_row(level =col , .before = pos)
      if(all){
        uBaselineTable[uBaselineTable$level==col,2]<-uBaselineTable[1,2]
      }else{
        uBaselineTable[uBaselineTable$level==col,2]<-length(aa[[col]][aa[[col]]!="missing"])
      }
    }
    uBaselineTable<-rbind(c("Characteristics","Number(%)"),
                          uBaselineTable[1:nrow(uBaselineTable),])
    empty_idx <- which(uBaselineTable$level == "")
    empty_idx <- empty_idx[-1]
    numeric_col<-c()
    for(col in univar_predictors){
      if (is.numeric(aa[[col]])) {
        numeric_col<-c(numeric_col,col)
      }
    }
    uBaselineTable$level[empty_idx] <- numeric_col
    rownames(uBaselineTable)<-NULL
    colnames(uBaselineTable)<-uBaselineTable[1,]
    uBaselineTable<-uBaselineTable[-2,]
    if(all){
      uBaselineTable<-uBaselineTable
      Uni_cox<-Uni_cox
      Uni_cox$Characteristics<-uBaselineTable$Characteristics
    }else{
      Uni_cox$Characteristics<-uBaselineTable$Characteristics
      uBaselineTable<-uBaselineTable[uBaselineTable$Characteristics!="missing",]
      Uni_cox<-Uni_cox[Uni_cox$Characteristics!="missing",]
    }
    #summary<-grepl(paste(factor_col, collapse = "|"),Uni_cox$Characteristics)
    summary<- Uni_cox$Characteristics %in% factor_col
    summary[[1]]<-TRUE
    #合并
    Uni_cox<-cbind(Uni_cox,uBaselineTable[,2])
    colnames(Uni_cox)[[7]]<-"Number(%)"
    Uni_cox<-Uni_cox[,c(1,7,2:6)]
    if(all){
      Uni_cox<-Uni_cox
    }else{
      Uni_cox[nrow(Uni_cox)+1,1]<-"Total"
      Uni_cox[nrow(Uni_cox),2]<-nrow(aa)
    }

    for(i in 3:5) {Uni_cox[, i] = as.numeric(Uni_cox[, i])}

    #summary1<-grepl(paste(univar_predictors, collapse = "|"),Uni_cox$Characteristics)
    summary1<- Uni_cox$Characteristics %in% univar_predictors
    summary1[[1]]<-TRUE
    Uni_cox$Characteristics[!summary1] <- paste0("  ", Uni_cox$Characteristics[!summary1])
    if (all){
      last_row <- nrow(Uni_cox)+1
    }else{
      last_row <- nrow(Uni_cox)
    }
    hrzl_lines <- list("1" = gpar(lty = 1, lwd = 2),
                       "2" = gpar(lty = 2))
    if (last_row >= 3) {
      hrzl_lines[[as.character(last_row)]] <- gpar(lwd = 2, lty = 1)
    }

    if(is.null(xticks1)){
      xticks<-c(Uni_cox[,4],Uni_cox[,5])
      xticks <- na.omit(xticks)

      # 2) 去掉离群值（IQR 法）
      qnt <- quantile(xticks, probs = c(.25, .75), na.rm = TRUE)
      iqr <- qnt[2] - qnt[1]
      lower <- qnt[1] - 1.5 * iqr
      upper <- qnt[2] + 1.5 * iqr
      xticks_clean <- xticks[xticks >= lower & xticks <= upper]

      # 3) 均分为 6 段（包括 0 与最大值）
      min_val <- if (log2) floor(min(xticks_clean)) else 0 #向下取整，确保整数
      max_val <- ceiling(max(xticks_clean))  # 向上取整，确保整数
      xticks1 <- seq(min_val, max_val, length.out = 7)
    }else{
      xticks1<-xticks1
    }

    rows1  <- which(Uni_cox$P < 0.05 & Uni_cox$P >= 0.01)
    Uni_cox[rows1,"P"]<-paste0(Uni_cox[rows1,"P"],"*")
    rows2  <- which(Uni_cox$P < 0.01 & Uni_cox$P >= 0.001)
    Uni_cox[rows2,"P"]<-paste0(Uni_cox[rows2,"P"],"**")
    rows3  <- which(Uni_cox$P =="<0.001")
    Uni_cox[rows3,"P"]<-paste0(Uni_cox[rows3,"P"],"***")

    col<-rep("white", nrow(Uni_cox))
    col[c(rows1, rows2, rows3)] <- "#f97"


    if (show_plots) {
      if (!is.null(Uni_cox) && nrow(Uni_cox) > 0) {
        if(forestplot){
          library(forestplot)
          Uni_cox$Characteristics[Uni_cox$Characteristics%in%univar_predictors]<-gsub("_|-"," ",Uni_cox$Characteristics[Uni_cox$Characteristics%in%univar_predictors])
          p<-forestplot( Uni_cox[,c(1,2,7,6)],
                         labeltext = as.character(Uni_cox$Characteristics),
                         mean = Uni_cox[,3],
                         lower = Uni_cox[,4],
                         upper = Uni_cox[,5],
                         zero = if (log2) 0 else 1,
                         boxsize = 0.4,
                         graph.pos = "right",
                         xlab=if (log2) "log2(Hazard Ratio)" else "HR",
                         hrzl_lines=hrzl_lines,
                         title = title_univar,
                         xticks=if (log2) seq(min(xticks1), ceiling(max(xticks1)) - max(0, floor(min(xticks1))), by = 1) else
                           seq(0, ceiling(max(xticks1)) - max(0, floor(min(xticks1))), by = 1),
                         is.summary = summary,
                         txt_gp = fpTxtGp(label = gpar(cex = 1),
                                          ticks = gpar(cex = 1.2),
                                          xlab = gpar(cex = 1.3),
                                          title = gpar(cex = 1.5)),
                         lwd.zero=1,
                         lwd.ci=1.6,
                         lwd.xaxis=1,
                         lty.ci=1,
                         ci.vertices =T,
                         ci.vertices.height=0.2,
                         clip=c(0.1,8),
                         #----------------#行间距、字间距/box形状
                         ineheight=unit(8, 'mm'),
                         line.margin=unit(8, 'mm'),
                         colgap=unit(6, 'mm'),
                         col=fpColors(zero = zero_col,
                                      box = ci_col,
                                      lines = ci_line),
                         fn.ci_norm="fpDrawCircleCI")|>
            fp_decorate_graph(grid=gpar(lty = 2, col = "black"),
                              graph.pos = 4)
          dev.off()
        }else{
          library(forestploter)
          dt<-Uni_cox
          dt$Characteristics[dt$Characteristics%in%multivar_predictors]<-gsub("_|-"," ",dt$Characteristics[dt$Characteristics%in%multivar_predictors])
          dt$` ` <- paste(rep(" ", 20), collapse = " ")
          colnames(dt)[c(1,2,6,7,8)]<-dt[1,c(1,2,6,7,8)]
          dt<-dt[-1,]
          colnames(dt)[[7]]<-"log2(HR(95%CI))"
          for (j in c(1,2,6,7,8)) {
            dt[[j]] <- ifelse(is.na(dt[[j]]), " ", dt[[j]])
          }
          rows <- which(grepl("*", dt$`P value`, fixed = TRUE))   # 找到带 * 的行号
          dt$HR[dt$HR == 0] <- 1e-5
          tm <- forest_theme(base_size = 10,
                             # Confidence interval point shape, line type/color/width
                             ci_pch = ci_pch,
                             ci_col = "grey20",
                             ci_fill = "grey20",
                             ci_alpha = 0.8,
                             ci_lty = 1,
                             ci_lwd = 1,
                             ci_Theight = 0.2, # Set a T end at the end of CI
                             # Reference line width/type/color
                             refline_gp = gpar(lwd = 1, lty = "dashed", col = zero_col),
                             # Vertical line width/type/color
                             vertline_lwd = 1,
                             vertline_lty = "dashed",
                             vertline_col = "grey20",
                             # Change summary color for filling and borders
                             summary_fill = "#4575b4",
                             summary_col = "#4575b4",
                             title_just = "left",
                             title_gp = gpar(cex = 1.2, fontface = "bold", col = "black"),
                             # Footnote font size/face/color
                             footnote_gp = gpar(cex = 0.7, fontface = "italic", col = "blue"),
                             core = list(bg_params = list(fill = "white")),
                             colhead = list(fg_params = list(hjust = 0.5, x = 0.5)),
                             arrow_gp = gpar(col = "blue"))

          p<-forest(dt[,c(1,2,7,8,6)],
                    est = dt$HR,       #效应值
                    lower = dt$CI5,     #可信区间下限
                    upper = dt$CI95,      #可信区间上限
                    sizes = if (log2) abs((dt$HR)/2) else (dt$HR)/2,     #黑框的大小
                    ci_column = 4,   #在那一列画森林图，要选空的那一列
                    ref_line = if (log2) 0 else 1,
                    title = if (is.null(footnote)) NULL else title_univar,
                    arrow_lab = c("Worser", "Better"),
                    xlim = if (log2) c(min(xticks1),ceiling(max(xticks1))) else
                      c(if (min(xticks1) < 1) 0 else min(xticks1),
                        ceiling(max(xticks1))),
                    ticks_at = if (log2) seq(min(xticks1), ceiling(max(xticks1)) - max(0, floor(min(xticks1))), by = 1) else
                      seq(0, ceiling(max(xticks1)) - max(0, floor(min(xticks1))), by = 1),
                    footnote = if (is.null(footnote)) paste0("\n",title_univar) else paste0("\n",footnote),
                    theme = tm)
          p<-add_border(p,part = "header",
                        row = c(0),
                        gp = gpar(lty = c(1),lwd=c(2)))
          p<-add_border(p,part = "header",
                        row = c(1),
                        gp = gpar(lty = c(2),lwd=c(1)))
          n=if(all){nrow(dt)+1} else {nrow(dt)}
          p<-add_border(p,part = "header",
                        row = c(n),
                        gp = gpar(lty = c(1),lwd=c(2)))
          p<-edit_plot(p, row = rows, gp = gpar(col = "red", fontface = "italic"))
          factor_col<-gsub("-|_"," ",factor_col)
          summary<-which(dt$Characteristics %in% factor_col)
          #summary<-which(grepl(paste(factor_col, collapse = "|"),dt$Characteristics))
          p<-edit_plot(p, row = summary, gp = gpar(col = "black", fontface = "bold"))
          if(all){
            p<-p
          }else{
            p<-edit_plot(p, row = nrow(dt), gp = gpar(col = "black", fontface = "bold"))
          }
          p<-edit_plot(p,
                       row = rows,
                       col= 4,
                       which = "ci",
                       gp=gpar(col=ci_line,fill =ci_col,
                               lwd=2))
        }
      }}
  }else{
    Uni_cox<-NULL
    p<-NULL
    uBaselineTable<-NULL
  }


  # 3) 多因素分析（若提供）: Multivariate
  if(Multivariate){
    Mult_cox <- NULL
    if (!is.null(multivar_predictors) && length(multivar_predictors) > 0) {
      mv_preds <- multivar_predictors[multivar_predictors %in% colnames(aa)]
      if (length(mv_preds) > 0) {
        mv_form <- as.formula(paste0("Surv(", time_col, ",", status_col, ") ~ ",
                                     paste(sprintf("`%s`", mv_preds), collapse = " + "))
        )
        mv_cox <- coxph(mv_form, data = aa)
        mv_sum <- summary(mv_cox)
        mv_table <- data.frame(round(mv_sum$conf.int, 3))
        mv_table$p <- round(mv_sum$coefficients[, 5], 3)
        mv_table<-mv_table[,-2]
        mv_table<-tibble::rownames_to_column(mv_table, var ="Characteristics")
        colnames(mv_table)<-c("Characteristics","HR","CI5","CI95","p")
        Mult_cox <- mv_table
      }
    }
    if(log2){
      Mult_cox$HR<-round(log2(Mult_cox$HR),3)
      Mult_cox$CI5<-round(log2(Mult_cox$CI5),3)
      Mult_cox$CI95<-round(log2(Mult_cox$CI95),3)
    }
    if (!is.null(Mult_cox ) && nrow(Mult_cox ) > 0) {
      if(log2){
        Mult_cox$HR.CI95 <- paste0(Mult_cox$HR, " (", Mult_cox$CI5, " ~ ", Mult_cox$CI95, ")")
      }else {
        Mult_cox$HR.CI95 <- paste0(Mult_cox$HR, " (", Mult_cox$CI5, " - ", Mult_cox$CI95, ")")
      }
      Mult_cox$p[Mult_cox $p == 0] <- "<0.001"
    }

    Mult_cox<-rbind(c("Characteristics",NA,NA,NA,"P value","HR (95%CI)"),
                    Mult_cox[1:nrow(Mult_cox),])
    factor_col<-c()
    for(col in multivar_predictors){
      if (!is.numeric(aa[[col]])) {
        factor_col<-c(factor_col,col)
      }
    }
    for (col in factor_col) {
      aa_col <- unique(aa[[col]])
      aa_col[is.na(aa_col)] <- "missing"
      Mult_cox_col <- str_remove(Mult_cox$Characteristics[grepl(col, Mult_cox$Characteristics)], col)
      t_col <- aa_col[!aa_col %in% Mult_cox_col]
      if (length(t_col) == 0) next   # 没有缺失值就跳过
      # 找到第一个出现 col 的行号
      #pos <-which(col %in% Mult_cox$Characteristics)[1]
      pos <- which(grepl(col, Mult_cox$Characteristics))[1]
      Mult_cox <- Mult_cox %>%
        add_row(Characteristics =t_col , .before = pos) %>%
        add_row(Characteristics = col, .before = pos)
      Mult_cox[Mult_cox$Characteristics == t_col,"HR.CI95"]<-"Reference"
    }
    str<-paste(factor_col,collapse ="|")
    Mult_cox$Characteristics<-str_remove(Mult_cox$Characteristics,str)
    Mult_cox$Characteristics[Mult_cox$Characteristics==""]<-factor_col

    # 4) 基线表（可选，使用 tableone）
    mBaselineTable <- NULL
    if (use_baseline_table) {
      # 以 predictor_vars 作为变量做基线表
      myVars <- multivar_predictors[multivar_predictors %in% colnames(aa)]
      catVars<-c()
      for (col in myVars) {
        if (!is.numeric(aa[[col]])) {
          catVars<-c(catVars,col)
        }
      }
      if (length(myVars) > 0) {
        tbl <- print(CreateTableOne(vars=myVars,
                                    data = aa,
                                    factorVars = catVars),
                     showAllLevels=TRUE)
        mBaselineTable <-as.data.frame(tbl)
      }
    }
    #基线表处理
    empty_idx <- which(mBaselineTable$level == "")
    empty_idx <- empty_idx[-1]
    numeric_col<-c()
    for(col in multivar_predictors){
      if (is.numeric(aa[[col]])) {
        numeric_col<-c(numeric_col,col)
      }
    }
    mBaselineTable$level[empty_idx] <- numeric_col
    for (col in factor_col) {
      pos <- which(grepl(col, rownames(mBaselineTable)))[1]
      mBaselineTable <- mBaselineTable %>%
        add_row(level =col , .before = pos)
      if(all){
        mBaselineTable[mBaselineTable$level==col,2]<-mBaselineTable[1,2]
      }else{
        mBaselineTable[mBaselineTable$level==col,2]<-length(aa[[col]][aa[[col]]!="missing"])
        mBaselineTable<-mBaselineTable[mBaselineTable$level!="missing",]
        Mult_cox<-Mult_cox[Mult_cox$Characteristics!="missing",]
      }
    }
    #summary<-grepl(paste(factor_col, collapse = "|"),Mult_cox$Characteristics)
    summary<- Mult_cox$Characteristics %in% factor_col
    summary[[1]]<-TRUE

    mBaselineTable<-rbind(c("Characteristics","Number(%)"),
                          mBaselineTable[1:nrow(mBaselineTable),])
    rownames(mBaselineTable)<-NULL
    colnames(mBaselineTable)<-mBaselineTable[1,]
    mBaselineTable<-mBaselineTable[-2,]

    #mBaselineTable$Characteristics<-Mult_cox$Characteristics
    Mult_cox$Characteristics<-mBaselineTable$Characteristics
    #合并
    Mult_cox<-cbind(Mult_cox,mBaselineTable[,2])
    colnames(Mult_cox)[[7]]<-"Number(%)"
    Mult_cox<-Mult_cox[,c(1,7,2:6)]

    if(all){
      Mult_cox<-Mult_cox
    }else{
      Mult_cox[nrow(Mult_cox)+1,1]<-"Total"
      Mult_cox[nrow(Mult_cox),2]<-nrow(aa)
    }
    for(i in 3:5) {Mult_cox[, i] = as.numeric(Mult_cox[, i])}


    #summary1<-grepl(paste(multivar_predictors, collapse = "|"),Mult_cox$Characteristics)
    summary1<-Mult_cox$Characteristics %in% multivar_predictors
    summary1[[1]]<-TRUE
    Mult_cox$Characteristics[!summary1] <- paste0("  ", Mult_cox$Characteristics[!summary1])

    if (all){
      last_row <- nrow(Mult_cox)+1
    }else{
      last_row <- nrow(Mult_cox)
    }
    hrzl_lines <- list("1" = gpar(lty = 1, lwd = 2),
                       "2" = gpar(lty = 2))
    if (last_row >= 3) {
      hrzl_lines[[as.character(last_row)]] <- gpar(lwd = 2, lty = 1)
    }

    if(is.null(xticks2)){
      xticks<-c(Mult_cox[,4],Mult_cox[,5])
      xticks <- na.omit(xticks)

      # 2) 去掉离群值（IQR 法）
      qnt <- quantile(xticks, probs = c(.25, .75), na.rm = TRUE)
      iqr <- qnt[2] - qnt[1]
      lower <- qnt[1] - 1.5 * iqr
      upper <- qnt[2] + 1.5 * iqr
      xticks_clean <- xticks[xticks >= lower & xticks <= upper]

      # 3) 均分为 6 段（包括 0 与最大值）
      min_val <- if (log2) floor(min(xticks_clean)) else 0 #向下取整，确保整数
      max_val <- ceiling(max(xticks_clean))  # 向上取整，确保整数
      xticks2 <- seq(min_val, max_val, length.out = 7)
    }else{
      xticks2<-xticks2
    }

    rows1  <- which(Mult_cox$p < 0.05 & Mult_cox$p >= 0.01)
    Mult_cox[rows1,"p"]<-paste0(Mult_cox[rows1,"p"],"*")
    rows2  <- which(Mult_cox$p < 0.01 & Mult_cox$p >= 0.001)
    Mult_cox[rows2,"p"]<-paste0(Mult_cox[rows2,"p"],"**")
    rows3  <- which(Mult_cox$p =="<0.001")
    Mult_cox[rows3,"p"]<-paste0(Mult_cox[rows3,"p"],"***")

    if (show_plots) {
      if (!is.null(Mult_cox) && nrow(Mult_cox) > 0) {
        #sum_value1<- as.character(length(rownames(cox$multivariate)))
        if(forestplot){
          library(forestplot)
          Mult_cox$Characteristics[Mult_cox$Characteristics%in%multivar_predictors]<-gsub("_|-"," ",Mult_cox$Characteristics[Mult_cox$Characteristics%in%multivar_predictors])
          p1<-forestplot(Mult_cox[,c(1,2,7,6)],
                         labeltext = as.character(Mult_cox$Characteristics),
                         mean = Mult_cox[,3],
                         lower = Mult_cox[,4],
                         upper = Mult_cox[,5],
                         zero = if (log2) 0 else 1,
                         boxsize = 0.4,
                         graph.pos = "right",
                         xlab = if (log2) "log2(Hazard Ratio)" else "HR",
                         hrzl_lines=hrzl_lines,
                         title = title_multivar,
                         xticks=if (log2) seq(min(xticks2), ceiling(max(xticks2)) - max(0, floor(min(xticks2))), by = 1) else
                           seq(0, ceiling(max(xticks2)) - max(0, floor(min(xticks2))), by = 1),
                         is.summary =summary, #rep(FALSE, nrow(Mult_cox)),
                         txt_gp = fpTxtGp(label = gpar(cex = 1),
                                          ticks = gpar(cex = 1.2),
                                          xlab = gpar(cex = 1.3),
                                          title = gpar(cex = 1.5)),
                         lwd.zero=1,
                         lwd.ci=1.6,
                         lwd.xaxis=1,
                         lty.ci=1,
                         ci.vertices =T,
                         ci.vertices.height=0.2,
                         clip=c(0.1,8),
                         #----------------#行间距、字间距/box形状
                         ineheight=unit(8, 'mm'),
                         line.margin=unit(8, 'mm'),
                         colgap=unit(6, 'mm'),
                         col=fpColors(zero = zero_col,
                                      box = ci_col,
                                      lines = ci_line),
                         fn.ci_norm="fpDrawCircleCI")|>
            fp_decorate_graph(grid=gpar(lty = 2, col = "black"),
                              graph.pos = 4)
          dev.off()
        }else{
          library(forestploter)
          dt<-Mult_cox
          dt$Characteristics[dt$Characteristics%in%multivar_predictors]<-gsub("_|-"," ",dt$Characteristics[dt$Characteristics%in%multivar_predictors])
          dt$` ` <- paste(rep(" ", 20), collapse = " ")
          colnames(dt)[c(1,2,6,7,8)]<-dt[1,c(1,2,6,7,8)]
          dt<-dt[-1,]
          colnames(dt)[[7]]<-"log2(HR(95%CI))"
          for (j in c(1,2,6,7,8)) {
            dt[[j]] <- ifelse(is.na(dt[[j]]), " ", dt[[j]])
          }

          rows <- which(grepl("*", dt$`P value`, fixed = TRUE))   # 找到带 * 的行号
          dt$HR[dt$HR == 0] <- 1e-5
          xlim <-if (log2) c(min(xticks2),ceiling(max(xticks2))) else
            c(if (min(xticks2) < 1) 0 else min(xticks2),
              ceiling(max(xticks2)))
          ticks_at <- if (log2) seq(min(xticks2), ceiling(max(xticks2)) - max(0, floor(min(xticks2))), by = 1) else
            seq(0, ceiling(max(xticks2)) - max(0, floor(min(xticks2))), by = 1)

          tm <- forest_theme(base_size = 10,
                             # Confidence interval point shape, line type/color/width
                             ci_pch = ci_pch,
                             ci_col = "grey20",
                             ci_fill = "grey20",
                             ci_alpha = 0.8,
                             ci_lty = 1,
                             ci_lwd = 1,
                             ci_Theight = 0.2, # Set a T end at the end of CI
                             # Reference line width/type/color
                             refline_gp = gpar(lwd = 1, lty = 1, col =  zero_col),
                             # Vertical line width/type/color
                             vertline_lwd = 1,
                             vertline_lty = 2,
                             vertline_col = "grey20",
                             # Change summary color for filling and borders
                             summary_fill = "#4575b4",
                             summary_col = "#4575b4",
                             title_just = "left",
                             title_gp = gpar(cex = 1.2, fontface = "bold", col = "black"),
                             # Footnote font size/face/color
                             footnote_gp = gpar(cex = 0.7, fontface = "italic", col = "blue"),
                             core = list(bg_params = list(fill = "white")),
                             colhead = list(fg_params = list(hjust = 0.5, x = 0.5)),
                             arrow_gp = gpar(col = "blue"))

          p1<-forest(dt[,c(1,2,7,8,6)],
                     est = dt$HR,       #效应值
                     lower = dt$CI5,     #可信区间下限
                     upper = dt$CI95,      #可信区间上限
                     sizes = if (log2) abs(dt$HR/2) else (dt$HR)/2,     #黑框的大小
                     ci_column = 4,   #在那一列画森林图，要选空的那一列
                     ref_line = if (log2) 0 else 1,
                     title = if (is.null(footnote))NULL else title_multivar,
                     is_summary = c(rep(FALSE, nrow(dt)-1), TRUE),
                     arrow_lab = c("Worser", "Better"),
                     xlim =xlim,
                     ticks_at = ticks_at,
                     footnote = if (is.null(footnote)) paste0("\n",title_multivar) else paste0("\n",footnote),
                     theme = tm)
          p1<-add_border(p1,part = "header",
                         row = c(0),
                         gp = gpar(lty = c(1),lwd=c(2)))
          p1<-add_border(p1,part = "header",
                         row = c(1),
                         gp = gpar(lty = c(2),lwd=c(1)))
          n=if(all){nrow(dt)+1} else {nrow(dt)}
          p1<-add_border(p1,part = "header",
                         row = c(n),
                         gp = gpar(lty = c(1),lwd=c(2)))
          p1<-edit_plot(p1, row = rows, gp = gpar(col = "red", fontface = "italic"))
          factor_col<-gsub("-|_"," ",factor_col)
          summary<-which(dt$Characteristics%in% factor_col)
          #summary<-which(grepl(paste(factor_col, collapse = "|"),dt$Characteristics))
          p1<-edit_plot(p1, row = summary, gp = gpar(col = "black", fontface = "bold"))
          if(all){
            p1<-p1
          }else{
            p1<-edit_plot(p1, row = nrow(dt), gp = gpar(col = "black", fontface = "bold"))
          }

          p1<-edit_plot(p1,
                        row = rows,
                        col= 4,
                        which = "ci",
                        gp=gpar(col=ci_line,fill =ci_col,
                                lwd=2))
        }
      }
    }
  }else{
    Mult_cox<-NULL
    p1<-NULL
    mBaselineTable<-NULL
  }
  # 6) 返回结果
  res <- list(univariate=list(univariate = Uni_cox,
                              uniforest=p,
                              baseline_table = uBaselineTable),
              multivariate=list(multivariate = Mult_cox,
                                multiforest=p1,
                                baseline_table = mBaselineTable),
              data_used = aa)
  return(res)
}


