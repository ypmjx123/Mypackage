#' Differential Gene Expression Data Visualization
#'
#' This function visualizes the Differential Gene Expression Data.
#'
#' @param exprdata    is the data of gene expression for samples
#' @param pdata       A data frame that contains groups of genes.例如data.frame(SMPLE_ID,Group)
#' @param Type        A character vector indicating the type.eg Type=c("Wild", "Mut")
#' @param tumor       the source of exprdata
#' @param diff_method 基因表达差异分析方法c("limma","DESeq2","edgeR")
#' @param contrastfml the Type for limma; eg contrastfml="Wild - Mut"
#' @param datatype    exprdata的类型是"counts"时可以用"DESeq2"/"edgeR"方法;"limma"不挑
#' @param P.Value     p阈值for差异分析
#' @param logFC       logFC阈值for差异分析
#' @param ann_colors  ComplexHeatmap::pheatmap 参数for sample annotation
#' @param color       ComplexHeatmap::pheatmap 参数gene expression
#' @param Regulate    select Regulate=c("Up","Down") for pathway enrichment
#' @param geom_text   defalut KEGG=FALSE;if TURE,Volcano plot will annotate with
#' @param tidyHeatmap defalut tidyHeatmap=TRUE
#' @param heatmap_col defalut heatmap_col=NULL,tidyHeatmap::heatmap参数palette_value Regulate
#' @param top50       defalut top50=TRUE 热图展示所用
#' @param GO          defalut GO=TRUE
#' @param GO.plot     GO富集图展现形式，GO.plot=c("dotplot","barplot","heatplot","cnetplot","gseaplot")
#' @param split       defalut split=FALSE;if TURE,GO with "BP", "MF" and "CC",else GO with "BP"
#' @param KEGG        defalut KEGG=TRUE
#' @param KEGG.plot   KEGG富集图展现形式，GO.plot=c("dotplot","barplot","heatplot","cnetplot","gseaplot")
#' @param rel_heights relative heights of subplots. defalut rel_heights = c(1.5, 0.5, 1)
#'
#' @return A list of differential analysis results,volcano,heatmap and enrichment plot of GO AND KEGG
#' @export
##' @author Pengmin Yang
#' @examples
#' data("exp_raw")
#' data("Gene_group_CRC1")
#' exp_CRC<-exp_geneIDtoSYMBOL(exp=exp_raw,genecoltype="ENTREZID")
#' exp_CRC<-exp_CRC$data
#' result<-limma.dif.visual(exprdata=exp_CRC[,-1],
#'                            pdata=Gene_group_CRC1,datatype="TPM",
#'                            Type=c("Wild", "Mut"),diff_method="limma",
#'                            contrastfml= "Wild - Mut",
#'                            tumor="CRC TCGA",
#'                            P.Value=0.05,
#'                            logFC=0.5,tidyHeatmap=TRUE,
#'                            color= NULL,
#'                            ann_colors = list(regulate = c(Down = "#1B9E77", Up = "#D95F02"),
#'                            PREX2 = c(Wild = "#757575", Mut = "#FF4040")),
#'                            Regulate=c("Up","Down"),GO=TRUE,GO.plot="dotplot",split=TRUE,
#'                            KEGG=TRUE,KEGG.plot="dotplot",rel_heights= c(1.5, 0.5, 1))
#' print(result)

limma.dif.visual<-function (exprdata,datatype="counts", pdata,Type,
                            tumor,diff_method=c("limma","DESeq2","edgeR"),
                            contrastfml,P.Value,logFC,ann_colors,color,Regulate,geom_text=FALSE,
                            tidyHeatmap=TRUE,heatmap_col=NULL,top50=TRUE,
                            GO=TRUE,GO.plot="dotplot",split=FALSE,
                            KEGG=TRUE,KEGG.plot="dotplot",rel_heights= c(1.5, 0.5, 1))
{
  library(limma)
  library(DESeq2)
  library(edgeR)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(dplyr)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(base)
  exprdata<-exprdata[!is.na(exprdata[,2]),]
  if(!is.na(as.numeric(rownames(exprdata)[[length(rownames(exprdata))]]))){
    if (any(duplicated(exprdata[[1]]))) {
      print(" expression data genes duplicated ")
      exprdata <- aggregate(as.formula(paste0(".", " ~ ", colnames(exprdata)[[1]])), data = exprdata, FUN = mean)
    }
    rownames(exprdata)<-exprdata[[1]]
    exprdata<-exprdata[,-1]
  }
  if (!all(colnames(exprdata) == pdata[, 1])) {
    print(" expression data do not match pdata, we conduct it")
    exprdata<- exprdata[, colnames(exprdata) %in% pdata[, 1], drop = FALSE]
    pdata1<-data.frame(SAMPLE_ID=colnames(exprdata))
    pdata1<-merge(pdata,pdata1,by.x=colnames(pdata)[[1]],by.y=colnames(pdata1)[[1]],all.y=T)
    pdata<-pdata1
  }
  pdata<-pdata
  if(diff_method=="limma"){
    pdata<-unique(pdata)
    group_list <- as.character(pdata[, 2])
    design <- model.matrix(~0 + factor(group_list))
    colnames(design) <- levels(as.factor(pdata[, 2]))
    rownames(design) <- colnames(exprdata)
    contrast.matrix <- makeContrasts(contrasts = contrastfml,
                                     levels = design)
    fit <- lmFit(exprdata, design)
    fit <- contrasts.fit(fit, contrast.matrix)
    fit <- eBayes(fit)
    dif <- topTable(fit, adjust.method = "BH", coef = contrastfml,
                    number = Inf)
    dif$regulate <- ifelse(dif$P.Value > P.Value, "Unchanged",
                           ifelse(dif$logFC > logFC, "Up",
                                  ifelse(dif$logFC < -logFC, "Down", "Unchanged")))
    minlogFC<-min(dif$logFC[dif$regulate == "Down"], na.rm = TRUE)
    maxlogFC<-max(dif$logFC[dif$regulate == "Up"], na.rm = TRUE)
    dif<-dif[dif$logFC <= maxlogFC & dif$logFC >= minlogFC, ]
    Down<-paste("Down",paste0("(",length(dif$regulate[dif$regulate=="Down"])," ","genes",")"))
    Up<-paste("Up",paste0("(",length(dif$regulate[dif$regulate=="Up"])," ","genes",")"))
    Unchanged<-paste("Unchanged",paste0("(",length(dif$regulate[dif$regulate=="Unchanged"])," ","genes",")"))
    dif1<-dif
    dif1$gene<-rownames(dif1)
    dif1$regulate<-ifelse(dif$regulate=="Down",Down,ifelse(dif$regulate=="Up",Up,ifelse(dif$regulate=="Unchanged",Unchanged,dif$regulate)))
  }else if(diff_method=="DESeq2"){
    pdata[[colnames(pdata)[[2]]]] <- factor(pdata[[colnames(pdata)[[2]]]], levels = Type)
    target<-colnames(pdata)[[2]]
    colnames(pdata)[[2]]<-"group"
    group <- factor(pdata[,2], levels = Type)
    exprdata <- round(exprdata)
    dds <- DESeqDataSetFromMatrix(
      countData = exprdata, #表达矩阵
      colData = pdata, #表达矩阵列名和分组的对应关系
      design = ~ group) #实验设计，condition对应着colData里的condition，也就是分组信息
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("group", rev(levels(group))))
    resOrdered <- res[order(res$padj), ]
    DEG <- as.data.frame(resOrdered)
    DEG_deseq2 <- na.omit(DEG)
    dif<-DEG_deseq2
    #logFC <- with(dif,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)) )
    colnames(pdata)[[2]]<-target
    dif$regulate <- ifelse(dif$padj > P.Value, "Unchanged",
                           ifelse(dif$log2FoldChange > logFC, "Up",
                                  ifelse(dif$log2FoldChange < -logFC, "Down", "Unchanged")))
    minlogFC<-min(dif$log2FoldChange[dif$regulate == "Down"], na.rm = TRUE)
    maxlogFC<-max(dif$log2FoldChange[dif$regulate == "Up"], na.rm = TRUE)
    dif<-dif[dif$log2FoldChange <= maxlogFC & dif$log2FoldChange >= minlogFC, ]
    Down<-paste("Down",paste0("(",length(dif$regulate[dif$regulate=="Down"])," ","genes",")"))
    Up<-paste("Up",paste0("(",length(dif$regulate[dif$regulate=="Up"])," ","genes",")"))
    Unchanged<-paste("Unchanged",paste0("(",length(dif$regulate[dif$regulate=="Unchanged"])," ","genes",")"))
    dif1<-dif
    dif1$gene<-rownames(dif1)
    dif1$regulate<-ifelse(dif$regulate=="Down",Down,ifelse(dif$regulate=="Up",Up,ifelse(dif$regulate=="Unchanged",Unchanged,dif$regulate)))
  }else if(diff_method=="edgeR"){
    pdata[[colnames(pdata)[[2]]]] <- factor(pdata[[colnames(pdata)[[2]]]], levels = Type)
    #colnames(pdata)[[2]]<-"group"
    group <- factor(pdata[,2], levels = Type)
    exprdata <- round(exprdata)
    d <- DGEList(counts = exprdata, group = group)
    keep <- rowSums(cpm(d) > 1) >= 2
    d <- d[keep, , keep.lib.sizes = FALSE]
    # 更新样本的库大小信息
    d$samples$lib.size <- colSums(d$counts)
    d <- calcNormFactors(d)
    dge = d
    # 创建设计矩阵，用于指定差异分析模型
    design <- model.matrix(~0 + factor(group))
    rownames(design) <- colnames(dge)
    colnames(design) <- levels(factor(group))
    # 估计数据的离散度 —— common离散度、trended离散度、tagwise离散度
    dge <- estimateGLMCommonDisp(dge, design)
    dge <- estimateGLMTrendedDisp(dge, design)
    dge <- estimateGLMTagwiseDisp(dge, design)
    # 在估计的模型基础上进行 广义线性模型 (GLM) 拟合
    fit <- glmFit(dge, design)
    lrt <- glmLRT(fit, contrast = c(-1, 1))
    nrDEG <- topTags(lrt, n = nrow(dge))
    # 将差异表达基因结果转换为数据框形式
    DEG_edgeR <- as.data.frame(nrDEG)
    dif<-DEG_edgeR
    #logFC <- with(dif,mean(abs(logFC)) + 2*sd(abs(logFC)) )
    dif$regulate <- ifelse(dif$FDR > P.Value, "Unchanged",
                           ifelse(dif$logFC  > logFC, "Up",
                                  ifelse(dif$logFC  < -logFC, "Down", "Unchanged")))
    minlogFC<-min(dif$logFC [dif$regulate == "Down"], na.rm = TRUE)
    maxlogFC<-max(dif$logFC [dif$regulate == "Up"], na.rm = TRUE)
    dif<-dif[dif$logFC  <= maxlogFC & dif$logFC  >= minlogFC, ]
    Down<-paste("Down",paste0("(",length(dif$regulate[dif$regulate=="Down"])," ","genes",")"))
    Up<-paste("Up",paste0("(",length(dif$regulate[dif$regulate=="Up"])," ","genes",")"))
    Unchanged<-paste("Unchanged",paste0("(",length(dif$regulate[dif$regulate=="Unchanged"])," ","genes",")"))
    dif1<-dif
    dif1$gene<-rownames(dif1)
    dif1$regulate<-ifelse(dif$regulate=="Down",Down,ifelse(dif$regulate=="Up",Up,ifelse(dif$regulate=="Unchanged",Unchanged,dif$regulate)))
  }
  if(diff_method=="limma"){
    top_20 <- bind_rows(
      dif1 %>%
        filter(regulate == Up) %>%
        arrange(P.Value, desc(abs(logFC))) %>%
        head(20),
      dif1 %>%
        filter(regulate == Down) %>%
        arrange(P.Value, desc(abs(logFC))) %>%
        head(20)
    )
    Exression<-dif1$regulate
    p<-ggplot(dif1, aes(logFC, -log10(P.Value)))+
      geom_point(aes(col=Exression,shape = Exression))+
      scale_shape_manual(values = c(25, 21,24))+
      #scale_shape_manual(values=c(Down=25, Up=24, Unchanged=21))+
      scale_color_manual(values=c(ann_colors[[1]][[1]],"grey",ann_colors[[1]][[2]]))+
      labs(title = " ")+
      geom_vline(xintercept=c(-logFC,logFC), colour="black", linetype="dashed")+
      geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+
      labs(x="log2(FoldChange)",y="-log10(Pvalue)",title = paste("Volcano plot of",tumor))+
      theme(axis.text=element_text(size=13,color = "black"),axis.title=element_text(size=13),
            plot.title = element_text(hjust = 0,size=13,face="plain"),
            legend.position ="bottom",
            panel.border = element_rect(color = "black", fill = NA),
            panel.background = element_rect(fill = "transparent"))+
      str(dif1)+#theme_bw()+
      geom_label_repel(data = top_20,max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                       aes(logFC, -log10(P.Value), label = gene),
                       size = 3, fill="#CCFFFF")
    if(geom_text){
      if(logFC!=0){
        min<-min(dif$logFC[dif$regulate == "Down"], na.rm = TRUE)
        max<-max(dif$logFC[dif$regulate == "Up"], na.rm = TRUE)
        p<-p + geom_text(x = (-logFC+min)/2, y = 0, label = paste(length(dif$regulate[dif$regulate=="Down"]),"genes"),
                         color = ann_colors[[1]][[1]], size = 4.5, vjust = 0)
        p<-p + geom_text(x = 0, y = 5, label = paste(length(dif$regulate[dif$regulate=="Unchanged"]),"genes"),
                         color = "grey", size = 3, vjust = 0)
        p<-p + geom_text(x = (logFC+max)/2, y = -log10(0.05), label = paste(length(dif$regulate[dif$regulate=="Up"]),"genes"),
                         color = ann_colors[[1]][[2]], size = 4.5, vjust = 0)
      }else{
        min<-min(dif$logFC[dif$regulate == "Down"], na.rm = TRUE)
        max<-max(dif$logFC[dif$regulate == "Up"], na.rm = TRUE)
        Pmax<-max(-log10(dif$P.Value), na.rm = TRUE)
        p<-p + geom_text(x = min/2, y = Pmax/2, label = paste(length(dif$regulate[dif$regulate=="Down"]),"genes"),
                         color = ann_colors[[1]][[1]], size = 4.5, vjust = 0)
        p<-p + geom_text(x = 0, y = -0.4, label = paste(length(dif$regulate[dif$regulate=="Unchanged"]),"genes"),
                         color = "grey", size = 3, vjust = 0)
        p<-p + geom_text(x = max/2, y = Pmax/2, label = paste(length(dif$regulate[dif$regulate=="Up"]),"genes"),
                         color = ann_colors[[1]][[2]], size = 4.5, vjust = 0)
      }

      #print(p)
    }
  }else if(diff_method=="DESeq2"){
    top_20 <- bind_rows(
      dif1 %>%
        filter(regulate == Up) %>%
        arrange( padj, desc(abs(log2FoldChange))) %>%
        head(20),
      dif1 %>%
        filter(regulate == Down) %>%
        arrange( padj, desc(abs(log2FoldChange))) %>%
        head(20)
    )
    Exression<-dif1$regulate
    p<-ggplot(dif1, aes(log2FoldChange, -log10(padj)))+
      geom_point(aes(col=Exression,shape = Exression))+
      scale_shape_manual(values = c(25, 21,24))+
      #scale_shape_manual(values=c(Down=25, Up=24, Unchanged=21))+
      scale_color_manual(values=c(ann_colors[[1]][[1]],"grey",ann_colors[[1]][[2]]))+
      labs(title = " ")+
      geom_vline(xintercept=c(-logFC,logFC), colour="black", linetype="dashed")+
      geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+
      labs(x="log2(FoldChange)",y="-log10(Pvalue.adjust)",title = paste("Volcano plot of",tumor))+
      theme(axis.text=element_text(size=13,color = "black"),axis.title=element_text(size=13),
            plot.title = element_text(hjust = 0,size=13,face="plain"),
            legend.position ="bottom",
            panel.border = element_rect(color = "black", fill = NA),
            panel.background = element_rect(fill = "transparent"))+
      str(dif1)+#theme_bw()+
      geom_label_repel(data = top_20,max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                       aes(log2FoldChange, -log10(padj), label = gene),
                       size = 3, fill="#CCFFFF")
    if(geom_text){
      if(logFC!=0){
        min<-min(dif$log2FoldChange[dif$regulate == "Down"], na.rm = TRUE)
        max<-max(dif$log2FoldChange[dif$regulate == "Up"], na.rm = TRUE)
        p<-p + geom_text(x = (-logFC+min)/2, y = 0, label = paste(length(dif$regulate[dif$regulate=="Down"]),"genes"),
                         color = ann_colors[[1]][[1]], size = 4.5, vjust = 0)
        p<-p + geom_text(x = 0, y = 5, label = paste(length(dif$regulate[dif$regulate=="Unchanged"]),"genes"),
                         color = "grey", size = 3, vjust = 0)
        p<-p + geom_text(x = (logFC+max)/2, y = -log10(0.05), label = paste(length(dif$regulate[dif$regulate=="Up"]),"genes"),
                         color = ann_colors[[1]][[2]], size = 4.5, vjust = 0)
      }else{
        min<-min(dif$log2FoldChange[dif$regulate == "Down"], na.rm = TRUE)
        max<-max(dif$log2FoldChange[dif$regulate == "Up"], na.rm = TRUE)
        Pmax<-max(-log10(dif$padj), na.rm = TRUE)
        p<-p + geom_text(x = min/2, y = Pmax/2, label = paste(length(dif$regulate[dif$regulate=="Down"]),"genes"),
                         color = ann_colors[[1]][[1]], size = 4.5, vjust = 0)
        p<-p + geom_text(x = 0, y = -0.4, label = paste(length(dif$regulate[dif$regulate=="Unchanged"]),"genes"),
                         color = "grey", size = 3, vjust = 0)
        p<-p + geom_text(x = max/2, y = Pmax/2, label = paste(length(dif$regulate[dif$regulate=="Up"]),"genes"),
                         color = ann_colors[[1]][[2]], size = 4.5, vjust = 0)
      }

      #print(p)
    }
  }else if(diff_method=="edgeR"){
    top_20 <- bind_rows(
      dif1 %>%
        filter(regulate == Up) %>%
        arrange( FDR, desc(abs(logFC))) %>%
        head(20),
      dif1 %>%
        filter(regulate == Down) %>%
        arrange( FDR, desc(abs(logFC))) %>%
        head(20)
    )
    Exression<-dif1$regulate
    p<-ggplot(dif1, aes(logFC, -log10(FDR)))+
      geom_point(aes(col=Exression,shape = Exression))+
      scale_shape_manual(values = c(25, 21,24))+
      #scale_shape_manual(values=c(Down=25, Up=24, Unchanged=21))+
      scale_color_manual(values=c(ann_colors[[1]][[1]],"grey",ann_colors[[1]][[2]]))+
      labs(title = " ")+
      geom_vline(xintercept=c(-logFC,logFC), colour="black", linetype="dashed")+
      geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+
      labs(x="log2(FoldChange)",y="-log10(FDR)",title = paste("Volcano plot of",tumor))+
      theme(axis.text=element_text(size=13,color = "black"),axis.title=element_text(size=13),
            plot.title = element_text(hjust = 0,size=13,face="plain"),
            legend.position ="bottom",
            panel.border = element_rect(color = "black", fill = NA),
            panel.background = element_rect(fill = "transparent"))+
      str(dif1)+#theme_bw()+
      geom_label_repel(data = top_20,max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                       aes(logFC, -log10(FDR), label = gene),
                       size = 3, fill="#CCFFFF")
    if(geom_text){
      if(logFC!=0){
        min<-min(dif$logFC[dif$regulate == "Down"], na.rm = TRUE)
        max<-max(dif$logFC[dif$regulate == "Up"], na.rm = TRUE)
        p<-p + geom_text(x = (-logFC+min)/2, y = 0, label = paste(length(dif$regulate[dif$regulate=="Down"]),"genes"),
                         color = ann_colors[[1]][[1]], size = 4.5, vjust = 0)
        p<-p + geom_text(x = 0, y = 5, label = paste(length(dif$regulate[dif$regulate=="Unchanged"]),"genes"),
                         color = "grey", size = 3, vjust = 0)
        p<-p + geom_text(x = (logFC+max)/2, y = -log10(0.05), label = paste(length(dif$regulate[dif$regulate=="Up"]),"genes"),
                         color = ann_colors[[1]][[2]], size = 4.5, vjust = 0)
      }else{
        min<-min(dif$logFC[dif$regulate == "Down"], na.rm = TRUE)
        max<-max(dif$logFC[dif$regulate == "Up"], na.rm = TRUE)
        Pmax<-max(-log10(dif$FDR), na.rm = TRUE)
        p<-p + geom_text(x = min/2, y = Pmax/2, label = paste(length(dif$regulate[dif$regulate=="Down"]),"genes"),
                         color = ann_colors[[1]][[1]], size = 4.5, vjust = 0)
        p<-p + geom_text(x = 0, y = -0.4, label = paste(length(dif$regulate[dif$regulate=="Unchanged"]),"genes"),
                         color = "grey", size = 3, vjust = 0)
        p<-p + geom_text(x = max/2, y = Pmax/2, label = paste(length(dif$regulate[dif$regulate=="Up"]),"genes"),
                         color = ann_colors[[1]][[2]], size = 4.5, vjust = 0)
      }

      #print(p)
    }

  }

  a<-length(colnames(dif))
  dif$gene<-rownames(dif)
  rownames(dif)<-NULL
  dif<-dif[,c(a+1,1:a)]
  dif_id <-dif[dif$regulate %in% c("Up", "Down"), ]
  if(diff_method%in%c("limma","edgeR")){
    dif_id <- dif_id[order(-dif_id$logFC), ]
  }else if(diff_method=="DESeq2"){
    dif_id <- dif_id[order(-dif_id$log2FoldChange), ]
  }
  if(top50){
    if(length(dif$gene[dif$regulate=="Up"])>=50&&length(dif$gene[dif$regulate=="Down"])>=50){
      dif_id <-rbind(head(dif_id, 50),
                     tail(dif_id, 50))
    }else if(length(dif$gene[dif$regulate=="Up"])>=50&&length(dif$gene[dif$regulate=="Down"])<50){
      dif_id <-rbind(head(dif_id, 50),
                     tail(dif_id, length(dif$gene[dif$regulate=="Down"])))
    }else if(length(dif$gene[dif$regulate=="Up"])<50&&length(dif$gene[dif$regulate=="Down"])>=50){
      dif_id <-rbind(head(dif_id, length(dif$gene[dif$regulate=="Up"])),
                     tail(dif_id, 50))
    }else if(length(dif$gene[dif$regulate=="Up"])<50&&length(dif$gene[dif$regulate=="Down"])<50){
      dif_id <-rbind(head(dif_id, length(dif$gene[dif$regulate=="Up"])),
                     tail(dif_id, length(dif$gene[dif$regulate=="Down"])))
    }
  }
  annotation_row<-dif_id[,c("gene","regulate")]
  dif_ids <-unique(dif_id[, "gene"])
  dif_exprdata <- exprdata[dif_ids,]
  hmexp <- na.omit(dif_exprdata)
  pdata[[colnames(pdata)[[2]]]] <- factor(pdata[[colnames(pdata)[[2]]]], levels = Type)
  hmexp <- hmexp[, match(as.character(pdata[[colnames(pdata)[[1]]]]), colnames(hmexp))]
  annotation_col <-  data.frame(Gene= factor(c(rep(Type[[1]], as.numeric(length(pdata[[colnames(pdata)[[2]]]][pdata[[colnames(pdata)[[2]]]]==Type[[1]]]))),
                                               rep(Type[[2]],as.numeric(length(pdata[[colnames(pdata)[[2]]]][pdata[[colnames(pdata)[[2]]]]==Type[[2]]]))))))
  colnames(annotation_col)[colnames(annotation_col) == "Gene"] <-colnames(pdata)[[2]]
  rownames(annotation_col) <- colnames(hmexp)
  hmexp_sorted <- hmexp[match(as.character(annotation_row$gene), rownames(hmexp)), ]
  annotation_row<-data.frame(regulate=annotation_row[,-1])
  rownames(annotation_row)<-dif_id[, "gene"]
  if(tidyHeatmap){
    input1<-hmexp_sorted
    input1<-t(input1)
    input1<-as.data.frame(input1)
    num1<-as.numeric(length(colnames(input1)))
    input1$SAMPLE_ID<-rownames(input1)
    num2<-as.numeric(length(colnames(input1)))
    input1<-input1[,c(num2,1:num1)]
    rownames(input1)<-NULL
    input1<-merge(pdata,input1,by.x=colnames(pdata)[[1]],by.y=colnames(input1)[[1]],all.x=T)
    feas1<-colnames(input1)[3:length(colnames(input1))]
    condiction<-dif_id[,c("gene","regulate")]
    p1<-Mypackage::sig_Heatmap(input = input1, features = feas1,ID =colnames(input1)[[1]],show_plot=F,
                               condiction=condiction,id_condiction=colnames(condiction)[[1]],col_condiction=colnames(condiction)[[2]],
                               cols_group=c(ann_colors[[2]][[1]],ann_colors[[2]][[2]]),row_group=c(ann_colors[[1]][[2]],ann_colors[[1]][[1]]),
                               legend_show=TRUE,column_title_size=10,row_title_size=8,
                               heatmap_col=color,
                               #heatmap_col=c("#0505FA", "#FFFFFF", "#FA050D"),
                               group =colnames(pdata)[[2]],row_title="Regulate", scale = TRUE,name="Expression")
  }
  else{p1<-pheatmap(hmexp_sorted,
                    annotation_col = annotation_col,cutree_rows=2,
                    annotation_row=annotation_row,
                    annotation_colors = ann_colors,
                    color = color ,clustering_distance_rows = "correlation",
                    name = "Exression",
                    cluster_rows = T,#行聚类
                    cluster_cols = F,#列聚类
                    show_rownames = T,#行名
                    show_colnames = F,#列名
                    scale = "row",
                    fontsize = 12,
                    fontsize_row = 7,
                    fontsize_col = 6,
                    border = F)
  }
  #print(p1)
  gene <- dif$gene[dif$regulate %in% Regulate]
  gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
  data_all <-merge(gene,dif,by.x="SYMBOL",by.y="gene",all.x=T)
  if(diff_method%in%c("limma","edgeR")){
    data_all_sort <- data_all %>%
      arrange(desc(logFC))
    geneList = data_all_sort$logFC #把foldchange按照从大到小提取出来
  }else if(diff_method=="DESeq2"){
    data_all_sort <- data_all %>%
      arrange(desc(log2FoldChange))
    geneList = data_all_sort$log2FoldChange #把foldchange按照从大到小提取出来
  }
  names(geneList) <- data_all_sort$ENTREZID
  if(GO){
    if(GO.plot!="barplot"){
      gse.GO <- tryCatch({
        gseGO(
          geneList, # geneList
          ont = if(split) "ALL" else "BP",  # Optional "BP", "MF" and "CC" or "ALL"
          OrgDb = org.Hs.eg.db, # Human gene annotation
          keyType = "ENTREZID",
          pvalueCutoff = 0.05,
          pAdjustMethod = "BH" # P-value adjustment method
        )
      }, error = function(e) {
        message("GO GSEA analysis failed: ", e$message)
        return(NULL)
      })
      if (nrow(gse.GO@result)!=0) {
        if(GO.plot=="dotplot"){
          if(split){
            GO1 <- dotplot(gse.GO, color="pvalue", title=paste("GO enrichment of", tumor),
                           showCategory =5,split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale='free')
          }else{
            GO1 <- dotplot(gse.GO, color="pvalue", title=paste("GO enrichment of", tumor))
          }
        }else if (GO.plot=="heatplot"){
          GO1 <- heatplot(gse.GO, foldChange=geneList)+labs(title = paste("GO enrichment of", tumor))
        }else if (GO.plot=="cnetplot"){
          GO1 <- cnetplot(gse.GO, categorySize="pvalue", foldChange=geneList,colorEdge = TRUE)+
            labs(title = paste("GO enrichment of", tumor))
        }else if (GO.plot=="gseaplot"){
          print("GO enrichment")
          GO1 <-GSEAplot2(gse.GO,
                          geneSetID=if(length(gse.GO@result[["ID"]])>10) c(1:10) else c(1:c(1:length(gse.GO@result[["ID"]]))),
                          ES_geom = "line",legend.position ="none",pvalue_table=FALSE,
                          title = paste("GO enrichment of", tumor),rel_heights = rel_heights,
                          base_size = 20,Type=if(diff_method=="limma") c(Type[[2]],Type[[1]]) else Type )
        }

      } else {
        gse.GO<-NULL
        GO1 <-NULL
        message("GO plot is NULL, skipping print.")
      }
    }else if(GO.plot=="barplot"){
      gse.GO<-enrichGO(gene = names(geneList), #需要分析的基因的EntrezID
                       OrgDb = org.Hs.eg.db, #人基因数据库
                       pvalueCutoff =0.05, #设置pvalue界值
                       qvalueCutoff = 0.05, #设置qvalue界值(FDR校正后的p值）
                       ont=if(split) "all" else "BP",, #选择功能富集的类型，可选BP、MF、CC，这里选择all。
                       readable =T)
      if (nrow(gse.GO@result)!=0) {
        if(split){
          GO1 <- barplot(height=gse.GO, x = "GeneRatio",color="pvalue", title=paste("GO enrichment of", tumor),
                         showCategory =5,split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale='free')
        }else{
          GO1 <- barplot(height=gse.GO, x = "GeneRatio", color="pvalue", title=paste("GO enrichment of", tumor))
        }
      }else{
        gse.GO<-NULL
        GO1 <-NULL
        message("GO plot is NULL, skipping print.")
      }
    }
  }else{gse.GO<-NULL
  GO1 <-NULL}
  #print(GO1)
  if(KEGG){
    if (KEGG.plot!="barplot"){
      gse.KEGG <- tryCatch({
        gseKEGG(geneList,
                organism = "hsa",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH")
      }, error = function(e) {
        message("KEGG GSEA analysis failed: ", e$message)
        return(NULL)  # Return NULL or appropriate value
      })
      if (nrow(gse.KEGG@result)!=0) {
        if(KEGG.plot=="dotplot"){
          KEGG1 <- dotplot(gse.KEGG, color="pvalue", title=paste("KEGG enrichment of", tumor))
        }else if (KEGG.plot=="heatplot"){
          KEGG1 <- heatplot(gse.KEGG, foldChange=geneList)+labs(title = paste("KEGG enrichment of", tumor))
        }else if (KEGG.plot=="cnetplot"){
          KEGG1 <- cnetplot(gse.KEGG, categorySize="pvalue", foldChange=geneList,colorEdge = TRUE)+
            labs(title = paste("KEGG enrichment of", tumor))
        }else if (KEGG.plot=="gseaplot"){
          print("KEGG enrichment")
          KEGG1 <-GSEAplot2(gse.KEGG,
                            geneSetID=if(length(gse.KEGG@result[["ID"]])>10) c(1:10) else c(1:c(1:length(gse.KEGG@result[["ID"]]))),
                            ES_geom = "line",legend.position ="none",pvalue_table=FALSE,
                            title = paste("KEGG enrichment of", tumor),rel_heights = rel_heights,
                            base_size = 20,Type=if(diff_method=="limma") c(Type[[2]],Type[[1]]) else Type)

        }
        #print(KEGG1)
      } else {
        message("KEGG plot is NULL, skipping print.")
        gse.KEGG<-NULL
        KEGG1 <-NULL
      }
    }else if (KEGG.plot=="barplot"){
      gse.KEGG<-enrichKEGG(gene = names(geneList), #需要分析的基因的EntrezID
                           organism = "hsa",
                           keyType = "kegg",
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           #universe,
                           minGSSize = 10,
                           maxGSSize = 500,
                           qvalueCutoff = 0.2,
                           use_internal_data = FALSE
      )
      if (nrow(gse.KEGG@result)!=0) {
        gse.KEGG1<-gse.KEGG
        gse.KEGG1@result<-gse.KEGG1@result[gse.KEGG1@result$category == "Organismal Systems", ]
        KEGG1 <- barplot(height=gse.KEGG1, x = "geneRatio",color="pvalue", title=paste("KEGG enrichment of", tumor),
                         showCategory =10)
      }else{
        message("KEGG plot is NULL, skipping print.")
        gse.KEGG<-NULL
        KEGG1 <-NULL
      }
    }
  }else{gse.KEGG<-NULL
  KEGG1 <-NULL}
  results <- list(
    method=diff_method,
    cutoff=list(logFC=logFC,P.Value=P.Value),
    dif = dif,
    dif_id = dif_id,
    volcano_plot = p,
    hmexp=hmexp_sorted,
    heatmap = p1,
    GO=gse.GO,
    GO_plot = if (!is.null(GO1)) GO1 else NULL,
    KEGG=gse.KEGG,
    KEGG_plot = if (!is.null(KEGG1)) KEGG1 else NULL
  )
  return(results)
}

