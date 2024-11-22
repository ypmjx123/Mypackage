#' Pathway Mutation Visualization
#'
#' This function visualizes the mutation rate of tumor pathways based on SNV and gene data.
#'
#' @param exprdata    is the data of gene expression for samples
#' @param pdata       A data frame that contains groups of genes.例如data.frame(SMPLE_ID,Group)
#' @param Type        A character vector indicating the type.eg Type=c("Wild", "Mut")
#' @param tumor       the source of exprdata
#' @param contrastfml the Type for limma; eg contrastfml="Wild - Mut"
#' @param P.Value     p阈值for差异分析
#' @param logFC       logFC阈值for差异分析
#' @param ann_colors  ComplexHeatmap::pheatmap 参数for sample annotation
#' @param color       ComplexHeatmap::pheatmap 参数gene expression
#' @param Regulate    select Regulate=c("Up","Down") for pathway enrichment
#'
#' @return A list of differential analysis results,volcano,heatmap and enrichment plot of GO AND KEGG
#' @export
#' @examples
#' data("exp_CRC")
#' data("Gene_group_CRC1")
#' result<-limma.dif.visual(exprdata=exp_CRC[,-2],
#'                            pdata=Gene_group_CRC1,
#'                            Type=c("Wild", "Mut"),
#'                            contrastfml= "Wild - Mut",
#'                            tumor="CRC TCGA",
#'                            P.Value=0.05,
#'                            logFC=0.5,
#'                            color= colorRampPalette(c("#0505FA", "#FFFFFF", "#FA050D"))(20),
#'                            ann_colors = list(regulate = c(Down = "#1B9E77", Up = "#D95F02"),
#'                            PREX2 = c(Wild = "#757575", Mut = "#FF4040")),
#'                            Regulate=c("Up","Down"))
#' print(result)
limma.dif.visual<-function (exprdata, pdata,Type, tumor,contrastfml,
                            P.Value,logFC,ann_colors,color,Regulate)
{
  library(limma)
  library(tidyverse)
  library(ggplot2)
  library(ComplexHeatmap)
  library(dplyr)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  #library(pathview)
  library(enrichplot)
  library(base)
  exprdata<-exprdata[!is.na(exprdata[,2]),]
  if (any(duplicated(exprdata[[1]]))) {
    print(" expression data genes duplicated ")
    exprdata <- aggregate(as.formula(paste0(".", " ~ ", colnames(exprdata)[[1]])), data = exprdata, FUN = mean)
  }
  rownames(exprdata)<-exprdata[[1]]
  exprdata<-exprdata[,-1]
  if (!all(colnames(exprdata) == pdata[, 1])) {
    print(" expression data do not match pdata, we conduct it")
    exprdata<- exprdata[, colnames(exprdata) %in% pdata[, 1], drop = FALSE]
    pdata1<-data.frame(SAMPLE_ID=colnames(exprdata))
    pdata1<-merge(pdata,pdata1,by="SAMPLE_ID",all.y=T)
    pdata<-pdata1
  }
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
  Exression<-dif$regulate
  p<-ggplot(dif, aes(logFC, -log10(P.Value)))+
    geom_point(aes(col=Exression))+
    scale_color_manual(values=c("green","grey","red"))+
    labs(title = " ")+
    geom_vline(xintercept=c(-logFC,logFC), colour="black", linetype="dashed")+
    geom_hline(yintercept = -log10(0.05),colour="black", linetype="dashed")+
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
    labs(x="log2(FoldChange)",y="-log10(Pvalue)",
         title = paste("Volcano plot of",tumor))+
    theme(axis.text=element_text(size=13),axis.title=element_text(size=13))+
    str(dif)+theme_bw()
  min<-min(dif$logFC[dif$regulate == "Down"], na.rm = TRUE)
  max<-max(dif$logFC[dif$regulate == "Up"], na.rm = TRUE)
  p<-p + geom_text(x = (-logFC+min)/2, y = 0, label = paste(length(dif$regulate[dif$regulate=="Down"]),"genes"),
                   color = "green", size = 4, vjust = 0)
  p<-p + geom_text(x = 0, y = 5, label = paste(length(dif$regulate[dif$regulate=="Unchanged"]),"genes"),
                   color = "grey", size = 4, vjust = 0)
  p<-p + geom_text(x = (logFC+max)/2, y = -log10(0.05), label = paste(length(dif$regulate[dif$regulate=="Up"]),"genes"),
                   color = "red", size = 4, vjust = 0)
  print(p)
  a<-length(colnames(dif))
  dif$gene<-rownames(dif)
  rownames(dif)<-NULL
  dif<-dif[,c(a+1,1:a)]
  dif_id <-dif[dif$regulate %in% c("Up", "Down"), ]
  dif_id <- dif_id[order(-dif_id$logFC), ]
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
  annotation_row<-dif_id[,c(1,8)]
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
  p1<-pheatmap(hmexp_sorted,
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
  print(p1)
  gene <- dif$gene[dif$regulate %in% Regulate]
  gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
  data_all <-merge(gene,dif,by.x="SYMBOL",by.y="gene",all.x=T)
  data_all_sort <- data_all %>%
    arrange(desc(logFC))
  geneList = data_all_sort$logFC #把foldchange按照从大到小提取出来
  names(geneList) <- data_all_sort$ENTREZID
  gse.GO <- tryCatch({
    gseGO(
      geneList, # geneList
      ont = "BP",  # Optional "BP", "MF" and "CC" or "ALL"
      OrgDb = org.Hs.eg.db, # Human gene annotation
      keyType = "ENTREZID",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH" # P-value adjustment method
    )
  }, error = function(e) {
    message("GO GSEA analysis failed: ", e$message)
    return(NULL)
  })

  if (!is.null(gse.GO)) {
    GO1 <- dotplot(gse.GO, color="pvalue", title=paste("GO enrichment of", tumor))
    print(GO1)
  } else {
    message("GO plot is NULL, skipping print.")
  }
  gse.KEGG <- tryCatch({
    gseKEGG(geneList,
            organism = "hsa",
            pvalueCutoff = 0.05,
            pAdjustMethod = "BH")
  }, error = function(e) {
    message("KEGG GSEA analysis failed: ", e$message)
    return(NULL)  # Return NULL or appropriate value
  })

  if (!is.null(gse.KEGG)) {
    KEGG1 <- dotplot(gse.KEGG, color="pvalue", title=paste("KEGG enrichment of", tumor))
    print(KEGG1)
  } else {
    message("KEGG plot is NULL, skipping print.")
  }
  results <- list(
    dif = dif,
    dif_id = dif_id,
    volcano_plot = p,
    heatmap = p1,
    GO_plot = if (!is.null(GO1)) GO1 else NULL,
    KEGG_plot = if (!is.null(KEGG1)) KEGG1 else NULL
  )
  return(results)
}
