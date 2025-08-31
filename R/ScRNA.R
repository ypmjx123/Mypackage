#' 单细胞RNA数据分析处理
#'
#' 该函数用于处理单细胞RNA测序数据，包含数据标准化、特征选择和降维分析。
#'
#' @param file.type 输入文件类型，可选"10X"、"h5"或"txt/csv/tsv"，默认"10X"
#' @param file 文件路径/文件夹路径，10X数据需传文件夹路径
#' @param multifile 是否处理多个文件，默认FALSE
#' @param file_folder 多文件模式下的父文件夹路径
#' @param project 项目名称标识
#' @param organism 物种类型，可选"Human"或"Mouse"，默认Human
#' @param clusters 是否进行差异表达分析，默认FALSE
#' @param reduction 降维可视化方法，可选"tsne"，"umap"，“pca”，默认tsne，If not specified, first searches for umap, then tsne, then pca
#' @param interaction.method 单细胞数据整合方法，c("interaction","harmony")
#'
#' @return  处理后的Seurat对象
#' @export
#' @author Pengmin Yang
ScRNA<-function(file.type="10X",file="GSM3729179/",multifile=FALSE,
                file_folder,project,organism="Human",clusters=FALSE,reduction = "tsne",
                interaction.method=c("interaction","harmony")){
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(SingleR)
  library(celldex)
  library(pheatmap)
  library(RColorBrewer)
  library(Matrix)
  library("scrapper")
  library(BiocParallel)
  library(harmony)
  if(multifile){
    read_txt_to_seurat <- function(file_path) {
      # 读取txt文件
      data <- read.table(file_path, header = TRUE, row.names = 1, sep = "\t")
      #data=fread(file_path,data.table = F)
      #rownames(data) <- data[,1]
      #data <- data[,-1]
      # 转换为稀疏矩阵
      data <- as(as.matrix(data), "dgCMatrix")
      # 创建Seurat对象
      seurat_obj <- CreateSeuratObject(counts = data,
                                       project = regmatches(file_path, regexpr("GSM[0-9]+", file_path)),
                                       min.cells = 5,
                                       min.features = 300)
      seurat_obj<-RenameCells(seurat_obj ,add.cell.id = regmatches(file_path, regexpr("GSM[0-9]+", file_path)))
      seurat_obj@meta.data[["orig.ident"]]<-regmatches(file_path, regexpr("GSM[0-9]+", file_path))
      Idents(seurat_obj) <- as.factor(seurat_obj$orig.ident)
      return(seurat_obj)
    }
    read_10x_to_seurat <- function(folder_path) {
      # 读取10X数据
      data <- Read10X(data.dir = folder_path)
      # 创建Seurat对象
      seurat_obj <- CreateSeuratObject(counts = data,
                                       project = gsub(paste0(file_folder,"/"),"",folder_path),
                                       min.cells = 5,
                                       min.features = 300)
      seurat_obj<-RenameCells(seurat_obj ,add.cell.id = gsub(paste0(file_folder,"/"),"",folder_path))
      return(seurat_obj)
    }
    read_h5_to_seurat <- function(file_path) {
      data =  Read10X_h5(file_path)
      seurat_obj <- CreateSeuratObject(counts = data,
                                       project = project,
                                       min.cells = 5,
                                       min.features = 300)
      return(seurat_obj)
    }
    file_folder <- file_folder
    # 获取文件夹中的所有文件和子文件夹
    files <- list.files(file_folder, full.names = T)
    # 初始化存储Seurat对象的列表
    seurat_list <- list()
    # 遍历文件夹中的文件
    for (file in files) {
      if (grepl("\\.(txt|csv|tsv|xls)$", file, ignore.case = TRUE)) {
        # 如果是txt文件
        print(paste("reading",file))
        seurat_obj <- read_txt_to_seurat(file)
        print("readed")
      } else if (grepl("\\.h5$", file)) {
        # 如果是h5文件
        print(paste("reading",file))
        seurat_obj <- read_h5_to_seurat(file)
        print("readed")
      }else if (file.info(file)$isdir) {
        # 如果是10X文件夹
        print(paste("reading",file))
        seurat_obj <- read_10x_to_seurat(file)
        print("readed")
      } else {
        next
      }
      # 将Seurat对象添加到列表
      #seurat_list[[basename(file)]] <- seurat_obj
      seurat_list[[regmatches(file, regexpr("GSM[0-9]+", file))]] <- seurat_obj
    }
    # 合并所有Seurat对象
    seurat_list <- lapply(X = seurat_list, FUN = function(x) {
      x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
      rb.genes <- rownames(x)[grep("^Rp[sl]",rownames(x))]
      C<-GetAssayData(object = x, slot = "counts")
      percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
      x <- AddMetaData(x, percent.ribo, col.name = "percent.ribo")
      QC_plot<-VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 2)
      print(QC_plot)
      if(organism=="Human"){
        x <- subset(x,
                    subset = nFeature_RNA > quantile(x@meta.data[["nFeature_RNA"]], 0.25) - 1.5 * IQR(x@meta.data[["nFeature_RNA"]])&nFeature_RNA < quantile(x@meta.data[["nFeature_RNA"]], 0.75) + 1.5 * IQR(x@meta.data[["nFeature_RNA"]])
                    & nCount_RNA < quantile(x@meta.data[["nCount_RNA"]], 0.75) + 1.5 * IQR(x@meta.data[["nCount_RNA"]]) &nCount_RNA > quantile(x@meta.data[["nCount_RNA"]], 0.25) - 1.5 * IQR(x@meta.data[["nCount_RNA"]])
                    & percent.mt < quantile(x@meta.data[["percent.mt"]], 0.75) + 1.5 * IQR(x@meta.data[["percent.mt"]])& percent.mt> quantile(x@meta.data[["percent.mt"]], 0.25) - 1.5 * IQR(x@meta.data[["percent.mt"]])
        )
      }else if (organism=="Mouse"){
        x <- subset(x,
                    subset = nFeature_RNA > quantile(x@meta.data[["nFeature_RNA"]], 0.25) - 1.5 * IQR(x@meta.data[["nFeature_RNA"]])&nFeature_RNA < quantile(x@meta.data[["nFeature_RNA"]], 0.75) + 1.5 * IQR(x@meta.data[["nFeature_RNA"]])
                    & nCount_RNA < quantile(x@meta.data[["nCount_RNA"]], 0.75) + 1.5 * IQR(x@meta.data[["nCount_RNA"]]) &nCount_RNA > quantile(x@meta.data[["nCount_RNA"]], 0.25) - 1.5 * IQR(x@meta.data[["nCount_RNA"]])
        )
      }
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })
    if(interaction.method=="interaction"){
      # select features that are repeatedly variable across datasets for integration
      features <- SelectIntegrationFeatures(object.list = seurat_list)
      immune.anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)
      # this command creates an 'integrated' data assay
      combined_seurat <- IntegrateData(anchorset = immune.anchors)
      DefaultAssay(combined_seurat) <- "integrated"
      combined_seurat <- ScaleData(combined_seurat, verbose = FALSE)
      combined_seurat <- RunPCA(combined_seurat, npcs = 30, verbose = FALSE)
      p<-ElbowPlot(combined_seurat)
      print(p)
      combined_seurat <- RunUMAP(combined_seurat, reduction = "pca", dims = 1:20)
      combined_seurat = RunTSNE(combined_seurat, dims = 1:20)
      combined_seurat<- FindNeighbors(combined_seurat, reduction = "pca", dims = 1:20)
      combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)
      p1 <- DimPlot(combined_seurat, reduction = reduction, group.by = "orig.ident")
      p2 <- DimPlot(combined_seurat, reduction = reduction, label = TRUE, repel = TRUE)
      print(p1 + p2)
      if(clusters){
        combined_seurat.markers <- FindAllMarkers(combined_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
      }else(combined_seurat.markers <-NULL)
      if(organism=="Human"){
        ref <- celldex::HumanPrimaryCellAtlasData()#human
      }else if(organism=="Mouse"){
        ref<-celldex::MouseRNAseqData()#Mouse
      }
      pred.scRNA <- SingleR(test = combined_seurat@assays$integrated@data,
                            ref = ref,
                            labels = ref$label.main,
                            clusters = combined_seurat@active.ident)
      plotScoreHeatmap(pred.scRNA, clusters=pred.scRNA@rownames, fontsize.row = 9,show_colnames = T)
      new.cluster.ids <- pred.scRNA$pruned.labels#也可手动注释
      names(new.cluster.ids) <- levels(combined_seurat)
      combined_seurat <- RenameIdents(combined_seurat,new.cluster.ids)
      combined_seurat@meta.data$label=combined_seurat@active.ident
    }else if (interaction.method=="harmony"){
      combined_seurat <- merge(x=seurat_list[[1]], y = seurat_list[-1])
      combined_seurat[["percent.mt"]] <- PercentageFeatureSet(combined_seurat, pattern = "^MT-")
      combined_seurat <- NormalizeData(combined_seurat) %>%
        FindVariableFeatures()
      gc()
      combined_seurat <- ScaleData(combined_seurat, features = VariableFeatures(combined_seurat))
      combined_seurat <- RunPCA(combined_seurat, features = VariableFeatures(combined_seurat))
      combined_seurat <- IntegrateLayers(combined_seurat,HarmonyIntegration)%>%
        FindNeighbors(reduction = 'harmony', dims = 1:15)%>%
        FindClusters(resolution = 0.5)%>%
        RunUMAP(reduction = "harmony", dims = 1:15)%>%
        RunTSNE(reduction = "harmony", dims = 1:15)
      combined_seurat = JoinLayers(combined_seurat)
      if(clusters){
        combined_seurat.markers <- FindAllMarkers(combined_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
      }else(combined_seurat.markers <-NULL)
      if(organism=="Human"){
        ref <- celldex::HumanPrimaryCellAtlasData()#human
      }else if(organism=="Mouse"){
        ref<-celldex::MouseRNAseqData()#Mouse
      }
      pred.scRNA <- SingleR(test = combined_seurat@assays$RNA$data,
                            ref = ref,
                            labels = ref$label.main,
                            clusters = combined_seurat@active.ident)
      if(any(is.na(pred.scRNA$pruned.labels))){
        new.cluster.ids <- pred.scRNA$labels
      }else{
        new.cluster.ids <- pred.scRNA$pruned.labels
      }
      names(new.cluster.ids) <- levels(combined_seurat)
      combined_seurat<- RenameIdents(combined_seurat,new.cluster.ids)
      combined_seurat@meta.data$label=combined_seurat@active.ident
    }
    result_all<-list(seurat_list=seurat_list,
                     SeuratObject=combined_seurat,
                     Clusters=combined_seurat.markers,
                     Cell_annotion=pred.scRNA)
  }else{
    if(file.type=="txt/csv/tsv"){
      #txt/csv/tsv数据读取
      counts=fread(file,data.table = F)
      rownames(counts) <- counts[,1]
      counts <- counts[,-1]
    }else if(file.type=="10X"){
      #标准10X
      counts<-Read10X(file)
    }else if(file.type=="h5"){
      counts =  Read10X_h5(file)
    }
    sce=CreateSeuratObject(counts =  counts ,
                           project =  project ,
                           min.cells = 5,
                           min.features = 300,)

    sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
    rb.genes <- rownames(sce)[grep("^Rp[sl]",rownames(sce))]
    C<-GetAssayData(object = sce, slot = "counts")
    percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
    sce <- AddMetaData(sce, percent.ribo, col.name = "percent.ribo")
    QC_plot<-VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 2)
    QC_plot
    plot1<-FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2<-FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    CombinePlots(plots = list(plot1, plot2))
    #QC
    if(organism=="Human"){
      sce <- subset(sce,
                    subset = nFeature_RNA > quantile(sce@meta.data[["nFeature_RNA"]], 0.25) - 1.5 * IQR(sce@meta.data[["nFeature_RNA"]])&nFeature_RNA < quantile(sce@meta.data[["nFeature_RNA"]], 0.75) + 1.5 * IQR(sce@meta.data[["nFeature_RNA"]])
                    & nCount_RNA < quantile(sce@meta.data[["nCount_RNA"]], 0.75) + 1.5 * IQR(sce@meta.data[["nCount_RNA"]]) &nCount_RNA > quantile(sce@meta.data[["nCount_RNA"]], 0.25) - 1.5 * IQR(sce@meta.data[["nCount_RNA"]])
                    & percent.mt < quantile(sce@meta.data[["percent.mt"]], 0.75) + 1.5 * IQR(sce@meta.data[["percent.mt"]])& percent.mt> quantile(sce@meta.data[["percent.mt"]], 0.25) - 1.5 * IQR(sce@meta.data[["percent.mt"]])
      )
    }else if (organism=="Mouse"){
      sce <- subset(sce,
                    subset = nFeature_RNA > quantile(sce@meta.data[["nFeature_RNA"]], 0.25) - 1.5 * IQR(sce@meta.data[["nFeature_RNA"]])&nFeature_RNA < quantile(sce@meta.data[["nFeature_RNA"]], 0.75) + 1.5 * IQR(sce@meta.data[["nFeature_RNA"]])
                    & nCount_RNA < quantile(sce@meta.data[["nCount_RNA"]], 0.75) + 1.5 * IQR(sce@meta.data[["nCount_RNA"]]) &nCount_RNA > quantile(sce@meta.data[["nCount_RNA"]], 0.25) - 1.5 * IQR(sce@meta.data[["nCount_RNA"]])
      )
    }

    #归一
    sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
    sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
    VariableFeaturePlot(sce) + LabelPoints(plot = VariableFeaturePlot(sce) , points =  head(VariableFeatures(sce), 10), repel = TRUE)
    sce <- ScaleData(sce , features = rownames(sce ))
    sce <- RunPCA(sce, features = VariableFeatures(object = sce))
    print(sce[["pca"]], dims = 1:5, nfeatures = 5)  #查看前五个高变基因
    VizDimLoadings(sce, dims = 1:2, reduction = "pca")
    DimHeatmap(sce, dims = 1:15, cells = 500, balanced = TRUE)
    sce<- FindNeighbors(sce, dims = 1:10)
    sce <- FindClusters(sce, resolution = 0.5)
    sce <- RunUMAP(sce, dims = 1:10)
    sce = RunTSNE(sce, dims = 1:30)
    if(clusters){
      sce.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    }else(sce.markers<-NULL)
    #sce.markers1<-merge(sce.markers1,data[data$Tissue_standard=="Urinary bladder",c("Tissue_standard","Cell_Marker","Cell_standard")],by.x="gene",by.y="Cell_Marker",all.x=T)
    #top10 <- sce.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
    #DoHeatmap(sce, features = top10$gene) + NoLegend()
    #细胞注释
    if(organism=="Human"){
      hpca.ref <- celldex::HumanPrimaryCellAtlasData()#human
    }else if(organism=="Mouse"){
      hpca.ref<-celldex::MouseRNAseqData()#Mouse
    }
    star <- as.SingleCellExperiment(DietSeurat(sce))
    hpca.main <- SingleR(test = star,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
    hpca.fine <- SingleR(test = star,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)
    hpca.main$pruned.labels<-ifelse(is.na(hpca.main$pruned.labels),"unknown",hpca.main$pruned.labels)
    hpca.fine$pruned.labels<-ifelse(is.na(hpca.fine$pruned.labels),"unknown",hpca.fine$pruned.labels)
    sce@meta.data$hpca.main   <- hpca.main$labels
    sce@meta.data$hpca.fine   <- hpca.fine$labels
    sce <- SetIdent(sce, value = "hpca.main")
    ncluster <- length(unique(sce[[]]$hpca.main))
    mycol <- colorRampPalette(brewer.pal(8, "Set2"))(ncluster)
    p<-DimPlot(sce,
               label = F , repel = T, reduction = reduction,
               label.size = 5,group.by = "hpca.main",
               cols = mycol)
    #p<-DimPlot(sce, group.by = "labels",reduction = reduction,label=F)#pca,umap,tsne
    #p[[1]][["labels"]][["title"]]<-NULL
    p
    result_all<-list(
      SeuratObject=sce,
      QC_plot= QC_plot,
      Clusters=sce.markers,
      Cell_annotion=list(main=hpca.main,
                         fine=hpca.fine),
      Cell_plot=p
    )
  }
  return(result_all)
}


