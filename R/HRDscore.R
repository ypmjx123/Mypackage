#' Scar score
#'
#' Determining genomic scar score (telomeric allelic imbalance, loss-off heterozigosity, large-scle transitions), signs of homologous recombination deficiency
#'
#' @param dirpath THE dir path of Allele specific Copy Number Segment files. sangerbox->数据下载->TCGA->cohort->Copy Number Variation(拷贝数)->Allele-specific Copy Number Segment
#' @param file_id_map file name match the sample ID. 一般在下载的文件最后一个,file_id_map_XX.TXT
#' @param ploidy Optional parameter, may be used if the ploidy of the sample is known. default is NULL, if have set ploidy=1
#' @param tumor  Allele specific Copy Number Segment files数据来源,当ploidy非NULL时,需提供; data("tumor_ploidy"),unique(tumor_ploidy$Tumor)
#' @param reference reference genome: either grch38 or grch37 or mouse. Default is grch38.

#'
#' @returns HRD scores
#' @export
#'
##' @author Pengmin Yang
##'
#' @examples
#' root_dir <-system.file("HRD", package = "Mypackage")
#' data("metadata")
#' HRD<-HRDscore(dirpath=root_dir,file_id_map=metadata,ploidy=NA,reference="grch38")
#' head(HRD)
#'
HRDscore<-function(dirpath="E:/R/TCGA/HRD/TCGA/OV-Allele-specific_Copy_Number_Segment/",
                   file_id_map=metadata,reference,
                   ploidy=NULL,tumor="OV"){
  library(data.table)
  library(scarHRD)
  cnv_file <- dir(dirpath)
  file <- paste0(dirpath,"/",cnv_file)
  file_list <- list()
  fred_cnv <- function(x){
    # 读取文件
    file_list[[x]] <- data.table::fread(file = x, data.table = FALSE)
    # 提取文件路径中第一个 '/' 之前的字符作为 file.id
    file.id_value <- sub("^.*/([^/]+)\\.TXT$", "\\1", x)
    # 添加到数据框中
    file_list[[x]]$file.id <- file.id_value
    # 返回数据框
    file_list[[x]]
  }
  library(dplyr)
  cnv_OV <- lapply(file, fred_cnv)
  #行合并
  cnv_OV <- do.call(dplyr::bind_rows, cnv_OV)
  cnv_OV$GDC_Aliquot<-cnv_OV$file.id
  cnv_OV$file.id<-NULL
  if(is.null(ploidy)){
    colnames(cnv_OV)<-c("SampleID", "Chromosome", "Start_position", "End_position","total_cn", "A_cn", "B_cn")
  }else{
    cnv_OV$ploidy<-NA
    colnames(cnv_OV)<-c("SampleID", "Chromosome", "Start_position", "End_position","total_cn", "A_cn", "B_cn","ploidy")
  }
  metadata<-file_id_map
  #cnv_OV$Sample <- metadata$cases.submitter_id[match(cnv_OV$SampleID,metadata$aliquots.aliquot_id)]
  cnv_OV$Sample <- metadata$submitter_id[match(cnv_OV$SampleID,metadata$file_id)]
  #rownames(cnv_OV)<-cnv_OV$SampleID
  cnv_OV$SampleID<-cnv_OV$Sample
  cnv_OV$Sample<-NULL
  if(is.null(ploidy)){
    cnv_OV$ploidy<-NA
    HRDs <- data.frame(matrix(ncol = 5, nrow = 0))
    colnames(HRDs) <- c("SampleID", "LOH", "TAI", "LST", "HRDsum")
    for(Sample in unique(cnv_OV$SampleID)){
      cnv_OV1<-cnv_OV[cnv_OV$SampleID==Sample,]
      rio::export(cnv_OV1,format = "txt")
      scar_score("cnv_OV1.tsv", reference = reference, chr.in.names=T, seqz = F, ploidy = ploidy,
                 sizelimitLOH = 1.5e+07, outputdir = NULL)
      HRD<-rio::import(paste0(Sample,"_HRDresults.txt"))
      file.remove(paste0(Sample, "_HRDresults.txt"))
      colnames(HRD)<-c("SampleID", "LOH", "TAI", "LST", "HRDsum")
      HRDs<-rbind(HRDs,HRD)
    }
    colnames(HRDs)<-c("SampleID","LOH","TAI","LST","HRDsum")
    HRDs$HRD<-ifelse(HRDs$HRDsum>=42,"HRD+","HRD-")
  }else{
    HRDs <- data.frame(matrix(ncol = 17, nrow = 0))
    colnames(HRDs) <- c("SampleID", "LOH", "TAI","Mean size", "Interstitial AI","Mean Size","Whole chr AI","Telomeric LOH","Mean size","Interstitial LOH", "Mean Size","Whole chr LOH","Ploidy","Aberrant cell fraction","LST", "HRDsum","adjustHRDsum")
    for(Sample in unique(cnv_OV$SampleID)){
      cnv_OV1<-cnv_OV[cnv_OV$SampleID==Sample,]
      data("tumor_ploidy")
      if (any(grepl(Sample, tumor_ploidy$SampleName) & grepl(tumor, tumor_ploidy$Tumor))) {
        ploidy <- tumor_ploidy$ploidy[which(grepl(Sample, tumor_ploidy$SampleName) & grepl(tumor, tumor_ploidy$Tumor))]
      } else {
        ploidy <- 0
      }
      cnv_OV1$ploidy<-ploidy
      rio::export(cnv_OV1,format = "txt")
      scar_score("cnv_OV1.tsv", reference = reference, chr.in.names=T, seqz = F, ploidy = ploidy,
                 sizelimitLOH = 1.5e+07, outputdir = NULL)
      HRD<-rio::import(paste0(Sample,"_HRDresults.txt"))
      file.remove(paste0(Sample, "_HRDresults.txt"))
      colnames(HRD)<-c("SampleID", "LOH", "TAI",colnames(HRD)[4:14], "LST", "HRDsum","adjusted-HRDsum")
      #sum_adjust<-LST-15.5*ploidy+LOH+TAI
      HRD$Ploidy <-ploidy
      HRDs<-rbind(HRDs,HRD)
    }
    colnames(HRDs)<-c("SampleID", "LOH", "TAI", colnames(HRD)[4:14],"LST","HRDsum", "adjusted-HRDsum")
    HRDs$HRD<-ifelse(HRDs$`adjusted-HRDsum`>=median(HRDs$`adjusted-HRDsum`),"HRD+","HRD-")
  }
  rio::export(HRDs,format = "xlsx")
  return(HRDs)
}
