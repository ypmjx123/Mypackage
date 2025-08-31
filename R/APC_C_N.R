#' Sample grouping with APC gene mutation status
#'
#' This function groups samples based on mutations at both ends of the APC gene according to MAF format data.
#'
#' @param SNV            A data frame containing SNV data.即MAF格式的数据
#' @param is.TCGA        is TCGA data or not
#' @param sample_id_data data.frame格式入组患者编号，表头名称需为Tumor_Sample_Barcode
#' @return A data frame of sample grouping with APC gene
#' @export
#' @references
#' 1. Mondaca S, Walch H, Nandakumar S, Chatila WK, Schultz N, Yaeger R. \bold{Specific Mutations in APC,
#'    but Not Alterations in DNA Damage Response, Associate With Outcomes of Patients With Metastatic Colorectal Cancer.}
#'    \emph{Gastroenterology}. 2020;159(5):1975-1978.e4. doi:10.1053/j.gastro.2020.07.041\cr
#' 2. Peng H, Ying J, Zang J, et al. \bold{Specific Mutations in APC,
#'    with Prognostic Implications in Metastatic Colorectal Cancer.}
#'    \emph{Cancer Res Treat}. 2023;55(4):1270-1280. doi:10.4143/crt.2023.415
##' @author Pengmin Yang
#' @examples
#' data("SNV_data")
#' data("gene_group_data")
#' colnames(gene_group_data)[colnames(gene_group_data) == "SAMPLE_ID"] <- "Tumor_Sample_Barcode"
#' result <- APC_C_N(SNV = SNV_data,
#'                   is.TCGA=TRUE,
#'                   sample_id_data=gene_group_data)
#' print(result)
APC_C_N<-function(SNV,is.TCGA=TRUE,sample_id_data){
  library(dplyr)
  library(stringr)
  sample_noAPC<-unique(SNV$Tumor_Sample_Barcode[SNV$Hugo_Symbol!= "APC"])
  data_APC_Type1 <- data.frame(Tumor_Sample_Barcode = sample_noAPC)
  data_APC<-SNV[SNV$Hugo_Symbol== "APC",]
  if(is.TCGA==FALSE){
    for (i in 1:nrow(data_APC)) {
      if(is.na(data_APC[i, "HGVSp"])){
        if ((grepl(">", data_APC[i, "HGVSc"])) &&
            (!grepl("del|ins|dup", data_APC[i, "HGVSc"])) &&
            (!grepl("_", data_APC[i, "HGVSc"]))) {
          num <- str_extract(data_APC[i, "HGVSc"], "(?<=\\+|-)(\\d+)(?=[A-Za-z])")
          num <- ifelse(is.na(num), 0, as.numeric(num))
          if (num <= 2) {
            if (data_APC[i, "Start_Position"]< 112175489) {
              data_APC[i,"APC_Type"] <- "APC_N"
            } else {
              data_APC[i,"APC_Type"] <- "APC_C"
            }
          } else {
            data_APC[i,"APC_Type"] <- "other mutation"
          }
        } else if ((!grepl(">", data_APC[i, "HGVSc"])) &&
                   (grepl("del|ins|dup", data_APC[i, "HGVSc"])) &&
                   (grepl("_", data_APC[i, "HGVSc"]))){
          num1 <- str_extract(data_APC[i, "HGVSc"], "(?<=\\+|-)(\\d+)(?=_)")
          num1 <- ifelse(is.na(num1), 0, as.numeric(num1))
          num2 <- str_extract(data_APC[i, "HGVSc"], "(?<=\\+|-)(\\d+)(?=del|ins|dup)")
          num2 <- ifelse(is.na(num2), 0, as.numeric(num2))
          if ((num1 <= 2 || num2 <= 2) || (num1 <= 2 && num2 <= 2)) {
            if (data_APC[i,"Start_Position"] < 112175489) {
              data_APC[i,"APC_Type"] <- "APC_N"
            } else {
              data_APC[i,"APC_Type"] <- "APC_C"
            }
          } else {
            data_APC[i,"APC_Type"] <- "other mutation"
          }
        }else if ((!grepl(">", data_APC[i, "HGVSc"])) &&
                  (grepl("del|ins|dup", data_APC[i, "HGVSc"])) &&
                  (!grepl("_", data_APC[i, "HGVSc"]))){
          num3 <- str_extract(data_APC[i, "HGVSc"], "(?<=\\+|-)(\\d+)(?=del|ins|dup)")
          num3 <- ifelse(is.na(num3), 0, as.numeric(num3))
          if (num3 <= 2) {
            if (data_APC[i,"Start_Position"] < 112175489) {
              data_APC[i,"APC_Type"] <- "APC_N"
            } else {
              data_APC[i,"APC_Type"] <- "APC_C"
            }
          } else {
            data_APC[i,"APC_Type"] <- "other mutation"
          }
        }
      } else if (!is.na(data_APC[i, "HGVSp"])){
        num4=as.numeric(gsub("\\D", "", data_APC[i, "HGVSp_Short"]))
        if (!grepl("\\*|fs", data_APC[i, "HGVSp"])) {
          data_APC[i,"APC_Type"] <- "other mutation"
        } else if (num4<1400) {
          data_APC[i,"APC_Type"] <- "APC_N"
        } else if (num4>= 1400){
          data_APC[i,"APC_Type"] <- "APC_C"
        }
      }
    }
  }else{
    for (i in 1:nrow(data_APC)) {
      if(is.na(data_APC[i, "HGVSp"])){
        if (is.na(data_APC[i, "HGVSp_Short"])){
          data_APC[i,"APC_Type"] <- "other mutation"
        }else {
          num=as.numeric(gsub("\\D", "", data_APC[i, "HGVSp_Short"]))
          if (num >= 1 && num < 1400){
            data_APC[i,"APC_Type"] <- "APC_N"
          }else if (num>=1400){
            data_APC[i,"APC_Type"] <- "APC_C"
          }
        }
      } else if (!is.na(data_APC[i, "HGVSp"])){
        if ((!grepl("\\*|fs", data_APC[i, "HGVSp"]))){
          data_APC[i,"APC_Type"] <- "other mutation"
        }else if (grepl("*", data_APC[i, "HGVSp_Short"]) & !grepl("fs", data_APC[i, "HGVSp_Short"])){
          num1=as.numeric(gsub("\\D", "", data_APC[i, "HGVSp_Short"]))
          if (num1 >= 1 && num1 < 1400){
            data_APC[i,"APC_Type"] <- "APC_N"
          }else if (num1>=1400){
            data_APC[i,"APC_Type"] <- "APC_C"
          }
        }
        else if (grepl("*", data_APC[i, "HGVSp_Short"]) & grepl("fs", data_APC[i, "HGVSp_Short"]))
        { short_string <- sub("\\*.*", "", data_APC[i, "HGVSp_Short"])
        num2=as.numeric(gsub("\\D", "", short_string))
        if (num2 >= 1 && num2 < 1400){
          data_APC[i,"APC_Type"] <- "APC_N"
        }else if (num2>=1400){
          data_APC[i,"APC_Type"] <- "APC_C"
        }
        }
      }
    }
  }
  data_APC_Type2 <- data_APC[, c("Tumor_Sample_Barcode", "APC_Type")]
  data_APC_Type2 <- aggregate(APC_Type ~ Tumor_Sample_Barcode, data_APC_Type2, function(x) {
    if (all(x == "other mutation")) {
      "other mutation"
    } else if (all(x == "APC_C")) {
      "APC_C"
    } else if (all(x == "APC_N")) {
      "APC_N"
    } else if ("APC_N" %in% x & "other mutation" %in% x & !("APC_C" %in% x)) {
      "APC_N"
    } else if ("APC_C" %in% x & "other mutation" %in% x & !("APC_N" %in% x)) {
      "APC_C"
    } else if ("APC_N" %in% x & "APC_C" %in% x & !("other mutation" %in% x)) {
      "APC_C"
    } else if ("APC_N" %in% x & "APC_C" %in% x & "other mutation" %in% x) {
      "APC_C"
    }
  })
  data_APC_Type2 <- unique(data_APC_Type2)
  data_APC_Type1 <- data_APC_Type1[!data_APC_Type1$Tumor_Sample_Barcode %in% data_APC_Type2$Tumor_Sample_Barcode, ]
  data_APC_Type1 <- data.frame(Tumor_Sample_Barcode = data_APC_Type1)
  data_APC_Type1$APC_Type <- "Wild"
  APC <- rbind(data_APC_Type1, data_APC_Type2)
  sample_id_data <- unique(sample_id_data$Tumor_Sample_Barcode)
  new_samples <- setdiff(sample_id_data, APC$Tumor_Sample_Barcode)
  APC_new <- data.frame(Tumor_Sample_Barcode = new_samples, APC_Type = "Wild")
  APC <- rbind(APC_new, APC)
}
