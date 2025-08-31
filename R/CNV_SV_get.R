#' Retrieve CNV and SV data from the NGS registration form
#'
#' This function Retrieve CNV and SV data from the NGS registration form.
#'
#' @param data    is the data of NGS registration form
#' @param Type    A list.eg.Type=c("重排/融合","扩增")
#' @return A data frame with gene amplification and fusion information with each sample
#' @export
##' @author Pengmin Yang
#' @examples
#' data("data")
#' result<-CNV_SV_get(data=data,Type=c("重排/融合","扩增"))
#' print(result)

CNV_SV_get<-function(data,Type=c("重排/融合","扩增")){
  library(tidyr)
  library(stringr)
  data$Amplification<-NA
  results <- data.frame(Gene = character(),
                        Sample_name = character(),
                        CN1 = numeric(),
                        CN = character(),
                        stringsAsFactors = FALSE)
  for(type in Type){
    data<- data[, c(1:2)]
    A<-colnames(data)[[2]]
    for (j in 1:nrow(data)){
      if (grepl(";", data[j, A])){
        split_data <- str_split(data[j, A], ";")
        # 遍历拆分后的每个子字符串
        for(k in 1:length(split_data[[1]])){
          if (grepl(type, split_data[[1]][[k]])){
            genes <- gsub("\\(.*", "", gsub(":.*", "", split_data[[1]][[k]]))
            genes <- trimws(genes)
            data[j,genes] <- sub(".*\\((.*?)\\).*", "\\1", split_data[[1]][[k]])
          }
          else {
            data[j, "Amplification"] <- 0}
        }
      } else {
        if (grepl(type, data[j, A])){
          gene <- gsub("\\(.*", "", gsub(":.*", "", data[j, A]))
          gene <- trimws(gene)
          data[j,gene] <- sub(".*\\((.*?)\\).*", "\\1", data[j, A])
        } else {
          data[j, "Amplification"] <- 0
        }
      }
    }
    genes <- unique(colnames(data)[-(1:3)])
    data<- data[, c(1:2)]
    for (gene in genes){
      for (j in 1:nrow(data)){
        data[j,gene]<-0
        if (grepl(";", data[j,A])){
          split_data <- str_split(data[j,A], ";")
          # 遍历拆分后的每个子字符串
          for(k in 1:length(split_data[[1]])){
            if (grepl(gene, split_data[[1]][[k]]) && grepl(type, split_data[[1]][[k]])){
              value1 <- sub(".*\\((.*?)\\).*", "\\1", split_data[[1]][[k]])
              data[j,gene]<-value1
            }
          }
        }else{
          if (grepl(gene, data[j,A]) && grepl(type, data[j,A])){
            value2 <- sub(".*\\((.*?)\\).*", "\\1", data[j,A])
            data[j,gene]<-value2
          }else{
            data[j,gene]<-0
          }
        }
      }
    }
    num<-length(colnames(data))
    data_long <- tidyr::pivot_longer(data, cols = 3:num, names_to = "Gene", values_to = "CN1")
    data_filtered<- data_long[data_long$CN1 != 0, ]
    colnames(data_filtered)[colnames(data_filtered) == colnames(data_filtered)[[1]]] <-"Sample_name"
    data_filtered<- data_filtered[, c("Gene", "Sample_name","CN1")]
    if(type=="重排/融合"){
      data_filtered$CN1 <- paste0(data_filtered$CN1, "%")
      data_filtered$CN<-"Fusion"
    }else{data_filtered$CN<-"Amplification"}
    results <- rbind(results,data_filtered)
  }
  return(results)
}
