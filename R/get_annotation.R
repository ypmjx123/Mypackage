#' Test annotation get
#'
#' This function performs statistical difference testing between every two groups
#'
#' @param data 长数据格式，最后一列为值
#' @param group 分组的名称
#' @param types A character vector indicating the type.例如c("Wild", "Mut")
#' @param test 参考ggpubr::compare_means，根据需要选择
#' @return A list with test annotation
#' @export
##' @author Pengmin Yang
get_annotation <- function(data, group,types,test="wilcox.test") {
  annotations <- c()
  data_subset_filtered <- data[data[[group]] %in% types, ]#
  num<-length(colnames(data_subset_filtered))
  test_result <- ggpubr::compare_means(as.formula(paste0(colnames(data[,num]) , " ~ ", group)), data_subset_filtered, #
                                       method = test)
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
