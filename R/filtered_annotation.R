#' Test annotation filtered
#'
#' This function performs filtered annotation
#'
#' @param annotations result of get_annotation
#' @param x_min ggpubr::geom_signif 参数xmin
#' @param x_max ggpubr::geom_signif 参数xmax
#' @return A list with filtered annotation,xmin and xmax
#' @export
##' @author Pengmin Yang

filtered_annotation <- function(annotations, x_min, x_max){
  valid_indices <- which(!annotations %in% c("NS.", "ns"))
  list(
    annotations = annotations[valid_indices],
    x_min = x_min[valid_indices],
    x_max = x_max[valid_indices]
  )
}
