#' Creates a 'InputHeatmap' object from 'tbl_df' on evaluation creates a 'ComplexHeatmap'
#'
#' heatmap() takes a tbl object and easily produces a ComplexHeatmap plot, with integration with tibble and dplyr frameworks.
#'
#' @param input data frame with 'ID', 'group' and features
#' @param ID This parameter specifies the column name in the input data that represents the unique identifier. The default value is "ID".
#' @param features 	This parameter is a vector of column names representing the features to be used in the heatmap.
#' @param group This parameter is a column name in the input data that represents the grouping variable.
#' @param condiction This parameter is an optional data frame that defines additional conditions to be applied in the heatmap. It should have two columns: one specifying the variables and another specifying the conditions. The default value is NULL.
#' @param id_condiction This parameter specifies the column name in the condiction data frame that represents the variables. The default value is "vars".
#' @param col_condiction This parameter specifies the column name in the condiction data frame that represents the conditions. The default value is "condiction".
#' @param scale This parameter specifies whether to scale the heatmap by column. If set to TRUE, the heatmap will be scaled by column. The default value is FALSE.
#' @param palette This parameter specifies the palette to be used in the heatmap. It can be a numeric value representing the palette index or a string specifying the palette name. The default value is 2.
#' @param palette_group This parameter specifies the palette group to be used in the heatmap. It should be a string representing the palette group name. The default value is "jama".
#' @param heatmap_col defalut heatmap_col=NULL. tidyHeatmap::heatmap参数palette_value
#' @param show_col This parameter specifies whether to show the color boxes in the heatmap. If set to TRUE, the color boxes will be shown. The default value is FALSE.
#' @param show_palettes This parameter specifies whether to show the available palettes. If set to TRUE, the available palettes will be shown. The default value is FALSE.
#' @param cols_group This parameter is a vector of colors to be used for the grouping cols variable. If not specified, the colors will be randomly assigned.
#' @param row_group This parameter is a vector of colors to be used for the grouping row variable. If not specified, the colors will be randomly assigned.
#' @param legend_show This parameter specifies whether to show the annotaion legend.The default value is TRUE.
#' @param show_plot This parameter specifies whether to show the resulting heatmap plot. If set to TRUE, the plot will be shown. The default value is TRUE.
#' @param width This parameter specifies the width of the resulting heatmap plot in inches. The default value is 8.
#' @param height This parameter specifies the height of the resulting heatmap plot in inches. If not specified, the height will be calculated based on the number of features.
#' @param size_col This parameter specifies the font size of the column names in the heatmap plot. The default value is 10.
#' @param size_row This parameter specifies the font size of the row names in the heatmap plot. The default value is 8.
#' @param angle_col This parameter specifies the rotation angle of the column names in the heatmap plot. The default value is 90.
#' @param column_title This parameter specifies the title of the column in the heatmap plot. If not specified, no title will be displayed.
#' @param column_title_size This parameter specifies the title size of the column in the heatmap plot.
#' @param row_title This parameter specifies the title of the row in the heatmap plot. If not specified, no title will be displayed.
#' @param row_title_size This parameter specifies the title size of the row in the heatmap plot.
#' @param show_heatmap_col_name This parameter specifies whether to show the column names in the heatmap plot. If set to TRUE, the column names will be shown. The default value is FALSE.
#' @param name This parameter is the title of heatmap_col.
#'
#' @return A heatmap plot object.
#' @export
##' @author Pengmin Yang
#' @examples
#' data("input1")
#' data("condiction")
#'
#' feas1<-colnames(input1)[3:102]
#' sig_Heatmap(input = input1, features = feas1,ID ="SAMPLE_ID",show_plot=F,
#'             condiction=condiction,id_condiction=colnames(condiction)[[1]],col_condiction=colnames(condiction)[[2]],
#'              cols_group=c("#757575","#FF4040"),row_group=c("red","green"),
#'              legend_show=TRUE,column_title_size=10,row_title_size=8,
#'              heatmap_col=NULL,
#'              #heatmap_col=c("#0505FA", "#FFFFFF", "#FA050D"),
#'              group = "PREX2",row_title="Regulate", scale = TRUE,name="Expression")
#'
sig_Heatmap<-function (input, ID = "ID", features, group, condiction = NULL,
                       id_condiction = "vars", col_condiction = "condiction", scale = FALSE,
                       palette = 2, palette_group = "jama",heatmap_col, show_col = F, show_palettes = F,
                       cols_group = NULL, row_group= NULL,legend_show=TRUE,show_plot = T, width = 8, height = NULL,
                       size_col = 10, size_row = 8, angle_col = 90, column_title = NULL,column_title_size=13.2,
                       row_title = NULL,row_title_size=13.2, show_heatmap_col_name = F, name)
{
  input <- as.data.frame(input)
  features <- features[features %in% colnames(input)]
  colnames(input)[which(colnames(input) == ID)] <- "idd"
  #input <- input[, c("idd", group, features)]
  colnames(input)[which(colnames(input) == group)] <- "target_group"
  input <- input[!is.na(input[, "target_group"]), ]
  num<-as.numeric(length(colnames(input)[!colnames(input) %in% features]))+1
  pf_long_group <- tidyr::pivot_longer(input, num:ncol(input),
                                       names_to = "variables", values_to = "value")
  if (!is.null(condiction)) {
    #condiction <- condiction[, c(id_condiction, col_condiction)]
    colnames(condiction)[which(colnames(condiction) == id_condiction)] <- "vars"
    colnames(condiction)[which(colnames(condiction) == col_condiction)] <- "condiction"
    pf_long_group <- merge(pf_long_group, condiction, by.x = "variables",
                           by.y = "vars", all.x = TRUE, all.y = FALSE)
    pf_long_group$condiction <- ifelse(is.na(pf_long_group$condiction),
                                       "NE", pf_long_group$condiction)
    print(head(pf_long_group))
  }
  if (is.null(height)) {
    height_heatmap <- length(features) * 0.1 + 3
  }
  else {
    height_heatmap <- height
  }
  if(is.null(heatmap_col)){
    heatmap_col <- IOBR::palettes(category = "tidyheatmap", palette = palette,show_col = show_col,
                                  show_message = show_palettes)
  }
  else{
    heatmap_col<-heatmap_col
  }
  if (!is.null(cols_group)) {
    color_box <- cols_group
  }
  else {
    if (is.null(palette)) {
      color_box <- IOBR::palettes(category = "random", show_col = show_col,
                                  show_message = show_palettes)
    }
    else {
      color_box <- IOBR::palettes(category = "box", palette = palette_group,
                                  show_col = show_col, show_message = show_palettes)
    }
  }
  if (!is.null(condiction)) {
    target_level1 <- unique(as.character(pf_long_group$condiction))
    target_level1 <- target_level1[!is.na(target_level1)]
    n <- length(target_level1)
    color_box1 <- color_box[1:n]
  }
  target_level2 <- unique(as.character(pf_long_group$target_group))
  target_level2 <- target_level2[!is.na(target_level2)]
  n <- length(target_level2)
  color_box2 <- color_box[1:n]
  if (scale) {
    scale = "row"
  }
  else {
    pf_long_group$value[pf_long_group$value > 2.5] = 2.5
    pf_long_group$value[pf_long_group$value < -2.5] = -2.5
    scale = "none"
  }
  if (is.null(condiction)) {
    pp <- pf_long_group %>% dplyr::group_by(target_group) %>%
      tidyHeatmap::heatmap(.column = idd, .row = variables,
                           .value = value, palette_grouping = list(c(cols_group)),
                           scale = scale, column_title = column_title,
                           column_title_gp = gpar(fontsize = column_title_size),
                           row_title = row_title, row_title_gp = gpar(fontsize = row_title_size),
                           palette_value = heatmap_col, show_column_names = show_heatmap_col_name,
                           column_names_gp = grid::gpar(fontsize = size_col),
                           row_names_gp = grid::gpar(fontsize = size_row),
                           column_names_rot = angle_col,heatmap_legend_param = list(title = name))
    if(legend_show){
      pp<-pp%>%annotation_tile(target_group,palette=cols_group,size = unit(0.3, "cm"),
                               show_annotation_name = T,annotation_name_align = T,
                               annotation_name_rot = 0,annotation_name_gp = gpar(fontsize =column_title_size),
                               annotation_name_side="right")
      pp@group_top_annotation<-list()
    }
  }
  else {
    pp <- pf_long_group %>%
      dplyr::group_by(target_group,condiction) %>%
      tidyHeatmap::heatmap(.column = idd, .row = variables, .value = value,
                           palette_grouping = list(c(row_group),c(cols_group)), scale = scale,
                           column_title = NULL, column_title_gp = gpar(fontsize = column_title_size),
                           row_title = NULL, row_title_gp = gpar(fontsize = row_title_size),
                           palette_value = heatmap_col,show_column_names = show_heatmap_col_name,
                           column_names_gp = grid::gpar(fontsize = size_col),
                           row_names_gp = grid::gpar(fontsize = size_row), column_names_rot = angle_col,
                           heatmap_legend_param = list(title = name))
    if(legend_show){
      pp<-pp%>%annotation_tile(target_group,palette=cols_group,size = unit(0.3, "cm"),
                               show_annotation_name = T,annotation_name_align = T,
                               annotation_name_rot = 0,annotation_name_gp = gpar(fontsize =column_title_size),
                               annotation_name_side="right")%>%
        annotation_tile(condiction,palette=row_group,size = unit(0.3, "cm"),
                        show_annotation_name = T,annotation_name_align = T,
                        annotation_name_rot = 270,annotation_name_gp = gpar(fontsize =row_title_size),
                        annotation_name_side="top")
      pp@group_top_annotation<-list()
      pp@group_left_annotation<-list()
    }
  }
  if(is.null(condiction)){
    colnames(pp@data)[which(colnames(pp@data) == "target_group")]<- paste(group,"Group")
    pp@top_annotation[["col_name"]]<- paste(group,"Group")
  }else{
    colnames(pp@data)[which(colnames(pp@data) == "target_group")]<- paste(group,"Group")
    colnames(pp@data)[which(colnames(pp@data) == "condiction")]<- row_title
    pp@top_annotation[["col_name"]]<- paste(group,"Group")
    pp@left_annotation[["col_name"]]<- row_title
  }
  if (show_plot)
    print(pp)
  return(pp)
}
