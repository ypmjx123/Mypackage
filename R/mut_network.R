#' Gene Mutation Co_Occurence/Mutually_Exclusive analysis and visualization
#'
#' @param SNV A data frame containing SNV data.即MAF格式的数据
#' @param top 突变率前x的基因, top=10
#' @param pValue pValue
#' @param customdata 基因注释数据,第一列为基因名，第二列为注释信息参考data("Gene",package = "Mypackage")，默认为NULL
#'
#' @returns A list of table of Gene Mutation Co_Occurence/Mutually_Exclusive analysis and the network plot
#' @export
#'
##' @author Pengmin Yang
##'
#' @examples
#' data("mutation_CRC")
#' result<-mut_network(SNV=mutation_CRC,
#'                     top=20,pValue=0.01,customdata=NULL)
#'
mut_network<-function(SNV,top,pValue,customdata=NULL){
  library(igraph)
  library(ggraph)
  library(maftools)
  library(tidyr)
  library(dplyr)
  maf<-read.maf(maf=SNV)
  p<-maftools::somaticInteractions(maf = maf, top = top)
  edges <- p[, c(1,2,3,10)]
  pValue <- as.numeric(pValue)
  edges <- edges[edges$pValue < pValue,]
  edges_list <- as.data.frame(edges)
  nodes <- data.frame(name = unique(c(edges_list$gene1,edges_list$gene2)))
  for(gene in nodes$name){
    count<-length(unique(maf@data$Tumor_Sample_Barcode[maf@data$Hugo_Symbol==gene]))
    vaf<-length(unique(maf@data$Tumor_Sample_Barcode[maf@data$Hugo_Symbol==gene]))/length(unique(maf@data$Tumor_Sample_Barcode))
    nodes[nodes$name==gene,"count"]<-count
    nodes[nodes$name==gene,"vaf"]<-vaf
  }
  if(is.null(customdata)){
    data("Gene",package = "Mypackage")
  }else{
    Gene<-customdata
    colnames(Gene)[1]<-"Hugo Symbol"
    colnames(Gene)[2]<-"Gene Type"
  }
  nodes1<-merge(nodes,Gene,by.x="name",by.y="Hugo Symbol",all.x=T)
  nodes1$`Gene Type`<-ifelse(is.na(nodes1$`Gene Type`),"Unknown",nodes1$`Gene Type`)
  nodes1 <- nodes1 %>%
    distinct(name, .keep_all = TRUE)
  g <- graph_from_data_frame(edges_list, vertices = nodes1, directed = FALSE)
  p1<-ggraph(g, layout = 'igraph', algorithm = 'kk') +
    geom_node_point(aes(size=vaf,color=`Gene Type`)) +
    geom_node_text(aes(label = name), vjust = 1.5) +
    geom_edge_link(aes(color = Event), alpha = 0.8) +
    theme_void()+
    theme(legend.position ="right")+
    geom_edge_density(aes(fill = Event),show.legend=F)
  result=list(network=p,
              plot=p1)
  return(result)
}
