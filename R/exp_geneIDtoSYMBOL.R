#' Epression data gene ID to SYMBOL
#'
#' This function convert the gene name of expression data to SYMBOL.
#'
#' @param exp          is the data of gene expression for samples
#' @param genecoltype  exp gene的类型。 c("ENTREZID","ENSEMBL","MAP");SYMBOL:"CEBPG",ENTREZID:"1054",ENSEMBL:"ENSG00000153879",MAP:"19q13.11"
#' @export
##' @author Pengmin Yang
#' @examples
#' data("exp_raw")
#' result<-exp_geneIDtoSYMBOL(exp=exp_raw,genecoltype="ENTREZID")
#' print(result[["gene"]][1:10,1:2])
#' print(result[["data"]][1:10,1:5])
exp_geneIDtoSYMBOL<-function(exp,genecoltype="ENTREZID"){
  if(genecoltype=="ENTREZID"){
    A <- clusterProfiler::bitr(exp[[colnames(exp)[[1]]]],fromType = 'ENTREZID',toType = 'SYMBOL',OrgDb = "org.Hs.eg.db")
    exp<-merge(A,exp,by.x= 'ENTREZID',by.y=colnames(exp)[[1]],all.x=T)
    exp<-exp[!is.na(exp[[colnames(exp)[[3]]]]),]
    if (any(duplicated(exp[[1]]))) {
      exp<-exp[,-1]
      exp <- aggregate(as.formula(paste0(".", " ~ ", colnames(exp)[[1]])), data = exp, FUN = mean)
      exp<-merge(A,exp,by.x= 'SYMBOL',by.y=colnames(exp)[[1]],all.x=T)
    }
  }else if(genecoltype=="ENSEMBL"){
    A <- clusterProfiler::bitr(exp[[colnames(exp)[[1]]]],fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = "org.Hs.eg.db")
    exp<-merge(A,exp,by.x= 'ENSEMBL',by.y=colnames(exp)[[1]],all.x=T)
    exp<-exp[!is.na(exp[[colnames(exp)[[3]]]]),]
    if (any(duplicated(exp[[1]]))) {
      exp<-exp[,-1]
      exp <- aggregate(as.formula(paste0(".", " ~ ", colnames(exp)[[1]])), data = exp, FUN = mean)
      exp<-merge(A,exp,by.x= 'SYMBOL',by.y=colnames(exp)[[1]],all.x=T)
    }
  }else if(genecoltype=="MAP"){
    A <- clusterProfiler::bitr(exp[[colnames(exp)[[1]]]],fromType = 'MAP',toType = 'SYMBOL',OrgDb = "org.Hs.eg.db")
    exp<-merge(A,exp,by.x= 'MAP',by.y=colnames(exp)[[1]],all.x=T)
    exp<-exp[!is.na(exp[[colnames(exp)[[3]]]]),]
    if (any(duplicated(exp[[1]]))) {
      exp<-exp[,-1]
      exp <- aggregate(as.formula(paste0(".", " ~ ", colnames(exp)[[1]])), data = exp, FUN = mean)
      exp<-merge(A,exp,by.x= 'SYMBOL',by.y=colnames(exp)[[1]],all.x=T)
    }
  }
  results<-list(gene=A,
                data=exp)
  return(results)
}
