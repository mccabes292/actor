#' Create actor Object
#'
#' To scale and flip contributions from a cross validation multi-omics analysis in order to be align with the contributions from the full analysis.
#' @param countData Sample isoform data. First column should be titled "gene_id" and second "feature_id"
#' @param reduceTissue TRUE or FALSE indicating whether or not the model should reduce to the top 5 tissues.
#' @param tissueList List of tissues that will be used for fitting the model. If left NULL all tissues will be used.
#' @param refPanel Reference panel user can add. Default is GTEx which is provided internally.
#' @return directorObj - Object to be used for model fitting
#' @importFrom dplyr group_by summarize summarize_at vars
#' @importFrom pheatmap pheatmap
#' @importFrom matrixStats logSumExp
#' @importFrom magrittr %>%
#' @importFrom usedist dist_make
#' @importFrom MCMCpack rdirichlet
#' @importFrom stats cutree hclust model.matrix runif var
#' @importFrom grDevices dev.off pdf
#'
#' @export




actor=function(countData,reduceTissue=TRUE,tissueList=NULL,refPanel=NULL){
  #Value to be plugged in for NA Alpha Values in the DM
  minAlpha=0.00005
  minPhiZero=5e-324
  minLike=-2000
  tissueFilterProp=0.025
  if(!all(colnames(countData)[1:2]==c("gene_id","feature_id") )   ){
    stop("First 2 colnames are not labeled 'gene_id','feature_id'")
  }
  if(is.null(refPanel)){
    alphaDirSet=gtexRefPanel

  }else{
    alphaDirSet=refPanel
  }

  message("Filtering Genes")
  #filter down to genes in the reference panel
  isoData=countData[countData$gene_id%in%unique(alphaDirSet$gene_id),]


  #Calculate Gene Expression for the isoform data
  geneExp=getGeneExpression2(isoData[,-(1:2)],isoData$gene_id )
  fullGeneCount=function(vec){
    all(vec>10)
  }

  geneExpFull=apply(geneExp,1,fullGeneCount)

  isoData=isoData[isoData$gene_id%in%names(geneExpFull[geneExpFull==TRUE]),]
  geneExp=geneExp[unique(isoData$gene_id),]




  likeList=calcLikelihood(isoData=isoData,geneExp=geneExp,minAlpha=minAlpha,minLike = minLike,tissueFilterProp=tissueFilterProp,reduceTissue = reduceTissue,tissueList=tissueList,alphaDirSet=alphaDirSet)
  class(likeList)="actor"
  likeList$numIter=10
  likeList$numClass=10
  likeList$vbStop=0.1
  likeList$newtStop=0.1

  return(likeList)

}
