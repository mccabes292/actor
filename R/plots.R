#' Create plots for visualization of output
#'
#'
#'
#' @param x Object of class `actorFit`
#' @param plotType Specifies the type of plot to be generated
#' \code{"TissueMem"} - Plots posterior probabilities for the tissue membership  \cr\cr
#' \code{"ClassMem"} - Plots posterior probabilities for gene class membership \cr\cr
#' \code{"GeneClass"} - Plots gene class posterior estimates of each tissue.
#' @param annotMat_Col Dataframe for annotation bars to be added to the columns of the heatmap
#' @param annotMat_Row Dataframe for annotation bars to be added to the rows of the heatmap
#' @param fontSize Font size for text in plot
#' @param addAnnot TRUE or FALSE for adding annotation bars to the heatmap.
#' @param ... other parameters to be passed through to plotting functions.
#' @export
#'


plot.actorFit=function(x,plotType,addAnnot=FALSE, annotMat_Col=NULL,annotMat_Row=NULL,fontSize=10,...){
  if(!(plotType%in%c("TissueMem","ClassMem","GeneClass"))){
    stop("Not a valid plotType")
  }
  # if(addAnnot==FALSE & (!is.na(annotMat_Col) | !is.na(annotMat_Row)  )){
  #   addAnnot=TRUE
  # }


  if(plotType=="TissueMem"){



    distMat=dist_make((x$tissuePhi),distFunc)

    gc=factor(apply(x$classLambda,1,which.max))
    if(is.null(annotMat_Col)){
      annotCol=data.frame("Gene_Class_Model"=gc)
    }else{
      if(!all(rownames(annotMat_Col)==names(gc))){
        warning("Rownames of annotMat_Col do not match rownames of model object.")
      }
      annotCol=data.frame('Gene_Class_Model'=gc,annotMat_Col)
    }
    if(addAnnot==TRUE){
      p1=pheatmap(t(x$tissuePhi),cluster_cols = hclust(distMat),fontsize = fontSize,annotation_col = annotCol,annotation_row=annotMat_Row,show_colnames = FALSE)
      return(p1)
    }
    p1=pheatmap(t(x$tissuePhi),cluster_cols = hclust(distMat),fontsize = fontSize,show_colnames = FALSE)
    return(p1)

  }
  if(plotType=="ClassMem"){
    gc=factor(apply(x$tissuePhi,1,which.max))
    if(is.null(annotMat_Col)){
      annotCol=data.frame("Tissue"=gc)
    }else{
      annotCol=data.frame('Tissue'=gc,annotMat_Col)
    }


    distMat=dist_make((x$classLambda),distFunc)
    if(addAnnot==TRUE){
      p1=pheatmap(t(x$classLambda),cluster_cols = hclust(distMat),fontsize = fontSize,annotation_col = annotCol,annotation_row=annotMat_Row,show_colnames = FALSE)

      return(p1)
    }
    p1=pheatmap(t(x$classLambda),cluster_cols = hclust(distMat),fontsize = fontSize,show_colnames = FALSE)

  }
  if(plotType=="GeneClass"){
    classProb=x$classRho/sum(x$classRho)
    if(is.null(annotMat_Col)){
      annotCol=data.frame("Class_Probability"=classProb)
    }else{
      annotCol=data.frame("Class_Probability"=classProb,annotMat_Col)
    }
    if(addAnnot==TRUE){
      p1=pheatmap(x$dirEta,cluster_cols = FALSE,cluster_rows = FALSE,fontsize = fontSize,annotation_col = annotCol,annotation_row=annotMat_Row,show_colnames = FALSE)
      return(p1)
    }
    p1=pheatmap(x$dirEta,cluster_cols = FALSE,cluster_rows = FALSE,fontsize = fontSize,show_colnames = FALSE)
    return(p1)
}



}
