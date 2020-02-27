#' Create Gene Level Plots
#'
#' To scale and flip contributions from a cross validation multi-omics analysis in order to be align with the contributions from the full analysis.
#' @param actorObj Output from actor model
#' @param isoData Sample data to be included in heatmap. First column must be "gene_id" and second column "feature_id"
#' @param geneList List of genes to plot.
#' @param tissueList List of tissues to plot. Must be a subset of the tissues used in `actorObj`. If left blank then the columns of `likelihoodSum` from `actorObj` will be used.
#' @param outFile PDF file path to write gene plots to. Must end in .pdf.
#' @param refPanel Reference panel user can add. Default is GTEx which is provided internally.
#' @export
#'


makeGenePlots=function(actorObj,isoData,geneList,tissueList=NULL,outFile=NULL,refPanel=NULL){
  if(is.null(refPanel)){
    alphaDirSet=gtexRefPanel

  }else{
    alphaDirSet=refPanel
  }
  if(is.null(tissueList)){
    tissueList=colnames(actorObj$posteriorEst)
  }
  if(!all(tissueList%in%colnames(actorObj$posteriorEst))){
    stop("tissueList contains tissues which are not contained in the actor object.")
  }
  if(!all(geneList%in%rownames(actorObj$posteriorEst))){
    stop("geneList contains genes which are not contained in the actor object.")
  }

  alphaDirSetRed=alphaDirSet[alphaDirSet$gene_id%in%geneList,c("gene_id","feature_id",tissueList)]
  gtexProb=data.frame("gene_id"=alphaDirSetRed$gene_id,"feature_id"=alphaDirSetRed$feature_id,alphaDirSetRed[,-(1:2)]/getGeneExpression(alphaDirSetRed[,-(1:2)],alphaDirSetRed$gene_id))
  tissuePrec=getGeneExpression2(alphaDirSetRed[,-(1:2)],alphaDirSetRed$gene_id)

  isoDataRed=isoData[isoData$gene_id%in%geneList,]
  isoProb=data.frame("gene_id"=isoDataRed$gene_id,"feature_id"=isoDataRed$feature_id, isoDataRed[,-(1:2)]/getGeneExpression(isoDataRed[,-(1:2)],isoDataRed$gene_id  ) )
  numSamp=ncol(isoProb)-2
  gtexSampMerge=merge(gtexProb,isoProb,by=c("gene_id","feature_id"),all.x=TRUE,all.y=TRUE)
  gtexSampMerge[is.na(gtexSampMerge)]=0
  if(!is.null(outFile)){

      pdf(file=outFile)


  }
  for(i in 1:length(geneList)){
    message(paste("Making Plot for ",geneList[i]))
    plotDf=gtexSampMerge[gtexSampMerge$gene_id==geneList[i],]

    annot=data.frame("Inv. Precision"=c(1/unlist(tissuePrec[geneList[i],]),rep(NA,numSamp) ),"Phi"=c( unlist(actorObj$tissuePhi[geneList[i],tissueList ]),rep(NA,numSamp) ) ,"Data"=factor(c(rep(0,length(tissueList)),rep(1,numSamp)  ) ,labels=c("Reference Panel","Experimental")))
    if(is.na(var(annot[,1],na.rm=TRUE))){
      annot=data.frame(annot[,-1])
    }
    rownames(annot)=colnames(plotDf)[-(1:2)]
    rownames(plotDf)=gtexSampMerge[gtexSampMerge$gene_id==geneList[i],2]
    p1=pheatmap(t(plotDf[,-(1:2)]) ,annotation_row=annot,main=geneList[i],breaks=seq(0,1,0.01),cluster_rows = FALSE)

    if(length(geneList)==1){

      return(p1)
    }
  }

  if(!is.null(outFile)){
    dev.off()

  }

}
