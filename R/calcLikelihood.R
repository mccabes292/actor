#Likelihood Function


calcLikelihood=function(isoData,geneExp,minAlpha,minLike,tissueFilterProp,reduceTissue=TRUE,tissueList=NULL){
  if(is.null(tissueList)){
    keepTissue=colnames(alphaDirSet)[-(1:2)]
  }else{
    keepTissue=tissueList
  }
  alphaDirSet=alphaDirSet[alphaDirSet$gene_id%in%rownames(geneExp),colnames(alphaDirSet)%in%c("gene_id","feature_id",keepTissue)]

  likeFilter=minLike
  mergeSet=data.frame(merge(alphaDirSet,isoData[,(1:2)],by=c("gene_id","feature_id"),all.x=TRUE,all.y=TRUE))
  mergeSet[is.na(mergeSet)]=minAlpha


  ##Remove Genes which are not present at all in GTEx
  rowMax=data.frame("gene_id"=mergeSet$gene_id, "feature_id"=mergeSet$feature_id,"rowMax"=apply(mergeSet[,-(1:2)],1,max))
  maxRowMax=rowMax%>%group_by(gene_id)%>%summarize(max(rowMax))
  keepGene=maxRowMax[!(maxRowMax[,2]<=minAlpha),1]
  mergeSetFull=mergeSet[mergeSet$gene_id%in%unlist(keepGene),]
  isoData=merge(isoData,mergeSetFull[,1:2],by=c("gene_id","feature_id"),all.y=TRUE)
  isoData[is.na(isoData)]=0
  geneCount=table(isoData$gene_id)


  geneExp=geneExp[unique(isoData$gene_id),]

  # #Exclude Genes based on an average log count
  # avgLogCt=log(apply(geneExp,1,mean)+1)
  # numGenes=nrow(geneExp)
  # #numKeepExpr=1000
  # numKeepExpr=numGenes
  # keepGenes=names(avgLogCt[order(avgLogCt,decreasing=TRUE)[1:numKeepExpr]])
  #
  # isoData=isoData[isoData$gene_id%in%keepGenes,]
  # geneExp=data.matrix(geneExp[unique(isoData$gene_id),])


  #Subset dirichlet matrix to final set of genes
  alphaDirSet=mergeSetFull[match(isoData$feature_id,mergeSetFull$feature_id),]

  tissuePrec=getGeneExpression(alphaDirSet[,-(1:2)],alphaDirSet$gene_id)
  tissueEst=data.frame(alphaDirSet[,1:2],alphaDirSet[,-(1:2)]/tissuePrec)
  tissueEst[alphaDirSet<=minAlpha]=NA

  estVar=data.frame("gene_id"=tissueEst$gene_id,"varRow"=apply(tissueEst[,-(1:2)],1,var,na.rm=TRUE))
  estVar%>%group_by(gene_id)%>%summarize(varSum=sum(varRow,na.rm=TRUE))->geneVar
  keepGeneFinal=geneVar$gene_id[geneVar$varSum>=0.001]
  tissueEst=tissueEst[tissueEst$gene_id%in%keepGeneFinal,]
  alphaDirSet=alphaDirSet[alphaDirSet$gene_id%in%keepGeneFinal,]
  tissuePrec=getGeneExpression(alphaDirSet[,-(1:2)],alphaDirSet$gene_id)
  isoData=isoData[isoData$gene_id%in%keepGeneFinal,]
  geneExp=geneExp[keepGeneFinal,]


  t1=alphaDirSet[,-(1:2)]
  t1[tissuePrec>50]=t1[tissuePrec>50]*50/(tissuePrec[tissuePrec>50])
  alphaDirSet[,-(1:2)]=t1




  maxAlphaDir=alphaDirSet%>%group_by(gene_id)%>%summarize_at(max,.vars=vars(-feature_id))
  maxAlphaDir2=data.matrix(maxAlphaDir[,-1])
  rownames(maxAlphaDir2)=maxAlphaDir$gene_id
  colnames(maxAlphaDir2)=colnames(alphaDirSet)[-(1:2)]
  ############Currently matched isoforms with data but need to adjust if not#####
  numGenes=nrow(geneExp)
  numSamp=ncol(geneExp)
  numTissue=ncol(alphaDirSet)-2
  #Pre Calculate values not dependent on model parameters
  sumDirichlet=rowsum(alphaDirSet[,-(1:2)],alphaDirSet$gene_id)
  #Likelihood Function
  lgSumDir=lgamma(sumDirichlet)
  geneLFact=lfactorial(geneExp)
  lg3=rowsum(lgamma(alphaDirSet[,-(1:2)]),alphaDirSet$gene_id  )
  isoLFact=rowsum(lfactorial(isoData[,-(1:2)]),isoData$gene_id)

  likelihoodSum=matrix(rep(0,numGenes*numTissue),ncol=numTissue)
  likelihoodSamp=array(rep(NA,numGenes*numTissue*numSamp),dim=c(numGenes,numTissue,numSamp),dimnames = list(rownames(geneExp),colnames(alphaDirSet)[-(1:2)],colnames(geneExp)))

  posteriorList=array(rep(NA,numGenes*numTissue*numSamp),dim=c(numGenes,numTissue,numSamp),dimnames = list(rownames(geneExp),colnames(alphaDirSet)[-(1:2)],colnames(geneExp)))

  message("Calculating Likelihood")
  for(i in 1:numSamp){
    lg1=lgamma(geneExp[,i]+sumDirichlet )
    lg2=rowsum(lgamma(unlist(isoData[,i+2])+alphaDirSet[,-(1:2)] ),alphaDirSet$gene_id)
    likelihoodSamp[,,i]=data.matrix(lgSumDir-lg1+lg2-lg3+geneLFact[,i]-isoLFact[,i])
    likelihoodSamp[,,i][maxAlphaDir2<=minAlpha]<-minLike
    temp1=apply(likelihoodSamp[,,i],1,logSumExp)
    posteriorList[,,i]=data.matrix(exp(likelihoodSamp[,,i]-temp1))
    rownames(posteriorList[,,i])=rownames(likelihoodSum)
    colnames(posteriorList[,,i])=colnames(likelihoodSum)
    likelihoodSum=likelihoodSum+likelihoodSamp[,,i]

  }

  logSumLike=apply(likelihoodSum,1,logSumExp)
  posteriorEst=exp(likelihoodSum-logSumLike)
  colnames(posteriorEst)=colnames(likelihoodSum)




  if(reduceTissue==TRUE & ncol(likelihoodSum)>5){
    message("Reducing to top 5 tissues")
    ###Test Function
    flag=FALSE
    i=0
    estTemp=tissueEst[tissueEst$gene_id%in%rownames(likelihoodSum),c("gene_id","feature_id",colnames(likelihoodSum))]
    while(flag==FALSE){
      print(i)
      #Remove genes where max falls below our threshold
      likeMax=apply(likelihoodSum,1,max)
      likelihoodSum=likelihoodSum[likeMax>likeFilter*numSamp,]
      posteriorEst=posteriorEst[likeMax>likeFilter*numSamp,]

      #Remove least similar tissue if below threshold
      postTissueSum=apply(posteriorEst,2,sum)
      tissFiltNum=nrow(likelihoodSum)*tissueFilterProp
      likelihoodSum=likelihoodSum[,-which.min(postTissueSum)]
      logSumLike=apply(likelihoodSum,1,logSumExp)
      posteriorEst=exp(likelihoodSum-logSumLike)
      postTissueSum=apply(posteriorEst,2,sum)

      tissFiltNum=nrow(likelihoodSum)*tissueFilterProp
      likeMax=apply(likelihoodSum,1,max)
      i=i+1

      #Remove genes where max falls below our threshold
      likeMax=apply(likelihoodSum,1,max)
      estTemp=estTemp[estTemp$gene_id%in%rownames(likelihoodSum),c("gene_id","feature_id",colnames(likelihoodSum))]
      rm=data.frame("gene_id"=estTemp$gene_id,"feature_id"=estTemp$feature_id,"maxVal"=apply(estTemp[,-(1:2)],1,max,na.rm=TRUE),"varRow"=apply(estTemp[,-(1:2)],1,var,na.rm=TRUE))
      rm%>%group_by(gene_id)%>%summarize("secondVal"=maxVal[order(maxVal,decreasing=TRUE)[2]],"varSum"=sum(varRow,na.rm=TRUE))->secondVal
      secondValKeep=secondVal[secondVal$secondVal>=0.1,]
      varSumKeep=secondValKeep[secondValKeep$varSum>=0.001,]
      likelihoodSum=likelihoodSum[(likeMax>likeFilter*numSamp)&(rownames(likelihoodSum)%in%varSumKeep$gene_id),]
      posteriorEst=posteriorEst[rownames(likelihoodSum),]
      flag=(ncol(likelihoodSum)==5)
    }
  }
  #Initialize Parameter Values
  numGenes=nrow(likelihoodSum)
  numSamp=ncol(isoData)-2
  numTissue=ncol(likelihoodSum)

  returnList=list("likelihoodSum"=likelihoodSum,"posteriorEst"=posteriorEst,"numGenes"=numGenes,"numTissue"=ncol(likelihoodSum),"numSamp"=numSamp)

  return(returnList)
}







