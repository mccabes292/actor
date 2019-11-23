initializeParams=function(actorObj){
  #Parameters from original model
  actorObj$classRho<-data.frame(rep(0.5,actorObj$numClass))
  rownames(actorObj$classRho)<-paste("geneClass",1:actorObj$numClass,sep="")
  actorObj$classOmega<-actorObj$classRho
  actorObj$tissueBeta=matrix(rep(3,actorObj$numTissue),ncol=1)
  rownames(actorObj$tissueBeta)=colnames(actorObj$likelihoodSum)
  colnames(actorObj$tissueBeta)="beta"
  #Parameters from VB Model
  #classLambda=matrix(rep(1/numClass,numClass*numGenes),ncol=numClass)
  actorObj$classLambda=rdirichlet(actorObj$numGenes,rep(0.25,actorObj$numClass))
  colnames(actorObj$classLambda)=paste("geneClass",1:actorObj$numClass,sep="")
  rownames(actorObj$classLambda)=rownames(actorObj$likelihoodSum)
  logSumPhi=apply(actorObj$likelihoodSum,1,logSumExp)
  actorObj$tissuePhi=exp(actorObj$likelihoodSum-logSumPhi)
  rownames(actorObj$tissuePhi)=rownames(actorObj$likelihoodSum)
  colnames(actorObj$tissuePhi)=colnames(actorObj$likelihoodSum)



  actorObj$dirEta=matrix(runif(actorObj$numTissue*actorObj$numClass,0,2),ncol=actorObj$numClass)

  rownames(actorObj$dirEta)=colnames(actorObj$likelihoodSum)
  colnames(actorObj$dirEta)=rownames(actorObj$classOmega)


  return(actorObj)
}

#####Update Functions for model fitting
#Phi
updatePhiFunction=function(classLambda,dirEta,likelihoodSum){
  update1=(likelihoodSum+(classLambda%*%(t(digamma(dirEta))-digamma(colSums(dirEta)))))
  logSum=apply(update1,1,logSumExp)

  return(exp(update1-logSum))
}


#Lambda
updateLambdaFunction=function(tissuePhi,dirEta,classRho){
  diGamRho=as.vector(digamma(classRho[,1])-digamma(sum(classRho[,1])))
  t1=exp(t((t(digamma(dirEta))-digamma(colSums(dirEta)))%*%t(tissuePhi)+diGamRho  ))
  rs=rowSums(t1)
  return(t1/rs)
}

#Eta Update
updateEtaFunction=function(tissuePhi,classLambda,tissueBeta){
  return((t(tissuePhi)%*%classLambda+tissueBeta[,1]))
}



####Model Parameters
#Gamma
updateGammaFunction=function(classLambda){
  numGenes=nrow(classLambda)
  cs=colSums(classLambda)
  return(cs/numGenes)
}


updateRhoFunction=function(classLambda,classOmega){
  return(colSums(classLambda)+classOmega)
}






updateOmegaFunction=function(classOmega,classRho,newtStop){
  numClass=nrow(classOmega)
  diGamRhoSum=sum(as.vector(digamma(classRho[,1])-digamma(sum(classRho[,1]))))
  omegaDiff=10
  omegaOld=(classOmega[1,1])
  while(omegaDiff>newtStop){
    score=numClass*digamma(numClass*omegaOld)-numClass*digamma(omegaOld)+diGamRhoSum
    info=numClass*numClass*trigamma(numClass*omegaOld)-numClass*trigamma(omegaOld)
    omegaNew=(log(omegaOld)-1/omegaOld*score/info   )
    omegaDiff=(omegaNew-log(omegaOld))^2
    omegaOld=exp(omegaNew)

  }
  return(omegaOld)
}

##################Check on this function. Exiting if info is NA
updateBetaFunction=function(tissueBeta,dirEta,newtStop){
  numClass=ncol(dirEta)
  numTissue=nrow(tissueBeta)
  diGamDiffSum=sum(t(digamma(dirEta))-digamma(colSums(dirEta)))
  betaOld=tissueBeta[1,1]
  betaDiff=10
  while(betaDiff>newtStop){
    score=numClass*numTissue*digamma(numTissue*betaOld)-numClass*numTissue*digamma(betaOld)+diGamDiffSum
    info=numClass*numTissue*numTissue*trigamma(numTissue*betaOld)-numClass*numTissue*trigamma(betaOld)
    if(is.na(info)){
      warning("betaFail")
      return(5e-2)
    }
    betaNew=log(betaOld)-1/betaOld*score/info
    betaDiff=(betaNew-log(betaOld))^2
    betaOld=exp(betaNew)
  }
  return(betaOld)

}

