

calcELBOFull=function(actorObj){
  diGamDiff=(t(digamma(actorObj$dirEta))-digamma(colSums(actorObj$dirEta)))
  diGamRho=digamma(actorObj$classRho[,1])-digamma(sum(actorObj$classRho[,1]))
  diGamOmega=digamma(actorObj$classOmega[,1])-digamma(sum(actorObj$classOmega[,1]))
  eX=sum(actorObj$tissuePhi*actorObj$likelihoodSum,na.rm=TRUE)
  eTC=sum(diag((diGamDiff%*%t(actorObj$tissuePhi)+diGamRho)%*%actorObj$classLambda))
  eGamma=lgamma(sum(actorObj$classOmega[,1]))-sum(lgamma(actorObj$classOmega))+sum((actorObj$classOmega[,1]-1)*diGamRho)



  eTheta=actorObj$numClass*(lgamma(sum(actorObj$tissueBeta))-sum(lgamma(actorObj$tissueBeta)))+sum(diGamDiff%*%(actorObj$tissueBeta-1))
  eqTheta=sum(lgamma(colSums(actorObj$dirEta)),-lgamma(actorObj$dirEta))+sum(diag(diGamDiff%*%(actorObj$dirEta-1)))
  eqTC=sum(actorObj$tissuePhi*log(actorObj$tissuePhi),actorObj$classLambda*log(actorObj$classLambda),na.rm=TRUE)
  eqGamma=lgamma(sum(actorObj$classRho[,1]))-sum(lgamma(actorObj$classRho)) + sum((actorObj$classRho[,1]-1)*diGamRho )
  return(eX+eTC+eTheta+eGamma-eqTheta-eqTC-eqGamma)
}








