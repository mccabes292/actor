#' Fit model
#'
#' To scale and flip contributions from a cross validation multi-omics analysis in order to be align with the contributions from the full analysis.
#' @param actorObj actor object from `actor`. Model specifications can be set within the object.
#' @return Object of class `actorFit` which can be used for post model fitting interpretations and analyses.
#'
#' @export



actorFit=function(actorObj){
  distMat=dist_make((actorObj$posteriorEst),distFunc)
  hc1=hclust(distMat)
  start=proc.time()
  actorObjList=NULL
  likeVec=rep(NA,actorObj$numIter)
  for(j in 1:actorObj$numIter){
    message(paste("Beggining iteration",j))
    tempObj=initializeParams(actorObj)
    if(tempObj$numClass>1){
      treeSplit=cutree(hc1,k=actorObj$numClass)
      tempObj$classLambda=model.matrix(~factor(treeSplit)-1)
      colnames(tempObj$classLambda)=paste("geneClass",1:tempObj$numClass,sep="")
      rownames(tempObj$classLambda)=rownames(tempObj$tissuePhi)
    }


    tempObj$tissueBeta=data.matrix(rep(updateBetaFunction(tempObj$tissueBeta,tempObj$dirEta,tempObj$newtStop),tempObj$numTissue))

    if(tempObj$numClass>1){
      tempObj$classRho=updateRhoFunction(tempObj$classLambda,tempObj$classOmega)
      tempObj$classOmega=data.matrix(rep(updateOmegaFunction(tempObj$classOmega,tempObj$classRho,tempObj$newtStop ),tempObj$numClass))


    }
    tempObj$dirEta=updateEtaFunction(tempObj$tissuePhi,tempObj$classLambda,tempObj$tissueBeta)
    elboFullOld=10
    elboFullDiff=10
    while(elboFullDiff>tempObj$vbStop){
      elboDiff=10
      elboOld=calcELBOFull(tempObj)
      while(elboDiff>tempObj$vbStop){
        tempObj$tissueBeta=data.matrix(rep(updateBetaFunction(tempObj$tissueBeta,tempObj$dirEta,tempObj$newtStop),tempObj$numTissue))
        if(tempObj$numClass>1){
          tempObj$classOmega=data.matrix(rep(updateOmegaFunction(tempObj$classOmega,tempObj$classRho,tempObj$newtStop ),tempObj$numClass))

        }


        elboNew=calcELBOFull(tempObj)
        elboDiff=abs(elboNew-elboOld)
        elboOld=elboNew

      }
      elboOld=10
      elboDiff=10
      i=0

      while(elboDiff>tempObj$vbStop){
        if(tempObj$numClass>1){
          tempObj$classLambda=updateLambdaFunction(tempObj$tissuePhi,tempObj$dirEta,tempObj$classRho )
          tempObj$classRho=updateRhoFunction(tempObj$classLambda,tempObj$classOmega)
        }

        tempObj$tissuePhi=updatePhiFunction(tempObj$classLambda,tempObj$dirEta,tempObj$likelihoodSum)
        tempObj$dirEta=updateEtaFunction(tempObj$tissuePhi,tempObj$classLambda,tempObj$tissueBeta)

        elboNew=calcELBOFull(tempObj)
        elboDiff=abs(elboNew-elboOld)
        elboOld=elboNew

      }



      elboNewFull=calcELBOFull(tempObj)
      elboFullDiff=abs(elboNewFull-elboFullOld)
      elboFullOld=elboNewFull
    }

    actorObjList[[j]]=tempObj
    likeVec[j]=calcELBOFull(tempObj)
    print(paste("ELBO:",likeVec[j]))
  }
  end=proc.time()

  bestIter=which.max(likeVec)

  modelList=actorObjList[[bestIter]]
  class(modelList)="actorFit"
  modelList$ELBO=likeVec[bestIter]
  modelList$Time=end-start
  return(modelList)

}
