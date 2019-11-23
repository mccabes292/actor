#####Helper Functions####

#Function to calculate proportions of each iso
getSampleProportions <- function(cts, gene_id) {
  gene.cts <-rowsum(cts,  gene_id)
  total.cts <- gene.cts[match(gene_id, rownames(gene.cts)),]
  return(cts/total.cts)
}
#Function to calculate Gene Expression
getGeneExpression <- function(cts, gene_id) {
  gene.cts <-rowsum(cts,  gene_id)
  total.cts <- gene.cts[match(gene_id, rownames(gene.cts)),]
  return(total.cts)
}
#Function to calculate Gene Expression
getGeneExpression2 <- function(cts, gene_id) {
  gene.cts <-rowsum(cts,  gene_id)
  return(gene.cts)

}
##Function to calculate the parameters of the BB distribution. estvec is the proportion estimates and precVec is the precision
calculateDirichletAlpha=function(estVec,precVec  ){
  cn=colnames(estVec)
  estTemp=precVec[match(estVec$gene_id,precVec$gene_id),]
  dTempAlpha1=data.frame("gene_id"=estVec$gene_id,"feature_id"=estVec$feature_id,  estTemp[,-1]*(estVec[,-(1:2)]))
  colnames(dTempAlpha1)=cn
  return(dTempAlpha1)
}
#Function to subset gene names
strp=function(x){
  substr(x,1,15)
}

distFunc=function(x,y){
  sqrt(2*sum((sqrt(x)-sqrt(y))^2))
}
