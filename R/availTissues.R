#' Print the available list of GTEx tissues. If no object is provided, the function provides all available GTEx tissues.
#'
#' @param actorObj Output from actor model
#' @return A vector of all available tissues to plot
#'
#' @export

availableTissues=function(actorObj=NULL){
  if(is.null(actorObj)){

    return(colnames(gtexRefPanel)[-(1:2)])
  }
  return(colnames(actorObj$tissuePhi))

}
