#' Print the available list of GTEx tissues.
#'
#'
#' @return A vector of all available tissues
#'
#' @export

availableTissues=function(){
  return(colnames(alphaDirSet)[-(1:2)])

}
