
#' Identify regulators enriched in each single-cell regulatory network
#'
#' @param regulatory_network
#' @param Kmeans_clustering_ENS
#' @param TFFDR1
#' @param TFFDR2
#' 
#' @return
#' @export
#'
#' @examples

Identify_enriched_scRegulators <- function(regulatory_network, Kmeans_clustering_ENS, TFFDR1 = 10, TFFDR2 = 10){
  library(IReNA)
  
  TFs_list <- network_analysis(regulatory_network, Kmeans_clustering_ENS, TFFDR1 = TFFDR1, TFFDR2 = TFFDR2)

  return(TFs_list)
}
