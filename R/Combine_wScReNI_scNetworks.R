
#' combine single-cell regulatory networks from wScReNI
#'
#' @param sub.scatac.top
#' @param network.path
#' @param cell.index
#' 
#' @return
#' @export
#'
#' @examples

Combine_wScReNI_scNetworks <- function(sub.scatac.top, network.path, cell.index=NULL){

  setwd(paste0(network.path, '/wScReNI'))
  cellnames <- colnames(sub.scatac.top)
  
  if(is.null(cell.index)){
     cell.index <- 1:length(cellnames)
  }
  
  wScReNI_scNetworks <- list()
  for(i in cell.index){
    print(c(i, cellnames[i]))
    wScReNI_scNetworks0 <- read.table(paste0(i, '.', cellnames[i], '.network.txt'))
    wScReNI_scNetworks0 <- as.matrix(wScReNI_scNetworks0)
    colnames(wScReNI_scNetworks0) <- rownames(wScReNI_scNetworks0)
    wScReNI_scNetworks[[i]] <- wScReNI_scNetworks0
  }
  names(wScReNI_scNetworks) <- cellnames
  
  return(wScReNI_scNetworks)
}

