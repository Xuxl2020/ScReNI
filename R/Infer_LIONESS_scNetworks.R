
#' Infer single-cell specific networks using LIONESS
#'
#' @param exprMatrix
#' @param nCores
#' @param nTrees
#' 
#' @return
#' @export
#'
#' @examples

Infer_LIONESS_scNetworks <- function(exprMatrix, nCores = 60, nTrees = 100){
  
  library(GENIE3)
  library(foreach)
  library(doParallel)
  
  set.seed(100)
  genie3_res <- GENIE3::GENIE3(exprMatrix, nCores = nCores, nTrees = nTrees, verbose = TRUE)
  genie3_res <- rbind(genie3_res[nrow(genie3_res), ], genie3_res[-nrow(genie3_res), ])
  rownames(genie3_res) <- colnames(genie3_res)
  
  ncell = ncol(exprMatrix)
  
  system.time({
    cl <- parallel::makeCluster(nCores)
    doParallel::registerDoParallel(cl)
    genie3_scNetwork <- foreach(i = 1:ncell, .combine = 'c', .multicombine = TRUE)  %dopar% {
      sub_res <- GENIE3::GENIE3(exprMatrix[, -i], nCores = 1, nTrees = nTrees, verbose = TRUE)
      sc_res <- as.matrix(ncell*(genie3_res - sub_res) + sub_res)
      # sc_res <- ifelse(sc_res>0, sc_res, 0)
      list(sc_res)
    }
    stopCluster(cl)
  })
  names(genie3_scNetwork) <- colnames(exprMatrix)
  
  return(genie3_scNetwork)
  
}

