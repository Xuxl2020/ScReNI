

#' Infer single-cell specific networks using kScReNI
#'
#' @param exprMatrix
#' @param nfeatures
#' @param knn
#' @param nthread
#' @param nTrees
#'
#' @return
#' @export
#'
#' @examples


Infer_kScReNI_scNetworks <- function(exprMatrix, nfeatures=4000, knn=20, nthread=20, nTrees=100){
  
  library(Seurat) ###remotes::install_version("Seurat", "5.0.1", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
  library(doParallel) ###1.0.17
  library(GENIE3) ###1.22.0
  
  pbmc <- CreateSeuratObject(counts = exprMatrix)
  all.genes <- rownames(pbmc)  
  pbmc <- NormalizeData(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = nfeatures)
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  pbmc <- FindNeighbors(pbmc, k.param = knn, features = VariableFeatures(object = pbmc))
  mat <- as.matrix(pbmc@graphs$RNA_snn)
  
  ncell = ncol(exprMatrix) 
  
  system.time({
    cl <- parallel::makeCluster(nthread)
    doParallel::registerDoParallel(cl)
    scNet_list <- foreach(i = 1:ncell, .combine='c', .multicombine=TRUE)  %dopar% {
      snn.expr <- exprMatrix[, order(mat[i,],decreasing = T)[1:(knn+1)]] 
      set.seed(100)
      sub_res <- GENIE3::GENIE3(snn.expr, nCores = 1, nTrees = 100, verbose = TRUE)
      sub_res <- ifelse(sub_res=='NaN', 0, sub_res)
      list(sub_res)
    }  
    stopCluster(cl) 
  })
  names(scNet_list) <- colnames(exprMatrix)
  
  scNet_list
}
