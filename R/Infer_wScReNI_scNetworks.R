###Construct single-cell regulatory networks using wScReNI

Infer_wScReNI_scNetworks <- function(exprMatrix, gene_peak_overlap_matrix, gene_peak_overlap_labs, nearest.neighbors.idx, 
                                       network.path, data.name, cell.index=NULL, nthread = 50, max.cell.per.batch = 10){

  # the function to calcualte regulatory weights among genes and peaks using the random forest
  gene_peak_randomForest <- function(exprMatrix, gene_peak_overlap_matrix, gene_peak_overlap_labs, K="sqrt", nb.trees=100, nthread=20, 
                                    importance.measure="IncNodePurity", seed=NULL, trace=TRUE, ...) {
    library(foreach)
    library(doParallel)
  
    # set random number generator seed if seed is given
    if (!is.null(seed)) { set.seed(seed) }
    
    # to be nice, report when parameter importance.measure is not correctly spelled
    if (importance.measure != "IncNodePurity" && importance.measure != "%IncMSE") {
      stop("Parameter importance.measure must be \"IncNodePurity\" or \"%IncMSE\"")
    }
    # Check if nodesize parameter is in the input arguments
    args <- list(...)
    nodesize.in.args <- "nodesize" %in% names(args)
    # transpose expression matrix to (samples x genes)
    exprMatrix <- t(as.matrix(exprMatrix))
    gene_peak_overlap_matrix <- t(as.matrix(gene_peak_overlap_matrix))
    # setup weight matrix
    num.samples <- dim(exprMatrix)[1]
    num.genes <- dim(exprMatrix)[2]
    gene.names <- colnames(exprMatrix)
    
    num.peaks <- dim(gene_peak_overlap_matrix)[2]
    peak.names <- colnames(gene_peak_overlap_matrix)
    input.gene.names <- gene.names
    
    # compute importances for every target gene
    cl <- parallel::makeCluster(nthread)
    doParallel::registerDoParallel(cl)
        
    weights <- foreach(target.gene.idx = seq(from=1, to=num.genes), .combine='cbind', .multicombine=TRUE)  %dopar% {
      target.gene.name <- gene.names[target.gene.idx]
      # remove target gene from input genes
      these.input.gene.names <- setdiff(input.gene.names, target.gene.name)
      x <- exprMatrix[,these.input.gene.names]
      x.new <- cbind(x, gene_peak_overlap_matrix)
      num.input.vars <- length(these.input.gene.names) + num.peaks
      y <- exprMatrix[,target.gene.name]
      
      if(sum(y)>0){
        # set mtry
        if (class(K) == "numeric") {
          mtry <- K
        } else if (K == "sqrt") {
          mtry <- round(sqrt(num.input.vars))
        } else if (K == "all") {
          mtry <- num.input.vars
        } else {
          stop("Parameter K must be \"sqrt\", or \"all\", or an integer")
        }
        if (trace) {
          cat(paste("K = ", mtry,", ", nb.trees, " trees\n\n", sep=""))
          flush.console()
        }
        if (importance.measure == "%IncMSE") {
          if (nodesize.in.args) {
            rf <- randomForest::randomForest(x.new, y, mtry=mtry, ntree=nb.trees, importance=TRUE,...)
          } else {
            # By default, grow fully developed trees
            rf <- randomForest::randomForest(x.new, y, mtry=mtry, ntree=nb.trees, importance=TRUE, nodesize=1,...)
          }
          
        } else {
          # Normalize output
          y <- y / sd(y)
          if (nodesize.in.args) {
            rf <- randomForest::randomForest(x.new, y, mtry=mtry, ntree=nb.trees, importance=FALSE,...)
          } else {
            # By default, grow fully developed trees
            rf <- randomForest::randomForest(x.new, y, mtry=mtry, ntree=nb.trees, importance=FALSE, nodesize=1,...)
          }
        }
        
        im <- randomForest::importance(rf)[,importance.measure]
        im.names <- names(im)
        
        weight.matrix = c()
        for(j in (1:num.genes)[-target.gene.idx]){
          y.peak.coef <- sum(im[im.names %in% gene_peak_overlap_labs@peaks[gene_peak_overlap_labs@genes %in% target.gene.name]])
          peakj.coef <- sum(im[im.names %in% gene_peak_overlap_labs@peaks[gene_peak_overlap_labs@genes %in% gene.names[j]]])
          
          gene.coef <- im[im.names %in% gene.names[j]]
          TFs <- gene_peak_overlap_labs@TFs[gene_peak_overlap_labs@genes %in% target.gene.name]
          TFs <- unique(unlist(strsplit(TFs, ';')))
          
          # genej is not TF
          if(gene.names[j]  %in% TFs){
            weight.matrix[j] <- gene.coef + y.peak.coef + peakj.coef
          } else {
            weight.matrix[j] <-  gene.coef + peakj.coef
          }
          
          # y is TF
          if(length(gene_peak_overlap_labs@labels[gene_peak_overlap_labs@genes %in% target.gene.name])==0){
            weight.matrix[target.gene.idx] <- 0
          }else if(unique(gene_peak_overlap_labs@labels[gene_peak_overlap_labs@genes %in% target.gene.name])=='TF'){
            weight.matrix[target.gene.idx] <-  y.peak.coef
          }else{
            weight.matrix[target.gene.idx] <- 0
          }
        }
      }else{
        weight.matrix=rep(0,num.genes)
      }
            
      weight.matrix
    }
    
    rownames(weights) <- gene.names
    colnames(weights) <- gene.names
    weights / num.samples
  }


  library(dplyr) ###1.1.2
  library(foreach)
  library(doParallel)
  library(Seurat) ###remotes::install_version("Seurat", "5.0.1", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
    
  colnames(gene_peak_overlap_matrix) <- colnames(exprMatrix)
  ncell = ncol(exprMatrix)
  
  dir.create(network.path)
  dir.create(paste0(network.path, "wScReNI"))
  
  print(paste("Total number of cells:", ncell))
  if(is.null(cell.index)){
    cell.start=1
    cell.end=ncell
    loop.num=ceiling(ncell/max.cell.per.batch)
  }else{
    cell.start=cell.index[1]
    cell.end=cell.index[length(cell.index)]
    loop.num=ceiling(length(cell.index)/max.cell.per.batch)
  }
  
  scNet_list2 <- list()
  system.time({
    for(loop in 1:loop.num){
      if(loop==loop.num){ ncell2 <- (max.cell.per.batch*(loop-1)+cell.start):cell.end
      }else{ ncell2 <- (max.cell.per.batch*(loop-1)+cell.start):(max.cell.per.batch*loop+cell.start-1) }
      print(paste('Cell', ncell2[1], 'to cell', ncell2[length(ncell2)]))
      
      nthread2 <- min(nthread, floor(length(ncell2)*1.5))
      cl <- parallel::makeCluster(nthread2)
      doParallel::registerDoParallel(cl)      
      
      scNet_list <- foreach(i = ncell2, .combine='c', .multicombine=TRUE)  %dopar% {
        wnn.expr <- exprMatrix[, c(i, nearest.neighbors.idx[i,]) ] 
        wnn.peak <- gene_peak_overlap_matrix[, c(i, nearest.neighbors.idx[i,]) ]
        # new computer
        sub_res <- gene_peak_randomForest(wnn.expr, wnn.peak, gene_peak_overlap_labs, nthread = 1)
        sc_res <- as.matrix(sub_res)
        exprMatrix.df <- as.data.frame(exprMatrix)
        tmp_name = colnames(exprMatrix.df)[i]
        write.table(sc_res, paste0(network.path, "wScReNI/", i, ".", tmp_name, ".network.txt"), sep="\t")
        list(sc_res)
      }
      
      stopCluster(cl)
      scNet_list2 <- c(scNet_list2, scNet_list)
    }
  })
  
  return(scNet_list2)
}
