#' Title
#'
#' @param scNetworks
#' @param top
#' @param cell_type_annotation
#' @param ntype
#' @param column_name
#' 
#' @return
#' @export
#'
#' @examples

calculate_scNetwork_degree <- function(scNetworks, top, cell_type_annotation, ntype=NULL, column_name='undup.sub.celltype') {
  
  # Function to count nonzero elements
  count_nonzero <- function(scNet) {
    sum(scNet != 0)
  }
  
  # Extract network types
  Network_types <- names(scNetworks)
  
  # Initialize list to store results
  degree_all <- list()
  
  for (i in seq_along(Network_types)) {
    scNet_list <- scNetworks[[i]]
    cell_num <- length(scNet_list)
    weights <- scNet_list
    
    # Calculate number of nonzero pairs
    nonzero_num <- sapply(scNetworks[[1]], count_nonzero)
    
    # Initialize matrices
    indegree <- matrix(0, nrow(weights[[1]]), cell_num)
    outdegree <- matrix(0, nrow(weights[[1]]), cell_num)
    rownames(indegree) <- rownames(outdegree) <- rownames(weights[[1]])
    colnames(indegree) <- colnames(outdegree) <- seq_len(cell_num)
    
    for (j in seq_len(cell_num)) {
      if(Network_types[i] != 'CSN') {
        weights1 <- weights[[j]]
        idx <- order(weights1, decreasing = TRUE)[1:top[j]]
        weights1[idx] <- 1
        weights1[-idx] <- 0
        indegree[, j] <- rowSums(weights1)
        outdegree[, j] <- colSums(weights1)
      } else {
        indegree[, j] <- rowSums(weights[[j]])
        outdegree[, j] <- colSums(weights[[j]])
      }  
    }
    
    # Store degree matrices
    degree <- list(indegree = indegree, outdegree = outdegree)
    
    # Iterate over in/out degree
    for (j in 1:2) {
      # Create Seurat object and perform dimensionality reduction
      degree_data <- if (j == 1) indegree else outdegree
      degree_umap <- CreateSeuratObject(counts = degree_data + 1)
      degree_umap <- NormalizeData(degree_umap)
      degree_umap <- ScaleData(degree_umap)
      degree_umap <- FindVariableFeatures(degree_umap, selection.method = "vst", nfeatures = 4000)
      degree_umap <- RunPCA(degree_umap, features = VariableFeatures(object = degree_umap))
      degree_umap <- RunUMAP(degree_umap, dims = 1:10)
      degree_umap <- FindNeighbors(degree_umap, dim=1:20) %>% FindClusters(resolution = 0.5)
      
      # Cluster data
      Cluster_umap <- Idents(degree_umap)
      degree_hclust <- dist(cor(log(degree_data + 1)))
      hc1 <- hclust(degree_hclust)
      Cluster_hclust <- cutree(hc1, k=ntype)
      
      # Calculate ARI
      True_label <- cell_type_annotation
      degree_umap_ARI <- flexclust::randIndex(table(Cluster_umap, True_label), correct = TRUE)
      degree_hclust_ARI <- flexclust::randIndex(table(Cluster_hclust, True_label), correct = TRUE)
      
      # Store results
      if(j == 1) {
        degree[['in.degree.umap']] <- degree_umap
        degree[['in.degree.umap.ARI']] <- degree_umap_ARI
        degree[['in.degree.hclust.ARI']] <- degree_hclust_ARI
      } else {
        degree[['out.degree.umap']] <- degree_umap
        degree[['out.degree.umap.ARI']] <- degree_umap_ARI
        degree[['out.degree.hclust.ARI']] <- degree_hclust_ARI
      }
    }
    
    degree_all[[Network_types[i]]] <- degree
  }
  
  return(degree_all)
}

