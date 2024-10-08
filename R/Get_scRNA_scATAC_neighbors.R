
#' Title
#'
#' @param coembed
#' @param celltypes
#' @param datatypes
#' 
#' @return
#' @export
#'
#' @examples

Get_scRNA_scATAC_neighbors <- function(coembed, celltypes='celltypes', datatypes='datatypes'){
  umap.embeddings <- coembed@reductions$umap@cell.embeddings
  cell.inf <- coembed@meta.data[, c(celltypes, datatypes)]
  
  umap.embeddings1 <- cbind(umap.embeddings, cell.inf[match(rownames(cell.inf), rownames(umap.embeddings)), ])
  rna.points <- umap.embeddings1[umap.embeddings1[, datatypes]=='scRNAseq', 1:ncol(umap.embeddings)]
  atac.points <- umap.embeddings1[umap.embeddings1[, datatypes]=='scATACseq', 1:ncol(umap.embeddings)]
  print(paste("Number of scRNAseq cells:", nrow(rna.points)))
  print(paste("Number of scATACseq cells:", nrow(atac.points)))
  
  ##calculate the distances from one point to multiple points
  distances <- apply(rna.points, 1, function(x)which.min(sqrt((atac.points[, 1] - x[1])^2 + (atac.points[, 2] - x[2])^2)))
  
  rna.points.id <- rownames(rna.points)
  atac.points.id <- rownames(atac.points)[distances]
  rna.points.celltype = cell.inf[match(rna.points.id, rownames(cell.inf)), celltypes]
  rna.atac.id <- data.frame(rna.points.id=rna.points.id, atac.points.id=atac.points.id, celltype=rna.points.celltype)
  
  ###remove scRNAseq cells which are matched to multiple scATACseq cells
  undup.rna.points.id <- rna.points.id[!duplicated(atac.points.id)]
  undup.atac.points.id <- atac.points.id[!duplicated(atac.points.id)] 
  undup.cell.types <- rna.points.celltype[!duplicated(atac.points.id)] 
  undup.rna.atac.id <- data.frame(undup.rna.points.id=undup.rna.points.id, undup.atac.points.id=undup.atac.points.id, undup.cell.types=undup.cell.types)

  return(list(rna.atac.id, undup.rna.atac.id))
}

