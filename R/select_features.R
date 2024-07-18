###select genes by rowSums>0 in expression and vst in seurat

select_features <- function(scrna, nfeatures = 4000, datatype = 'RNA', SeuratLayer=T){
  library(Seurat) ###remotes::install_version("Seurat", "5.0.1", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

  scrna <- scrna[rowSums(scrna)>0, ]
  scrna_seurat <- CreateSeuratObject(scrna, assay = datatype)
  scrna_seurat <- NormalizeData(scrna_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  scrna_seurat <- FindVariableFeatures(scrna_seurat, selection.method = "vst", nfeatures=nfeatures)
  var.genes <- VariableFeatures(scrna_seurat)
  scrna_seurat_top <- subset(scrna_seurat, features= var.genes)
  if(SeuratLayer==T){
    scrna_seurat_data <- as.matrix(scrna_seurat_top@assays[[datatype]]@layers$data)
    scrna_seurat_counts <- as.matrix(scrna_seurat_top@assays[[datatype]]@layers$counts)
  }else{
    scrna_seurat_data <- as.matrix(scrna_seurat_top@assays[[datatype]]@data)
    scrna_seurat_counts <- as.matrix(scrna_seurat_top@assays[[datatype]]@counts)
  }
  rownames(scrna_seurat_data) <- rownames(scrna_seurat_top)
  colnames(scrna_seurat_data) <- colnames(scrna_seurat_top)
  scrna_seurat_data <- scrna_seurat_data[var.genes, ]
  rownames(scrna_seurat_counts) <- rownames(scrna_seurat_top)
  colnames(scrna_seurat_counts) <- colnames(scrna_seurat_top)
  scrna_seurat_counts <- scrna_seurat_counts[var.genes, ]

  return(list(data=scrna_seurat_data, counts=scrna_seurat_counts))
}

