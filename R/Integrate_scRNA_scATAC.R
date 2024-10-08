
#' @export
#'

Integrate_scRNA_scATAC <- function(scRNAseq, scATACseq, AnchorsDim=50, IntegratedDimensions, KNN=20, data.type='unpaired', species='mouse'){
  
  if(data.type=='unpaired'){
    print("Step 1.1: integrate unpaired scRNA-seq and scATAC-seq")
    HarmonyDim = IntegratedDimensions
    library(Signac) ## 1.10.0
    library(SeuratObject) ## remotes::install_version("SeuratObject", "5.0.1", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
    library(Seurat) ## remotes::install_version("Seurat", "5.0.1", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
    library(Matrix) ## remotes::install_version("Matrix", version = "1.6-5")
    library(harmony) ## remotes::install_version("harmony", version = "0.1.1")
    library(sctransform) ## remotes::install_version("sctransform", version = "0.4.1")
    library(GenomeInfoDb)
    library(EnsDb.Mmusculus.v79)
    library(ggplot2)
    library(patchwork)
    library(future)
    plan("multisession", workers = 4)
    options(future.globals.maxSize = 100000 * 1024^2)    
    
    scRNAseq <- ScaleData(scRNAseq, features = rownames(scRNAseq))
    DefaultAssay(scATACseq) = "ACTIVITY"
    scATACseq <- ScaleData(scATACseq, features = rownames(scATACseq))
    
    transfer.anchors1 <- FindTransferAnchors(reference = scRNAseq, query = scATACseq, 
                                             features = VariableFeatures(object = scRNAseq),
                                             reference.assay = "RNA", query.assay = "ACTIVITY", 
                                             dims = 1:AnchorsDim, reduction = "cca") 
    
    #the imputation was restricted to variable genes from scRNA-seq
    genes.use <- VariableFeatures(scRNAseq)
    DefaultAssay(scRNAseq) <- "RNA"
    scRNA_object <- NormalizeData(scRNAseq)
    refdata <- GetAssayData(scRNA_object, assay = "RNA", slot = "data")[genes.use, ]
    
    imputation <- TransferData(anchorset = transfer.anchors1, refdata = refdata,
                               dims = 1:AnchorsDim, weight.reduction = scATACseq[["lsi"]])
    
    scATACseq[["RNA"]] <- imputation
    DefaultAssay(scATACseq) <- "RNA"
    
    scRNA_object <- AddMetaData(scRNAseq, metadata = rep("scRNAseq", nrow(scRNAseq@meta.data)), col.name = "datatypes")
    scATAC_object <- AddMetaData(scATACseq, metadata = rep("scATACseq", nrow(scATACseq@meta.data)), col.name = "datatypes")
    
    #run Harmony to integrate scRNAseq and scATACseq
    coembed <- merge(x = scRNA_object, y = scATAC_object)  
    DefaultAssay(coembed) <- "RNA"
    coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
    coembed <- RunPCA(coembed, features = genes.use)
    coembed <- RunHarmony(coembed, group.by = "datatypes", lambda = 0.5, 
                          dims.use = 2:HarmonyDim, reduction.use = 'lsi')
    coembed <- RunUMAP(coembed, reduction = "harmony", dims = 1:HarmonyDim)
    coembed <- FindNeighbors(coembed, k.param = KNN, features = VariableFeatures(object = coembed))
    
    return(coembed)
    
  }else if(data.type=='paired'){
    
    print("Step 1.1: integrate paired scRNA-seq and scATAC-seq")
    UmapDims = IntegratedDimensions
    
    library(dplyr) 
    library(Seurat) 
    
    # Create Seurat object
    if(class(scRNAseq)[1]=='Seurat'){
      scRNAseq = as.matrix(scRNAseq[["RNA"]]$counts)
    }

    if(class(scATACseq)[1]=='Seurat'){
      scATACseq = as.matrix(scATACseq[["ATAC"]]$counts)
    }
    colnames(scATACseq)  <- colnames(scRNAseq)

    pbmc <- CreateSeuratObject(counts = scRNAseq)
    
    if(species=='mouse'){ 
      genomeVersion = 'mm10'
    }else if(species=='human'){
      genomeVersion = 'hg38'
    }else{ 
      print('Unknown species')
    }
    
    chrom_assay <- CreateChromatinAssay(counts = scATACseq, 
                                        sep = c(":", "-"), genome = genomeVersion, min.cells = 0
    )
    
    pbmc[["ATAC"]] <- chrom_assay
    
    # RNA analysis
    DefaultAssay(pbmc) <- "RNA"
    pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:UmapDims, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
    
    # ATAC analysis
    DefaultAssay(pbmc) <- "ATAC"
    pbmc <- RunTFIDF(pbmc)
    pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
    pbmc <- RunSVD(pbmc)
    pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:UmapDims, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
    pbmc <- FindMultiModalNeighbors(pbmc, k.nn = KNN, reduction.list = list("pca", "lsi"), dims.list = list(1:UmapDims, 1:UmapDims), knn.range = 100)
    pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
    pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
    
    return(pbmc)
  }
  
}

