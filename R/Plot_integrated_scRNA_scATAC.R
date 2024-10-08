#' Title
#'
#' @param coembed
#' @param figure.folder
#' @param data.name
#' @param data.type
#' @param groups
#' @param genes
#' 
#' @return
#' @export
#'
#' @examples

Plot_integrated_scRNA_scATAC <- function(coembed, figure.folder=result.path, data.name, data.type='unpaired', groups=c('datatypes', 'samples', 'celltypes'), genes='Rlbp1'){
  
  if(data.type=='unpaired'){
    
    setwd(figure.folder)
    
    print(DimPlot(coembed, group.by = groups[1], reduction = "umap", pt.size = 0.8, label = F) 
          | DimPlot(coembed, group.by = groups[2], reduction = "umap", pt.size = 0.8, label = F)
          | DimPlot(coembed, group.by = groups[3], reduction = "umap", pt.size = 0.8, label = T))

    pdf(paste0(data.name, "_integrated_scRNA_scATAC_", paste(groups, collapse='_'), ".pdf"), height=10, width=20)
    print(DimPlot(coembed, group.by = groups[1], reduction = "umap", pt.size = 0.8, label = F) 
          | DimPlot(coembed, group.by = groups[2], reduction = "umap", pt.size = 0.8, label = F)
          | DimPlot(coembed, group.by = groups[3], reduction = "umap", pt.size = 0.8, label = T))
    dev.off()
    
    pdf(paste0(data.name, "_integrated_scRNA_scATAC_", groups[1], ".pdf"), height=10, width=10)
    print(DimPlot(coembed, group.by = groups[1], pt.size = 0.8, reduction = "umap", label = F))
    dev.off()
    
    pdf(paste0(data.name, "_integrated_scRNA_scATAC_", groups[2], ".pdf"), height=10, width=10)
    print(DimPlot(coembed, group.by = groups[2], pt.size = 0.8, reduction = "umap", label = F))
    dev.off()
    
    pdf(paste0(data.name, "_integrated_scRNA_scATAC_", groups[1], "_", groups[2], ".pdf"), height=10, width=20)
    print(DimPlot(coembed, group.by = groups[1], pt.size = 0.8, reduction = "umap", split.by = groups[2], ncol=3, label = F))
    dev.off()
    
    pdf(paste0(data.name, "_integrated_scRNA_scATAC_", groups[2], "_", groups[1], ".pdf"), height=10, width=20)
    print(DimPlot(coembed, group.by = groups[2], pt.size = 0.8, reduction = "umap", split.by = groups[1], ncol=3, label = F))
    dev.off()
    
    pdf(paste0(data.name, "_integrated_scRNA_scATAC_", groups[3], "_", groups[1], ".pdf"), height=10, width=20)
    print(DimPlot(coembed, group.by = groups[3], pt.size = 0.8, reduction = "umap", split.by = groups[1], ncol=2, label = F))
    dev.off()
    
    ##plot gene expression
    pdf(paste0(data.name, "_integrated_scRNA_scATAC_D30_", paste(genes, collapse='_'), ".pdf"), height=10, width=20)
    print(FeaturePlot(coembed, features = genes, reduction = 'umap', max.cutoff = 3, ncol = 1))
    dev.off()
    
  }else if(data.type=='paired'){
    setwd(figure.folder) 
    pdf(paste0(data.name, "_wnnCelltype.pdf"), width = 16, height = 6)
    DimPlot(coembed, reduction = "wnn.umap", group.by = "celltypes")
    dev.off()
    
    pdf(paste0(data.name, "_RNACelltype.pdf"), width = 8, height = 6)
    DimPlot(coembed, reduction = "umap.rna", group.by = "celltypes")
    dev.off()
    
    pdf(paste0(data.name, "_ATACCelltype.pdf"), width = 8, height = 6)
    DimPlot(coembed, reduction = "umap.atac", group.by = "celltypes")
    dev.off()
  }
  
}


