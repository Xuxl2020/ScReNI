###source('/data2/jiewang/Bioinformatics/Networks/SingleCnetwork/Programs/Plot_scNetwork_clustering.R')
###Perform clustering based on the gene degrees of cell-specific networks

Plot_scNetwork_clustering <- function(degree_all, path, data.name, cell_type_annotation, col,
                                      cell_type_column){
  
  library(ComplexHeatmap)
  library(Seurat) ###remotes::install_version("Seurat", "5.0.1", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
  
  Network_types <- names(degree_all)
  
  for(i in 1:length(degree_all)){
    print(Network_types[i])
    
    degree <- degree_all[[i]]
    indegree <- degree[['indegree']]
    outdegree <- degree[['outdegree']]
    in.degree <- degree[['in.degree.umap']]
    out.degree <- degree[['out.degree.umap']]
    
    ###plot the heatmap based on the similarity of the degree matrix
    column_ha = HeatmapAnnotation("Cell_type" = cell_type_annotation,
                                  col = col)
    indegree_heatmap <- Heatmap(cor(log(indegree+1)), top_annotation = column_ha, show_row_names = F, show_column_names = F)
    pdf(paste0(path, data.name, '_', Network_types[i], "_indegree_correlation_heatmap.pdf"), width = 6.5, height = 5)
    print(indegree_heatmap)
    dev.off()
    
    outdegree_heatmap <- Heatmap(cor(log(outdegree+1)), top_annotation = column_ha, show_row_names = F, show_column_names = F)
    pdf(paste0(path, data.name, '_', Network_types[i], "_outdegree_correlation_heatmap.pdf"), width = 6.5, height = 5)
    print(outdegree_heatmap)
    dev.off()
    
    ###plot UMAP based on the degree matrix
    pdf(paste0(path, data.name, '_', Network_types[i], "_indegree_umap.pdf"), height=10, width=10)
    in.degree[['lab']] <- sub.celltype
    print(DimPlot(in.degree, reduction = "umap", group.by = "lab", label = F))
    dev.off()
    
    pdf(paste0(path, data.name, '_', Network_types[i], "_outdegree_umap.pdf"), height=10, width=10)
    out.degree[['lab']] <- sub.celltype
    print(DimPlot(out.degree, reduction = "umap", group.by = "lab", label = F))
    dev.off()
  }
}


