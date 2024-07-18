Select_partial_cells_for_scNewtorks <- function(coembed, Partial_cell, cell.num, data.type='unpaired', scRNA_scATAC_neighbor_undup=NULL, cell.name=NULL, Celltypename='undup.cell.types'){
  
  if(data.type=='unpaired'){
    RNAname = 'rna' 
    ATACname = 'atac' 
    RNAindex <- grep(RNAname, colnames(scRNA_scATAC_neighbor_undup))
    ATACindex <- grep(ATACname, colnames(scRNA_scATAC_neighbor_undup))
    Celltypeindex <- grep(Celltypename, colnames(scRNA_scATAC_neighbor_undup))
    sub.coembed <- subset(coembed, celltypes %in% Partial_cell)
    scrna <- as.matrix(sub.coembed@assays$RNA@counts[, scRNA_scATAC_neighbor_undup[, RNAindex]])
    scatac <- as.matrix(sub.coembed@assays$ATAC@counts[, scRNA_scATAC_neighbor_undup[, ATACindex]])
    colnames(scatac) = colnames(scrna)
    scatac <- scatac[grep("chr", rownames(scatac)), ]
    # annotations of cell types
    cell.type <- scRNA_scATAC_neighbor_undup[, Celltypeindex]
    cell.type.num <- table(cell.type)
    cell.type.name <- names(cell.type.num)
    
    if(is.null(cell.name)){
      target_cell <- c()
      for(i in 1:length(cell.type.name)){
        target_cell <- c(target_cell, sample(which(cell.type==cell.type.name[i]), min(cell.type.num[i], cell.num)))
      }
      sub.scrna <- scrna[, target_cell]
      sub.scatac <- scatac[, target_cell]
    }else{
      sub.scrna <- scrna[, cell.name]
      sub.scatac <- scatac[, cell.name]
    }
    
    sub.celltype <- scRNA_scATAC_neighbor_undup[match(colnames(sub.scrna), scRNA_scATAC_neighbor_undup[, RNAindex]), ]
    
  }else if(data.type=='paired'){
    
    sub.coembed <- subset(coembed, celltypes %in% Partial_cell)
    
    sampled_cells <- lapply(unique(sub.coembed$celltypes), function(group) {
      cells_in_group <- which(coembed$celltypes == group)
      sample(cells_in_group, cell.num)
    })
    
    sampled_cells <- unlist(sampled_cells)
    sampled_cells_names <- colnames(sub.coembed)[sampled_cells]
    sub.coembed1 <- sub.coembed[, sampled_cells_names]
    sub.celltype <- as.data.frame(sub.coembed1@meta.data[,c("celltypes")], row.names = rownames(sub.coembed1@meta.data))
    colnames(sub.celltype) <- "sub.celltype"
    
    sub.scrna <- as.matrix(sub.coembed1@assays$RNA$counts)
    sub.scatac <- as.matrix(sub.coembed1@assays$ATAC$counts)
    colnames(sub.scatac) = colnames(sub.scrna)
  }
  
    sub.scatac <- sub.scatac[grep("chr", rownames(sub.scatac)), ]
    
    # format rownames and colnames in scatac
    colnames(sub.scatac) <- gsub('_', '-', colnames(sub.scatac))
    
    # standardized rownames for scCAT-seq data to follow 'chrX:start-end'
    sc_peaks_name <- strsplit(rownames(sub.scatac), '-')
    sc_peaks_name <- t(as.data.frame(sc_peaks_name))
    sc_peaks_name <- paste0(sc_peaks_name[,1], ':', sc_peaks_name[,2], '-', sc_peaks_name[,3])
    rownames(sub.scatac) <- sc_peaks_name
    
    return(list(sub.scrna, sub.scatac, sub.celltype))
  
}

