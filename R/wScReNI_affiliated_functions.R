#' Title
#'
#' @param gtf_data
#' @param scrna
#' @param gene_name_type
#' @param scatac 
#' @param upstream_len
#' @param downstream_len
#' 
#' @return
#' @export
#'
#' @examples
# overlaps between genes (target genes and TFs) and peaks by TSS
peak_gene_overlap_GR1 <- function(gtf_data, scrna, gene_name_type, scatac, upstream_len, downstream_len){
  scrna_gene_id <- rownames(scrna)
  gene_tss <- get_tss_region(gtf_data, scrna_gene_id, gene_name_type, upstream_len, downstream_len)
  gene_tss_GR <- GRanges(paste0(gene_tss[,2], ':', gene_tss[,3], '-', gene_tss[,4]))
  gene_tss_GR$gene_id = gene_tss[,1]
  peak_GR <- GRanges(rownames(scatac))
  peak_GR$peak_id <- rownames(scatac)
  overlap_ind <- findOverlaps(gene_tss_GR, peak_GR)
  overlap_gene <- gene_tss_GR[overlap_ind@from,]$gene_id
  overlap_peak_GR <- peak_GR[overlap_ind@to,]
  overlap_peak_GR$overlap_gene <- overlap_gene
  return(overlap_peak_GR)
}

#' Title
#'
#' @param exprMatrix
#' @param scatac
#' @param peak_gene_overlap_GR
#' @param threshold
#' @param nthread
#' @param downstream_len
#' 
#' @return
#' @export
#'
#' @examples

# correlation between gene and peak
gene_peak_corr1 <- function(exprMatrix, scatac, peak_gene_overlap_GR, threshold = 0.1, nthread=20){
  library(foreach)
  
  overlap_gene = peak_gene_overlap_GR$overlap_gene
  overlap_peak = peak_gene_overlap_GR$peak_id
  overlap_gene_peak = cbind(overlap_gene, overlap_peak)
  len <- length(overlap_gene)
  system.time({
    cl <- parallel::makeCluster(nthread)
    doParallel::registerDoParallel(cl)
    cors <- foreach(i = 1:len, .combine='c', .multicombine=TRUE)  %dopar% {
      
      id1 <- Matrix::which(exprMatrix[overlap_gene_peak[i,1],]!=0)
      id2 <- Matrix::which(scatac[overlap_gene_peak[i,2],]!=0)
      id <- intersect(id1, id2)
      cor(as.numeric(exprMatrix[overlap_gene_peak[i,1], id]), 
          as.numeric(scatac[overlap_gene_peak[i,2], id]), 
          method = "spearman")
    }
  })
  
  overlap_gene_peak1 <- overlap_gene_peak[which(abs(cors) > threshold),]
  return(overlap_gene_peak1)
}

#' Title
#'
#' @param peak_gene_overlap_GR2
#' @param gene_name_type
#' @param motif_database
#' @param motif_pwm
#' @param genome_database
#' @param pvalue=5^(-4)
#' @param nthread=10
#' @return
#' @export
#'
#' @examples

Get_motif_peak_pair_df0 <- function(peak_gene_overlap_GR2, gene_name_type, 
                                    motif_database, motif_pwm, 
                                    genome_database, pvalue=5^(-4), nthread=10){
  library(foreach)
  library(doParallel)
  
  overlap_gene <- peak_gene_overlap_GR2$overlap_gene
  overlap_gene_motif <- motifs_select(motif_database, overlap_gene, gene_name_type)
  overlap_gene_motif_GR <- motif_pwm[names(motif_pwm@listData) %in% overlap_gene_motif$Accession, ]
  peak_gene_overlap_GR.reduced <- reduce(peak_gene_overlap_GR2)
  len <- length(peak_gene_overlap_GR.reduced)
  
  system.time({
    cl <- parallel::makeCluster(nthread)
    doParallel::registerDoParallel(cl)
    motif_ix_list <- foreach(i = 1:len, .combine=rbind, .multicombine=TRUE)  %dopar% {
      # 在peak鉴定motif
      motif_ix <- motifmatchr::matchMotifs(overlap_gene_motif_GR, peak_gene_overlap_GR.reduced[i],
                                           genome = genome_database, out = 'positions', p.cutoff = pvalue)
      # 筛除未在peak区域鉴定到的motif
      length_motif_ix <- lapply(motif_ix, length)
      motif_id <- names(length_motif_ix[length_motif_ix>0])
      # 将结果整合为motif-peak的数据框
      tmp_peak.df <- as.data.frame(peak_gene_overlap_GR.reduced[i])
      tmp_peak_name <- paste0(tmp_peak.df$seqnames, ":", tmp_peak.df$start, "-", tmp_peak.df$end)
      tmp_peak_motif.df <- data.frame(motif_id=motif_id, peak_name=rep(tmp_peak_name, length(motif_id)))
      colnames(tmp_peak_motif.df) <- c("motif_id", "peak_name")
      return(tmp_peak_motif.df)
    }
    stopCluster(cl)
  })
  
  return(motif_ix_list)
}

#' Title
#'
#' @param motif_database
#' @param motif_peak_pair_df
#' @param peak_gene_overlap_GR2
#' 
#' @return
#' @export
#'
#' @examples

peak_gene_TF_match <- function(motif_database, motif_peak_pair_df, peak_gene_overlap_GR2){
  # Matching tfs and motifs
  TF_motif_pair.list <- lapply(1:nrow(motif_database), get_TF_motif_pair, ref=motif_database)
  TF_motif_pair.df <- do.call("rbind", TF_motif_pair.list)
  # 根据ref（tf-motif）关系，关联TF-gene之间的关系
  # 获取TF-motif-peak之间的关系
  TF_motif_peak_pair.df <- merge(motif_peak_pair_df, TF_motif_pair.df, by = "motif_id")
  # 获取TF-motif-peak-gene之间的关系
  peak_gene_pair.df <- as.data.frame(peak_gene_overlap_GR2)[, c("peak_id", "overlap_gene")]
  colnames(peak_gene_pair.df) <- c("peak_name", "overlap_gene")
  TF_motif_peak_gene.pair.df <- merge(TF_motif_peak_pair.df, peak_gene_pair.df, by = "peak_name")
  tmp_TF_peak_gene.df <- unique(TF_motif_peak_gene.pair.df[, c("peak_name", "TFs", "overlap_gene")])
  peak_gene_TF <- aggregate(tmp_TF_peak_gene.df$TFs, 
                            by = list(tmp_TF_peak_gene.df$peak_name, tmp_TF_peak_gene.df$overlap_gene), 
                            FUN = paste0, collapse=";")
  colnames(peak_gene_TF) <- c( "peak.name", "gene.name", "TFs")
  return(peak_gene_TF)
}

#' Title
#'
#' @param motif_database
#' @param overlap_gene
#' @param gene_name_type
#' 
#' @return
#' @export
#'
#' @examples
# find motifs
motifs_select <- function(motif_database, overlap_gene, gene_name_type){
  if(gene_name_type == 'symbol'){
    motif_mat <- as.matrix(motif_database[[4]])
  } else if (gene_name_type == 'id'){
    motif_mat <- as.matrix(motif_database[[5]])
  } else {
    'please input correct gene_name_type.'
  }
  b <- apply(motif_mat, 1, function(x){
    a = sum(unlist(strsplit(x,';')) %in% overlap_gene)
    ifelse(a>0, 1, 0)
  })
  motif <- motif_database[which(b==1),]
  return(motif)
}

#' Title
#'
#' @param scatac
#' @param peak_gene_TF
#' @param pvalue
#' 
#' @return
#' @export
#'
#' @examples
peakMat <- function(scatac, peak_gene_TF, pvalue=10^(-5)){
  peak.matrix_gene_peak_overlap <- as.matrix(scatac[rownames(scatac) %in% peak_gene_TF$peak.name, ])
  e_atac <- matrix(rnorm(nrow(peak.matrix_gene_peak_overlap)*ncol(peak.matrix_gene_peak_overlap), 0, pvalue), 
                   nrow(peak.matrix_gene_peak_overlap), ncol(peak.matrix_gene_peak_overlap))
  peak.matrix_gene_peak_overlap <- peak.matrix_gene_peak_overlap + e_atac

  return(peak.matrix_gene_peak_overlap)
}

#' Title
#'
#' @param peak_gene_TF
#' 
#' @return
#' @export
#'
#' @examples
peak_gene_TF_labs <- function(peak_gene_TF){
  
  labels <- rep('target', length(peak_gene_TF$gene.name))
  
  TFs <- unique(unlist(strsplit(peak_gene_TF$TFs, ';')))
  TF_match <- TFs[TFs %in% peak_gene_TF$gene.name]
  idx <- which(peak_gene_TF$gene.name %in% TF_match)
  labels[idx] <- 'TF'
  
  # Define the corresponding matrix of genes, peaks, and transcription factors
  setClass("information match", slots = list(labels="character", genes="character", peaks="character", TFs="character"))
  labs <- new("information match", labels = labels, genes = peak_gene_TF$gene.name, peaks = peak_gene_TF$peak.name, TFs = peak_gene_TF$TF)

  return(labs)
}

#' Title
#'
#' @param gtf_data
#' @param scrna_gene_id
#' @param gene_name_type
#' @param upstream_len
#' @param downstream_len
#' @return
#' @export
#'
#' @examples

# get promoter region of gene
get_tss_region <- function(gtf_data, scrna_gene_id, gene_name_type, upstream_len, downstream_len){
  gtf_gene <- gtf_data[gtf_data$type == 'gene',]
  if(gene_name_type == 'symbol'){
    gtf_gene_id <- gtf_gene$gene_name
  } else if (gene_name_type == 'id'){
    gtf_gene_id <- gtf_gene$gene_id
  } else {
    'please input correct gene_name_type.'
  }
  gene_match_inf <- gtf_gene[gtf_gene_id %in% scrna_gene_id,]
  gene_id <- gtf_gene_id[gtf_gene_id %in% scrna_gene_id]
  start <- ifelse(gene_match_inf$strand == '-', gene_match_inf$end - downstream_len, gene_match_inf$start - upstream_len )
  end <- ifelse(gene_match_inf$strand == '-', gene_match_inf$end + upstream_len, gene_match_inf$start + downstream_len)
  if(gene_name_type == 'symbol'){
    chr <- gene_match_inf$seqnames
  } else if (gene_name_type == 'id'){
    chr <- paste0('chr', gene_match_inf$seqnames)
  } else {
    'please input correct gene_name_type.'
  }
  tss_info <- data.frame(gene_name = gene_id, chr = chr, start = start, end = end)
  return(tss_info)
}


get_TF_motif_pair <- function(i, ref, gene_name_type='symbol'){
  
  MotifTFs_i <- ref[i,]
  
  if(gene_name_type == 'symbol'){
    TF_motif_pair <- data.frame(TFs = unlist(strsplit(MotifTFs_i$TFs, ";")),
                                motif_id = rep(MotifTFs_i$Accession, length(unlist(strsplit(MotifTFs_i$TFs,";")))))
  } else if (gene_name_type == 'id'){
    TF_motif_pair <- data.frame(TFs = unlist(strsplit(MotifTFs_i$EnsemblID,";")),
                                motif_id = rep(MotifTFs_i$Accession, length(unlist(strsplit(MotifTFs_i$EnsemblID,";")))))
  } else {
    'please input correct gene_name_type.'
  }
}




