#' Title
#'
#' @param gtf_data
#' @param scrna
#' @param scatac
#' @param motif_database
#' @param motif_pwm
#' @param genome_database
#' @param gene_name_type
#' @param upstream_len
#' @param downstream_len
#' @param threshold
#' @param nthread
#' @param pvalue
#' 
#' @return
#' @export
#'
#' @examples

Infer_gene_peak_relationships <- function(gtf_data, scrna, scatac, motif_database, 
         motif_pwm, genome_database, gene_name_type = 'symbol', upstream_len = 1000*250, 
         downstream_len = 1000*250, threshold = 0.1, nthread = 20, pvalue = 5^(-4) ) {
  
  # find overlaps between peaks and genes based on TSS
  peak_gene_overlap_GR <- peak_gene_overlap_GR1(gtf_data = gtf_data, scrna = scrna, 
                                                gene_name_type = gene_name_type, scatac = scatac,
                                                upstream_len = upstream_len, downstream_len = downstream_len)
  
  # Relationship links are assigned based on correlation between gene and peak
  overlap_gene_peak1 <- gene_peak_corr1(exprMatrix = as.matrix(scrna), scatac = scatac, 
                                        peak_gene_overlap_GR = peak_gene_overlap_GR, 
                                        threshold = threshold, nthread = nthread)
  
  # Organising overlap_gene_peak1 into GRanges format
  peak_gene_overlap <- data.frame(gene_id = overlap_gene_peak1[,1], 
                                  peak_id = overlap_gene_peak1[,2])
  peak_gene_overlap <- na.omit(peak_gene_overlap)
  peak_GR <- GRanges(peak_gene_overlap$peak_id)
  peak_GR$peak_id <- peak_gene_overlap$peak_id
  peak_gene_overlap_GR2 <- peak_GR
  peak_gene_overlap_GR2$overlap_gene <- peak_gene_overlap$gene_id
  
  motif_peak_pair_df <- Get_motif_peak_pair_df0(peak_gene_overlap_GR2, gene_name_type = gene_name_type, 
                                                motif_database, motif_pwm = motif_pwm, 
                                                genome_database, pvalue = pvalue, nthread = nthread)
  
  # match peaks, genes and TFs
  peak_gene_TF <- peak_gene_TF_match(motif_database, motif_peak_pair_df, peak_gene_overlap_GR2)
  
  return(peak_gene_TF)
}
