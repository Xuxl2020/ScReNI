###source('/data2/jwang/Retina/RetReg/CellType/Programs/CACIMAR_R_NewFunctions.R')


enrich_module <- function(Kmeans_result=Kmeans_clustering_ENS, org.db=org.Gg.eg.db, enrich.db,fun_num = 5,
                          pvalueCutoff = 0.05, use_internal_data = TRUE, organism = NULL) {
  all_gene <- Kmeans_result
  all_gene<-all_gene[order(all_gene$KmeansGroup),]
  le<-levels(as.factor(all_gene$KmeansGroup))
  for (i in le) {
    acc11 <- all_gene[all_gene$KmeansGroup == i,]$NMDA.Symbol
    gene1 <- clusterProfiler::bitr(acc11, fromType = "SYMBOL",
                                   toType = c("ENTREZID"),
                                   OrgDb = org.db)
    if (enrich.db =='KEGG') {
      k1 <- clusterProfiler::enrichKEGG(gene = gene1$ENTREZID,
                                        pvalueCutoff = pvalueCutoff
                                        ,organism = organism
                                        ,use_internal_data = use_internal_data)
    }else if(enrich.db =='GO'){
      k1 = clusterProfiler::enrichGO(gene = gene1$ENTREZID,
                                     OrgDb = org.db,
                                     keyType = "ENTREZID",
                                     ont = "BP",
                                     pvalueCutoff = pvalueCutoff)
    }
    acc2 <- k1@result
    acc2$'-log10(q-value)' <- -log10(acc2$qvalue)
    acc2 <- acc2[order(-acc2$`-log10(q-value)`),]
    if (i=='1' | i == 1) {
      acc21 <- acc2[1:fun_num,]
      acc21$module <- rep(i,fun_num)
    }else{
      acc22 <- acc2[1:fun_num,]
      acc22$module<-rep(i,fun_num)
      acc21 <- rbind(acc21,acc22)
    }
  }
  acc21 <- acc21[,c(1,2,11,10,3:9)]
  return(acc21)
}




Str_to_GR <- function(x){
  sp = strsplit(x,split='-')
  chr = sapply(sp,function(x) x[[1]])
  s = as.numeric(sapply(sp,function(x) x[[2]]))
  e = as.numeric(sapply(sp,function(x) x[[3]]))
  GR_out = GRanges(chr,IRanges(s,e))
  return(GR_out)
}

Seqnames <- function(x){
  out= as.character(seqnames(x))
  return(out)
}

Start <- function(x){
  out= as.numeric(start(x))
  return(out)
}

End <- function(x){
  out= as.numeric(end(x))
  return(out)
}

Must_to_GR <- function(x){
  library('pbapply')
  chr_all = pblapply(x,Seqnames)
  start_all = pblapply(x,Start)
  end_all = pblapply(x,End)
  len_all = pblapply(x,function(x) length(x))
  #####
  chr_all = as.character(unlist(chr_all))
  start_all = as.numeric(unlist(start_all))
  end_all = as.numeric(unlist(end_all))
  ##### #####
  names_all = rep(names(x),len_all)
  GR_out = GRanges(chr_all,IRanges(start_all,end_all),motifs=names_all)
  return(GR_out)
}

motifs_select <- function(motif,gene){
  index <- c()
  if (stringr::str_sub(gene[1],1,3)=='ENS') {
    col_idx = 5
  }else{col_idx = 4}
  for (i in 1:nrow(motif)) {
    judge <- c()
    gene1 <- strsplit(motif[i,col_idx],';')[[1]]
    for (j in gene1) {
      if (j %in% gene) {
        judge <- c(judge,'YSE')
      }
    }
    if ('YSE' %in% judge) {
      index <- c(index,i)
    }
  }
  motif1 <- motif[index,]
  return(motif1)
}

#' Title
#'
#' @param GR
#' @param gene.use
#' @param motifdb
#' @param pvalue.cutoff
#' @param BSdb
#'
#' @return
#' @export
#'
#' @examples
identify_region_tfs <- function(GR,gene.use,motifdb,
                                pvalue.cutoff = 5e-05,BSdb){
  motif_use = motifs_select(motifdb,gene.use)
  PWM = Transfac_PWMatrixList
  PWM = PWM[motif_use$Accession]
  matched_motif <- motifmatchr::matchMotifs(PWM,
                                            GR,genome = BSdb,
                                            out='positions',p.cutoff = pvalue.cutoff)
  matched_motif <- Must_to_GR(matched_motif)
  overlapped_region <- findOverlaps(GR,matched_motif)
  enriched_tf <- c()
  for (i in unique(overlapped_region@from)) {
    all_motif <- matched_motif$motifs[overlapped_region[overlapped_region@from==i]@to]
    all_motif <- motifdb[motifdb$Accession%in%all_motif,]
    all_tf <- paste(unique(unlist(strsplit(all_motif$EnsemblID,';')))
                    ,collapse = ';')
    names(all_tf) <- gene.use[i]
    enriched_tf <- c(enriched_tf,all_tf)
  }
  regulation <- data.frame('TF'=enriched_tf,'Target'=names(enriched_tf))
  return(regulation)
}

overlap_peak_motif <- function(peak,motif,motifdb){
  overlaped = findOverlaps(peak,motif)
  peak_motif = cbind(as.data.frame(peak[overlaped@from]),as.data.frame(motif[overlaped@to]))
  peak_motif$TF = motifdb[match(peak_motif$motifs,motifdb$Accession),4]
  return(peak_motif)
}

make_tf_target <- function(atac_out){
  tf = atac_out$TF
  tf = paste0(tf,'#',atac_out$symbol)
  tf_target = unlist(map(tf,~paste_gene(.x)))
  return(tf_target)
}

paste_gene <- function(gene){
  tf = strsplit(gene,'#')[[1]][1]
  target = strsplit(gene,'#')[[1]][2]
  tf = unlist(strsplit(tf,';'))
  return(paste0(tf,'-',target))
}


filter_regulation_fimo <- function(fimo_regulation,regulatory_relationships){
  if (!'TF' %in% colnames(regulatory_relationships)) {
    stop('regulatory_relationships should contain "TF" column')
  }
  if (!'Target' %in% colnames(regulatory_relationships)) {
    stop('regulatory_relationships should contain "Target" column')
  }
  fimo_pair <- apply(fimo_regulation, 1, function(x1){
    TFs <- strsplit(x1[1],';')[[1]]
    Target <- x1[2]
    regulation <- paste(TFs,Target)
    return(regulation)
  })
  fimo_pair <- unlist(fimo_pair)
  regulation_pair <- paste(regulatory_relationships[,1],regulatory_relationships[,4])
  regulation1 <- regulatory_relationships[regulation_pair %in% fimo_pair,]
  return(regulation1)
}

split_motif <- function(motif){
  gene1 <- strsplit(motif[5],';')[[1]]
}

enrich_module <- function(Kmeans_result=Kmeans_clustering_ENS, org.db=org.Gg.eg.db, enrich.db ,fun_num = 5,
                          pvalueCutoff = 0.05, use_internal_data = TRUE, organism = NULL) {
  all_gene <- Kmeans_result
  all_gene<-all_gene[order(all_gene$KmeansGroup),]
  le<-levels(as.factor(all_gene$KmeansGroup))
  for (i in le) {
    i=1
    acc11<-rownames(all_gene[all_gene$KmeansGroup == i,])
    gene1 <- clusterProfiler::bitr(acc11, fromType = "ENSEMBL",
                                   toType = c("SYMBOL", "ENTREZID"),
                                   OrgDb = org.db)
    if (enrich.db =='KEGG') {
      k1 <- clusterProfiler::enrichKEGG(gene = gene1$ENTREZID,
                                        pvalueCutoff = pvalueCutoff
                                        ,organism = organism
                                        ,use_internal_data = use_internal_data)
    }else if(enrich.db =='GO'){
      k1 = clusterProfiler::enrichGO(gene = gene1$ENTREZID,
                                     OrgDb = org.db,
                                     keyType = "ENTREZID",
                                     ont = "BP",
                                     pvalueCutoff = pvalueCutoff)
    }
    acc2 <- k1@result
    acc2$'-log10(q-value)' <- -log10(acc2$qvalue)
    acc2 <- acc2[order(-acc2$`-log10(q-value)`),]
    if (i=='1' | i == 1) {
      acc21 <- acc2[1:fun_num,]
      acc21$module <- rep(i,fun_num)
    }else{
      acc22 <- acc2[1:fun_num,]
      acc22$module<-rep(i,fun_num)
      acc21 <- rbind(acc21,acc22)
    }
  }
  acc21 <- acc21[,c(1,2,11,10,3:9)]
  return(acc21)
}

