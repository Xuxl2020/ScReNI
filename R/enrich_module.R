#' Title
#'
#' @param Kmeans_result
#' @param org.db
#' @param enrich.db
#' @param fun_num
#' @param pvalueCutoff
#' @param use_internal_data
#' @param organism
#' 
#' @return
#' @export
#'
#' @examples

enrich_module <- function(Kmeans_result, org.db, enrich.db,fun_num = 5,
                          pvalueCutoff = 0.05, use_internal_data = TRUE, organism = NULL) {
  all_gene <- Kmeans_result
  all_gene<-all_gene[order(all_gene$KmeansGroup),]
  le<-levels(as.factor(all_gene$KmeansGroup))
  for (i in le) {
    acc11 <- all_gene[all_gene$KmeansGroup == i,]$Symbol
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

