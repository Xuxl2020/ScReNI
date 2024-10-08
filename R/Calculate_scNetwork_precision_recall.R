
#' Calculate the precision and recall of single-cell specific network
#'
#' @param scNetworks
#' @param TF_target_pair
#' @param top_number
#' @param gene_id_gene_name_pair
#' @param gene_name_type
#' 
#' 
#' @return
#' @export
#'
#' @examples


Calculate_scNetwork_precision_recall <- function(scNetworks, TF_target_pair, top_number=c(1000, seq(2000, 10000, by=2000), 20000), gene_id_gene_name_pair=NULL, gene_name_type=NULL){
       
  cell_num <- length(scNetworks[[1]])
  nonzero_num <- unlist(lapply(scNetworks[['CSN']], function(scNet){
    num1 <- length(scNet[scNet!=0])
    } ) )
  nonzero_num1 <- as.matrix(nonzero_num)
  colnames(nonzero_num1) <- 'Number.of.nonzero.pairs'
  
  Network_types <- names(scNetworks)  
  top_number <- unique(c(0, top_number))
  scNet_precision_all <- c(); scNet_recall_all <- c()  
  scNet_precision_recall_all <- list()
  
  Index1 <- setdiff(1:length(scNetworks), grep('CSN', names(scNetworks)))
  scNetworks1 <- list()
  for(i in 1:length(Index1)){
    scNetworks1[[i]] <- scNetworks[[Index1[i]]]
  }
  
  for(k in 1:length(top_number)){
    if(top_number[k]==0){
      ###calculate precision and recall using the same number of gene pairs with CSN
      Network_types1 <- Network_types
      scNetworks2 <- scNetworks
      top_number1 <- nonzero_num1
    }else{
      ###calculate precision and recall using the top gene pairs
      Network_types1 <- Network_types[Network_types!='CSN']
      scNetworks2 <- scNetworks1
      top_number1 <- matrix(rep(top_number[k], cell_num), ncol=1)
    }
    
    scNet_precision_recall_all[[k]] <- c()
    scNet_precision_recall2 <- c()
    for(i in 1:length(scNetworks2)){
      scNet_list <- scNetworks2[[i]]
      
      scNet_precision_recall <- c()
      for(j in 1:length(scNet_list)){
        print(c(top_number[k], Network_types1[i], j))
        scNet_precision_recall0 <- calculate_precision_recall(scNet_list[[j]], TF_target_pair, top_number1[j, ], gene_id_gene_name_pair, gene_name_type)
        scNet_precision_recall <- rbind(scNet_precision_recall, c(Network_types1[i], scNet_precision_recall0))
      }
      rownames(scNet_precision_recall) <- names(scNet_list)
      scNet_precision_recall2 <- rbind(scNet_precision_recall2, scNet_precision_recall)
    }
    scNet_precision_recall2 <- data.frame(scNet_precision_recall2)
    colnames(scNet_precision_recall2) <- c('scNetwork_type', 'precision', 'recall')
    scNet_precision_recall2[, 1] <- as.factor(scNet_precision_recall2[, 1])
    for(i in 2:ncol(scNet_precision_recall2)){ 
      scNet_precision_recall2[, i] <- as.numeric(scNet_precision_recall2[, i])
    }
    scNet_precision_recall_all[[k]] <- scNet_precision_recall2
  }
  names(scNet_precision_recall_all) <- top_number
  
  return(scNet_precision_recall_all)
}

