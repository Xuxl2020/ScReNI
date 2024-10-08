
#' Calculate the precision and recall of top gene pairs from single-cell specific network
#'
#' @param scNet_precision_recall
#' @param top_number
#' 
#' @return
#' @export
#'
#' @examples

Calculate_scNetwork_precision_recall_top <- function(scNet_precision_recall, top_number){
  
  ###get the precision and recall of top gene pairs for all methods except of CSN
  scNet_precision_recall_top1 <- c()
  for(i in 1:length(top_number)){
   scNet_precision_recall_top <- scNet_precision_recall[[i]]
   scNet_precision_recall_top1 <- rbind(scNet_precision_recall_top1, cbind(top_number[i], scNet_precision_recall_top))
  }
  colnames(scNet_precision_recall_top1)[1] <- 'top_number'  
  scNet_precision_recall_top1$top_number <- as.factor(scNet_precision_recall_top1$top_number)
  for(i in 3:ncol(scNet_precision_recall_top1)){ 
    scNet_precision_recall_top1[, i] <- as.numeric(scNet_precision_recall_top1[, i])
  }
      
  scNet_precision_top <- summarySE(scNet_precision_recall_top1, measurevar="precision", groupvars=c("scNetwork_type", "top_number"))
  scNet_recall_top <- summarySE(scNet_precision_recall_top1, measurevar="recall", groupvars=c("scNetwork_type", "top_number"))
  
  return(list(scNet_precision_top, scNet_recall_top))
}

