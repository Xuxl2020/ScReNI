
#' Title
#'
#' @param CSNall_precision_recall
#' @param scNet_precision_recall_top
#' @param width
#' @param height
#' @param values
#' @param path
#' @param data.name
#' 
#' @return
#' @export
#'
#' @examples

Plot_scNetwork_precision_recall <- function(CSNall_precision_recall, scNet_precision_recall_top, width=7, height=5, values, path, data.name){
  
  library(ggplot2)
  
  ###plot the precision and recall of CSN and all other methods
  pdf(paste0(path, data.name, '_CSNall_precision.pdf'), width = width, height = height)
  p <- ggplot(CSNall_precision_recall, aes(x=scNetwork_type, y=precision, color=scNetwork_type)) + 
    geom_boxplot() + scale_color_manual(values=values)
  print(p + theme_classic() + theme(axis.text=element_text(size=14), axis.title=element_text(size=16)))
  dev.off()
  
  pdf(paste0(path, data.name, '_CSNall_recall.pdf'), width = width, height = height)
  p <- ggplot(CSNall_precision_recall, aes(x=scNetwork_type, y=recall, color=scNetwork_type)) + 
    geom_boxplot() + scale_color_manual(values=values)
  print(p + theme_classic() + theme(axis.text=element_text(size=14), axis.title=element_text(size=16)))
  dev.off()
  
  
  ###plot the precision and recall of top gene pairs for all methods except CSN
  pdf(paste0(path, data.name, '_Alltop_precision.pdf'), width = width, height = height)
  p <- ggplot(scNet_precision_recall_top[[1]], aes(x=top_number, y=precision, fill=scNetwork_type)) + 
    geom_bar(stat="identity", position = position_dodge()) + scale_fill_manual(values=values) + 
    geom_errorbar(aes(ymin=precision-se, ymax=precision+se), width=.2, position=position_dodge(.9))
  print(p + theme_classic() + theme(axis.text=element_text(size=14), axis.title=element_text(size=16)))
  dev.off()
  
  pdf(paste0(path, data.name, '_Alltop_recall.pdf'), width = width, height = height)
  p <- ggplot(scNet_precision_recall_top[[2]], aes(x=top_number, y=recall, fill=scNetwork_type)) + 
    geom_bar(stat="identity", position = position_dodge()) + scale_fill_manual(values=values) + 
    geom_errorbar(aes(ymin=recall-se, ymax=recall+se), width=.2, position=position_dodge(.9))
  print(p + theme_classic() + theme(axis.text=element_text(size=14), axis.title=element_text(size=16)))
  dev.off()
  
}


