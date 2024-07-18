###source('/data2/jiewang/Bioinformatics/Networks/SingleCnetwork/Programs/Deal_gene_information.R')

Deal_gene_information <- function(gtf_regions, gene_name_type = 'symbol'){
  #remove non-protein coding genes
  gtf_regions <- gtf_regions[gtf_regions$gene_type %in% "protein_coding",]
  gene_id_gene_name_pair <- unique(gtf_regions[,c("gene_id", "gene_name")])
  gene_id_gene_name_pair$gene_id <- gsub("[.][0-9]+", "", gene_id_gene_name_pair$gene_id)
  gene_id_gene_name_pair <- gene_id_gene_name_pair[grep("_PAR_Y", gene_id_gene_name_pair$gene_id, invert=T),]
  
  #remove one-to-multiple genes
  duplicated_gene_names <- gene_id_gene_name_pair$gene_name[duplicated(gene_id_gene_name_pair$gene_name)]
  duplicated_gene.df <- gene_id_gene_name_pair[gene_id_gene_name_pair$gene_name%in%duplicated_gene_names,]
  
  duplicated_gene_ID_for_remove <- c()
  
  for(i in unique(duplicated_gene.df$gene_name)){
    if(gene_name_type == 'symbol'){
      duplicated_gene_ID_for_remove <- c(duplicated_gene.df[duplicated_gene.df$gene_name%in%i, "gene_name"][1], duplicated_gene_ID_for_remove)
    } else if (gene_name_type == 'id'){
      duplicated_gene_ID_for_remove <- c(duplicated_gene.df[duplicated_gene.df$gene_name%in%i, "gene_id"][1], duplicated_gene_ID_for_remove)
    } else {
      'please input correct gene_name_type.'
    }
  }
  
  if(gene_name_type == 'symbol'){
    gene_id_gene_name_pair <- gene_id_gene_name_pair[!gene_id_gene_name_pair$gene_name%in%duplicated_gene_ID_for_remove,]
    rownames(gene_id_gene_name_pair) <- gene_id_gene_name_pair$gene_name
  } else if (gene_name_type == 'id'){
    gene_id_gene_name_pair <- gene_id_gene_name_pair[!gene_id_gene_name_pair$gene_id%in%duplicated_gene_ID_for_remove,]
    rownames(gene_id_gene_name_pair) <- gene_id_gene_name_pair$gene_id
  } else {
    'please input correct gene_name_type.'
  }

  return(gene_id_gene_name_pair)
}


calculate_precision_recall <- function(scNetwork_weights, TF_target_pair, top_number=1000, gene_id_gene_name_pair=NULL, gene_name_type=NULL){

  link <- reshape2::melt(as.matrix(scNetwork_weights), na.rm=T)
  link <- link[order(-link[, 3]), ]
  link <- data.frame(from.gene = link[,1], to.gene = link[,2], im = link[,3])
  
  if(is.null(gene_id_gene_name_pair)){
    rownames(link) <- paste(link$from.gene, link$to.gene, sep="_")
  }else if(!is.null(gene_name_type)){
    if(gene_name_type == 'symbol'){
      link$from.gene.id <- gene_id_gene_name_pair[link$from.gene, "gene_id"]
      link$to.gene.id <- gene_id_gene_name_pair[link$to.gene, "gene_id"]
      # delete NA values in the pairs
      link <- na.omit(link)
      # use the pairs as the row names
      rownames(link) <- paste(link$from.gene, link$to.gene, sep="_")
    } else if (gene_name_type == 'id'){
      link$from.gene.symbol <- gene_id_gene_name_pair[link$from.gene, "gene_name"]
      link$to.gene.symbol <- gene_id_gene_name_pair[link$to.gene, "gene_name"]
      # delete NA values in the pairs
      link <- na.omit(link)
      # use the pairs as the row names
      rownames(link) <- paste(link$from.gene.symbol, link$to.gene.symbol, sep="_")
    } else {
      'please input correct gene_name_type.'
    }
  }
  link[rownames(link)%in%TF_target_pair, "true_state"] <- 1
  numerator = length(which(link[1:top_number, ]$true_state==1))
  precision_denominator = length(link[1:top_number, 1])
  precision = numerator/precision_denominator
  
  recall_denominator = length(which(link$true_state==1))
  recall = numerator/recall_denominator

  return(c(precision, recall))
}


## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

