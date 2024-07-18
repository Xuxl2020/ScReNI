###source('/data2/jwang/Retina/RetReg/CellType/Programs/reconstruct_network_part.R')

get_Enriched_TFs <- function(GeneCor1, Kmeans_result, TFFdrThr1 = 2) {
  GeneCor1$TF <- factor(GeneCor1$TF)
  GeneCorP1 <- GeneCor1[GeneCor1$Correlation > 0, ]
  GeneCorN1 <- GeneCor1[GeneCor1$Correlation < 0, ]
  Module1 <- Kmeans_result
  TF0 <- table(GeneCor1$TF)
  TF01 <- rep(0, length(TF0))
  names(TF01) <- names(TF0)
  TFP01 <- table(GeneCorP1$TF)
  TFN01 <- table(GeneCorN1$TF)
  TFP1 <- TF01
  TFP1[match(names(TFP01), names(TF01))] <- TFP01
  TFN1 <- TF01
  TFN1[match(names(TFN01), names(TF01))] <- TFN01

  uGroup1 <- sort(unique(Module1$KmeansGroup))
  pTF4 <- c()
  for (i in 1:length(uGroup1)) {
    Module2 <- rownames(Module1)[Module1$KmeansGroup == uGroup1[i]]
    for (j in 1:2) {
      if (j == 1) {
        GeneCorPN <- GeneCorP1
        TF1 <- TFP1
        GeneCor2 <- GeneCorP1[is.element(GeneCorP1$Target, Module2), ]
      } else {
        GeneCorPN <- GeneCorN1
        TF1 <- TFN1
        GeneCor2 <- GeneCorN1[is.element(GeneCorN1$Target, Module2), ]
      }
      TF02 <- table(GeneCor2$TF)
      TF2 <- TF01
      TF2[match(names(TF02), names(TF01))] <- TF02
      TF3 <- cbind(TF2, TF1, nrow(GeneCorPN) - TF1, nrow(GeneCor2))

      ### perform hypergeometric test
      pTF3 <- apply(TF3, 1, function(X1) {
        X1 <- as.numeric(X1)
        if (X1[1] == 0 | X1[1] < 5 & X1[1] < 0.02 * X1[4]) {
          P1 <- 1
        } else {
          P1 <- phyper(X1[1], X1[2], X1[3], X1[4], lower.tail = FALSE)
        }
      })
      pTF31 <- p.adjust(pTF3, method = "fdr")
      pTF32 <- cbind(apply(TF3, 1, function(x1) {
        paste(x1, collapse = ";")
      }), pTF3, pTF31)
      pTF32 <- as.data.frame(pTF32)
      pTF32[,2] <- as.numeric(pTF32[,2])
      pTF32[,3] <- as.numeric(pTF32[,3])
      if (i == 1 & j == 1) {
        pTF4 <- pTF32
      } else {
        pTF4 <- cbind(pTF4, pTF32)
      }
    }
    colnames(pTF4)[(ncol(pTF4) - 5):ncol(pTF4)] <- paste0(c("Pnum", "Pp", "Pfdr"
                                                            , "Nnum", "Np",
                                                            "Nfdr"), uGroup1[i])
  }
  pTF4Mod <- Module1[match(rownames(pTF4), rownames(Module1)), ]
  Ind1 <- 1:length(uGroup1)
  pTF4Min <- t(apply(pTF4[, grep("fdr", colnames(pTF4))], 1, function(x1) {
    x12 <- as.numeric(as.character(x1))
    x2 <- -log10(x12)
    x21 <- x2[seq(1, length(x2), 2)]
    x22 <- x2[seq(2, length(x2), 2)]
    x212 <- paste(Ind1[x21 > TFFdrThr1], collapse = ";")
    x222 <- paste(Ind1[x22 > TFFdrThr1], collapse = ";")
    if (x212 == "") {
      x212 <- "NA"
    }
    if (x222 == "") {
      x222 <- "NA"
    }

    x3 <- max(x2)
    ix2 <- which.max(x2)
    if (ix2 %% 2 == 1) {
      ix3 <- paste0("P", floor((ix2 + 1) / 2))
    } else {
      ix3 <- paste0("N", floor((ix2 + 1) / 2))
    }

    return(c(x3, ix3, x212, x222))
  }))
  colnames(pTF4Min) <- c("TFMinNlogfdr", "TFMinGroup", "SigActModules", "SigRepModules")
  pTF4Min <- as.data.frame(pTF4Min)
  pTF4Min[,1] <- as.numeric(pTF4Min[,1])
  pTF42 <- cbind(pTF4Mod, pTF4Min, pTF4)

  #### Enriched TFs
  Ind42 <- sort(unique(c(1:nrow(pTF42))[as.numeric(as.character(pTF42[, "TFMinNlogfdr"])) > TFFdrThr1]))
  EnrichTF1 <- pTF42[Ind42, ]
  print(paste("Total TFs:", nrow(pTF42)))
  print(paste("Enriched TFs:", nrow(EnrichTF1)))

  GeneCorTFfdr <- pTF4Min[match(as.character(GeneCor1$TF), rownames(pTF4Min)), ]
  Regulation <- apply(GeneCor1, 1, function(x1) {
    x2 <- "Positive"
    if (as.numeric(x1["Correlation"]) < 0) {
      x2 <- "Negative"
    }
    return(x2)
  })

  GeneCor3 <- cbind(GeneCor1[, grep("TF", colnames(GeneCor1))], GeneCorTFfdr,
                    GeneCor1[, c(grep("Target", colnames(GeneCor1))[1]:ncol(GeneCor1))],
                    Regulation)
  #### Edges' regulator belongs to enriched TFs
  EnTFReg1 <- GeneCor3[is.element(GeneCor3$TF, rownames(EnrichTF1)), ]
  #### Edges' regulator and target belong to enriched TFs
  EnTFTarg1 <- GeneCor3[is.element(GeneCor3$TF, rownames(EnrichTF1)) & is.element(GeneCor3$Target,
                                                                                  rownames(EnrichTF1)), ]

  #### Get edges within each module
  EnTFGroup1 <- sort(unique(c(unique(EnTFTarg1$TFGroup), unique(EnTFTarg1$TargetGroup))))
  EnTFTarg2 <- c()
  for (i in 1:length(EnTFGroup1)) {
    EnTFTarg2 <- rbind(EnTFTarg2, EnTFTarg1[EnTFTarg1$TFGroup == EnTFGroup1[i] &
                                              EnTFTarg1$TargetGroup == EnTFGroup1[i], ])
  }
  list1 <- list(pTF42, EnrichTF1, EnTFReg1, EnTFTarg1, EnTFTarg2)
  names(list1) <- c("Cor_TFs", "Cor_EnTFs", "FOSF_RegMTF_Cor_EnTFs",
                    "FOSF_RegMTF_Cor_EnTFsTarg", "FOSF_RegMTF_Cor_EnTFsTargM")
  return(list1)
}


get_regulation_of_TFs_to_modules <- function(TFs_list, Thr = 2) {
  con1 <- TFs_list[['Cor_EnTFs']]
  Thr1 <- Thr
  name1 <- colnames(con1)[grep("^[PN]fdr\\d+$", colnames(con1))]
  ind1 <- grep("^[PN]fdr\\d+$", colnames(con1))
  col1 <- c()
  col1 <- c(paste("TF", "TFSymbol", "TFGroup", "TargetModule", "TargetGroup",
                  "Regulation", "Nlogfdr", sep = "\t"))
  TF <- c()
  for (i in 1:nrow(con1)) {
    for (j in ind1) {
      if (con1[i, ][j] == 0) {
        fdr1 <- Inf
      }
      else {
        fdr1 <- -log(con1[i, ][j]) / log(10)
      }
      if (fdr1 > Thr1) {
        name2 <- colnames(con1)[j]
        var1 <- paste(con1[i, ][1:2], collapse = "\t")
        if (grepl("^P", name2) == TRUE) {
          name3 <- "Positive"
        } else if (grepl("^N", name2) == TRUE) {
          name3 <- "Negative"
        }
        TG <- stringr::str_extract(name2, "\\d")
        var2 <- paste(rownames(con1)[i], var1, paste0("Group", TG), TG,
                      name3, fdr1, sep = "\t")
        col1 <- c(col1, var2)
        TF <- c(TF, con1[i, ][1]$Symbol)
      }
    }
  }
  col2 <- as.data.frame(col1)
  TF_list <- TF[!duplicated(TF)]
  TF_module_regulation <- strsplit(col2[1,], '\t')[[1]]
  TF_module_regulation <- matrix(TF_module_regulation, ncol =
                                   length(TF_module_regulation))
  for (i in 2:nrow(col2)) {
    out2 <- strsplit(col2[i,], '\t')[[1]]
    TF_module_regulation <- rbind(TF_module_regulation,out2)
  }
  TF_module_regulation <- as.data.frame(TF_module_regulation)
  colnames(TF_module_regulation) <- TF_module_regulation[1,]
  TF_module_regulation <- TF_module_regulation[-1,]
  TF_module_regulation[,7] <- as.numeric(TF_module_regulation[,7])
  TFs_list[['TF_list']] <- TF_list
  TFs_list[['TF_module_regulation']] <- TF_module_regulation
  return(TFs_list)
}


get_partial_regulations <- function(TFs_list) {
  con1 <- TFs_list[['FOSF_RegMTF_Cor_EnTFs']]
  hash2 <- TFs_list[['TF_list']]
  con1$TFSymbol <- as.character(con1$TFSymbol)
  con1$TargetSymbol <- as.character(con1$TargetSymbol)
  rowcount <- c()
  for (i in 1:nrow(con1)) {
    if (con1[i, ][2] %in% hash2 & con1[i, ][9] %in% hash2 == TRUE) {
      rowcount <- c(rowcount, i)
    }
  }
  col1 <- con1[rowcount, ]
  TFs_list[['TF_network']] <- col1
  return(TFs_list)
}


merge_Module_Regulations <- function(TFs_list, Kmeans_result, ModuleThr1 = 0.05) {
  TF1 <- Kmeans_result
  Regulation1 <- TFs_list[['FOSF_RegMTF_Cor_EnTFsTarg']]
  Regulation1$Correlation <- as.numeric(Regulation1$Correlation)
  Module1 <- sort(unique(TF1$KmeansGroup))

  RegulationNum1 <- c()
  RegulationP <- Regulation1[Regulation1$Regulation == "Positive", ]
  RegulationN <- Regulation1[Regulation1$Regulation == "Negative", ]
  for (i in 1:length(Module1)) {
    Regulation1A <- Regulation1[Regulation1$TFGroup == Module1[i], ]
    Regulation1AP <- Regulation1A[Regulation1A$Regulation == "Positive", ]
    Regulation1AN <- Regulation1A[Regulation1A$Regulation == "Negative", ]
    Regulation2P <- RegulationP[RegulationP$TargetGroup == Module1[i], ]
    Regulation2N <- RegulationN[RegulationN$TargetGroup == Module1[i], ]
    for (j in i:length(Module1)) {
      Regulation2A <- Regulation1[Regulation1$TFGroup == Module1[j], ]
      Regulation2AP <- Regulation2A[Regulation2A$Regulation == "Positive", ]
      Regulation2AN <- Regulation2A[Regulation2A$Regulation == "Negative", ]
      Regulation1P <- RegulationP[RegulationP$TargetGroup == Module1[j], ]
      Regulation1N <- RegulationN[RegulationN$TargetGroup == Module1[j], ]

      Regulation12P <- Regulation1AP[Regulation1AP$TargetGroup == Module1[j], ]
      Regulation12N <- Regulation1AN[Regulation1AN$TargetGroup == Module1[j], ]
      Regulation21P <- Regulation2AP[Regulation2AP$TargetGroup == Module1[i], ]
      Regulation21N <- Regulation2AN[Regulation2AN$TargetGroup == Module1[i], ]

      Regulation12Pnum <- c(Module1[i], Module1[j], "Positive",
                            mean(Regulation12P$Correlation), nrow(Regulation12P),
                            nrow(Regulation1P), nrow(RegulationP) - nrow(Regulation1P), nrow(Regulation1AP))
      Regulation12Nnum <- c(Module1[i], Module1[j], "Negative",
                            mean(Regulation12N$Correlation), nrow(Regulation12N),
                            nrow(Regulation1N), nrow(RegulationN) - nrow(Regulation1N), nrow(Regulation1AN))
      Regulation21Pnum <- c(Module1[j], Module1[i], "Positive",
                            mean(Regulation21P$Correlation), nrow(Regulation21P),
                            nrow(Regulation2P), nrow(RegulationP) - nrow(Regulation2P), nrow(Regulation2AP))
      Regulation21Nnum <- c(Module1[j], Module1[i], "Negative",
                            mean(Regulation21N$Correlation), nrow(Regulation21N),
                            nrow(Regulation2N), nrow(RegulationN) - nrow(Regulation2N), nrow(Regulation2AN))
      if (i == j) {
        RegulationNum1 <- rbind(RegulationNum1, Regulation12Pnum, Regulation12Nnum)
      } else {
        RegulationNum1 <- rbind(RegulationNum1, Regulation12Pnum, Regulation12Nnum,
                                Regulation21Pnum, Regulation21Nnum)
      }
    }
  }

  ### perform hypergeometric test
  RegulationNum1 <- as.data.frame(RegulationNum1)
  RegulationNum1 <- RegulationNum1[RegulationNum1$V4 != "NaN", ]
  RegulationP1 <- apply(RegulationNum1[, 5:ncol(RegulationNum1)], 1, function(X1) {
    X1 <- as.numeric(X1)
    if (X1[1] < 4) {
      P1 <- 1
    } else {
      P1 <- phyper(X1[1], X1[2], X1[3], X1[4], lower.tail = FALSE)
    }
  })
  RegulationP2 <- -log10(p.adjust(RegulationP1, method = "fdr"))
  RegulationP3 <- cbind(RegulationNum1[, 1:4], apply(RegulationNum1[, 5:ncol(RegulationNum1)], 1, function(x1) {
    paste(x1, collapse = ";")
  }), RegulationP1, RegulationP2)
  colnames(RegulationP3)[1:7] <- c("TFGroup", "TargetGroup", "Regulation",
                                   "Correlation", "NumberRegulation", "Pvalue", "NlogFdr")

  RegulationP4 <- RegulationP3[as.numeric(as.character(RegulationP3[, "NlogFdr"])) >
                                 -log10(ModuleThr1), ]
  RegulationP4 <- RegulationP4[order(RegulationP4[, "TargetGroup"]), ]
  RegulationP4 <- RegulationP4[order(RegulationP4[, "TFGroup"]), ]
  print(paste("Significant regulations:", nrow(RegulationP4)))
  TFs_list[['intramodular_network']] <- RegulationP4
  return(TFs_list)
}


