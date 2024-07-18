function (RNA1, K1 = 1, ColumnGroup1 = NULL, Scale1 = "row",
    Range1 = c(-Inf, Inf), Reorder1 = TRUE, RevOrder1 = -1, NAcolum1 = NULL)
{
    validInput(K1, "K1", "numeric")
    validInput(Range1, "Range1", "numeric")
    RowGroup1 = NULL
    NumColumnBlank1 = NULL
    if (is.numeric(K1)) {
        K2 <- K1
    }
    else {
        K2 <- nrow(K1)
    }
    if (Scale1 == "row") {
        RNA1 <- t(scale(t(RNA1)))
    }
    RNA1[is.na(RNA1)] = 0
    RNA1[is.nan(RNA1)] = 0
    if (!is.null(ColumnGroup1)) {
        print("sepate columns according to varible ColumnGroup1")
        ColumnBlank1 <- array(0, dim = c(nrow(RNA1), NumColumnBlank1))
        uColumnGroup1 <- unique(ColumnGroup1)
        if (length(uColumnGroup1) == 1) {
        }
        else {
            for (i in 1:length(uColumnGroup1)) {
                RNA12 <- RNA1[, ColumnGroup1 == uColumnGroup1[i]]
                if (i == 1) {
                  RNA21 <- cbind(RNA12, ColumnBlank1)
                }
                else if (i < length(uColumnGroup1)) {
                  RNA21 <- cbind(RNA21, RNA12, ColumnBlank1)
                }
                else {
                  RNA21 <- cbind(RNA21, RNA12)
                }
            }
        }
    }
    else {
        RNA21 <- RNA1
    }
    print("Perform k-means")
    if (is.null(RowGroup1)) {
        if (!is.null(NAcolum1)) {
            cRNA1 <- kmeans(RNA1[, -NAcolum1], K1)
            RNA02 <- cbind(cRNA1$cluster, RNA1[, -NAcolum1])
            Cluster1 <- cRNA1$cluster
        }
        else {
            if (is.numeric(K1)) {
                if (K1 == 1) {
                  Cluster1 <- rep(1, nrow(RNA1))
                }
                else {
                  cRNA1 <- kmeans(RNA1, K1)
                  Cluster1 <- cRNA1$cluster
                }
            }
            else {
                cRNA1 <- kmeans(RNA1, K1)
                Cluster1 <- cRNA1$cluster
            }
            RNA02 <- cbind(Cluster1, RNA1)
        }
        RNA2 <- cbind(Cluster1, RNA21)
        NameInd1 <- 1:K2
        colnames(RNA2)[1] <- c("KmeansGroup")
        print(table(RNA2[, "KmeansGroup"]))
    }
    else {
        if (is.factor(RowGroup1)) {
            Name1 <- levels(RowGroup1)
            Name2 <- sort(levels(RowGroup1))
            NameInd1 <- match(Name1, Name2)
            RowGroup1 <- as.numeric(RowGroup1)
            uRowGroup1 <- NameInd1
        }
        else {
            uRowGroup1 <- sort(unique(RowGroup1))
        }
        if (!is.null(NAcolum1)) {
            RNA02 <- cbind(RowGroup1, RNA1[, -NAcolum1])
        }
        else {
            RNA02 <- cbind(RowGroup1, RNA1)
        }
        RNA2 <- cbind(RowGroup1, RNA21)
        K1 <- length(uRowGroup1)
        colnames(RNA2)[1] <- c("KmeansGroup")
        tRNA2 <- table(RNA2[, "KmeansGroup"])
        if (is.factor(RowGroup1)) {
            names(tRNA2) <- Name2
            tRNA2 <- tRNA2[NameInd1]
        }
        print(tRNA2)
    }
    print("Sort genes")
    RNA22 <- c()
    for (i in 1:K2) {
        if (is.null(RowGroup1)) {
            RNA20 <- RNA2[RNA2[, "KmeansGroup"] == i, ]
            RNA03 <- RNA02[RNA02[, 1] == i, ]
        }
        else {
            RNA20 <- RNA2[RNA2[, "KmeansGroup"] == uRowGroup1[i],
                ]
            RNA03 <- RNA02[RNA02[, 1] == uRowGroup1[i], ]
        }
        if (Reorder1 == TRUE) {
            if (nrow(RNA03) > 1) {
                Hier1 <- hclust(as.dist((1 - cor(t(RNA03[, 2:ncol(RNA03)])))/2))
                if (RevOrder1[1] != -1) {
                  RevOrder2 <- FALSE
                  for (j in 1:length(RevOrder1)) {
                    if (RevOrder1[j] == i) {
                      RevOrder2 <- TRUE
                      break
                    }
                  }
                  if (RevOrder2 == TRUE) {
                    Ind1 <- rev(Hier1$order)
                  }
                  else {
                    Ind1 <- Hier1$order
                  }
                }
                else {
                  Ind1 <- Hier1$order
                }
            }
            else {
                Ind1 <- 1:nrow(RNA20)
            }
        }
        else {
            Ind1 <- 1:nrow(RNA20)
        }
        RNA22 <- rbind(RNA22, RNA20[Ind1, ])
    }
    print("Revise outlier")
    RNA23 <- RNA22[, 2:ncol(RNA22)]
    print(paste("Number of outlier:", c(length(RNA23[RNA23 <
        Range1[1]]), length(RNA23[RNA23 > Range1[2]]))))
    RNA23[RNA23 < Range1[1]] <- Range1[1]
    RNA23[RNA23 > Range1[2]] <- Range1[2]
    RNA3 <- cbind(RNA22[, 1], RNA23)
    colnames(RNA3)[1] <- "Module"
    RNA4 <- RNA22[RNA22[, 1] != 0, ]
    RNA5 <- cbind(RNA4[match(rownames(RNA1), rownames(RNA4)),
        "KmeansGroup"], RNA1)
    colnames(RNA5)[1] <- "KmeansGroup"
    RNA6 <- RNA5[order(RNA5[, "KmeansGroup"]), ]
    RNA6 <- as.data.frame(RNA6)
    return(RNA6)
}
