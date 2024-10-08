#' Title
#'
#' @param data
#' @param alpha
#' @param boxsize
#' @param weighted
#' 
#' @return
#' @export
#'
#' @examples

Infer_CSN_scNetworks <- function(data, alpha = 0.01, boxsize = 0.1, weighted = 0) {

  if (is.null(weighted)) {
    weighted <- 0
  }
  if (is.null(boxsize)) {
    boxsize <- 0.1
  }
  if (is.null(alpha)) {
    alpha <- 0.01
  }

  n <- nrow(data)
  m <- ncol(data)

  upper <- matrix(0, nrow = n, ncol = m)
  lower <- matrix(0, nrow = n, ncol = m)

  for (i in 1:n) {
    sorted <- sort(data[i, ])  #rank expression values of gene i in all cells
    s1 <- sorted
    s2 <- order(data[i, ])
    n3 <- m - sum(sign(s1))  #calculate the number of cells which have zero expression
    h <- round(boxsize/2*sum(sign(s1)))  #get the threshold according to the boxsize
    k <- 1
    while (k <= m) {
      s <- 0
      while (k + s + 1 <= m && s1[k + s + 1] == s1[k]) {
        s <- s + 1  #determine the number of continuously equal values in gene expression s1
      }
      if (s >= h) {
        upper[i, s2[k:(k + s)]] <- data[i, s2[k]]
        lower[i, s2[k:(k + s)]] <- data[i, s2[k]]  #when the number of continuously equal values in gene expression s1 is more than a threshold, set the lower and upper to the expression of gene k
      } else {
        upper[i, s2[k:(k + s)]] <- data[i, s2[min(m, k + s + h)]]
        lower[i, s2[k:(k + s)]] <- data[i, s2[max(n3 * (n3 > h) + 1, k - h)]]
      }
      k <- k + s + 1
    }
  }

  csn <- list()
  B <- matrix(0, nrow = n, ncol = m)
  p <- -qnorm(alpha, 0, 1)

  for (k in 1:m) {
    for (j in 1:m) {
      B[, j] <- data[, j] <= upper[, k] & data[, j] >= lower[, k]
    }
    a <-  apply(B, 1, sum)
    d <- (B%*%t(B)*m-a%*%t(a))/sqrt((a%*%t(a))*((m-a)%*%t(m-a))/(m-1)+.Machine$double.eps)  #calculate ??
    diag(d) <- 0
    if (weighted == 1) {
      csn[[k]] <- d * (d > 0)
    } else {
      csn[[k]] <- d > p
      csn[[k]][csn[[k]]] <- 1
    }
    gene_names <- rownames(data)
    rownames(csn[[k]]) <- gene_names
    colnames(csn[[k]]) <- gene_names
    print(paste("Cell", k, "is completed"))
  }
  names(csn) <- colnames(data)

  return(csn)
}
