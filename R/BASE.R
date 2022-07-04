#' Souce code for the BASE algorithm.
#' This code is adapted from Varn et al., DOI: 10.1158/0008-5472.CAN-16-2490
#'
#' @param data: numeric matrix; Patient gene expression data, with genes as rows and patients as columns
#' @param reg: numeric matrix; Reference immune cell weights created using SuppSoftware1
#' @param perm: numeric; Number of permutations to perform during the normalization step
#' @param median.norm: logical;  parameter stating whether to median center the expression values for each gene
#'
#'
base_algorithm <- function(data, reg, perm = 100, median.norm = T) {
  ## quantile normalization
  myrk <- matrix(0, nrow(data), ncol(data))
  xx <- myrk
  for (k in 1:ncol(data))
  {
    myrk[, k] <- rank(data[, k])
    xx[, k] <- sort(data[, k])
  }
  mymed <- apply(xx, 1, median, na.rm = T)
  for (k in 1:ncol(data))
  {
    data[, k] <- mymed[myrk[, k]]
  }

  comGene <- intersect(row.names(data), row.names(reg))
  data <- data[comGene, ]
  reg <- reg[comGene, ]

  if (median.norm) {
    mymed <- apply(data, 1, median)
    data <- data - mymed
  }

  cnum <- ncol(data)
  rnum <- nrow(data)
  es <- matrix(0, cnum, ncol(reg))
  cat("\ncalculate ES\n")
  for (k in 1:cnum)
  {
    cat("\r", k)
    myorder <- order(data[, k], decreasing = T)
    cur.exp <- data[myorder, k]
    cur.reg <- reg[myorder, ]
    fg1 <- as.matrix(abs(cur.reg * cur.exp))
    bg1 <- as.matrix(abs((1 - cur.reg) * cur.exp))
    for (i in 2:nrow(fg1))
    {
      fg1[i, ] <- fg1[i, ] + fg1[i - 1, ]
      bg1[i, ] <- bg1[i, ] + bg1[i - 1, ]
    }
    for (i in 1:ncol(fg1))
    {
      fg1[, i] <- fg1[, i] / fg1[rnum, i]
      bg1[, i] <- bg1[, i] / bg1[rnum, i]
    }
    xx <- fg1 - bg1
    tmp <- apply(xx, 2, max)
    pos.es <- ifelse(tmp > 0, tmp, 0)
    tmp <- apply(xx, 2, min)
    neg.es <- ifelse(tmp < 0, tmp, 0)
    es[k, ] <- ifelse(pos.es > abs(neg.es), pos.es, neg.es)
  }

  ## perm
  cat("\n Permutation \n")
  cur.reg <- reg

  pos.es <- neg.es <- matrix(0, ncol(reg), perm)
  for (k in 1:perm)
  {
    cat("\r", k)
    se <- sample(1:cnum, 1)
    cur.exp <- sample(data[, se])
    fg1 <- as.matrix(abs(cur.reg * cur.exp))
    bg1 <- as.matrix(abs((1 - cur.reg) * cur.exp))
    for (i in 2:nrow(fg1))
    {
      fg1[i, ] <- fg1[i, ] + fg1[i - 1, ]
      bg1[i, ] <- bg1[i, ] + bg1[i - 1, ]
    }
    for (i in 1:ncol(fg1))
    {
      fg1[, i] <- fg1[, i] / fg1[rnum, i]
      bg1[, i] <- bg1[, i] / bg1[rnum, i]
    }
    xx <- fg1 - bg1
    tmp <- apply(xx, 2, max)
    pos.es[, k] <- ifelse(tmp > 0, tmp, 0)
    tmp <- apply(xx, 2, min)
    neg.es[, k] <- ifelse(tmp < 0, tmp, 0)
  }

  ## normalize
  pavg <- apply(pos.es, 1, mean)
  navg <- abs(apply(neg.es, 1, mean))
  pos.npes <- pos.es / pavg
  neg.npes <- neg.es / navg
  for (k in 1:nrow(es))
  {
    tmp <- es[k, ]
    es[k, ] <- ifelse(tmp > 0, tmp / pavg, tmp / navg)
  }
  res <- es
  colnames(res) <- paste(colnames(reg), ".ES", sep = "")
  row.names(res) <- colnames(data)

  # Now we have to adjust the results
  ncol.half <- ncol(res) / 2

  CLP.up <- res[, 1:ncol.half]
  CLP.down <- res[, (ncol.half + 1):ncol(res)]

  CLP.scores <- CLP.up - CLP.down
  colnames(CLP.scores) <- str_replace_all(
    colnames(CLP.scores),
    "\\_UP\\.ES", ""
  )
  colnames(CLP.scores) <- make.names(colnames(CLP.scores), allow_ = FALSE)

  CLP.scores <- CLP.scores[, is.finite(colSums(CLP.scores))]

  return(CLP.scores)
}



#' Souce code to create the compendium used kin the BASE algorithm, containing
#'    up- and down-regulated weight sets that specify the
# 	    specificity by which each gene is expressed by a given cell.
#' This code is adapted from Varn et al., DOI: 10.1158/0008-5472.CAN-16-2490
#'
#' @param signature_matrix: numeric matrix; The signature matrix from which the compendium will be built.
#'    Must contain genes on rows and cell on columns
#'
create_base_compendium <- function(signature_matrix) {
  myinf1 <- signature_matrix

  med <- apply(myinf1, 1, median)
  myinf1 <- myinf1 - med

  avg <- apply(myinf1, 2, mean)
  std <- apply(myinf1, 2, sd)
  for (k in 1:ncol(myinf1))
  {
    myinf1[, k] <- (myinf1[, k] - avg[k]) / std[k]
  }


  res1 <- myinf1
  for (k in 1:ncol(res1))
  {
    tmp <- myinf1[, k]
    tmp[tmp < 0] <- 0
    tmp <- -log10(pnorm(-tmp) * 2)
    tmp[tmp > 10] <- 10
    res1[, k] <- tmp
  }
  colnames(res1) <- paste(colnames(myinf1), "_up", sep = "")

  res2 <- myinf1
  for (k in 1:ncol(res2))
  {
    tmp <- myinf1[, k]
    tmp[tmp > 0] <- 0
    tmp <- -log10(pnorm(tmp) * 2)
    tmp[tmp > 10] <- 10
    res2[, k] <- tmp
  }
  colnames(res2) <- paste(colnames(myinf1), "_dn", sep = "")

  res <- cbind(res1, res2)

  minv <- min(res)
  maxv <- max(res)
  res <- (res - minv) / (maxv - minv)
  colnames(res) <- gsub(" ", "", colnames(res))


  colnames(res) <- toupper(colnames(res))

  return(res)
}
