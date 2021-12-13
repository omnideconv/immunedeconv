#' Source code for the Linear Least Square Regression performed in seqImmuCC.
#'
#' This code is provided by the authors and the method available at:
#' http://218.4.234.74:3200/immune/
#'
#' The method is described in Chen et al. doi: 10.3389/fimmu.2018.01286
#' 
#' @importFrom preprocessCore normalize.quantiles
#' 
#' @param signature signature matrix
#' @param SampleData sample expression profile
#' @param w vector to perform weighted least squares
#' @param QN logical. If TRUE expression matrix is quantile normalized
#' @param sig.scale logical. If TRUE, scales the signature matrix
#' @param sig.stand logical.If TRUE, standardizes the signature matrix
#' @param sample.scale logical. If TRUE, scales the sample expression 
#' @param log logical. If TRUE, log transforms signature and expression data
LLSR <- function(signature, SampleData, w=NA, QN=T, sig.scale=F, sig.stand=T, sample.scale=T, log=T){
 
  # Expression profile format standarlization
  signature <- data.matrix(signature)
  SampleData <- data.matrix(SampleData)
  
  # Make sure the data is non log transformed
  if(log==F){
    if (max(signature) < 50) {signature <- 2^signature}   
    if (max(SampleData) < 50) {SampleData <- 2^SampleData}
  }else{
    if (max(signature) > 50) {signature <- log2(signature+1)}   
    if (max(SampleData) > 50) {SampleData <- log2(SampleData+1)}
  }
  
  
  #standardize sig matrix
  if (sig.stand == T) signature <- (signature - mean(signature)) / sd(as.vector(signature)) 
  
  # scale of the signature data
  if (sig.scale == T) {
    signature <- t(scale(t(signature)))
  } else {
    signature <- signature
  }
  
  
  # Quantile normalization of the sample data
  library(preprocessCore)
  if (QN == T) {
    tmpc <- colnames(SampleData)
    tmpr <- rownames(SampleData)
    data <- normalize.quantiles(SampleData)
    colnames(SampleData) <- tmpc
    rownames(SampleData) <- tmpr
  } 
  
  common.gene <- intersect(rownames(signature), rownames(SampleData))
  signature <- signature[common.gene, ]
  SampleData <- SampleData[common.gene, ]
  
  N1 <- ncol(SampleData)
  fraction <- c()
  rsem <- c()
  corrv <- c()
  for (i in seq(N1)) {
    y <- SampleData[, i]
    if(sample.scale==T) y <- (y - mean(y)) / sd(y)
    if (is.na(w[1])) tmp <- lsfit(signature, y, intercept=F) else tmp <- lsfit(signature, y, w, intercept=F)
    tmp.res <- tmp$residuals
    tmp.coef <- tmp$coefficients
    tmp.coef[which(tmp.coef <0)] <- 0
    cat(length(tmp.coef), "\n")
    
    u <- sweep(signature, MARGIN=2, tmp$coefficients, '*')
    v <- apply(u, 1, sum)
    rsem <- c(rsem, sqrt(mean((v-SampleData[, i])^2)))
    cat ("The RSEM value is", rsem[i], "\n")
    corrv <- c(corrv, cor(v, SampleData[, i]))
    
    tmp.coef <- tmp.coef/sum(tmp.coef) *100
    fraction <- rbind(fraction, tmp.coef)
  }
  fraction <- cbind(fraction, corrv, rsem)
  
  colnames(fraction) <- c(colnames(signature), "Correlation", "RSEM")
  rownames(fraction) <- colnames(SampleData)
  data.matrix(fraction)
}