deconvolute_abis_default <- function(expression, array = FALSE){
  
  # 2 sets of operations for RNA/Microarray
  
  signature.path <- system.file('extdata', 'abis', package = 'immunedeconv', mustWork=TRUE)
  
  if (array == FALSE){
    signature.path.rna <- paste0(signature.path, '/sigmatrixRNAseq.txt')
    signature <- read.delim(signature.path.rna, check.names = F)
    genes <- intersect(rownames(expression), rownames(signature))
    Dec <- (apply(expression[genes, , drop=F], 2, function(x) coef(rlm( as.matrix(signature[genes,]), x, maxit =100 ))))*100
    
  } else {
    signature.path.micro <- paste0(signature.path, '/sigmatrixMicro.txt')
    signature <- read.delim(signature.path.micro, check.names = F)
    genes <- intersect(rownames(expression), rownames(signature))
    expression2 <- normalize.quantiles.use.target(as.matrix(expression[genes,]), target[,1],copy=TRUE,subset=NULL)
    colnames(expression2) <-colnames(expression)
    rownames(expression2) <- genes
    Dec <- (apply(expression2, 2, function(x) coef(rlm( as.matrix(signature[genes,]), x, maxit =100 )))) *100
  }
  
  Dec <- signif(Dec, 3)
  
  Dec
}