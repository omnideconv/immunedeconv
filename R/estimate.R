#'  Source code for the ESTIMATE algorithm:
#'  Estimate of Stromal and Immune Cells in Malignant Tumor Tissues from Expression Data
#'  (https://doi.org/10.1038/ncomms3612)
#'  Source: http://r-forge.r-project.org/projects/estimate/
#'  Copyright: 2013-2022, MD Anderson Cancer Center (MDACC)
#'  License: GPL-2
#'
#'  Source code adapted from: http://download.r-forge.r-project.org/bin/windows/contrib/3.3/estimate_1.0.13.zip
#'



#'  ESTIMATE algortihm: This function computes stromal, immune, and ESTIMATE scores
#'  per sample using gene expression data, through GSEA. The ESTIMATE score is used
#'  to compute an estimate of the tumor purity
#'  *Warning*: The tumor purity estimation was originally intended only for microarray
#'  affymetrix data
#'
#'  @param gene_expression_matrix a m x n matrix with m genes and n samples
#'
#'  @export
deconvolute_estimate <- function(gene_expression_matrix){

  estimate.genes.path <- system.file('extdata', 'estimate', 'estimate_genes.rds',
                                         package = 'immunedeconv', mustWork=TRUE)

  estimate.files <- readRDS(estimate.genes.path)

  common_genes <- estimate.files$common_genes
  signature_geneset <- estimate.files$signature_geneset


  merged.df <- merge(common_genes, gene_expression_matrix, by.x = "GeneSymbol", by.y = "row.names")
  rownames(merged.df) <- merged.df$GeneSymbol
  filtered_matrix <- merged.df[, -1:-(ncol(common_genes))]
  print(sprintf("Merged dataset includes %d genes (%d mismatched).",
                nrow(filtered_matrix), nrow(common_genes) - nrow(filtered_matrix)))


  m <- filtered_matrix
  gene.names <- rownames(m)
  sample.names <- colnames(m)

  Ns <- length(m[1, ])
  Ng <- length(m[, 1])


  for (j in 1:Ns) {
    m[, j] <- rank(m[, j], ties.method = "average")
  }
  m <- 10000 * m/Ng
  gs <- as.matrix(signature_geneset[, -1], dimnames = NULL)
  N.gs <- 2
  gs.names <- row.names(signature_geneset)
  score.matrix <- matrix(0, nrow = N.gs, ncol = Ns)
  for (gs.i in 1:N.gs) {
    gene.set <- gs[gs.i, ]
    gene.overlap <- intersect(gene.set, gene.names)
    print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=",
                length(gene.overlap)))
    if (length(gene.overlap) == 0) {
      score.matrix[gs.i, ] <- rep(NA, Ns)
      next
    }
    else {
      ES.vector <- vector(length = Ns)
      for (S.index in 1:Ns) {
        gene.list <- order(m[, S.index], decreasing = TRUE)
        gene.set2 <- match(gene.overlap, gene.names)
        correl.vector <- m[gene.list, S.index]
        TAG <- sign(match(gene.list, gene.set2, nomatch = 0))
        no.TAG <- 1 - TAG
        N <- length(gene.list)
        Nh <- length(gene.set2)
        Nm <- N - Nh
        correl.vector <- abs(correl.vector)^0.25
        sum.correl <- sum(correl.vector[TAG == 1])
        P0 <- no.TAG/Nm
        F0 <- cumsum(P0)
        Pn <- TAG * correl.vector/sum.correl
        Fn <- cumsum(Pn)
        RES <- Fn - F0
        max.ES <- max(RES)
        min.ES <- min(RES)
        if (max.ES > -min.ES) {
          arg.ES <- which.max(RES)
        }
        else {
          arg.ES <- which.min(RES)
        }
        ES <- sum(RES)
        EnrichmentScore <- list(ES = ES, arg.ES = arg.ES,
                                RES = RES, indicator = TAG)
        ES.vector[S.index] <- EnrichmentScore$ES
      }
      score.matrix[gs.i, ] <- ES.vector
    }
  }
  score.data <- data.frame(score.matrix)
  names(score.data) <- sample.names
  row.names(score.data) <- gs.names
  estimate.score <- apply(score.data, 2, sum)

  convert_row_estimate_score_to_tumor_purity <- function(x) {
    stopifnot(is.numeric(x))
    cos(0.6049872018 + 0.0001467884 * x)
  }
  est.new <- NULL
  for (i in 1:length(estimate.score)){
    est_i <- convert_row_estimate_score_to_tumor_purity(estimate.score[i])
    est.new <- rbind(est.new, est_i)
    if (est_i >= 0) {
      next
    } else {
     message(paste(sample.names[i], ": out of bounds",
                      sep = ""))}
  }

  colnames(est.new) <- c("TumorPurity")
  estimate.t1 <- cbind(estimate.score, est.new)
  x.bad.tumor.purities <- estimate.t1[, "TumorPurity"] < 0
  estimate.t1[x.bad.tumor.purities, "TumorPurity"] <- NA
  score.data <- rbind(score.data, t(estimate.t1))
  rownames(score.data) <- c("StromalScore", "ImmuneScore",
                            "ESTIMATEScore", "TumorPurity")

  score.data
}
