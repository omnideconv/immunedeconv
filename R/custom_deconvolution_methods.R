#' Collection of deconvolution methods that allow custom signature matrices.
#'
#' @description Aliases to some of the deconvolution functions with extra arguments.
#' @import methods
#' @import dplyr
#' @importFrom testit assert
#' @import readr
#' @importFrom tibble as_tibble
#' @importFrom EPIC EPIC
#' @importFrom rlang dots_list
#' @importFrom stats aggregate lm lsfit median qqplot
#' @importFrom utils capture.output read.csv read.table tail write.table
#'
#' @name custom_deconvolution
NULL



#' List of methods that support the use of a custom signature
#'
#' The available methods are
#' `epic`, `cibersort`, `cibersort_abs`, `consensus_tme`, `seqimmucc`
#'
#' The object is a named vector. The names correspond to the display name of the method,
#' the values to the internal name.
#'
#' @export
custom_deconvolution_methods <- c(
  "EPIC" = "epic",
  "CIBERSORT" = "cibersort",
  "ConsensusTME" = "consensus_tme",
  "seqImmuCC" = "seqimmucc"
)


###########################################################################
# Deconvolution with custom signature matrix
#
# These functions let the users specify a custom signature matrix for the analysis
###########################################################################



#' Deconvolute using CIBERSORT or CIBERSORT abs and a custom signature matrix.
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples
#' @param signature_matrix a m x l matrix with m genes and l cell types. The
#'   matrix should contain only a subset of the genes useful for the analysis.
#' @param QN boolean. Wheter to quantile normalize the data. Data should be normalized
#'   when the signature matrix is derived from different studies/sample batches
#' @param absolute Set to TRUE for CIBERSORT absolute mode.
#' @param abs_method Choose method to compute absolute score (only if `absolute=TRUE`).
#' @param ... passed through to the original CIBERSORT function. A native argument takes precedence
#'   over an immunedeconv argument (e.g. `QN` takes precedence over `arrays`). Documentation
#'   is not publicly available. Log in to the CIBERSORT website for details.
#'
#' @note the gene expression and the signature matrix should be provided in the same normalization
#' @export
deconvolute_cibersort_custom <- function(gene_expression_matrix, signature_matrix, QN = FALSE,
                                         absolute = FALSE, abs_method = "sig.score",
                                         ...) {
  assert("CIBERSORT.R is provided", exists("cibersort_binary", envir = config_env))
  source(get("cibersort_binary", envir = config_env))

  temp.expression.file <- tempfile()
  temp.signature.file <- tempfile()
  write_tsv(as_tibble(gene_expression_matrix, rownames = "gene_symbol"), path = temp.expression.file)
  write_tsv(as_tibble(signature_matrix, rownames = "gene_symbol"), path = temp.signature.file)


  arguments <- dots_list(temp.signature.file, temp.expression.file,
    perm = 0,
    QN = QN, absolute = absolute, abs_method = abs_method, ..., .homonyms = "last"
  )

  call <- rlang::call2(CIBERSORT, !!!arguments)

  results <- eval(call)
  results <- results %>%
    .[, !colnames(.) %in% c("RMSE", "P-value", "Correlation")]


  return(t(results))
}

#' Deconvolute using EPIC and a custom signature matrix.
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples
#' @param signature_matrix a m x l matrix with m genes and l cell types. This matrix
#'    should contain the whole set of genes
#' @param signature_genes a character vector of the gene names to use as signature
#'    needs to be smaller than the genes in the signature matrix
#' @param genes_var (optional) a m x l matrix with m genes and l cell types, with
#'    the variability of each gene expression for each cell type.
#'    This will be used in the optimization
#' @param mrna_quantities (optional) A named numeric vector with
#'    the amount of mRNA in arbitrary units for each of the
#'    reference cells and of the other uncharacterized cells.
#' @param ... passed through to EPIC. A native argument takes precedence
#'   over an immunedeconv argument.
#'   See [EPIC](https://rdrr.io/github/GfellerLab/EPIC/man/EPIC.html)
#' @export
deconvolute_epic_custom <- function(gene_expression_matrix, signature_matrix,
                                    signature_genes, genes_var = NULL, mrna_quantities = NULL,
                                    ...) {
  ref <- list()
  ref$refProfiles <- signature_matrix
  ref$sigGenes <- signature_genes
  if (!is.null(genes_var)) {
    ref$refProfiles.var <- genes_var
  }

  mrna_cell <- mrna_quantities
  if (is.null(mrna_quantities)) {
    mRNA_cell <- c("default" = 1.)
  }

  arguments <- dots_list(
    bulk = gene_expression_matrix,
    reference = ref, mRNA_cell = mRNA_cell, ..., .homonyms = "last"
  )

  call <- rlang::call2(EPIC::EPIC, !!!arguments)
  epic_res_raw <- eval(call)

  t(epic_res_raw$cellFractions)
}




#' Deconvolute using ConsesnusTME and a custom signature matrix
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples. Data
#'    should be TPM normalized and log10 scaled.
#' @param signature_genes a list with each element containing genes to represent a cell type. The cell types
#'    should be the names of each element of the list.
#' @param stat_method Choose statistical framework to generate the entichment scores.
#'     Default: 'ssgsea'. Available methods: 'ssgsea', 'gsva', 'plage', 'zscore', 'singScore'.
#'     These mirror the parameter options of \code{GSVA::gsva()} with the exception of \code{singScore}
#'     which leverages \code{singscore::multiScore()}
#' @note ConsensusTME uses tumor-specific consensus built gene signatures. In this case
#'    only the user-provided signature will be used
#' @export
#'
deconvolute_consensus_tme_custom <- function(gene_expression_matrix, signature_genes, stat_method = "ssgsea") {
  results <- ConsensusTME::geneSetEnrichment(
    gene_expression_matrix, signature_genes,
    stat_method
  )

  return(results)
}


#' Deconvolute using seqImmuCC (LLSR regression) and a custom signature matrix
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples. Data
#'    should be TPM normalized and log10 scaled.
#' @param signature_matrix a m x l matrix with m genes and l cell types. The
#'   matrix should contain only a subset of the genes useful for the analysis.
#' @param ... passed to the original seqImmuCC_LLSR fnuction
#' @export
#'
deconvolute_seqimmucc_custom <- function(gene_expression_matrix,
                                         signature_matrix,
                                         ...) {
  results <- seqImmuCC_LLSR(signature_matrix, gene_expression_matrix, ...)
  results <- results[, !colnames(results) %in% c("Correlation", "RMSE")]

  return(t(results))
}

