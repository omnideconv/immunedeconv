#' Collection of deconvolution methods that allow custom signature matrices.
#' 
#' @import methods
#' @import dplyr
#' @importFrom testit assert
#' @import readr
#' @importFrom tibble as_tibble
#' @importFrom EPIC EPIC
#' @importFrom rlang dots_list
#' @importFrom stats aggregate lm lsfit median qqplot
#' @importFrom utils capture.output read.csv read.table tail write.table
NULL



#' List of methods that support the use of a custom signature
#'
#' The available methods are
#' `epic`, `xcell`, `cibersort`, `cibersort_abs`, `consensus_tme`, `base`
#'
#' The object is a named vector. The names correspond to the display name of the method,
#' the values to the internal name.
#'
#' @export
deconvolution_methods = c("EPIC"="epic",
                          "xCell"="xcell",
                          "CIBERSORT"="cibersort",
                          "CIBERSORT (abs.)"="cibersort_abs",
                          "ConsensusTME"="consensus_tme", 
                          "BASE"="base")


###########################################################################
# Deconvolution with custom signature matrix
# 
# These functions let the users specify a custom signature matrix for the analysis 
###########################################################################



#' Deconvolute using CIBERSORT or CIBERSORT abs and a custom signature matrix.
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples
#' @param signature_matrix a m x l matrix with m genes and l cell types
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
deconvolute_cibersort_custom = function(gene_expression_matrix, signature_matrix, QN = FALSE, 
                                        absolute = FALSE, abs_method="sig.score",
                                        ...){
  

  assert("CIBERSORT.R is provided", exists("cibersort_binary", envir=config_env))
  source(get("cibersort_binary", envir=config_env))
  
  temp.expression.file = tempfile()
  temp.signature.file = tempfile()
  write_tsv(as_tibble(gene_expression_matrix, rownames = 'gene_symbol'), path = temp.expression.file)
  write_tsv(as_tibble(signature_matrix, rownames = 'gene_symbol'), path = temp.signature.file)
  
  
  arguments = dots_list(temp.signature.file, temp.expression.file, perm=perm,
                        QN=QN, absolute=absolute, abs_method=abs_method, ..., .homonyms="last")
  
  call = rlang::call2(CIBERSORT, !!!arguments)

  results = eval(call)
  results = results %>%
    .[, !colnames(.) %in% c("RMSE", "P-value", "Correlation")]
  
  
  return(t(results))
  
}

#' Deconvolute using EPIC and a custom signature matrix.
#' 
#' @param gene_expression_matrix
#' @param signature_matrix
#' @param signature_genes
#' @param mrna_cells
#' @param ...
deconvolute_epic_custom = function(gene_expression_matrix, signature_matrix, 
                                   signature_genes, mrna_cells = NULL, 
                                   ...){}




#' Deconvolute using ConsesnusTME and a custom signature matrix
#' 
#' @param gene_expression_matrix a m x n matrix with m genes and n samples. Data
#'    should be TPM normalized and log10 scaled.
#' @param signature_genes a list with each element containing genes to represent a cell type. The cell types
#' should be the names of each element of the list.
#' @param stat_method Choose statistical framework to generate the entichment scores. 
#'     Default: 'ssgsea'. Available methods: 'ssgsea', 'gsva', 'plage', 'zscore', 'singScore'.
#'     These mirror the parameter options of \code{GSVA::gsva()} with the exception of \code{singScore}
#'     which leverages \code{singscore::multiScore()}
#' @note ConsensusTME uses tumor-specific consensus built gene signatures. In this case
#'    only the user-provided signature will be used
#' @export
#' 
deconvolute_consensus_tme_custom = function(gene_expression_matrix, signature_genes, stat_method = 'ssgsea'){
  
  results = ConsensusTME::geneSetEnrichment(gene_expression_matrix, signature_genes, 
                                            stat_method)
  
  return(results)
  
}



#' Deconvolute using BASE and a custom signature matrix
#' 
#' @param gene_expression_matrix a m x n matrix with m genes and n samples. Data
#'    should be TPM normalized and log10 scaled.
#' @param signature_matrix a m x l matrix with m genes and l cell types. Data 
#'    should be non normalized, as the normalization wil be done in the construction 
#'    of the compendium (internal structure)
#' @param n_permutations the number of permutations of each sample expression
#'    to generate. These are used to normalize the results. 
#' @param log10 logical. if TRUE, log10 transforms the expression matrix. 
#' @export
#' 
deconvolute_base_custom = function(gene_expression_matrix, 
                                   signature_matrix, 
                                   n_permutations = 100, 
                                   log10 = TRUE){
  
  new.cell.compendium <- create_base_compendium(signature_matrix)
  results = base_algorithm(gene_expression_matrix, new.cell.compendium, perm = n_permutations)
  
  return(t(results))
}
  
  
  
  



deconvolute_with_custom_signature = function(gene_expression_matrix, method, ...){
  
}
