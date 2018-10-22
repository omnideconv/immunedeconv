#' Collection of immune cell deconvolution methods.
#'
#' @docType package
#' @name immune_deconvolution_methods
#' @import methods
#' @import dplyr
#' @importFrom testit assert
#' @import readr
#' @importFrom tibble as_tibble
NULL


#' List of supported immune deconvolution methods
#'
#' The methods currently supported are
#' `mcp_counter`, `epic`, `quantiseq`, `xcell`, `cibersort`, `cibersort_abs`, `timer`
#'
#' The object is a named vector. The names correspond to the display name of the method,
#' the values to the internal name.
#'
#' @export
deconvolution_methods = c("MCPcounter"="mcp_counter",
                          "EPIC"="epic",
                          "quanTIseq"="quantiseq",
                          "xCell"="xcell",
                          "CIBERSORT"="cibersort",
                          "CIBERSORT (abs.)"="cibersort_abs",
                          "TIMER"="timer",
                          "random pred."="random")



#' Deconvolute using the awseome RANDOM technique
#'
#' Here is a good place to add some documentation.
deconvolute_random = function(gene_expression_matrix) {
  # list of the cell types we want to 'predict'
  cell_types = c("CD4+ Tcell", "CD8+ Tcell", "NK cell", "Macrophage",
                 "Monocyte")
  n_samples = ncol(gene_expression_matrix)

  # generate random values
  results = matrix(runif(length(cell_types) * n_samples), ncol=n_samples)

  # rescale the values to sum to 1 for each sample
  results = apply(results, 2, function(x) {x/sum(x)})
  rownames(results) = cell_types
  colnames(results) = colnames(gene_expression_matrix)

  results
}

#' Data object from xCell.
#'
#' For some reason, this object is not properly exported from the xCell namespace.
#' This is a workaround, that `xCellAnalysis` can be properly called from this package.
#'
#' @export
xCell.data = xCell::xCell.data

#' Set Path to CIBERSORT R script (`CIBERSORT.R`)
#'
#' CIBERSORT is only freely available to academic users.
#' A license an the binary can be obtained from https://cibersort.stanford.edu.
#'
#' @export
set_cibersort_binary = function(path) {
  assign("cibersort_binary", path, envir=config_env)
}

#' Set Path to CIBERSORT matrix file (`LM22.txt`)
#'
#' CIBERSORT is only freely available to academic users.
#' A license an the binary can be obtained from https://cibersort.stanford.edu.
#'
#' @export
set_cibersort_mat = function(path) {
  assign("cibersort_mat", path, envir=config_env)
}



###########################################################################
# Deconvolution functions for consistenctly accessing each method
#
# These functions are called from the generic `deconvolute()` function
###########################################################################


#' Deconvolute using the TIMER technique
#'
#' Unlike the other methods, TIMER needs the specification of the
#' cancer type for each sample.
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples
#' @param indications a n-vector giving and indication string (e.g. 'brca') for each sample.
#'     Accepted indications are 'kich', 'blca', 'brca', 'cesc', 'gbm', 'hnsc', 'kirp', 'lgg',
#'     'lihc', 'luad', 'lusc', 'prad', 'sarc', 'pcpg', 'paad', 'tgct',
#'     'ucec', 'ov', 'skcm', 'dlbc', 'kirc', 'acc', 'meso', 'thca',
#'     'uvm', 'ucs', 'thym', 'esca', 'stad', 'read', 'coad', 'chol'
#' @param tumor ignored for this method
#' @param arrays ignored for this method
deconvolute_timer = function(gene_expression_matrix, indications=NULL) {
  indications = tolower(indications)
  assert("indications fit to mixture matrix", length(indications) == ncol(gene_expression_matrix))
  args = new.env()
  args$outdir = tempdir()
  args$batch = tempfile()
  lapply(unique(indications), function(ind) {
    tmp_file = tempfile()
    tmp_mat = gene_expression_matrix[, indications == ind, drop=FALSE] %>% as_tibble(rownames = "gene_symbol")
    write_tsv(tmp_mat, tmp_file)
    cat(paste0(tmp_file, ",", ind, "\n"), file=args$batch, append=TRUE)
  })
  # reorder results to be consistent with input matrix
  results = deconvolute_timer.default(args)[, make.names(colnames(gene_expression_matrix))]
  colnames(results) = colnames(gene_expression_matrix)
  results
}


deconvolute_xcell = function(gene_expression_matrix, arrays, ...) {
  rnaseq = !arrays
  invisible(capture.output(res <- xCell::xCellAnalysis(gene_expression_matrix, rnaseq=rnaseq, ...)))
  res
}


deconvolute_mcp_counter = function(gene_expression_matrix, feature_types="HUGO_symbols", ...) {
  MCPcounter::MCPcounter.estimate(gene_expression_matrix, featuresType=feature_types, ...)
}


deconvolute_epic = function(gene_expression_matrix, tumor, scale_mrna, ...) {
  ref = ifelse(tumor, "TRef", "BRef")
  mRNA_cell = NULL
  if(!scale_mrna) mRNA_cell = c("default"=1.)
  epic_res_raw = EPIC::EPIC(bulk=gene_expression_matrix,
                            reference=ref, mRNA_cell = mRNA_cell, ...)
  t(epic_res_raw$cellFractions)
}


deconvolute_quantiseq = function(gene_expresssion_matrix, tumor, arrays, scale_mrna) {
  deconvolute_quantiseq.default(gene_expresssion_matrix, tumor=tumor, arrays=arrays, mRNAscale = scale_mrna) %>%
    as_tibble() %>%
    select(-Sample) %>%
    as.matrix() %>%
    t()
}

deconvolute_cibersort = function(gene_expression_matrix,
                                 arrays,
                                 absolute=FALSE,
                                 abs_method="sig.score") {
  # the authors reccomend to disable quantile normalizeation for RNA seq.
  # (see CIBERSORT website).
  quantile_norm = arrays
  assert("CIBERSORT.R is provided", exists("cibersort_binary", envir=config_env))
  assert("CIBERSORT signature matrix is provided", exists("cibersort_mat", envir=config_env))
  source(get("cibersort_binary", envir=config_env))

  tmp_mat = tempfile()
  write_tsv(as_tibble(gene_expression_matrix, rownames="gene_symbol"), path=tmp_mat)
  res = CIBERSORT(get("cibersort_mat", envir=config_env), tmp_mat, perm=0,
                  QN=quantile_norm, absolute=absolute, abs_method=abs_method)

  res = res %>%
    t() %>%
    .[!rownames(.) %in% c("RMSE", "P-value", "Correlation"), ]
  return(res)
}


#' Annotate unified cell_type names
#'
#' (map the cell_types of the different methods to a common name)
annotate_cell_type = function(result_table, method) {
  cell_type_map %>%
    filter(method_dataset == !!method) %>%
    inner_join(result_table, by="method_cell_type") %>%
    select(-method_cell_type, -method_dataset)
}


#' Convert a `Biobase::ExpressionSet` to a gene-expression matrix.
#'
#' @param eset `ExpressionSet`
#' @param column column name of the `fData()` table, which contains the HGNC gene symbols.
#' @return matrix with gene symbols as rownames and sample identifiers as colnames.
#'
#' @export
eset_to_matrix = function(eset, column) {
  expr_mat = exprs(eset)
  rownames(expr_mat) = fData(eset) %>% pull(!!column)
  return(expr_mat)
}


#' Perform an immune cell deconvolution on a dataset.
#'
#' @param gene_expression A gene expression matrix or a Biobase ExpressionSet.
#'   Either: A numeric matrix or data.frame with HGNC gene symbols as rownames and sample identifiers as colnames.
#'   Or: A Biobase ExpressionSet with HGNC symbols in an fData column (see `column` parameter)
#'   In both cases, data must be on non-log scale.
#' @param column Only releveant if `gene_expression` is an ExpressionSet. Defines in which column
#'   of fData the HGNC symbol can be found.
#' @param method a string specifying the method.
#'   Supported methods are `xcell`, `...`
#' @param indications a character vector with one indication per
#'   sample for TIMER. Argument is ignored for all other methods.
#' @param tumor use a signature matrix/procedure optimized for tumor samples,
#'   if supported by the method. Currently affects EPIC and quanTIseq.
#' @param arrays Runs methods in a mode optimized for microarray data.
#'   Currently affects quanTIseq and CIBERSORT.
#' @param rmgenes a character vector of gene symbols. Exclude these genes from the analysis.
#'   Use this to exclude e.g. noisy genes.
#' @param scale_mrna logical. If FALSE, disable correction for mRNA content of different cell types.
#'   This is supported by methods that compute an absolute score (EPIC and quanTIseq)
#' @param ... arguments passed to the respective method
#' @return `data.frame` with `cell_type` as first column and a column with the
#'     calculated cell fractions for each sample.
#'
#' @examples
#' # Not run: deconvolute(gene_expression_matrix, "xcell")
#'
#' @name deconvolute
#' @export deconvolute
deconvolute = function(gene_expression, method=deconvolution_methods,
                       indications=NULL, tumor=TRUE,
                       arrays=FALSE, column="gene_symbol",
                       rmgenes=NULL, scale_mrna=TRUE,
                       ...) {
  message(paste0("\n", ">>> Running ", method))

  # convert expression set to matrix, if required.
  if(is(gene_expression, "ExpressionSet")) {
    gene_expression = gene_expression %>% eset_to_matrix(column)
  }

  if(!is.null(rmgenes)) {
    gene_expression = gene_expression[!rownames(gene_expression) %in% rmgenes,]
  }

  # run selected method
  res = switch(method,
         xcell = deconvolute_xcell(gene_expression, arrays=arrays, ...),
         mcp_counter = deconvolute_mcp_counter(gene_expression, ...),
         epic = deconvolute_epic(gene_expression, tumor=tumor, scale_mrna=scale_mrna, ...),
         quantiseq = deconvolute_quantiseq(gene_expression,
                                           tumor=tumor, arrays=arrays, scale_mrna=scale_mrna, ...),
         cibersort = deconvolute_cibersort(gene_expression, absolute = FALSE,
                                           arrays=arrays, ...),
         cibersort_abs = deconvolute_cibersort(gene_expression, absolute = TRUE,
                                               arrays=arrays, ...),
         random = deconvolute_random(gene_expression),
         timer = deconvolute_timer(gene_expression, indications=indications, ...))

  # convert to tibble and annotate unified cell_type names
  res = res %>%
    as_tibble(rownames="method_cell_type") %>%
    annotate_cell_type(method=method)

  return(res)
}

