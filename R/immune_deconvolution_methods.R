#' Collection of immune cell deconvolution methods.
#'
#' @description Immunedeconv is an an R package for unified access to computational methods for
#'  estimating immune cell fractions from bulk RNA sequencing data.
#' @docType package
#' @name immunedeconv
#' @import methods
#' @import dplyr
#' @importFrom testit assert
#' @import readr
#' @importFrom tibble as_tibble
#' @importFrom EPIC EPIC
#' @importFrom rlang dots_list
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline barplot box mtext
#' @importFrom stats aggregate lm lsfit median qqplot
#' @importFrom utils capture.output read.csv read.table tail write.table
NULL


#' List of supported immune deconvolution methods
#'
#' The methods currently supported are
#' `mcp_counter`, `epic`, `quantiseq`, `xcell`, `cibersort`, `cibersort_abs`, `timer`, `abis`, `consensus_tme`, `estimate`
#'
#' The object is a named vector. The names correspond to the display name of the method,
#' the values to the internal name.
#'
#' @export
deconvolution_methods <- c(
  "MCPcounter" = "mcp_counter",
  "EPIC" = "epic",
  "quanTIseq" = "quantiseq",
  "xCell" = "xcell",
  "CIBERSORT" = "cibersort",
  "CIBERSORT (abs.)" = "cibersort_abs",
  "TIMER" = "timer",
  "ConsensusTME" = "consensus_tme",
  "ABIS" = "abis",
  "ESTIMATE" = "estimate"
)

#' Data object from xCell.
#'
#' For some reason, this object is not properly exported from the xCell namespace.
#' This is a workaround, that `xCellAnalysis` can be properly called from this package.
#'
#' @export
xCell.data <- xCell::xCell.data

#' Set Path to CIBERSORT R script (`CIBERSORT.R`)
#'
#' CIBERSORT is only freely available to academic users.
#' A license an the binary can be obtained from https://cibersort.stanford.edu.
#'
#' @param path path to cibersort R script.
#'
#' @export
set_cibersort_binary <- function(path) {
  assign("cibersort_binary", path, envir = config_env)
}

#' Set Path to CIBERSORT matrix file (`LM22.txt`)
#'
#' CIBERSORT is only freely available to academic users.
#' A license an the binary can be obtained from https://cibersort.stanford.edu.
#'
#' @param path path to cibersort matrix.
#'
#' @export
set_cibersort_mat <- function(path) {
  assign("cibersort_mat", path, envir = config_env)
}



###########################################################################
# Deconvolution functions for consistently accessing each method
#
# These functions are called from the generic `deconvolute()` function.
# They can also be used by the end-user to access method-specific
# arguments.
###########################################################################


#' Deconvolute using TIMER
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
#' @export
deconvolute_timer <- function(gene_expression_matrix, indications = NULL) {
  indications <- tolower(indications)
  assert("indications fit to mixture matrix", length(indications) == ncol(gene_expression_matrix))
  args <- new.env()
  args$outdir <- tempdir()
  args$batch <- tempfile()
  lapply(unique(indications), function(ind) {
    tmp_file <- tempfile()
    tmp_mat <- gene_expression_matrix[, indications == ind, drop = FALSE] %>% as_tibble(rownames = "gene_symbol")
    write_tsv(tmp_mat, tmp_file)
    cat(paste0(tmp_file, ",", ind, "\n"), file = args$batch, append = TRUE)
  })
  # reorder results to be consistent with input matrix
  results <- deconvolute_timer.default(args)[, make.names(colnames(gene_expression_matrix))]
  colnames(results) <- colnames(gene_expression_matrix)
  results
}


#' Deconvolute using xCell
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples
#' @param arrays Set to TRUE if microarray data, to FALSE for RNASeq (`rnaseq` parameter in xCell)
#' @param expected_cell_types a character list of the cell types to use
#'   in the analysis. If NULL runs xCell with all cell types.
#'   The spillover compensation step may over compensate, thus it is always better
#'   to run xCell with a list of cell types that are expected to be in the mixture.
#'   The names of cell types in this list must be a subset of the cell types that are inferred by xCell.
#'   (`cell.types.use` parameter in xCell)
#' @param ... Passed through to original xCell function. A native argument takes precedence
#'   over an immunedeconv argument (e.g. `rnaseq` takes precedence over `arrays`)
#'   See [xCellAnalysis](https://rdrr.io/github/dviraran/xCell/man/xCellAnalysis.html)
#'
#' @export
deconvolute_xcell <- function(gene_expression_matrix, arrays, expected_cell_types = NULL, ...) {
  rnaseq <- !arrays

  # map the 'expected cell types' to their xCell counterpart.
  if (!is.null(expected_cell_types)) {
    get_children_xcell <- function(cell_type) get_all_children(cell_type, "xcell")
    cell_types_xcell <- lapply(expected_cell_types, get_children_xcell) %>%
      unlist() %>%
      unique()
  } else {
    cell_types_xcell <- NULL
  }

  arguments <- dots_list(gene_expression_matrix,
    rnaseq = rnaseq,
    cell.types.use = cell_types_xcell,
    parallel.sz = config_env$xcell_cores, ..., .homonyms = "last"
  )
  call <- rlang::call2(xCell::xCellAnalysis, !!!arguments)
  invisible(capture.output(res <- eval(call)))

  res
}


#' Deconvolute using MCP-counter
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples
#' @param feature_types type of identifiers used for expression features. May be
#'  one of `"affy133P2_probesets","HUGO_symbols","ENTREZ_ID"`
#' @param ... passed through to original MCP-counter function. A native argument takes precedence
#'   over an immunedeconv argument (e.g. `featureType` takes precedence over `feature_types`)
#'  See [MCPcounter.estimate](https://github.com/ebecht/MCPcounter/blob/master/Source/R/MCPcounter.R#L19).
#'
#' @export
deconvolute_mcp_counter <- function(gene_expression_matrix, feature_types = "HUGO_symbols", ...) {
  arguments <- dots_list(gene_expression_matrix, featuresType = feature_types, ..., .homonyms = "last")
  call <- rlang::call2(MCPcounter::MCPcounter.estimate, !!!arguments)
  eval(call)
}


#' Deconvolute using EPIC
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples
#' @param tumor Set to TRUE if working with tumor data. Will choose the `TRef`
#'  signature matrix in that case, `BRef` otherwise (through EPIC's `reference` parameter)
#' @param scale_mrna Set to FALSE to disable correction for cell type-specific differences
#'  in mRNA content (through EPIC's `mRNA_cell` parameter)
#' @param ... passed through to EPIC. A native argument takes precedence
#'   over an immunedeconv argument (e.g. `ref` takes precedence over `tumor`)
#'  See [EPIC](https://rdrr.io/github/GfellerLab/EPIC/man/EPIC.html)
#'
#' @export
deconvolute_epic <- function(gene_expression_matrix, tumor, scale_mrna, ...) {
  ref <- ifelse(tumor, "TRef", "BRef")
  mRNA_cell <- NULL
  if (!scale_mrna) mRNA_cell <- c("default" = 1.)

  arguments <- dots_list(
    bulk = gene_expression_matrix,
    reference = ref, mRNA_cell = mRNA_cell, ..., .homonyms = "last"
  )
  call <- rlang::call2(EPIC::EPIC, !!!arguments)
  epic_res_raw <- eval(call)

  t(epic_res_raw$cellFractions)
}


#' Deconvolute using quanTIseq
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples
#' @param tumor Set to TRUE if dealing with tumor samples. if TRUE, signature genes with
#'   high expression in tumor samples are removed.
#' @param arrays Set to TRUE if working with Microarray data instead of RNA-seq
#' @param scale_mrna Set to FALSE to disable correction for cell type-specific differences
#'  in mRNA content
#' @param ... passed through to original quantiseq method. A native argument takes precedence
#'   over an immunedeconv argument (e.g. `mRNAscale` takes precedence over `scale_mrna`)
#'   See `deconvolute_quantiseq.default()`.
#'
#' @export
deconvolute_quantiseq <- function(gene_expression_matrix, tumor, arrays, scale_mrna, ...) {
  res <- quantiseqr::run_quantiseq(
    expression_data = gene_expression_matrix,
    is_arraydata = arrays,
    is_tumordata = tumor,
    scale_mRNA = scale_mrna,
    ...
  )

  sample_names <- res$Sample
  res_mat <- res %>%
    as_tibble() %>%
    select(-Sample) %>%
    as.matrix()
  rownames(res_mat) <- sample_names
  t(res_mat)
}


#' Deconvolute using CIBERSORT or CIBERSORT abs.
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples
#' @param arrays Set to TRUE if working with Microarray data instead of RNA-seq.
#'   As recommended by the authors, quantile normalization will be enabled
#'   for Microarrays and disabled for RNAseq.
#' @param absolute Set to TRUE for CIBERSORT absolute mode.
#' @param abs_method Choose method to compute absolute score (only if `absolute=TRUE`).
#' @param ... passed through to the original CIBERSORT function. A native argument takes precedence
#'   over an immunedeconv argument (e.g. `QN` takes precedence over `arrays`). Documentation
#'   is not publicly available. Log in to the CIBERSORT website for details.
#'
#' @export
deconvolute_cibersort <- function(gene_expression_matrix,
                                  arrays,
                                  absolute = FALSE,
                                  abs_method = "sig.score",
                                  ...) {
  # the authors reccomend to disable quantile normalizeation for RNA seq.
  # (see CIBERSORT website).
  quantile_norm <- arrays
  assert("CIBERSORT.R is provided", exists("cibersort_binary", envir = config_env))
  assert("CIBERSORT signature matrix is provided", exists("cibersort_mat", envir = config_env))
  source(get("cibersort_binary", envir = config_env))

  tmp_mat <- tempfile()
  write_tsv(as_tibble(gene_expression_matrix, rownames = "gene_symbol"), file = tmp_mat)

  arguments <- dots_list(get("cibersort_mat", envir = config_env), tmp_mat,
    perm = 0,
    QN = quantile_norm, absolute = absolute, abs_method = abs_method, ..., .homonyms = "last"
  )
  call <- rlang::call2(CIBERSORT, !!!arguments)
  res <- eval(call)

  res <- res %>%
    t() %>%
    .[!rownames(.) %in% c("RMSE", "P-value", "Correlation"), ]
  return(res)
}

#' Deconvolute using ABIS.
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples. Data
#'    must be TPM-normalized, since the matrix accounts for mRNA bias
#' @param arrays Set to TRUE if working with Microarray data instead of RNA-seq.
#'    A different signature matrix will be used.
#' @param ... additional arguments
#'
#' @export
deconvolute_abis <- function(gene_expression_matrix,
                             arrays = FALSE, ...) {
  results <- deconvolute_abis_default(gene_expression_matrix, arrays)
  return(results)
}

#' Deconvolute using ConsensusTME.
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples
#' @param indications a n-vector giving and indication string (e.g. 'brca') for each sample.
#'     The method requires at least 2 samples of a certain cancer type.
#'     Accepted indications are 'kich', 'blca', 'brca', 'cesc', 'gbm', 'hnsc', 'kirp', 'lgg',
#'     'lihc', 'luad', 'lusc', 'prad', 'sarc', 'pcpg', 'paad', 'tgct',
#'     'ucec', 'ov', 'skcm', 'dlbc', 'kirc', 'acc', 'meso', 'thca',
#'     'uvm', 'ucs', 'thym', 'esca', 'stad', 'read', 'coad', 'chol'
#' @param method Choose statistical framework to generate the entichment scores.
#'     Default: 'ssgsea'. Available methods: 'ssgsea', 'gsva', 'plage', 'zscore', 'singScore'.
#'     These mirror the parameter options of \code{GSVA::gsva()} with the exception of \code{singScore}
#'     which leverages \code{singscore::multiScore()}
#' @param ... passed through to the original ConsensusTME function. A native argument takes precedence
#'   over an immunedeconv argument. Documentation can be found at http://consensusTME.org
#'
#' @export
deconvolute_consensus_tme <- function(gene_expression_matrix,
                                      indications = NULL,
                                      method = "ssgsea",
                                      ...) {
  indications <- toupper(indications)
  indications[indications == "UNFILTERED"] <- "Unfiltered"
  assert("indications fit to mixture matrix", length(indications) == ncol(gene_expression_matrix))

  gene_expression_matrix <- as.matrix(gene_expression_matrix)

  sorting <- order(indications)
  indications <- indications[sorting]
  gene_expression_matrix <- gene_expression_matrix[, sorting]


  tumor.types <- unique(indications)
  list.results <- list()
  for (t in tumor.types) {
    cur.samples <- indications == t
    cur.results <- ConsensusTME::consensusTMEAnalysis(
      as.matrix(gene_expression_matrix[, cur.samples]),
      t, method
    )

    list.results[[t]] <- cur.results
  }

  rownames <- Reduce(intersect, lapply(list.results, rownames))
  list.results <- lapply(list.results, function(x) x[rownames, ])
  results <- Reduce(cbind, list.results)

  colnames(results) <- colnames(gene_expression_matrix)

  return(results)
}



#' Annotate unified cell_type names
#'
#' map the cell_types of the different methods to a common name
#'
#' @param result_table output of `deconvolute`
#' @param method one of `immune_deconvolution_methods`.
annotate_cell_type <- function(result_table, method) {
  cell_type_map %>%
    filter(method_dataset == !!method) %>%
    inner_join(result_table, by = "method_cell_type") %>%
    select(-method_cell_type, -method_dataset)
}


#' Convert a `Biobase::ExpressionSet` to a gene-expression matrix.
#'
#' @param eset `ExpressionSet`
#' @param column column name of the `fData()` table, which contains the HGNC gene symbols.
#' @return matrix with gene symbols as rownames and sample identifiers as colnames.
#'
#' @importFrom Biobase exprs fData pData ExpressionSet
#'
#' @export
eset_to_matrix <- function(eset, column) {
  expr_mat <- exprs(eset)
  rownames(expr_mat) <- fData(eset) %>% pull(!!column)
  return(expr_mat)
}


#' Perform an immune cell deconvolution on a dataset.
#'
#' @param gene_expression A gene expression matrix or a Biobase ExpressionSet.
#'   Either: A numeric matrix or data.frame with HGNC gene symbols as rownames and sample identifiers as colnames.
#'   Or: A Biobase ExpressionSet with HGNC symbols in an fData column (see `column` parameter)
#'   In both cases, data must be on non-log scale.
#' @param column Only relevant if `gene_expression` is an ExpressionSet. Defines in which column
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
#' @param expected_cell_types Limit the analysis to the cell types given in this list. If the cell
#'   types present in the sample are known *a priori*, setting this can improve results for
#'   xCell (see https://github.com/grst/immunedeconv/issues/1).
#' @param ... arguments passed to the respective method
#' @return `data.frame` with `cell_type` as first column and a column with the
#'     calculated cell fractions for each sample.
#'
#' @examples
#' # Not run: deconvolute(gene_expression_matrix, "xcell")
#'
#' @name deconvolute
#' @export deconvolute
deconvolute <- function(gene_expression, method = deconvolution_methods,
                        indications = NULL, tumor = TRUE,
                        arrays = FALSE, column = "gene_symbol",
                        rmgenes = NULL, scale_mrna = TRUE,
                        expected_cell_types = NULL,
                        ...) {
  message(paste0("\n", ">>> Running ", method))

  # convert expression set to matrix, if required.
  if (is(gene_expression, "ExpressionSet")) {
    gene_expression <- gene_expression %>% eset_to_matrix(column)
  }

  if (!is.null(rmgenes)) {
    gene_expression <- gene_expression[!rownames(gene_expression) %in% rmgenes, ]
  }

  # run selected method
  res <- switch(method,
    xcell = deconvolute_xcell(gene_expression, arrays = arrays, expected_cell_types = expected_cell_types, ...),
    mcp_counter = deconvolute_mcp_counter(gene_expression, ...),
    epic = deconvolute_epic(gene_expression, tumor = tumor, scale_mrna = scale_mrna, ...),
    quantiseq = deconvolute_quantiseq(gene_expression,
      tumor = tumor, arrays = arrays, scale_mrna = scale_mrna, ...
    ),
    cibersort = deconvolute_cibersort(gene_expression,
      absolute = FALSE,
      arrays = arrays, ...
    ),
    cibersort_abs = deconvolute_cibersort(gene_expression,
      absolute = TRUE,
      arrays = arrays, ...
    ),
    timer = deconvolute_timer(gene_expression, indications = indications, ...),
    abis = deconvolute_abis(gene_expression, arrays = arrays),
    consensus_tme = deconvolute_consensus_tme(gene_expression, indications = indications, ...),
    estimate = deconvolute_estimate(gene_expression)
  )

  # convert to tibble and annotate unified cell_type names
  res <- res %>%
    as_tibble(rownames = "method_cell_type") %>%
    annotate_cell_type(method = method)

  return(res)
}
