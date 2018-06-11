#' Collection of immune cell deconvolution methods. 
#'
#' @docType package
#' @name immune_deconvolution_methods
#' @import methods
#' @import dplyr
#' @importFrom testit assert
#' @import readr
NULL


#' list of supported immune deconvolution methods
#' @export
deconvolution_methods = c("mcp_counter", "epic", "quantiseq", "xcell", "cibersort", "cibersort_abs")


#' a mapping of the various cell_types between the different methods to a common standard
#' @export
celltype2method_mapping = readxl::read_xlsx(system.file("extdata", "cell_type_mapping.xlsx", package="immunedeconv", mustWork=TRUE),
                                            sheet="celltype2method")


#' Data object from xCell. 
#' 
#' For some reason, this object is not properly exported from the xCell namespace. 
#' This is a workaround, that `xCellAnalysis` can be properly called from this package. 
#' 
#' @export
xCell.data = xCell::xCell.data

#' set Path to CIBERSORT R script
#'
#' @export
set_cibersort_binary = function(path) {
  assign("cibersort_binary", path, envir=config_env)
}

#' set Path to CIBERSORT matrix file
#' 
#' @export
set_cibersort_mat = function(path) {
  assign("cibersort_mat", path, envir=config_env)
}


#' TODO: documentation
deconvolute_xcell = function(gene_expression_matrix, ...) {
  invisible(capture.output(res <- xCell::xCellAnalysis(gene_expression_matrix, ...)))
  res
}


deconvolute_mcp_counter = function(gene_expression_matrix, feature_types="HUGO_symbols", ...) {
  MCPcounter::MCPcounter.estimate(gene_expression_matrix, featuresType=feature_types, ...)
}


deconvolute_epic = function(gene_expression_matrix, ...) {
  epic_res_raw = EPIC::EPIC(bulk=gene_expression_matrix, ...)
  t(epic_res_raw$cellFractions)
}


deconvolute_quantiseq = function(gene_expresssion_matrix) {
  deconvolute_quantiseq.default(gene_expresssion_matrix) %>% 
    as_tibble() %>%
    select(-Sample) %>%
    as.matrix() %>%
    t()
}

deconvolute_cibersort = function(gene_expression_matrix,
                                 absolute=FALSE,
                                 abs_method="sig.score") {
  assert("CIBERSORT.R is provided", exists("cibersort_binary", envir=config_env))
  assert("CIBERSORT signature matrix is provided", exists("cibersort_mat", envir=config_env))
  source(get("cibersort_binary", envir=config_env))
  
  tmp_mat = tempfile()
  write_tsv(as_tibble(gene_expression_matrix, rownames="gene_symbol"), path=tmp_mat)
  res = CIBERSORT(get("cibersort_mat", envir=config_env), tmp_mat, perm=0, QN=TRUE, absolute=absolute, abs_method=abs_method)
  
  res = res %>% 
    t() %>%
    .[!rownames(.) %in% c("RMSE", "P-value", "Correlation"), ]
  return(res)
}


export_for_timer = function(gene_expression_matrix, path=stop("Specify output path. ")) {
  write_tsv(as_tibble(gene_expression_matrix, rownames="gene_symbol"), path=path)
}


import_from_timer = function(input_file) {
  read_csv(input_file) %>% 
    as.data.frame() %>%
    column_to_rownames("sampleID") %>%
    t() %>%
    as_tibble(rownames="method_cell_type") %>%
    annotate_cell_type("timer")
}

#' Annotate unified cell_type names 
#' 
#' (map the cell_types of the different methods to a common name)
annotate_cell_type = function(result_table, method) {
  celltype2method_mapping %>%
    select(method_cell_type=!!method, cell_type) %>% 
    inner_join(result_table, by="method_cell_type") %>%
    select(-method_cell_type)
}


#' convert an `Biobase::ExpressionSet` to a gene-expression matrix. 
#' 
#' @param eset `ExpressionSet`
#' @param column column name of the `fData()` table, which contains the HGNC gene symbols.  
#' @return matrix with gene symbols as rownames and sample identifiers as colnames. 
eset_to_matrix = function(eset, column) {
  expr_mat = exprs(eset)
  rownames(expr_mat) = fData(eset) %>% pull(!!column)
  return(expr_mat)
}


#' Perform an immune cell deconvolution on a dataset. 
#' 
#' @param gene_expression A numeric matrix with HGNC gene symbols as rownames and sample identifiers as colnames. 
#'   Data must be on non-log scale. 
#' @param method a string specifying the method. Supported methods are `xcell`, `...`
#' @param ... arguments passed to the respective method
#' @return `data.frame` with `cell_type` as first column and a column with the 
#'     calculated cell fractions for each sample. 
#'     
#' @examples
#' # Not run: deconvolute(gene_expression_matrix, "xcell")
#' 
#' @name deconvolute
deconvolute.default = function(gene_expression, method=deconvolution_methods, ...) {
  # run selected method
  res = switch(method,
         xcell = deconvolute_xcell(gene_expression),
         mcp_counter = deconvolute_mcp_counter(gene_expression),
         epic = deconvolute_epic(gene_expression),
         quantiseq = deconvolute_quantiseq(gene_expression),
         cibersort = deconvolute_cibersort(gene_expression, absolute = FALSE),
         cibersort_abs = deconvolute_cibersort(gene_expression, absolute = TRUE))
  
  # convert to tibble and annotate unified cell_type names
  res = res %>%
    as_tibble(rownames="method_cell_type") %>%
    annotate_cell_type(method=method)
  
  return(res)
}


#' @rdname deconvolute
setGeneric("deconvolute", function(gene_expression, method=deconvolution_methods, column="gene_symbol") {
  standardGeneric("deconvolute")
})


#' @describeIn deconvolute `eset` is a `Biobase::ExpressionSet`. 
#' `fData` contains a column with HGNC gene symbols (specify column name)
#' 
#' @param eset Expression set 
#' @param col column name of the `fData` column that contains HGNC gene symbols. 
setMethod("deconvolute", representation(gene_expression="eSet", method="character"),
          function(gene_expression, method, column="gene_symbol") {
            gene_expression %>%
              eset_to_matrix(column) %>% 
              deconvolute.default(method)
          })

#' @describeIn deconvolute `matrix` is matrix with HGNC gene symbols as row names. 
setMethod("deconvolute", representation(gene_expression="matrix", method="character"),
          function(gene_expression, method) {
            deconvolute.default(gene_expression, method)
          })




