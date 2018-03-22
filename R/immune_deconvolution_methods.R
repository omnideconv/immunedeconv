#' Collection of immune cell deconvolution methods. 
#'
#' @docType package
#' @name immune_deconvolution_methods
#' @import methods

NULL

#' list of supported immune deconvolution methods
#' @export
deconvolution_methods = c("mcp_counter", "epic", "quantiseq", "xcell")

#' a mapping of the various cell_types between the different methods to a common standard
cell_type_mapping = readxl::read_xlsx(system.file("extdata", "cell_type_mapping.xlsx", package="immunedeconv", mustWork=TRUE))

xCell.data = xCell::xCell.data

#'
deconvolute_xcell = function(gene_expression_matrix, ...) {
  # for some reason, this object cannot be exported from the namespace. 

  xCellAnalysis(gene_expression_matrix, ...)
}

deconvolute_mcp_counter = function(gene_expression_matrix, feature_types="HUGO_symbols", ...) {
  MCPcounter.estimate(gene_expression, feature_types=feature_types, ...)
}

deconvolute_epic = function(gene_expression_matrix, ...) {
  epic_res_raw = EPIC(bulk=gene_expression, ...)
  t(epic_res_raw$cellFractions)
}

#' Add together two numbers.
#' 
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
#' @name deconvolute
deconvolute.default = function(gene_expression, method=deconvolution_methods, ...) {
  # run selected method
  res = switch(method,
         xcell = deconvolute_xcell(gene_expression),
         mcp_counter = deconvolute_mcp_counter(gene_expression),
         epic = deconvolute_epic(gene_expression))
  
  # convert to tibble
  res = res %>% as_tibble(rownames="method_cell_type")
  
  # annotate unified cell_type names
  res = cell_type_mapping %>%
    select(method_cell_type=!!method, cell_type) %>% 
    inner_join(res, by="method_cell_type") %>%
    select(-method_cell_type)
  
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
            expr_mat = exprs(gene_expression)
            rownames(expr_mat) = fData(gene_expression) %>% pull(!!column)
            deconvolute.default(expr_mat, method)
          })

#' @describeIn deconvolute `matrix` is matrix with HGNC gene symbols as row names. 
setMethod("deconvolute", representation(gene_expression="matrix", method="character"),
          function(gene_expression, method) {
            deconvolute.default(gene_expression, method)
          })




