#' Collection of immune cell deconvolution methods. 
#'
#' @docType package
#' @name immune_deconvolution_methods
NULL



#' Add together two numbers.
#' 
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
#' @export
deconvolute = function(gene_expression_matrix, method=c("mcp_counter", "epic", "quantiseq", "xcell")) {
  res = switch(method,
         xcell=xCellAnalysis(gene_expression_matrix),
         mcp_counter=MCPcounter.estimate(gene_expression_matrix, "HUGO_symbols"),
         epic=function() {
            epic_res_raw = EPIC(bulk=eset_mat)
            t(epic_res_raw$cellFractions)
         })
  res
}