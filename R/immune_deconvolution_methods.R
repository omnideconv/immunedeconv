#' Collection of immune cell deconvolution methods. 
#'
#' @docType package
#' @name immune_deconvolution_methods
NULL

#' list of supported immune deconvolution methods
#' @export
deconvolution_methods = c("mcp_counter", "epic", "quantiseq", "xcell")

#' a mapping of the various cell_types between the different methods to a common standard
cell_type_mapping = read_xlsx(system.file("extdata", "cell_type_mapping.xlsx", package="immunedeconv", mustWork=TRUE))

#' Add together two numbers.
#' 
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
#' @export
deconvolute = function(gene_expression_matrix, method=deconvolution_methods) {
  res = switch(method,
         xcell=xCellAnalysis(gene_expression_matrix),
         mcp_counter=MCPcounter.estimate(gene_expression_matrix, "HUGO_symbols"),
         epic=function() {
            epic_res_raw = EPIC(bulk=eset_mat)
            t(epic_res_raw$cellFractions)
         }) %>% 
    as_tibble(rownames="method_cell_type")
  
  res = cell_type_mapping %>%
    select(method_cell_type=!!method, cell_type) %>% 
    inner_join(res, by="method_cell_type") %>%
    select(-method_cell_type)
  
  res
}