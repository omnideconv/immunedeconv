available_datasets = c("single_cell_schelker")

#' @export
download_dataset = function(which_dataset=available_datasets) {
  
}

#' @export
map_results_to_dataset = function(results, which_dataset=available_datasets) {
  
}

#' make a random bulk sample from a single cell dataset
#' 
#' @export
make_random_bulk = function(eset, n_cells=500) {
  cell_ids = base::sample(dim(eset)["Samples"], size=n_cells)
  reduced_eset = eset[,cell_ids]
  # simulated bulk tissue as sum of all selected single cells
  expr = apply(exprs(reduced_eset), 1, sum) %>% as_tibble()
  proportions = pData(reduced_eset) %>% 
    select(cell_type) %>% 
    group_by(cell_type) %>%
    count() %>%
    mutate(fraction=n/!!n_cells) %>%
    select(-n) %>%
    right_join(cell_names[,"cell_type"], by="cell_type") %>%
    mutate(fraction=ifelse(is.na(fraction),0,fraction))
  return(list(expr=expr, proportions=proportions))
}


#' make a random expression set from a single cell dataset
#' 
#' @export
make_bulk_eset = function(eset, n=100, n_cells=500) {
  bulk_samples = lapply(1:n, function(x){ make_random_bulk(eset) })
  expr = bind_cols(lapply(bulk_samples, function(x){x$expr})) %>% as.matrix()
  pdata = bind_rows(lapply(bulk_samples, function(x) {
    # setting row names on a tibble is deprecated
    suppressWarnings({
      x$proportions %>% column_to_rownames(var="cell_type") %>% t() %>% as_tibble()
    })
  }))
  # setting row names on a tibble is deprecated
  suppressWarnings({
    rownames(pdata) = colnames(expr)
  })
  rownames(expr) = fData(eset)$gene_symbol
  ExpressionSet(expr, 
                phenoData = new("AnnotatedDataFrame", data=pdata),
                featureData = new("AnnotatedDataFrame", fData(eset)))
}


