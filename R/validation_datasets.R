#' @importFrom assertthat assert_that
#' @import dplyr
NULL

#' @export
download_dataset = function(which_dataset=available_datasets) {
  
}

#' a mapping of the various cell_types between the different methods to a common standard
#' @export
celltype2dataset_mapping = readxl::read_xlsx(system.file("extdata", "cell_type_mapping.xlsx", package="immunedeconv", mustWork=TRUE),
                                             sheet="celltype2dataset")

#' List of available validation datasets
#' 
#' @export
available_datasets = celltype2dataset_mapping %>% select(-cell_type) %>% colnames()

#' annotate cell types with the cell types available in the dataset
#' @export
map_results_to_dataset = function(results, which_dataset=available_datasets) {
  assert_that(which_dataset %in% available_datasets, msg="the chosen dataset is not available. ")
  celltype2dataset_mapping %>% 
    select(cell_type, !!which_dataset) %>% 
    drop_na() %>% 
    inner_join(results, by=c("cell_type")) %>% 
    select(-cell_type)
}

#' make a random bulk sample from a single-cell dataset
#' 
#' @param eset `Biobase::ExpressionSet` with a `cell_type` column in `pData`.
#' 
#' @export
make_random_bulk = function(eset, n_cells=500) {
  # list of available cell_types
  available_cell_types = pData(eset) %>%
    as_tibble() %>%
    select(cell_type) %>%
    distinct()
  
  # create random subsample of available cells
  cell_ids = base::sample(dim(eset)["Samples"], size=n_cells)
  reduced_eset = eset[,cell_ids]
  
  # simulated bulk tissue as sum of all selected single cells
  expr = apply(exprs(reduced_eset), 1, sum) %>% as_tibble()
  
  # calculate known proportions of the sample
  proportions = pData(reduced_eset) %>% 
    select(cell_type) %>% 
    group_by(cell_type) %>%
    count() %>%
    mutate(fraction=n/!!n_cells) %>%
    select(-n) %>%
    right_join(available_cell_types, by="cell_type") %>%
    mutate(fraction=ifelse(is.na(fraction),0,fraction))
  
  return(list(expr=expr, proportions=proportions))
}


#' make a random expression set from a single-cell dataset
#' 
#' @export
make_bulk_eset = function(eset, n=100, n_cells=500) {
  bulk_samples = lapply(1:n, function(x){ make_random_bulk(eset) })
  expr = bind_cols(lapply(bulk_samples, function(x){x$expr})) %>% as.matrix()
  pdata = bind_rows(lapply(bulk_samples, function(x) {
    # suprressWarning: setting row names on a tibble is deprecated
    suppressWarnings({
      x$proportions %>% column_to_rownames(var="cell_type") %>% t() %>% as_tibble()
    })
  }))
  
  # suprressWarning: setting row names on a tibble is deprecated
  suppressWarnings({
    rownames(pdata) = colnames(expr)
  })
  rownames(expr) = fData(eset)$gene_symbol
  
  ExpressionSet(expr, 
                phenoData = new("AnnotatedDataFrame", data=pdata),
                featureData = new("AnnotatedDataFrame", fData(eset)))
}


