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
available_datasets = colnames(dplyr::select(celltype2dataset_mapping, -cell_type))


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
#' @param cell_fractions named list indicating the fraction of each cell type 
#'    which will be in the sample. 
#' 
#' @export
make_random_bulk = function(eset, cell_fractions, n_cells=500) {
  cell_numbers = round(cell_fractions * n_cells)
  
  # sample n cells from each cell type as specified in `cell_numbers`
  cell_ids = lapply(names(cell_numbers), function(cell_type) {
    n = cell_numbers[[cell_type]]
    pData(eset) %>%
      as_tibble() %>%
      rowid_to_column("id") %>%
      filter(cell_type==!!cell_type) %>% 
      dplyr::sample_n(n, replace=TRUE)
  }) %>% 
    bind_rows() %>%
    pull(id)
  
  # create random subsample of available cells
  reduced_eset = eset[,cell_ids]
  
  # simulated bulk tissue as sum of all selected single cells
  expr = apply(exprs(reduced_eset), 1, sum) %>% as_tibble()

  return(expr)
}


#' make a random expression set from a single-cell dataset
#' 
#' @param eset `Biobase::ExpressionSet` with a `cell_type` column in `pData`.
#' @param cell_fractions n x n_cell_types dataframe with the fraction for each 
#'     sample in each row. 
#' @param n_cell number of single cells to use in each sample
#' 
#' @export
make_bulk_eset = function(eset, cell_fractions, n_cells=500) {
  expr = lapply(1:nrow(cell_fractions), function(i) { 
    make_random_bulk(eset, cell_fractions[i,], n_cells=n_cells) 
  }) %>% 
    bind_cols() %>% 
    as.matrix()
  
  rownames(expr) = fData(eset)$gene_symbol
  
  pdata = cell_fractions %>% as.data.frame()
  rownames(pdata) = colnames(expr)
  
  ExpressionSet(expr, 
                phenoData = new("AnnotatedDataFrame", data=pdata),
                featureData = new("AnnotatedDataFrame", fData(eset)))
}


