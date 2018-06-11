#' @importFrom testit assert
#' @importFrom dplyr select
#' @import magrittr
NULL

#' Table mapping the cell types from methods/datasets to a single, controlled vocabulary
#' 
#' @export
cell_type_mapping = readxl::read_xlsx(system.file("extdata", "cell_type_mapping.xlsx", package="immunedeconv", mustWork=TRUE), 
                                     sheet="mapping") %>% select(method_dataset, method_cell_type, cell_type)

#' Available methods and datasets
#' 
#' @export
available_datasets = cell_type_mapping %>% pull(method_dataset) %>% unique()
  
cell_type_list =  readxl::read_excel(system.file("extdata", "cell_type_mapping.xlsx", package="immunedeconv", mustWork=TRUE),
                                    sheet = "controlled_vocabulary") %>% select(parent, cell_type)

#' Available cell types in the controlled vocabulary organized as a lineage tree. 
#' 
#' @details a `data.tree` object
#' @export
cell_type_tree = cell_type_list %>% as.data.frame() %>% data.tree::FromDataFrameNetwork()


# Access nodes by name in O(1). Node names are unique in our tree. 
assert(length(cell_type_list$cell_type) == length(unique(cell_type_list$cell_type)))
node_by_name = cell_type_tree$Get(function(node){node})


#' Use a tree-hierarchy to map cell types among different methods. 
#' 
#' `cell_type`` refers to a cell type from the controlled vocabulary (CV). 
#' `method_cell_type` refers to a cell type from a method or dataset. 
#' 
#' @param use_cell_types list of cell types from the CV to map to
#' @param fractions named vector or list of cell type fractions. Names corresponds to the method_cell_types. 
#' @param method1 method or dataset with which `fractions` was generated
#' 
#' @return numeric vector with CV cell types as names
#' 
#' @export
map_cell_types = function(use_cell_types, fractions, method) { 
  fractions = as.list(fractions)
  lapply(use_cell_types, function(cell_type) {
    find_children(node_by_name[[cell_type]], fractions1, method1)
  })
}


#' recursively traverse the cell type tree to resolve the mapping. 
#' 
#' @param node data.tree::Node corresponding to a controlled vocabulary cell type
#' @param fractions a named list of fractions for each method_cell_type
#' @param method character identifing the method in the celltype_mapping
#' 
#' @return numeric Either (1) the value of the method_cell_type mapped to cell_type, 
#'                  (2) the sum of all child nodes (recursively) of cell_type
#'                  (3) NA, if the mapping cannot be resolved, i.e. at least one of the child nodes is missing. 
find_children = function(node, fractions, method) {
  cell_type = node$name
  tmp_method_celltype = cell_type_mapping %>% filter(cell_type == !!cell_type, method_dataset == method) %>% pull(method_cell_type)
  assert("Method cell type is uniquely mapped to a cell type", length(tmp_method_celltype) <= 1)
  if(length(tmp_method_celltype) == 1) {
    assert("method_cell_type is available in the given fractions vector", tmp_method_celltype %in% names(fractions))
    fractions[[tmp_method_celltype]]
  } else {
    if(!node$isLeaf) {
      lapply(node$children, function(child) {
        find_children(child, fractions, method)
      }) %>% purrr::reduce(sum)
    } else {
      NA
    }
  }
}


#' #' annotate cell types with the cell types available in the dataset
#' #' @export
#' map_results_to_dataset = function(results, which_dataset=available_datasets) {
#'   assert("the chosen dataset is not available. ", which_dataset %in% available_datasets)
#'   celltype2dataset_mapping %>% 
#'     select(cell_type, !!which_dataset) %>% 
#'     drop_na() %>% 
#'     inner_join(results, by=c("cell_type")) %>% 
#'     select(-cell_type)
#' }
