#' @importFrom testit assert
NULL


celltype_mapping = readxl::read_xlsx(system.file("extdata", "cell_type_mapping.xlsx", package="immunedeconv", mustWork=TRUE), 
                                     sheet="mapping") %>% select(method_dataset, method_cell_type, cell_type)
celltype_list =  readxl::read_excel(system.file("extdata", "cell_type_mapping.xlsx", package="immunedeconv", mustWork=TRUE),
                                    sheet = "controlled_vocabulary") %>% select(parent, cell_type)
celltype_tree = celltype_list %>% as.data.frame() %>% data.tree::FromDataFrameNetwork()

#' List of available validation datasets
#' 
#' @export
available_datasets = colnames(dplyr::select(celltype2dataset_mapping, -cell_type))


# Access nodes by name in O(1). Node names are unique in our tree. 
assert(length(celltype_list$cell_type) == length(unique(celltype_list$cell_type)))
node_by_name = celltype_tree$Get(function(node){node})


#' Use a tree-hierarchy to map cell types among different methods. 
#' 
#' @param use_cell_types list of cell types from the controlled vocabulary to map
#' @param fractions1 named vector or list of cell type fractions. Name corresponds to the cell types
#' @param fractions2 
#' @param method1 method or dataset to use for mapping `fractions1`
#' @param method2
#' 
#' @return named vector. Name of the vector corresponds to the controlled vocabulary cell types. 
map_cell_types = function(use_cell_types, fractions1, fractions2, method1, method2) { 
  fractions1 = as.list(fractions1)
  fractions2 = as.list(fractions2)
  # build reference
  ref = lapply(use_cell_types, function(cell_type) {
    find_children(node_by_name[[cell_type]], fractions1, method1)
  })
  # build estimate
  estimate = lapply(use_cell_types, function(cell_type) {
    find_children(node_by_name[[cell_type]], fractions2, method2)
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
  tmp_method_celltype = celltype_mapping %>% filter(cell_type == !!cell_type, method_dataset == method) %>% pull(method_cell_type)
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


#' annotate cell types with the cell types available in the dataset
#' @export
map_results_to_dataset = function(results, which_dataset=available_datasets) {
  assert("the chosen dataset is not available. ", which_dataset %in% available_datasets)
  celltype2dataset_mapping %>% 
    select(cell_type, !!which_dataset) %>% 
    drop_na() %>% 
    inner_join(results, by=c("cell_type")) %>% 
    select(-cell_type)
}
