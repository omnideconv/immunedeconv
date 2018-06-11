#@importFrom testit assert

# celltype_mapping = readxl::read_xlsx(system.file("extdata", "cell_type_mapping2.xlsx", package="immunedeconv", mustWork=TRUE),
# sheet="mapping")
celltype_mapping = readxl::read_xlsx("/storage/home/sturm/projects/immune_deconvolution_methods/inst/extdata/cell_type_mapping2.xlsx", 
                                     sheet="mapping") %>% select(method_dataset, method_cell_type, cell_type)
celltype_list =  readxl::read_excel("/storage/home/sturm/projects/immune_deconvolution_methods/inst/extdata/cell_type_mapping2.xlsx",
                                    sheet = "controlled_vocabulary") %>% select(parent, cell_type)
celltype_tree = celltype_list %>% as.data.frame() %>% data.tree::FromDataFrameNetwork()

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
  lapply(use_cell_types, function(cell_type) {
    # build reference
    ref = find_children(node_by_name[[cell_type]], fractions1, method1)
    # build estimate
  })
}


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
