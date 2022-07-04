#' Functions for mapping cell types between deconvolution methods and datasets.
#'
#' makes use of the manually curated mapping in `inst/extdata/cell_type_mapping.xlsx`.
#'
#' @importFrom testit assert
#' @importFrom dplyr select
#' @importFrom stats na.omit
#' @import magrittr
#'
#' @name cell_type_mapping
NULL

#' Table mapping the cell types from methods/datasets to a single, controlled vocabulary.
#'
#' Columns: `method_dataset`, `method_cell_type`, `cell_type`.
#'
#' See `inst/extdata/cell_type_mapping.xlsx` for more details.
#'
#' @export
cell_type_map <- NULL
# gets attached on .onLoad, see zzz.R
.get_cell_type_map <- function() {
  readxl::read_xlsx(system.file("extdata", "cell_type_mapping.xlsx",
    package = "immunedeconv", mustWork = TRUE
  ),
  sheet = "mapping"
  ) %>%
    select(method_dataset, method_cell_type, cell_type) %>%
    na.omit()
}


#' Available methods and datasets.
#'
#' A list of all methods (e.g. `cibersort`) and datasets (e.g. `schelker_ovarian`) for
#' that the cell types are mapped to the controlled vocabulary.
#'
#' @export
available_datasets <- NULL
# gets attached on .onLoad, see zzz.R
.get_available_datasets <- function() {
  cell_type_map %>%
    pull(method_dataset) %>%
    unique()
}


#' List with controlled cell-type vocabulary
cell_type_list <- NULL
.get_cell_type_list <- function() {
  suppressWarnings({
    tmp_sheet <- readxl::read_excel(system.file("extdata", "cell_type_mapping.xlsx",
      package = "immunedeconv", mustWork = TRUE
    ),
    sheet = "controlled_vocabulary"
    )
  })
  tmp_list <- tmp_sheet %>%
    select(parent, cell_type, optional) %>%
    mutate(optional = optional %in% TRUE)
  assert("Node names are unique", length(tmp_list$cell_type) == length(unique(tmp_list$cell_type)))
  tmp_list
}


#' Available cell types in the controlled vocabulary organized as a lineage tree.
#'
#' @details a `data.tree` object
#' @name cell_type_tree
#' @export
cell_type_tree <- NULL
# gets attached on .onLoad, see zzz.R
.get_cell_type_tree <- function() {
  cell_type_list %>%
    as.data.frame() %>%
    data.tree::FromDataFrameNetwork()
}


#' Lookup dictionary for cell-type nodes
#'
#' Access nodes by name in O(1). Node names are unique in our tree.
#' gets attached on .onLoad, see zzz.R
node_by_name <- NULL
# gets attached on .onLoad, see zzz.R
.get_node_by_name <- function() {
  cell_type_tree$Get(function(node) {
    node
  })
}


#' Use a tree-hierarchy to map cell types among different methods.
#'
#' `cell_type` refers to a cell type from the controlled vocabulary (CV).
#' `method_cell_type` refers to a cell type from a method or dataset.
#'
#' @param use_cell_types list of cell types from the CV to map to
#' @param fractions Dataframe with cell types as rows and samples as columns. Rownames corresponds to the method_cell_types.
#'                  A named vector will be coerced into a one-column dataframe.
#' @param method_dataset method or dataset with which `fractions` was generated. If NULL, the input data already uses the CV.
#'
#' @return numeric vector with CV cell types as names
#'
#' @export
map_cell_types <- function(use_cell_types, fractions, method_dataset = NULL) {
  fractions <- as.data.frame(fractions)
  tmp_res <- lapply(use_cell_types, function(cell_type) {
    assert(paste("cell type is in node list:", cell_type), cell_type %in% names(node_by_name))
    find_children(node_by_name[[cell_type]], fractions, method_dataset)
  })
  do.call("rbind", tmp_res)
}

#' sum up the fraction of all child nodes to a single value.
#' Take into account which cell types are optional.
#'
#' @param df with samples in columns and cell types in rows. Contains a logical column `optional` which states which
#'   cell type may be missing
summarise_children <- function(df) {
  # store the 'optional' column in a vector and remove it from the dataframe.
  optional <- df$optional
  df$optional <- NULL

  # generate a vector which indicates samples which have *only* NA cell types.
  # we do not want to set them to zero, despite all cell types being optional
  all_na <- apply(is.na(df), 2, all)

  # create a mask which cell types we set to zero.
  # -> we want to set all NAs to zero that are optional, unless there are only NAs
  set_to_zero <- (optional %*% t(!all_na)) & is.na(df)
  df[set_to_zero] <- 0

  # finally, sum up the children.
  df_sum <- summarise_all(df, funs(sum))
  df_sum
}

#' Recursive helper function for map_cell_types.
#'
#' Traverses the cell type hierarchy, summing up cell type fractions when necessary.
#'
#' @param node data.tree::Node corresponding to a controlled vocabulary cell type
#' @param fractions a named list of fractions for each method_cell_type
#' @param method_dataset character identifying the method or dataset in the celltype_mapping. If NULL, the input
#'                         vector already uses controlled vocabulary.
#'
#' @return numeric Either (1) the value of the method_cell_type mapped to cell_type,
#'                  (2) the sum of all child nodes (recursively) of cell_type
#'                  (3) NA, if the mapping cannot be resolved, i.e. at least one of the child nodes is missing.
find_children <- function(node, fractions, method_dataset = NULL) {
  cell_type <- node$name
  if (is.null(method_dataset)) {
    tmp_cell_type <- rownames(fractions)[rownames(fractions) == cell_type]
  } else {
    tmp_cell_type <- cell_type_map %>%
      filter(cell_type == !!cell_type, method_dataset == !!method_dataset) %>%
      pull(cell_type)
  }
  assert("Method cell type is uniquely mapped to a cell type", length(tmp_cell_type) <= 1)
  if (length(tmp_cell_type) == 1 && tmp_cell_type %in% rownames(fractions)) {
    # assert(paste("tmp_cell_type is available in the given fractions vector:", tmp_cell_type), tmp_cell_type %in% rownames(fractions))
    fractions[tmp_cell_type, , drop = FALSE]
  } else {
    if (!node$isLeaf) {
      # recursively sum up the child nodes
      tmp_sum <- lapply(node$children, function(child) {
        tmp_row <- find_children(child, fractions, method_dataset)
        tmp_row$optional <- child$optional
        tmp_row
      }) %>%
        bind_rows() %>%
        summarise_children()
      rownames(tmp_sum) <- cell_type
      tmp_sum
    } else {
      # return NA row
      tmp_na <- as.data.frame(matrix(nrow = 1, ncol = ncol(fractions)))
      rownames(tmp_na) <- cell_type
      colnames(tmp_na) <- colnames(fractions)
      tmp_na
    }
  }
}


#' Map a result table as generated by `deconvolute` to a user-defined list of cell types.
#'
#' @seealso map_cell_types
#'
#' @param result result data.frame generated by the deconvolution method
#' @param use_cell_types list of cell_types to map to
#' @param method method or cell type to used. If method is NULL, it is expected to already use controlled vocabulary terms
#'
#' @importFrom tibble column_to_rownames
#'
#' @export
map_result_to_celltypes <- function(result, use_cell_types, method = NULL) {
  result_mat <- result %>%
    as.data.frame() %>%
    column_to_rownames("cell_type")
  map_cell_types(use_cell_types, result_mat, method)
}


#' Get all children of a certain cell type.
#'
#' If method is NULL, this function will return ALL children of the given
#' cell type in the cell type hierarchy.
#'
#' If a method is given, this function will return the cell type names
#' as used by `method` and stop as soon as a cell type maps to a cell type
#' provided by `method`.
#'
#' @param cell_type cell type name that appears in the cell type hierarchy.
#' @param method method or dataset from the cell type mapping.
#'
#' @return character vector of cell type names.
#'
#' @export
get_all_children <- function(cell_type, method = NULL) {
  if (!cell_type %in% names(node_by_name)) {
    stop(sprintf("unknown cell type: %s", cell_type))
  }
  if (is.null(method)) {
    names(node_by_name[[cell_type]]$Get("name"))
  } else {
    .get_all_children(node_by_name[[cell_type]], method)
  }
}

#' Recursive helper function for get_all_children.
#'
#' @param node Node in the cell type hierarchy. Will look for children of this node.
#' @param method method or dataset from the cell type mapping.
.get_all_children <- function(node, method) {
  tmp_method_cell_type <- cell_type_map %>%
    filter(method_dataset == method, cell_type == node$name) %>%
    pull(method_cell_type)
  assert("Method cell type is uniquely mapped to a cell type", length(tmp_method_cell_type) <= 1)
  if (length(tmp_method_cell_type) == 1) {
    # stop at the highest level the cell type is mapped to a method
    return(tmp_method_cell_type)
  } else if (!node$isLeaf) {
    # otherwise, continue searching children
    lapply(node$children, function(child) {
      .get_all_children(child, method)
    }) %>%
      unlist() %>%
      unname()
  }
}
