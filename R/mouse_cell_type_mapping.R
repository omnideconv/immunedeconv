#' Functions for mapping cell types between deconvolution methods.
#'
#' @importFrom testit assert
#' @importFrom dplyr select
#' @import ComICS
#' @importFrom matrixStats rowMedians
#' @import magrittr
#' @import stringr
#'
#' @name mouse_cell_type_mapping
NULL

#' Since DCQ and BASE provide estimates for several cell types, this function
#'    combines the results to align them with the rest of the methods
#'
#' @param deconvolution.results a dataframe containing the results of the
#'    deconvolution experiment (DCQ or BASE). Cell types are in columns, samples in rows.
#' @param method the method used to combine the estimates: sum, for regression-based
#'    methods as DCQ, or median for the scoring methods such as BASE.
#' @return a data frame with the aggregated cell types
#' @export
reduce_mouse_cell_types <- function(deconvolution.results,
                                    method = c("sum", "median")) {
  table.mapping.path <- system.file("extdata", "mouse_deconvolution", "cell_correspondance.txt",
    package = "immunedeconv", mustWork = TRUE
  )
  cell.types.table <- read.table(table.mapping.path, sep = "\t", header = TRUE)
  cell.types <- unique(cell.types.table$Subtype)

  if (method == "median") {
    deconvolution.results <- scale(deconvolution.results)
  }
  results.transformed <- as.data.frame(deconvolution.results)

  for (cell in cell.types) {
    current.columns <- colnames(results.transformed) %in% cell.types.table[which(cell.types.table$Subtype == cell), 1]

    # Some cell types might not be present for both DCQ and BASE (ex. Treg)
    if (sum(current.columns) > 0) {
      if (method == "median") {
        results.transformed[cell] <- results.transformed[current.columns] %>%
          as.matrix(.) %>%
          rowMedians(.)
      } else {
        results.transformed[cell] <- results.transformed[current.columns] %>%
          rowSums(.)
      }
    }
  }

  results.transformed <- results.transformed[colnames(results.transformed) %in% cell.types]

  return(results.transformed)
}




#' This function returns the list of all cell types in BASE/DCQ results,
#'    along with the cell type they are mapped to
#' @param method the method for which we want to get the cell types
#' @return a dataframe with two columns: the exact cell type and the
#'    correspondent, more general one
#' @export
get_dcq_cell_types_list <- function(method = c("dcq", "base")) {
  table.mapping.path <- system.file("extdata", "mouse_deconvolution", "cell_correspondance.txt",
    package = "immunedeconv", mustWork = TRUE
  )
  cell.types.table <- read.table(table.mapping.path, sep = "\t", header = TRUE)

  if (method == "dcq") {
    all.cell.types <- immgen_dat %>%
      colnames(.)
  } else if (method == "base") {
    base.compendium.path <- system.file("extdata", "mouse_deconvolution", "BASE_cell_compendium.rds",
      package = "immunedeconv", mustWork = TRUE
    )

    BASE.compendium <- readRDS(base.compendium.path)
    all.cell.types <- colnames(BASE.compendium) %>%
      str_split_fixed(., "\\_", 2) %>%
      .[, 1] %>%
      unique(.)
  }

  method.output <- data.frame(
    "Name" = all.cell.types
  )

  method.output <- left_join(method.output, cell.types.table) %>%
    select(., c("Name", "Subtype"))

  method.output
}
