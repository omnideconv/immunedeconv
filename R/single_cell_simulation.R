#' Functions for generating simulated bulk samples from a single cell dataset.
#'
#' @importFrom testit assert
#' @import dplyr
#' @importFrom tibble rowid_to_column
#'
#' @name single_cell_simulation
NULL


#' Scale sample to TPM
#'
#' a TPM sample always sums to 1,000,000.
#'
#' This function takes an input vector and
#' rescales it, so that the sum of the vector equals 1e6.
#'
#' @param sample numeric vector of gene expression values
#'
#' @export
scale_to_million <- function(sample) {
  (sample / sum(sample)) * 1e6
}

#' Make a random bulk sample from a single-cell dataset
#'
#' Take random single cells from the input expression set, given
#' the cell_fractions vector and combine them by taking the mean.
#'
#' @param eset `Biobase::ExpressionSet` with a `cell_type` column in `pData`.
#'    Each sample is the gene expression of a single cell. The cell
#'    type of each cell needs to be denoted in the pData table in the
#'    `cell_type` column.
#' @param cell_fractions named list indicating the fraction of each cell type
#'    which will be in the sample. The names of the list need to correspond
#'    to the `cell_type` column in the ExpressionSet.
#' @param n_cells number of single cells to integrate into a sample
#' @param combine callback function used to aggregate the counts.
#'
#' @return a single-column tibble with simulated expression for each gene.
#'   The column sum is scaled to 1 million (TPM)
#'
#' @examples
#' suppressPackageStartupMessages(library(Biobase))
#' suppressPackageStartupMessages(library(tibble))
#'
#' # generate toy matrix with three marker genes and three cell types
#' expr <- matrix(c(
#'   rep(c(1, 0, 0), 300),
#'   rep(c(0, 1, 0), 300),
#'   rep(c(0, 0, 1), 300)
#' ), nrow = 3)
#'
#' # generate a featureData and phenoData data-frame.
#' # row names must be consistent between expr and featureData.
#' gene_names <- c("CD8A", "CD4", "CD19")
#' rownames(expr) <- gene_names
#' cell_types <- c(rep("T cell CD8+", 300), rep("T cell CD4+", 300), rep("B cell", 300))
#' pdata <- data.frame(cell_type = cell_types)
#' fdata <- data.frame(gene_symbol = gene_names)
#' rownames(fdata) <- gene_names
#'
#' # tie expr, fdata and pdata together in expression set
#' eset <- ExpressionSet(expr,
#'   phenoData = as(pdata, "AnnotatedDataFrame"),
#'   featureData = as(fdata, "AnnotatedDataFrame")
#' )
#'
#' # make a random bulk sample.
#' make_random_bulk(eset,
#'   c("T cell CD8+" = 0.3, "B cell" = 0.4, "T cell CD4+" = 0.3),
#'   n_cells = 1000
#' )
#'
#' @export
make_random_bulk <- function(eset, cell_fractions, n_cells = 500, combine = mean) {
  cell_numbers <- round(cell_fractions * n_cells)

  # sample n cells from each cell type as specified in `cell_numbers`
  cell_ids <- lapply(names(cell_numbers), function(cell_type) {
    n <- cell_numbers[[cell_type]]
    pData(eset) %>%
      as_tibble() %>%
      rowid_to_column("id") %>%
      filter(cell_type == !!cell_type) %>%
      dplyr::sample_n(n, replace = TRUE)
  }) %>%
    bind_rows() %>%
    pull(id)

  # create random subsample of available cells
  reduced_eset <- eset[, cell_ids]

  # simulated bulk tissue as mean of all selected single cells
  expr <- apply(exprs(reduced_eset), 1, combine) %>%
    scale_to_million() %>%
    as_tibble()

  return(expr)
}


#' Make a random expression set from a single-cell dataset
#'
#' @param eset `Biobase::ExpressionSet` with a `cell_type` column in `pData`.
#' @param cell_fractions n x n_cell_types dataframe with the fraction for each
#'     sample in each row.
#' @param n_cells number of single cells to use in each sample
#' @param combine callback function used to aggregate the counts.
#'
#' @seealso make_random_bulk
#'
#' @return matrix with genes in rows and samples in columns. The column
#'   sum is scaled to 1 million (TPM).
#'
#' @examples
#' suppressPackageStartupMessages(library(Biobase))
#' suppressPackageStartupMessages(library(tibble))
#'
#' # generate toy matrix with three marker genes and three cell types
#' expr <- matrix(c(
#'   rep(c(1, 0, 0), 300),
#'   rep(c(0, 1, 0), 300),
#'   rep(c(0, 0, 1), 300)
#' ), nrow = 3)
#'
#' # generate a featureData and phenoData data-frame.
#' # row names must be consistent between expr and featureData.
#' gene_names <- c("CD8A", "CD4", "CD19")
#' rownames(expr) <- gene_names
#' cell_types <- c(rep("T cell CD8+", 300), rep("T cell CD4+", 300), rep("B cell", 300))
#' pdata <- data.frame(cell_type = cell_types)
#' fdata <- data.frame(gene_symbol = gene_names)
#' rownames(fdata) <- gene_names
#'
#' # tie expr, fdata and pdata together in expression set
#' eset <- ExpressionSet(expr,
#'   phenoData = as(pdata, "AnnotatedDataFrame"),
#'   featureData = as(fdata, "AnnotatedDataFrame")
#' )
#'
#' # create data frame with cell fractions (one row per sample)
#' desired_cell_fractions <- tibble(
#'   "T cell CD8+" = c(0.1, 0.2, 0.3),
#'   "T cell CD4+" = c(0.9, 0.7, 0.5),
#'   "B cell" = c(0, 0.1, 0.2)
#' )
#' new_eset <- make_bulk_eset(eset, desired_cell_fractions, n_cells = 500)
#' exprs(new_eset)
#'
#' @export
make_bulk_eset <- function(eset, cell_fractions, n_cells = 500, combine = mean) {
  expr <- lapply(1:nrow(cell_fractions), function(i) {
    make_random_bulk(eset, cell_fractions[i, ], n_cells = n_cells, combine = combine)
  }) %>%
    bind_cols() %>%
    as.matrix()

  rownames(expr) <- fData(eset)$gene_symbol

  pdata <- cell_fractions %>% as.data.frame()
  rownames(pdata) <- colnames(expr)

  ExpressionSet(expr,
    phenoData = new("AnnotatedDataFrame", data = pdata),
    featureData = new("AnnotatedDataFrame", fData(eset))
  )
}
