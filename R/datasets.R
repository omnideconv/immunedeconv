#' Example RNA-seq dataset from the EPIC publication.
#
#' A dataset of four RNA-seq samples from patients with metastatic melanoma.
#'
#' @format a gene expression matrix with gene_symbols as rownames and sample identifiers as colnames.
#'
#' @source Racle et al. (2017), eLIFE, @url{https://doi.org/10.7554/eLife.26476.029}
"dataset_racle$expr_mat"


#' Immune cell proportions from the Racle-dataset measured with FACS.
#
#' Immune cell ccontents of four samples from patients with metastatic melanoma profiled
#' with Fluorescence activated cell sorting (FACS).
#'
#' @format a data.frame with three columns: sample, cell_type, true_fraction
#'
#' @source Racle et al. (2017), eLIFE, @url{https://doi.org/10.7554/eLife.26476.029}
"dataset_racle$ref"
