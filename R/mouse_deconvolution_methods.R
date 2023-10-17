#' Collection of deconvolution methods for mouse data.
#'
#' @import methods
#' @import dplyr
#' @importFrom testit assert
#' @import readr
#' @importFrom matrixStats rowSds
#' @import stringr
#' @import mMCPcounter
#' @import ComICS
#' @import ImmuCellAImouse
#' @importFrom biomaRt useEnsembl getLDS
#' @importFrom tibble as_tibble
#' @importFrom rlang dots_list
#' @importFrom utils capture.output read.csv read.table tail write.table
#'
#' @name mouse_deconvolution
NULL


#' List of supported mouse deconvolution methods

#' The methods currently supported are
#' `mmcp_counter`, `seqimmucc`, `dcq`, `base`, `immucellai_mouse`
#' The object is a named vector. The names correspond to the display name of the method,
#' the values to the internal name.
#'
#' @export
deconvolution_methods_mouse <- c(
  "mMCPcounter" = "mmcp_counter",
  "seqImmuCC" = "seqimmucc",
  "DCQ" = "dcq",
  "BASE" = "base",
  "ImmuCellAI_mouse" = "immucellai_mouse"
)




############################################################################
# Deconvolution functions for consistently accessing each method
#
# These functions are called from the generic `deconvolute_mouse()` function.
# They can also be used by the end-user to access method-specific
# arguments.
############################################################################



#' Deconvolute using mMCP-counter
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples. Should be TPM normalized
#' @param log2 logical. If TRUE, log2 transforms the expression matrix
#' @param gene_id specifies the format of the gene IDs. Possible choices are 'Gene.Symbol'
#'    (default), 'ENSEMBL.ID' or 'Probes' (Affymetrix array probes)
#' @param genome specifies the mouse genome version to use, GCRm39 (default) or GCRm38
#' @export
#'
deconvolute_mmcp_counter <- function(gene_expression_matrix, log2 = TRUE,
                                     gene_id = "Gene.Symbol", genome = "GCRm39") {
  if (log2 == TRUE) {
    gene_expression_matrix <- log2(gene_expression_matrix + 1)
  }
  arguments <- dots_list(
    exp = gene_expression_matrix,
    features = gene_id,
    genomeVersion = genome
  )

  call <- rlang::call2(mMCPcounter::mMCPcounter.estimate, !!!arguments)
  results <- eval(call)
  results <- as.matrix(results)
  return(results)
}

#' Deconvolute using seqImmuCC
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples. Data
#'    should NOT be normalized (raw counts)
#' @param algorithm the method to use for deconvolution. Possible choices are
#'    'SVR' for CIBERSORT or 'LLSR' for least squares regression
#' @param ... passed through to original quantiseq method. A native argument takes precedence
#'   over an immunedeconv argument
#' @export
#'
deconvolute_seqimmucc <- function(gene_expression_matrix,
                                  algorithm = c("SVR", "LLSR"), ...) {
  signature.path <- system.file("extdata", "mouse_deconvolution", "sig_matr_seqImmuCC.txt",
    package = "immunedeconv", mustWork = TRUE
  )

  if (algorithm == "SVR") {
    assert("CIBERSORT.R is provided", exists("cibersort_binary", envir = config_env))
    source(get("cibersort_binary", envir = config_env))

    temp.expression.file <- tempfile()

    write_tsv(as_tibble(gene_expression_matrix, rownames = "gene_symbol"), path = temp.expression.file)
    perm <- 100
    # results = CIBERSORT(signature.path, temp.expression.file, perm, QN=FALSE, ...)
    # results = results %>%
    #  .[, !colnames(.) %in% c("RMSE", "P-value", "Correlation")]


    arguments <- dots_list(signature.path, temp.expression.file,
      perm = perm,
      QN = FALSE, ..., .homonyms = "last"
    )

    call <- rlang::call2(CIBERSORT, !!!arguments)
    results <- eval(call)
    results <- results %>%
      .[, !colnames(.) %in% c("RMSE", "P-value", "Correlation")]
  } else if (algorithm == "LLSR") {
    signature <- read.table(signature.path, header = T, sep = "\t", row.names = 1, check.names = F)
    results <- seqImmuCC_LLSR(signature, gene_expression_matrix, ...)
  }


  return(t(results))
}

#' Deconvolute using DCQ
#'
#'
#' Cell types obtained with this method can be further reduced with the
#'    map_result_to_celltypes function, using 'DCQ' as the method
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples
#' @param ref_samples a vector containing the names of the samples to be
#'    used to standardize the expressions. If NULL, all the samples are used.
#' @param n_repeats the number of models to generate. Predicted cell quantities
#'    will then be averaged.
#' @param combine_cells logical. This method estimates several cell types (~200).
#'    If TRUE (default), these are combined into 19 major immune-stromal cell types.
#' @param ... passed through to original quantiseq method. A native argument takes precedence
#'   over an immunedeconv argument
#' @export
#'
deconvolute_dcq <- function(gene_expression_matrix,
                            ref_samples = NULL, n_repeats = 10,
                            combine_cells = TRUE, ...) {
  # We need to scale the counts
  gene_expression_matrix <- as.matrix(gene_expression_matrix)

  if (is.null(ref_samples)) {
    rows.means <- rowMeans(gene_expression_matrix)
    rows.SD <- rowSds(gene_expression_matrix)
  } else {
    rows.means <- rowMeans(gene_expression_matrix[, ref_samples])
    rows.SD <- rowSds(gene_expression_matrix[, ref_samples])
  }
  rows.SD[rows.SD == 0] <- 1

  gene_expression_matrix <- (gene_expression_matrix - rows.means) / rows.SD

  arguments <- dots_list(
    reference_data = ComICS::immgen_dat,
    mix_data = gene_expression_matrix,
    marker_set = ComICS::DCQ_mar,
    number_of_repeats = n_repeats,
    ..., .homonyms = "last"
  )
  call <- rlang::call2(dcq, !!!arguments)
  results <- eval(call)
  results <- results$average

  if (combine_cells == TRUE) {
    results <- reduce_mouse_cell_types(results, "sum")
  }

  return(t(results))
}

#' Deconvolute using BASE
#'
#' Cell types obtained with this method can be further reduced with the
#'    map_result_to_celltypes function, using 'DCQ' as the method
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples. Data
#'    should be TPM normalized and log10 scaled.
#' @param n_permutations the number of permutations of each sample expression
#'    to generate. These are used to normalize the results.
#' @param log10 logical. if TRUE, log10 transforms the expression matrix.
#' @param combine_cells logical. This method estimates several cell types (~180).
#'    If TRUE (default), these are combined into 19 major immune-stromal cell types.
#' @export
deconvolute_base_algorithm <- function(gene_expression_matrix, n_permutations = 10,
                                       log10 = TRUE, combine_cells = TRUE) {
  base.compendium.path <- system.file("extdata", "mouse_deconvolution", "BASE_cell_compendium.rds",
    package = "immunedeconv", mustWork = TRUE
  )
  cell.compendium <- readRDS(base.compendium.path)

  results <- base_algorithm(gene_expression_matrix, cell.compendium, perm = n_permutations)

  if (combine_cells == TRUE) {
    results <- reduce_mouse_cell_types(results, "median")
  }

  return(t(results))
}


#' List of cells in the different immucellai layers (mouse method)

#' @param layer one of '1', '2', '3'
#' @param print logical, wether to print the cell types included in the layer
#' @export
#'
get_immucellai_mouse_layers <- function(layer = c('1', '2', '3'), print = FALSE) {

  layer <- match.arg(layer);

  cell_types <- switch(layer,
                       '1' = c('B_cell', 'Dendritic_cells', 'T_cell', 'Granulocytes', 'Monocytes', 'NK', 'Macrophage'),
                       '2' = c('B1_cell', 'Follicular_B', 'Germinal_center_B', 'Marginal_Zone_B',
                               'Memory_B', 'Plasma_cell', 'cDC1', 'cDC2', 'MoDC', 'pDC', 'Basophil',
                               'Eosinophil', 'mast_cell', 'Neutrophils', 'M1_macrophage', 'M2_macrophage',
                               'CD4_T_cell', 'CD8_T_cell', 'NKT', 'Tgd'),
                       '3' = c('CD4_Tm', 'Naive_CD4_T', 'T_helper_cell', 'Treg', 'CD8_Tc', 'CD8_Tcm', 'CD8_Tem',
                               'CD8_Tex', 'Naive_CD8_T')
  )

  if(print){

    data <- tibble("method_cell_type" = cell_types) %>%
      annotate_cell_type(method = 'immucellai_mouse')
    message(paste0("Cell types in layer ", layer, ": \n \n", data[['cell_type']]))
  }

  return(cell_types)
}


#' Deconvolute using ImmuCellAI-mouse
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples. Should be TPM normalized
#' @param arrays Set to TRUE if microarray data, to FALSE for RNASeq (`rnaseq` parameter in ImmuCellAI_mouse)
#' @param layer2 logical, wether to return the fractions for cell types in layer 2
#' @param layer3 logical, wether to return the fractions for cell types in layer 3 (if TRUE, this will NOT override layer2)
#' @export
#'
deconvolute_immucellai_mouse <- function(gene_expression_matrix, arrays,
                                         layer2 = FALSE, layer3 = FALSE, group = FALSE,
                                         customer = FALSE, ...) {

  data.type <- ifelse(arrays, 'microarray', 'rnaseq')
  group.tag <- ifelse(group, 1, 0)

  arguments <- dots_list(
    sample = gene_expression_matrix,
    data_type = data.type,
    group_tag = group.tag,
    customer = 0, ...,
    .homonyms = "last"
  )

  call <- rlang::call2(ImmuCellAI_mouse, !!!arguments)
  results <- eval(call)
  results <- as.data.frame(results$abundance)

  if(!layer2 && !layer3){
    cell_types <- get_immucellai_mouse_layers('1', FALSE)

  } else if(layer2){

    cell_types_1 <- get_immucellai_mouse_layers('1', FALSE)
    cell_types_2 <- get_immucellai_mouse_layers('2', FALSE)

    cell_types <- c(cell_types_1, cell_types_2)
    cell_types <- cell_types[!cell_types %in% c('T_cell', 'Dendritic_cells', 'Granulocytes', 'B_cell', 'Macrophage')]

    if(layer3){
      cell_types_3 <- get_immucellai_mouse_layers('3', FALSE)
      cell_types <- c(cell_types, cell_types_3)
      cell_types <- cell_types[!cell_types %in% c('CD4_T_cell', 'CD8_T_cell')]
    }

  }

  results <- results[, c(cell_types, 'Infiltration_score')]
  results <- t(results)

  return(results)
}



#' Perform deconvolution on a mouse RNAseq dataset
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples.
#'    Data should be TPM normalized, except when using seqImmuCC where raw
#'    counts are prefereable.
#' @param method string specifying the method
#' @param rmgenes noisy genes to be removed from the analysis
#' @param column Only relevant if `gene_expression` is an ExpressionSet. Defines in which column
#'   of fData the HGNC symbol can be found.
#' @param arrays Set to TRUE if microarray data, to FALSE for RNASeq. Default FALSE
#' @param algorithm statistcal algorithm for SeqImmuCC (ignored by all other methods)
#' @param ... additional arguments passed to the respective method
#' @return a dataset with the estimated fractions/scores, with samples in rows and cell types in column
#'
#' @export
deconvolute_mouse <- function(gene_expression_matrix,
                              method = deconvolution_methods_mouse,
                              rmgenes = NULL, column = "gene_symbol",
                              arrays = FALSE, algorithm = NULL, ...) {
  message(paste0("\n", ">>> Running ", method))


  if (!is.null(rmgenes)) {
    gene_expression_matrix <- gene_expression_matrix[!rownames(gene_expression_matrix) %in% rmgenes, ]
  }

  results <- switch(method,
    mmcp_counter = deconvolute_mmcp_counter(gene_expression_matrix, ...),
    seqimmucc = deconvolute_seqimmucc(gene_expression_matrix, algorithm, ...),
    dcq = deconvolute_dcq(gene_expression_matrix, ...),
    base = deconvolute_base_algorithm(gene_expression_matrix, ...),
    immucellai_mouse = deconvolute_immucellai_mouse(gene_expression_matrix, arrays, layer2=TRUE, layer3=TRUE, ...)
  )

  results <- results %>%
    as_tibble(., rownames = "method_cell_type") %>%
    annotate_cell_type(method = method)

  return(results)
}


#' This function converts the mouse gene symbols into corresponding human ones.
#'
#' This function relies on the `biomaRt`` package.
#'
#' @param gene_expression_matrix a m x n matrix with m genes and n samples.
#'    Gene symbols must be the rownames of the matrix.
#' @param mirror the ensembl mirror to use. Possible choices are 'www' (default),
#'    'uswest', 'useast', 'asia'
#' @return the same matrix, with the counts for the corresponding human genes.
#'    This matrix can directly be used with the immunedeconv methods. A message
#'    will display the ratio of original genes which were converted.
#'
#' @export
mouse_genes_to_human <- function(gene_expression_matrix, mirror = "www") {
  gene.names.mouse <- rownames(gene_expression_matrix)
  gene_expression_matrix$gene_name <- gene.names.mouse

  # human = useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
  # mouse = useMart('ensembl', dataset = 'mmusculus_gene_ensembl')
  human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = mirror)
  mouse <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", mirror = mirror)
  genes.retrieved <- getLDS(
    attributes = c("mgi_symbol"), filters = "mgi_symbol", values = gene.names.mouse,
    mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows = T
  )


  newGenes.counts <- gene_expression_matrix %>%
    left_join(., genes.retrieved, by = c("gene_name" = "MGI.symbol")) %>%
    select(., -c("gene_name")) %>%
    select(., c("HGNC.symbol", everything())) %>%
    .[!(is.na(.$HGNC.symbol)), ]

  colnames(newGenes.counts)[1] <- "gene_name"
  newGenes.counts <- newGenes.counts[!(duplicated(newGenes.counts$gene_name)), ] %>%
    as.data.frame(.)
  rownames(newGenes.counts) <- newGenes.counts$gene_name
  newGenes.counts <- select(newGenes.counts, -c("gene_name"))

  fraction <- 100 * (nrow(newGenes.counts) / nrow(gene_expression_matrix)) %>%
    round(., 1)

  message(paste0("ATTENTION: Only the ", fraction, "% of genes was maintained"))

  return(newGenes.counts)
}
