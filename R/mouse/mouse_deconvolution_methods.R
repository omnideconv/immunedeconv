#' Collection of deconvolution methods for mouse data.
require(tidyverse)
require(magrittr)
require(rlang)
require(ComICS)
require(mMCPcounter)
require(biomaRt)

#' List of supported mouse deconvolution methods
#'
#' The methods currently supported are
#' `mMcp_counter`, `seqimmucc`, `dcq`, `base`
#'
#' The object is a named vector. The names correspond to the display name of the method,
#' the values to the internal name.
#'
#' @export
deconvolution_methods_mouse = c("mMCPcounter"="mMcp_counter",
                                "seqImmuCC"="seqimmucc", 
                                "DCQ"="dcq", 
                                "BASE"="base")




############################################################################
# Deconvolution functions for consistently accessing each method
#
# These functions are called from the generic `deconvolute_mouse()` function.
# They can also be used by the end-user to access method-specific
# arguments.
############################################################################



#' Deconvolute using mMCP-counter
#' 
#' @param gene.expression.matrix a m x n matrix with m genes and n samples. Should be TPM normalized
#' @param log2 logical. If TRUE, log2 transforms the expression matrix
#' @export
#' 
deconvolute_mMCP = function(gene.expression.matrix, log2 = TRUE){
  if(log2 == TRUE){gene.expression.matrix = log2(gene.expression.matrix + 1)}
  call = rlang::call2(mMCPcounter::mMCPcounter.estimate)
  results = eval(call)
  return(t(results))
}


#' Deconvolute using seqImmuCC
#' 
#' @param gene.expression.matrix a m x n matrix with m genes and n samples. Data 
#'    should NOT be normalized (raw counts)
#' @param method the method to use for deconvolution. Possible choices are
#'    'SVR' for CIBERSORT or 'LLSR' for least squares regression
#' @export
#'     
deconvolute_seqImmuCC = function(gene.expression.matrix, method = c('SVR', 'LLSR')){
  
  signature.path <- system.file('extdata', 'mouse_deconvolution', 'sig_matr',
                                   package = 'immunedeconv', mustWork=TRUE)
  source('./seqImmuCC_LLSR.R')
  
  if (method == 'SVR'){
    
    assert("CIBERSORT.R is provided", exists("cibersort_binary", envir=config_env))
    source(get("cibersort_binary", envir=config_env))
    
    temp.expression.file = tempfile()
    
    write_tsv(as_tibble(gene.expression.matrix, rownames = 'gene_symbol'), path = temp.expression.file)
    perm = 100
    results = CIBERSORT(signature.path, temp.expression.file, perm, QN=FALSE)
    results = results %>%
     .[, !colnames(.) %in% c("RMSE", "P-value", "Correlation")]
    
  } else if (method == 'LLSR'){
    
    signature = read.table(signature.path, header=T, sep="\t", row.names=1, check.names=F)
    results = LLSR(signature, gene.expression.matrix)
  }
  
  
  return(results)
}

#' Deconvolute using DCQ
#' 
#' @param gene.expression.matrix a m x n matrix with m genes and n samples
#' @param ref.samples a vector containing the names of the samples to be 
#'    used to standardize the expressions. If NULL, all the samples are used. 
#' @param n.repeats the number of models to generate. Predicted cell quantities
#'    will then be averaged.    
#' @export
#'      
deconvolute_DCQ = function(gene.expression.matrix, 
                           ref.samples = NULL, n.repeats = 10){
  
  # We need to scale the counts
  gene.expression.matrix = as.matrix(gene.expression.matrix)
  
  if (is.null(ref.samples)){
    rows.means = rowMeans(gene.expression.matrix)
    rows.SD = rowSds(gene.expression.matrix)
  } else {
    rows.means = rowMeans(gene.expression.matrix[, samples])
    rows.SD = rowSds(gene.expression.matrix[, samples])
  }
  rows.SD[rows.SD==0] = 1
  
  gene.expression.matrix = (gene.expression.matrix - rows.means)/rows.SD
  
  arguments = dots_list(reference_data = immgen_dat, 
                         mix_data = gene.expression.matrix, 
                         marker_set = DCQ_mar, 
                         number_of_repeats = n.repeats)
  call = rlang::call2(dcq, !!!arguments)
  results = eval(call)
  return(results$average)
}

#' Deconvolute using BASE
#' 
#' @param gene.expression.matrix a m x n matrix with m genes and n samples. Data
#'    should be TPM normalized and log10 scaled.
#' @param n.permutations the number of permutations of each sample expression
#'    to generate. These are used to normalize the results. 
#' @param log10 logical. if TRUE, log10 transforms the expression matrix. 
#' 
#' @export
#' 
deconvolute_BASE = function(gene.expression.matrix, n.permutations = 100, log10 = TRUE){
  
  base.compendium.path <- system.file('extdata', 'mouse_deconvolution', 'BASE_immune_compendium.rds',
                                      package = 'immunedeconv', mustWork=TRUE)
  immune.compendium = load(base.compendium.path)
  
  source('./BASE.R')
  
  results = base_full_pipeline(gene.expression.matrix, immune.compendium, perm = n.permutations, median.norm = T)
  
  return(results)
}  



#' Perform deconvolution on a mouse RNAseq dataset
#' 
#' @param gene.expression.matrix a m x n matrix with m genes and n samples. 
#'    Data should be TP normalized, except when using seqImmuCC where raw 
#'    counts are prefereable.
#' @param method string specifying the method
#' @param rmgenes noisy genes to be removed from the analysis
#' @param ... additional arguments passed to the respective method
#' @return a dataset with the estimated fractions/scores, with samples in rows and cell types in column
#' 
#' @export

deconvolute_mouse = function(gene.expression.matrix, 
                              method = deconvolution_methods_mouse, 
                              rmgenes = NULL, ...){
  message(paste0("\n", ">>> Running ", method))
  
  if(!is.null(rmgenes)) {
    gene.expression.matrix = gene.expression.matrix[!rownames(gene.expression.matrix) %in% rmgenes, ]
  }
  
  results = switch(method, 
                   mMcp_counter = deconvolute_mMCP(gene.expression.matrix, ...),
                   seqimmucc = deconvolute_seqImmuCC(gene.expression.matrix, ...),
                   dcq = deconvolute_DCQ(gene.expression.matrix, ...),
                   base = deconvolute_BASE(gene.expression.matrix, ...))
  
  # Whenever BASE/DCQ are used, the cell types are reduced
  if(method == 'dcq'){
    results = reduce_cell_types(results, 'dcq', 'sum')
  } else if(method == 'base'){
    results = reduce_cell_types(results, 'base', 'median')
  }
  
  results = t(results) %>%
    as_tibble(., rownames = 'method_cell_type')
  
  return(results)
}


#' This function converts the mouse gene symbols into corresponding human ones.
#' 
#' @param gene.expression.matrix a m x n matrix with m genes and n samples. 
#'    Gene symbols must be the rownames of the matrix. 
#' @return the same matrix, with the counts for the corresponding human genes    
#'    
#' @export
mouse_genes_to_human = function(gene.expression.matrix){
  
  gene.names.mouse = rownames(gene.expression.matrix)
  gene.expression.matrix$gene_name = gene.names.mouse
  
  library(biomaRt)
  human = useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
  mouse = useMart('ensembl', dataset = 'mmusculus_gene_ensembl')
  genes.retrieved = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = gene.names.mouse, 
                           mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  
  
  newGenes.counts <- gene.expression.matrix %>%
    left_join(., genes.retrieved, by = c('gene_name' = 'MGI.symbol')) %>%
    select(., -c('gene_name')) %>%
    select(., c('HGNC.symbol', everything())) %>%
    .[!(is.na(.$HGNC.symbol)), ]
  
  colnames(newGenes.counts)[1] <- 'gene_name'
  newGenes.counts <- newGenes.counts[!(duplicated(newGenes.counts$gene_name)), ] %>%
    as.data.frame(.)
  rownames(newGenes.counts) <- newGenes.counts$gene_name
  newGenes.counts <- select(newGenes.counts, -c('gene_name'))
  
  fraction = 100*(nrow(newGenes.counts)/nrow(gene.expression.matrix)) %>%
    round(., 1)
  
  message(paste0('ATTENTION: Only the ', fraction, '% of genes was maintained'))
  
  return(newGenes.counts)
  
}