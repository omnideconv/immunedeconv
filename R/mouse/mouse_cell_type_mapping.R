#' Functions for mapping cell types between deconvolution methods.

require(tidyverse)
require(magrittr)
require(rlang)
require(ComICS)


#' Since DCQ and BASE provide estimates for several cell types, this function 
#'    combines the results to align them with the rest of the methods
#'
#' @param deconvolution.results a dataframe containing the results of the 
#'    deconvolution experiment (DCQ or BASE). Cell types are in columns, samples in rows.
#' @param annotation which annotation has to be used. Possible values are 'dcq', 'base' or 
#'    'user' for user provided annotations
#' @param cell.table a dataframe with two columns: prefixes of cell names and 
#'    the corresponding cell type. Users can submit their own as well. 
#' @param method the method used to combine the estimates: sum, for regression-based 
#'    methods as DCQ, or median for the scoring methods such as BASE.   
#' @return a data frame with the aggregated cell types     
#' @export

reduce_mouse_cell_types <- function(deconvolution.results, 
                                    annotation = c('dcq', 'base', 'user'),
                                    cell.table = NULL, 
                                    method = c('sum', 'median')){
  
  dcq.mapping.path <- system.file('extdata', 'mouse_deconvolution', 'cell_types_table_DCQ.rds',
                                  package = 'immunedeconv', mustWork=TRUE)
  base.mapping.path <- system.file('extdata', 'mouse_deconvolution', 'cell_types_table_BASE.rds',
                                   package = 'immunedeconv', mustWork=TRUE)
  
  
  cell.types.table = switch(annotation,
      dcq = load(dcq.mapping.path), 
      base = load(base.mapping.path),
      user = cell.table)
  
  if (method == 'median'){deconvolution.results = scale(deconvolution.results)}
  
  results.transformed = deconvolution.results
  
  j = 1
  while (j <= nrow(cell.types.table)){
    substr <- cell.types.table[j, 1]
    cell <- cell.types.table[j, 2]
    
    str <- paste0('^', substr)
    cols <- colnames(deconv.results) %>%
      .[str_detect(., str)]
    
    if (j < nrow(cell.types.table)){
      next.cell.type <- cell.types.table[j+1, 2]
      while (next.cell.type == cell && j < nrow(cell.types.table)){
        j <- j+1
        substr <- cell.types.table[j, 1]
        str <- paste0('^', substr)
        
        cols.new <- colnames(deconv.results) %>%
          .[str_detect(., str)]
        
        cols <- c(cols, cols.new)
        next.cell.type <- cell.types.table[j+1, 2] 
      }
    }
    
    if (method == 'median'){
      results.transformed[cell] <- select(results.transformed, cols) %>%
        as.matrix(.) %>%
        rowMedians(.)
    } else {
      results.transformed[cell] <- select(results.transformed, cols) %>%
        as.matrix(.) %>%
        rowSums(.)
    }
    
    deconv.results <- select(deconv.results, -cols)
    j <- j+1
  }
  
  results.transformed <- select(results.transformed, unique(cell.types.table$cell_type))
  results.transformed
}




#' This function returns the list of all cell types in BASE/DCQ results, 
#'    along with the cell type they are mapped to
#' @param method the method for which we want to get the cell types
#' @return a dataframe with two columns: the exact cell type and the 
#'    correspondent, more general one  
#' @export    

get_cell_types_list <- function(method = c('dcq', 'base')){
  
  if (method == 'dcq'){
    
    dcq.mapping.path <- system.file('extdata', 'mouse_deconvolution', 'cell_types_table_DCQ.rds',
                package = 'immunedeconv', mustWork=TRUE)
    mapping <- load(dcq.mapping.path)
    
    all.cell.types <- ComICS::immgen_dat %>%
      colnames(.)
  } else if (method == 'base'){
    
    base.mapping.path <- system.file('extdata', 'mouse_deconvolution', 'cell_types_table_BASE.rds',
                                    package = 'immunedeconv', mustWork=TRUE)
    base.compendium.path <- system.file('extdata', 'mouse_deconvolution', 'BASE_immune_compendium.rds',
                                        package = 'immunedeconv', mustWork=TRUE)
    
    mapping <- load(base.mapping.path)
    BASE.compendium <- load(base.compendium.path)
    all.cell.types <- colnames(BASE.compendium) %>%
      str_split_fixed(., '\\_', 2) %>%
      .[, 1] %>%
      unique(.)
  }
  
  cell.types.df <- data.frame(
    cell_types = all.cell.types, 
    correspondance = rep('', length(all.cell.types))
  )
  
  for (j in 1:nrow(mapping)){
    cur.prefix <- paste0('^', mapping[j, 1])
    cell.types.df[which(str_detect(cell.types.df$cell_types, cur.prefix)), 2] <- mapping[j, 2]
  }
  
  cell.types.df
}