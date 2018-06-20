library(sva)
library(crayon)
library(sqldf)

TimerINFO <- function(string) {
  cat(green(sprintf('## %s\n', string)))
}

TimerINFO('Loading Timer Utilities')

immuneCuratedData <- paste(baseDir, '/data/precalculated/immune.expression.curated.RData', sep='')

ConvertImmuneProbeToRefgene <- function(curated.ref){
  ##----- function to preprocess the reference dataset, not necessary if the processed data "curated.ref.genes.Rdata" is available -----##

  tmpDD <- data.frame(curated.ref)
  tmpDD <- tmpDD[order(rownames(tmpDD)), ]
  ## sort the immune expression data by rownames

  colnames(tmpDD) <- gsub('\\.', '_', colnames(tmpDD))
  genes <- sapply(strsplit(rownames(tmpDD), ';'), function(x) x[[1]])
  ## remove the probe ID, only keep gene names

  tmpDD <- cbind(genes, tmpDD)
  tmpDD <- tmpDD[order(genes), ]
  ## sort by gene names

  tmp0 <- c()
  cnt <- 0

  TimerINFO('Aggregating immune expression data')

  for(i in colnames(tmpDD)[2:ncol(tmpDD)]){
    ## start from the second column (the first column is gene information)
    cat(sprintf("(%d of %d) %s\n", cnt, ncol(tmpDD) - 1 ,i))
    tmp <- sqldf(paste('select max(', i, ') from tmpDD group by genes', sep=''))
    ## select the maximum probe expression level when a gene has multiple probes

    if(length(tmp0) == 0) tmp0 <- tmp else tmp0 <- cbind(tmp0, tmp)
    cnt <- cnt + 1
  }
  colnames(tmp0) <- colnames(tmpDD)[2:ncol(tmpDD)]
  rownames(tmp0) <- unique(tmpDD[, 1])
  curated.ref.genes <- tmp0
  return(curated.ref.genes)
}


LoadImmuneGeneExpression <- function() {
  ## Load gene expression data for immune cells

  ## Returns:
  ##   A data frame of expression data for immune cells
  ##   (cols for immune cell sample, rows for "gene name;probe ID")
  if (file.exists(immuneCuratedData)) {
    ## See below: list(genes=curated.ref.genes, celltypes=curated.cell.types)
    curated.data <- get(load(immuneCuratedData))
    return(curated.data)
  }

  exp <- get(load(paste(baseDir,'/data/immune_datasets/HPCTimmune.Rdata',sep='')))

  ##----- Select single reference samples of pre-selected immune cell types -----##
  B_cell <- 362:385
  T_cell.CD4 <- grep('T_cell.CD4',colnames(exp))
  T_cell.CD8 <- grep('T_cell.CD8',colnames(exp))
  NK <- 328:331
  Neutrophil <- 344:361
  Macrophage <- 66:80
  DC <- 151:238

  curated.ref <- exp[, c(B_cell, T_cell.CD4, T_cell.CD8,
                         NK, Neutrophil, Macrophage, DC)]

  curated.cell.types <- colnames(curated.ref)
  names(curated.cell.types) <- c(rep('B_cell', length(B_cell)),
                                 rep('T_cell.CD4', length(T_cell.CD4)),
                                 rep('T_cell.CD8', length(T_cell.CD8)),
                                 rep('NK', length(NK)),
                                 rep('Neutrophil', length(Neutrophil)),
                                 rep('Macrophage', length(Macrophage)),
                                 rep('DC', length(DC)))
  curated.ref.genes <- ConvertImmuneProbeToRefgene(curated.ref)
  ret <- list(genes=curated.ref.genes, celltypes=curated.cell.types)
  save(ret,file=immuneCuratedData)
  return(ret)
}

RemoveBatchEffect <- function(cancer.exp, immune.exp, immune.cellType) {
  ## intersect the gene names of cancer.exp and immune.exp
  tmp.dd <- as.matrix(cancer.exp)
  tmp <- sapply(strsplit(rownames(cancer.exp), '\\|'),
                function(x) x[[1]])
  rownames(tmp.dd) <- tmp
  tmp.dd <- as.matrix(tmp.dd[which(nchar(tmp)>1), ])
  tmp.ss <- intersect(rownames(tmp.dd), rownames(immune.exp))

  ## bind cancer and immune expression data into one dataframe
  N1 <- ncol(tmp.dd)

  tmp.dd <- cbind(tmp.dd[tmp.ss, ], immune.exp[tmp.ss, ])
  tmp.dd <- as.matrix(tmp.dd)
  mode(tmp.dd) <- 'numeric'

  ## remove batch effects
  N2 <- ncol(immune.exp)
  tmp.batch <- c(rep(1, N1), rep(2, N2))
  tmp.dd0 <- ComBat(tmp.dd, tmp.batch, c())

  ## separate cancer and immune expression data after batch effect removing
  dd.br <- tmp.dd0[, 1:N1]
  immune.exp.br <- tmp.dd0[, (N1+1):(N1+N2)]

  ## a immune category has multiple samples, use the median expression level for a gene
  tmp0 <- c()
  for(kk in unique(names(immune.cellType))){
    tmp.vv <- which(names(immune.cellType)==kk)
    tmp0 <- cbind(tmp0, apply(immune.exp.br[, tmp.vv], 1, median, na.rm=T))
  }


  immune.exp.agg.br <- tmp0
  colnames(immune.exp.agg.br) <- unique(names(immune.cellType))
  return(list(as.matrix(dd.br), immune.exp.br, immune.exp.agg.br))
}

