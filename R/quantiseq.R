#' Use quanTIseq to deconvolute a gene expression matrix.
#'
#' Source code from https://github.com/FFinotello/quanTIseq
#'
#' F. Finotello, C. Mayer, C. Plattner, G. Laschober, D. Rieder,
#' H. Hackl, A. Krogsdam, W. Posch, D. Wilflingseder, S. Sopper, M. Jsselsteijn,
#' D. Johnsons, Y. Xu, Y. Wang, M. E. Sanders, M. V. Estrada, P. Ericsson-Gonzalez,
#' J. Balko, N. F. de Miranda, Z. Trajanoski. "quanTIseq: quantifying immune contexture of human tumors".
#' bioRxiv 223180. https://doi.org/10.1101/223180.
#'
#' @param mix.mat table with the gene TPM (or microarray expression values) for all samples to be deconvoluted
#'     (Gene symbols on the first column and sample IDs on the first row). Expression data must be on non-log scale
#' @param arrays specifies whether expression data are from microarrays (instead of RNA-seq).
#'     If TRUE, the "--rmgenes" parameter is set to "none".
#' @param signame name of the signature matrix. Currently only `TIL10` is available.
#' @param tumor 	specifies whether expression data are from tumor samples. If TRUE, signature genes
#'     with high expression in tumor samples are removed.
#'     Default: FALSE.
#' @param mRNAscale specifies whether cell fractions must be scaled to account for cell-type-specific
#'     mRNA content.
#'     Default: TRUE.
#' @param method deconvolution method to be used: "hampel", "huber", or "bisquare" for robust regression
#'     with Huber, Hampel, or Tukey bisquare estimators, respectively, or "lsei" for constrained
#'     least squares regression. The fraction of uncharacterized cells ("other") is computed only
#'     by the "lsei" method.
#'     Default: "lsei".
#' @param btotalcells compute cell densities instead of fractions
#'     Default: FALSE
#' @param rmgenes Default: "default" for RNAseq, "none" for microArray data
#'
#' @import preprocessCore
#' @export
deconvolute_quantiseq.default = function(mix.mat,
                                         arrays=FALSE,
                                         signame="TIL10",
                                         tumor=FALSE,
                                         mRNAscale=TRUE,
                                         method="lsei",
                                         btotalcells=FALSE,
                                         rmgenes="unassigned") {

  message("\nRunning quanTIseq deconvolution module\n")

  # List of genes to be discarded
  if (rmgenes=="unassigned" && arrays==TRUE) { # For Microarrays
    rmgenes<-"none"

  } else if (rmgenes=="unassigned" && arrays==FALSE) { # For RNA-seq
    rmgenes<-"default"

  }

  # Files
  listsig<-c("TIL10")
  if (signame %in% listsig) {

    sig.mat.file <- system.file("extdata", "quantiseq", paste0(signame, "_signature.txt"),
                                package="immunedeconv", mustWork=TRUE)

    mRNA.file <- system.file("extdata", "quantiseq", paste0(signame, "_mRNA_scaling.txt"),
                             package="immunedeconv", mustWork=TRUE)

    fileab <- system.file("extdata", "quantiseq", paste0(signame,"_TCGA_aberrant_immune_genes.txt"),
                          package="immunedeconv", mustWork=TRUE)

    if (rmgenes=="default") {
      filerm <- system.file("extdata", "quantiseq", paste0(signame,"_rmgenes.txt"),
                            package="immunedeconv", mustWork=TRUE)
    } else if (rmgenes=="path") {
      filerm <- system.file("extdata", "quantiseq", paste0(signame,"rmgenes.txt"),
                            package="immunedeconv", mustWork=TRUE)
    }


  } else {

    sig.mat.file<-paste0(signame, "_signature.txt")
    mRNA.file<-paste0(signame, "_mRNA_scaling.txt")

  }

  if(is.numeric(mix.mat[[1,1]])!= TRUE){
    stop("Wrong input format for the mixture matrix! Please follow the instructions of the documentation.")
  }

  # Load signature
  sig.mat<-read.table(sig.mat.file, header=TRUE, sep="\t", row.names=1)

  # Load normalization factors (set all to 1 if mRNAscale==FALSE)
  if (mRNAscale) {

    mRNA<-read.table(mRNA.file,
                     sep="\t",
                     header=FALSE,
                     stringsAsFactors=FALSE)
    colnames(mRNA)<-c("celltype", "scaling")
    mRNA<-as.vector(as.matrix(mRNA$scaling[match(colnames(sig.mat), mRNA$celltype)]))

  } else {

    mRNA<-rep(1, ncol(sig.mat))

  }

  # Preprocess mixture matrix
  message(paste0("Gene expression normalization and re-annotation (arrays: ",
    arrays, ")\n"))
  mix.mat<-fixMixture(mix.mat, arrays=arrays)

  # Remove noisy genes
  if (rmgenes!="none") {

    if (signame %in% listsig) {

      lrmgenes<-as.vector(read.table(filerm, header=FALSE, sep="\t")[,1])
      n1<-nrow(sig.mat)
      sig.mat<-sig.mat[!rownames(sig.mat) %in% lrmgenes,, drop=FALSE]
      n2<-nrow(sig.mat)
      message(paste0("Removing ", n1-n2, " noisy genes\n"))

    }
  }

  # Fix tumor data
  if (tumor) {

    if (signame %in% listsig) {

    abgenes<-as.vector(read.table(fileab, header=FALSE, sep="\t")[,1])
    n1<-nrow(sig.mat)
    sig.mat<-sig.mat[!rownames(sig.mat) %in% abgenes,, drop=FALSE]
    n2<-nrow(sig.mat)
    message(paste0("Removing ", n1-n2, " genes with high expression in tumors\n"))

    }
  }

  # Signature genes present in the mixture
  ns<-nrow(sig.mat)
  us<-length(intersect(rownames(sig.mat), rownames(mix.mat)))
  perc<-round(us*100/ns,digits=2)
  message(paste0("Signature genes found in data set: ",
    us, "/", ns, " (", perc, "%)\n"))

  # Run deconvolution
  message(paste0("Mixture deconvolution (method: ", method, ")\n"))
  results1<-quanTIseq(sig.mat,
                      mix.mat,
                      scaling=mRNA,
                      method=method)
  if ("Tregs" %in% colnames(sig.mat) && "T.cells.CD4" %in% colnames(sig.mat) && method %in% c("lsei")) {

    minTregs<-0.02
    i<-which(colnames(sig.mat)=="T.cells.CD4")
    results2<-quanTIseq(sig.mat[,-i],
      mix.mat,
      scaling=mRNA[-i],
      method=method)

    ind<-which(results1[,"Tregs"]<minTregs)
    if (length(ind)>0) {

      results1[ind,"Tregs"]<-(results2[ind,"Tregs"]+results1[ind,"Tregs"])/2
      results1[ind,"T.cells.CD4"]<-pmax(0,results1[ind,"T.cells.CD4"]-(results2[ind,"Tregs"]+results1[ind,"Tregs"])/2)

    }

  }
  results<-results1
  results<-results/apply(results,1,sum)

  message("Deconvolution sucessful!")



  # Save results using user's output ID
  DCres <- results


  if (btotalcells == TRUE) {
    celldens <- data.frame(celldensities(DCres))
    # fileout2 <- paste0(output, prefix, "_cell_densities.txt")
    celldens<-cbind(rownames(celldens), celldens)
    colnames(celldens)[1]<-"Sample"
    return(celldens)
    # write.table(celldens,
    #             sep="\t",
    #             row.names=FALSE,
    #             quote=FALSE,
    #             file=fileout2)
  } else {
    # fileout<-paste0(output, prefix, "_cell_fractions.txt")
    # cast to dataframe, otherwise cbind will cast to character.
    results = data.frame(results)
    results<-cbind(rownames(results), results)
    colnames(results)[1]<-"Sample"
    return(results)
    # write.table(results,
    #             sep="\t",
    #             row.names=FALSE,
    #             quote=FALSE,
    #             file=fileout)
  }

}
