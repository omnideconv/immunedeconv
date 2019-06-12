# `immunedeconv` - an R package for unified access to computational methods for estimating immune cell fractions from bulk RNA sequencing data.

[![travis](https://travis-ci.com/grst/immunedeconv.svg?branch=master)](https://travis-ci.com/grst/immunedeconv) ![license](https://img.shields.io/badge/license-BSD-green.svg) [![docs](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://grst.github.io/immunedeconv)

## Basic usage
```R
immunedeconv::deconvolute(gene_expression_matrix, "quantiseq")
```

where `gene_expression_matrix` is a matrix with genes in rows and samples in columns. The rownames must be
[HGNC](https://www.genenames.org/) symbols and the colnames must be sample names. The method can be one of
```
quantiseq
timer
cibersort
cibersort_abs
mcp_counter
xcell
epic
```

For more detailed usage instructions, see the
[Documentation](https://grst.github.io/immunedeconv/articles/immunedeconv.html).


## Available methods

| method | citation |
|--------|----------|
| [quanTIseq](http://icbi.at/software/quantiseq/doc/index.html) | Finotello, F., Mayer, C., Plattner, C., Laschober, G., Rieder, D., Hackl, H., … Trajanoski, Z. (2017). quanTIseq: quantifying immune contexture of human tumors. BioRxiv, 223180. https://doi.org/10.1101/223180 |
| [TIMER](http://cistrome.org/TIMER/) | Li, B., Severson, E., Pignon, J.-C., Zhao, H., Li, T., Novak, J., … Liu, X. S. (2016). Comprehensive analyses of tumor immunity: implications for cancer immunotherapy. Genome Biology, 17(1), 174.  https://doi.org/10.1186/s13059-016-1028-7 |
| [CIBERSORT](https://cibersort.stanford.edu/) | Newman, A. M., Liu, C. L., Green, M. R., Gentles, A. J., Feng, W., Xu, Y., … Alizadeh, A. A. (2015). Robust enumeration of cell subsets from tissue expression profiles. Nature Methods, 12(5), 453–457.  https://doi.org/10.1038/nmeth.3337 |
| [MCPCounter](https://github.com/ebecht/MCPcounter) | Becht, E., Giraldo, N. A., Lacroix, L., Buttard, B., Elarouci, N., Petitprez, F., … de Reyniès, A. (2016). Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. Genome Biology, 17(1), 218. https://doi.org/10.1186/s13059-016-1070-5 |
| [xCell](http://xcell.ucsf.edu/) | Aran, D., Hu, Z., & Butte, A. J. (2017). xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome Biology, 18(1), 220. https://doi.org/10.1186/s13059-017-1349-1 |
| [EPIC](https://gfellerlab.shinyapps.io/EPIC_1-1/) | Racle, J., de Jonge, K., Baumgaertner, P., Speiser, D. E., & Gfeller, D. (2017). Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. ELife, 6, e26476. https://doi.org/10.7554/eLife.26476 |


### Comparison of the methods
For a benchmark comparison of the methods, please see our [publication](https://doi.org/10.1101/463828).
If you would like to benchmark additional methods, please see our [benchmark
pipeline](https://github.com/grst/immune_deconvolution_benchmark).


## Installation
System requirements: linux and R >= 3.4.1. The package has been tested on CentOS Linux 7.5 with R 3.4.1. 

### Conda
The easiest way to retrieve this package and all its dependencies is to use [Anaconda](https://conda.io/miniconda.html).
Install typically completes within minutes. 

1. Download [Miniconda](https://conda.io/miniconda.html), if you don't have a conda installation already.

2. (Optional) create and activate an environment for deconvolution:
```
conda create -n deconvolution
conda activate deconvolution
```

3. Install the `immunedeconv` package
```
conda install --override -c grst -c bioconda -c conda-forge/label/cf201901 r-immunedeconv
```
**Note:** due to a recent [conda compiler update](https://github.com/conda/conda/issues/8413), immunedeconv needs
to be installed using the `cf201901` label of conda-forge. I'm working on making it work with the default version...

`conda` will automatically install the package and all dependencies.
You can then open an `R` instance within the environment and use the package.


### Standard R Package
We highly recommend using `conda`, as it will avoid incompatibilities between
different package versions. That being said, you can also install `immunedeconv`
as a regular R package in your default R installation. Installation typically completes within 30 minutes, depending 
on how many dependency packages need to be compiled. 

The easiest way to do so is to use the `remotes` package, which will automatically download all CRAN, Bioconductor and GitHub dependencies:

```R
install.packages("remotes")
remotes::install_github("grst/immunedeconv")
```

## Citation
If you use this package, please cite

Sturm, G., Finotello F. et al. "Comprehensive evaluation of cell-type quantification methods for immuno-oncology",
bioRxiv, https://doi.org/10.1101/463828
