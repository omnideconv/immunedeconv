# `immunedeconv` - an R package for unified access to computational methods for estimating immune cell fractions from bulk RNA sequencing data.
[![travis](https://travis-ci.com/icbi-lab/immunedeconv.svg?branch=master)](https://travis-ci.com/icbi-lab/immunedeconv)
[![appveyor](https://ci.appveyor.com/api/projects/status/j2fb3fd097kqahg5/branch/master?svg=true)](https://ci.appveyor.com/project/icbi-lab/immunedeconv/branch/master)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/r-immunedeconv/README.html)
[![license](https://img.shields.io/badge/license-BSD-green.svg)](https://github.com/icbi-lab/immunedeconv/blob/master/LICENSE.md)
[![docs](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://icbi-lab.github.io/immunedeconv)

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

For more detailed usage instructions, see the Documentation:
* [Getting started](https://icbi-lab.github.io/immunedeconv/articles/immunedeconv.html).
* [Detailed example](https://icbi-lab.github.io/immunedeconv/articles/detailed_example.html).


## Available methods, Licenses, Citations
Note that, while *immunedeconv* itself is free ([BSD](https://github.com/icbi-lab/immunedeconv/blob/master/LICENSE.md)), you may need to obtain a license to use the individual methods. See the table below for more information. If you use this package in your work, please cite both our package and the method(s) you are using.

> Sturm, G., Finotello, F., Petitprez, F., Zhang, J. D., Baumbach, J., Fridman, W. H., ..., List, M., Aneichyk, T. (2019). Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology. Bioinformatics, 35(14), i436-i445. https://doi.org/10.1093/bioinformatics/btz363



| method | license | citation |
|--------|---------|----------|
| [quanTIseq](http://icbi.at/software/quantiseq/doc/index.html) | free ([BSD](https://github.com/icbi-lab/immunedeconv/blob/master/LICENSE.md)) | Finotello, F., Mayer, C., Plattner, C., Laschober, G., Rieder, D., Hackl, H., ..., Sopper, S. (2019). Molecular and pharmacological modulators of the tumor immune contexture revealed by deconvolution of RNA-seq data. Genome medicine, 11(1), 34. https://doi.org/10.1186/s13073-019-0638-6 |
| [TIMER](http://cistrome.org/TIMER/) | free ([GPL 2.0](http://cistrome.org/TIMER/download.html)) | Li, B., Severson, E., Pignon, J.-C., Zhao, H., Li, T., Novak, J., … Liu, X. S. (2016). Comprehensive analyses of tumor immunity: implications for cancer immunotherapy. Genome Biology, 17(1), 174.  https://doi.org/10.1186/s13059-016-1028-7 |
| [CIBERSORT](https://cibersort.stanford.edu/) | free for non-commerical use only | Newman, A. M., Liu, C. L., Green, M. R., Gentles, A. J., Feng, W., Xu, Y., … Alizadeh, A. A. (2015). Robust enumeration of cell subsets from tissue expression profiles. Nature Methods, 12(5), 453–457.  https://doi.org/10.1038/nmeth.3337 |
| [MCPCounter](https://github.com/ebecht/MCPcounter) | free ([GPL 3.0](https://github.com/ebecht/MCPcounter/blob/master/Source/License)) | Becht, E., Giraldo, N. A., Lacroix, L., Buttard, B., Elarouci, N., Petitprez, F., … de Reyniès, A. (2016). Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. Genome Biology, 17(1), 218. https://doi.org/10.1186/s13059-016-1070-5 |
| [xCell](http://xcell.ucsf.edu/) | free ([GPL 3.0](https://github.com/dviraran/xCell/blob/master/DESCRIPTION)) | Aran, D., Hu, Z., & Butte, A. J. (2017). xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome Biology, 18(1), 220. https://doi.org/10.1186/s13059-017-1349-1 |
| [EPIC](https://gfellerlab.shinyapps.io/EPIC_1-1/) | free for non-commercial use only ([Academic License](https://github.com/GfellerLab/EPIC/blob/master/LICENSE)) | Racle, J., de Jonge, K., Baumgaertner, P., Speiser, D. E., & Gfeller, D. (2017). Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. ELife, 6, e26476. https://doi.org/10.7554/eLife.26476 |


### Comparison of the methods
For a benchmark comparison of the methods, please see our [publication](https://doi.org/10.1101/463828).
If you would like to benchmark additional methods, please see our [benchmark
pipeline](https://github.com/icbi-lab/immune_deconvolution_benchmark).


## Installation
System requirements: R >= 3.4.1. Only linux is officially supported, but Mac/Windows should work, too.

### Bioconda (Linux/MacOS only)
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
conda install -c bioconda -c conda-forge r-immunedeconv
```

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
remotes::install_github("icbi-lab/immunedeconv")
```


