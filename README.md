![immunedeconv](man/figures/immunedeconv_logo_sm.png)

**an R package for unified access to computational methods for estimating immune cell fractions from bulk RNA sequencing data.**

![tests](https://github.com/icbi-lab/immunedeconv/workflows/tests/badge.svg)
![test-conda](https://github.com/icbi-lab/immunedeconv/workflows/conda/badge.svg)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/r-immunedeconv/README.html)
[![license](https://img.shields.io/badge/license-BSD-green.svg)](https://github.com/icbi-lab/immunedeconv/blob/master/LICENSE.md)
[![docs](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://icbi-lab.github.io/immunedeconv)

## Basic usage

Deconvolution of human data:

```R
immunedeconv::deconvolute(gene_expression_matrix, "quantiseq")
```

Deconvolution of mouse data:

```R
immunedeconv::deconvolute_mouse(gene_expression_matrix, "mmcp_counter")
```

where `gene_expression_matrix` is a matrix with genes in rows and samples in columns. The rownames must be
[HGNC](https://www.genenames.org/) symbols for human data, or [MGI](http://www.informatics.jax.org/mgihome/nomen/) gene symbols for mouse data. 
The colnames must be sample names. For human data, the method can be one of

```
quantiseq
timer
cibersort
cibersort_abs
mcp_counter
xcell
epic
abis
consensus_tme
estimate
```

The [ESTIMATE](https://bioinformatics.mdanderson.org/public-software/estimate/) algorithm, which computes a score for the tumoral, immune and stromal components and the fraction of tumor purity of a sample, has been implemented. 

```R
immunedeconv::deconvolute_estimate(gene_expression_matrix)
```

The methods available for the deconvolution of mouse data are 

```
mmcp_counter
seqimmucc
dcq
base
```

In addition, human-based methods can be used to deconvolute mouse data through the conversion to orthologous gene names


```r
gene_expression_matrix <- immunedeconv::mouse_genes_to_human(gene_expression_matrix)
immunedeconv::deconvolute(gene_expression_matrix, "quantiseq")
```

Finally, certain methods can be used with custom signatures, consisting of either a signature matrix or signature genes 
for the cell types of interest. Since the information used to deconvolute the bulk is user-provided, these functions can be 
used for different tissues and organisms. 
The functions may require different input data formats, related to the requirements of each method. Please refer to their documentation. 
The available methods are


```r
base:  deconvolute_base_custom()
cibersort norm/abs:  deconvolute_cibersort_custom()
epic: deconvolute_epic_custom()
consensus_tme: deconvolute_consensus_tme_custom()
```



For more detailed usage instructions, see the Documentation:
* [Getting started](https://omnideconv.org/immunedeconv/articles/immunedeconv.html).
* [Detailed example](https://omnideconv.org/immunedeconv/articles/detailed_example.html).
* [Detailed example - mouse](https://omnideconv.org/immunedeconv/articles/detailed_example_mouse.html).


## Available methods, Licenses, Citations
Note that, while *immunedeconv* itself is free ([BSD](https://github.com/icbi-lab/immunedeconv/blob/master/LICENSE)), you may need to obtain a license to use the individual methods. See the table below for more information. If you use this package in your work, please cite both our package and the method(s) you are using.

> Sturm, G., Finotello, F., Petitprez, F., Zhang, J. D., Baumbach, J., Fridman, W. H., ..., List, M., Aneichyk, T. (2019). Comprehensive evaluation of transcriptome-based cell-type quantification methods for immuno-oncology. Bioinformatics, 35(14), i436-i445. https://doi.org/10.1093/bioinformatics/btz363



| method | organism |license | citation |
|--------|---------|----------|----------|
| [quanTIseq](http://icbi.at/software/quantiseq/doc/index.html) | human | free ([BSD](https://github.com/icbi-lab/immunedeconv/blob/master/LICENSE.md)) | Finotello, F., Mayer, C., Plattner, C., Laschober, G., Rieder, D., Hackl, H., ..., Sopper, S. (2019). Molecular and pharmacological modulators of the tumor immune contexture revealed by deconvolution of RNA-seq data. Genome medicine, 11(1), 34. https://doi.org/10.1186/s13073-019-0638-6 |
| [TIMER](http://cistrome.org/TIMER/) | human | free ([GPL 2.0](http://cistrome.org/TIMER/download.html)) | Li, B., Severson, E., Pignon, J.-C., Zhao, H., Li, T., Novak, J., … Liu, X. S. (2016). Comprehensive analyses of tumor immunity: implications for cancer immunotherapy. Genome Biology, 17(1), 174.  https://doi.org/10.1186/s13059-016-1028-7 |
| [CIBERSORT](https://cibersort.stanford.edu/) | human | free for non-commerical use only | Newman, A. M., Liu, C. L., Green, M. R., Gentles, A. J., Feng, W., Xu, Y., … Alizadeh, A. A. (2015). Robust enumeration of cell subsets from tissue expression profiles. Nature Methods, 12(5), 453–457.  https://doi.org/10.1038/nmeth.3337 |
| [MCPCounter](https://github.com/ebecht/MCPcounter) | human | free ([GPL 3.0](https://github.com/ebecht/MCPcounter/blob/master/Source/License)) | Becht, E., Giraldo, N. A., Lacroix, L., Buttard, B., Elarouci, N., Petitprez, F., … de Reyniès, A. (2016). Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression. Genome Biology, 17(1), 218. https://doi.org/10.1186/s13059-016-1070-5 |
| [xCell](http://xcell.ucsf.edu/) | human | free ([GPL 3.0](https://github.com/dviraran/xCell/blob/master/DESCRIPTION)) | Aran, D., Hu, Z., & Butte, A. J. (2017). xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome Biology, 18(1), 220. https://doi.org/10.1186/s13059-017-1349-1 |
| [EPIC](https://gfellerlab.shinyapps.io/EPIC_1-1/) | human | free for non-commercial use only ([Academic License](https://github.com/GfellerLab/EPIC/blob/master/LICENSE)) | Racle, J., de Jonge, K., Baumgaertner, P., Speiser, D. E., & Gfeller, D. (2017). Simultaneous enumeration of cancer and immune cell types from bulk tumor gene expression data. ELife, 6, e26476. https://doi.org/10.7554/eLife.26476 |
| [ESTIMATE](https://gfellerlab.shinyapps.io/EPIC_1-1/) | human | free ([GPL 2.0](https://bioinformatics.mdanderson.org/public-software/estimate/)) | Yoshihara, K., Shahmoradgoli, M., Martínez, E., Vegesna, R., Kim, H., Torres-Garcia, W., Treviño, V., Shen, H., Laird, P. W., Levine, D. A., Carter, S. L., Getz, G., Stemke-Hale, K., Mills, G. B., & Verhaak, R. G. (2013). Inferring tumour purity and stromal and immune cell admixture from expression data. Nature communications, 4, 2612. https://doi.org/10.1038/ncomms3612 |
| [ABIS](https://giannimonaco.shinyapps.io/ABIS/) | human | free ([GPL 2.0](https://github.com/giannimonaco/ABIS)) | Monaco, G., Lee, B., Xu, W., Mustafah, S., Hwang, Y. Y., ..., Larbi, A. (2019). RNA-Seq Signatures Normalized by mRNA Abundance Allow Absolute Deconvolution of Human Immune Cell Types. Cell reports, 26(6), 1627–1640.e7. https://doi.org/10.1016/j.celrep.2019.01.041|
| [ConsensusTME](https://olliecast.shinyapps.io/Deconvolution_Benchmarking/) | human | free ([GPL 3.0](https://github.com/cansysbio/ConsensusTME/blob/master/LICENSE.md)) | Jiménez-Sánchez, A., Cast, O., & Miller, M. L. (2019). Comprehensive Benchmarking and Integration of Tumor Microenvironment Cell Estimation Methods. Cancer research, 79(24), 6238–6246. https://doi.org/10.1158/0008-5472.CAN-18-3560 |
|[mMCPCounter](https://github.com/cit-bioinfo/mMCP-counter)| mouse | free ([GPL 3.0](https://github.com/cit-bioinfo/mMCP-counter/blob/master/LICENSE.md))| Petitprez, F., Levy, S., Sun, C. M., Meylan, M., ..., de Reyniès, A. (2020). The murine Microenvironment Cell Population counter method to estimate abundance of tissue-infiltrating immune and stromal cell populations in murine samples using gene expression. Genome medicine, 12(1), 86. https://doi.org/10.1186/s13073-020-00783-w |
|[seqImmuCC](218.4.234.74:3200/immune/)| mouse | free for non-commerical use only | Chen, Z., Quan, L., Huang, A., Zhao, Q., Yuan, Y., Yuan, X., ..., Wu, A. (2018). seq-ImmuCC: Cell-Centric View of Tissue Transcriptome Measuring Cellular Compositions of Immune Microenvironment From Mouse RNA-Seq Data. Frontiers in immunology, 9, 1286. https://doi.org/10.3389/fimmu.2018.01286 |
|[DCQ](http://dcq.tau.ac.il/)| mouse | free ([GPL 2.0](https://cran.r-project.org/web/packages/ComICS/index.html))| Altboum, Z., Steuerman, Y., David, E., Barnett-Itzhaki, Z., Valadarsky, L., ..., Amit, I. (2014). Digital cell quantification identifies global immune cell dynamics during influenza infection. Molecular systems biology, 10(2), 720. https://doi.org/10.1002/msb.134947 |
|BASE| mouse | free | Varn, F. S., Andrews, E. H., Mullins, D. W., & Cheng, C. (2016). Integrative analysis of breast cancer reveals prognostic haematopoietic activity and patient-specific immune response profiles. Nature communications, 7, 10248. https://doi.org/10.1038/ncomms10248  |

### Comparison of the methods
For a benchmark comparison of the human-based methods, please see our [publication](https://doi.org/10.1101/463828).
If you would like to benchmark additional methods, please see our [benchmark
pipeline](https://github.com/icbi-lab/immune_deconvolution_benchmark).


## Installation
System requirements: R >= 4.1. Only linux is officially supported, but Mac/Windows should work, too.

### Bioconda (Linux/MacOS only)
The easiest way to retrieve this package and all its dependencies is to use [Anaconda](https://conda.io/miniconda.html).
The installation typically completes within minutes.

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
as a regular R package in your default R installation. The installation typically completes within 30 minutes, depending
on how many dependency packages need to be compiled.

The easiest way to do so is to use the `remotes` package, which will automatically download all CRAN, Bioconductor and GitHub dependencies:

```R
install.packages("remotes")
remotes::install_github("omnideconv/immunedeconv")
```

## Credits

This package was originally developed by [Gregor Sturm](https://github.com/grst) in 2018 at [Pieris Pharmaceuticals GmbH](https://www.pieris.com/) in collaboration with [Markus List](https://biomedical-big-data.de/), [Tatsiana Aneichyk](https://www.independentdatalab.com/team), and [Francesca Finotello](https://computationalbiomedicinegroup.github.io/). Gregor Sturm continued to support this package while at [icbi-lab](https://icbi.at). In 2022, this repository moved to the [omnideconv](https://omnideconv.org) organization, a joint effort of the [List Lab](https://biomedical-big-data.de/) and [Finotello Lab](https://computationalbiomedicinegroup.github.io/) dedicated to improve accessibility of deconvolution methods. At this point [Lorenzo Merotto](https://github.com/LorenzoMerotto) became primary maintainer of the immunedeconv package. 

