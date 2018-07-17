# immune_deconvolution_methods (`immunedeconv`)

This package provides a unified interface to computational methods for estimating immune cell fractions from bulk RNA sequencing data.

```
immunedeconv::deconvolute(gene_expression_matrix, "quanTIseq")
```

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
For a benchmark comparison of the methods, please see our [publication](TODO).
If you would like to benchmark additional methods, please see our [benchmark
pipeline](https://github.com/grst/immune_deconvolution_benchmark).


## Installation
### Conda
The easiest way to retrieve this package and all its dependencies is to use [Anaconda](https://conda.io/miniconda.html).

1. If you have not already, download [Miniconda](https://conda.io/miniconda.html).

2. Create an environment for deconvolution: `conda create -n deconvolution`

3. Activate the environment `conda activate deconvolution`

4. Add additional conda channels:
```
conda config --add channels grst
conda config --add channels bioconda
conda config --add channels r
```

5. Install the `immunedeconv` package
```
conda install -c grst r-immunedeconv
```


### Standard R Package
You can also install `immunedeconv` as a regular R package in your environment.
You need install the following non-CRAN dependencies. If you use a very recent version of
[devtools](https://github.com/r-lib/devtools), it will also resolve these dependencies automatically.

**Bioconductor**
```R
biocLite("proprocessCore")
biocLite("Biobase")
biocLite("GSVA")
biocLite("sva")
biocLite("GSEABase")
```

**GitHub**
```R
library(devtools)
install_github('dviranran/xCell')
install_github('GfellerLab/EPIC')
install_github('ebecht/MCPcounter/Source')
```

Finally, install the `immunedeconv` package itself by running
```R
devtools::install_github('grst/immune_deconvolution_methods')
```

## Usage


## Developing
* add your own method

## Benchmark
see immune_deconvolution_benchmarks

## Citation
If you use this package, please cite

Sturm, G. and Aneichyk, T. "Benchmarking methods for estimating immune cell abundance from bulk RNA-sequencing data",
*manuscript in preparation*
