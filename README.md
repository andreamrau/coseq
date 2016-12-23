# coseq: Co-Expression Analysis of Sequencing Data

Authors: Andrea Rau and Cathy Maugis-Rabusseau

Gaussian and Poisson mixture models are implemented to cluster gene expression profiles from high-throughput sequencing data. 
Parameter estimation is performed using the EM algorithm and model selection criteria (to choose the number of clusters and
data transformation) are provided. 

A typical call to coseq to fit a Gaussian mixture model on arcsin- or logit-transformed normalized
RNA-seq profiles takes the following form:
```
library(coseq)
run_arcsin <- coseq(counts, K=2:10, model="Normal", transformation="arcsin")
run_logit <- coseq(counts, K=2:10, model="Normal", transformation="logit")
```
where `counts` represents a (n×q) matrix or data frame of read counts for n genes in q samples
and `K=2:10` provides the desired range of numbers of clusters (here, 2 to 10). We note that
this function directly calls the [Rmixmod](https://cran.r-project.org/web/packages/Rmixmod/index.html) R package to fit Gaussian mixture models. 
The output of the `coseq` function is an
S4 object of class `coseqResults` on which standard `plot` and `summary` functions can be directly applied; the former
uses functionalities from the [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) package. The option of parallelization
via the [BiocParallel](https://bioconductor.org/packages/release/bioc/html/BiocParallel.html) Bioconductor package is also provided.

### Reference

Rau, A. and Maugis-Rabusseau, C. (2016) Transformation and model choice for co-expression analayis of RNA-seq data. 
*Briefings in Bioinformatics* (to appear), doi: http://dx.doi.org/10.1101/065607.

### License

The coseq package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public 
License, version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied 
warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

A copy of the GNU General Public License, version 3, is available at http://www.r-project.org/Licenses/GPL-3.