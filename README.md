# coseq: Co-Expression Analysis of Sequencing Data

Authors: Andrea Rau, Cathy Maugis-Rabusseau, and Antoine Godichon-Baggioni

Co-expression analysis for expression profiles arising from high-throughput sequencing data. Feature (e.g., gene) profiles 
are clustered using adapted transformations and mixture models or a K-means algorithm, and model selection criteria
(to choose an appropriate number of clusters) are provided.

A typical call to coseq to apply the K-means algorithm to logCLR-transformed normalized
RNA-seq profiles takes the following form:
```
library(coseq)
## The following two lines are equivalent:
run_kmeans <- coseq(counts, K=2:10)
run_kmeans <- coseq(counts, K=2:10, model="kmeans", transformation="logclr")
```
where `counts` represents a (n×q) matrix or data frame of read counts for n genes in q samples
and `K=2:10` provides the desired range of numbers of clusters (here, 2 to 10). 

A typical call to coseq to fit a Gaussian mixture model on arcsin- or logit-transformed normalized
RNA-seq profiles takes the following form:
```
run_arcsin <- coseq(counts, K=2:10, model="Normal", transformation="arcsin")
run_logit <- coseq(counts, K=2:10, model="Normal", transformation="logit")
```
We note that this function directly calls the [Rmixmod](https://cran.r-project.org/web/packages/Rmixmod/index.html) 
R package to fit Gaussian mixture models. 

The output of the `coseq` function is an
S4 object of class `coseqResults` on which standard `plot` and `summary` functions can be directly applied; the former
uses functionalities from the [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) package. The option of parallelization
via the [BiocParallel](https://bioconductor.org/packages/release/bioc/html/BiocParallel.html) Bioconductor package is also provided.

### Reference

Rau, A. and Maugis-Rabusseau, C. (2018) Transformation and model choice for co-expression analysis of RNA-seq data. *Briefings in Bioinformatics*, 19(3)-425-436.

Godichon-Baggioni, A., Maugis-Rabusseau, C. and Rau, A. (2018) Clustering transformed compositional data using K-means, with applications in gene expression and bicycle sharing system data. *Journal of Applied Statistics*, doi:10.1080/02664763.2018.1454894.

### License

The coseq package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public 
License, version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but without any warranty; without even the implied 
warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.

A copy of the GNU General Public License, version 3, is available at http://www.r-project.org/Licenses/GPL-3.