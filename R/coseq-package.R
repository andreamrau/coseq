#' Co-expression and co-abundance analysis of high-throughput sequencing data
#'
#' Mixture models are implemented to cluster genes from high-throughput
#' transcriptome sequencing (RNA-seq) data. Parameter estimation is performed
#' using the EM algorithm, and model selection is performed using
#' either the slope heuristics or the integrated completed likelihood (ICL)
#' criterion.
#'
#' \tabular{ll}{ Package: \tab coseq\cr Type: \tab Package\cr Version:
#' \tab 0.99.7\cr Date: \tab 2016-12-20\cr License: \tab GPL (>=3)\cr LazyLoad:
#' \tab yes\cr }
#'
#' @name coseq-package
#' @aliases coseq-package
#' @docType package
#' @author Andrea Rau, Cathy Maugis-Rabusseau
#'
#' Maintainer: Andrea Rau <\url{andrea.rau@@jouy.inra.fr}>
#' @references
#' Rau, A. and Maugis-Rabusseau, C. (2016) Transformation and model choice for
#' co-expression analayis of RNA-seq data. bioRxiv, doi: http://dx.doi.org/10.1101/065607.
#'
#' Rau, A., Maugis-Rabusseau, C., Martin-Magniette, M.-L., Celeux,
#' G. (2015) Co-expression analysis of high-throughput transcriptome sequencing
#' data with Poisson mixture models. Bioinformatics, doi:
#' 10.1093/bioinformatics/btu845.
#'
#' Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C. (2011)
#' Clustering high-throughput sequencing data with Poisson mixture models.
#' Inria Research Report 7786. Available at
#' \url{http://hal.inria.fr/inria-00638082}.
#' @keywords models cluster
#' @example /inst/examples/coseq-package.R
NULL