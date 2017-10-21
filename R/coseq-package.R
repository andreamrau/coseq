#' Co-expression and co-abundance analysis of high-throughput sequencing data
#'
#' Co-expression analysis for expression profiles arising from high-throughput sequencing
#' data. Feature (e.g., gene) profiles are clustered using adapted transformations and
#' mixture models or a K-means algorithm, and model selection criteria
#' (to choose an appropriate number of clusters) are provided.
#'
#' \tabular{ll}{ Package: \tab coseq\cr Type: \tab Package\cr Version:
#' \tab 1.1.3\cr Date: \tab 2017-10-20\cr License: \tab GPL (>=3)\cr LazyLoad:
#' \tab yes\cr }
#'
#' @name coseq-package
#' @aliases coseq-package
#' @docType package
#' @author Andrea Rau, Cathy Maugis-Rabusseau, Antoine Godichon-Baggioni
#'
#' Maintainer: Andrea Rau <\url{andrea.rau@@jouy.inra.fr}>
#' @references
#' Godichon-Baggioni, A., Maugis-Rabusseau, C. and Rau, A. (2017) Clustering
#' transformed compositional data using K-means, with applications in gene
#' expression and bicycle sharing system data. arXiv:1704.06150.
#'
#' Rau, A. and Maugis-Rabusseau, C. (2017) Transformation and model choice for
#' co-expression analayis of RNA-seq data. Briefings in Bioinformatics,
#' doi: http://dx.doi.org/10.1101/065607.
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



#' RNA-seq data from the mouse neocortex in Fietz et al. (2012)
#'
#' This dataset represents RNA-seq data from mouse neocortex RNA-seq data in five
#' embryonic (day 14.5) mice by analyzing the transcriptome of three regions: the
#' ventricular zone (VZ), subventricular zone (SVZ) and cortical place (CP).
#'
#' @name fietz
#' @docType data
#' @references \url{https://perso.math.univ-toulouse.fr/maugis/mixstatseq/packages}
#' @usage data(fietz)
#' @keywords datasets
#' @format An ExpressionSet named \code{fietz.eset} containing the phenotype data and
#' expression data for the Fietz et al. (2012) experiment. Phenotype data may be
#' accessed using the \code{pData} function, and expression data may be accessed
#' using the \code{exprs} function.
#' @return Object of class \sQuote{ExpressionSet}. Matrix of counts can be accessed after
#' loading the \sQuote{Biobase} package and calling \code{exprs(fietz))}.
#' @source Digital Expression Explorer (http://dee.bakeridi.edu.au/).
#' @references
#' Fietz, S. A., et al. (2012). Transcriptomes
#' of germinal zones of human and mouse fetal neocortex suggest a role of extracellular
#' matrix in progenitor self-renewal. Proceedings of the National Academy of Sciences,
#' 109(29):11836-11841.
NULL