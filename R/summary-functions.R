#' Summarize results from coseq clustering
#'
#' A function to summarize the clustering results obtained from a Poisson or
#' Gaussian mixture model estimated using \code{coseq}. In particular,
#' the function provides the number of clusters selected for the ICL
#' model selection approach (or alternatively, for the capushe non-asymptotic approach
#' if K-means clustering is used), number of genes assigned to each cluster, and
#' if desired the per-gene cluster means.
#'
#' Provides the following summary of results:
#'
#' 1) Number of clusters and model selection criterion used, if applicable.
#'
#' 2) Number of observations across all clusters with a maximum conditional
#' probability greater than 90% (and corresponding percentage of total
#' observations) for the selected model.
#'
#' 3) Number of observations per cluster with a maximum conditional probability
#' greater than 90% (and corresponding percentage of total observations per
#' cluster) for the selected model.
#'
#' 4) If desired, the \eqn{\ensuremath\boldsymbol{\mu}}{\mu} values and
#' \eqn{\ensuremath\boldsymbol{\pi}}{\pi} values for the selected
#' model in the case of a Gaussian mixture model.
#'
#' @rdname summary
#' @aliases
#' summary
#' summary-methods
#' summary,coseqResults-method
#'
#' @param object An object of class \code{"coseqResults"}
#' @param y_profiles Data used for clustering if per-cluster means are desired
#' @param digits Integer indicating the number of decimal places to be used
#' for mixture model parameters
#' @param ... Additional arguments
#' @author Andrea Rau
#' @seealso \code{\link{coseq}}
#' @references
#' Rau, A. and Maugis-Rabusseau, C. (2017) Transformation and model choice for
#' co-expression analayis of RNA-seq data. Briefings in Bioinformatics,
#' doi: http://dx.doi.org/10.1101/065607.
#'
#' Godichon-Baggioni, A., Maugis-Rabusseau, C. and Rau, A. (2017) Clustering
#' transformed compositional data using K-means, with applications in gene
#' expression and bicycle sharing system data. arXiv:1704.06150.
#'
#' @return Summary of the \code{coseqResults} object.
#' @keywords methods
#' @example /inst/examples/coseq-package.R
#' @export
setMethod("summary", signature(object="coseqResults"), function(object, y_profiles, digits=3, ...) {
  x <- object
  if(!is(x, "coseqResults")) stop(paste0(sQuote("object")), " must be of class ",
         paste0(dQuote("coseqResults")), sep="")
  cat("*************************************************\n")
  if(model(x) == "Poisson") cat("Model: ", model(x), "\n", sep = "")
  if(model(x) == "Normal") cat("Model: ", metadata(x)$GaussianModel, "\n", sep="")
  if(model(x) == "kmeans") cat("K-means algorithm\n", sep="")
  cat("Transformation: ", transformationType(x), "\n", sep = "")

  cat("*************************************************\n")
  clustNum <-  paste(metadata(x)$nbCluster, collapse=",")
  clustErr <- paste(metadata(x)$nbClusterError, collapse=",")
  clustErr <- ifelse(clustErr == "", "---", clustErr)
  cat("Clusters fit: ", clustNum, "\n", sep = "")
  cat("Clusters with errors: ", clustErr, "\n", sep= "")
  if(!is.null(ICL(x))) {
    cat("Selected number of clusters via ICL: ", ncol(assay(x)), "\n", sep = "")
    cat("ICL of selected model: ", min(ICL(x)), "\n", sep = "")
  }
  if(is.null(ICL(x))) {
    cat("Selected number of clusters via capushe: ", ncol(assay(x)), "\n", sep = "")
  }
  cat("*************************************************\n")

  probaPost <- assay(x)
  labels <- apply(probaPost, 1, which.max)
  map <- apply(probaPost, 1, max)
  length(which(map > 0.9))/length(map)
  g <- ncol(probaPost)
  cat("Number of clusters = ", ncol(probaPost), "\n", sep = "")
  if(model(x) != "kmeans") {
    cat("ICL = ", min(ICL(x)), "\n", sep = "")
  }
  cat("*************************************************\n")
  tab <- table(labels)
  names(tab) <- paste("Cluster", names(tab))
  cat("Cluster sizes:\n")
  print(tab)
  cat("\n")
  if(model(x) != "kmeans") {
    cat("Number of observations with MAP > 0.90 (% of total):\n")
    cat(length(which(map > 0.9)), " (", round(length(which(map > 0.9))/length(map)*100,2),
        "%)\n\n", sep = "")
    cat("Number of observations with MAP > 0.90 per cluster (% of total per cluster):\n")

    tab2 <- matrix(NA, nrow = 2, ncol = g)
    colnames(tab2) <- paste("Cluster", 1:g)
    rownames(tab2) <- rep("", 2)
    for(i in seq_len(g)) {
      if(sum(labels == i) > 1) {
        map.clust <- apply(matrix(probaPost[labels == i,], ncol=g), 1, max)
        tab2[1,i] <- length(which(map.clust > 0.9))
        tab2[2,i] <- paste("(", round(100*length(which(map.clust > 0.9))/length(map.clust),2),
                           "%)", sep = "")
      }
      if(sum(labels == i) == 1) {
        map.clust <- max(probaPost[labels == i,])
        tab2[1,i] <- length(which(map.clust > 0.9))
        tab2[2,i] <- paste("(", round(100*length(which(map.clust > 0.9))/length(map.clust),2),
                           "%)", sep = "")
      }
      if(sum(labels == i) == 0) {
        tab2[1,i] <- "---"
        tab2[2,i] <- "---"
      }
    }
    print(tab2, quote = FALSE)
    cat("\n")
  }

  ## Print out parameter estimates if Gaussian model and y_profiles provided
  if(!is.null(metadata(x)$GaussianModel) & !missing(y_profiles)) {
    param <- NormMixParam(x, y_profiles, digits=digits)
    mu <- param$mu
    pi <- param$pi
    g <- length(pi)
    rownames(mu) <- names(pi) <- names(tab)
    colnames(mu) <- colnames(y_profiles)
    cat("Mu:\n")
    print(round(mu,digits=digits))
    cat("\n")
    cat("Pi:\n")
    print(round(pi,digits=digits))
    cat("\n")
  }

})