#' Co-expression or co-abudance analysis of high-throughput sequencing data
#'
#' This is the primary user interface for the \code{coseq} package.
#' Generic S4 methods are implemented to perform co-expression or co-abudance analysis of
#' high-throughput sequencing data, with or without data transformation, using K-means or mixture models.
#' The supported classes are \code{matrix}, \code{data.frame}, and \code{DESeqDataSet}.
#' The output of \code{coseq} is an S4 object of class \code{coseqResults}.
#'
#' @inheritParams coseqRun
#' @param object Data to be clustered. May be provided as a y (\emph{n} x \emph{q})
#' matrix or data.frame of observed counts for \emph{n}
#' observations and \emph{q} variables, or an object of class \code{DESeqDataSet}
#' arising from a differential analysis via DESeq2.
#'
#' @return
#' An S4 object of class \code{coseqResults}, where conditional
#' probabilities of cluster membership for each gene in each model is stored as a SimpleList of assay
#' data, and the corresponding {log likelihood, ICL value, number of
#' clusters, and form of Gaussian model} for each model are stored as metadata.
#'
#' @aliases
#' coseq
#' coseq-methods
#' coseq,matrix-method
#' coseq,data.frame-method
#' coseq,DESeqDataSet-method
#'
#' @author Andrea Rau
#' @export
#' @example inst/examples/coseq-package.R
#' @keywords methods
#' @rdname coseq
#' @docType methods
#' @import methods
#' @importFrom stats p.adjust
setMethod("coseq",
          signature=signature(object="matrix"),
          definition=function(object, K, subset=NULL, model="kmeans", transformation="logclr",
                              normFactors="TMM", meanFilterCutoff=NULL,
                              modelChoice=ifelse(model=="kmeans", "DDSE", "ICL"),
                              parallel=FALSE,  BPPARAM=bpparam(), ...)
          {
            y <- object
            arg.user <- list(...)

            if(is.null(subset) | is.numeric(subset)) {
              run <- coseqRun(y=y, K=K, subset=subset, model=model, transformation=transformation,
                               normFactors=normFactors, meanFilterCutoff=meanFilterCutoff, modelChoice=modelChoice,
                               parallel=parallel,
                               BPPARAM=BPPARAM, ...)
            }
            if(is(subset, "DESeqResults")) {
              if(nrow(subset) != nrow(y)) stop("y and subset must have the same number of rows")
              res <- subset
              subset.index <- which(res$padj < metadata(res)$alpha)
              cat("****************************************\n")
              cat("Co-expression analysis on DESeq2 output:\n")
              cat(paste(length(subset.index), "DE genes at p-adj <", metadata(res)$alpha, "\n"))
              run <- coseqRun(y=y, K=K, subset=subset.index, model=model, transformation=transformation,
                              normFactors="DESeq", meanFilterCutoff=NULL, modelChoice=modelChoice,
                               parallel=parallel,
                               BPPARAM=BPPARAM, ...)
            }
            if(is(subset, "DGELRT")) {
              res <- subset
              if(is.null(arg.user$alpha)) {
                alpha <- 0.05
              } else {
                alpha <- arg.user$alpha
              }
              subset.index <- which(p.adjust(subset$table$PValue, method="BH") < alpha)
              cat("****************************************\n")
              cat("Co-expression analysis on edgeR output:\n")
              cat(paste(length(subset.index), "DE genes at p-adj <", alpha, "\n"))
              run <- coseqRun(y=y, K=K, subset=subset.index, model=model, transformation=transformation,
                              normFactors="TMM", meanFilterCutoff=NULL, modelChoice=modelChoice,
                               parallel=parallel,
                               BPPARAM=BPPARAM, ...)
            }
            return(run)
          })



#########################################################################################
#' @rdname coseq
#' @export
setMethod("coseq", signature=signature(object="data.frame"),
          definition=function(object, K, subset=NULL, model="kmeans", transformation="logclr",
                              normFactors="TMM", meanFilterCutoff=NULL,
                              modelChoice=ifelse(model=="kmeans", "DDSE", "ICL"),
                              parallel=FALSE,  BPPARAM=bpparam(), ...)
          {
            y <- as.matrix(object)
            run <- coseq(y, K=K, subset=subset, model=model, transformation=transformation, normFactors=normFactors,
                         meanFilterCutoff=meanFilterCutoff, modelChoice=modelChoice, parallel=parallel,
                         BPPPARAM=BPPARAM, ...)
            return(run)
          })



#########################################################################################
#' @rdname coseq
#' @export
#' @importMethodsFrom DESeq2 counts sizeFactors
#' @importFrom DESeq2 results
setMethod("coseq", signature=signature(object="DESeqDataSet"),
          definition=function(object, K, model="kmeans", transformation="logclr",
                              normFactors="TMM", meanFilterCutoff=NULL,
                              modelChoice=ifelse(model=="kmeans", "DDSE", "ICL"),
                              parallel=FALSE, BPPARAM=bpparam(), ...)
          {
            ## Parse ellipsis function separately for DESeq and coseq
            dots <- dots_DESeq <- dots_coseq <- list(...)
            DESeq_argindex <- setdiff(names(dots_DESeq), c("alpha"))
            coseq_argindex <-
              setdiff(names(dots_coseq), c("conds", "geneIDs",
                                           "parallel", "BPPARAM", "alg.type", "init.runs",
                                           "init.type", "GaussianModel", "init.iter",
                                           "cutoff", "verbose", "digits", "fixed.lambda",
                                           "equal.proportions", "prev.labels",
                                           "prev.probaPost", "interpretation", "EM.verbose",
                                           "wrapper", "modelChoice"))
            dots_DESeq[DESeq_argindex] <- NULL
            dots_coseq[coseq_argindex] <- NULL

            y <- object
            count_matrix <- counts(y)
            if(is.null(dots_DESeq$alpha)) res <- results(y)
            if(!is.null(dots_DESeq$alpha)) res <- results(y, alpha=dots_DESeq$alpha)
            normFactors <- sizeFactors(y)
            subset.index <- which(res$padj < metadata(res)$alpha)
            cat("****************************************\n")
            cat("Co-expression analysis on DESeq2 output:\n")
            cat(paste(length(subset.index), "DE genes at p-adj <", metadata(res)$alpha, "\n"))
            meanFilterCutoff <- NULL

            run <- do.call(coseqRun, list(y=count_matrix, K=K, subset=subset.index, model=model,
                                          transformation=transformation,
                                          meanFilterCutoff=meanFilterCutoff, normFactors=normFactors, dots_coseq))
            return(run)

          })





#########################################################################################
#' Accessors for the assigned cluster labels of a coseqResults object.
#'
#' The counts slot holds the count data as a matrix of non-negative integer
#' count values, one row for each observational unit (gene or the like), and one
#' column for each sample.
#'
#' @docType methods
#' @rdname coseqHelpers
#' @aliases
#' clusters
#' clusters,coseqResults-method
#' clusters,RangedSummarizedExperiment-method
#' likelihood
#' likelihood,MixmodCluster-method
#' likelihood,RangedSummarizedExperiment-method
#' likelihood,coseqResults-method
#' likelihood,NULL-method
#' nbCluster
#' nbCluster,MixmodCluster-method
#' nbCluster,RangedSummarizedExperiment-method
#' nbCluster,coseqResults-method
#' nbCluster,NULL-method
#' ICL
#' ICL,RangedSummarizedExperiment-method
#' ICL,mixmodCluster-method
#' ICL,coseqResults-method
#' ICL,NULL-method
#' profiles
#' profiles,coseqResults-method
#' tcounts
#' tcounts,coseqResults-method
#' transformationType
#' transformationType,coseqResults-method
#' show
#' show,coseqResults-method
#' model
#' model,coseqResults-method
#' proba
#' proba,MixmodCluster-method
#' DDSEextract
#' DDSEextract,Capushe-method
#' Djumpextract
#' Djumpextract,Capushe-method
#' @param object a \code{coseqResults}, \code{RangedSummarizedExperiment}, or
#' \code{MixmodCluster} object.
#' @param K numeric indicating the model to be used (if NULL of missing, the model chosen
#' by ICL is used by default)
#' @param ... Additional optional parameters
#' @return Output varies depending on the method. \code{clusters} returns a vector of cluster
#' labels for each gene for the desired model.
#' @author Andrea Rau
#' @export
#' @example inst/examples/coseq-package.R
setMethod("clusters", signature(object="coseqResults"),
          function(object, K) {
            if(missing(K))
              pp <- assay(object)
            if(!missing(K))
              pp <- coseqFullResults(object)[[paste0("K=", K)]]
            labels <- apply(pp, 1, which.max)
            names(labels) <- rownames(pp)
            return(labels)
})

#' @rdname coseqHelpers
#' @export
setMethod("clusters", signature(object="RangedSummarizedExperiment"),
          function(object, ...) {
            pp <- assay(object)
            labels <- apply(pp, 1, which.max)
            return(labels)
})

#' @rdname coseqHelpers
#' @export
setMethod("clusters", signature(object="matrix"),
          function(object, ...) {
            labels <- apply(object, 1, which.max)
            return(labels)
})

#' @rdname coseqHelpers
#' @export
setMethod("clusters", signature(object="data.frame"),
          function(object, ...) {
            labels <- apply(object, 1, which.max)
            return(labels)
          })

#########################################################################################

#' @rdname coseqHelpers
#' @export
setMethod("likelihood", "MixmodCluster", function(object) object["bestResult"]@likelihood)

#' @rdname coseqHelpers
#' @export
setMethod("likelihood", "RangedSummarizedExperiment", function(object) metadata(object)$logLike)

#' @rdname coseqHelpers
#' @export
setMethod("likelihood", "coseqResults", function(object) metadata(object)$logLike)

#' @rdname coseqHelpers
#' @export
setMethod("likelihood", "NULL", function(object) NA)

#########################################################################################

#' @rdname coseqHelpers
#' @export
setMethod("nbCluster", "MixmodCluster", function(object) object["bestResult"]@nbCluster)

#' @rdname coseqHelpers
#' @export
setMethod("nbCluster", "RangedSummarizedExperiment", function(object) metadata(object)$nbCluster)

#' @rdname coseqHelpers
#' @export
setMethod("nbCluster", "coseqResults", function(object) metadata(object)$nbCluster)

#' @rdname coseqHelpers
#' @export
setMethod("nbCluster", "NULL", function(object) NA)

#########################################################################################

#' @rdname coseqHelpers
#' @export
setMethod("ICL", "RangedSummarizedExperiment", function(object) metadata(object)$ICL)

#' @rdname coseqHelpers
#' @export
setMethod("ICL", "MixmodCluster", function(object) object["bestResult"]@criterionValue)

#' @rdname coseqHelpers
#' @export
setMethod("ICL", "coseqResults", function(object) metadata(object)$ICL)

#' @rdname coseqHelpers
#' @export
setMethod("ICL", "NULL", function(object) NA)


#########################################################################################

#' @rdname coseqHelpers
#' @export
setMethod("profiles", "coseqResults", function(object) object@y_profiles)

#' @rdname coseqHelpers
#' @export
setMethod("tcounts", "coseqResults", function(object) object@tcounts)

#' @rdname coseqHelpers
#' @export
setMethod("transformationType", "coseqResults", function(object) object@transformation)

#' @rdname coseqHelpers
#' @export
setMethod("model", "coseqResults", function(object) object@model)

#' @rdname coseqHelpers
#' @export
setMethod("coseqFullResults", "coseqResults", function(object) object@allResults)

#########################################################################################

#' @rdname coseqHelpers
#' @export
setMethod("show", signature = "coseqResults", definition = function(object) {
  cat("An object of class ", class(object), "\n", sep="")
  if(nrow(profiles(object)))
    cat(" ", nrow(profiles(object)), " features by ", ncol(profiles(object)), " samples. \n", sep = "")
  if(!nrow(profiles(object)))
    cat(" ", nrow(object), " features by ", ncol(object), " clusters \n", sep = "")
  mods <- as.numeric(unlist(lapply(strsplit(names(coseqFullResults(object)), split="=",
                                            fixed=TRUE),
                        function(xx) xx[2])))
  cat(" ", "Models fit: K = ", min(mods), " ... ", max(mods), "\n", sep="")
  cat(" ", "Chosen clustering model: K = ", ncol(assay(object)), sep="")
  invisible(NULL)
})

#########################################################################################

#' @rdname coseqHelpers
#' @export
setMethod("proba", "MixmodCluster", function(object) object["bestResult"]@proba)
#' @rdname coseqHelpers
#' @export
#' @importClassesFrom capushe Capushe
setMethod("DDSEextract", "Capushe", function(object) object@DDSE@model)
#' @rdname coseqHelpers
#' @export
setMethod("Djumpextract", "Capushe", function(object) object@Djump@model)

#########################################################################################
#' Pairwise comparisons of ARI values among a set of clustering partitions
#'
#' Provides the adjusted rand index (ARI) between pairs of clustering paritions.
#'
#' @name compareARI
#' @rdname compareARI
#' @aliases
#' compareARI
#' compareARI-methods
#' compareARI,coseqResults-method
#' compareARI,RangedSummarizedExperiment-method
#' compareARI,matrix-method
#' compareARI,data.frame-method
#' @param object Object of class \code{coseqResults} or \code{RangedSummarizedExperiment},
#' or alternatively a \emph{n} x \emph{M} \code{data.frame} or \code{matrix}
#' containing the clustering partitions for \emph{M} different models
#' @param K If \code{NULL}, pairwise ARI values will be calculated among every model in
#' object \code{x}. Otherwise, \code{K} provides a vector of cluster numbers identifying
#' a subset of models in \code{x}.
#' @param parallel If \code{FALSE}, no parallelization. If \code{TRUE}, parallel
#' execution using BiocParallel (see next argument \code{BPPARAM}).
#' Note that parallelization is unlikely to be helpful unless the number of
#' observations \emph{n} in the clustering partitions or the number of
#' models \emph{M} are very large.
#' @param BPPARAM Optional parameter object passed internally to \code{bplapply}
#' when \code{parallel=TRUE}. If not specified, the parameters last registered
#' with \code{register} will be used.
#' @param plot If \code{TRUE}, provide a heatmap using corrplot to visualize
#' the calculated pairwise ARI values.
#' @param ... Additional optional parameters for corrplot
#'
#' @return Matrix of adjusted rand index values calculated between each pair
#' of models.
#'
#' @author Andrea Rau
#'
#' @export
#' @importFrom HTSCluster highDimensionARI
#' @importFrom corrplot corrplot
#' @example inst/examples/coseq-package.R
setMethod("compareARI", "coseqResults", function(object, K=NULL,
                                                 parallel=FALSE, BPPARAM=bpparam(),
                                                 plot=TRUE, ...) {
  compareARI.coseqResults(object, K=K, parallel=parallel, BPPARAM=BPPARAM, plot=plot, ...)
})


#' @rdname compareARI
#' @export
setMethod("compareARI", "matrix", function(object,
                                          parallel=FALSE, BPPARAM=bpparam(),
                                          plot=TRUE, ...) {
  compareARI.matrix(object, parallel=parallel, BPPARAM=BPPARAM,  plot=plot, ...)
})


#' @rdname compareARI
#' @export
setMethod("compareARI", "data.frame", function(object,
                                               parallel=FALSE, BPPARAM=bpparam(),
                                               plot=TRUE, ...) {
  compareARI.matrix(object, parallel=parallel, BPPARAM=BPPARAM,  plot=plot, ...)
})






#######################
## Unexported functions
#######################

compareARI.labels <- function(labels, parallel=FALSE, BPPARAM=bpparam(), plot=TRUE, ...) {
  arg.user <- list(...)
  if(is.null(arg.user$digits)) arg.user$digits<-2
  full_labels <- labels

  ARI <- matrix(0, nrow=ncol(full_labels), ncol=ncol(full_labels))
  rownames(ARI) <- colnames(ARI) <- colnames(full_labels)
  index <- ARI
  index[upper.tri(index)] <- seq(1, (ncol(ARI) * nrow(ARI) - ncol(ARI))/2)

  if(!parallel) {
    tmp <- lapply(seq_len(max(index, na.rm=TRUE)), function(ii) {
      index2 <- which(index == ii, arr.ind=TRUE)
      ARItmp <- highDimensionARI(full_labels[,index2[1]], full_labels[,index2[2]])
      return(list(ARItmp=ARItmp, ii=ii))
    })
  } else if(parallel) {
    tmp <- bplapply(seq_len(max(index, na.rm=TRUE)), function(ii) {
      index2 <- which(index == ii, arr.ind=TRUE)
      ARItmp <- highDimensionARI(full_labels[,index2[1]], full_labels[,index2[2]])
      return(list(ARItmp=ARItmp, ii=ii))
    }, BPPARAM=BPPARAM)
  }

  for(i in seq_len(length(tmp))) {
    new_index <- which(index == tmp[[i]]$ii, arr.ind=TRUE)
    ARI[new_index] <- tmp[[i]]$ARItmp
  }

  diag(ARI) <- 1

  if(plot) {
    corrplot(ARI, is.corr=FALSE, method="color", type="upper", p.mat=ARI, insig="p-value",
             sig.level=-1, tl.pos="d", addgrid.col="white",
             tl.col="black", ...)
  }

  ARI <- round(ARI, digits = arg.user$digits)
  ARI[lower.tri(ARI)] <- ""
  ARI <- data.frame(ARI, check.names=FALSE)
  return(ARI)
}


compareARI.coseqResults <- function(object, K=NULL, parallel=FALSE, BPPARAM=bpparam(), plot=TRUE, ...) {
  x <- object
  arg.user <- list(...)
  if(is.null(arg.user$digits)) arg.user$digits<-2
  ## For class coseqResults
  full_labels <- do.call("cbind", lapply(coseqFullResults(x), clusters))
  if(!is.null(K)) {
    if(!length(K)) {
      stop("K must be a vector of length at least 2.")
    }
    index <- which(substr(colnames(full_labels), 3, 10) %in% K)
    if(!length(index)) stop("None of the indicated models are included in argument x")
    full_labels <- full_labels[,index]
  }
  compareARI.labels(labels=full_labels, parallel=parallel, BPPARAM=BPPARAM, plot=plot, ...)
}


compareARI.matrix <- function(object, parallel=FALSE, BPPARAM=bpparam(), plot=TRUE, ...) {
  x <- object
  arg.user <- list(...)
  if(is.null(arg.user$digits)) arg.user$digits<-2

  ## For class data.frame or matrix
  full_labels <- x
  if(!length(colnames(full_labels))) {
    colnames(full_labels) <- paste("Model", seq_len(ncol(full_labels)))
  }
  compareARI.labels(labels=full_labels, parallel=parallel, BPPARAM=BPPARAM, plot=plot, ...)
}