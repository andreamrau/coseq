#' @rdname coseqResults
#' @export
#' @import SummarizedExperiment
setClass("coseqResults",
         contains = "RangedSummarizedExperiment",
         representation = representation(
           allResults="list",
           model="character",
           transformation="character",
           tcounts="DataFrame",
           y_profiles="DataFrame",
           normFactors="numeric"
         )
         )

## TODO:
## setValidity( ... )


#' coseqResults object and constructor
#'
#' \code{coseqResults} is a subclass of \code{RangedSummarizedExperiment},
#' used to store the co-expression results as well as some additional
#' information useful for plotting (\code{tcounts}, \code{y_profiles}) and
#' meta-information about the co-expression analysis (\code{transformation},
#' \code{normFactors}).
#'
#' This constructor function would not typically be used by "end users".
#' This simple class extends the \code{RangedSummarizedExperiment} class of the
#' SummarizedExperiment package
#' to allow other packages to write methods for results
#' objects from the coseq package. It is used by \code{\link{coseqRun}}
#' to wrap up the results table.
#'
#' @param SummarizedExperiment a \code{RangedSummarizedExperiment} of \code{coseq} results
#' @param allResults List of conditional probabilities of cluster membership for each gene,
#' in all models fit
#' @param model \code{"Normal"} or \code{"Poisson"}, the mixture model used for co-expression
#' @param transformation Transformation applied to counts to obtain \code{tcounts}
#' @param tcounts Transformed counts used for mixture model fitting
#' @param y_profiles y profiles used for \code{coseq} plotting
#' @param normFactors Scaling factors used for normalization
#'
#'
#' @return a coseqResults object
#' @docType class
#' @rdname coseqResults
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom S4Vectors DataFrame
#' @export
coseqResults <- function(SummarizedExperiment,
                         allResults, model=NULL,
                         transformation=NULL,
                         tcounts=NULL,
                         y_profiles=NULL,
                         normFactors=NULL) {
  se <- SummarizedExperiment
  if (!is(se, "RangedSummarizedExperiment")) {
      stop("'SummarizedExperiment' must be a RangedSummarizedExperiment object")
  }
  if(is.null(model)) model <- ""
  if(is.null(transformation)) transformation <- ""
  if(is.null(tcounts)) tcounts <- DataFrame(matrix(0, nrow=0, ncol=0))
  if(is.null(y_profiles)) y_profiles <- DataFrame(matrix(0, nrow=0, ncol=0))
  if(is.null(normFactors)) normFactors <- as.numeric()
  object <- new("coseqResults", se, allResults=allResults, model=model, transformation=transformation,
                tcounts=tcounts, y_profiles=y_profiles, normFactors=normFactors)
  return(object)
}