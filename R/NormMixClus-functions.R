#' Normal mixture model estimation and selection for a series of cluster numbers
#'
#' Perform co-expression and co-abudance analysis of high-throughput
#' sequencing data, with or without data transformation, using a Normal
#' mixture models. The output of \code{NormMixClus} is an S4 object of
#' class \code{RangedSummarizedExperiment}.
#'
#' @param y_profiles (\emph{n} x \emph{q}) matrix of observed profiles for \emph{n}
#' observations and \emph{q} variables
#' @param K Number of clusters (a single value or a sequence of values).
#' @param subset Optional vector providing the indices of a subset of
#' genes that should be used for the co-expression analysis (i.e., row indices
#' of the data matrix \code{y}.
#' @param parallel If \code{FALSE}, no parallelization. If \code{TRUE}, parallel
#' execution using BiocParallel (see next argument \code{BPPARAM}). A note on running
#' in parallel using BiocParallel: it may be advantageous to remove large, unneeded objects
#' from the current R environment before calling the function, as it is possible that R's
#' internal garbage collection will copy these files while running on worker nodes.
#' @param BPPARAM Optional parameter object passed internally to \code{bplapply} when
#' \code{parallel=TRUE}. If not specified, the parameters last registered with \code{register}
#' will be used.
#' @param ... Additional optional parameters to be passed to \code{\link{NormMixClusK}}.
#'
#' @return
#' An S4 object of class \code{coseqResults}, with conditional
#' probabilities of cluster membership for each gene in each model stored as a list of assay
#' data, and corresponding {log likelihood, ICL value, number of
#' clusters, and form of Gaussian model} for each model stored as metadata.
#'
#' @author Andrea Rau, Cathy Maugis-Rabusseau
#'
#' @example inst/examples/NormMixClus.R
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel bpparam
#' @export

NormMixClus <- function(y_profiles, K, subset=NULL, parallel=TRUE, BPPARAM=bpparam(), ...){

  subset.index <- subset

  ## Parse ellipsis function
  providedArgs <- list(...)
  arg.user <- list(alg.type="EM", init.runs=50, init.type="small-em", init.iter=20,
                       iter=1000, cutoff=0.001, GaussianModel="Gaussian_pk_Lk_Ck",
                       verbose=TRUE, digits=3)
  arg.user[names(providedArgs)] <- providedArgs

  y_profiles <- as.data.frame(y_profiles)
  ## In case only a subset of data are to be used for analysis
  if(!is.null(subset.index)) {
    y_profiles <- y_profiles[subset.index,]
  }

  all.results <- vector("list", length = length(K))
  names(all.results) <- paste0("K=", K)
  if(arg.user$verbose) {
    cat("Running K =", min(K), "...\n")
  }
  all.results[[1]] <- suppressWarnings(NormMixClusK(y_profiles, K=min(K), alg.type=arg.user$alg.type,
                                                     init.runs=arg.user$init.runs,
                                                     init.type=arg.user$init.type,
                                                     init.iter=arg.user$init.iter,
                                                     iter=arg.user$iter,
                                                     cutoff=arg.user$cutoff,
                                                     GaussianModel=arg.user$GaussianModel,
                                                     digits=arg.user$digits))
  index <- 2
  remainingK <- K[-which.min(K)]
  if (length(remainingK)) {
    ## In the case where parallelization is NOT used
    if(!parallel) {
      for (k in remainingK) {
        if(arg.user$verbose) {
          cat("Running K =", k, "...\n")
        }
        all.results[[index]] <- suppressWarnings(NormMixClusK(y_profiles=y_profiles, K=k,
                                                               alg.type=arg.user$alg.type,
                                                               init.runs=arg.user$init.runs,
                                                               init.type=arg.user$init.type,
                                                               init.iter=arg.user$init.iter,
                                                               iter=arg.user$iter,
                                                               cutoff=arg.user$cutoff,
                                                               GaussianModel=arg.user$GaussianModel,
                                                               digits=arg.user$digits))
        index <- index + 1
      }
      ## In the case where parallelization IS used
    } else if(parallel) {
      tmp <- bplapply(remainingK, function(ii) {
        if(arg.user$verbose == TRUE) {
          cat("Running K =", ii, "...\n")
        }
        res <- suppressWarnings(NormMixClusK(y_profiles=y_profiles, K=as.numeric(ii),
                                              alg.type=arg.user$alg.type, init.runs=arg.user$init.runs,
                                              init.type=arg.user$init.type, init.iter=arg.user$init.iter,
                                              iter=arg.user$iter, cutoff=arg.user$cutoff,
                                              GaussianModel=arg.user$GaussianModel, digits=arg.user$digits))

        return(res)}, BPPARAM=BPPARAM)
      if(!sum(unlist(lapply(tmp, nrow)))) {
        stop(paste("All models of form", arg.user$GaussianModel, "resulted in estimation errors.
This is likely due to singular covariance matrices when this form of Gaussian mixture is used.
Try rerunning coseq with either a spherical (pk_Lk_Bk) or diagonal (pk_Lk_I)
covariance matrix form instead:

coseq(..., GaussianModel = \"Gaussian_pk_Lk_Bk\")
coseq(..., GaussianModel = \"Gaussian_pk_Lk_I\")"))
      }
      Kmods <- paste0("K=", unlist(lapply(tmp, function(x) nbCluster(x))))
      all.results[match(Kmods, names(all.results))] <- tmp[which(unlist(lapply(tmp, function(xx)
        length(nbCluster(xx)) != 0)) == TRUE)]
    }
  }

  if(!sum(unlist(lapply(all.results, nrow)))) {
    stop(paste("All models of form", arg.user$GaussianModel, "resulted in estimation errors.
This is likely due to singular covariance matrices when this form of Gaussian mixture is used.
Try rerunning coseq with either a spherical (pk_Lk_Bk) or diagonal (pk_Lk_I)
covariance matrix form instead:

coseq(..., GaussianModel = \"Gaussian_pk_Lk_Bk\")
coseq(..., GaussianModel = \"Gaussian_pk_Lk_I\")"))
  }

  nbClust.all <- unlist(lapply(all.results, function(x) nbCluster(x)))
  logLike.all <- unlist(lapply(all.results, function(x) likelihood(x)))
  ICL.all <- unlist(lapply(all.results, function(x) ICL(x)))
  ICL.choose <- names(ICL.all)[which.min(ICL.all)]
  select.results <- all.results[[ICL.choose]]

  pp <- lapply(all.results, function(x) assay(x))
  names(pp) <- names(all.results)

  select.results <- SummarizedExperiment(assay(select.results),
                                         metadata=list(nbCluster=nbClust.all,
                                                       logLike=logLike.all,
                                                       ICL=ICL.all,
                                                       nbClusterError=K[!K %in% nbClust.all],
                                                       GaussianModel=arg.user$GaussianModel))
  RES <- coseqResults(as(select.results, "RangedSummarizedExperiment"), allResults=pp)
  return(RES)
}


#' Normal mixture model estimation
#'
#' Perform co-expression and co-abudance analysis of high-throughput
#' sequencing data, with or without data transformation, using a Normal
#' mixture models for single number of clusters \emph{K}.
#' The output of \code{NormMixClusK} is an S4 object of
#' class \code{RangedSummarizedExperiment}.
#'
#' @param y_profiles y (\emph{n} x \emph{q}) matrix of observed profiles for \emph{n}
#' observations and \emph{q} variables
#' @param K Number of clusters (a single value).
#' @param init.type Type of initialization strategy to be used:
#' \dQuote{\code{small-em}} for the Small-EM strategy, \dQuote{\code{random}}, \dQuote{\code{CEM}},
#' or \dQuote{\code{SEMMax}}
#' @param GaussianModel One of the 28 forms of Gaussian models defined in Rmixmod,
#' by default equal to the \code{"Gaussian_pk_Lk_Ck"} (i.e., a general family model with free
#' proportions, free volume, free shape, and free orientation)
#' @param init.runs Number of runs to be used for the Small-EM strategy, with a default value of 50
#' @param init.iter Number of iterations to be used within each run for the
#' Small-EM strategry, with a default value of 20
#' @param alg.type Algorithm to be used for parameter estimation:
#' \dQuote{\code{EM}}, \dQuote{\code{CEM}}, \dQuote{\code{SEM}}
#' @param cutoff Cutoff to declare algorithm convergence
#' @param iter Maximum number of iterations to be run for the chosen algorithm
#' @param verbose If \code{TRUE}, verbose output is created
#' @param digits Integer indicating the number of decimal places to be used for the
#' \code{probaPost} output
#'
#' @return
#' An S4 object of class \code{RangedSummarizedExperiment}, with conditional
#' probabilities of cluster membership for each gene stored as assay data, and
#' {log likelihood, ICL value, number of
#' clusters, and form of Gaussian model} stored as metadata.
#'
#' @author Cathy Maugis-Rabusseau, Andrea Rau
#'
#' @importFrom Rmixmod mixmodGaussianModel
#' @importFrom Rmixmod mixmodStrategy
#' @importFrom Rmixmod mixmodCluster
#' @example inst/examples/NormMixClus.R
#'
#' @export

NormMixClusK <- function(y_profiles, K, alg.type="EM", init.runs=50,
                          init.type="small-em", GaussianModel="Gaussian_pk_Lk_Ck",
                          init.iter=20, iter=1000, cutoff=0.001, verbose=TRUE, digits=3) {

  if(!is.data.frame(y_profiles)) y_profiles <- as.data.frame(y_profiles)

  models <- mixmodGaussianModel(listModels=c(GaussianModel))
  # strategy
  if(init.type == "small-em") init.type <- "smallEM"
  strats<-mixmodStrategy(algo=alg.type, nbTry = 1,
                         initMethod = init.type, nbTryInInit = init.runs,
                         nbIterationInInit = init.iter, nbIterationInAlgo = iter,
                         epsilonInInit = cutoff, epsilonInAlgo = cutoff,
                         seed = NULL)
  xem <- mixmodCluster(y_profiles, nbCluster=K, models=models, strategy=strats,
                       criterion="ICL")
  pp <- round(proba(xem), digits=digits )
  if(ncol(pp)) {
    colnames(pp) <- paste0("Cluster_", seq_len(ncol(pp)))
    if(!is.null(rownames(y_profiles))) rownames(pp) <- rownames(y_profiles)
    if(is.null(rownames(y_profiles))) rownames(pp) <- seq_len(nrow(pp))
  }

  res <- SummarizedExperiment(assays=pp,
                              metadata=list(logLike=likelihood(xem),
                              ICL=ICL(xem),
                              nbCluster=nbCluster(xem),
                              GaussianModel=GaussianModel))
  return(as(res, "RangedSummarizedExperiment"))
}


#' Calculate the mean and covariance for a Normal mixture model
#'
#' Calculates the mean and covariance parameters for a normal mixture model
#' of the form pK_Lk_Ck
#'
#' @param coseqResults Object of class \code{coseqResults} or \code{RangedSummarizedExperiment}
#' (as output from the \code{NormMixClus} or \code{NormMixClusK} functions)
#' @param y_profiles y (\emph{n} x \emph{q}) matrix of observed profiles for \emph{n}
#' observations and \emph{q} variables, required for \code{x} of class \code{RangedSummarizedExperiment}
#' @param K The model used for parameter estimation for objects \code{x} of
#' class \code{coseq} or \code{NormMixClus}. When \code{NULL}, the model selected
#' by the ICL criterion is used; otherwise, \code{K} should designate the number
#' of clusters in the desired model
#' @param digits Integer indicating the number of decimal places to be used for output
#' @param plot If \code{true}, produce heatmaps to visualize the estimated per-cluster
#' correlation matrices
#' @param ... Additional optional parameters to pass to \code{corrplot}, if desired
#'
#' @return
#' \item{pi }{ Vector of dimension \emph{K} with the estimated cluster proportions from
#' the Gaussian mixture model, where \emph{K} is the number of clusters}
#' \item{mu }{ Matrix of dimension \emph{K} x \emph{d} containing the estimated mean
#' vector from the Gaussian mixture model, where \emph{d} is the
#' number of samples in the data \code{y_profiles} and \emph{K} is the number of clusters}
#' \item{Sigma }{ Array of dimension \emph{d} x \emph{d} x \emph{K} containing the
#' estimated covariance matrices from the Gaussian mixture model, where \emph{d} is the
#' number of samples in the data \code{y_profiles} and \emph{K} is the number of clusters}
#' \item{rho }{ Array of dimension \emph{d} x \emph{d} x \emph{K} containing the
#' estimated correlation matrices from the Gaussian mixture model, where \emph{d} is the
#' number of samples in the data \code{y_profiles} and \emph{K} is the number of clusters}
#'
#' @author Andrea Rau, Cathy Maugis-Rabusseau
#'
#' @example inst/examples/NormMixClus.R
#' @export
#' @importFrom stats cov2cor
#' @importFrom grDevices n2mfrow
#' @importFrom graphics mtext par

NormMixParam <- function(coseqResults, y_profiles=NULL, K=NULL, digits=3, plot=FALSE, ...) {
  x <- coseqResults
  ## TODO: object could be data.frame or matrix
  if (!is(x, "coseqResults") & !is(x, "RangedSummarizedExperiment")) {
    stop(paste0(sQuote("coseqResults")), " must be of class ",
         paste0(dQuote("coseqResults")), " or ", paste0(dQuote("RangedSummarizedExperiment")))
  }
  if(is(x, "coseqResults")) {
    if(nrow(tcounts(x))>0) y_profiles <- as.matrix(as.data.frame(tcounts(x)))
    if(!nrow(tcounts(x)) & is.null(y_profiles)) stop("y_profiles argument needed.")
    if(is.null(K)) {
      probaPost <- assay(x)
    } else {
      if(K == "ICL") probaPost <- assay(x)
      if(K != "ICL") probaPost <- coseqFullResults(x)[[paste0("K=", K)]]
    }
    GaussianModel <- metadata(x)$GaussianModel
  }
  if(is(x, "RangedSummarizedExperiment")) {
    if(is.null(y_profiles)) stop(paste0(dQuote("y_profiles")), " must
                                         be included for class ",
                                         paste0(dQuote("RangedSummarizedExperiment")))
    probaPost <- assay(x)
    GaussianModel <- metadata(x)$GaussianModel
  }

  if(GaussianModel != "Gaussian_pk_Lk_Ck")
    stop("Recalculation of Gaussian parameters currently only supported for",
         paste(dQuote("Gaussian_pk_Lk_Ck")))

  pi <- apply(probaPost,2,sum) / nrow(y_profiles)
  mu <- matrix(0, nrow=ncol(probaPost), ncol=ncol(y_profiles))
  Sigma <- array(0, dim=c(ncol(y_profiles), ncol(y_profiles), ncol(probaPost)))
  for (k in seq_len(ncol(probaPost))) {
    mu[k,]<-apply(probaPost[,k] * y_profiles,2,sum) / sum(probaPost[,k])
    Sigma[,,k] <- (t(y_profiles) - mu[k,]) %*% ( t(t(y_profiles) - mu[k,]) * probaPost[,k])
    Sigma[,,k] <- Sigma[,,k] / sum(probaPost[,k])
  }
  rho <- lapply(seq_len(ncol(probaPost)), function(xx) cov2cor(Sigma[,,xx]))
  rho <- array(unlist(rho), dim = c(dim(rho[[1]]), length(rho)))

  dimnames(rho) <- list(colnames(y_profiles), colnames(y_profiles),
                        paste("Cluster", seq_len(ncol(probaPost))))
  dimnames(Sigma) <- list(colnames(y_profiles), colnames(y_profiles),
                          paste("Cluster", seq_len(ncol(probaPost))))
  colnames(mu) <- colnames(y_profiles)
  rownames(mu) <- paste("Cluster", seq_len(ncol(probaPost)))

  if(plot) {
    par(mfrow=n2mfrow(ncol(probaPost)))
    for(kk in seq_len(ncol(probaPost))) {
      corrplot(rho[,,kk],  method="ellipse", type="upper",
               tl.pos="d", ...)
      mtext(paste("K =", kk), side=2, las=1, line=-1)
    }
  }

  param <- list(pi=round(pi, digits), mu=round(mu, digits), Sigma=round(Sigma, digits),
                rho=round(rho, digits))
  return(param)
}