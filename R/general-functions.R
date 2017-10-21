#' Co-expression analysis
#'
#' Function for primary code to perform co-expression analysis, with or without data transformation,
#' using mixture models. The output of \code{coseqRun} is an S4 object of class \code{coseqResults}.
#'
#' @param y (\emph{n} x \emph{q}) matrix of observed counts for \emph{n}
#' observations and \emph{q} variables
#' @param K Number of clusters (a single value or a vector of values)
#' @param conds Vector of length \emph{q} defining the condition (treatment
#' group) for each variable (column) in \code{y}
#' @param normFactors The type of estimator to be used to normalize for differences in
#' library size: (\dQuote{\code{TC}} for total count, \dQuote{\code{UQ}} for
#' upper quantile, \dQuote{\code{Med}} for median, \dQuote{\code{DESeq}} for
#' the normalization method in the DESeq package, and \dQuote{\code{TMM}} for
#' the TMM normalization method (Robinson and Oshlack, 2010). Can also be a
#' vector (of length \emph{q}) containing pre-estimated library size estimates
#' for each sample.
#' @param model Type of mixture model to use (\dQuote{\code{Poisson}} or \dQuote{\code{Normal}}), or alternatively
#' \dQuote{\code{kmeans}} for a K-means algorithm
#' @param transformation Transformation type to be used: \dQuote{\code{voom}}, \dQuote{\code{logRPKM}}
#' (if \code{geneLength} is provided by user), \dQuote{\code{arcsin}}, \dQuote{\code{logit}},
#' \dQuote{\code{logMedianRef}}, \dQuote{\code{profile}}, \dQuote{\code{logclr}}, \dQuote{\code{clr}},
#' \dQuote{\code{alr}}, \dQuote{\code{ilr}}, or \dQuote{\code{none}}
#' @param subset Optional vector providing the indices of a subset of
#' genes that should be used for the co-expression analysis (i.e., row indices
#' of the data matrix \code{y}. For the generic function \code{coseq}, the results of a previously
#' run differential analysis may be used to select a subset of genes on which to perform the
#' co-expression analysis. If this is desired, \code{subset.index} can also be an object of class
#' DESeqResults (from the \code{results} function in \code{DESeq2}).
#' @param meanFilterCutoff Value used to filter low mean normalized counts if desired (by default,
#' set to a value of 50)
#' @param modelChoice Criterion used to select the best model. For Gaussian mixture models,
#' \dQuote{\code{ICL}} (integrated completed likelihood criterion) is currently supported. For Poisson
#' mixture models, \dQuote{\code{ICL}}, \dQuote{\code{BIC}} (Bayesian information criterion), and a
#' non-asymptotic criterion calibrated via the slope heuristics  using either the \dQuote{\code{DDSE}}
#' (data-driven slope estimation) or \dQuote{\code{Djump}} (dimension jump) approaches may be used.
#' See the \code{HTSCluster} package documentation for more details about the slope heuristics approaches.
#' @param parallel If \code{FALSE}, no parallelization. If \code{TRUE}, parallel
#' execution using BiocParallel (see next argument \code{BPPARAM}). A note on running
#' in parallel using BiocParallel: it may be advantageous to remove large, unneeded objects
#' from the current R environment before calling the function, as it is possible that R's
#' internal garbage collection will copy these files while running on worker nodes.
#' @param BPPARAM Optional parameter object passed internally to \code{bplapply} when
#' \code{parallel=TRUE}. If not specified, the parameters last registered with \code{register}
#' will be used.
#' @param ... Additional optional parameters.
#'
#' @return
#' An S4 object of class \code{coseqResults} whose assays contain a \code{SimpleList}
#' object, where each element in the list corresponds to the conditional probabilities of cluster membership
#' for each gene in each model. Meta data (accessible via \code{metatdata} include the \code{model} used
#' (either \code{Normal} or \code{Poisson}), the \code{transformation} used on the data, the
#' transformed data using to estimate model (\code{tcounts}), the normalized profiles for use in plotting
#' (\code{y_profiles}), and the normalization factors used in the analysis (\code{normFactors}).
#'
#' @author Andrea Rau
#'
#' @export
#' @importFrom HTSCluster PoisMixClus
#' @importFrom HTSCluster PoisMixClusWrapper
#' @importFrom stats na.omit kmeans
#' @importFrom capushe capushe
#' @importClassesFrom S4Vectors DataFrame
#' @importMethodsFrom S4Vectors metadata
#'
#' @examples
#' ## Simulate toy data, n = 300 observations
#' set.seed(12345)
#' countmat <- matrix(runif(300*4, min=0, max=500), nrow=300, ncol=4)
#' countmat <- countmat[which(rowSums(countmat) > 0),]
#' conds <- rep(c("A","B","C","D"), each=2)
#'
#' ## Run the K-means for K = 2,3,4 with logCLR transformation
#' ## The following are equivalent:
#' run <- coseqRun(y=countmat, K=2:15)
#' run <- coseq(object=countmat, K=2:15, transformation="logclr", model="kmeans")
#'
#' ## Run the Normal mixture model for K = 2,3,4 with arcsine transformation
#' ## The following are equivalent:
#' run <- coseqRun(y=countmat, K=2:4, iter=5, transformation="arcsin", model="Normal")
#' run <- coseq(object=countmat, K=2:4, iter=5, transformation="arcsin", model="Normal")
#'
coseqRun <- function(y, K, conds=NULL, normFactors="TMM", model="kmeans", transformation="logclr",
                      subset=NULL, meanFilterCutoff=50,
                      modelChoice=ifelse(model=="kmeans", "DDSE", "ICL"),
                      parallel=FALSE, BPPARAM=bpparam(), ...) {

  if(!is.null(subset)) y <- y[subset,]
  if(!all.equal(round(abs(K)), K)) stop("K should be a vector of cluster numbers.")
  y_profiles <- transformRNAseq(y=y, normFactors=normFactors, transformation="profile",
                                meanFilterCutoff=meanFilterCutoff, verbose=FALSE)$tcounts

  ## Parse ellipsis function
  providedArgs <- list(...)
  arg.user <- list(alg.type="EM", init.runs=50, init.type="small-em", init.iter=20,
                   iter=1000, cutoff=0.001, GaussianModel="Gaussian_pk_Lk_Ck",
                   verbose=TRUE, digits=3,
                   Kmin.init="small-em", split.init=FALSE, fixed.lambda=NA,
                   equal.proportions=FALSE, EM.verbose=FALSE, interpretation="sum",
                   geneLength=NA, iter.max=50, nstart=50, algorithm="MacQueen", trace=FALSE)

  if(model == "Poisson") {
    arg.user$init.runs <- 1
    arg.user$init.iter <- 10
    arg.user$cutoff <- 1e-05
    if(is.null(subset)) arg.user$subset.index <- NA
  }
  arg.user[names(providedArgs)] <- providedArgs

  cat("****************************************\n")
  cat("coseq analysis:", model, "approach &", transformation, "transformation\n")
  cat("K =", min(K), "to", max(K), "\n")
  cat("****************************************\n")

  ########################
  ## POISSON MIXTURE MODEL
  ########################
  if(length(model) & model == "Poisson") {
    if(transformation != "none") stop("Poisson mixture model may only be applied on raw counts.")
    if(is.null(conds)) {
      message("Poisson mixture model fit assuming each sample is an independent condition.")
      conds <- seq_len(ncol(y))
    }

    ## Grouping columns of y in order of condition (all replicates put together)
    o.ycols <- order(conds)
    y <- y[,o.ycols]
    conds <- conds[o.ycols]
    conds.names <- unique(conds)
    d <- length(unique(conds))
    r <- as.vector(table(conds))
    if(!length(rownames(y))) rn <- 1:nrow(y)
    if(length(rownames(y))) rn <- rownames(y)
    y <- as.matrix(y, nrow = nrow(y), ncol = ncol(y))
    rownames(y) <- rn


    n <- dim(y)[1]
    cols <- dim(y)[2]
    w <- rowSums(y)


    tcounts <- transformRNAseq(y=y, normFactors=normFactors, transformation="none",
                                geneLength=arg.user$geneLength,
                                meanFilterCutoff=meanFilterCutoff, verbose=FALSE)
    conds <- as.vector(conds)

    if(!parallel) {
      run <- suppressWarnings(PoisMixClusWrapper(y=tcounts$tcounts, gmin=min(K),
                                                 gmax=max(K), conds=conds,
                                                 norm=tcounts$ellnorm / sum(tcounts$ellnorm),
                                                 gmin.init.type=arg.user$Kmin.init,
                                                 split.init=arg.user$split.init, subset.index=NA,
                                                 init.runs=arg.user$init.runs, init.iter=arg.user$init.iter,
                                                 alg.type=arg.user$alg.type, cutoff=arg.user$cutoff,
                                                 iter=arg.user$iter,
                                                 fixed.lambda=arg.user$fixed.lambda,
                                                 equal.proportions=arg.user$equal.proportions,
                                                 verbose=arg.user$verbose, EM.verbose=arg.user$EM.verbose,
                                                 interpretation=arg.user$interpretation))
      names(run$all.results) <- paste0("K=", K)
      run$nbCluster.all <- K
    }


    if(parallel) {
      all.results <- vector("list", length = length(K))
      names(all.results) <- paste0("K=", K)
      cat("Running K =", min(K), "...\n")
      if(arg.user$split.init == TRUE) {
        warning("Splitting initialization is not compatible with parallelization.")
      }
      run <- PoisMixClus(y=tcounts$tcounts, g=min(K), conds=conds,
                         norm=tcounts$ellnorm / sum(tcounts$ellnorm),
                         init.type=arg.user$Kmin.init,
                         subset.index=NA,
                         wrapper=TRUE, init.runs=arg.user$init.runs,
                         init.iter=arg.user$init.iter, alg.type=arg.user$alg.type,
                         cutoff=arg.user$cutoff, iter=arg.user$iter,
                         fixed.lambda=arg.user$fixed.lambda,
                         equal.proportions=arg.user$equal.proportions,
                         prev.labels=NA, prev.probaPost = NA,
                         verbose=arg.user$verbose,
                         EM.verbose=arg.user$EM.verbose,
                         interpretation=arg.user$interpretation)
      all.results[[1]] <- run

      index <- 2
      remainingK <- K[-which(K == min(K))]
      if(length(remainingK) > 0) {
        tmp <- bplapply(remainingK, function(ii, P_y, P_conds, P_norm,
                                             P_init.type,
                                             P_init.runs, P_init.iter,
                                             P_alg.type, P_cutoff,
                                             P_iter, P_fixed.lambda,
                                             P_equal.proportions, P_verbose,
                                             P_interpretation, P_EM.verbose) {
          cat("Running K =", ii, "...\n")
          res <- PoisMixClus(g=as.numeric(ii),
                             y=P_y,
                             conds=P_conds,
                             norm=P_norm,
                             init.type=P_init.type,
                             subset.index=NA,
                             wrapper=TRUE,
                             prev.probaPost=NA,
                             prev.labels=NA,
                             init.runs = P_init.runs,
                             init.iter = P_init.iter,
                             alg.type = P_alg.type,
                             cutoff = P_cutoff,
                             iter = P_iter,
                             fixed.lambda = P_fixed.lambda,
                             equal.proportions = P_equal.proportions,
                             verbose = P_verbose,
                             interpretation = P_interpretation,
                             EM.verbose = P_EM.verbose)
          return(res)},
          P_y=tcounts$tcounts, P_conds=conds, P_norm=tcounts$ellnorm / sum(tcounts$ellnorm),
          P_init.type=arg.user$Kmin.init,
          P_init.runs=arg.user$init.runs, P_init.iter=arg.user$init.iter,
          P_alg.type=arg.user$alg.type, P_cutoff=arg.user$cutoff,
          P_iter=arg.user$iter, P_fixed.lambda=arg.user$fixed.lambda,
          P_equal.proportions=arg.user$equal.proportions, P_verbose=arg.user$verbose,
          P_interpretation=arg.user$interpretation, P_EM.verbose=arg.user$EM.verbose,
          BPPARAM=BPPARAM)
        Kmods <- paste0("K=", unlist(lapply(tmp, function(x) ncol(x$lambda))))
        all.results[-1] <- tmp[na.omit(match(names(all.results), Kmods))]
      }

      logLike.all <- unlist(lapply(all.results, function(x) x$log.like))
      ICL.all <- unlist(lapply(all.results, function(x) x$ICL))
      ICL.choose <- which.max(ICL.all)
      select.results <- all.results[[ICL.choose]]
      select.results$model.selection <- "ICL"

      BIC.all <- unlist(lapply(all.results, function(x) x$BIC))
      BIC.choose <- which.max(BIC.all)
      select.results2 <- all.results[[BIC.choose]]
      select.results2$model.selection <- "BIC"

      # Apply capushe: only if at least 10 models are considered
      if(c(max(K) - min(K) + 1) <= 10) {
        message("Note: slope heuristics for model selection only applied if > 10 models are fit.")
        DDSE.results <- NA
        Djump.results <- NA
        capushe <- NA
        ResCapushe <- NA
      }
      if(c(max(K) - min(K) + 1) > 10) {
        message("Note: diagnostic plots for slope heuristics (Djump and DDSE) should be examined to \n
                ensure that sufficiently complex models have been considered.")
        Kchoice <- K
        np <- (Kchoice-1) + (length(unique(conds))-1)*(Kchoice)
        mat <- cbind(Kchoice, np/n, np/n, -logLike.all/n)
        ResCapushe <- suppressWarnings(capushe(mat, n))
        DDSE <- DDSEextract(ResCapushe)
        Djump <- Djumpextract(ResCapushe)
        DDSE.results <- all.results[[paste0("K=", DDSE)]]
        Djump.results <- all.results[[paste0("K=", Djump)]]
        DDSE.results$model.selection <- "DDSE"
        Djump.results$model.selection <- "Djump"
      }
      run <- list(nbCluster.all=K, logLike.all = logLike.all,
                  ICL.all = ICL.all,
                  capushe = ResCapushe,
                  all.results = all.results,
                  DDSE.results = DDSE.results,
                  Djump.results = Djump.results,
                  BIC.results = select.results2,
                  ICL.results = select.results)
    }


    if(modelChoice == "BIC") final.results <- run$BIC.results
    if(modelChoice == "ICL") final.results <- run$ICL.results
    if(modelChoice == "DDSE") final.results <- run$DDSE.results
    if(modelChoice == "Djump") final.results <- run$Djump.results

    final.results$probaPost <- round(final.results$probaPost, arg.user$digits)
    colnames(final.results$probaPost) <- paste0("Cluster_", seq_len(ncol(final.results$probaPost)))
    rownames(final.results$probaPost) <- rownames(tcounts$tcounts)
    for(jj in seq_len(length(run$all.results))) {
      run$all.results[[jj]]$probaPost <- round(run$all.results[[jj]]$probaPost,
                                               arg.user$digits)
    }

    ## CONVERT RESULTS TO SUMMARIZEDEXPERIMENT0 CLASS
    nbClust.all <- run$nbCluster.all
    logLike.all <- run$logLike.all
    ICL.all <- run$ICL.all

    pp <- lapply(run$all.results, function(x) {
      tmp <- round(x$probaPost, digits=arg.user$digits)
      colnames(tmp) <- paste0("Cluster_", seq_len(ncol(tmp)))
      rownames(tmp) <- rownames(tcounts$tcounts)
      return(tmp)
    })
    names(pp) <- names(run$all.results)


    ## Format HTSCluster results
    ICL.results <- SummarizedExperiment(final.results$probaPost, metadata=list(nbCluster=nbClust.all,
                                                                               logLike=logLike.all,ICL=ICL.all))

    tcountsDF = as(tcounts$tcounts, "DataFrame")
    y_profilesDF <- y_profiles <- as(tcounts$tcounts/rowSums(tcounts$tcounts), "DataFrame")
    run <- coseqResults(as(ICL.results, "RangedSummarizedExperiment"),
                        allResults = pp,
                        model="Poisson",
                        transformation="none",
                        tcounts=tcountsDF,
                        y_profiles=y_profilesDF,
                        normFactors=tcounts$ellnorm / sum(tcounts$ellnorm))
  }



  ########################
  ## NORMAL MIXTURE MODEL
  ########################
  if(length(model) & model == "Normal") {

    if(modelChoice != "ICL") message("Note: only ICL is currently supported for model choice for Normal mixture models.")
    tcounts <- transformRNAseq(y=y, normFactors=normFactors, transformation=transformation,
                                geneLength=arg.user$geneLength, meanFilterCutoff=meanFilterCutoff, verbose=FALSE)
    run <- NormMixClus(y_profiles=tcounts$tcounts, K=K, subset=NULL,
                       parallel=parallel,
                       BPPARAM=BPPARAM, alg.type=arg.user$alg.type,
                       init.runs=arg.user$init.runs,
                       init.type=arg.user$init.type, init.iter=arg.user$init.iter,
                       iter=arg.user$iter, cutoff=arg.user$cutoff,
                       verbose=arg.user$verbose, digits=arg.user$digits,
                       GaussianModel=arg.user$GaussianModel)
  }


  ########################
  ## K-means
  ########################

  if(length(model) & model == "kmeans") {

    if(modelChoice != "DDSE") message("Note: only DDSE is currently supported for model choice for K-means.")
    if(!transformation %in% c("clr", "alr", "logclr", "ilr")) {
      message("Transformation used is: ", transformation, "\n
              Typically one of the following profile transformations is used with K-means: clr, alr, ilr, logclr")
    }
    tcounts <- transformRNAseq(y=y, normFactors=normFactors, transformation=transformation,
                               geneLength=arg.user$geneLength,
                               meanFilterCutoff=meanFilterCutoff, verbose=FALSE)

    n <- nrow(tcounts$tcounts)
    d <- ncol(tcounts$tcounts)
    km_cluster <- vector("list", length(K)) ## Initialisation pour les vecteurs de classification
    tot_withinss <- c() ## Initialisation pour l'inertie intra
    names(km_cluster) <- paste0("K=", K)

    Kerr <- c()
    if(!parallel) {
      for(k in K) {
        if(arg.user$verbose) cat("Running K =", k, "...\n")
        km <- suppressWarnings(kmeans(tcounts$tcounts, centers=k,
                                      iter.max=arg.user$iter.max,
                                      nstart=arg.user$nstart,
                                      algorithm=arg.user$algorithm,
                                      trace=arg.user$trace))
        km_cluster[[paste0("K=",k)]] <- km$cluster
        tot_withinss <- c(tot_withinss, km$tot.withinss)
        if(!is.null(km$ifault)) Kerr <- c(Kerr, k)
      }
      names(tot_withinss) <- paste0("K=", K)
    }
    if(parallel) {
      tot_withinss <- rep(NA, length(K)) ## Initialisation pour l'inertie intra
      names(tot_withinss) <- paste0("K=", K)
      tmp <- bplapply(K, function(ii, km_tcounts, km_iter.max, km_nstart, km_algorithm,
                                  km_trace, km_verbose) {
        if(km_verbose) cat("Running K =", ii, "...\n")
        km <- suppressWarnings(kmeans(km_tcounts, centers=as.numeric(ii),
                                      iter.max=km_iter.max,
                                      nstart=km_nstart,
                                      algorithm=km_algorithm,
                                      trace=km_trace))
        return(km)}, km_tcounts=tcounts$tcounts, km_iter.max=arg.user$iter.max,
        km_nstart=arg.user$nstart, km_algorithm=arg.user$algorithm, km_trace=arg.user$trace,
        km_verbose=arg.user$verbose, BPPARAM=BPPARAM)
      Kmods <- paste0("K=", unlist(lapply(tmp, function(x) length(x$size))))
      km_cluster <- lapply(tmp[na.omit(match(names(km_cluster), Kmods))], function(x) x$cluster)
      names(km_cluster) <- paste0("K=", K)
      tot_withinss <- unlist(lapply(tmp[na.omit(match(names(km_cluster), Kmods))],
                                    function(x) x$tot.withinss))
      names(tot_withinss) <- Kmods
      Kerr <- unlist(lapply(tmp[na.omit(match(names(km_cluster), Kmods))],
                                    function(x) ifelse(is.null(x$ifault),
                                                       NA, length(x$size))))
      Kerr <- Kerr[which(!is.na(Kerr))]
    }

    if(length(K) < 10)
      warning("Be careful: for model selection via capushe, at least 10 models should be estimated.")
    cap <- suppressWarnings(capushe(matrix(c(K, sqrt(n*d*K), sqrt(n*d*K), tot_withinss), ncol=4)))
    K_select <- paste0("K=",DDSEextract(cap)[1]) # nombre de classes s?l?ctionn? par capushe
    cluster_select <- km_cluster[[K_select]]

#   ## Keep the estimated posterior probabilities for K-means rather than just the cluster labels
#    pp_select <- round(kmeansProbaPost(clusters=cluster_select, tcounts=tcounts$tcounts), arg.user$digits)
    pp_select <- matrix(0, nrow=nrow(tcounts$tcounts), ncol=as.numeric(DDSEextract(cap)[1]))
    pp_select[cbind(seq_len(length(cluster_select)), cluster_select)] <- 1
    colnames(pp_select) <- paste0("Cluster_", seq_len(ncol(pp_select)))
#    rownames(pp_select) <- rownames(tcounts$tcounts)
    rownames(pp_select) <- rownames(y_profiles)


    nbClust.all <- K
    names(nbClust.all) <- names(km_cluster)

    select.results <- SummarizedExperiment(pp_select,
                                           metadata=list(nbCluster=nbClust.all,
                                                         DDSE=paste0("K=", DDSEextract(cap)[1]),
                                                         Djump=paste0("K=", Djumpextract(cap)[1]),
                                                         capushe=cap,
                                                         tot_withinss=tot_withinss,
                                                         nbClusterError=Kerr))

    all.results <- vector("list", length(K))
    names(all.results) <- paste0("K=", K)
    for(k in K) {
      # all.results[[paste0("K=",k)]] <- round(kmeansProbaPost(clusters=km_cluster[[paste0("K=", k)]],
      #                                                 tcounts=tcounts$tcounts), arg.user$digits)
      all.results[[paste0("K=",k)]] <- matrix(0, nrow=nrow(tcounts$tcounts), ncol=k)
      all.results[[paste0("K=",k)]][cbind(seq_len(length(km_cluster[[paste0("K=", k)]])),
                                          km_cluster[[paste0("K=", k)]])] <- 1
      colnames(all.results[[paste0("K=",k)]]) <- paste0("Cluster_", seq_len(ncol(all.results[[paste0("K=",k)]])))
#      rownames(all.results[[paste0("K=",k)]]) <- rownames(tcounts$tcounts)
      rownames(all.results[[paste0("K=",k)]]) <- rownames(y_profiles)
    }

    run <- coseqResults(as(select.results, "RangedSummarizedExperiment"), allResults=all.results)
  }

  ####################################
  ## RETURN RESULTS
  ####################################

  ICL.results <- SummarizedExperiment(assay(run), metadata=metadata(run))
  all.results <- coseqFullResults(run)
  tcountsDF <- as(tcounts$tcounts, "DataFrame")
  y_profilesDF <- as(y_profiles, "DataFrame")

  colnames(tcountsDF) <- colnames(y_profilesDF)
  rownames(tcountsDF) <- rownames(y_profilesDF)

  RESULTS <- coseqResults(as(ICL.results, "RangedSummarizedExperiment"),
                          allResults=all.results,
                          model=model,
                          transformation=transformation,
                          tcounts=tcountsDF,
                          y_profiles=y_profilesDF,
                          normFactors=tcounts$snorm)
  return(RESULTS)
}



#' Calculate conditional probabilities of cluster membership for K-means clustering
#'
#' @param clusters Cluster labels arising from K-means clustering
#' @param tcounts Transformed counts clustered using K-means
#'
#' @return Conditional probabilities of cluster membership for each observation in each cluster
#' @export
#' @importFrom mvtnorm dmvnorm
#'
#' @examples
#' ## Example of K-means taken from ?kmeans help page
#' x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
#'            matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
#'            colnames(x) <- c("x", "y")
#' cl <- kmeans(x, 5)
#' probaPost <- kmeansProbaPost(cl$cluster, x)
#' head(probaPost)
kmeansProbaPost <- function(clusters,tcounts)
{
  results <- clusters
  K <- max(results)
  centers <- matrix(0,nrow=K,ncol=ncol(tcounts))
  dens <- matrix(0,nrow=nrow(tcounts),ncol=K)
  withinss <- c() ###initialisation pour les inerties inter
  size <- c() ### initialisation pour les tailles de classes
  atyp <- c() ### indices des classes a un element
  for(i in seq_len(K)) {
    I <- which(results==i)
    size <- c(size,length(I))
    if(length(I)>1) {
      centers[i,] <- colSums(tcounts[I,])/length(I)
      var <- sum((tcounts[I,] - matrix(rep(colSums(tcounts[I,])/length(I),length(I)),
                                       ncol=ncol(tcounts),byrow=TRUE))^2)/length(I)
      withinss <- c(withinss,var)
      dens[,i] <- 1/length(I)*exp(dmvnorm(tcounts,centers[i,], var/ncol(tcounts)*diag(ncol(tcounts)),
                                          log=TRUE))
    }
    if(length(I)==1) {
      centers[i,] <- tcounts[I,]
      var <- 0
      dens[,i] <- rep(0,nrow(tcounts))
      atyp <- c(atyp,i)
    }
  }
  # dens <- dens/rowSums(dens)
  # epsilon <- 1e-10
  # maxcut <- 1 - epsilon
  # mincut <- epsilon
  # dens <- apply(dens, 2, pmax, mincut)
  # dens <- apply(dens, 2, pmin, maxcut)

  tauik <- dens/matrix(rep(rowSums(dens),K),nrow=nrow(tcounts))
  if(length(atyp)>0) {
    for (i in atyp) {
      I <- which(results==i)
      tauik[I,i] <- 1
    }
  }
  rownames(tauik) <- rownames(tcounts)
  colnames(tauik) <- paste0("Cluster_", seq_len(K))
  return(tauik)
}






#' Compare corrected ICL values after data transformation
#'
#' Compare the corrected ICL values after applying the arcsin, logit, and logMedianRef
#' transformations in a coseq analysis
#'
#' @param x A list made up of \code{coseqResults} objects. At the current time, this function
#' only supports the comparison of \code{coseqResults} objects using \code{model="Normal"} and
#' \code{transformation = c("arcsin", "logit", "logMedianRef")}
#'
#' @return A plot of corrected ICL values for the models included in \code{x} (the list
#' of \code{coseqResults} objects)
#'
#' @author Andrea Rau, Cathy Maugis-Rabusseau
#'
#' @export
#' @example inst/examples/coseq-package.R
compareICL <- function(x) {

  coseq_list <- x

  if(length(coseq_list) > 3) stop("Function only currently supported for <= 3 coseq objects")
  transf <- unlist(lapply(coseq_list, function(xx) transformationType(xx)))
  check <- which(transf != "arcsin" & transf != "logit" & transf != "logMedianRef")
  if(length(check))
    stop("Function only currently supported for arcsin, logit, and logMedianRef transformations")
  ## TODO: CHECK THAT EACH TRANSFORMATION IS ONLY PRESENT ONCE
  resarcsin <- reslogit <- reslogMedianRef <- NULL
  if(length(which(transf == "arcsin"))) {
    resarcsin <- coseq_list[[which(transf == "arcsin")]]
    resarcsin_profiles <- as.matrix(as.data.frame(profiles(resarcsin)))
  }
  if(length(which(transf == "logit"))) {
    reslogit <- coseq_list[[which(transf == "logit")]]
    reslogit_profiles <- as.matrix(as.data.frame(profiles(reslogit)))
  }
  if(length(which(transf == "logMedianRef"))) {
    reslogMedianRef <- coseq_list[[which(transf == "logMedianRef")]]
    reslogMedianRef_profiles <- as.matrix(as.data.frame(profiles(reslogMedianRef)))
  }

  ## TODO: check that same data are used for all transformations
  #if(sum(profiles(resarcsin) != profiles(reslogit)) > 0 |
  #   sum(profiles(resarcsin) != profiles(reslogMedianRef)) > 0)
  #  stop("y_profiles in coseq objects not equal -- are models estimated on same data?")

  ## NOTE: REPLACE 0's with smallest value > 0, 1's with largest value < 1
  PP <- resarcsin_profiles
  PP[which(PP == 0)] <- min(PP[which(PP > 0)])
  PP[which(PP == 1)] <- max(PP[which(PP < 1)])

  n <- dim(PP)[1]
  p <- dim(PP)[2]
  qarcsin <- (n*p*log(2)) + (0.5*sum(sum(log(PP*(1-PP)))))
  qlogit <- (n*p*log(log(2))) + (sum(sum(log(PP*(1-PP)))))
  qlogmedianref <- (n*p*log(log(2))) + sum(sum(log(PP)))

  plotmat <- matrix(NA, nrow=0, ncol=3)
  colnames(plotmat) <- c("K", "ICL", "Transformation")

  if(!is.null(resarcsin)) {
    plotmat <- rbind(plotmat, cbind(metadata(resarcsin)$nbCluster,
                                    metadata(resarcsin)$ICL + (2*qarcsin),
                                    rep("arcsin", length(metadata(resarcsin)$ICL))))
  }
  if(!is.null(reslogit)) {
    plotmat <- rbind(plotmat, cbind(metadata(reslogit)$nbCluster,
                                    metadata(reslogit)$ICL + (2*qlogit),
                                    rep("logit", length(metadata(reslogit)$ICL))))
  }
  if(!is.null(reslogMedianRef)) {
    plotmat <- rbind(plotmat, cbind(metadata(reslogMedianRef)$nbCluster,
                                    metadata(reslogMedianRef)$ICL + (2*qlogmedianref),
                                    rep("logMedianRef",
                                        length(metadata(reslogMedianRef)$ICL))))
  }
  plotdf <- data.frame(plotmat, row.names=seq_len(nrow(plotmat)),
                       stringsAsFactors = FALSE)
  plotdf$Transformation <- factor(plotdf$Transformation)
  plotdf$ICL <-  as.numeric(plotdf$ICL)
  plotdf$K <-  as.numeric(plotdf$K)

  g1 <- ggplot(plotdf, aes_string(x="K", y="ICL")) +
    geom_line(aes_string(color="Transformation")) +
    geom_point(aes_string(color="Transformation")) +
    scale_y_continuous(name="Corrected ICL") + theme_bw()

  print(g1)
}


#' Calculation of per-cluster entropy
#'
#' Provides the calculation of per-cluster entropy, equivalent to
#' \deqn{Entropy(k) = \sum_{i \in C_k} \log (\tau_{ik})}
#' where \eqn{\tau_{ik}} is the conditional probability of gene \emph{i} belonging
#' to cluster \emph{k} and \eqn{C_k} corresponds to the set of indices of genes
#' attributed to cluster \emph{k}.
#'
#' @param probaPost Matrix containing the conditional probabilities of belonging
#' to each cluster for all observations
#'
#' @return Entropy per cluster
#'
#' @author Cathy Maugis-Rabusseau
#'
#' @examples
#' ## Generate artificial matrix of conditional probabilities for K=5 clusters
#' tmp <- matrix(runif(100*5), nrow=100, ncol=5)
#' probaPost <- tmp / rowSums(tmp)
#' clusterEntropy(probaPost)
#'
#' @export
clusterEntropy <- function(probaPost) {
  label <- apply(probaPost,1,which.max)
  entrop <- NULL
  for (k in seq_len(max(label))){
    II <- which(label==k)
    entrop <- c(entrop,sum(log(probaPost[II,k])))
  }
  return(entrop)
}


#' Calculation of within-cluster inertia
#'
#' Provides the calculation of within-cluster inertia, equivalent to
#' \deqn{Inertia(k) = \sum_{i \in C_k} (y_{ik} - \mu_k)^2}
#' where \eqn{\mu_k} is the mean of cluster \emph{k} and \eqn{C_k} corresponds to the set of indices of genes
#' attributed to cluster \emph{k}.
#'
#' @param profiles Matrix, data.frame, or DataFrame containing the (transformed) profiles used for the clustering
#' @param clusters Vector of cluster labels corresponding to the observations in \code{profiles}
#'
#' @return Within cluster inertia
#'
#' @author Andrea Rau, Antoine Godichon-Baggioni
#'
#' @export
#' @examples
#' ## Simulate toy data, n = 300 observations
#' set.seed(12345)
#' countmat <- matrix(runif(300*4, min=0, max=500), nrow=300, ncol=4)
#' countmat <- countmat[which(rowSums(countmat) > 0),]
#' conds <- rep(c("A","B","C","D"), each=2)
#'
#' ## Run the K-means algorithm for logclr profiles for K = 2,..., 20
#' run_kmeans <- coseq(object=countmat, K=2:20, transformation="logclr",
#' model="kmeans")
#' clusterInertia(profiles=tcounts(run_kmeans), clusters=clusters(run_kmeans))
#'
clusterInertia <- function(profiles, clusters) {
  inertia <- NULL
  ## TODO: check that the length of clusters equals profiles
  for(k in seq_len(max(clusters))){
    II <- which(clusters==k)
    if(length(II) > 1) {
      tmp <- as.matrix(profiles[II,])
      inertia <- c(inertia, sum(t(t(tmp) - colMeans(tmp))^2))
    } else inertia <- c(inertia, 0)
  }
  return(inertia)
}


#' Permute columns of a contingency table
#'
#' Permute the columns of a contingency table comparing two clusterings
#' to load the diagonal as much as possible.
#'
#' @param table_1 Partition from a first data clustering
#' @param table_2 Partition from a second data clustering
#'
#' @return Permuted contingency table
#' @export
#'
#' @examples
#' ## Generate arbitrary labels from two separate clustering results
#' labels_1 <- sample(1:10, 1000, replace=TRUE)  ## K=10 clusters
#' labels_2 <- sample(1:8, 1000, replace=TRUE)   ## K=8 clusters
#' matchContTable(labels_1, labels_2)
#'
#' @importFrom e1071 matchClasses
matchContTable <- function(table_1, table_2){
  tab <- table(table_1, table_2)
  ## Put larger clustering in rows if needed, nrow(tab) >= ncol(tab)
  transpose <- FALSE
  if(nrow(tab) < ncol(tab)) transpose <- TRUE
  if(transpose) tab <- t(tab)
  ## Order rows according to largest clusters
  tab <- tab[order(apply(tab,1,max), decreasing=TRUE),]
  ## Match best column with each row of tab
  ## Use unique indices as some columns might map to multiple rows
  index <- matchClasses(tab, method=ifelse(nrow(tab)==ncol(tab), "exact", "rowmax"))
  tabord <- tab[,unique(index)]
  if(transpose) tabord <- t(tabord)
  return(tabord)
}

#' Transform RNA-seq data using common transformations
#'
#' Application of common transformations for RNA-seq data prior to fitting a normal mixture model
#'
#' @param y (\emph{n} x \emph{q}) \code{matrix} or \code{data.frame} of observed counts
#' for \emph{n} observations and \emph{q} variables
#' @param normFactors The type of estimator to be used to normalize for differences in
#' library size: \dQuote{\code{TC}} for total count, \dQuote{\code{DESeq}} for
#' the normalization method in the DESeq package, and \dQuote{\code{TMM}} for
#' the TMM normalization method (Robinson and Oshlack, 2010). Can also be a
#' vector (of length \emph{q}) containing pre-estimated library size estimates
#' for each sample.
#' @param transformation Transformation type to be used: \dQuote{\code{arcsin}},
#' \dQuote{\code{logit}}, \dQuote{\code{logMedianRef}}, \dQuote{\code{profile}},
#' \dQuote{\code{voom}}, \dQuote{\code{logRPKM}} (if \code{geneLength} is provided by user),
#' \dQuote{\code{logclr}}, \dQuote{\code{clr}}, \dQuote{\code{alr}}, \dQuote{\code{ilr}},
#' \dQuote{\code{none}},
#' @param geneLength Vector of length equal to the number of rows in \dQuote{\code{y}} providing
#' the gene length (bp) for RPKM calculation
#' @param meanFilterCutoff Value used to filter low mean normalized counts
#' @param verbose If \code{TRUE}, include verbose output
#'
#' @return
#' \item{tcounts }{Transformed counts}
#' \item{normCounts }{Normalized counts}
#' \item{snorm }{Per-sample normalization factors divided by mean normalization factor}
#' \item{ellnorm }{Per-sample normalization factors}
#'
#' @export
#'
#' @examples
#' set.seed(12345)
#' countmat <- matrix(runif(300*4, min=0, max=500), nrow=300, ncol=4)
#' countmat <- countmat[which(rowSums(countmat) > 0),]
#' conds <- rep(c("A","B","C","D"), each=2)
#'
#' ## Arcsin transformation, TMM normalization
#' arcsin <- transformRNAseq(countmat, normFactors="TMM", transformation="arcsin")$tcounts
#' ## Logit transformation, TMM normalization
#' logit <- transformRNAseq(countmat, normFactors="TMM", transformation="logit")$tcounts
#' ## logCLR transformation, TMM normalization
#' logclr <- transformRNAseq(countmat, normFactors="TMM", transformation="logclr")$tcounts
#'
#' @importFrom edgeR calcNormFactors
#' @importFrom edgeR cpm
#' @importFrom HTSFilter HTSBasicFilter
#' @importFrom DESeq2 varianceStabilizingTransformation
#' @importFrom stats median
#' @importFrom compositions alr clr ilr
transformRNAseq <- function(y, normFactors="TMM", transformation="arcsin",
                             geneLength=NA, meanFilterCutoff=NULL, verbose=TRUE) {

  ##################################
  ## Calculate normalization factors
  ##################################
  ## Only calculate s values if they are not provided
  if(length(normFactors) != 1 & length(normFactors) != ncol(y)) stop(paste(sQuote("normFactors"), "must be one of
     the following: none, TC, DESeq, TMM, or a vector oflength equal to the number of columns
     in", sQuote("y")))
  ## If estimated from data, all genes should be used
  if(length(normFactors) == 1) {
    if(normFactors == "none") {
      libsize <- rep(1, ncol(y))
    }
    if(normFactors == "TMM") {
      f <- calcNormFactors(as.matrix(y), method="TMM")
      libsize <- colSums(y) * f
    }
    if(normFactors == "TC") {
      libsize <- colSums(y)
    }
    if(normFactors == "DESeq") {
      loggeomeans <- rowMeans(log(y))
      f <- apply(y, 2, function(x)
        exp(median((log(x) - loggeomeans)[is.finite(loggeomeans)])))
      libsize <- colSums(y) * f
    }
  }
  if(length(normFactors) > 1) {
    libsize <- normFactors
  }

  snorm <- libsize / (mean(libsize))
  ellnorm <- libsize

  ##################################
  ## Filter data based on mean count if desired
  ##################################
  if(!is.null(meanFilterCutoff)) {

    filter <- HTSBasicFilter(y, method="mean", cutoff.type="value",
                             cutoff=meanFilterCutoff, normalization=normFactors)
    if(verbose) {
      cat("Filter applied: remove observations with normalized mean < ", meanFilterCutoff, "\n")
      cat("               ", nrow(filter$filteredData), "observations retained for analysis\n")
    }
    newy <- filter$filteredData
    colnames(newy) <- colnames(y)
    rownames(newy) <- rownames(y)[which(filter$on==1)]
    y <- newy
    geneLength <- unlist(geneLength)[which(filter$on == 1)]
  }

  ##################################
  ## Transform data and calculate the first derivative of transformations
  ##################################

  normCounts <- t(t(y)/snorm +1)
  if(transformation == "none") {
    tcounts <- y
  }
  if(transformation == "profile") {
    tcounts <- normCounts / rowSums(normCounts)
  }
  if(transformation == "clrProfile") {
    profiles <- normCounts / rowSums(normCounts)
    gm <- apply(profiles, 1, function(x) exp(mean(log(x))))
    tcounts <- log(profiles / gm)
  }
  if(transformation == "voom") {
    tcounts <- log2(t((t(y + 0.5)/(ellnorm + 1)))*1e+06)
  }
  if(transformation == "logRPKM") {
    if(is.na(geneLength)[1]) stop("RPKM transformation requires input for gene lengths")
    if(length(unlist(geneLength))!=nrow(y)) stop("Gene length vector of different dimension than counts")
    Lcounts <- y/geneLength
    tcounts <- log2(t(t(Lcounts)/ellnorm) * 10^9 + 0.5)
    elltmp <- matrix(rep(ellnorm, each = nrow(y)), nrow=nrow(y))
    lentmp <- matrix(rep(geneLength, times=ncol(y)), nrow=nrow(y))
  }
  if(transformation == "arcsin") {
    props <- normCounts / rowSums(normCounts)
    tcounts <- asin(sqrt(props))
    w <- matrix(rep(rowSums(normCounts), times=ncol(y)), nrow=nrow(y))
  }
  if(transformation == "logit") {
    props <- normCounts / rowSums(normCounts)
    tcounts <- log2(props / (1-props))
    w <- matrix(rep(rowSums(normCounts), times=ncol(y)), nrow=nrow(y))
    stemp <- matrix(rep(snorm, each = nrow(y)), nrow=nrow(y))
  }
  if(transformation == "logMedianRef") {
    m <- apply((t(t(y)/snorm))+1, 1, median)
    tcounts <- log2((t(t(y)/snorm) + 1)/(m+1))
    stemp <- matrix(rep(snorm, each = nrow(y)), nrow=nrow(y))
  }
  if(transformation == "logclr") {
    profiles <- normCounts / rowSums(normCounts)
    tcounts <- logclr(profiles)
  }
  if(transformation == "clr") {
    profiles <- normCounts / rowSums(normCounts)
    tmp <- clr(profiles)
    ## Remove the rmult attributes to keep only a matrix
    attributes(tmp) <- NULL
    tmp <- matrix(tmp, nrow=nrow(normCounts), ncol=ncol(normCounts))
    rownames(tmp) <- rownames(normCounts)
    colnames(tmp) <- colnames(normCounts)
    tcounts <- tmp
  }
  if(transformation == "alr") {
    profiles <- normCounts / rowSums(normCounts)
    tmp <- alr(profiles)
    ## Remove the rmult attributes to keep only a matrix
    attributes(tmp) <- NULL
    tmp <- matrix(tmp, nrow=nrow(normCounts), ncol=ncol(normCounts))
    rownames(tmp) <- rownames(normCounts)
    colnames(tmp) <- colnames(normCounts)
    tcounts <- tmp
  }
  if(transformation == "ilr") {
    profiles <- normCounts / rowSums(normCounts)
    tmp <- ilr(profiles)
    ## Remove the rmult attributes to keep only a matrix
    attributes(tmp) <- NULL
    tmp <- matrix(tmp, nrow=nrow(normCounts), ncol=ncol(normCounts))
    rownames(tmp) <- rownames(normCounts)
    colnames(tmp) <- colnames(normCounts)
    tcounts <- tmp
  }
  ##################################
  ## Old transformations (kept for reference)
  ##################################
  if(transformation == "log") {
    tcounts <- log2(y + 1)
  }
  if(transformation == "normlog") {
    tcounts <- log2(t(t(y)/snorm) + 1)
  }
  if(transformation == "logMeanRef") {
    m <- apply((t(t(y)/snorm))+1, 1, mean)
    tcounts <- log2((t(t(y)/snorm) + 1)/(m+1))
  }
  if(transformation == "logGMeanRef") {
    m <- apply((t(t(y)/snorm))+1, 1, function(x) (prod(x+1))^(1/length(x)))
    tcounts <- log2((t(t(y)/snorm) + 1)/(m))
  }
  if(transformation == "vst") {
    tcounts <- varianceStabilizingTransformation(as.matrix(y), blind=TRUE,
                                                 fitType="parametric")
  }
  if(transformation == "moderatedCPM") {
    tcounts <- cpm(as.matrix(y), normalized.lib.sizes=TRUE, log=TRUE,
                   prior.count=0.25)
  }
  rownames(tcounts) <- rownames(y)
  colnames(tcounts) <- colnames(y)
  return(list(tcounts=tcounts, normCounts=normCounts, snorm=snorm, ellnorm=ellnorm))
}



#' Convert legacy coseq objects
#'
#' Convert legacy coseq S3 class objects to coseqResults S4 class objects
#'
#' @param object Object of S3 class \code{coseq} arising from a call to previous versions of
#' coseq (< 0.99.1)
#' @param digits integer indicating the number of decimal places (round) to retain in results.
#'
#' @return Converted object of S4 class \code{coseqResults} compatible with
#' recent versions of coseq (>= 0.99.1)
#' @export
convertLegacyCoseq <- function(object, digits=3) {
  if(!is(object, "coseq"))
    stop("This function is intended to convert legacy coseq S3 objects \n
          from coseq versions < 0.99.1 into coseqResults S4 objects from \n
          coseq versions >= 0.99.2.")
  if(object$model == "Poisson")
    pp_ICL <- round(object$results$selected.results$probaPost, digits=digits)
  if(object$model == "Normal")
    pp_ICL <- round(object$results$ICL.results$probaPost, digits=digits)
  colnames(pp_ICL) <- paste0("Cluster_", seq_len(ncol(pp_ICL)))
  if(!is.null(rownames(object$y_profiles)))
    rownames(pp_ICL) <- rownames(object$y_profiles)
  if(is.null(rownames(object$y_profiles)))
    rownames(pp_ICL) <- seq_len(nrow(pp_ICL))
  ICL.results <- SummarizedExperiment(pp_ICL,
                                      metadata=list(nbCluster=object$results$nbCluster.all,
                                                    logLike=object$results$logLike.all,
                                                    ICL=object$results$ICL.all,
                                      GaussianModel=object$results$ICL.results$GaussianModel))
  tcountsDF <- as(object$tcounts, "DataFrame")
  y_profilesDF = as(object$y_profiles, "DataFrame")
  pp <- lapply(object$results$all.results, function(x) {
    tmp <- round(x$probaPost, digits=digits)
    if(ncol(tmp)) {
      colnames(tmp) <- paste0("Cluster_", seq_len(ncol(tmp)))
      if(is.null(rownames(tmp))) {
        if(!is.null(rownames(object$y_profiles))) rownames(tmp) <- rownames(object$y_profiles)
        if(is.null(rownames(object$y_profiles))) rownames(tmp) <- seq_len(nrow(tmp))
      }
    }
    return(tmp)
  })
  newobj <- coseqResults(as(ICL.results, "RangedSummarizedExperiment"),
                         allResults = pp,
                         model=object$model,
                         transformation=object$transformation,
                         tcounts=tcountsDF,
                         y_profiles=y_profilesDF,
                         normFactors=numeric())
  return(newobj)
}



#' Calculate the Log Centered Log Ratio (logCLR) transformation
#'
#' @param profiles Matrix of profiles. Note that the presence of 0 values causes an error message to be produced.
#'
#' @return logCLR-transformed profiles
#' @export
logclr <- function(profiles)
{
  d <- ncol(profiles)
  n <- nrow(profiles)
  if(length(which(profiles == 0)))
    stop("log-CLR transformation cannot be performed when there are 0's in the profiles matrix")
  tprofiles <- profiles / exp(rowMeans(log(profiles)))
  tprofiles2 <- matrix(NA, nrow=nrow(tprofiles), ncol=ncol(tprofiles))
  tprofiles2[which(tprofiles <= 1)] <- -(log(1-log(tprofiles[which(tprofiles <= 1)]))^2)
  tprofiles2[which(tprofiles > 1)] <- (log(tprofiles[which(tprofiles > 1)]))^2
  return(tprofiles2)
}
