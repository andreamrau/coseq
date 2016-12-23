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
#' @param model Type of mixture model to use (\dQuote{\code{Poisson}} or \dQuote{\code{Normal}})
#' @param transformation Transformation type to be used: \dQuote{\code{voom}}, \dQuote{\code{logRPKM}}
#' (if \code{geneLength} is provided by user), \dQuote{\code{arcsin}}, \dQuote{\code{logit}},
#' \dQuote{\code{logMedianRef}}, \dQuote{\code{profile}}, \dQuote{\code{none}}
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
#' @importFrom stats na.omit
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
#' ## Run the Normal mixture model for K = 2,3,4
#' ## The following are equivalent:
#' run <- coseqRun(y=countmat, K=2:4, iter=5, transformation="arcsin")
#' run <- coseq(object=countmat, K=2:4, iter=5, transformation="arcsin")
#'
coseqRun <- function(y, K, conds=NULL, normFactors="TMM", model="Normal", transformation="arcsin",
                      subset=NULL, meanFilterCutoff=50, modelChoice="ICL",
                      parallel=FALSE, BPPARAM=bpparam(), ...) {

  subset.index <- subset

  ## Parse ellipsis function
  providedArgs <- list(...)
  arg.user <- list(alg.type="EM", init.runs=50, init.type="small-em", init.iter=20,
                   iter=1000, cutoff=0.001, GaussianModel="Gaussian_pk_Lk_Ck",
                   verbose=TRUE, digits=3,
                   Kmin.init="small-em", split.init=FALSE, fixed.lambda=NA,
                   equal.proportions=FALSE, EM.verbose=FALSE, interpretation="sum",
                   geneLength=NA)
  if(model == "Poisson") {
    arg.user$init.runs <- 1
    arg.user$inititer <- 10
    arg.user$cutoff <- 1e-05
    if(is.null(subset)) arg.user$subset.index <- NA
  }
  arg.user[names(providedArgs)] <- providedArgs

  y_profiles <- round(transformRNAseq(y=y, normFactors=normFactors, transformation="profile",
                                       meanFilterCutoff=meanFilterCutoff, verbose=TRUE)$tcounts,
                      digits=arg.user$digits)


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

    ## In case only a subset of data are to be used for analysis
    if(!is.null(subset)) {
      y <- y[subset.index,]
      n <- dim(y)[1]
      cols <- dim(y)[2]
      w <- rowSums(y)
    }
    if(is.null(subset)) {
      n <- dim(y)[1]
      cols <- dim(y)[2]
      w <- rowSums(y)
      subset.index <- NA
    }

    tcounts <- transformRNAseq(y=y, normFactors=normFactors, transformation="none",
                                geneLength=arg.user$geneLength,
                                meanFilterCutoff=meanFilterCutoff, verbose=FALSE)
    conds <- as.vector(conds)

    if(!parallel) {
      run <- suppressWarnings(PoisMixClusWrapper(y=tcounts$tcounts, gmin=min(K),
                                                 gmax=max(K), conds=conds,
                                                 norm=tcounts$ellnorm / sum(tcounts$ellnorm),
                                                 gmin.init.type=arg.user$Kmin.init,
                                                 split.init=arg.user$split.init, subset.index=subset.index,
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
                         subset.index=subset.index,
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
        tmp <- bplapply(remainingK, function(ii) {
          cat("Running K =", ii, "...\n")
          res <- PoisMixClus(g=as.numeric(ii), y=tcounts$tcounts,
                             conds=conds, norm=tcounts$ellnorm / sum(tcounts$ellnorm),
                             init.type=arg.user$Kmin.init,
                             subset.index=subset.index, wrapper=TRUE,
                             prev.probaPost=NA, prev.labels=NA,
                             init.runs = arg.user$init.runs, init.iter = arg.user$init.iter,
                             alg.type = arg.user$alg.type, cutoff = arg.user$cutoff,
                             iter = arg.user$iter, fixed.lambda = arg.user$fixed.lambda,
                             equal.proportions = arg.user$equal.proportions,
                             verbose = arg.user$verbose,
                             interpretation = arg.user$interpretation,
                             EM.verbose = arg.user$EM.verbose)
          return(res)}, BPPARAM=BPPARAM)
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
      tmp <- x$probaPost
      colnames(tmp) <- paste0("Cluster_", seq_len(ncol(tmp)))
      rownames(tmp) <- rownames(tcounts$tcounts)
      return(tmp)
    })
    names(pp) <- names(run$all.results)


    ## Format HTSCluster results
    ICL.results <- SummarizedExperiment(final.results$probaPost, metadata=list(nbCluster=nbClust.all,
                                                                               logLike=logLike.all,ICL=ICL.all))

    tcountsDF = as(tcounts$tcounts, "DataFrame")
    y_profilesDF = as(tcounts$tcounts/rowSums(tcounts$tcounts), "DataFrame")
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

    tcounts <- transformRNAseq(y=y, normFactors=normFactors, transformation=transformation,
                                geneLength=arg.user$geneLength,
                                meanFilterCutoff=meanFilterCutoff, verbose=FALSE)
    run <- NormMixClus(y_profiles=tcounts$tcounts, K=K, subset=subset.index,
                       parallel=parallel,
                       BPPARAM=BPPARAM, alg.type=arg.user$alg.type,
                       init.runs=arg.user$init.runs,
                       init.type=arg.user$init.type, init.iter=arg.user$init.iter,
                       iter=arg.user$iter, cutoff=arg.user$cutoff,
                       verbose=arg.user$verbose, digits=arg.user$digits,
                       GaussianModel=arg.user$GaussianModel)
  }


  if(!is.null(subset)) {
    tcounts$tcounts <- tcounts$tcounts[subset.index,]
    y_profiles <- y_profiles[subset.index,]
  }

  ####################################
  ## RETURN RESULTS
  ####################################

  ICL.results <- SummarizedExperiment(assay(run), metadata=metadata(run))
  all.results <- coseqFullResults(run)
  tcountsDF <- as(tcounts$tcounts, "DataFrame")
  y_profilesDF <- as(y_profiles, "DataFrame")

  RESULTS <- coseqResults(as(ICL.results, "RangedSummarizedExperiment"),
                          allResults=all.results,
                          model=model,
                          transformation=transformation,
                          tcounts=tcountsDF,
                          y_profiles=y_profilesDF,
                          normFactors=tcounts$snorm)
  return(RESULTS)
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
#' \dQuote{\code{none}}
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
#'
#' @importFrom edgeR calcNormFactors
#' @importFrom edgeR cpm
#' @importFrom HTSFilter HTSBasicFilter
#' @importFrom DESeq2 varianceStabilizingTransformation
#' @importFrom stats median


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