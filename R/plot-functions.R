#' Visualize results from coseq clustering
#'
#' Plot a coseqResults object.
#'
#' @rdname plot
#' @aliases
#' plot
#' plot-methods
#' plot,coseqResults-method
#'
#' @param x An object of class \code{"coseqResults"}
#' @param y_profiles y (\emph{n} x \emph{q}) matrix of observed profiles for \emph{n}
#' observations and \emph{q} variables to be used for graphing results (optional for
#' \code{logLike}, \code{ICL}, \code{probapost_boxplots}, and \code{probapost_barplots},
#' and by default takes value \code{x$tcounts} if \code{NULL})
#' @param K If desired, the specific model to use for plotting (or the specific cluster number(s)
#' to use for plotting in the case of \code{coseqModelPlots}). If \code{NULL},
#' all clusters will be visualized, and the model chosen by ICL will be plotted
#' @param threshold Threshold used for maximum conditional probability; only observations
#' with maximum conditional probability greater than this threshold are visualized
#' @param conds Condition labels, if desired
#' @param average_over_conds If \code{TRUE}, average values of \code{y_profiles} within
#' each condition identified by \code{conds} for the \code{profiles} and \code{boxplots}
#' plots. This argument is redundant to \code{collapse_reps = "sum"}, and \code{collapse_reps}
#' should be used instead.
#' @param collapse_reps If \code{"none"}, display all replicates. If \code{"sum"}, collapse replicates
#' within each condition by summing their profiles If \code{"average"}, collapse replicates within
#' each condition by averaging their profiles. For highly unbalanced experimental designs, using
#' \code{"average"} will likely provide more easily interpretable plots.
#' @param graphs Graphs to be produced, one (or more) of the following:
#' \code{"logLike"} (log-likelihood plotted versus number of clusters),
#' \code{"ICL"} (ICL plotted versus number of clusters),
#' \code{"profiles"} (line plots of profiles in each cluster), \code{"boxplots"}
#' (boxplots of profiles in each cluster), \code{"probapost_boxplots"} (boxplots of
#' maximum conditional probabilities per cluster), \code{"probapost_barplots"}
#' (number of observations with a maximum conditional probability greater than
#' \code{threshold} per cluster), \code{"probapost_histogram"} (histogram of maximum
#' conditional probabilities over all clusters)
#' @param order If \code{TRUE}, order clusters in \code{probapost_boxplot} by median and
#' \code{probapost_barplot} by number of observations with maximum conditional probability
#' greater than \code{threshold}
#' @param profiles_order If \code{NULL} or \code{FALSE}, line plots and boxplots of profiles are
#' plotted sequentially by cluster number (K=1, K=2, ...). If \code{TRUE}, line plots and boxplots of
#' profiles are plotted in an automatically calculated order (according to the Euclidean distance
#' between cluster means) to plot clusters with similar mean profiles next to one another.
#' Otherwise, the user may provide a vector (of length equal to the number of clusters in the
#' given model) providing the desired order of plots.
#' @param n_row Number of rows for plotting layout of line plots and boxplots of profiles.
#' @param n_col Number of columns for plotting layout of line plots and boxplots of profiles.
#' @param ...  Additional optional plotting arguments (e.g., xlab, ylab, use_sample_names, facet_labels)
#' @param object An object of class \code{"RangedSummarizedExperiment"} arising from a call to
#' \code{NormMixClus}
#' @param probaPost Matrix or data.frame of dimension (\emph{n} x \emph{K}) containing the
#' conditional probilities of cluster membership for \emph{n} genes in \emph{K} clusters
#' arising from a mixture model
#' @param add_lines If \code{TRUE}, add red lines representing means to boxplots; if \code{FALSE},
#' these will be suppressed.
#'
#' @return Named list of plots of the \code{coseqResults} object.
#'
#' @author Andrea Rau, Cathy Maugis-Rabusseau
#' @example inst/examples/coseq-package.R
#'
#' @export
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics boxplot
#' @importFrom graphics axis
#' @importFrom graphics barplot
#' @importFrom graphics matplot boxplot
#' @importFrom grDevices heat.colors
#' @importFrom scales alpha
#' @importFrom stats hclust dist na.omit
#' @import ggplot2
setMethod(f="plot", signature(x="coseqResults"),
          definition=function(x, y_profiles=NULL, K=NULL, threshold=0.8, conds=NULL,
                              average_over_conds=FALSE,
                              collapse_reps = "none",
                              graphs=c("logLike", "ICL",
                                      "profiles", "boxplots", "probapost_boxplots",
                                       "probapost_barplots", "probapost_histogram"),
                              order=FALSE, profiles_order=NULL, n_row=NULL, n_col=NULL,
                              add_lines = TRUE, ...) {
            # x <- object
            graph_objects <- c()
            ## Parse ellipsis function
            arg.user <- list(...)
            if(is.null(arg.user$alpha)) arg.user$alpha<-0.3

            object <- x
            if(is.null(y_profiles)) y_profiles <- profiles(object)
            ## Check on the length of conds, if provided
            if(!is.null(conds)) {
              if(length(conds) != ncol(y_profiles))
                stop("conds should be a vector of the same length as the number of samples used to create the coseq object.")
            }
            if("logLike" %in% graphs | "ICL" %in% graphs) {
              if(model(object) != "kmeans") {
                globalPlots <- coseqGlobalPlots(object, K=K, threshold=threshold, conds=conds,
                                                graphs=graphs,
                                                order=order, profiles_order=profiles_order, n_row=n_row,
                                                n_col=n_col, ...)
                graph_objects <- c(graph_objects, globalPlots)
              }
              if(model(object) == "kmeans")
                 print("Note: Log-likelihood and ICL plots are not available for coseqResults objects using K-means.")
            }

            ## Model-specific plots
            if(sum(graphs %in% c("profiles", "boxplots", "probapost_boxplots", "probapost_barplots",
                                 "probapost_histogram"))) {
              if(sum(graphs %in% c("profiles", "boxplots"))) {
                if(is.null(y_profiles)) stop("y_profiles needed to plot selected graphs")
                if(!nrow(y_profiles)) stop("y_profiles needed to plot selected graphs")
              }
              if(is.null(K) ) xx <- assay(object)
              if(!is.null(K)) {
                if(K == "ICL") xx <- assay(object)
                if(K != "ICL") {
                  if(!length(which(names(coseqFullResults(object)) == paste0("K=",K))))
                    stop("Selected model was not estimated by coseq")
                  if(length(K) > 1)
                    stop("K must be NULL, a single value, or ICL")
                  xx <- coseqFullResults(object)[[which(names(coseqFullResults(object)) == paste0("K=",K))]]
                }
              }

              if(model(object) == "kmeans") {
                probacalc <- kmeansProbaPost(clusters=apply(xx, 1, which.max),
                                             tcounts=as.matrix(as.data.frame(tcounts(object))))
                xx <- xx * probacalc
              }
              if(average_over_conds) {
                message("The average_over_conds argument is deprecated, and collapse_reps has been set to 'sum'.

(Note that this was the default behavior of average_over_conds in previous versions.
In applications where the number of replicates per condition is unbalanced, we suggest
using collapse_reps = 'average' instead.)")
                collapse_reps <- "sum"
              }
              modelPlots <- coseqModelPlots(probaPost=xx, y_profiles=y_profiles, K=NULL, threshold=threshold, conds=conds,
                   collapse_reps=collapse_reps,
                   graphs=graphs, order = order, alpha=arg.user$alpha,
                   profiles_order=profiles_order, n_row = n_row, n_col = n_col, add_lines = add_lines, ...)
              graph_objects <- c(graph_objects, modelPlots)
            }
        return(graph_objects)
})


#' @export
#' @rdname plot
coseqGlobalPlots <- function(object, graphs=c("logLike", "ICL"), ...) {
  if(is.null(metadata(object)$ICL))
    stop("This function only defined for objects resulting from a coseq analysis.")

  graph_objects <- c()

  ## Parse ellipsis function
  arg.user <- list(...)
  if(is.null(arg.user$alpha)) arg.user$alpha <- 0.3

  pl_data <- na.omit(data.frame(Cluster = nbCluster(object),
                        logLike = likelihood(object),
                        ICL = ICL(object)))

  ## Likelihood plot
  if("logLike" %in% graphs) {
    gg <- ggplot(pl_data, aes_string(x="Cluster", y="logLike")) +
      geom_point() + geom_line() +
      scale_y_continuous(name = ifelse(is.null(arg.user$ylab), "Log-likelihood", arg.user$ylab))
#    print(gg)
    graph_objects$logLike <- gg
  }

  ## ICL plot
  if("ICL" %in% graphs) {
    gg <- ggplot(pl_data, aes_string(x="Cluster", y="ICL")) +
      geom_point() + geom_line()
#    print(gg)
    graph_objects$ICL <- gg
  }
  return(graph_objects)
}


#' @export
#' @rdname plot
coseqModelPlots <- function(probaPost, y_profiles, K=NULL, threshold=0.8, conds=NULL,
                               collapse_reps="none",
                               graphs=c("profiles", "boxplots",
                                        "probapost_boxplots",
                                        "probapost_barplots",
                                        "probapost_histogram"),
                               order = FALSE, profiles_order=NULL,
                               n_row=NULL, n_col=NULL, add_lines = TRUE, ...) {

  graph_objects <- c()

  object <- probaPost
  if(nrow(object) != nrow(y_profiles)) stop("Something is wrong: number of rows in object do not match
                                            number of rows in y_profiles")

  labels <- apply(object, 1, which.max)
  proba <- apply(object, 1, max)

  if(length(conds)) {
    conds <- as.factor(conds)
    conds_vec <- rep(conds, each=nrow(y_profiles))
  }
  if(!length(conds)) conds_vec <- rep(NA, nrow(y_profiles)*ncol(y_profiles))

  ## Parse ellipsis function
  arg.user <- list(...)
  if(is.null(arg.user$alpha)) arg.user$alpha<-0.3

  rn <- rownames(y_profiles)
  cn <- colnames(y_profiles)
  if(is.null(rn)) rn <- seq_len(nrow(y_profiles))
  if(is.null(cn)) cn <- seq_len(ncol(y_profiles))

  #####################################################
  ## SET UP PLOTTING DATA.FRAME
  #####################################################
  if(collapse_reps == "none") {

    pl_data <- data.frame(ID=rep(rn, times=ncol(y_profiles)),
                          y_prof=matrix(as.matrix(as.data.frame(y_profiles)), ncol=1),
                          col_num=rep(seq_len(ncol(y_profiles)), each=nrow(y_profiles)),
                          col_nam=rep(cn, each=nrow(y_profiles)),
                          conds=conds_vec,
                          labels=rep(labels, times=ncol(y_profiles)),
                          proba=rep(proba, times=ncol(y_profiles)))
  }

  if(collapse_reps == "sum") {
    if(!length(conds)) stop("Conds argument needed when collapse_reps == 'sum'")
    y_profiles_c <- t(rowsum(t(as.matrix(as.data.frame(y_profiles))), conds))
    conds_vec <- factor(rep(colnames(y_profiles_c), each=nrow(y_profiles_c)),
                        levels=levels(conds))
    pl_data <- data.frame(ID=ifelse(rep(length(rownames(y_profiles_c))==0, nrow(y_profiles_c)),
                                    rep(seq_len(nrow(y_profiles_c)), times=ncol(y_profiles_c)),
                                    rownames(y_profiles_c)),
                          y_prof=matrix(as.matrix(y_profiles_c), ncol=1),
                          col_num=rep(seq_len(ncol(y_profiles_c)), each=nrow(y_profiles_c)),
                          col_nam=rep(colnames(y_profiles_c), each=nrow(y_profiles_c)),
                          conds=conds_vec,
                          labels=rep(labels, times=ncol(y_profiles_c)),
                          proba=rep(proba, times=ncol(y_profiles_c)))
  }

  if(collapse_reps == "average") {
    if(!length(conds)) stop("Conds argument needed when collapse_reps == 'average")
    y_profiles_c <- t(rowsum(t(as.matrix(as.data.frame(y_profiles))), conds)/as.numeric(table(conds)))
    conds_vec <- factor(rep(colnames(y_profiles_c), each=nrow(y_profiles_c)),
                        levels=levels(conds))
    pl_data <- data.frame(ID=ifelse(rep(length(rownames(y_profiles_c))==0, nrow(y_profiles_c)),
                                    rep(seq_len(nrow(y_profiles_c)), times=ncol(y_profiles_c)),
                                    rownames(y_profiles_c)),
                          y_prof=matrix(as.matrix(y_profiles_c), ncol=1),
                          col_num=rep(seq_len(ncol(y_profiles_c)), each=nrow(y_profiles_c)),
                          col_nam=rep(colnames(y_profiles_c), each=nrow(y_profiles_c)),
                          conds=conds_vec,
                          labels=rep(labels, times=ncol(y_profiles_c)),
                          proba=rep(proba, times=ncol(y_profiles_c)))
  }

  y_profiles <- as.data.frame(y_profiles)
  ## Reorder clusters if desired
  pl_data$labels <- factor(pl_data$labels)
  if(!is.null(profiles_order)) {
    if(length(profiles_order) == 1 & is.logical(profiles_order)) {
      ## If actual normal mixture model
      meanmat <- matrix(0, nrow=ncol(object), ncol=ncol(y_profiles))
      for (k in seq_len(ncol(probaPost))) {
        meanmat[k,]<- apply(probaPost[,k] * y_profiles,2,sum) / sum(probaPost[,k])
      }
      ord <- hclust(dist(meanmat))$order
    } else if(length(profiles_order) == length(unique(labels))) {
      ord <- profiles_order
    } else {
      ord <- sort(unique(labels))
    }
    pl_data$labels <- factor(pl_data$labels, levels=ord)
  }

  pl_data_tmp <- pl_data

  #####################################################
  ## PROFILE PLOTS
  #####################################################
  if("profiles" %in% graphs) {
    ## For one specific value of K
    if(!is.null(K) & length(K) == 1) pl_data_tmp <- pl_data[which(pl_data$labels == K),]
    ## For a subset of values of K
    if(!is.null(K) & length(K) > 1) {
      pl_data_tmp <- pl_data[which(pl_data$labels %in% K),]
      pl_data_tmp$labels <- droplevels(pl_data_tmp$labels)
    }

    ## Print all on the same page
      g1 <- ggplot(pl_data_tmp[which(pl_data_tmp$proba > threshold),]) +
        geom_line(colour=alpha("black", arg.user$alpha),
                  aes_string(x=ifelse(collapse_reps != "none", "conds", "col_num"), y="y_prof", group="ID")) +
        geom_line(data=pl_data_tmp[which(pl_data_tmp$proba < threshold),],
                  colour=alpha("red", arg.user$alpha),
                  aes_string(x=ifelse(collapse_reps != "none", "conds", "col_num"), y="y_prof", group="ID")) +
        theme_bw()
      if(collapse_reps == "none") g1 <- g1 +
          scale_y_continuous(name=ifelse(is.null(arg.user$ylab), "Expression profiles", arg.user$ylab)) +
          scale_x_continuous(name=ifelse(is.null(arg.user$xlab), "Sample number", arg.user$xlab))
      if(collapse_reps == "sum") g1 <- g1 +
          scale_y_continuous(name=ifelse(is.null(arg.user$ylab), "Summed expression profiles", arg.user$ylab)) +
          scale_x_discrete(name=ifelse(is.null(arg.user$xlab), "Conditions", arg.user$xlab))
      if(collapse_reps == "average") g1 <- g1 +
          scale_y_continuous(name=ifelse(is.null(arg.user$ylab), "Average expression profiles", arg.user$ylab)) +
          scale_x_discrete(name=ifelse(is.null(arg.user$xlab), "Conditions", arg.user$xlab))

      if(!is.null(K) & length(K) == 1) g1 <- g1 + ggtitle(paste("Cluster", K))
      if(is.null(K)) {
        if(is.null(arg.user$facet_labels)) g1 <- g1 + facet_wrap(~labels,
                                                                 nrow=n_row, ncol=n_col)
        if(!is.null(arg.user$facet_labels))
          g1 <- g1 + facet_wrap(~labels, labeller=labeller(labels = arg.user$facet_labels),
                                nrow=n_row, ncol=n_col)

      graph_objects$profiles <- g1
    }

    # ## Print on different pages
    # if(!is.null(n_row)) {
    #   g2_list <- lapply(levels(pl_data_tmp$labels), function(.x) {
    #     pl_data2 <- pl_data_tmp[which(pl_data_tmp$labels == .x),]
    #     pl_data2$labels <- factor(pl_data2$labels, levels=.x)
    #     g2 <- ggplot(pl_data2[which(pl_data2$proba > threshold),]) +
    #       geom_line(colour=alpha("black", arg.user$alpha),
    #                 aes_string(x=ifelse(collapse_reps != "none", "conds", "col_num"), y="y_prof", group="ID")) +
    #       geom_line(data=pl_data2[which(pl_data2$proba < threshold),],
    #                 colour=alpha("red", arg.user$alpha),
    #                 aes_string(x=ifelse(collapse_reps != "none", "conds", "col_num"), y="y_prof", group="ID")) +
    #       theme_bw()
    #     if(is.null(arg.user$facet_labels)) g2 <- g2 + facet_wrap(~labels)
    #     if(!is.null(arg.user$facet_labels))
    #       g2 <- g2 + facet_wrap(~labels, labeller=labeller(labels = arg.user$facet_labels))
    #     if(collapse_reps == "none") g2 <- g2 +
    #       scale_y_continuous(name=ifelse(is.null(arg.user$ylab), "Expression profiles", arg.user$ylab)) +
    #       scale_x_continuous(name=ifelse(is.null(arg.user$xlab), "Sample number", arg.user$xlab))
    #     if(collapse_reps == "sum") g2 <- g2 +
    #       scale_y_continuous(name=ifelse(is.null(arg.user$ylab), "Summed expression profiles", arg.user$ylab)) +
    #       scale_x_discrete(name=ifelse(is.null(arg.user$xlab), "Conditions", arg.user$xlab))
    #     if(collapse_reps == "average") g2 <- g2 +
    #       scale_y_continuous(name=ifelse(is.null(arg.user$ylab), "Average expression profiles", arg.user$ylab)) +
    #       scale_x_discrete(name=ifelse(is.null(arg.user$xlab), "Conditions", arg.user$xlab))
    #     return(g2)
    #   })
    #   g2 <- marrangeGrob(g2_list, ncol=n_col, nrow=n_row)
    #   graph_objects$profiles <- g2
    # }
  }

  #####################################################
  ## PROFILE BOXPLOTS
  #####################################################
  if("boxplots" %in% graphs) {
    ## Use sample number or column names?
    if(!is.null(arg.user$use_sample_names)) {
      levels(pl_data_tmp$col_nam) <- levels(pl_data_tmp$col_nam)[match(cn, levels(pl_data_tmp$col_nam))]
      pl_data_tmp$col_num <- pl_data_tmp$col_nam
    }

    ## For one specific value of K
    if(!is.null(K) & length(K) == 1) pl_data_tmp <- pl_data[which(pl_data$labels == K),]
    ## For a subset of values of K
    if(!is.null(K) & length(K) > 1) {
      pl_data_tmp <- pl_data[which(pl_data$labels %in% K),]
      pl_data_tmp$labels <- droplevels(pl_data_tmp$labels)
    }
    pl_data_tmp$col_num <- factor(pl_data_tmp$col_num)
    pl_data_tmp$conds <- factor(pl_data_tmp$conds)
    pl_data$col_num <- factor(pl_data$col_num)

    ## Print all on the same page
 #   if(is.null(n_row)) {
      g3 <- ggplot(pl_data_tmp,
                   aes_string(x=ifelse(collapse_reps == "none", "col_num", "conds"), y="y_prof"))
      if(!length(conds)) {
        g3 <- g3 +  geom_boxplot()
      }
      if(length(conds)) {
        g3 <- g3 + geom_boxplot(aes_string(fill="conds")) +
          scale_fill_discrete(name="Conditions")
      }
      if(add_lines == TRUE) {
        g3 <- g3 + stat_summary(fun=mean, geom="line", aes(group=1), colour="red")  +
          stat_summary(fun=mean, geom="point", colour="red")
      }
      if(!is.null(K) & length(K) > 1) g3 <- g3 +  ggtitle(paste("Cluster", K))
      if(collapse_reps == "none") g3 <- g3 +
        scale_y_continuous(name=ifelse(is.null(arg.user$ylab), "Expression profiles", arg.user$ylab)) +
        scale_x_discrete(name=ifelse(is.null(arg.user$xlab), "Sample number", arg.user$xlab))
      if(collapse_reps == "sum") g3 <- g3 +
        scale_y_continuous(name=ifelse(is.null(arg.user$ylab), "Summed expression profiles", arg.user$ylab)) +
        scale_x_discrete(name=ifelse(is.null(arg.user$xlab), "Conditions", arg.user$xlab))
      if(collapse_reps == "average") g3 <- g3 +
        scale_y_continuous(name=ifelse(is.null(arg.user$ylab), "Average expression profiles", arg.user$ylab)) +
        scale_x_discrete(name=ifelse(is.null(arg.user$xlab), "Conditions", arg.user$xlab))
      if(is.null(K)) {
        if(is.null(arg.user$facet_labels)) g3 <- g3 + facet_wrap(~labels, nrow=n_row, ncol=n_col)
        if(!is.null(arg.user$facet_labels))
          g3 <- g3 + facet_wrap(~labels, labeller=labeller(labels = arg.user$facet_labels),
                                nrow=n_row, ncol = n_col)
      }
 #     print(g3)
      graph_objects$boxplots <- g3

 #    ## Print on different pages
 #    if(!is.null(n_row)) {
 #      g4_list <- lapply(levels(pl_data_tmp$labels), function(.x) {
 #        pl_data2 <- pl_data_tmp[which(pl_data_tmp$labels == .x),]
 #        pl_data2$labels <- factor(pl_data2$labels, levels=.x)
 #        g4 <- ggplot(pl_data2, aes_string(x=ifelse(collapse_reps == "none", "conds", "col_num"), y="y_prof"))
 #        if(!length(conds)) {
 #          g4 <- g4 + geom_boxplot()
 #        }
 #        if(length(conds)) {
 #          g4 <- g4 + geom_boxplot(aes_string(fill="conds")) +
 #            scale_fill_discrete(name=ifelse(is.null(arg.user$xlab), "Conditions", arg.user$xlab))
 #        }
 #        g4 <- g4 + stat_summary(fun=mean, geom="line", aes(group=1), colour="red")  +
 #          stat_summary(fun=mean, geom="point", colour="red")
 #
 #        if(is.null(arg.user$facet_labels)) g4 <- g4 + facet_wrap(~labels)
 #        if(!is.null(arg.user$facet_labels))
 #          g4 <- g4 + facet_wrap(~labels, labeller=labeller(labels = arg.user$facet_labels))
 #        if(collapse_reps == "none") g4 <- g4 +
 #          scale_y_continuous(name=ifelse(is.null(arg.user$ylab), "Expression profiles", arg.user$ylab)) +
 #          scale_x_discrete(name=ifelse(is.null(arg.user$xlab), "Sample number", arg.user$xlab))
 #        if(collapse_reps == "sum") g4 <- g4 +
 #          scale_y_continuous(name=ifelse(is.null(arg.user$ylab), "Summed expression profiles", arg.user$ylab)) +
 #          scale_x_discrete(name=ifelse(is.null(arg.user$xlab), "Conditions", arg.user$xlab))
 #        if(collapse_reps == "average") g4 <- g4 +
 #          scale_y_continuous(name=ifelse(is.null(arg.user$ylab), "Average expression profiles", arg.user$ylab)) +
 #          scale_x_discrete(name=ifelse(is.null(arg.user$xlab), "Conditions", arg.user$xlab))
 #        return(g4)
 #      })
 #      g4 <- marrangeGrob(g4_list, ncol=n_col, nrow=n_row)
 # #     print(g4)
 #      graph_objects$boxplots <- g4
 #    }
  }


  #####################################################
  ## PROBAPOST BOXPLOTS
  #####################################################

  if("probapost_boxplots" %in% graphs) {
    pl_data <- data.frame(ID=rep(rn, times=ncol(y_profiles)),
                          y_prof=matrix(as.matrix(as.data.frame(y_profiles)), ncol=1),
                          col_num=rep(seq_len(ncol(y_profiles)), each=nrow(y_profiles)),
                          col_nam=rep(cn, each=nrow(y_profiles)),
                          labels=rep(labels, times=ncol(y_profiles)),
                          proba=rep(proba, times=ncol(y_profiles)))
    pl_data$labels <- factor(pl_data$labels)

    ## Single value of K or subset of values
    if(!is.null(K)) pl_data <- pl_data[which(pl_data$labels %in% K),]
    if(is.null(K) & order) {
      A <- boxplot(pl_data$proba~pl_data$label, plot=FALSE)
      J <- sort.int(A$stat[3,],index.return=TRUE,decreasing=TRUE)$ix
      pl_data$labels <- factor(pl_data$labels, levels=J)
    }
    gg <- ggplot(pl_data, aes_string(x="labels", y="proba")) +
        geom_boxplot() +  scale_x_discrete(name="Cluster") +
        scale_y_continuous(name="Max conditional probability")
 #   print(gg)
    graph_objects$probapost_boxplots <- gg
  }

  #####################################################
  ## PROBAPOST BARPLOTS
  #####################################################
  if("probapost_barplots" %in% graphs) {
    pl_data <- data.frame(ID=rep(rn, times=ncol(y_profiles)),
                          y_prof=matrix(as.matrix(as.data.frame(y_profiles)), ncol=1),
                          col_num=rep(seq_len(ncol(y_profiles)), each=nrow(y_profiles)),
                          col_nam=rep(cn, each=nrow(y_profiles)),
                          labels=rep(labels, times=ncol(y_profiles)),
                          proba=rep(proba, times=ncol(y_profiles)))
    pl_data$labels <- factor(pl_data$labels)

    pl_data2 <- pl_data[which(pl_data$col_num==1),]
    pl_data2$goodproba <- factor(ifelse(pl_data2$proba > threshold,
                                        paste(">", threshold), paste("<", threshold)),
                                 levels=c(paste(">", threshold),
                                          paste("<", threshold)))
    ## Single value of K or subset of values
    if(!is.null(K)) pl_data2 <- pl_data2[which(pl_data2$labels %in% K),]
    if(order) {
      pl_data2$labels <- factor(pl_data2$labels,
                               levels=names(sort(table(pl_data2$labels[which(
                                 pl_data2$goodproba == paste(">", threshold))]),
                                 decreasing=TRUE)))
    }

    gg <- ggplot(pl_data2, aes_string(x="labels", fill="goodproba")) +
      geom_bar() +
      scale_fill_brewer(direction=-1, palette="Accent",
                        name="Max\nconditional\nprobability") +
      scale_x_discrete(name="Cluster") +
      scale_y_continuous(name="Number of observations")
 #   print(gg)
    graph_objects$probapost_barplots <- gg
  }


  #####################################################
  ## PROBAPOST HISTOGRAM
  #####################################################
  if("probapost_histogram" %in% graphs) {
    pl_data <- data.frame(ID=rep(rn, times=ncol(y_profiles)),
                          y_prof=matrix(as.matrix(as.data.frame(y_profiles)), ncol=1),
                          col_num=rep(seq_len(ncol(y_profiles)), each=nrow(y_profiles)),
                          col_nam=rep(cn, each=nrow(y_profiles)),
                          labels=rep(labels, times=ncol(y_profiles)),
                          proba=rep(proba, times=ncol(y_profiles)))
    pl_data$labels <- factor(pl_data$labels)

    pl_data_tmp <- pl_data[which(pl_data$col_num == 1),]
    gg <- ggplot(pl_data_tmp, aes_string(x="proba")) +
      geom_histogram(binwidth = 0.01) +
      scale_x_continuous(name = "Maximum conditional probability") + theme_bw()
#    print(gg)
    graph_objects$probapost_histogram <- gg
  }
  return(graph_objects)
}
