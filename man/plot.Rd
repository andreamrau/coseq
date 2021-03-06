% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aaa-AllGenerics.R, R/plot-functions.R
\name{plot}
\alias{plot}
\alias{plot,coseqResults-method}
\alias{plot-methods}
\alias{coseqGlobalPlots}
\alias{coseqModelPlots}
\title{Visualize results from coseq clustering}
\usage{
plot(x, ...)

\S4method{plot}{coseqResults}(
  x,
  y_profiles = NULL,
  K = NULL,
  threshold = 0.8,
  conds = NULL,
  average_over_conds = FALSE,
  collapse_reps = "none",
  graphs = c("logLike", "ICL", "profiles", "boxplots", "probapost_boxplots",
    "probapost_barplots", "probapost_histogram"),
  order = FALSE,
  profiles_order = NULL,
  n_row = NULL,
  n_col = NULL,
  add_lines = TRUE,
  ...
)

coseqGlobalPlots(object, graphs = c("logLike", "ICL"), ...)

coseqModelPlots(
  probaPost,
  y_profiles,
  K = NULL,
  threshold = 0.8,
  conds = NULL,
  collapse_reps = "none",
  graphs = c("profiles", "boxplots", "probapost_boxplots", "probapost_barplots",
    "probapost_histogram"),
  order = FALSE,
  profiles_order = NULL,
  n_row = NULL,
  n_col = NULL,
  add_lines = TRUE,
  ...
)
}
\arguments{
\item{x}{An object of class \code{"coseqResults"}}

\item{...}{Additional optional plotting arguments (e.g., xlab, ylab, use_sample_names, facet_labels)}

\item{y_profiles}{y (\emph{n} x \emph{q}) matrix of observed profiles for \emph{n}
observations and \emph{q} variables to be used for graphing results (optional for
\code{logLike}, \code{ICL}, \code{probapost_boxplots}, and \code{probapost_barplots},
and by default takes value \code{x$tcounts} if \code{NULL})}

\item{K}{If desired, the specific model to use for plotting (or the specific cluster number(s)
to use for plotting in the case of \code{coseqModelPlots}). If \code{NULL},
all clusters will be visualized, and the model chosen by ICL will be plotted}

\item{threshold}{Threshold used for maximum conditional probability; only observations
with maximum conditional probability greater than this threshold are visualized}

\item{conds}{Condition labels, if desired}

\item{average_over_conds}{If \code{TRUE}, average values of \code{y_profiles} within
each condition identified by \code{conds} for the \code{profiles} and \code{boxplots}
plots. This argument is redundant to \code{collapse_reps = "sum"}, and \code{collapse_reps}
should be used instead.}

\item{collapse_reps}{If \code{"none"}, display all replicates. If \code{"sum"}, collapse replicates
within each condition by summing their profiles If \code{"average"}, collapse replicates within
each condition by averaging their profiles. For highly unbalanced experimental designs, using
\code{"average"} will likely provide more easily interpretable plots.}

\item{graphs}{Graphs to be produced, one (or more) of the following:
\code{"logLike"} (log-likelihood plotted versus number of clusters),
\code{"ICL"} (ICL plotted versus number of clusters),
\code{"profiles"} (line plots of profiles in each cluster), \code{"boxplots"}
(boxplots of profiles in each cluster), \code{"probapost_boxplots"} (boxplots of
maximum conditional probabilities per cluster), \code{"probapost_barplots"}
(number of observations with a maximum conditional probability greater than
\code{threshold} per cluster), \code{"probapost_histogram"} (histogram of maximum
conditional probabilities over all clusters)}

\item{order}{If \code{TRUE}, order clusters in \code{probapost_boxplot} by median and
\code{probapost_barplot} by number of observations with maximum conditional probability
greater than \code{threshold}}

\item{profiles_order}{If \code{NULL} or \code{FALSE}, line plots and boxplots of profiles are
plotted sequentially by cluster number (K=1, K=2, ...). If \code{TRUE}, line plots and boxplots of
profiles are plotted in an automatically calculated order (according to the Euclidean distance
between cluster means) to plot clusters with similar mean profiles next to one another.
Otherwise, the user may provide a vector (of length equal to the number of clusters in the
given model) providing the desired order of plots.}

\item{n_row}{Number of rows for plotting layout of line plots and boxplots of profiles.}

\item{n_col}{Number of columns for plotting layout of line plots and boxplots of profiles.}

\item{add_lines}{If \code{TRUE}, add red lines representing means to boxplots; if \code{FALSE},
these will be suppressed.}

\item{object}{An object of class \code{"RangedSummarizedExperiment"} arising from a call to
\code{NormMixClus}}

\item{probaPost}{Matrix or data.frame of dimension (\emph{n} x \emph{K}) containing the
conditional probilities of cluster membership for \emph{n} genes in \emph{K} clusters
arising from a mixture model}
}
\value{
Named list of plots of the \code{coseqResults} object.
}
\description{
Plot a coseqResults object.
}
\examples{
## Simulate toy data, n = 300 observations
set.seed(12345)
countmat <- matrix(runif(300*4, min=0, max=500), nrow=300, ncol=4)
countmat <- countmat[which(rowSums(countmat) > 0),]
conds <- rep(c("A","B","C","D"), each=2)

## Run the Normal mixture model for K = 2,3,4
run_arcsin <- coseq(object=countmat, K=2:4, iter=5, transformation="arcsin",
                    model="Normal", seed=12345)
run_arcsin

## Plot and summarize results
plot(run_arcsin)
summary(run_arcsin)

## Compare ARI values for all models (no plot generated here)
ARI <- compareARI(run_arcsin, plot=FALSE)

## Compare ICL values for models with arcsin and logit transformations
run_logit <- coseq(object=countmat, K=2:4, iter=5, transformation="logit",
                   model="Normal")
compareICL(list(run_arcsin, run_logit))

## Use accessor functions to explore results
clusters(run_arcsin)
likelihood(run_arcsin)
nbCluster(run_arcsin)
ICL(run_arcsin)

## Examine transformed counts and profiles used for graphing
tcounts(run_arcsin)
profiles(run_arcsin)

## Run the K-means algorithm for logclr profiles for K = 2,..., 20
run_kmeans <- coseq(object=countmat, K=2:20, transformation="logclr",
                    model="kmeans")
run_kmeans
}
\author{
Andrea Rau, Cathy Maugis-Rabusseau
}
