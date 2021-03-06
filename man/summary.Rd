% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary-functions.R
\name{summary,coseqResults-method}
\alias{summary,coseqResults-method}
\alias{summary}
\alias{summary-methods}
\title{Summarize results from coseq clustering}
\usage{
\S4method{summary}{coseqResults}(object, y_profiles, digits = 3, ...)
}
\arguments{
\item{object}{An object of class \code{"coseqResults"}}

\item{y_profiles}{Data used for clustering if per-cluster means are desired}

\item{digits}{Integer indicating the number of decimal places to be used
for mixture model parameters}

\item{...}{Additional arguments}
}
\value{
Summary of the \code{coseqResults} object.
}
\description{
A function to summarize the clustering results obtained from a Poisson or
Gaussian mixture model estimated using \code{coseq}. In particular,
the function provides the number of clusters selected for the ICL
model selection approach (or alternatively, for the capushe non-asymptotic approach
if K-means clustering is used), number of genes assigned to each cluster, and
if desired the per-gene cluster means.
}
\details{
Provides the following summary of results:

1) Number of clusters and model selection criterion used, if applicable.

2) Number of observations across all clusters with a maximum conditional
probability greater than 90% (and corresponding percentage of total
observations) for the selected model.

3) Number of observations per cluster with a maximum conditional probability
greater than 90% (and corresponding percentage of total observations per
cluster) for the selected model.

4) If desired, the \eqn{\ensuremath\boldsymbol{\mu}}{\mu} values and
\eqn{\ensuremath\boldsymbol{\pi}}{\pi} values for the selected
model in the case of a Gaussian mixture model.
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
\references{
Rau, A. and Maugis-Rabusseau, C. (2017) Transformation and model choice for
co-expression analayis of RNA-seq data. Briefings in Bioinformatics,
doi: http://dx.doi.org/10.1101/065607.

Godichon-Baggioni, A., Maugis-Rabusseau, C. and Rau, A. (2017) Clustering
transformed compositional data using K-means, with applications in gene
expression and bicycle sharing system data. arXiv:1704.06150.
}
\seealso{
\code{\link{coseq}}
}
\author{
Andrea Rau
}
\keyword{methods}
