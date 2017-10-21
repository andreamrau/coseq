##-------------------------------------------------------------------------------
context("coseq function tests")
set.seed(12345)
obj <- matrix(round(runif(50, 1, 100)), ncol=5, nrow=10)

test_that("K is a vector of positive integers", {
  expect_error(coseq(obj, K=-3, verbose=FALSE))
  expect_error(coseq(obj, K=matrix(rnorm(6), nrow=3), verbose=FALSE))
  expect_error(coseq(obj, K=c(-3, -2, -1, 0, 1, 2, 3), verbose=FALSE))
  expect_error(coseq(obj, K=c(0, 1, 2, 3), verbose=FALSE))
  expect_error(coseq(obj, K=c(1.5, 3.3, 5.2), verbose=FALSE))
})

obj <- data.frame(obj)

test_that("model, transformation, normFactors, ICL", {
  expect_error(coseq(obj, K=2:4, model="Gamma", verbose=FALSE))
  expect_error(coseq(obj, K=2:4, transformation="other", verbose=FALSE))
  expect_error(coseq(obj, K=2:4, transformation="other", verbose=FALSE))
})

coseq_res <- coseq(obj, K=2:4, norm="none", GaussianModel="Gaussian_pk_Lk_I",
                   verbose=FALSE, model="Normal", transformation="none")

test_that("coseq output", {
  expect_true(length(clusters(coseq_res)) == nrow(obj))
  expect_equal(max(clusters(coseq_res)),
               as.numeric(substr(names(which.min(ICL(coseq_res))), 3, 10)))
  expect_equal(min(clusters(coseq_res)), 1)
  expect_equal(as.numeric(substr(names(which.min(ICL(coseq_res))), 3, 10)),
               max(clusters(coseq_res)))
  expect_equal(length(nbCluster(coseq_res)), length(likelihood(coseq_res)))
  expect_equal(length(ICL(coseq_res)), length(likelihood(coseq_res)))
  expect_equal(nrow(profiles(coseq_res)), nrow(tcounts(coseq_res)))
  expect_equal(ncol(profiles(coseq_res)), ncol(tcounts(coseq_res)))
  expect_true(model(coseq_res) == "Normal")
  expect_true(transformationType(coseq_res) == "none")
  expect_is(coseqFullResults(coseq_res), "list")
})

test_that("compareARI", {
  tmp <- compareARI(coseq_res, plot=FALSE)
  expect_equal(nrow(tmp), length(nbCluster(coseq_res)))
  obj2 <- matrix(rnorm(100), ncol=4)
  expect_error(compareARI(obj2))
})

#-------------------------------------------------------------------------------
context("general function tests")
test_that("compareICL, clusterEntropy, transformRNAseq", {
  coseq_res2 <- coseq(obj, K=2:4, transformation="logMedianRef",
                      GaussianModel="Gaussian_pk_Lk_I", model="Normal",
                     verbose=FALSE, norm="none")
  expect_error(compareICL(list(coseq_res, coseq_res2)))
  expect_equal(length(clusterEntropy(assay(coseq_res))),
               as.numeric(substr(names(which.min(ICL(coseq_res))), 3, 10)))

  objtr <- transformRNAseq(obj, normFactors="none", transformation="arcsin",
                              geneLength=NA, meanFilterCutoff=NULL, verbose=FALSE)
  expect_equal(nrow(objtr$tcounts), nrow(obj))
  expect_equal(ncol(objtr$tcounts), ncol(obj))

})



#-------------------------------------------------------------------------------
context("NormMixClus tests")
test_that("NormMixClus working", {
  obj2 <- transformRNAseq(obj, normFactors="none", transformation="arcsin",
                          geneLength=NA, meanFilterCutoff=NULL, verbose=FALSE)$tcounts
  expect_error(NormMixClusK(obj2, K=-2))
  expect_error(NormMixParam(coseq_res))
  expect_error(coseqModelPlots(coseq_res, graphs="boxplots"))
})


#-------------------------------------------------------------------------------
context("plot and class tests")
test_that("plot and class construction", {
  expect_error(plot(coseq_res, K=-2))
  newobj <- SummarizedExperiment(list(foo=as.matrix(obj)),
                                 colData=DataFrame(group=factor(c(1,1,2,2,2))))
  expect_error(coseqResults(newobj, allResults=list(1,2,3)))
})

