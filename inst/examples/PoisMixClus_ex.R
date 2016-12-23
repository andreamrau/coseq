## Simulate toy data, n = 300 observations
set.seed(12345)
countmat <- matrix(round(runif(300*4, min=0, max=500)), nrow=300, ncol=4)
countmat <- countmat[which(rowSums(countmat) > 0),]
conds <- rep(c("A","B"), each=2)
rownames(countmat) <- paste0("Gene_", seq_len(nrow(countmat)))

## Run the Poisson mixture model for K = 2,3
run <- coseq(object=countmat, K=2:3, conds=conds, iter=5, model="Poisson")
run
summary(run)
plot(run, graphs="profiles")
