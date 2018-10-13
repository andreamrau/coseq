## ---- eval=FALSE---------------------------------------------------------
#  run <- coseq(counts, K=2:25)
#  summary(run)
#  plot(run)

## ---- eval=FALSE---------------------------------------------------------
#  clusters(run)

## ---- message=FALSE------------------------------------------------------
library(coseq)
library(Biobase)
library(corrplot)

data("fietz")
counts <- exprs(fietz)
conds <- pData(fietz)$tissue

## Equivalent to the following:
## counts <- read.table("http://www.math.univ-toulouse.fr/~maugis/coseq/Fietz_mouse_counts.txt",
##                       header=TRUE)
## conds <- c(rep("CP",5),rep("SVZ",5),rep("VZ",5))

## ----eval=FALSE----------------------------------------------------------
#  run <- coseq(..., parallel=TRUE)

## ------------------------------------------------------------------------
runLogCLR <- coseq(counts, K=2:25, transformation="logclr",norm="TMM", 
                      meanFilterCutoff=50, model="kmeans",
                   nstart=1, iter.max=10)

## ------------------------------------------------------------------------
runArcsin <- coseq(counts, K=2:20, model="Normal", transformation="arcsin", 
                   meanFilterCutoff=200, iter=10)
runLogit <- coseq(counts, K=2:20, model="Normal", transformation="logit", 
                  meanFilterCutoff=200, verbose=FALSE, iter=10)

## ------------------------------------------------------------------------
class(runArcsin)
runArcsin

## ---- figure=TRUE, fig.width=6, fig.height=4-----------------------------
compareICL(list(runArcsin, runLogit))

## ---- figure=TRUE, fig.width=7, fig.height=7-----------------------------
compareARI(runArcsin, K=8:12)

## ------------------------------------------------------------------------
summary(runArcsin)

## ---- fig.width=6, fig.height=4------------------------------------------
## To obtain all plots
## plot(runArcsin)
plot(runArcsin, graphs="boxplots")
plot(runArcsin, graphs="boxplots", conds=conds)
plot(runArcsin, graphs="boxplots", conds=conds, collapse_reps = "sum")
plot(runArcsin, graphs="profiles")

plot(runArcsin, graphs="probapost_boxplots")
plot(runArcsin, graphs="probapost_histogram")

## ---- figure=TRUE, fig.width=4, fig.height=4-----------------------------
rho <- NormMixParam(runArcsin)$rho
## Covariance matrix for cluster 1
rho1 <- rho[,,1]
colnames(rho1) <- rownames(rho1) <- paste0(colnames(rho1), "\n", conds)
corrplot(rho1)

## ------------------------------------------------------------------------
labels <- clusters(runArcsin)
table(labels)

metadata(runArcsin)
likelihood(runArcsin)
nbCluster(runArcsin)
ICL(runArcsin)
model(runArcsin)
transformationType(runArcsin)

## ------------------------------------------------------------------------
## arcsine-transformed normalized profiles
tcounts(runArcsin)

## normalized profiles
profiles(runArcsin)

## ------------------------------------------------------------------------
fullres <- coseqFullResults(runArcsin)
class(fullres)
length(fullres)
names(fullres)

## ---- message=FALSE, warning=FALSE---------------------------------------
library(DESeq2)
dds <- DESeqDataSetFromMatrix(counts, DataFrame(group=factor(conds)), ~group)
dds <- DESeq(dds, test="LRT", reduced = ~1)
res <- results(dds)
summary(res)

## By default, alpha = 0.10
run <- coseq(dds, K=2:15, verbose=FALSE)

## The following two lines provide identical results
run <- coseq(dds, K=2:15, alpha=0.05)
run <- coseq(dds, K=2:15, subset=results(dds, alpha=0.05))

## ---- message=FALSE, warning=FALSE---------------------------------------
library(edgeR)
y <- DGEList(counts=counts, group=factor(conds))
y <- calcNormFactors(y)
design <- model.matrix(~conds)
y <- estimateDisp(y, design)

## edgeR: QLF test
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef=2)

## edgeR: LRT test
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)

run <- coseq(counts, K=2:15, subset=lrt, alpha=0.1)
run <- coseq(counts, K=2:15, subset=qlf, alpha=0.1)

## ---- fig.width=6, fig.height=4------------------------------------------
## Simulate toy data, n = 300 observations
set.seed(12345)
countmat <- matrix(runif(300*10, min=0, max=500), nrow=300, ncol=10)
countmat <- countmat[which(rowSums(countmat) > 0),]
conds <- c(rep("A", 4), rep("B", 3), rep("D", 3))

## Run coseq
coexp <- coseq(object=countmat, K=2:15, iter=5, transformation="logclr",
                    model="kmeans")

## Original boxplot
p <- plot(coexp, graphs = "boxplots", conds = conds, 
     profiles_order = 10:1, collapse_reps = "average",
     n_row = 3, n_col = 4)
p$boxplots

## ---- fig.width=6, fig.height=4------------------------------------------
library(ggplot2)
p$boxplot + 
  theme_bw() + 
  coord_fixed(ratio = 20) +
  theme(axis.ticks = element_line(color="black", size=1.5),
        axis.line = element_line(color="black", size=1.5),
        text = element_text(size = 15)) +
  scale_fill_brewer(type = "qual")

