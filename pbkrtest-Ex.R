pkgname <- "pbkrtest"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('pbkrtest')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("KenRog")
### * KenRog

flush(stderr()); flush(stdout())

### Name: KenwardRoger
### Title: Ftest and degrees of freedom based on Kenward-Roger
###   approximation
### Aliases: KRmodcomp KRmodcomp.mer
### Keywords: function

### ** Examples

(fmLarge <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
## removing Day
(fmSmall <- lmer(Reaction ~ 1 + (Days|Subject), sleepstudy))
(KRmodcomp(fmLarge,fmSmall))

## The same test using a restriction matrix
L<-cbind(0,1)
(KRmodcomp(fmLarge,L))



cleanEx()
nameEx("PB_PBmodcomp")
### * PB_PBmodcomp

flush(stderr()); flush(stdout())

### Name: PBmodcomp
### Title: Model comparison of mixed models using parametric bootstrap
###   methods.
### Aliases: PBmodcomp PBmodcomp.mer makePBcluster stopPBcluster
###   plot.XXmodcomp
### Keywords: utilities models inference

### ** Examples

data(beets)
head(beets)
beet0<-lmer(sugpct~block+sow+harvest+(1|block:harvest), data=beets, REML=FALSE)
beet_no.harv <- update(beet0, .~.-harvest)
PBmodcomp(beet0, beet_no.harv, nsim=20)

## Not run: 
##D ## Vanilla
##D PBmodcomp(beet0, beet_no.harv)
##D ## Simulate reference distribution separately:
##D rr <- PBrefdist(beet0, beet_no.harv, nsim=20)
##D PBmodcomp(beet0, beet_no.harv, ref=rr)
## End(Not run)




cleanEx()
nameEx("PB_PBrefdist")
### * PB_PBrefdist

flush(stderr()); flush(stdout())

### Name: PBrefdist
### Title: Calculate reference distribution using parametric bootstrap
### Aliases: PBrefdist PBrefdist.mer
### Keywords: utilities models

### ** Examples

data(beets)
head(beets)
beet0<-lmer(sugpct~block+sow+harvest+(1|block:harvest), data=beets, REML=FALSE)
beet_no.harv <- update(beet0, .~.-harvest)
rr <- PBrefdist(beet0, beet_no.harv, nsim=20)
rr

## Note clearly many more than 10 simulations must be made in practice.

## Computations can be made in parallel using several processors:
## Not run: 
##D library(snow)
##D cl <- makeSOCKcluster(rep("localhost", 4))
##D clusterEvalQ(cl, library(lme4))
##D clusterSetupSPRNG(cl)
##D rr <- PBrefdist(beet0, beet_no.harv, nsim=20)
##D stopCluster(cl)
## End(Not run)
## Above, 4 cpu's are used and 5 simulations are made on each cpu.



cleanEx()
nameEx("beets")
### * beets

flush(stderr()); flush(stdout())

### Name: beets
### Title: Yield and sugar percentage in sugar beets from a split plot
###   experiment.
### Aliases: beets
### Keywords: datasets

### ** Examples

data(beets)
## maybe str(beets) ; plot(beets) ...

beets$bh <- with(beets, interaction(block, harvest))
summary(aov(yield~block+sow+harvest+Error(bh), beets))
summary(aov(sugpct~block+sow+harvest+Error(bh), beets))



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
