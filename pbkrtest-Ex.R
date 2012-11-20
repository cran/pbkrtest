pkgname <- "pbkrtest"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('pbkrtest')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("DATA-beets")
### * DATA-beets

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



cleanEx()
nameEx("DATA-budworm")
### * DATA-budworm

flush(stderr()); flush(stdout())

### Name: budworm
### Title: Effect of Insecticide on survivial of tobacco budworms
### Aliases: budworm
### Keywords: datasets

### ** Examples

data(budworm)

#function to caclulate the empirical logits
empirical.logit<- function(nevent,ntotal) {
y<-log ((nevent+0.5)/(ntotal-nevent+0.5))
y
}


#plot the empirical logits against log-dose

log.dose<-log(budworm$dose)
emp.logit<-empirical.logit(budworm$ndead,budworm$ntotal)
plot(log.dose,emp.logit,type='n',xlab='log-dose',ylab='emprirical logit')
title('budworm: emprirical logits of probability to die ')
male<-budworm$sex=='male'
female<-budworm$sex=='female'
lines(log.dose[male],emp.logit[male],type='b',lty=1,col=1)
lines(log.dose[female],emp.logit[female],type='b',lty=2,col=2)
legend(0.5,2,legend=c('male','female'),lty=c(1,2),col=c(1,2))





## Not run: 
##D * SAS example;
##D data budworm;
##D infile 'budworm.txt' firstobs=2;
##D input sex dose ndead ntotal;
##D run;
## End(Not run)




cleanEx()
nameEx("KR-modcomp")
### * KR-modcomp

flush(stderr()); flush(stdout())

### Name: KenwardRoger
### Title: Ftest and degrees of freedom based on Kenward-Roger
###   approximation
### Aliases: KRmodcomp KRmodcomp.mer KRmodcomp.lmerMod
### Keywords: function

### ** Examples

	(fmLarge <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
	## removing Days
	(fmSmall <- lmer(Reaction ~ 1 + (Days|Subject), sleepstudy))
	anova(fmLarge,fmSmall)
	KRmodcomp(fmLarge,fmSmall)

## The same test using a restriction matrix
L<-cbind(0,1)
KRmodcomp(fmLarge, L)



cleanEx()
nameEx("PB-modcomp")
### * PB-modcomp

flush(stderr()); flush(stdout())

### Name: PBmodcomp
### Title: Model comparison of mixed models using parametric bootstrap
###   methods.
### Aliases: PBmodcomp PBmodcomp.lm PBmodcomp.mer PBmodcomp.lmerMod getLRT
###   getLRT.lm getLRT.mer getLRT.lmerMod plot.XXmodcomp
### Keywords: utilities models inference

### ** Examples

data(beets)
head(beets)
beet0<-lmer(sugpct~block+sow+harvest+(1|block:harvest), data=beets, REML=FALSE)
beet_no.harv <- update(beet0, .~.-harvest)
PBmodcomp(beet0, beet_no.harv, nsim=20)

## Not run: 
##D (fmLarge <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
##D ## removing Days
##D (fmSmall <- lmer(Reaction ~ 1 + (Days|Subject), sleepstudy))
##D anova(fmLarge, fmSmall)
##D PBmodcomp(fmLarge, fmSmall)
##D 
##D ## The same test using a restriction matrix
##D L<-cbind(0,1)
##D PBmodcomp(fmLarge, L)
##D 
##D 
##D ## Vanilla
##D PBmodcomp(beet0, beet_no.harv, nsim=200)
##D 
##D ## Simulate reference distribution separately:
##D refdist <- PBrefdist(beet0, beet_no.harv, nsim=200)
##D PBmodcomp(beet0, beet_no.harv, ref=refdist)
##D 
##D ## Do computations with multiple processors:
##D ## Number of cores:
##D (nc <- detectCores())
##D ## Create clusters
##D cl <- makeCluster(rep("localhost", nc))
##D 
##D ## Then do:
##D PBmodcomp(beet0, beet_no.harv, cl=cl)
##D 
##D ## Or in two steps:
##D refdist <- PBrefdist(beet0, beet_no.harv, nsim=200, cl=cl)
##D PBmodcomp(beet0, beet_no.harv, ref=refdist)
##D 
##D ## It is recommended to stop the clusters before quitting R:
##D stopCluster(cl)
## End(Not run)




cleanEx()
nameEx("PB-refdist")
### * PB-refdist

flush(stderr()); flush(stdout())

### Name: PBrefdist
### Title: Calculate reference distribution using parametric bootstrap
### Aliases: PBrefdist PBrefdist.mer PBrefdist.lmerMod PBrefdist.lm
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
##D cl <- makeSOCKcluster(rep("localhost", 4))
##D clusterEvalQ(cl, library(lme4))
##D clusterSetupSPRNG(cl)
##D rr <- PBrefdist(beet0, beet_no.harv, nsim=20)
##D stopCluster(cl)
## End(Not run)
## Above, 4 cpu's are used and 5 simulations are made on each cpu.



cleanEx()
nameEx("getKR")
### * getKR

flush(stderr()); flush(stdout())

### Name: getKR
### Title: Extract (or "get") components from a 'KRmodcomp' object.
### Aliases: getKR
### Keywords: models

### ** Examples

data(beets, package='pbkrtest')
lg <- lmer(sugpct ~ block + sow + harvest + (1|block:harvest), 
              data=beets, REML=FALSE)
sm <- update(lg, .~. - harvest)
xx<-KRmodcomp(lg, sm)
getKR(xx, "ddf") # get denominator degrees of freedom.




cleanEx()
nameEx("vcovAdj")
### * vcovAdj

flush(stderr()); flush(stdout())

### Name: vcovAdj
### Title: Ajusted covariance matrix for linear mixed models according to
###   Kenward and Roger
### Aliases: vcovAdj LMM_Sigma_G
### Keywords: function

### ** Examples

(fmLarge <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy))
## removing Day
(vcovAdj(fmLarge,detail=0))



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
