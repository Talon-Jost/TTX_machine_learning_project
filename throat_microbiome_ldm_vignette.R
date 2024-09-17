#throat microbiome vignette

library(LDM)
devtools::install_github("yijuanhu/LDM", build_vignettes=TRUE)
library(vegan)


otu_presence = which(colSums(throat.otu.tab>0)>=5)
throat.otu.tab5 = throat.otu.tab[,otu_presence]
dim(throat.otu.tab5)

meta <- throat.meta

fit <- ldm(formula = throat.otu.tab5 | (Sex+AntibioticUse) ~ SmokingStatus + PackYears, 
           data = meta, 
           n.perm.max = 0)

fit$VE.global.freq.submodels
fit$global.tests.stopped
fit$otu.tests.stopped
fit$p.global.omni
fit$detected.otu.omni


w1 <- match(fit$detected.otu.omni[[1]], colnames(fit$q.otu.omni))
order <- order(TTX_ldm$p.otu.omni[1,])[1:10]
order <- as.data.frame(order)