
library(plyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(caret)

rm(list=ls()) 
gc()

setwd('/Users/ben/Library/CloudStorage/GoogleDrive-brharrison2000@gmail.com/My Drive/Documents/DogAgingProject')

#############################################
load('dog metabolome paper/data/blood.covariates.RData') # load CBC covariate names (bloodvars), and the data (b)
cellCovars <- bloodvars[grepl('cbc', bloodvars)] 

load("data/metabolome/ProcessedData/P1.technicalCovarsRemovedData") # the metabolome and dog meta data

head(dat)


rownames(dat) <- dat$dog_id
dat <- dat[b$dog_id, ]
dat <- join(dat, b)

load('data/genetic/grm_Mar_2024.RData')
dat <- dat[dat$dog_id %in% rownames(grm), ] # remove any dogs not in the GRM
grm <- grm[dat$dog_id, dat$dog_id] # align and chop GRM

dim(grm)

# normalize age and weight
dat$sqrtWT <- sqrt(dat$weight_at_DOC)
dat$sqrtAge <- sqrt(dat$age)

## scale all variables
vars <- c('sqrtAge', 'sqrtWT', 'sex', 'sterilization_status', 'hours_fasting', cellCovars)
scDat <- dat

par(mfrow=c(1,2))
boxplot(scDat[ ,vars], pch=20, cex=0.5) # pre-scaling
scDat[ ,vars] <- scale(scDat[ ,vars])
boxplot(scDat[ ,vars], pch=20, cex=0.5) # post-scaling

save(scDat, grm, dogmzs, cellCovars, vars, file='data/scaled.data_for_CBC_mixedModel')

rm(list=ls())
###################################
load('data/scaled.data_for_CBC_mixedModel')

###################################
# fit a mixed model, with fixed effects of age, sex, weight, etc, in the context of random effects of relatedness (as represented in the GRM)

Zmat <- diag(nrow(scDat)) # an empty design matrix for the random effects

# build a design matrix for the fixed effects; specify interaction terms here
vars[!vars%in%dogmzs]

tmp <- ' ~ sqrtAge + sqrtWT + sex * sterilization_status + hours_fasting'
cat(c(tmp, paste0('+ ', cellCovars)))

X <- model.matrix( ~ sqrtAge + sqrtWT + sex * sterilization_status + hours_fasting + krt_cbc_total_wbcs + krt_cbc_abs_bands + krt_cbc_abs_neutrophils + krt_cbc_abs_lymphocytes + krt_cbc_abs_monocytes + krt_cbc_abs_eosinophils + krt_cbc_abs_basophils + krt_cbc_rbc + krt_cbc_hgb + krt_cbc_hct + krt_cbc_mcv + krt_cbc_mch + krt_cbc_mchc + krt_cbc_rdw + krt_cbc_mpv + krt_cbc_pct + krt_cbc_retic_abs, data=scDat) 

mzdat <- scDat[ ,dogmzs]

###################################
# the next steps can be slow, use parallel processing
library(parallel)
library(doParallel)
library(EMMREML)

n.cores <- detectCores(all.tests = F, logical = T) # detects the number of available cores
clus <- makeCluster(n.cores) 
registerDoParallel(cores=n.cores)  

emma <- emmreml(y=mzdat[ ,1], X=X, Z=Zmat, K=grm, varbetahat = T, varuhat=T, PEVuhat=T, test=T) # the model fitting step, run an example to see what it does and to record the output rownames etc. (the code below will call the output)
emma$betahat

clusterExport(clus, varlist=c("mzdat", 'dogmzs', 'X', 'Zmat', 'grm'), envir = environment()) # specifies what things in the R environment that will be accessible to parallel processing function

###################################
# this code fits the same model twice, the first time it extracts the fixed effects, and the second time, it takes the random effects (they are called BLUPs):

# get fixed effects
clus <- makeCluster(n.cores) 
registerDoParallel(cores=n.cores)  
clusterExport(clus, varlist=c("mzdat", 'dogmzs', 'X', 'Zmat', 'grm'), envir = environment()) 

fixedEffects <- t(parApply(clus, mzdat, 2, function(Y) {
  library(EMMREML)
  emma <- emmreml(y=Y, X=X, Z=Zmat, K=grm, varbetahat = T, varuhat=T, PEVuhat=T, test=T) 
  p=emma$pvalbeta
  b=emma$betahat
  return(c(b, p[ ,"none"])) } ))

fixedEffects <- as.data.frame(fixedEffects)

colnames(fixedEffects) <- c(paste0('beta_', rownames(emma$betahat)), paste0('P_', rownames(emma$betahat)))


# get BLUEs (subtract these from y to get the residual Y after the fixed effects):
bhats <- as.matrix(t(fixedEffects[ ,grepl('beta', colnames(fixedEffects))]))
BLUEs <- X %*% bhats

# get random effects (BLUPs)
randomEffects <- parApply(clus, mzdat, 2, function(Y) {
  library(EMMREML)
  emma <- emmreml(y=Y, X=X, Z=Zmat, K=grm, varbetahat = T, varuhat=T, PEVuhat=T, test=T) 
  blup <- emma$uhat
  varblup <- emma$varuhat
  return(list(blup, varblup)) } )

save(file='CBC_CELLcovars_mixed_model_results', fixedEffects, randomEffects, BLUEs, scDat, vars, dogmzs, bhats, cellCovars, X)
#######################################


rm(list=ls())
#######################################
load('CBC_CELLcovars_mixed_model_results')


head(fixedEffects)
fdr <- as.data.frame(apply(fixedEffects[ ,grepl('P_', colnames(fixedEffects))], 2, function(x) p.adjust(x, 'fdr')))
colnames(fdr) <- gsub('P_', 'FDR_', colnames(fdr))

colSums(fdr<=0.05)

plot(fixedEffects$beta_sqrtAge, -log(fixedEffects$P_sqrtAge, 10))


plot(scDat$`N-Ac-Phenylalanine` ~ scDat$sqrtAge)

