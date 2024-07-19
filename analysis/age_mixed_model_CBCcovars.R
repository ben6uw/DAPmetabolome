
library(plyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(caret)

rm(list=ls()) 
gc()

setwd('/Users/ben/Library/CloudStorage/GoogleDrive-brharrison2000@gmail.com/My Drive/Documents/DogAgingProject/dog metabolome paper/')

#############################################
load('data/data_for_CBC_mixedModel')

# which CBC measures should be covariates in the mz ~ age modeling?
cbcCovars <- read.csv('cbcCovars.csv') # from a mediation analysis, the nubmer of either CBC~age mediated by mzs, or the number of mz~age mdeiated by CBCs.  The later are the CBCs we should add as coviarates to the mz ~ age model
head(cbcCovars)
CBC <- gsub('_value', '', cbcCovars$X) #simplify labels for the CBC variables
CBC <- gsub('krt_cbc_', '', CBC)
CBC <- gsub('krt_cp_', '', CBC)
cbcCovars$CBC <- CBC 
rm(CBC)

cbcCovars$X

cbcCovars$bloodTraitType <- factor(ifelse(grepl('cbc', cbcCovars$X), 'cell', 'chem'))
cbcCovars$bloodTraitType

ggplot(cbcCovars, aes(y=bloodvar..MZ..age, x=MZ..bloodvar..age, label=CBC, color=bloodTraitType))+
  geom_point()+
  ylab('CBC ~ [MZ] x age (how often a CBC trait is mediated)') +
  xlab('MZ ~ [CBC] x age (number of metabolites mediated)') +
  theme_bw(base_size = 16) +
  ggtitle('outcome of mediation analysis')+
  ggrepel::geom_text_repel()+
  scale_color_manual(values=c('purple', 'darkgreen'))+
  geom_abline(intercept = 0, slope=1, color='grey40')


table(cbcCovars$MZ..bloodvar..age)

head(cbcCovars)
cellCovars <- cbcCovars$X[grepl('cbc', cbcCovars$X)] # one option, is to consider 'CBC' blood cell counts as covars, and not the blood chemistry traits as covars (i.e. why correct for creatinine, when it is one of the mzs on the pannel?)

cbcMediatorCovs <- cbcCovars[cbcCovars$MZ..bloodvar..age > cbcCovars$bloodvar..MZ..age, ] # select cbc variables, may increase stringency to avoid convergence problems in the mixed model 
cbcMediatorCovs  <- cbcMediatorCovs$X

# normalize age and weight
b$sqrtWT <- sqrt(b$weight_at_DOC)
b$sqrtAge <- sqrt(b$age)

## scale all variables
vars <- c('sqrtAge', 'sqrtWT', 'sex', 'sterilization_status', 'hours_fasting', dogmzs, cellCovars)
scDat <- b

par(mfrow=c(1,2))
boxplot(scDat[ ,vars], pch=20, cex=0.5) # pre-scaling
scDat[ ,vars] <- scale(scDat[ ,vars])
boxplot(scDat[ ,vars], pch=20, cex=0.5) # post-scaling

save(scDat, grm, dogmzs, cellCovars, cbcMediatorCovs, vars, file='data/scaled.data_for_CBC_mixedModel')

rm(list=ls())
###################################
load('data/scaled.data_for_CBC_mixedModel')

table(scDat$dog_id %in% rownames(grm))
g <- scDat[scDat$dog_id %in% rownames(grm), ] # g is now composed of variables that have been scaled
grm <- grm[g$dog_id, g$dog_id] # subset and order the GRM by the dogs in the mz data

###################################
# fit a mixed model, with fixed effects of age, sex, weight, etc, in the context of random effects of relatedness (as represented in the GRM)

Zmat <- diag(nrow(g)) # an empty design matrix for the random effects

# build a design matrix for the fixed effects; specify interaction terms here
vars[!vars%in%dogmzs]

tmp <- ' ~ g$sqrtAge + g$sqrtWT + g$sex * g$sterilization_status + g$hours_fasting'
cat(c(tmp, paste0('+ g$', cellCovars)))

X <- model.matrix( ~ g$sqrtAge + g$sqrtWT + g$sex * g$sterilization_status + g$hours_fasting + g$krt_cbc_abs_bands + g$krt_cbc_abs_eosinophils + g$krt_cbc_abs_lymphocytes + g$krt_cbc_abs_monocytes + g$krt_cbc_abs_neutrophils + g$krt_cbc_hct + g$krt_cbc_hgb + g$krt_cbc_mch + g$krt_cbc_mcv + g$krt_cbc_mpv + g$krt_cbc_pct + g$krt_cbc_rbc + g$krt_cbc_rdw + g$krt_cbc_rel_bands + g$krt_cbc_rel_lymphocytes + g$krt_cbc_rel_monocytes + g$krt_cbc_rel_neutrophils + g$krt_cbc_retic_abs + g$krt_cbc_retic_per + g$krt_cbc_total_wbcs) 

mzdat <- g[ ,dogmzs]

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

colnames(fixedEffects) <- c(paste0('beta_', gsub('g\\$', '', rownames(emma$betahat))), paste0('P_', gsub('g\\$', '', rownames(emma$betahat))))

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


save(file='CBC_CELLcovars_mixed_model_results', fixedEffects, randomEffects, BLUEs, g, vars, dogmzs, bhats, cellCovars, cbcMediatorCovs, X)
#######################################


rm(list=ls())
#######################################
load('CBC_CELLcovars_mixed_model_results')


head(fixedEffects)
fdr <- as.data.frame(apply(fixedEffects[ ,grepl('P_', colnames(fixedEffects))], 2, function(x) p.adjust(x, 'fdr')))
colnames(fdr) <- gsub('P_', 'FDR_', colnames(fdr))

colSums(fdr<=0.05)


plot(fixedEffects$beta_sqrtAge, -log(fixedEffects$P_sqrtAge, 10))

head(fixedEffects)
mm <- data.frame('cbcCovars_ageBeta' =fixedEffects$beta_sqrtAge, 'cbcCovars_ageFDR' =fdr$FDR_sqrtAge)

load('mixed_model_results')
original.fdr <- as.data.frame(apply(fixedEffects[ ,grepl('P_', colnames(fixedEffects))], 2, function(x) p.adjust(x, 'fdr')))
colnames(original.fdr) <- gsub('P_', 'FDR_', colnames(original.fdr))
colSums(original.fdr<=0.05)

mm$original_ageBeta = fixedEffects$`beta_scale(sqrtAge)`
mm$original_ageFDR = original.fdr$`FDR_scale(sqrtAge)`

head(mm)
table('original'=ifelse(mm$original_ageFDR<=0.05, '5%FDR', 'NS'), 'cell covars model'=ifelse(mm$cbcCovars_ageFDR<=0.05, '5%FDR', 'NS'))

mm$mz <- dogmzs

ggplot(mm, aes(y=cbcCovars_ageBeta, x=original_ageBeta, label=mz))+
  geom_abline(intercept = 0, slope=1, color='grey50')+
  geom_point()+
  theme_bw(base_size = 16) +
  ylab('cell trait covariates model')+
  ggtitle('age effects with CELL-ONLY covariates')+
  ggrepel::geom_text_repel()


mm[mm$cbcCovars_ageFDR<=0.05 & mm$original_ageFDR<=0.05, ]




