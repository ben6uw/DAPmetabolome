
library(plyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(caret)

rm(list=ls()) 
gc()

setwd('/Users/ben/Library/CloudStorage/GoogleDrive-ben6@uw.edu/My Drive/Documents/DogAgingProject/dog metabolome paper/')


#############################################
### merge mz data and dog meta data

## targeted LCMS data for P1 dogs
dat<- readRDS('data/P1metabolites.RDS')
dat <- as.data.frame(t(dat))
dogmzs <- colnames(dat) # grab the names of the metabolites
dat$dog_id <- rownames(dat)

## sample info
samps <- read.csv('data/P1sampleinfo.csv')
samps$dog_id <- as.character(samps$dog_id)
head(samps)
table(samps$dog_id %in% dat$dog_id)
dat <- merge(samps, dat, by='dog_id')


load('data/2023_release/DAP_2023_DogOverview_v1.0.RData') # dog meta data; grab estimated date of birth (DOB), and the categorical regarding the nature of the estimate
DogOverview <- haven::zap_labels(DogOverview)
DogOverview$dog_id <-  as.character(DogOverview$dog_id)

dat <- merge(DogOverview[ ,c('dog_id', 'Estimated_DOB', 'AgeDOB_Estimation_Method')], dat, all.x=F)

prior.age <- dat$age

dat$age <- (as.numeric(as.Date(dat$DOC) - as.Date(dat$Estimated_DOB)))/365

ggplot(dat, aes(y=age, x=prior.age, color=AgeDOB_Estimation_Method))+
  geom_abline(intercept=0, slope=1, color='grey')+
  geom_point()+
  theme_classic(base_size = 14)+
  theme(legend.position = 'none')+
  facet_wrap(~ AgeDOB_Estimation_Method)


save(file='data/p1_update_March2024', dat, dogmzs)



rm(list=ls())
###################################
load('data/p1_update_March2024.Rdata')

## GRM: genome-relatedness matrix from PLINK
load('data/grm_Mar_2024.RData')

table(dat$dog_id %in% rownames(grm))
g <- dat[dat$dog_id %in% rownames(grm), ]
grm <- grm[g$dog_id, g$dog_id] # subset and order the GRM by the dogs in the mz data

###################################
# fit a mixed model, with fixed effects of age, sex, weight, etc, in the context of random effects of relatedness (as represented in the GRM)

g$sqrtAge <- sqrt(g$age) # normalize age by taking its square root
g$sqrtWT <- sqrt(g$weight_at_DOC) # normalize weight by its square root

Zmat <- diag(nrow(g)) # an empty design matrix for the random effects
X <- model.matrix( ~  scale(g$sqrtAge) * scale(g$sqrtWT) + scale(g$hours_fasting) + scale(g$sex) * scale(g$sterilization_status)) # builds a design matrix for the fixed effects; specify interaction terms here

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
colnames(fixedEffects) <- gsub('scale\\(', 'sc_', colnames(fixedEffects))
colnames(fixedEffects) <- gsub('\\)$', '', colnames(fixedEffects))
colnames(fixedEffects) <- gsub('\\)', '', colnames(fixedEffects))
colnames(fixedEffects) <- gsub('\\(', '', colnames(fixedEffects))

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

save(file='mixed_model_results', fixedEffects, randomEffects, BLUEs, g, dogmzs, bhats, X)
#######################################






