# blood chem covariates:
library(plyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(haven)
library(ordinal)

rm(list=ls())

setwd('/Users/ben/Library/CloudStorage/GoogleDrive-brharrison2000@gmail.com/My Drive/Documents/DogAgingProject/dog metabolome paper/')

# there are several types of data about the blood samples, 'CBC' (complete blood count, counts of blood cell types etc.), and 'ChemistryPanel' blood chemistry measurements.   
dir('data/2023_release')

cbc <- read.csv('data/2023_release/DAP_2023_SamplesResults_CBC_v1.0.csv', stringsAsFactors = T) # cbc data as of 2023 release (feb 2024)
colnames(cbc)[colnames(cbc)=='Sample_Year'] <- 'cohort'
cbc$cohort <- gsub('precision_', 'Precision ', cbc$cohort)
head(cbc)
cbc <- cbc[ ,!grepl('comments', colnames(cbc))]# this step omits the "krt_cbc_test_comments"
cbc <- cbc[ ,!grepl('DAP_Sample_ID', colnames(cbc))]


chem <- read.csv('data/2023_release/DAP_2023_SamplesResults_ChemistryPanel_v1.0.csv', stringsAsFactors = T) # cbc data as of 2023 release (feb 2024)
head(chem)
colnames(chem)[colnames(chem)=='Sample_Year'] <- 'cohort'
chem$cohort <- gsub('precision_', 'Precision ', chem$cohort)
head(chem)
chem<- chem[ ,!grepl('comments', colnames(chem))] # this step omits the comments
chem<- chem[ ,!grepl('DAP_Sample_ID', colnames(chem))] 

b <- merge(cbc, chem, by=c('dog_id', 'cohort')) # make b (blood)

bloodvars <- colnames(b)[-c(1:3)]
bloodvars

table(grepl('cbc', bloodvars)) # 45 cbc cell traits, and 42 chem traits


b <- b[b$cohort=="Precision 1", ]

# to simplify, remove invariant bloodvars
invars <- apply(b[ ,bloodvars], 2, function(x) var(as.numeric(x), na.rm=T))
invars <- invars[!is.na(invars)] 
table(invars==0) # 4 invariants
invars <- names(invars)[invars==0] 
bloodvars <- bloodvars[!bloodvars %in% invars] # remove 4 invariants

str(b[ ,bloodvars])
bloodvars <- bloodvars[bloodvars!='krt_cbc_hemoparasites'] # another invariant

colSums(is.na(b[ ,bloodvars])) # some of these could drop out due to NAs

# any covariate with NA won't enable 'adjustment' for its effect, unless the NAs are imputed
table(colSums(is.na(b[ ,bloodvars])) < 10)
# no chance to imput NAs when too common, chop by number of NAS:
bloodvars <- bloodvars[colSums(is.na(b[ ,bloodvars])) < 10] # this step leaves 46 blood covars
table(grepl('cbc', bloodvars))

colSums(is.na(b[ ,bloodvars])) # for the continuous vars, transform, then impute
str(b[ ,bloodvars])

# consider removing zero-inflated covars
par(mfrow=c(2,2))
plot(b$krt_cbc_abs_bands) 
plot(b$krt_cbc_abs_basophils)
plot(b$krt_cbc_rel_bands) # these are discrete (counts), in which zero inflation is fine
plot(b$krt_cbc_rel_basophils) # these are discrete (counts), in which zero inflation is fine

str(b[ ,bloodvars])

## remove all 'rel' measures (I worked with these and they are all redudant to the abs variables, and, due to dependence on other measures, are somewhat ambiguous)
bloodvars <- bloodvars[!grepl('rel', bloodvars)]
bloodvars <- bloodvars[!grepl('_per', bloodvars)]

table(grepl('cbc', bloodvars))


# before imputation:
# take the numerics
bloodnums <- select_if(b[ ,bloodvars], is.numeric)

norms <- apply(bloodnums, 2, shapiro.test) # test for normal dist
stat <- sapply(norms, function(x) x$statistic)
pval <- sapply(norms, function(x) x$p.value)

?shapiro.test

dim(bloodnums)

par(mfrow=c(4,6))
for(i in 1:ncol(bloodnums)){
  hist(bloodnums[ ,i], main=round(stat[i], 3), xlab=colnames(bloodnums)[i], ylab='') }

table(stat<0.96, ifelse(grepl('cbc', colnames(bloodnums)), 'cell', 'chem')) # a shapiro test stat of <0.96 seems a good threshold for log transformation

tmp <- data.frame(colnames(bloodnums), stat)
colnames(tmp)[1] <- c('var')

gsub('krt_', '', tmp$var[tmp$stat<0.96 & grepl('cbc', tmp$var)])


loggedBloodvars <- colnames(bloodnums)[stat<0.96]
loggedBloodvars <- gsub('krt_cbc_', '', loggedBloodvars)
loggedBloodvars <-  gsub('krt_cp_', '', loggedBloodvars)
loggedBloodvars <-  gsub('_value', '', loggedBloodvars)
loggedBloodvars

bloodnums[ ,stat<0.96] <- log(bloodnums[ ,stat<0.96]+1)
str(b[ ,bloodvars])

# after log-normalization
norms <- apply(bloodnums, 2, shapiro.test) # test for normal dist
stat <- sapply(norms, function(x) x$statistic)
pval <- sapply(norms, function(x) x$p.value)

par(mfrow=c(4,6))
for(i in 1:ncol(bloodnums)){
  hist(bloodnums[ ,i], main=round(stat[i], 3), xlab=colnames(bloodnums)[i]) }
# check for stragglers

# straggler
par(mfrow=c(2,2))
hist(b$krt_cp_ggt_value)
hist(log(b$krt_cp_ggt_value))
table(b$krt_cp_ggt_value)
table(log(b$krt_cp_ggt_value)) # ggt has a tail, but is digital

colSums(is.na(bloodnums)) # now for imputation
imp <- as.matrix(bloodnums)
imp <- impute::impute.knn(imp)
imp <- imp$data
b[ ,colnames(imp)] <- imp # replace with log-transformed and imputed data

colSums(is.na(b))
b <- b[ ,colSums(is.na(b)) == 0] # remove the last vars with NAs, these should just be non-numerics
str(b)

bloodvars <- colnames(b)[-c(1,2)] # 45 bloodvars

table(grepl('cbc', bloodvars)) # 17 cell traits, 21 chem traits

save(b, bloodvars, loggedBloodvars, file='data/blood.covariates.RData')



rm(list=ls())
###########################################################
# remove techcnical effects; these blood samples, like those for metabolomics, were 'Shipped to Texas' and so transit time and arrival temperature should be considered
load('data/blood.covariates.RData')
samp <- read.csv('data/P1sampleinfo.csv')

samp[1:4, 1:4]
samp$dog_id <- as.character(samp$dog_id)

str(samp)
hist(samp$elapse_DOC_LOG) # transit time
hist(samp$serum_temp) # arrival temperature

b <- merge(samp[ ,c('dog_id', 'elapse_DOC_LOG', 'serum_temp')], b)

median(b$serum_temp)
range(b$serum_temp)


pList <- list()

for(i in 1:length(bloodvars)){
s <- summary(lm(b[ ,bloodvars[i]] ~ elapse_DOC_LOG + serum_temp, b))
pList[[i]] <- s$coefficients[2:3, 4]}

names(pList) <- bloodvars
p <- as.data.frame(do.call(rbind, pList))

p[order(p$elapse_DOC_LOG), ][1:6, ]

ggplot(b, aes(y=krt_cbc_mchc, x=elapse_DOC_LOG))+
  geom_point()+
  theme_bw(base_size=18)+
  geom_smooth(method='lm', size=0)

ggplot(b, aes(y=krt_cbc_mcv, x=elapse_DOC_LOG))+
  geom_point()+
  theme_bw(base_size=18)

ggplot(b, aes(y=krt_cbc_rel_monocytes, x=elapse_DOC_LOG))+
  geom_point(alpha=0.5)+
  theme_bw(base_size=18)+
  geom_smooth(method='lm', size=0)

mean(b$elapse_DOC_LOG)
median(b$elapse_DOC_LOG)
range(b$elapse_DOC_LOG)

p[order(p$serum_temp), ][1:6, ]

ggplot(b, aes(y=krt_cbc_mpv, x=elapse_DOC_LOG))+
  geom_point(alpha=0.5)+
  theme_bw(base_size=18)

ggplot(b, aes(y=krt_cp_magnesium_value, x=elapse_DOC_LOG))+
  geom_point(alpha=0.5)+
  theme_bw(base_size=18)

ggplot(b, aes(y=krt_cbc_rdw, x=elapse_DOC_LOG))+
  geom_point(alpha=0.5)+
  theme_bw(base_size=18)+
  geom_smooth(method='lm', size=0)

## linear correction seems appropriate (ie, none are particularly non-linear)
# adjust for travel time and arrival temp
for(i in 1:length(bloodvars)){
  b[ ,bloodvars[i]] <- residuals(lm(b[ ,bloodvars[i]] ~ elapse_DOC_LOG + serum_temp, b))}

save(b, bloodvars, file='data/blood.covariates.RData')


rm(list=ls())
###########################################################
load('data/blood.covariates.RData')
dir()

## targeted LCMS data for P1 dogs
load('data/p1_update_March2024')

## GRM: genome-relatedness matrix from PLINK
load('data/grm_Mar_2024.RData')

table(dat$dog_id %in% rownames(grm))
g <- dat[dat$dog_id %in% rownames(grm), ]
grm <- grm[b$dog_id, b$dog_id] # subset and order the GRM by the dogs in the mz data

b <- merge(g, b[b$cohort == 'Precision 1', ], by='dog_id', all.y=F, all.x=T)

colSums(is.na(b[ ,bloodvars])) # 34 dogs with missing cbc data
rowSums(is.na(b[ ,bloodvars])) # 
b <- b[rowSums(is.na(b[ ,bloodvars]))==0, ] # omit them

grm <- grm[b$dog_id, b$dog_id] # subset and order the GRM by the dogs in the mz data
dim(grm)

save(b, grm, dogmzs, bloodvars, file='data/dataforCBC_ChemMixedModel.RData')


rm(list=ls())
######################################################################
# fit a mixed model for effects on CBC
######################################################################
load('data/dataforCBC_ChemMixedModel.RData')

# get fixed effects
library(parallel)
library(doParallel)
library(EMMREML)

n.cores <- detectCores(all.tests = F, logical = T) # detects the number of available cores
clus <- makeCluster(n.cores) 
registerDoParallel(cores=n.cores)  

b$sqrtAge <- sqrt(b$age) # normalize age by taking its square root
b$sqrtWT <- sqrt(b$weight_at_DOC) # normalize weight by its square root

Zmat <- diag(nrow(b)) # an empty design matrix for the random effects

vars<- c('sqrtAge', 'sqrtWT', 'hours_fasting', 'sex', 'sterilization_status', bloodvars, dogmzs)

## scale the variables:
scDat <- b
scDat[ ,vars] <- scale(scDat[ ,vars])

# run, just to get emma ##############
X <- model.matrix( ~ b$sqrtAge + b$sqrtWT + b$hours_fasting + b$sex * b$sterilization_status)
emma <- emmreml(y=scDat$Adenosine, X=X, Z=Zmat, K=grm, varbetahat = T, varuhat=T, PEVuhat=T, test=T) 
######################################

cbcdat <- scDat[ ,bloodvars]

X <- model.matrix( ~ b$sqrtAge + b$sqrtWT + b$hours_fasting + b$sex * b$sterilization_status)
clusterExport(clus, varlist=c('cbcdat', 'X', 'Zmat', 'grm', 'emma'), envir = environment()) 
  
tmp <- t(parApply(clus, cbcdat, 2, function(Y) {
    library(EMMREML)
    emma <- emmreml(y=Y, X=X, Z=Zmat, K=grm, varbetahat = T, varuhat=T, PEVuhat=T, test=T) 
    p=emma$pvalbeta
    b=emma$betahat
    return(c(b, p[ ,"none"])) } ))
  
ageOnCBCmm <- as.data.frame(tmp)
colnames(ageOnCBCmm) <- c(paste0('beta_', gsub('b\\$', '', colnames(X))), paste0('P_', gsub('b\\$', '', colnames(X))))

save(file='ageOnCBCmm', ageOnCBCmm, bloodvars)
######################################################################

rm(list=ls())
######################################################################
load('ageOnCBCmm')
head(ageOnCBCmm )

fdrs <- apply(ageOnCBCmm[ ,grepl('P_', colnames(ageOnCBCmm))], 2, function(x) p.adjust(x, 'fdr'))
colnames(fdrs) <- gsub('P_', 'FDR_', colnames(fdrs))
colSums(fdrs<=0.05)







rm(list=ls())
######################################################################
# fit a mixed model for effects of MZ on CBC
######################################################################
load('data/dataforCBC_ChemMixedModel.RData')

# get fixed effects
library(parallel)
library(doParallel)
library(EMMREML)

n.cores <- detectCores(all.tests = F, logical = T) # detects the number of available cores
clus <- makeCluster(n.cores) 
registerDoParallel(cores=n.cores)  

b$sqrtAge <- sqrt(b$age) # normalize age by taking its square root
b$sqrtWT <- sqrt(b$weight_at_DOC) # normalize weight by its square root

Zmat <- diag(nrow(b)) # an empty design matrix for the random effects

vars<- c('sqrtAge', 'sqrtWT', 'hours_fasting', 'sex', 'sterilization_status', bloodvars, dogmzs)

## scale the variables:
scDat <- b
scDat[ ,vars] <- scale(scDat[ ,vars])

# run, just to get emma ##############
X <- model.matrix( ~ b$sqrtAge + b$sqrtWT + b$hours_fasting + b$sex * b$sterilization_status)
emma <- emmreml(y=scDat$Adenosine, X=X, Z=Zmat, K=grm, varbetahat = T, varuhat=T, PEVuhat=T, test=T) 
######################################

CBCbyMZmmList <- list()
cbcdat <- scDat[ ,bloodvars]

for(k in 1:length(dogmzs)) {
  X <- model.matrix( ~ b[ ,dogmzs[k]] + b$sqrtAge + b$sqrtWT + b$hours_fasting + b$sex * b$sterilization_status)
  clusterExport(clus, varlist=c('cbcdat', 'dogmzs', 'X', 'Zmat', 'grm', 'emma'), envir = environment()) 
  
  tmp <- t(parApply(clus, cbcdat, 2, function(Y) {
    library(EMMREML)
    emma <- emmreml(y=Y, X=X, Z=Zmat, K=grm, varbetahat = T, varuhat=T, PEVuhat=T, test=T) 
    p=emma$pvalbeta
    b=emma$betahat
    return(c(b, p[ ,"none"])) } ))

  tmp <- as.data.frame(tmp)
  colnames(X)[grepl('blood', colnames(X))] <- dogmzs[k]
  colnames(tmp) <- c(paste0('beta_', gsub('b\\$', '', colnames(X))), paste0('P_', gsub('b\\$', '', colnames(X))))
  CBCbyMZmmList[[k]] <- tmp  }

names(CBCbyMZmmList) <- dogmzs


save(file='CBCbyMZmmList', CBCbyMZmmList, dogmzs, bloodvars)
######################################################################


rm(list=ls())
######################################################################
load('CBCbyMZmmList')

colnames(CBCbyMZmmList[[1]])

betas <- lapply(CBCbyMZmmList, function(x) x[ ,grepl('beta', colnames(x))])
betas <- do.call(cbind, betas)
betas <- betas[ ,grepl('dogmzs', colnames(betas))]
colnames(betas) <- dogmzs
betas

p <- lapply(CBCbyMZmmList, function(x) x[ ,grepl('P', colnames(x))])
p <- do.call(cbind,p)
p <- p[ ,grepl('dogmzs', colnames(p))]
colnames(p) <- dogmzs

fdr <- apply(p, 1, function(x) p.adjust(x, 'fdr'))
colSums(fdr <=0.05)
rowSums(fdr <=0.05)

df <- as.data.frame(betas)
df$bloodvar <- bloodvars
long <- pivot_longer(df, -bloodvar, values_to='beta', names_to = 'mz')
head(long)

df2 <- as.data.frame(fdr)
df2$mz <- dogmzs
long2 <- pivot_longer(df2, -mz, values_to='FDR', names_to = 'bloodvar')
long <- merge(long, long2)
head(long)
rm(long2)



rm(list=ls())
######################################################################
# with MZ as outcome
######################################################################
load('data/dataforCBC_ChemMixedModel.RData')

# get fixed effects
library(parallel)
library(doParallel)
library(EMMREML)

n.cores <- detectCores(all.tests = F, logical = T) # detects the number of available cores
clus <- makeCluster(n.cores) 
registerDoParallel(cores=n.cores)  

b$sqrtAge <- sqrt(b$age) # normalize age by taking its square root
b$sqrtWT <- sqrt(b$weight_at_DOC) # normalize weight by its square root

Zmat <- diag(nrow(b)) # an empty design matrix for the random effects

vars<- c('sqrtAge', 'sqrtWT', 'hours_fasting', 'sex', 'sterilization_status', bloodvars, dogmzs)

## scale the variables:
scDat <- b
scDat[ ,vars] <- scale(scDat[ ,vars])

# run, just to get emma ##############
X <- model.matrix( ~ b$sqrtAge + b$sqrtWT + b$hours_fasting + b$sex * b$sterilization_status)
emma <- emmreml(y=scDat$Adenosine, X=X, Z=Zmat, K=grm, varbetahat = T, varuhat=T, PEVuhat=T, test=T)  
#########################################

bloodEffectsList <- list()

mzdat <- scDat[ ,dogmzs]

for(k in 1:length(bloodvars)) {
  cbcX <- model.matrix( ~ b[ ,bloodvars[k]] + b$sqrtAge + b$sqrtWT + b$hours_fasting + b$sex * b$sterilization_status)
  clusterExport(clus, varlist=c('mzdat', 'dogmzs', 'cbcX', 'Zmat', 'grm', 'emma'), envir = environment()) 
  
  tmp <- t(parApply(clus, mzdat, 2, function(Y) {
    library(EMMREML)
    emma <- emmreml(y=Y, X=cbcX, Z=Zmat, K=grm, varbetahat = T, varuhat=T, PEVuhat=T, test=T) 
    p=emma$pvalbeta
    b=emma$betahat
    return(c(b, p[ ,"none"])) } ))
  

tmp <- as.data.frame(tmp)
colnames(cbcX)[grepl('blood', colnames(cbcX))] <- bloodvars[k]
colnames(tmp) <- c(paste0('beta_', gsub('b\\$', '', colnames(cbcX))), paste0('P_', gsub('b\\$', '', colnames(cbcX))))
bloodEffectsList[[k]] <- tmp  }

names(bloodEffectsList) <- bloodvars

save(file='CBCbyMZmmList', bloodEffectsList)


#####################################################
rm(list=ls())
load('mixed_model_results')
load('bloodEffectsList')

e <- bloodEffectsList
betas <- lapply(e, function(x) x[ ,grepl('beta_krt', colnames(x))])
betas <- do.call(cbind, betas)
rownames(betas) <- dogmzs
betas

p <- lapply(e, function(x) x[ ,grepl('P_krt', colnames(x))])
p <- do.call(cbind,p)
rownames(p) <- dogmzs

fdr <- apply(p, 2, function(x) p.adjust(x, 'fdr'))
colSums(fdr <=0.05)
rowSums(fdr <=0.05)
rm(e)

hits <- betas[rowSums(fdr<=0.05)>0, ]
pheatmap(t(hits))

df <- as.data.frame(betas)
df$mz <- dogmzs
long <- pivot_longer(df, -mz, values_to='beta', names_to = 'bloodvar')
head(long)

df2 <- as.data.frame(fdr)
df2$mz <- dogmzs
long2 <- pivot_longer(df2, -mz, values_to='FDR', names_to = 'bloodvar')
long <- merge(long, long2)
head(long)

rm(long2)


#######################################################################
# [for plotting the cbc/chem effects] remove covariate fixed effects, and BLUPs, :
load('BLUPsubtracted.RData') 
effectsSub <- matrix(nr=nrow(g), nc=length(dogmzs))

for(i in 1:length(dogmzs)) {
  effectsSub[ ,i] <- BLUPsubtracted[ ,dogmzs[i]] - BLUEs[ ,i] }
colnames(effectsSub) <- dogmzs
rownames(effectsSub) <- g$dog_id

load('data/dataforCBC_ChemMixedModel.RData')
effectsSub <- effectsSub[b$dog_id, ] # not all dog shave cbc chem data

b[ ,dogmzs] <- effectsSub 

b$sqrtAge <- sqrt(b$age)
b$sqrtWt <- sqrt(b$weight_at_DOC)

long[order(long$FDR), ][1:14, ]

dev.off()

ggplot(b, aes(y= Creatinine, x=krt_cp_creatinine_value))+
  geom_point(size=0.5)+
  geom_smooth(method='lm', se=F, color='purple', size=0.5)+
  theme_classic()

ggplot(b, aes(y= Creatinine, x=sqrtAge))+
  geom_point(size=0.5)+
  geom_smooth(method='lm', se=F, color='purple', size=0.5)+
  theme_classic()


p1 <- ggplot(b, aes(y= Allantoin, x=krt_cp_creatinine_value))+
  geom_point(size=0.5)+
  geom_smooth(method='lm', se=F, color='purple', size=0.5)+
  theme_classic()

p2 <- ggplot(b, aes(y= Allantoin, x= krt_cp_bun_value))+
  geom_point(size=0.5)+
  geom_smooth(method='lm', se=F, color='purple', size=0.5)+
  theme_classic()

p3 <- ggplot(b, aes(x= krt_cp_bun_value, y=krt_cp_creatinine_value))+
  geom_point(size=0.5)+
  geom_smooth(method='lm', se=F, color='purple', size=0.5)+
  theme_classic()


resAllantoin_Creat <- residuals(lm(Allantoin  ~ krt_cp_creatinine_value, b))

summary(lm(Allantoin ~ krt_cp_creatinine_value * krt_cp_bun_value, b))


p4 <- ggplot(b, aes(y= resAllantoin_Creat, x=krt_cp_bun_value))+
  geom_point(size=0.5)+
  ylab('residual allantoin ~ creatinine')+
  geom_smooth(method='lm', se=F, color='purple', size=0.5)+
  theme_classic()

ggarrange(p1, p2, p3, p4)





