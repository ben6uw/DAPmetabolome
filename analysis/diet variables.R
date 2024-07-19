
rm(list=ls())

# diet info from the 2023 release:
load('data/p1_update_March2024')
dat$dog_id

dir('data/2023_release')
cb <-read.csv('data/2023_release/DAP_2023_CODEBOOK_v1.0.csv')
head(cb)

load('data/2023_release/DAP_2023_HLES_dog_owner_v1.0.RData')
h <- HLES_dog_owner
rm(HLES_dog_owner)
h <- h[h$dog_id%in%dat$dog_id, ]

colnames(h)[grepl('diet', colnames(h))]

question <- cb$SurveyText[cb$Variable=='df_primary_diet_component']
primarydiet <- as.factor(h$df_primary_diet_component)
vals <- cb$ValueLabels[cb$Variable=='df_primary_diet_component']
vals <- gsub(' ', '', unlist(strsplit(vals, split=' \\| ')))

levels(primarydiet) <- vals

table(primarydiet)

diet <- data.frame('dog_id'=h$dog_id, primarydiet)
head(diet)
str(diet)


# to test for assoc, btw residual mz and diet, add diet to the fill mixed model:

###################################
load('data/scaled.data_for_CBC_mixedModel')

scDat <- merge(diet, scDat, all.x=F)
scDat <- scDat[!is.na(scDat$primarydiet), ]

scDat$source <- 'Other'
scDat$source[grepl('Commercial', scDat$primarydiet)] = 'Commercial'
scDat$source[grepl('Home', scDat$primarydiet)] = 'Home'
scDat$source <- as.factor(scDat$source)
table(scDat$source)

table(scDat$primarydiet)

scDat$dietnames <- scDat$primarydiet
scDat$dietnames <- gsub('Commerciallyprepared' ,'', scDat$dietnames)
scDat$dietnames <- gsub('Homeprepared' ,'', scDat$dietnames)
scDat$dietnames <- gsub('food' ,'', scDat$dietnames)
scDat$dietnames <- gsub('diet' ,'', scDat$dietnames)

table(scDat$dietnames)
scDat$dietnames <- factor(scDat$dietnames, levels=c("dry(kibble)", "refrigeratedorfrozenraw", "canned", "freeze-dried", "semi-dryorsemi-moist", "cooked", "raw", "Other"))

table(scDat$dog_id %in% rownames(grm))
grm <- grm[scDat$dog_id, scDat$dog_id] # subset and order the GRM by the dogs in the mz data

###################################
# fit a mixed model, with fixed effects of age, sex, weight, etc, in the context of random effects of relatedness (as represented in the GRM)

Zmat <- diag(nrow(scDat)) # an empty design matrix for the random effects

# build a design matrix for the fixed effects; specify interaction terms here

vars[!vars%in%dogmzs]

tmp <- ' ~ primarydiet + sqrtAge + sqrtWT + sex * sterilization_status + hours_fasting'
cat(c(tmp, paste('+', cellCovars)))

## NOTE: becuase diet is a factor the first level (set to kibble, the most popular) as the ref, all effects compare to it!
X <- model.matrix(  ~ dietnames + sqrtAge + sqrtWT + sex * sterilization_status + hours_fasting + krt_cbc_abs_bands + krt_cbc_abs_eosinophils + krt_cbc_abs_lymphocytes + krt_cbc_abs_monocytes + krt_cbc_abs_neutrophils + krt_cbc_hct + krt_cbc_hgb + krt_cbc_mch + krt_cbc_mcv + krt_cbc_mpv + krt_cbc_pct + krt_cbc_rbc + krt_cbc_rdw + krt_cbc_rel_bands + krt_cbc_rel_lymphocytes + krt_cbc_rel_monocytes + krt_cbc_rel_neutrophils + krt_cbc_retic_abs + krt_cbc_retic_per + krt_cbc_total_wbcs, data=scDat) 
X

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

betas <- fixedEffects[ ,1:length(rownames(emma$betaha))]
colnames(betas) <- gsub('dietnames', '', rownames(emma$betahat) )

Ps <- fixedEffects[ ,c((1+length(rownames(emma$betaha))):ncol(fixedEffects))]
colnames(Ps) <- gsub('dietnames', '', rownames(emma$betahat) )

fdr <- as.data.frame(apply(Ps, 2, function(x) p.adjust(x, 'fdr')))
colSums(fdr<=0.05) 

# get random effects (BLUPs)
randomEffects<- parApply(clus, mzdat, 2, function(Y) {
  library(EMMREML)
  emma <- emmreml(y=Y, X=X, Z=Zmat, K=grm, varbetahat = T, varuhat=T, PEVuhat=T, test=T) 
  blup <- emma$uhat
  varblup <- emma$varuhat
  return(list(blup, varblup)) } )


save(file='diet_mixed_model_results', fixedEffects, randomEffects, scDat, dogmzs, X, fdr, betas, Ps)




rm(list=ls())
##########################################################################
load('diet_mixed_model_results')

mp <- read.csv('modAAs.csv') # read the ptmAAs
head(mp)

colSums(fdr<=0.05)
dogmzs[fdr$canned <=0.05] %in% mp$mz[mp$type=='PTM']
dogmzs[fdr$refrigeratedorfrozenraw <=0.05] %in% mp$mz[mp$type=='PTM']
dogmzs[fdr$cooked <=0.05]%in% mp$mz[mp$type=='PTM']
dogmzs[fdr$raw <=0.05]%in% mp$mz[mp$type=='PTM']

dogmzs[fdr$refrigeratedorfrozenraw <=0.05] 
mp

## make a table for paper:
fdr[ ,dietnames]

levels(scDat$dietnames)
types <- levels(scDat$dietnames)
types[6:7] <- paste('home-made,', types[6:7])

dietFDR <- fdr[ ,c("(Intercept)", levels(scDat$dietnames)[-c(1)])]
colnames(dietFDR) <- levels(scDat$dietnames)
dietBetas <- betas[ ,c("(Intercept)", levels(scDat$dietnames)[-c(1)])]
colnames(dietBetas) <- types
dietFDRmat <- as.matrix(dietFDR[rowSums(dietFDR <=0.05)>0, colSums(dietFDR <=0.05)>0])
dietBetasMat <- as.matrix(dietBetas[rowSums(dietFDR <=0.05)>0, colSums(dietFDR <=0.05)>0])
pheatmap(ifelse(dietFDRmat <=0.05, 1, 0))
pheatmap(dietBetasMat)

ifelse(dietFDRmat <=0.05, 1, 0)


test_labels <-dietFDRmat
test_labels[dietFDRmat >0.05] <- ''
test_labels[dietFDRmat <=0.05] <- "*"
test_labels[ ,1] <- ''
pheatmap(dietBetasMat, display_numbers = test_labels, fontsize_number=14)



## I should take model residuals for plotting:

# get BLUEs (subtract these from y to get the residual Y after the fixed effects):
bhats <- as.matrix(t(betas))
rownames(bhats)

subbhats <- bhats[!rownames(bhats) %in% unique(scDat$dietnames), ] # remove variable of interest, [note]: as the refernce level, teh mean(main) effect of kiddle is 0, so its not among the BLUEs, the other diets simply modify from 0
subX <- X[ ,!colnames(X) %in% paste0('dietnames', unique(scDat$dietnames))]
BLUEs <- subX %*% subbhats
BLUEs

DietResiduals <- matrix(nr=nrow(scDat), nc=length(dogmzs))

for(i in 1:length(dogmzs)) {
  BLUPS <- unlist(randomEffects[[dogmzs[i]]][1]) # the BLUPS for a given mz
  DietResiduals [ ,i] <- scDat[ ,dogmzs[i]] - BLUEs[ ,i] - BLUPS }
colnames(DietResiduals) <- dogmzs

tmp <- scDat
tmp[ ,dogmzs] <- DietResiduals
DietResiduals <- tmp



ggplot(scDat, aes(y=`1/3-Methylhistidine`, x=dietnames, fill=source, color=source))+
  geom_boxplot(alpha=0.5)+
  geom_jitter(width=0.1)+
  coord_flip()+
  theme_classic(base_size = 16)+
  xlab('')+
  scale_x_discrete(limits=rev)

ggplot(scDat, aes(y=`N-Acetyl-Aspartate (NAA)`, x=dietnames, fill=source, color=source))+
  geom_boxplot(alpha=0.5)+
  geom_jitter(width=0.1)+
  coord_flip()+
  theme_classic(base_size = 16)+
  xlab('')+
  scale_x_discrete(limits=rev)

library(ggsignif)

p1 <-ggplot(scDat, aes(x=dietnames, fill=source, color=source))+
  geom_bar(stat='count')+  
  theme_bw(base_size = 14)+
  xlab('')+
  ylab('dogs (log scale)')+
  ggtitle('primary diet')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_y_log10() 


p2 <- ggplot(scDat, aes(y=`N-Acetyl-Aspartate (NAA)`, x=dietnames, fill=source, color=source))+
  geom_boxplot(alpha=0.5, outlier.shape=NA)+
  geom_jitter(width=0.1)+
  theme_bw(base_size = 14)+
  xlab('')+
  ylab('N-Ac-aspartate')+
  ylim(c(-3, 3.5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = 'none')+
  geom_signif(comparisons = list(c("dry(kibble)", "refrigeratedorfrozenraw")), 
              map_signif_level=T, color=1)

p3 <- ggplot(scDat, aes(y=`N6-Acetyl-Lysine`, x=dietnames, fill=source, color=source))+
  geom_boxplot(alpha=0.5, outlier.shape=NA)+
  geom_jitter(width=0.1)+
  theme_bw(base_size = 14)+
  xlab('')+
  ylab('N6-Ac-lysine')+
  ylim(c(-3.5, 4)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = 'none')+
  geom_signif(comparisons = list(c("dry(kibble)", "refrigeratedorfrozenraw")), 
              map_signif_level=T, color=1)

ggarrange(p1, ggarrange(p2, p3, ncol=1), nrow=1, widths=c(1, 0.7))


ggplot(scDat, aes(y=`Creatine`, x=dietnames, fill=source, color=source))+
  geom_boxplot(alpha=0.5)+
  geom_jitter(width=0.1)+
  theme_classic(base_size = 14)+
  xlab('')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))




table('age'=fdr$sqrtAge<=0.05, 'diet'=fdr$canned  <=0.05)
table('age'=fdr$sqrtAge<=0.05, 'diet'=fdr$refrigeratedorfrozenraw <=0.05)
table('age'=fdr$sqrtAge<=0.05, 'diet'=fdr$cooked <=0.05)
table('age'=fdr$sqrtAge<=0.05, 'diet'=fdr$raw <=0.05)

# plot the betas of diet by the betas of age...

plot(betas$cooked , betas$sqrtAge, pch=19)

dietCor <- cor(betas[ ,2:10])
diag(dietCor) = NA

pheatmap(dietCor)

plot(betas$cooked , betas$raw, pch=19)
plot(betas$refrigeratedorfrozenraw , betas$raw, pch=19)
pairs(betas[ ,2:10], pch=20)
## the betas among the food types are correlated, do not analyze beta assoc by food.



################################################################################
#### as suggested by Katie Tolbert, Vet nutritionist: to look withing a grouo with more conisistent diet, run models in only the kibble dogs
# also consider omitting those kibble dogs on vetereary therapy food.
################################################################################


###################################
# fit a mixed model, with fixed effects of age, sex, weight, etc, in the context of random effects of relatedness (as represented in the GRM)
table(scDat$dietnames)
table(scDat$dietnames == 'dry(kibble)')
kibblers <- scDat[scDat$dietnames == 'dry(kibble)', ]  # 653 of these
 
Zmat <- diag(nrow(kibblers)) # an empty design matrix for the random effects
grm <- grm[kibblers$dog_id, kibblers$dog_id]
# build a design matrix for the fixed effects; specify interaction terms here

vars[!vars%in%dogmzs]

tmp <- ' ~ sqrtAge + sqrtWT + sex * sterilization_status + hours_fasting'
cat(c(tmp, paste('+', cellCovars)))

## NOTE: becuase diet is a factor the first level (set to kibble, the most popular) as the ref, all effects compare to it!
X <- model.matrix(  ~ sqrtAge + sqrtWT + sex * sterilization_status + hours_fasting + krt_cbc_abs_bands + krt_cbc_abs_eosinophils + krt_cbc_abs_lymphocytes + krt_cbc_abs_monocytes + krt_cbc_abs_neutrophils + krt_cbc_hct + krt_cbc_hgb + krt_cbc_mch + krt_cbc_mcv + krt_cbc_mpv + krt_cbc_pct + krt_cbc_rbc + krt_cbc_rdw + krt_cbc_rel_bands + krt_cbc_rel_lymphocytes + krt_cbc_rel_monocytes + krt_cbc_rel_neutrophils + krt_cbc_retic_abs + krt_cbc_retic_per + krt_cbc_total_wbcs, data=kibblers) 
X

mzdat <- kibblers[ ,dogmzs]

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

betas <- fixedEffects[ ,1:length(rownames(emma$betaha))]
colnames(betas) <- gsub('dietnames', '', rownames(emma$betahat) )

Ps <- fixedEffects[ ,c((1+length(rownames(emma$betaha))):ncol(fixedEffects))]
colnames(Ps) <- gsub('dietnames', '', rownames(emma$betahat) )

fdr <- as.data.frame(apply(Ps, 2, function(x) p.adjust(x, 'fdr')))
colSums(fdr<=0.05)

mods <- read.csv('modAAs.csv')
ptm <- mods$mz[mods$type=='PTM']

table(fdr[mods$mz[mods$type=='PTM'], ]$sqrtAge <=0.05)

a$ageFDR_kibbers <- fdr$sqrtAge

rownames(a) <- a$mz
 
table('all dogs'=a[ptm, ]$age_FDR<0.05, 'kibblers'=a[ptm, ]$ageFDR_kibbers<=0.05)
a$mz[a$age_FDR<0.05 & a$ageFDR_kibbers>0.05] # "N-Acetyl-Aspartate (NAA)" 'dropped out, however the rest reamin, including all of the 5 that are up in old dogs

