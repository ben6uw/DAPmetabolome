library(ggplot2)
library(RColorBrewer)
library(plyr)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(usmap)
library(ggridges)

rm(list=ls())

setwd('/Users/ben/Library/CloudStorage/GoogleDrive-brharrison2000@gmail.com/My Drive/Documents/DogAgingProject')

#################################################################################################
## DAP, Precision and other Cohorts:
meta <- read.csv('data/2023_release/DAP_2023_SamplesResults_Metadata_v1.0.csv')
table(meta$Sample_Year, meta$Sample_Type)
m <- meta[meta$Sample_Type=='Metabolome' & meta$Sample_Year=='precision_1', ]
table(table(m$dog_id)) #1034 unique dogs

dogs <- read.csv('data/2023_release/SampledDogCohortKey_20240719.csv')
head(dogs)
table(dogs$Cohort_Label)

p <- dogs$dog_id[dogs$Cohort_Label=="Precision"]
bp <- dogs$dog_id[dogs$Cohort_Label=="BHS Precision"]
table(bp %in% p)
table(bp %in% m$dog_id)
table(p %in% m$dog_id)
#################################################################################################
rm(list=ls())


#################################################################################################
# signalment in the Precision Cohort baseline metabolome

rm(list=ls())
#################################################################################################
load("data/metabolome/ProcessedData/technicalCovarsRemovedData") # metabolome adjusted for technical covariates (shipping and LC-MS)

head(dat)
table(dat$cohort) # remove P2, and P3, merge EDTA with P!

dat$cohort[dat$cohort=='Precision 1 EDTA'] <- 'Precision 1'
dat <- dat[dat$cohort %in% 'Precision 1', ] #865 P1 dogs

save(dat, dogmzs, file='data/metabolome/ProcessedData/P1.technicalCovarsRemovedData')

## add CBC data
load('dog metabolome paper/data/blood.covariates.RData') # loads 'b', normalized CBC data for the P1 dogs
d <- merge(b, dat) # trims to 785 dogs
   
cellVars <- bloodvars[grepl('cbc', bloodvars)]
cellVars <- cellVars[!grepl('rel', cellVars)] # remove partly-redundant/colinear 'rel' relative counts, leave 'abs', absolute counts
cellVars

## add genetic breed data
load('dog metabolome paper/data/p1_update_March2024_with_Genetic_Breeds.RData') # loads new 'dat'
d <- merge(dat[ ,1:4], d) # trims to 784 dogs

d$sex <- as.factor(d$sex)
levels(d$sex) <- c('female', 'male')
d$sterilization_status <- as.factor(d$sterilization_status)
levels(d$sterilization_status) <- c('intact', 'sterilie')
table(d$sterilization_status, d$sex)


## define common breeds:
table(d$pct>0.6)
breedTable <- table(d$geneticTopBreed, '85%+ ancestry'=ifelse(d$pct>=0.85, "85%+", "<85%"))
write.csv(breedTable, quote=F, file='breedDistributionTable.csv') # distribution of breeds at >60% ancestry

ggplot(d, aes(x=100*pct))+
  geom_histogram(alpha=0.5, position = 'identity')+
  theme_classic(base_size = 16)+
  ylab('dogs')+
  xlab('top-breed ancestry (%)')

# representation of the top breeds:
head(breedTable)
topbreeds <- breedTable[ ,'85%+'][rev(order(breedTable[ ,'85%+']))][1:20] # the top 20 breeds
topbreeds <- topbreeds[topbreeds >=8] # at least 8 dogs
head(topbreeds)

cnt <- d[d$pct>=0.85, ] %>% group_by(geneticTopBreed, sex) %>%  dplyr::summarise(n=n())
head(cnt)

cnt <- cnt[cnt$geneticTopBreed %in% names(topbreeds), ]
cnt$geneticTopBreed <- factor(cnt$geneticTopBreed, levels=names(topbreeds))

ggplot(cnt, aes(geneticTopBreed, n, fill=sex)) + 
  geom_col() +
  theme_classic() +
  coord_flip()+
  scale_x_discrete(limits=rev)+
  ylab('dogs')+
  xlab('')


sum(cnt$n[cnt$geneticTopBreed != "Newfoundland"])
nrow(d)-sum(cnt$n[cnt$geneticTopBreed != "Newfoundland"]) # dogs not among the top breeds (excluding Newfoundlands)
cnt$geneticTopBreed

# exclude Newfoundland from breed-level analysis as it is confounded with sex
d$commonPurebred <- d$geneticTopBreed
d$commonPurebred[!d$geneticTopBreed %in% names(topbreeds)[names(topbreeds)!= "Newfoundland"] | d$pct <0.85] ='remaining dogs' # threshold for data among the common purebreds

table(d$commonPurebred)

d$commonPurebred <- as.factor(d$commonPurebred)
levels(d$commonPurebred)

levels(d$commonPurebred)[levels(d$commonPurebred)=="Cavalier King Charles Spaniel"] <- "Cav. King Ch. Spaniel"
levels(d$commonPurebred)[levels(d$commonPurebred)=="German Shepherd Dog"] <- "German Shepherd"
levels(d$commonPurebred)[levels(d$commonPurebred)=="American Pitbull Terrier"] <- "Am. Pitbull Terrier"

table(d$commonPurebred)
table(d$commonPurebred != 'remaining dogs')

cnt <- d %>% group_by(commonPurebred, sex) %>%  dplyr::summarise(n=n())
head(cnt)

x <- d %>% group_by(commonPurebred) %>% dplyr::summarise(n=n())
cnt$commonPurebred <- factor(cnt$commonPurebred, levels=rev(c(x$commonPurebred[order(x$n)])))

ggplot(subset(cnt, commonPurebred!= 'remaining dogs'), aes(commonPurebred, n, fill=sex)) + 
  geom_col() +
  theme_classic(base_size = 14) +
  coord_flip()+
  scale_x_discrete(limits=rev)+
  ylab('dogs')+
  xlab('')


#################################################################################################
#################################################################################################
d$lifestage_at_DOC <- factor(d$lifestage_at_DOC, levels=c('Puppy', 'Young', 'Mature', 'Senior'))

d$sex
ggplot(d, aes(x=age, fill=sex))+
  geom_histogram(alpha=0.5, position = 'identity')+
  theme_classic(base_size = 14)+
  ylab('dogs')+
  xlab('age (years)')

summary(d$age)
range(d$age)
table(d$sex)/nrow(d)
table(d$sterilization_status)/nrow(d)
table(d$sterilization_status, d$sex)
median(d$age[d$sex=='female']) #sexmedian(d$age[d$sex=='female']) # median female age
median(d$age[d$sex=='male']) # median male age


## weight distribution
breedmeans <- d %>% group_by(commonPurebred) %>% summarise_at(vars(weight_at_DOC), mean)
d$commonPurebred <- factor(d$commonPurebred, levels=c(breedmeans$commonPurebred[order(breedmeans$weight_at_DOC)]))

levels(d$commonPurebred)[levels(d$commonPurebred)=="Cavalier King Charles Spaniel"] <- "Cav. King Ch. Spaniel"
levels(d$commonPurebred)[levels(d$commonPurebred)=="German Shepherd Dog"] <- "German Shepherd"
levels(d$commonPurebred)[levels(d$commonPurebred)=="American Pitbull Terrier"] <- "Am. Pitbull Terrier"
levels(d$commonPurebred)[levels(d$commonPurebred)=="mixed_breed"] <- "remaining dogs"

mutts <- subset(d, commonPurebred == 'remaining dogs')
mutts$commonPurebred <- droplevels(mutts$commonPurebred)
mutts$commonPurebred <- 'remaining dogs'
range(d$weight_at_DOC)

p1 <- ggplot(subset(d, commonPurebred != 'remaining dogs'), aes(x=weight_at_DOC, y=commonPurebred, fill=commonPurebred, color = commonPurebred)) +
  geom_density_ridges(alpha = 0.5, color=NA)+
  theme_bw(base_size = 16)+
  theme(legend.position = 'none')+
  scale_y_discrete(limits=rev)+
  xlab('')+
  ylab('')+
  xlim(0, 100)


p2 <- ggplot(mutts, aes(x=weight_at_DOC, y=commonPurebred, fill=commonPurebred, color = commonPurebred))+
  geom_density_ridges(alpha = 0.5, color=NA, scale=18)+
  theme_bw(base_size = 16)+
  scale_fill_manual(values=c('grey')) +
  theme(legend.position = 'none')+
  scale_y_discrete(expand = expansion(add = c(0.15, 0.6)))+
  xlab('weight (kg)')+
  ylab('')+
  xlim(0, 100)

dev.off()
ggarrange(p1, p2, ncol=1, align='v', heights=c(1, 0.5))


####################################################################
### US map states

table(d$state)
mean(table(d$state))
range(table(d$state))
median(table(d$state))

states <- as.data.frame(table(d$state))
states
colnames(states) <- c('state', 'dogs')

usmap::plot_usmap(data=states, values='dogs', regions = "states") +
  scale_fill_continuous(states, low = "darkgoldenrod1", high = "red", name = "dogs", label = scales::comma, na.value="white") +
  labs(title = "Precision Cohort") +
  theme(legend.position = 'right', text = element_text(size = 16))+
  theme(panel.background = element_rect(color = "blue"))
  
dev.off()
hist(states$dogs, border=0, col='grey40', 20, main='', xlab='dogs', ylab='n states')
####################################################################

save(d, dogmzs, bloodvars, cellVars, file='data/metabolome/ProcessedData/P1.metabolome.withBreedData')





#################################################################################################
# PCA
#################################################################################################
rm(list=ls())
##############################################################################
load('data/metabolome/ProcessedData/P1.metabolome.withBreedData')

pca <- prcomp(d[ ,dogmzs], scale=T)

eigs <- pca$sdev^2
o <- AssocTests::tw(eigs, eigenL=length(eigs), criticalpoint=0.9793) # alpha at 5%
sigEigs <- print(o$SigntEigenL)
PCs <- as.data.frame(pca$x[ ,1:sigEigs])
pcs <- colnames(PCs)
PCs$dog_id <- d$dog_id

propVar <- print(summary(pca)$importance[2, 1:sigEigs]) # prop. var exp
sum(summary(pca)$importance[2, ][1:sigEigs]) # cumulative proportion of variance explained

d$sqrtAge <- sqrt(d$age)
d$sqrtWT <- sqrt(d$weight_at_DOC)

pcs <- names(PCs)
d <- merge(d, PCs)
d$commonPurebred <- as.factor(d$commonPurebred)
d$lifestage_at_DOC <- factor(d$lifestage_at_DOC, levels=c('Puppy', 'Young', 'Mature', 'Senior'))

cat(paste(cellVars, "+")) # trick to get model terms

library(car) # enables type III anova/ancova

aList <- list()

for(k in 1:sigEigs) {
aList[[k]] <- Anova(aov(lm(d[ ,pcs[k]] ~ sqrtAge + sqrtWT + lifestage_at_DOC + sex + sterilization_status + hours_fasting + commonPurebred + krt_cbc_total_wbcs + krt_cbc_abs_bands + krt_cbc_abs_neutrophils + krt_cbc_abs_lymphocytes + krt_cbc_abs_monocytes + krt_cbc_abs_eosinophils + krt_cbc_abs_basophils + krt_cbc_rbc + krt_cbc_hgb + krt_cbc_hct + krt_cbc_mcv + krt_cbc_mch + krt_cbc_mchc + krt_cbc_rdw + krt_cbc_mpv + krt_cbc_pct + krt_cbc_retic_abs, d)), type="III")  }

terms <- rownames(aList[[1]])

p <- as.data.frame(t(sapply(aList, function(x) x$`Pr(>F)`))) # pvalues
rownames(p) <- pcs[1:sigEigs]
colnames(p) <- terms
table(p$sqrtAge <0.05)

ss <- as.data.frame(t(sapply(aList, function(x) x$`Sum Sq`))) # sum of squares
rownames(ss) <- pcs[1:sigEigs]
colnames(ss) <- terms
ss$PC <- rownames(ss)

## only look at SS among the CBC variables:
ss <- reshape2::melt(ss[ ,c('PC', cellVars)])
head(ss)
colnames(ss)[2:3] <- c('term', 'SS') 
ss$term <- gsub('krt_cbc_', '', ss$term)
ss$PC <- factor(ss$PC, levels=pcs[1:sigEigs])
levels(ss$PC) <- 1:sigEigs
ss$term <- factor(ss$term)
levels(ss$term)
levels(ss$term) <- c("Abs Bands", "Abs Basophils", "Abs Eosinophils", "Abs Lymphocytes", "Abs Monocytes", "Abs Neutrophils", "HCT", "HGB", "MCH", "MCHC", "MCV", "MPV", "PCT", "RBC", "RDW", "Retic Abs Count", "Total WBCs") # simpler names - match names on Table S2

colours=pals::glasbey (nlevels(ss$term))

cbcANCOVAplot <- ggplot(ss, aes(y=SS, x=PC, fill=term))+
  geom_bar(stat='identity')+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(size = 12))+
  scale_fill_manual(values=colours, name = "CBC term")

cbcANCOVAplot



# combine the var by CBC variables
ss <- as.data.frame(t(sapply(aList, function(x) x$`Sum Sq`))) # sum of squares
rownames(ss) <- pcs[1:sigEigs]
colnames(ss) <- terms

ss$CBC <- apply(ss[ ,cellVars], 1, sum) 
ss <- ss[ ,!colnames(ss) %in% cellVars]
ss$PC <- rownames(ss)

ss <- reshape2::melt(ss)
head(ss)
colnames(ss)[2:3] <- c('term', 'SS') 
ss <- ss[ss$term!='(Intercept)', ]
ss$PC <- factor(ss$PC, levels=pcs[1:sigEigs])
levels(ss$PC) <- 1:sigEigs
ss$term <- factor(ss$term)
levels(ss$term)
ss$term <- factor(ss$term, levels=c("sqrtAge", "sqrtWT", "lifestage_at_DOC", "sex", "sterilization_status", "commonPurebred", "hours_fasting", "CBC", "Residuals"))
levels(ss$term) <- c('age', 'weight', 'lifestage', 'sex', 'sterilization', 'breed', 'fasting', 'cbc', 'residuals') # simpler names

nlevels(ss$term)
colors <- c(RColorBrewer::brewer.pal(nlevels(ss$term)-1, 'Set1'), 'grey80')

ggplot(ss, aes(y=SS, x=PC, fill=term))+
  geom_bar(stat='identity')+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(size = 12))+
  scale_fill_manual(values=colors)

ggplot(subset(ss, term!= 'residuals'), aes(y=SS, x=PC, fill=term))+
  geom_bar(stat='identity')+
  theme_bw(base_size = 16)+
  scale_fill_manual(values=colors, '')+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())


residualANCOVAplot <-  ggplot(ss, aes(y=SS, x=PC, fill=term))+
  geom_bar(stat='identity')+
  theme_classic(base_size = 16)+
  theme(axis.text.x = element_text(size = 12))+
  scale_fill_manual(values=colors, name = "ANCOVA term")


ggarrange(residualANCOVAplot, cbcANCOVAplot)



# PCS ~ common breed to see if PCs show breed-specificity:
pcVars <- paste0('PC', 1:sigEigs)

ANOVAbreedP <- sapply(aList, function(x) x$`Pr(>F)`[7])
breedPCs <- pcVars[ANOVAbreedP<=0.05] # PCS for which breed is significant at P<=5%

breedPCs

l <- pivot_longer(d[ ,colnames(d) %in% c(terms, breedPCs)], cols=all_of(breedPCs), names_to='PC')
head(l)
l$PC <- factor(l$PC, levels=breedPCs)
l$commonPurebred <- factor(l$commonPurebred, levels=unique(l$commonPurebred))
l$commonPurebred

ggplot(l, aes(y=value, x=commonPurebred, color=commonPurebred, fill=commonPurebred))+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  geom_jitter(width=0.01, size=0.2)+
  theme_classic(base_size = 12)+
  xlab('breed')+
  theme(axis.text.x=element_blank())+
  facet_wrap(~PC, nr=4, scales='free')+
  scale_fill_discrete('breed')+
  scale_color_discrete('breed')
    
# order: 
# breeds by size

l$sqrtWT
meanWT <- aggregate(sqrtWT ~ commonPurebred, l, mean)
meanWT
brwtOrder <- as.character(meanWT$commonPurebred[order(meanWT$sqrtWT)])
l$commonPurebred <- factor(l$commonPurebred, levels=brwtOrder)

ggplot(l, aes(y=value, x=commonPurebred, color=commonPurebred, fill=commonPurebred))+
  geom_boxplot(alpha=0.5, outlier.shape = NA, size=0.2)+
  geom_jitter(width=0.01, size=0.2)+
  theme_classic(base_size = 12)+
  theme(axis.text.x=element_blank())+
  facet_wrap(~PC, nr=1, scales='free')+
  ylab('PC value')+
  xlab('breeds, ordered by mean cohort dog weight')+
  scale_fill_discrete(name = "")+
  scale_color_discrete(name = "")

meanAge <- aggregate(sqrtAge ~ commonPurebred, l, mean)
plot(meanAge$sqrtAge, meanWT$sqrtWT, pch=19)


ggplot(d, aes(y=PC4, x=sqrtAge))+
  geom_point()+
  theme_classic(base_size = 16)+
  geom_smooth(method='lm', se=T, size=0) +
  xlab(expression("age" ~ sqrt(years)))

## the effect of age on PC4 in a full model:
##NOTE: ran all PCs in model mixed model with GRM below, beyond the ANOVA

summary(lm(PC4 ~ sqrtAge + sqrtWT + lifestage_at_DOC + sex + sterilization_status + hours_fasting + commonPurebred + krt_cbc_total_wbcs + krt_cbc_abs_bands + krt_cbc_abs_neutrophils + krt_cbc_abs_lymphocytes + krt_cbc_abs_monocytes + krt_cbc_abs_eosinophils + krt_cbc_abs_basophils + krt_cbc_rbc + krt_cbc_hgb + krt_cbc_hct + krt_cbc_mcv + krt_cbc_mch + krt_cbc_mchc + krt_cbc_rdw + krt_cbc_mpv + krt_cbc_pct + krt_cbc_retic_abs, d))

p1 <- ggplot(subset(ss, term!= 'residuals'), aes(y=SS, x=PC, fill=term))+
  geom_bar(stat='identity')+
  ylab('ANOVA (SS)')+
  theme_classic(base_size = 16)+
  scale_fill_manual(values=colors)+
  theme(strip.text.x = element_blank(),
      strip.background = element_rect(colour="white", fill="white"),
      legend.position=c(0.8,.7))

p2 <- ggplot(d, aes(y=PC4, x=sqrtAge))+
  geom_point(size=0.8)+
  theme_classic(base_size = 16)+
  geom_smooth(method='lm', se=T, size=0) +
  xlab(expression("age" ~ sqrt(years)))

ggarrange(p1, p2, widths=c(1, 0.9))


save(ss, terms, aList, file='dog metabolome paper/PCANOVA_result')

load('dog metabolome paper/PCANOVA_result')
ss[ss$term=='breed', ]


### how much of the variance can be explaiend by the signalment variables (1-residual SS) per PCs:
ss <- as.data.frame(t(sapply(aList, function(x) x$`Sum Sq`)))
rownames(ss) <- pcs[1:sigEigs]
colnames(ss) <- terms
ss
varExp <-  t(ss %>% ungroup() %>%
  mutate(across()/rowSums(across())))
round(varExp, 5)
head(varExp)

max(varExp['lifestage_at_DOC', ])

ss$totalvariance <- apply(ss, 1, sum) 

propVar <- apply(ss, 2, function(x) x/ss$totalvariance)
propVar[1, ]

apply(propVar, 2, max) # the maximum proportion of variance (of PCs) explained by each model term/level 

table(d$commonPurebred)

ss$CBC <- apply(ss[ ,cellVars], 1, sum) 
ss <- ss[ ,!colnames(ss) %in% cellVars]
head(ss)
hist(ss$CBC/ss$totalvariance)
range(ss$CBC/ss$totalvariance)
mean(ss$CBC/ss$totalvariance)

ss$PC <- rownames(ss)

range(1-(ss$Residuals/ss$totalvariance)) # range of the variance explained by all covariates

ss$PC <- factor(ss$PC, levels=pcs[1:sigEigs])
levels(ss$PC) <- 1:sigEigs

ggplot(ss, aes(y=1-(Residuals/totalvariance), PC))+
  geom_bar(stat='identity')+
  theme_classic(base_size = 16)+
  ylab('proportion of variance')

cat(paste0("'",levels(d$commonPurebred), "',"))

d$commonPurebred <- factor(d$commonPurebred, levels=c('Cav. King Ch. Spaniel', 'French Bulldog', 'Border Collie', 'Poodle', 'Golden Retriever', 'Labrador Retriever', 'German Shepherd', 'Great Dane', 'remaining dogs'))

breedcolors <- c(RColorBrewer::brewer.pal(12, name='Paired'), 'grey70')
breedcolors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#f542d4", "#B15928", "grey70") 

ggplot(d, aes(y=PC2, x=sqrtWT, color=commonPurebred))+
  geom_point()+
  theme_classic(base_size = 16)+
  geom_smooth(method='lm', se=F) +
  scale_color_manual(values=breedcolors)+
  labs(color="breed")+
  xlab(expression("weight" ~ sqrt(kg)))


save(PCs, sigEigs, file='mzPCs')



## PC snp heritability: variance in metabolome PCs explained by the relatedness in the grm:
rm(list=ls())
###################################

load('mzPCs')
load('dog metabolome paper/data/scaled.data_for_CBC_mixedModel') # as of Aug 27th 2024, the mz data have been replaced with the newly-normalized data

PCs <- merge(scDat, PCs)
table(PCs$dog_id %in% rownames(grm))
PCs <- PCs[PCs$dog_id %in% rownames(grm), ] # g is now composed of variables that have been scaled
grm <- grm[PCs$dog_id, PCs$dog_id] # subset and order the GRM by the dogs in the mz data

###################################
# fit a mixed model, with fixed effects of age, sex, weight, etc, in the context of random effects of relatedness (as represented in the GRM)

Zmat <- diag(nrow(PCs)) # an empty design matrix for the random effects

# build a design matrix for the fixed effects; specify interaction terms here
vars[!vars%in%dogmzs]

tmp <- ' ~ sqrtAge + sqrtWT + sex * sterilization_status + hours_fasting'
cat(c(tmp, paste0('+ ', cellVars))) # manually copy and paste this into the line below.

X <- model.matrix( ~ sqrtAge + sqrtWT + sex * sterilization_status + hours_fasting + krt_cbc_total_wbcs + krt_cbc_abs_bands + krt_cbc_abs_neutrophils + krt_cbc_abs_lymphocytes + krt_cbc_abs_monocytes + krt_cbc_abs_eosinophils + krt_cbc_abs_basophils + krt_cbc_rbc + krt_cbc_hgb + krt_cbc_hct + krt_cbc_mcv + krt_cbc_mch + krt_cbc_mchc + krt_cbc_rdw + krt_cbc_mpv + krt_cbc_pct + krt_cbc_retic_abs, data=PCs) 

pcVars <- paste0('PC', 1:sigEigs)
pcDat <- PCs[ ,pcVars]

###################################
# the next steps can be slow, use parallel processing
library(parallel)
library(doParallel)
library(EMMREML)

n.cores <- detectCores(all.tests = F, logical = T) # detects the number of available cores
clus <- makeCluster(n.cores) 
registerDoParallel(cores=n.cores)  

emma <- emmreml(y=pcDat[ ,1], X=X, Z=Zmat, K=grm, varbetahat = T, varuhat=T, PEVuhat=T, test=T) # the model fitting step, run an example to see what it does and to record the output rownames etc. (the code below will call the output)
emma$Vu/(emma$Vu+emma$Ve) # variance in Y explained by relatedness
 

###################################
# this code fits the same model twice, the first time it extracts the fixed effects, and the second time, it takes the random effects (they are called BLUPs):

# get fixed effects
clus <- makeCluster(n.cores) 
registerDoParallel(cores=n.cores)  
clusterExport(clus, varlist=c("pcDat", 'pcVars', 'X', 'Zmat', 'grm'), envir = environment()) 

mList <- list()

for(i in 1:length(pcVars)) {
  clusterExport(clus, varlist=c("pcDat", 'pcVars', 'X', 'Zmat', 'grm'), envir = environment()) 
  Y <- pcDat[ ,i]
  mList[[i]] <- emmreml(y=Y, X=X, Z=Zmat, K=grm, varbetahat = T, varuhat=T, PEVuhat=T, test=T) }


save(mList, PCs, pcVars, X, file='dog metabolome paper/PCMixedModel_results')




rm(list=ls())
##########################################################
# heritability of PCs
##########################################################
load('dog metabolome paper/PCMixedModel_results')

h <- sapply(mList, function(x) x$Vu/sum(x$Vu+x$Ve))
PCherit <- data.frame(PC=pcVars, heritability=h)
head(PCherit)

PCherit$PC <- factor(PCherit$PC, levels=pcVars)
levels(PCherit$PC) <- 1:sigEigs

ggplot(PCherit, aes(x = PC, y=heritability)) +
  geom_bar(stat="identity") +
  theme_bw(base_size = 14) +
  labs(y= expression(H["SNP"]))+
  xlab('Principal Component')
  

load('dog metabolome paper/PCANOVA_result')
head(ss)

nlevels(ss$term)
colors <- c(RColorBrewer::brewer.pal(nlevels(ss$term)-1, 'Set1'), 'grey80')

anovaPlot  <- ggplot(subset(ss, term!= 'residuals'), aes(y=SS, x=PC, fill=term))+
  geom_bar(stat='identity')+
  ylab('ANOVA (SS)')+
  xlab('Principal Component')+
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )+
  scale_fill_manual(values=colors, name = "")+
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour="black", fill="white"),
        legend.position=c(0.8,.63))

hsnpPlot <- ggplot(PCherit, aes(x = PC, y=heritability)) +
  geom_bar(stat="identity") +
  theme_bw(base_size = 14) +
 theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )+
  labs(y= expression(H["SNP"]))+
  xlab('Principal Component')


ggarrange(anovaPlot, NULL, hsnpPlot,   ncol = 1, heights=c(1, 0.05, 0.7))
  
ggarrange(ggarrange(anovaPlot, NULL, hsnpPlot,   ncol = 1, heights=c(1, 0.05, 0.7)), pc4dot)

save(file='Fig2_panels', anovaPlot, hsnpPlot, pc4dot)


mean(PCherit$heritability)
range(PCherit$heritability)

head(ss)
par(mfrow=c(1,1))
plot(PCherit$heritability ~ ss$SS[ss$term=='breed'], pch=19) # breed effect is not a good predictor of heritability
summary(lm(PCherit$heritability ~ ss$SS[ss$term=='breed']))

p <- numeric()
for(i in 1:nlevels(ss$term)) {
s <- summary(lm(PCherit$heritability ~ ss$SS[ss$term==levels(ss$term)[i]]))
p[i] <- s$coefficients[2,4] }
names(p) <- levels(ss$term)
p # fasting and sterilization show some evidence of assoc with heritability of the PC

plot(PCherit$heritability ~ ss$SS[ss$term=='fasting'], pch=19) # Ah! its negative, ie, the more a PC is affected by fasting, the less heritable it appears, that is sensible
plot(PCherit$heritability ~ ss$SS[ss$term=='sterilization'], pch=19) # Ah! same for sterilization










