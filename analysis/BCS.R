
rm(list=ls())

dir('data')
load('data/p1_update_March2024') # P1 metabolome data

v <- read.csv('data/PrecisionVEMRExtraction_20240326.csv')
colnames(v)
colnames(v) <- gsub('vemr_', '', colnames(v))

v$dog_id <- as.character(v$dog_id)
table(table(v$dog_id)) # 1002 dogs, each with one row
rownames(v) <- v$dog_id

table(dat$dog_id %in% v$dog_id) # 712 mz dogs in the VEMR data
rownames(dat) <- dat$dog_id
inter <- intersect(dat$dog_id , v$dog_id)
v <- v[inter, ]
dat <- dat[inter, ]

# [note]: these data contain a muscle condition score (mcs)
table(v$bcs)
table(v$mcs) # 99 = not found/unclear # only 26 informative muscle condition scores
table(v$mcs, v$bcs)


long <- pivot_longer(v[ ,c('dog_id', 'bcs', 'bcs_2', 'bcs_3')], -dog_id, values_to='bcs', names_to='record')

ggplot(long, aes(y=bcs, x=record, group=dog_id))+
  geom_jitter(width=0.01)+
  geom_line(alpha=0.1)+
  theme_classic(base_size = 18)+
  ggtitle('change in BCS\nnote "bcs" is the most contemporary to blood draw')

hist(v$bcs_3-v$bcs_2)
hist(v$bcs_2-v$bcs)
hist(c(v$bcs_3-v$bcs_2, v$bcs_2-v$bcs), border=0, col='grey30') # histogram of all delta-bcs

### will measure time btw BCS and blood draw after merging

## harmonize the 5-pt and 9-pt bcs: 
# [note]: as far as I can tell, each dog was only measured on one scale, even if they have >1 BCS measurement
v$bcs_choice
v$bcs_choice[v$bcs_choice==99] <- NA
v$bcs_choice
v$bcs_choice <- as.factor(v$bcs_choice)
levels(v$bcs_choice)
levels(v$bcs_choice) <- c('5-point', '9-point', 'unknown')

ggplot(subset(v, bcs_choice %in% c('5-point', '9-point')), aes(y=as.factor(bcs), fill=bcs_choice))+
  geom_bar(stat='count')+
  facet_wrap(~bcs_choice, ncol=1, scales="free_y")+
  theme_bw(base_size = 14)+
  coord_flip()+
  ylab('Body Condition Score')+
  xlab('P1 dogs')+
  ggtitle('BCS for 215 P1 dogs')

v$adjBCS <- v$bcs
v$adjBCS[v$bcs_choice == '5-point' & !is.na(v$bcs_choice)] <- 2*v$bcs[v$bcs_choice == '5-point' & !is.na(v$bcs_choice)]-1

p1 <- ggplot(subset(v, bcs_choice %in% c('5-point', '9-point')), aes(y=bcs, fill=bcs_choice))+
  geom_bar(stat='count')+
  facet_wrap(~bcs_choice, ncol=1, scales="free_y")+
  theme_bw(base_size = 14)+
  coord_flip()+
  ylab('Body Condition Score')+
  xlab('P1 dogs')+
  ggtitle('pre alignment')

p2 <- ggplot(subset(v, bcs_choice %in% c('5-point', '9-point')), aes(y=adjBCS, fill=bcs_choice))+
  geom_bar(stat='count')+
  facet_wrap(~bcs_choice, ncol=1, scales="free_y")+
  theme_bw(base_size = 14)+
  coord_flip()+
  ylab('Body Condition Score')+
  xlab('P1 dogs')+
  ggtitle('post alignment')

ggpubr::ggarrange(p1, p2, legend = F)

table(v$bcs, v$bcs_choice)
table(v$adjBCS, v$bcs_choice)
table(round(v$adjBCS), v$bcs_choice)

v$adjBCS_2 <- v$bcs_2
v$adjBCS_2[v$bcs_choice == '5-point' & !is.na(v$bcs_choice)] <- 2 * v$bcs_2[v$bcs_choice == '5-point' & !is.na(v$bcs_choice)]-1

v$adjBCS_3 <- v$bcs_3
v$adjBCS_3[v$bcs_choice == '5-point' & !is.na(v$bcs_choice)] <- 2 * v$bcs_3[v$bcs_choice == '5-point' & !is.na(v$bcs_choice)]-1

table(v$sex, dat$sex) # incompatible sex assignment in VEMR and mzdat, default to mz dat:

v <- merge(v[ ,!colnames(v)%in% 'sex'], dat)
v[1:10,1:4]
dat[1:10, 1:4]


## find the BCS closest to the DOC:
hist(as.numeric(difftime(as.Date(v$DOC), as.Date(v$bcs_date)))/365, border=0, col='grey30', 30, xlab='time between blood draw and (first) BCS record')

v$timelapse <- as.numeric(difftime(as.Date(v$DOC), as.Date(v$bcs_date)))/365
v$timelapse_2 <- as.numeric(difftime(as.Date(v$DOC), as.Date(v$bcs_date_2)))/365
v$timelapse_3 <- as.numeric(difftime(as.Date(v$DOC), as.Date(v$bcs_date_3)))/365

hist(v$timelapse, 30)
hist(c(v$timelapse, v$timelapse_2, v$timelapse_3), 30, border=0, col='grey30')

par(mfrow=c(1,2))
plot(v$timelapse, v$timelapse_2)
plot(v$timelapse, v$timelapse_3)
# the 'bcs' (not bcs_2 or bcs_3) are closest to the DOC

long <- pivot_longer(v[ c('dog_id', 'timelapse', 'timelapse_2', 'timelapse_3')], -dog_id)
head(long)

ggplot(long, aes(y=value, x=name, group=dog_id))+
  geom_jitter(width=0.01)+
  geom_line(alpha=0.1)+
  ylab('time until blood draw (years)')+
  theme_classic(base_size = 16)+
  ggtitle("time btw blood draw and bcs\nis shorest for 'bcs', not bcs_2, or bcs_3")


# the "SamplesResults_Metadata" file contains data taken on the same vet visit as the blood draw, when available, take their BCS:
d <- read.csv('data/SamplesResults_Metadata_20240112.csv')
d$dog_id <- as.character(d$dog_id)
colnames(d)
d[1:4,1:4]
table(d$Sample_Year)
d <- d[d$Sample_Year =='precision_1', ]
d <- d[d$Sample_Type =='Metabolome', ]
table(table(d$dog_id))
rownames(d) <- d$dog_id
d <- d[inter,]

table(d$Sample_Dog_Body_Condition_Score, d$Sample_Dog_Body_Condition_Score_Scale)
d$DOCadjBCS <- d$Sample_Dog_Body_Condition_Score
hist(d$DOCadjBCS)
table(is.na(d$Sample_Dog_Body_Condition_Score_Scale))

## convert 5pt to 9pt scale:
d$DOCadjBCS[d$Sample_Dog_Body_Condition_Score_Scale == 'Out of 5' & !is.na(d$Sample_Dog_Body_Condition_Score_Scale)] <- 2 * d$DOCadjBCS[d$Sample_Dog_Body_Condition_Score_Scale == 'Out of 5' & !is.na(d$Sample_Dog_Body_Condition_Score_Scale)]-1
hist(d$DOCadjBCS)

v <- merge(v, d[ ,c('dog_id', 'DOCadjBCS', 'Sample_Dog_Body_Condition_Score_Scale')], all.x=T)


table(v$adjBCS)
table(round(v$adjBCS)) # rounded with 0 digits (simplest)
v$rdBCS <- round(v$adjBCS)

par(mfrow=c(1,3))

plot(DOCadjBCS ~ rdBCS, v, pch=19, col=rgb(0,0,0,0.5), ylab='BCS at blood draw', xlab='BCS from VEMR', cex.lab=1.5)
abline(0,1, lty=2)

mean(v$DOCadjBCS-v$rdBCS, na.rm=T)
hist(v$DOCadjBCS-v$rdBCS, main='', xlab='difference in BCS (DOC-VEMR)', las=1, 8, cex.lab=1.5) # histogram of difference in BCS from DOC and from previous (but most recent) VEMR

plot(DOCadjBCS-rdBCS ~ timelapse, v[!is.na(v$DOCadjBCS) & !is.na(v$rdBCS), ], ylab='difference in BCS (DOC-VEMR)', xlab='time laspse (years)', pch=19, col=rgb(0,0,0,0.5), cex.lab=1.5)
abline(0, 0, lty=2)
abline(lm(DOCadjBCS-rdBCS ~ timelapse, v[!is.na(v$DOCadjBCS) & !is.na(v$rdBCS), ]))

table('BCS at blood draw'=v$DOCadjBCS[!is.na(v$DOCadjBCS) & !is.na(v$rdBCS)], 'VEMR BCS'=v$rdBCS[!is.na(v$DOCadjBCS) & !is.na(v$rdBCS)])

table(v$DOCadjBCS)

table(is.na(v$DOCadjBCS))
table(is.na(v$adjBCS))
table(is.na(v$DOCadjBCS), is.na(v$adjBCS))



## combine bcs, keep bcs at blood draw, and fill in most-resent past bcs of bcs missing from blood draw
v$agBCS <- v$DOCadjBCS
table(is.na(v$agBCS)) # 185 bcs records
v$agBCS[is.na(v$agBCS)] <- v$adjBCS[is.na(v$agBCS)]
table(is.na(v$agBCS)) # to 583 bcs records after filling in with VEMR records

v$rdBCS[v$rdBCS==10 & !is.na(v$rdBCS)] <- 9 # bump score of 10 back to the max of 9

###########################################
### proceed analyzing bcs
v$sqrtAge <- sqrt(v$age)

ggplot(v, aes(y=agBCS, x=sqrtAge))+
  geom_smooth(method='lm', size=0)+
  geom_point()+
  theme_bw(base_size = 14)

summary(lm(agBCS ~ sqrtAge, v))
summary(ordinal::clm(as.factor(agBCS) ~ sqrtAge, data=v)) 

ggplot(v, aes(y=rdBCS, x=sqrtAge))+
  geom_smooth(method='lm', size=0)+
  geom_point()+
  theme_bw(base_size = 14)

summary(lm(rdBCS ~ sqrtAge, v))
summary(ordinal::clm(as.factor(rdBCS) ~ sqrtAge, data=v))

ggplot(subset(v, !is.na(rdBCS)), aes(x=as.factor(rdBCS), y=sqrtAge))+
  geom_violin(trim=F, scale='count')+
  geom_jitter(alpha=0.5, width=0.05)+
  theme_bw(base_size = 14)+
  xlab('BCS')+
  ggtitle('n=583 P1 dogs\nordinal BCS ~ age P<0.00062')

v$sqrtWT <- sqrt(v$weight_at_DOC)

summary(ordinal::clm(as.factor(rdBCS) ~ sqrtAge, data=v))
summary(ordinal::clm(as.factor(rdBCS) ~ sqrtWT, data=v))
summary(ordinal::clm(as.factor(rdBCS) ~ sqrtAge + sqrtWT, data=v)) # WT is signficant on BCS. 
summary(ordinal::clm(as.factor(rdBCS) ~ sqrtAge * sqrtWT, data=v))  #  age x wt makes age affect go away.
# this implies that the effect of age on BCS is on weight ...could it also be that the effect of wt on BCS is through age!

psych::mediate(rdBCS ~ sqrtAge + (sqrtWT), data=v)
psych::mediate(rdBCS ~ sqrtWT + (sqrtAge), data=v)

w <- ggplot(subset(v, !is.na(rdBCS)), aes(x=as.factor(rdBCS), y=sqrtWT))+
  geom_violin(trim=F, scale='count')+
  geom_jitter(alpha=0.5, width=0.05)+
  theme_bw(base_size = 14)

a <- ggplot(subset(v, !is.na(rdBCS)), aes(x=as.factor(rdBCS), y=sqrtAge))+
  geom_violin(trim=F, scale='count')+
  geom_jitter(alpha=0.5, width=0.05)+
  theme_bw(base_size = 14)

ggarrange(a, w)

ggplot(subset(v, !is.na(rdBCS)), aes(y=sqrtWT, x=sqrtAge))+
  geom_smooth(method='lm', alpha=0.2)+
  geom_point(alpha=0.5)+
  theme_bw(base_size = 14)+
  xlab('age (sqrt Yrs)')+
  facet_wrap(~rdBCS, nrow=1)

plot(sqrtWT ~ as.factor(size_at_DOC), v)

v$size_at_DOC <- factor(v$size_at_DOC, c('Toy', 'Medium', 'Standard', 'Large', 'Giant'))
plot(sqrtWT ~ size_at_DOC, v)

ggplot(subset(v, !is.na(rdBCS)), aes(y=rdBCS, x=sqrtAge))+
  geom_smooth(method='lm', alpha=0.2)+
  geom_point(alpha=0.5)+
  theme_bw(base_size = 14)+
  xlab('age (sqrt Yrs)')+
  facet_wrap(~size_at_DOC, nrow=1)

summary(ordinal::clm(as.factor(rdBCS) ~ sqrtAge * size_at_DOC, data=v))
anova(ordinal::clm(as.factor(rdBCS) ~ sqrtAge * size_at_DOC, data=v))

v$clumpedSizes <- v$size_at_DOC 
v$clumpedSizes[v$clumpedSizes == 'Giant'] <- 'Large'
v$clumpedSizes <- droplevels(v$clumpedSizes)

summary(ordinal::clm(as.factor(rdBCS) ~ sqrtAge * clumpedSizes, data=v))
anova(ordinal::clm(as.factor(rdBCS) ~ sqrtAge * clumpedSizes, data=v))

ggplot(subset(v, !is.na(rdBCS)), aes(y=rdBCS, x=sqrtAge))+
  geom_smooth(method='lm', alpha=0.2, size=0)+
  geom_point(alpha=0.5)+
  theme_bw(base_size = 14)+
  xlab(expression('age ' * sqrt(years)))+
  facet_wrap(~clumpedSizes, nrow=1)

summary(lm(rdBCS ~ sqrtAge * clumpedSizes, data=v)) # linear-response model 
anova(lm(rdBCS ~ sqrtAge * clumpedSizes, data=v)) # linear-response model 

# other signalments: sex:
summary(ordinal::clm(as.factor(rdBCS) ~ sqrtAge * sex + sqrt(weight_at_DOC) + hours_fasting, data=v)) # sex, sex x age and fasting are not significant.
save(v, file='data/BCS_added')



rm(list=ls())
##################################################################

library(parallel)
library(doParallel)
library(mediation)

# merge cbc variables:
load('data/BCS_added')
load('data/data_for_CBC_mixedModel')

head(v)
v <- merge(b, v[ ,c('dog_id', 'rdBCS')]) 
v <- v[!is.na(v$rdBCS), ] # 512 dogs have both cbc and vemr records w/bcs

rm(grm)


### does wt (or size) mediate the BCS <- age relationship?
v$sqrtWT <- sqrt(v$weight_at_DOC)
v$sqrtAge <- sqrt(v$age)

v$size_at_DOC <- factor(v$size_at_DOC, c('Toy', 'Medium', 'Standard', 'Large', 'Giant'))
v$clumpedSizes <- v$size_at_DOC 
v$clumpedSizes[v$clumpedSizes == 'Giant'] <- 'Large'
v$clumpedSizes <- droplevels(v$clumpedSizes)

summary(ordinal::clm(as.factor(rdBCS) ~ sqrtAge + sex * sterilization_status + sqrt(weight_at_DOC) + hours_fasting, data=v)) 

ggplot(subset(v, !is.na(rdBCS)), aes(y=rdBCS, x=sqrtAge))+
  geom_smooth(method='lm', alpha=0.2, size=0)+
  geom_point(alpha=0.5)+
  theme_bw(base_size = 14)+
  xlab(expression('age ' * sqrt(years)))+
  facet_wrap(~clumpedSizes, nrow=1)

ggplot(subset(v, !is.na(rdBCS)), aes(y=rdBCS, x=sqrtAge))+
  geom_point(alpha=0.5)+
  theme_bw(base_size = 16)+
  ylab('Body Condition Score')+
  xlab(expression('age ' * sqrt(years)))

ggplot(subset(v, !is.na(rdBCS)), aes(x=rdBCS, y=sqrtAge))+
  geom_jitter(alpha=0.5, width=0.05)+
  theme_bw(base_size = 16)+
  coord_flip()+
  xlab('Body Condition Score')+
  ylab(expression('age ' * sqrt(years)))



## scale all variables
vars <- c('sqrtAge', 'sqrtWT', 'sex', 'sterilization_status', 'hours_fasting', dogmzs)
scDat <- v

boxplot(scDat[ ,vars]) # pre-scaling
scDat[ ,vars] <- scale(scDat[ ,vars])
boxplot(scDat[ ,vars]) # post-scaling

summary(ordinal::clm(as.factor(rdBCS) ~ sqrtAge * sex + sqrt(weight_at_DOC) + hours_fasting, data=scDat)) 

ggplot(subset(scDat, !is.na(rdBCS)), aes(y=rdBCS, x=sqrtAge))+
  geom_smooth(method='lm', alpha=0.2, size=0)+
  geom_point(alpha=0.5)+
  theme_bw(base_size = 14)+
  xlab(expression('scaled age ' * sqrt(years)))+
  facet_wrap(~clumpedSizes, nrow=1)


## does dog weight mediate the effect of age on BCS?
f1 <- lm(rdBCS ~ sqrtAge + hours_fasting + sex * sterilization_status, data = scDat)
f2 <- lm(rdBCS ~ sqrtAge + sqrtWT + hours_fasting + sex * sterilization_status, data = scDat) 
m <- mediate(f1, f2, treat = "sqrtAge", mediator = "sqrtWT",  sims=1000, boot=T) 
summary(m) # ambiguous, bc no effect from all terms

## does age meditate the effect of weight on BCS
f1 <- lm(rdBCS ~ sqrtWT+ hours_fasting + sex * sterilization_status, data = scDat)
f2 <- lm(rdBCS ~ sqrtWT + sqrtAge + hours_fasting + sex * sterilization_status, data = scDat) 
m <- mediate(f1, f2, treat = "sqrtWT", mediator = "sqrtAge",  sims=1000, boot=T) 
summary(m)

## does age meditate the effect of size on BCS
f1 <- lm(rdBCS ~ clumpedSizes + hours_fasting + sex * sterilization_status, data = scDat)
f2 <- lm(rdBCS ~ clumpedSizes + sqrtAge + hours_fasting + sex * sterilization_status, data = scDat) 
m <- mediate(f1, f2, treat = "clumpedSizes", mediator = "sqrtAge",  sims=1000, boot=T) 
summary(m) # ambiguous, bc no effect from all terms, model wants 2 factor levels for 'treat'

## does size meditate the effect of age on BCS
f1 <- lm(rdBCS ~ clumpedSizes + hours_fasting + sex * sterilization_status, data = scDat)
f2 <- lm(rdBCS ~ sqrtAge + clumpedSizes + hours_fasting + sex * sterilization_status, data = scDat) 
m <- mediate(f1, f2, treat = "sqrtAge", mediator = "clumpedSizes",  sims=1000, boot=T) 
summary(m) # does not handle clumped size as a mediator


## scale all variables
vars <- c('sqrtAge', 'sqrtWT', 'sex', 'sterilization_status', 'hours_fasting', dogmzs, bloodvars)
scDat <- v

table(scDat$rdBCS)
table(v$rdBCS)

boxplot(scDat[ ,vars]) # pre-scaling
scDat[ ,vars] <- scale(scDat[ ,vars])
boxplot(scDat[ ,vars]) # post-scaling

table(is.na(scDat$rdBCS))
scDat <- scDat[!is.na(scDat$rdBCS), ]


save(file='data.for.BCS.mediation', scDat, dogmzs, bloodvars, v)


rm(list=ls())
######################################################################################
load('data.for.BCS.mediation')

dogmzs
ggplot(scDat, aes(y=`4-Guanidinobutanoate`, x=sqrtAge))+
  geom_point()+
  theme_bw(base_size = 14)+
  geom_smooth(method='lm', size=0)+
  facet_wrap(~rdBCS)

scDat$clumpedBCS <- scDat$rdBCS
scDat$clumpedBCS[scDat$clumpedBCS %in% c(3, 4)] <- '3_4'
scDat$clumpedBCS[scDat$clumpedBCS %in% c(5:6)] <- '5_6'
scDat$clumpedBCS[scDat$clumpedBCS %in% c(7, 8, 9)] <- '7_9'
scDat$clumpedBCS <- as.factor(scDat$clumpedBCS)
table(scDat$clumpedBCS)


ggplot(scDat, aes(y=`4-Guanidinobutanoate`, x=sqrtAge))+
  geom_point()+
  theme_bw(base_size = 14)+
  geom_smooth(method='lm', size=0)+
  facet_wrap(~clumpedBCS)

anova(lm(`4-Guanidinobutanoate` ~ sqrtAge * clumpedBCS, scDat))

summary(lm(`4-Guanidinobutanoate` ~ sqrtAge, scDat))
summary(lm(`4-Guanidinobutanoate` ~ sqrtAge * rdBCS, scDat))


p1 <- ggplot(scDat, aes(y=`4-Guanidinobutanoate`, x=sqrtAge))+
  geom_point()+
  theme_bw(base_size = 14)+
  geom_smooth(method='lm', size=0)+
  facet_wrap(~clumpedBCS)

p2 <- ggplot(scDat, aes(y=Deoxycarnitine, x=sqrtAge))+
  geom_point()+
  theme_bw(base_size = 14)+
  geom_smooth(method='lm', size=0)+
  facet_wrap(~clumpedBCS)

p3 <- ggplot(scDat, aes(y=`Glutaric Acid`, x=sqrtAge))+
  geom_point()+
  theme_bw(base_size = 14)+
  geom_smooth(method='lm', size=0)+
  facet_wrap(~clumpedBCS)

ggarrange(p1, p2, p3, ncol=1)



######################################################################################
## mediation analysis
######################################################################################
load('data.for.BCS.mediation')


sensitiveTest <- function(x) ifelse(x$ind.d0[x$rho==0]==0, 'passes', 'fails')
medOut<- function(x) data.frame('ADE' = x$z0, 'ADElowerCI' = x$z0.ci[1], 'ADEupperCI'= x$z0.ci[2], 'ADE_P' = x$z0.p, 'ACME' = x$d0, 'ACMElowerCI' = x$d0.ci[1], 'ACMEupperCI'= x$d0.ci[2], 'ACME_P' = x$d0.p)

n.cores <- detectCores(all.tests = F, logical = T)
clus <- makeCluster(n.cores) 
registerDoParallel(cores=n.cores)  


# forward mediation path 
  medModList <- list()
  for (i in 1:length(dogmzs)) {
    clusterExport(clus, varlist=c('v', 'dogmzs'), envir = environment()) 
    scDat$MZ <- scDat[ ,dogmzs[i]]
    med.fit <- lm(rdBCS ~ sqrtAge + sqrtWT + hours_fasting + sex * sterilization_status, data = scDat) 
    out.fit <- lm(rdBCS ~ MZ + sqrtAge + sqrtWT + hours_fasting + sex * sterilization_status, data = scDat)
    medModList[[i]] <- mediate(med.fit, out.fit, treat = "sqrtAge", mediator = "MZ",  sims=1000, boot=T) 
  }
  
  names(medModList) <- dogmzs
  sensList <- lapply(medModList, function(x) medsens(x))
  sumList <- lapply(medModList, function(x) summary(x))
  out <- lapply(sumList, medOut)
  out <- do.call(rbind, out)
  out$mz <- dogmzs
  out$ADE_FDR <- p.adjust(out$ADE_P, 'fdr')
  out$ACME_FDR <- p.adjust(out$ACME_P, 'fdr')
  out$sensitivity <- sapply(sensList, sensitiveTest)
  
  save(out, medModList, sensList, file='BCS_mediation/mediation.rdBCS.by.MZs.forward')

  rm(sumList)
  rm(sensList)
  
 # reverse mediation path
  revModList <- list()
  for (i in 1:length(dogmzs)) {
    scDat$MZ <- scDat[ ,dogmzs[i]]
    med.fit <- lm(MZ ~ sqrtAge + sqrtWT + hours_fasting + sex * sterilization_status, data = scDat) 
    out.fit <- lm(MZ ~ rdBCS + sqrtAge + sqrtWT + hours_fasting + sex * sterilization_status, data = scDat)
    revModList[[i]] <- mediate(med.fit, out.fit, treat = "sqrtAge", mediator = "rdBCS", sims = 1000, boot=T) }
  
  names(revModList) <- dogmzs
  
  revsensList <- lapply(revModList, function(x) medsens(x))
  revsumList <- lapply(revModList, function(x) summary(x))
  revout <- lapply(revsumList, medOut)
  revout <- do.call(rbind, revout)
  revout$mz <- dogmzs
  revout$ADE_FDR <- p.adjust(revout$ADE_P, 'fdr')
  revout$ACME_FDR <- p.adjust(revout$ACME_P, 'fdr')
  revout$sensitivity <- sapply(revsensList, sensitiveTest)
  
  save(revout, revModList, revsensList, file='BCS_mediation/mediation.rdBCS.by.MZs.reverse')
  rm(revsumList)
  rm(revsensList)
  
  revout$direction <- 'MZ<-BCS<-age'
  out$direction <- 'BCS<-MZ<-age'
  
  o <- rbind(out,revout)
  o$direction
  o$direction <- as.factor(o$direction)
  
  ggplot(o, aes(y=ACME, x=direction, color=ACME_FDR<=0.05, group=mz))+
    geom_line(alpha=0.5, color='grey')+
    geom_jitter(width=0.01, size=1)+
    theme_classic(base_size=16)  

  
# the FDR and P for 4-Guanidinobutanoate, Deoxycarnitine, Glutaric Acid, are 0, so run more perms to resolve them:
load('data.for.BCS.mediation')
load('BCS_mediation/mediation.rdBCS.by.MZs.reverse')
revout$ACME_FDR 

# deeper look
digTheseOut <- c('4-Guanidinobutanoate', 'Deoxycarnitine', 'Glutaric Acid')
deeper <- list()
for (i in 1:length(digTheseOut)) {
  scDat$MZ <- scDat[ ,digTheseOut[i]]
  med.fit <- lm(MZ ~ sqrtAge + sqrtWT + hours_fasting + sex * sterilization_status, data = scDat) 
  out.fit <- lm(MZ ~ rdBCS + sqrtAge + sqrtWT + hours_fasting + sex * sterilization_status, data = scDat)
  deeper[[i]] <- mediate(med.fit, out.fit, treat = "sqrtAge", mediator = "rdBCS", sims = 100000, boot=T) }

names(deeper) <- digTheseOut

medOut<- function(x) data.frame('ADE' = x$z0, 'ADElowerCI' = x$z0.ci[1], 'ADEupperCI'= x$z0.ci[2], 'ADE_P' = x$z0.p, 'ACME' = x$d0, 'ACMElowerCI' = x$d0.ci[1], 'ACMEupperCI'= x$d0.ci[2], 'ACME_P' = x$d0.p)

sumList <- lapply(deeper, function(x) summary(x))
revout <- lapply(sumList, medOut)
revout <- do.call(rbind, revout)
revout$mz <- digTheseOut
revout$ADE_FDR <- p.adjust(revout$ADE_P, 'fdr', n=137)
revout$ACME_FDR <- p.adjust(revout$ACME_P, 'fdr', n=137)

revout

scDat$MZ <- scDat[ ,'Glutaric Acid']
med.fit <- lm(MZ ~ sqrtAge + sqrtWT + hours_fasting + sex * sterilization_status, data = scDat) 
out.fit <- lm(MZ ~ rdBCS + sqrtAge + sqrtWT + hours_fasting + sex * sterilization_status, data = scDat)
m <- mediate(med.fit, out.fit, treat = "sqrtAge", mediator = "rdBCS", sims = 1000000, boot=T) 

revout['Glutaric Acid', ] <- medOut(summary(m))
revout$ADE_FDR <- p.adjust(revout$ADE_P, 'fdr', n=137)
revout$ACME_FDR <- p.adjust(revout$ACME_P, 'fdr', n=137)
revout
save(revout, file= 'deeperPermResults')

revout$ACMEupperCI- revout$ACME
revout$ACMElowerCI- revout$ACME

rm(list=ls())
################################### ################################### ###################################
## July 9th 2024, conversaion with Roger Fielding at Tufts, he thought about asking if the effect of age on ptmAAs was affected by 'what a dog is working with' (ie, its BCS).  Specifically:

# Is the effect of age on ptmAA mediated by BCS?? 
# before mediation analysis, test if BCS matters for mz while controlling for all other covars in the mixed model:  
###################################
  load('data/scaled.data_for_CBC_mixedModel')
  load('data.for.BCS.mediation')
  
  scDat$rdBCS[ scDat$rdBCS==10] <- 9 # join the BCS 10 with the 9s
  
  
  table(scDat$dog_id %in% rownames(grm))
  grm <- grm[scDat$dog_id, scDat$dog_id] # subset and order the GRM by the dogs in the mz data
  
  ###################################
  # fit a mixed model, with fixed effects of age, sex, weight, etc, in the context of random effects of relatedness (as represented in the GRM)
  
  Zmat <- diag(nrow(scDat)) # an empty design matrix for the random effects
  
  # build a design matrix for the fixed effects; specify interaction terms here
  vars[!vars%in%dogmzs]
  
  tmp <- ' ~ rdBCS + sqrtAge + sqrtWT + sex * sterilization_status + hours_fasting'
  cat(c(tmp, paste('+', cellCovars)))
  
  X <- model.matrix(  ~ rdBCS + sqrtAge + sqrtWT + sex * sterilization_status + hours_fasting + krt_cbc_abs_bands + krt_cbc_abs_eosinophils + krt_cbc_abs_lymphocytes + krt_cbc_abs_monocytes + krt_cbc_abs_neutrophils + krt_cbc_hct + krt_cbc_hgb + krt_cbc_mch + krt_cbc_mcv + krt_cbc_mpv + krt_cbc_pct + krt_cbc_rbc + krt_cbc_rdw + krt_cbc_rel_bands + krt_cbc_rel_lymphocytes + krt_cbc_rel_monocytes + krt_cbc_rel_neutrophils + krt_cbc_retic_abs + krt_cbc_retic_per + krt_cbc_total_wbcs, data=scDat) 
  
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
colnames(betas) <- rownames(emma$betahat)  
Ps <- fixedEffects[ ,c((1+length(rownames(emma$betaha))):ncol(fixedEffects))]
colnames(Ps) <- rownames(emma$betahat)  

fdr <- as.data.frame(apply(Ps, 2, function(x) p.adjust(x, 'fdr')))
colSums(fdr<=0.05)

dogmzs[fdr$rdBCS<=0.05] # four mzs assoc with BCS
betas['Citrulline',]

## plot MZ ~ BCS using model residuals (BLUPs, BLUEs)

# get BLUEs (subtract these from y to get the residual Y after the fixed effects):
bhats <- as.matrix(t(betas))
subbhats <- bhats[rownames(bhats) != 'rdBCS', ] # remove variable of interest
subX <- X[ ,colnames(X)!= 'rdBCS']
BLUEs <- subX %*% subbhats

BLUEs

# get random effects (BLUPs)
randomEffects<- parApply(clus, mzdat, 2, function(Y) {
  library(EMMREML)
  emma <- emmreml(y=Y, X=X, Z=Zmat, K=grm, varbetahat = T, varuhat=T, PEVuhat=T, test=T) 
  blup <- emma$uhat
  varblup <- emma$varuhat
  return(list(blup, varblup)) } )

save(file='BCSmixed_model_results', fixedEffects, randomEffects, BLUEs, scDat, dogmzs, X, v)




rm(list=ls())
######################################################################
# summarize and plot results
load('BCSmixed_model_results')


betas <- fixedEffects[ ,1:(ncol(fixedEffects)/2)]
colnames(betas) <- colnames(X)  
Ps <- fixedEffects[ ,c((1+(ncol(fixedEffects)/2)):ncol(fixedEffects))]
colnames(Ps) <- colnames(X)  

fdr <- as.data.frame(apply(Ps, 2, function(x) p.adjust(x, 'fdr')))
colSums(fdr<=0.05)

dogmzs[fdr$rdBCS<=0.05] # four mzs assoc with BCS

# get residuals
mmResiduals <- matrix(nr=nrow(scDat), nc=length(dogmzs))

for(i in 1:length(dogmzs)) {
  BLUPS <- unlist(randomEffects[[dogmzs[i]]][1]) # the BLUPS for a given mz
  mmResiduals [ ,i] <- scDat[ ,dogmzs[i]] - BLUEs[ ,i] - BLUPS }
colnames(mmResiduals) <- dogmzs

tmp <- scDat
tmp[ ,dogmzs] <- mmResiduals
mmResiduals <- tmp
rm(tmp)

dogmzs[fdr$rdBCS<=0.05] # four mzs assoc with BCS

bcsAgePlot <- ggplot(subset(v, !is.na(rdBCS)), aes(x=rdBCS, y=sqrtAge))+
  geom_jitter(alpha=0.5, width=0.05)+
  theme_bw(base_size = 12)+
  coord_flip()+
  xlab('Body Condition Score')+
  ylab(expression('age ' * sqrt(years)))

p1 <- ggplot(scDat, aes(y=`4-Guanidinobutanoate`, x=rdBCS))+
  geom_point()+
  geom_smooth(size=0, alpha=0.5, method='lm')+
  xlab('Body Condition Score')+
  theme_bw(base_size = 12)

p2 <-ggplot(scDat, aes(y=Deoxycarnitine, x=rdBCS))+
  geom_point()+
  geom_smooth(size=0, alpha=0.5, method='lm')+
  xlab('Body Condition Score')+
  theme_bw(base_size = 12)

p3 <- ggplot(scDat, aes(y=`Glutaric Acid`, x=rdBCS))+
  geom_point()+
  geom_smooth(size=0, alpha=0.5, method='lm')+
  xlab('Body Condition Score')+
  theme_bw(base_size = 12)

p4 <- ggplot(scDat, aes(y=Citrulline, x=rdBCS))+
  geom_point()+
  geom_smooth(size=0, alpha=0.5, method='lm')+
  xlab('Body Condition Score')+
  theme_bw(base_size = 12)

ggarrange(bcsAgePlot, p1, p2, p3, p4, nrow=1)

  
rm(list=ls())
#######################################
load('data/scaled.data_for_CBC_mixedModel')
load('data.for.BCS.mediation')
load('BCS_mediation/mediation.rdBCS.by.MZs.reverse')

rm(revsumList)
rm(revsensList)

revout$direction <- 'MZ<-BCS<-age'

head(revout)
dogmzs[revout$ACME_FDR <=0.05]

scDat$MZ <- scDat$`4-Guanidinobutanoate`

med.fit <- lm(MZ ~ sqrtAge + sqrtWT + hours_fasting + sex * sterilization_status+ krt_cbc_abs_bands + krt_cbc_abs_eosinophils + krt_cbc_abs_lymphocytes + krt_cbc_abs_monocytes + krt_cbc_abs_neutrophils + krt_cbc_hct + krt_cbc_hgb + krt_cbc_mch + krt_cbc_mcv + krt_cbc_mpv + krt_cbc_pct + krt_cbc_rbc + krt_cbc_rdw + krt_cbc_rel_bands + krt_cbc_rel_lymphocytes + krt_cbc_rel_monocytes + krt_cbc_rel_neutrophils + krt_cbc_retic_abs + krt_cbc_retic_per + krt_cbc_total_wbcs, data = scDat) 

out.fit <- lm(MZ ~ rdBCS + sqrtAge + sqrtWT + hours_fasting + sex * sterilization_status+ krt_cbc_abs_bands + krt_cbc_abs_eosinophils + krt_cbc_abs_lymphocytes + krt_cbc_abs_monocytes + krt_cbc_abs_neutrophils + krt_cbc_hct + krt_cbc_hgb + krt_cbc_mch + krt_cbc_mcv + krt_cbc_mpv + krt_cbc_pct + krt_cbc_rbc + krt_cbc_rdw + krt_cbc_rel_bands + krt_cbc_rel_lymphocytes + krt_cbc_rel_monocytes + krt_cbc_rel_neutrophils + krt_cbc_retic_abs + krt_cbc_retic_per + krt_cbc_total_wbcs, data = scDat)

m1 <- mediate(med.fit, out.fit, treat = "sqrtAge", mediator = "rdBCS", sims = 1000, boot=T) # may have convergece problems with so many mediators
summary(m1)

plot(m1, main='4-Guanidinobutanoate ~ age\n mediation by BCS')


scDat$MZ <- scDat$Deoxycarnitine

med.fit <- lm(MZ ~ sqrtAge + sqrtWT + hours_fasting + sex * sterilization_status+ krt_cbc_abs_bands + krt_cbc_abs_eosinophils + krt_cbc_abs_lymphocytes + krt_cbc_abs_monocytes + krt_cbc_abs_neutrophils + krt_cbc_hct + krt_cbc_hgb + krt_cbc_mch + krt_cbc_mcv + krt_cbc_mpv + krt_cbc_pct + krt_cbc_rbc + krt_cbc_rdw + krt_cbc_rel_bands + krt_cbc_rel_lymphocytes + krt_cbc_rel_monocytes + krt_cbc_rel_neutrophils + krt_cbc_retic_abs + krt_cbc_retic_per + krt_cbc_total_wbcs, data = scDat) 

out.fit <- lm(MZ ~ rdBCS + sqrtAge + sqrtWT + hours_fasting + sex * sterilization_status+ krt_cbc_abs_bands + krt_cbc_abs_eosinophils + krt_cbc_abs_lymphocytes + krt_cbc_abs_monocytes + krt_cbc_abs_neutrophils + krt_cbc_hct + krt_cbc_hgb + krt_cbc_mch + krt_cbc_mcv + krt_cbc_mpv + krt_cbc_pct + krt_cbc_rbc + krt_cbc_rdw + krt_cbc_rel_bands + krt_cbc_rel_lymphocytes + krt_cbc_rel_monocytes + krt_cbc_rel_neutrophils + krt_cbc_retic_abs + krt_cbc_retic_per + krt_cbc_total_wbcs, data = scDat)

m2 <- mediate(med.fit, out.fit, treat = "sqrtAge", mediator = "rdBCS", sims = 1000, boot=T)
summary(m2)

plot(m2, main='Deoxycarnitine ~ age\n mediation by BCS')


betas["Deoxycarnitine",]
fdr["Deoxycarnitine",]
summary(m2)

scDat$MZ <- scDat$`Glutaric Acid`

med.fit <- lm(MZ ~ sqrtAge + sqrtWT + hours_fasting + sex * sterilization_status+ krt_cbc_abs_bands + krt_cbc_abs_eosinophils + krt_cbc_abs_lymphocytes + krt_cbc_abs_monocytes + krt_cbc_abs_neutrophils + krt_cbc_hct + krt_cbc_hgb + krt_cbc_mch + krt_cbc_mcv + krt_cbc_mpv + krt_cbc_pct + krt_cbc_rbc + krt_cbc_rdw + krt_cbc_rel_bands + krt_cbc_rel_lymphocytes + krt_cbc_rel_monocytes + krt_cbc_rel_neutrophils + krt_cbc_retic_abs + krt_cbc_retic_per + krt_cbc_total_wbcs, data = scDat) 

out.fit <- lm(MZ ~ rdBCS + sqrtAge + sqrtWT + hours_fasting + sex * sterilization_status+ krt_cbc_abs_bands + krt_cbc_abs_eosinophils + krt_cbc_abs_lymphocytes + krt_cbc_abs_monocytes + krt_cbc_abs_neutrophils + krt_cbc_hct + krt_cbc_hgb + krt_cbc_mch + krt_cbc_mcv + krt_cbc_mpv + krt_cbc_pct + krt_cbc_rbc + krt_cbc_rdw + krt_cbc_rel_bands + krt_cbc_rel_lymphocytes + krt_cbc_rel_monocytes + krt_cbc_rel_neutrophils + krt_cbc_retic_abs + krt_cbc_retic_per + krt_cbc_total_wbcs, data = scDat)

m3 <- mediate(med.fit, out.fit, treat = "sqrtAge", mediator = "rdBCS", sims = 1000, boot=T)
summary(m3)

plot(m3, main='Glutaric Acid ~ age\n mediation by BCS')


par(mfrow=c(3,1))
summary(m1)
summary(m2)
summary(m3)
plot(m1, main='4-Guanidinobutanoate ~ age\n mediation by BCS')
plot(m2, main='Deoxycarnitine ~ age\n mediation by BCS')
plot(m3, main='Glutaric Acid ~ age\n mediation by BCS')




