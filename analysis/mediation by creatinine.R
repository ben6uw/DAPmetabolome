library(parallel)
library(doParallel)
library(mediation)

rm(list=ls())
###################################
load('data/scaled.data_for_CBC_mixedModel')


###################################
# can CREAT be a mediator? ...is it age associated?
###################################
plot(krt_cp_creatinine_value ~ sqrt(age), scDat)


Zmat <- diag(nrow(scDat))

X <- model.matrix( ~ sqrtAge + sqrtWT + sex * sterilization_status + hours_fasting + krt_cbc_total_wbcs + krt_cbc_abs_bands + krt_cbc_abs_neutrophils + krt_cbc_abs_lymphocytes + krt_cbc_abs_monocytes + krt_cbc_abs_eosinophils + krt_cbc_abs_basophils + krt_cbc_rbc + krt_cbc_hgb + krt_cbc_hct + krt_cbc_mcv + krt_cbc_mch + krt_cbc_mchc + krt_cbc_rdw + krt_cbc_mpv + krt_cbc_pct + krt_cbc_retic_abs, data=scDat) 

library(EMMREML)
emma <- emmreml(y=scDat$krt_cp_creatinine_value, X=X, Z=Zmat, K=grm, varbetahat = T, varuhat=T, PEVuhat=T, test=T)  
creatModel <- data.frame(P=emma$pvalbeta[ ,'none'], beta=emma$betahat)
head(creatModel) # creatinine ~ age P=0.84 
# creat does not fit the definition of a mediator as it is not age associated.
###################################
###################################



###################################
# can BUN be a mediator? ...is it age associated?
###################################
plot(krt_cp_bun_value ~ sqrt(age), scDat)

Zmat <- diag(nrow(scDat))

X <- model.matrix( ~ sqrtAge + sqrtWT + sex * sterilization_status + hours_fasting + krt_cbc_total_wbcs + krt_cbc_abs_bands + krt_cbc_abs_neutrophils + krt_cbc_abs_lymphocytes + krt_cbc_abs_monocytes + krt_cbc_abs_eosinophils + krt_cbc_abs_basophils + krt_cbc_rbc + krt_cbc_hgb + krt_cbc_hct + krt_cbc_mcv + krt_cbc_mch + krt_cbc_mchc + krt_cbc_rdw + krt_cbc_mpv + krt_cbc_pct + krt_cbc_retic_abs, data=scDat) 

library(EMMREML)
emma <- emmreml(y=scDat$krt_cp_bun_value, X=X, Z=Zmat, K=grm, varbetahat = T, varuhat=T, PEVuhat=T, test=T)  
bunModel <- data.frame(P=emma$pvalbeta[ ,'none'], beta=emma$betahat)
head(bunModel) # BUN ~ age P=0.86 
# BUN does not fit the definition of a mediator as it is not age associated.
###################################
###################################




###################################
# mediation analysis
###################################

sensitiveTest <- function(x) ifelse(x$ind.d0[x$rho==0]==0, 'passes', 'fails')
medOut<- function(x) data.frame('Total Effect'= x$tau.coef, 'ADE' = x$z0, 'ADElowerCI' = x$z0.ci[1], 'ADEupperCI'= x$z0.ci[2], 'ADE_P' = x$z0.p, 'ACME' = x$d0, 'ACMElowerCI' = x$d0.ci[1], 'ACMEupperCI'= x$d0.ci[2], 'ACME_P' = x$d0.p, 'ADEpropOfTotal Effect'=x$z0/x$tau.coef)

medModList <- list()
clusterExport(clus, varlist=c('scDat', 'dogmzs'), envir = environment()) 

for (i in 1:length(dogmzs)) {
  scDat$MZ <- scDat[ ,dogmzs[i]]
  med.fit <- lm(MZ ~ sqrtAge + sqrtWT + hours_fasting + sex * sterilization_status, data = scDat) 
  out.fit <- lm(MZ ~ krt_cp_creatinine_value + sqrtAge + sqrtWT + hours_fasting + sex * sterilization_status, data = scDat)
  medModList[[i]] <- mediate(med.fit, out.fit, treat = "sqrtAge", mediator = "krt_cp_creatinine_value",  sims=10000, boot=T) 
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
out$mediator <- 'krt_cp_creatinine_value'
save(out, medModList, sensList, file='mediationByCREAT')


load('mediationByCREAT')
out[out$ACME_FDR<=0.05, ]

plot(medModList[['N-Ac-Phenylalanine']])

# how often are the ptmAAs mediated by creat?
mods <- read.csv('dog metabolome paper/modAAs.csv')
table(mods$type)
out$PTM <- out$mz %in% mods$mz[mods$type=='PTM']

table('mediated'=out$ACME_FDR<=0.05, 'ptmAA'=out$PTM) # but, how many of these ptmAAs are age-associated????
fisher.test(table('mediated'=out$ACME_FDR<=0.05, 'ptmAA'=out$PTM))
a <- read.csv('ageMZsforDaniel.csv')

head(a)
head(out)
x <- merge(a, out)

table('mediated'=x$ACME_FDR<=0.05, 'ptmAA'=x$PTM, 'age.assoc'=x$age_FDR<=0.05) # the age effect of 6 of the 7 age-assoc ptmAAs are mediated by creatinine.  
x$mz[x$ACME_FDR<=0.05 & x$PTM & x$age_FDR<=0.05] # ptms mediated by creat
x$mz[x$ACME_FDR<=0.05 & x$PTM & x$age_FDR<=0.05 & x$B_age>0]
x$mz[x$ACME_FDR<=0.05 & x$PTM & x$age_FDR<=0.05 & x$B_age<0]

c("Dimethylarginine (A/SDMA)", "N-Ac-Alanine", "N-Ac-L-Glutamine", "N-Ac-Phenylalanine", "N-Ac-Tryptophan", "N-Acetyl-Aspartate (NAA)", 'Hydroxyproline')


barplot(out$ADEpropOfTotal.Effect[out$mz %in% x$mz[x$ACME_FDR<=0.05 & x$PTM & x$age_FDR<=0.05]])

hist(out$ADEpropOfTotal.Effect[out$ACME_FDR<=0.05])

s <- summary(medModList[['N-Ac-Phenylalanine']])
s$z0 / s$tau.coef

s

0.375/0.71 # ACME / Total Effect = prop mediated
plot(medModList[['N-Ac-Phenylalanine']])
s$n0
s


sumList <- lapply(medModList, function(x) summary(x))
p <- lapply(sumList, function(x) data.frame('propMed'=x$n0, 'LCI'=x$n.avg.ci[1], 'UCI'=x$n.avg.ci[2]))
p <- do.call(rbind, p)
p$mz <- rownames(p)

p <- p[p$mz %in% intersect(a$mz[a$age_FDR<=0.05], out$mz[out$ACME_FDR <=0.05]), ] # limit to only those age-assoc. mz, that are also mediated by creat

out$mz[out$ACME_FDR <=0.05]

ggplot(p, aes(y=propMed, x=mz))+
  geom_bar(stat='identity') + 
  geom_errorbar(aes(ymax = UCI, ymin=LCI, width=0.1))+
  theme_bw(base_size = 16)+
  theme(legend.position = 'none')+
  coord_flip(ylim = c(0, 1)) + 
  xlab('')+
  ylab('proportion mediated (+/- 95%CI)')+
  ggtitle('mediation by serum creatinine')


ggplot(p, aes(y=propMed, x=mz))+
  geom_bar(stat='identity') + 
  geom_errorbar(aes(ymax = UCI, ymin=LCI, width=0.1))+
  theme_bw(base_size = 16)+
  theme(legend.position = 'none')+
  coord_flip() + 
  xlab('')+
  ylab('proportion mediated (+/- 95%CI)')+
  ggtitle('mediation by serum creatinine')


p <- lapply(sumList, function(x) data.frame('propMed'=x$n0, 'LCI'=x$n.avg.ci[1], 'UCI'=x$n.avg.ci[2]))
p <- do.call(rbind, p)
p$mz <- rownames(p)

p <- p[c("Dimethylarginine (A/SDMA)", "N-Ac-Alanine", "N-Ac-L-Glutamine", "N-Ac-Phenylalanine", "N-Ac-Tryptophan", "N-Acetyl-Aspartate (NAA)", 'Hydroxyproline'), ]
p

p$mz <- factor(p$mz, levels=p$mz[order(p$propMed)])
levels(p$mz) <- c("hydroxyproline", "N-Ac-tryptophan", "N-Ac-glutamine", "N-Ac-phenylalanine", "N-Ac-aspartate", "N-Ac-alanine", "dimethylarginine")

ggplot(p, aes(y=propMed, x=mz, fill=as.factor(c(1,1,1,1,1,1,0))))+
  geom_bar(stat='identity', colour="black") + 
  geom_errorbar(aes(ymax = UCI, ymin=LCI, width=0.1))+
  theme_bw(base_size = 16)+
  theme(legend.position = 'none')+
  coord_flip(ylim = c(0, 1)) + 
  xlab('')+
  ylab('proportion mediated (+/- 95%CI)')+
  ggtitle('mediation by serum creatinine')+
  scale_fill_manual(values=c('white', 'grey40'))


range(p$propMed[1:6])
mean(p$propMed[1:6])
mean(p$propMed[c(2:4, 5)]) # up with age



# is creatinite (a proxy for GFR) associated with uSG??
plot(krt_cp_creatinine_value ~ krt_urine_sg, u)

load('dog metabolome paper/UrineMixedModelData')

tmp <- ' ~ sqrtAge + sqrtWT + sex * sterilization_status + hours_fasting'
tmp2 <- ', data=u)'
cat(c(tmp, paste0('+ ', cellCovars)), tmp2)

X <- model.matrix(   ~ krt_cp_creatinine_value + sqrtAge + sqrtWT + sex * sterilization_status + hours_fasting + krt_cbc_total_wbcs + krt_cbc_abs_bands + krt_cbc_abs_neutrophils + krt_cbc_abs_lymphocytes + krt_cbc_abs_monocytes + krt_cbc_abs_eosinophils + krt_cbc_abs_basophils + krt_cbc_rbc + krt_cbc_hgb + krt_cbc_hct + krt_cbc_mcv + krt_cbc_mch + krt_cbc_mchc + krt_cbc_rdw + krt_cbc_mpv + krt_cbc_pct + krt_cbc_retic_abs , data=u)

###################################
# the next steps can be slow, use parallel processing
library(parallel)
library(doParallel)
library(EMMREML)

n.cores <- detectCores(all.tests = F, logical = T) # detects the number of available cores
clus <- makeCluster(n.cores) 
registerDoParallel(cores=n.cores)  

emma <- emmreml(y=u$krt_urine_sg, X=X, Z=Zmat, K=grm, varbetahat = T, varuhat=T, PEVuhat=T, test=T) # the model fitting step, run an example to see what it does and to record the output rownames etc. (the code below will call the output)
head(emma)


data.frame('beta'=emma$betahat, 'P'=emma$pvalbeta[ ,'none'])

plot(krt_urine_sg ~ krt_cp_creatinine_value, u)






