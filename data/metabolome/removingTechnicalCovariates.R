library(tidyverse)
library(gtools)
library(impute)
library(dplyr)
library(outliers)
library(sva)
library(reshape)
library(ggbiplot)
library(gridExtra)
library(here)

rm(list=ls())
setwd('/Users/ben/Library/CloudStorage/GoogleDrive-brharrison2000@gmail.com/My Drive/Documents/DogAgingProject')

load("data/metabolome/ProcessedData/normalizedData")

par(mfrow=c(2,2))
## adjust for technical covariates (pre-LCMS variables)
hist(dat$elapse_DOC_LOG) # travel time from collection location to Texas A&M
hist(dat$serum_temp) # sample temperature at Texas A&M
plot( serum_temp ~ elapse_DOC_LOG, dat, las=1)

## get an idea of the issues:
p <- list()
for(i in 1:length(dogmzs)){
s <- summary(lm(dat[ ,dogmzs[i]] ~ elapse_DOC_LOG + serum_temp, dat, na.action=na.exclude))
p[[i]] <- s$coefficients[2:3,4] }
p <- as.data.frame(do.call(rbind, p))
rownames(p) <- dogmzs

p[order(p$elapse_DOC_LOG), ][1:10, ]

par(mfrow=c(2,2))
plot(Choline ~ elapse_DOC_LOG, dat, pch=19, cex=0.5) # plot the offended
plot(Choline ~ serum_temp, dat, pch=19, cex=0.5) # plot an offender


tech <- dat[ ,dogmzs]

for(i in 1:length(dogmzs)){
  tech[ ,dogmzs[i]] <- residuals(lm(tech[ ,dogmzs[i]] ~ elapse_DOC_LOG + serum_temp, dat, na.action=na.exclude)) } # 6 samples with NA for serum temp

dat[ ,dogmzs] <- tech

plot(Choline ~ elapse_DOC_LOG, dat, pch=19, cex=0.5) # plot the offended
plot(Choline ~ serum_temp, dat, pch=19, cex=0.5) # plot an offender

save(dat, dogmzs, file="data/metabolome/ProcessedData/technicalCovarsRemovedData")

rm(list=ls())




