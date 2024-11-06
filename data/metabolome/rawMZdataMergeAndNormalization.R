library(tidyverse)
library(gtools)
library(impute)
library(dplyr)
library(outliers)
library(sva)
library(reshape)
library(ggbiplot)
library(gridExtra)
library(ggpubr)

rm(list=ls())
setwd('/Users/ben/Library/CloudStorage/GoogleDrive-brharrison2000@gmail.com/My Drive/Documents/DogAgingProject')


# load merged data on line: 289

########################################################################################
# load and merge raw LC-MS data:
########################################################################################
#Run Order with QC
batch1.run <- read.csv("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#1-3_Run_Order.csv")
batch2.run <- read.csv("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#4-7_Run_Order.csv")
batch3.run <- read.csv("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#8-13_Run_Order.csv")
batch4.run <- read.csv("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#14-21_Run_Order.csv")
batch5.run <- read.csv("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#22-35_Run_Order.csv")

#QC Data
batch1.QC <- read.csv("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#1-3_QC.csv", check.names = F)
batch2.QC <- read.csv("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#4-7_QC.csv", check.names = F)
batch3.QC <- read.csv("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#8-13_QC.csv", check.names = F)
batch4.QC <- read.csv("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#14-21_QC.csv", check.names = F)
batch5.QC <- read.csv("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#22-35_QC.csv", check.names = F)

#Loading metabolite files
batch1.data <- read.csv("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#1-3_Data.csv")
batch2.data <- read.csv("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#4-7_Data.csv")
batch3.data <- read.csv("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#8-13_Data.csv")
batch4.data <- read.csv("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#14-21_Data.csv")
batch5.data <- read.csv("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#22-35_Data.csv")

#Sample Info file for renaming
sampleinfo <- read.csv("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/sampleinfo.csv")

#Subset sampleinfo into three batches run
sample1 <- sampleinfo[sampleinfo$run == 1, ]
sample2 <- sampleinfo[sampleinfo$run == 2, ]
sample3 <- sampleinfo[sampleinfo$run == 3, ]
sample4 <- sampleinfo[sampleinfo$run == 4, ]
sample5 <- sampleinfo[sampleinfo$run == 5, ]

#Renumber Run Order
batch1.run$run <- row(batch1.run)
batch2.run$run <- row(batch2.run)
batch3.run$run <- row(batch3.run)
batch4.run$run <- row(batch4.run)
batch5.run$run <- row(batch5.run)

#Rename colnames
colnames(batch1.run) <- c('run_order', 'new_run_order')
colnames(batch2.run) <- c('run_order', 'new_run_order')
colnames(batch3.run) <- c('run_order', 'new_run_order')
colnames(batch4.run) <- c('run_order', 'new_run_order')
colnames(batch5.run) <- c('run_order', 'new_run_order')

#Renumbering Run Order
batch2.run$new_run_order <- batch2.run$new_run_order + 95
batch3.run$new_run_order <- batch3.run$new_run_order + 292
batch4.run$new_run_order <- batch4.run$new_run_order + 592
batch5.run$new_run_order <- batch5.run$new_run_order + 992

#Reorder Run
sample1 <- merge(sample1, batch1.run, by = "run_order", all = T)
sample1$new_SampleID <- paste("X", sample1$new_run_order, sep = "")
sample1 <- sample1[mixedorder(sample1[, 33]), ]
sample1 <- sample1[-62, ]

sample2 <- merge(sample2, batch2.run, by = "run_order", all = T)
sample2$new_SampleID <- paste("X", sample2$new_run_order, sep = "")
sample2 <- sample2[mixedorder(sample2[, 33]), ]

sample3 <- merge(sample3, batch3.run, by = "run_order", all = T)
sample3$new_SampleID <- paste("X", sample3$new_run_order, sep = "")
sample3 <- sample3[mixedorder(sample3[, 33]), ]
sample3 <- sample3[-c(80, 235), ]

sample4 <- merge(sample4, batch4.run, by = "run_order", all = T)
sample4$new_SampleID <- paste("X", sample4$new_run_order, sep = "")
sample4 <- sample4[mixedorder(sample4[, 33]), ]
sample4 <- sample4[-c(366,377,391,398), ]

sample5 <- merge(sample5, batch5.run, by = "run_order", all = T)
sample5$new_SampleID <- paste("X", sample5$new_run_order, sep = "")
sample5 <- sample5[mixedorder(sample5[, 33]), ]
sample5 <- sample5[-344, ]

batch1.data <- batch1.data[, -c(53)]
batch1.colnames <- sample1$new_SampleID[match(names(batch1.data), sample1$SampleID)]
names(batch1.data)[!is.na(batch1.colnames)] <- batch1.colnames[!is.na(batch1.colnames)]
batch1.data <- batch1.data

QC1.colnames <- sample1$new_SampleID[match(names(batch1.QC), sample1$run_order)]
names(batch1.QC)[!is.na(QC1.colnames)] <- QC1.colnames[!is.na(QC1.colnames)]
batch1.QC <- batch1.QC

#Merging two data sets
batch1.data <- merge(batch1.data[, -c(2,3)], batch1.QC, by.x = "COMPOUND", sort = F)

batch1.data[batch1.data == "N/A"] <- NA
batch1.standards <- batch1.data[-c(1:361), ]
batch1.data <- batch1.data[-c(362:393), ]

#Matching sample ID with QC included Run order
batch2.colnames <- sample2$new_SampleID[match(names(batch2.data), sample2$SampleID)]
names(batch2.data)[!is.na(batch2.colnames)] <- batch2.colnames[!is.na(batch2.colnames)]
batch2.data <- batch2.data

QC2.colnames <- sample2$new_SampleID[match(names(batch2.QC), sample2$run_order)]
names(batch2.QC)[!is.na(QC2.colnames)] <- QC2.colnames[!is.na(QC2.colnames)]
batch2.QC <- batch2.QC

#Merging two data sets
batch2.data <- merge(batch2.data[, -c(2,3)], batch2.QC, by.x = "COMPOUND", sort = F)

batch2.data[batch2.data == "N/A"] <- NA
batch2.standards <- batch2.data[-c(1:361), ]
batch2.data <- batch2.data[-c(362:393), ]

#Matching sample ID with QC included Run order
batch3.data <- batch3.data[, -c(67, 192)]
batch3.colnames <- sample3$new_SampleID[match(names(batch3.data), sample3$SampleID)]
names(batch3.data)[!is.na(batch3.colnames)] <- batch3.colnames[!is.na(batch3.colnames)]
batch3.data <- batch3.data

QC3.colnames <- sample3$new_SampleID[match(names(batch3.QC), sample3$run_order)]
names(batch3.QC)[!is.na(QC3.colnames)] <- QC3.colnames[!is.na(QC3.colnames)]
batch3.QC <- batch3.QC

#Merging two data sets
batch3.data <- merge(batch3.data[, -c(2,3)], batch3.QC, by.x = "COMPOUND", sort = F)

batch3.data[batch3.data == "N/A"] <- NA
batch3.standards <- batch3.data[-c(1:361), ]
batch3.data <- batch3.data[-c(362:393), ]

#Matching sample ID with QC included Run order
batch4.data <- batch4.data[, -c(295,304,316,323)]
batch4.colnames <- sample4$new_SampleID[match(names(batch4.data), sample4$SampleID)]
names(batch4.data)[!is.na(batch4.colnames)] <- batch4.colnames[!is.na(batch4.colnames)]
batch4.data <- batch4.data

QC4.colnames <- sample4$new_SampleID[match(names(batch4.QC), sample4$run_order)]
names(batch4.QC)[!is.na(QC4.colnames)] <- QC4.colnames[!is.na(QC4.colnames)]
batch4.QC <- batch4.QC

#Merging two data sets
batch4.data <- merge(batch4.data[, -c(2,3)], batch4.QC, by.x = "COMPOUND", sort = F)

batch4.data[batch4.data == "N/A"] <- NA
batch4.standards <- batch4.data[-c(1:361), ]
batch4.data <- batch4.data[-c(362:393), ]

##Matching sample ID with QC included Run order
batch5.data <- batch5.data[, -279]
batch5.colnames <- sample5$new_SampleID[match(names(batch5.data), sample5$SampleID)]
names(batch5.data)[!is.na(batch5.colnames)] <- batch5.colnames[!is.na(batch5.colnames)]
batch5.data <- batch5.data

QC5.colnames <- sample5$new_SampleID[match(names(batch5.QC), sample5$run_order)]
names(batch5.QC)[!is.na(QC5.colnames)] <- QC5.colnames[!is.na(QC5.colnames)]
batch5.QC <- batch5.QC

##Merging two data sets
batch5.data <- merge(batch5.data[, -c(2,3)], batch5.QC, by.x = "COMPOUND", sort = F)
batch5.data[batch5.data == "N/A"] <- NA
batch5.standards <- batch5.data[-c(1:361), ]
batch5.data <- batch5.data[-c(362:393), ]

#Reorder each batch
#Reorder column names based on sample ID
batch1.data <- batch1.data[ , c(1, 1 + mixedorder(names(batch1.data[ , -1])))]
batch2.data <- batch2.data[ , c(1, 1 + mixedorder(names(batch2.data[ , -1])))]
batch3.data <- batch3.data[ , c(1, 1 + mixedorder(names(batch3.data[ , -1])))]
batch4.data <- batch4.data[ , c(1, 1 + mixedorder(names(batch4.data[ , -1])))]
batch5.data <- batch5.data[ , c(1, 1 + mixedorder(names(batch5.data[ , -1])))]

#Reorder rows by alphabet
batch1.data <- batch1.data[mixedorder(batch1.data[ , 1]), ]
batch2.data <- batch2.data[mixedorder(batch2.data[ , 1]), ]
batch3.data <- batch3.data[mixedorder(batch3.data[ , 1]), ]
batch4.data <- batch4.data[mixedorder(batch4.data[ , 1]), ]
batch5.data <- batch5.data[mixedorder(batch5.data[ , 1]), ]

for (i in 1:length(sample1$run_order)){
  if (grepl("*S", sample1$run_order[i])){
    sample1$group[i] <- 'QC'
  }
  else if (grepl("*I", sample1$run_order[i])){
    sample1$group[i] <- 'Standard'
  }
  else{sample1$group[i] <- 'Sample'}
}

for (i in 1:length(sample2$run_order)){
  if (grepl("*S", sample2$run_order[i])){
    sample2$group[i] <- 'QC'
  }
  else if (grepl("*I", sample2$run_order[i])){
    sample2$group[i] <- 'Standard'
  }
  else{sample2$group[i] <- 'Sample'}
}

for (i in 1:length(sample3$run_order)){
  if (grepl("*S", sample3$run_order[i])){
    sample3$group[i] <- 'QC'
  }
  else if (grepl("*I", sample3$run_order[i])){
    sample3$group[i] <- 'Standard'
  }
  else{sample3$group[i] <- 'Sample'}
}

for (i in 1:length(sample4$run_order)){
  if (grepl("*S", sample4$run_order[i])){
    sample4$group[i] <- 'QC'
  }
  else if (grepl("*I", sample4$run_order[i])){
    sample4$group[i] <- 'Standard'
  }
  else{sample4$group[i] <- 'Sample'}
}

for (i in 1:length(sample5$run_order)){
  if (grepl("*S", sample5$run_order[i])){
    sample5$group[i] <- 'QC'
  }
  else if (grepl("*I", sample5$run_order[i])){
    sample5$group[i] <- 'Standard'
  }
  else{sample5$group[i] <- 'Sample'}
}

sample1$run <- 1
sample2$run <- 2
sample3$run <- 3
sample4$run <- 4
sample5$run <- 5


#Merging all batches
dat <- merge(batch1.data, batch2.data, by = "COMPOUND", sort = F)
dat <- merge(dat, batch3.data, by = "COMPOUND", sort = F)
dat <- merge(dat, batch4.data, by = "COMPOUND", sort = F)
dat <- merge(dat, batch5.data, by = "COMPOUND", sort = F)

#Reorder column names based on sample ID
dat <- dat[ , c(1, 1 + mixedorder(names(dat[ , -1])))]
#Reorder rows based off alphabet
dat <- dat[mixedorder(dat[ , 1]), ]

#Convert data frame to numeric matrix
#Separate column and row headers since those are not being treated as numbers
metabolite.rows <- dat[ , 1]
metabolite.cols <- names(dat[ , -1])

dat <- sapply(dat[ , -1], function(x) as.numeric(as.character(x)))

#Reattaching column and row headers
rownames(dat) <- metabolite.rows
colnames(dat) <- metabolite.cols

samps <- rbind(sample1, sample2, sample3, sample4, sample5)

panel <- rownames(dat)  # names of all metabolites on panel

########################################################################################
save(file='data/metabolome/merged.data', dat, samps, panel) # save mz data and sample info
########################################################################################



rm(list=ls())
################################################################################################
## quality control, normalization, imputation
################################################################################################
load('data/metabolome/merged.data')

par(mfrow=c(2,2))
boxplot(t(dat), pch=19, cex=0.1, xlab='metabolites', ylab='raw peak area') 
boxplot(dat[ ,1:100], pch=19, cex=0.1, xlab='samples', ylab='raw peak area')

dat <- log(dat)
dat <- scale(dat, scale=F)
dat <- as.data.frame(t(dat))

boxplot(dat, pch=19, cex=0.1, xlab='metabolites', main='log mzs\nscaled by sample') # by mz
boxplot(t(dat[ ,1:100]), pch=19, cex=0.1, xlab='samples', main='log mzs\nscaled by sample') # by samples

dat <- data.frame(samps, dat)
colnames(dat)[!colnames(dat) %in% colnames(samps)] <- panel # trick to thwart the re-coding of colnames for the mzs
dat$group <- as.factor(dat$group)
dat$prep_batch <- as.factor(dat$prep_batch) 
dat$run_order <- as.numeric(dat$run_order)
table(is.na(dat$run_order), dat$group) # run orders recorded for the samples, the interspersed QCs are within the run order
dat$dog_id <- as.character(dat$dog_id)
  
table(dat$group)
table(dat$graded_hemolysis, dat$group) # hemolysis score for plasma samples
table(dat$prep)
table(dat$cohort)
table(dat$cohort, dat$graded_hemolysis) # hemolysis among the samples from each cohort

# Exclude Hemolyzed Samples (remove score =4)---------panel# Exclude Hemolyzed Samples (remove score =4)--------------------------------------------
dat <- dat[dat$group=='Sample', ]
dat <- dat[dat$grade != 4, ]

# Removing the most- missing metabolites ----------------------------------------------
mzdat <- dat[ ,panel]
table(colSums(is.na(mzdat))) # NA count per metabolite, 108 complete metabolites, 180 completely missing metabolites
mzdat
par(mfrow=c(1,1))
hist(colSums(is.na(mzdat)/nrow(mzdat)), border=0, col='grey40', xlab='missingness (0 to 1)', main='', 30)

K <-print(round(nrow(mzdat) * 0.1)) # 10% of 1274 samples = 127
mzdat <- mzdat[ ,colSums(is.na(mzdat)) <= K] # 10% missingness threshold, no more than 10% missing

dogmzs <- colnames(mzdat) # these are the dog metabolites that we will work with
dat <- dat[ ,!colnames(dat) %in% panel] # remove mzs from dat
dat <- data.frame(dat, mzdat)
colnames(dat)[!colnames(dat) %in% colnames(samps)] <- dogmzs # trick to thwart the re-coding of colnames


# plots to show batch (prep) and LCMS run order effects:
par(mfrow=c(5,4), mar=c(4,4,1,1)+0.5)

for (i in 1:length(dogmzs)) {
  if(!all(is.na(dat[,dogmzs[i]]))) {
    plot(dat[ ,dogmzs[i]], pch=20, cex=0.5,
         main = dogmzs[i],
         xlab = "", ylab = "Abundance", las = 2, col = dat$prep_batch) +
      abline(v = c(74,229,444,726)) }} # vertical lines delineate LCMS runs

head(dat)
p1 <- ggplot(dat, aes(y=`Trimethylamine-N-Oxide (TMAO)`, x=1:nrow(dat), color=prep_batch))+
  geom_point(size=0.5)+
  xlab('')+
  theme_bw(base_size = 12)+
  geom_smooth(method='lm', se=F)+
  theme(legend.position = 'none')

p2 <- ggplot(dat, aes(y=`Trimethylamine-N-Oxide (TMAO)`, x=1:nrow(dat), color=prep_batch))+
  geom_point(size=0.5)+
  xlab('')+
  theme_bw(base_size = 8)+
  geom_smooth(method='lm', se=F)+
  facet_wrap(~prep_batch, scales='free')+
  theme(legend.position = 'none')

ggarrange(p1, p2)


# correct batch and run order
# then impute missing

table(colSums(is.na(dat[ ,dogmzs]))) # missingness, up to 125 NAs per metabolite, which is below 10% of the samples

table(dat$prep_batch) # seems prep batches are large enough to estimate a linear effect of run order.  
## doing mz ~ run order x prep batch would also better accomdate non-linear run order effects (see valine ~ run order in the last LCMS run = non-linear) 
## this would also simultaneously adjust for batch and LCMS run effects:

batchRunCorrected <- dat[ ,dogmzs]
# this loop keeps the NAs in the residuals
for(i in 1:length(dogmzs)){
  batchRunCorrected[ ,dogmzs[i]] <- residuals(lm(dat[ ,dogmzs[i]] ~ run_order * prep_batch, dat, na.action=na.exclude)) }


# confirm that this saved the position of all NAs
table(colSums(is.na(dat[ ,dogmzs])))
table(colSums(is.na(batchRunCorrected[ ,dogmzs])))
dat[,dogmzs][is.na(batchRunCorrected)] # should all be NA

par(mfrow=c(2,2))
plot(batchRunCorrected$`1-Methylnicotinamide`, dat$`1-Methylnicotinamide`)
plot(batchRunCorrected$Xanthine, dat$Xanthine)

par(mfrow=c(5,4), mar=c(4,4,1,1)+0.5)

for (i in 1:length(dogmzs)) {
  if(!all(is.na(dat[,dogmzs[i]]))) {
    plot(batchRunCorrected[ ,dogmzs[i]], pch=20, cex=0.5,
         main = dogmzs[i],
         xlab = "", ylab = "Abundance", las = 2, col = dat$prep_batch) +
      abline(v = c(74,229,444,726)) }} # vertical lines delineate LCMS runs

tmp <- dat
tmp[ ,dogmzs] <- batchRunCorrected

p3 <- ggplot(tmp, aes(y=`Trimethylamine-N-Oxide (TMAO)`, x=1:nrow(dat), color=prep_batch))+
  geom_point(size=0.5)+
  xlab('')+
  theme_bw(base_size = 12)+
  geom_smooth(method='lm', se=F)+
  theme(legend.position = 'none')

ggarrange(ggarrange(p1, p3), p2, ncol=1, heights = c(0.65, 1))

dat[ ,dogmzs] <- batchRunCorrected

# it looks like there are differences in mz variance among the LCMS runs.
# example TMAO

# it looks like there are differences in mz variance among the prep batches too.
# example pyruvate

# if we scale (DO NOT ALSO CENTER, scale ONLY, see GLUCOSE) by prep batch, it would handle both LCMS run and prep batch
dev.off()

scat <- dat %>% group_by(prep_batch) %>%  mutate(across(all_of(dogmzs), function(x) scale(x, center=F)))

par(mfrow=c(5,4), mar=c(4,4,1,1)+0.5)

for (i in 1:length(dogmzs)) {
  if(!all(is.na(dat[,dogmzs[i]]))) {
    plot(scat[ ,dogmzs[i]], pch=20, cex=0.5,
         main = dogmzs[i],
         xlab = "", ylab = "Abundance", las = 2, col = dat$prep_batch) +
      abline(v = c(74,229,444,726)) }} # vertical lines delineate LCMS runs


p4 <- ggplot(scat, aes(y=`Trimethylamine-N-Oxide (TMAO)`, x=1:nrow(dat), color=prep_batch))+
  geom_point(size=0.5)+
  xlab('')+
  theme_bw(base_size = 16)+
  geom_smooth(method='lm', se=F)+
  theme(legend.position = 'none')

ggarrange(p1, p3, p4, nrow=1)


# no remaining secular trends (nearly an axiom at this point).
#

# Summarize missingness
miss <- colSums(is.na(scat[ ,dogmzs]))
max(miss)
mean(miss)
max(miss)/nrow(scat)

# Impute
##############################################################
imp <- impute.knn(t(scat[ ,dogmzs]))
imp <- t(imp$data)
imp[1:4,1:4]

for (i in 1:length(dogmzs)) {
  if(!all(is.na(dat[,dogmzs[i]]))) {
    plot(imp[ ,dogmzs[i]], pch=20, cex=0.5,
         main = dogmzs[i],
         xlab = "", ylab = "Abundance", las = 2, col = dat$prep_batch) +
      abline(v = c(74,229,444,726)) }} # vertical lines delineate LCMS runs

dat[ ,dogmzs] <- imp


save(dat, dogmzs, file="data/metabolome/ProcessedData/normalizedData")
###############################
rm(list=ls())


load('data/metabolome/ProcessedData/normalizedData')


