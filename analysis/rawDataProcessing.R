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

here::i_am("data/metabolome/DataCleanUp20231120.Rmd")

#Run Order with QC
batch1.run <- read.csv(here("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#1-3_Run_Order.csv"))
batch2.run <- read.csv(here("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#4-7_Run_Order.csv"))
batch3.run <- read.csv(here("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#8-13_Run_Order.csv"))
batch4.run <- read.csv(here("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#14-21_Run_Order.csv"))
batch5.run <- read.csv(here("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#22-35_Run_Order.csv"))

#QC Data
batch1.QC <- read.csv(here("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#1-3_QC.csv"), check.names = F)
batch2.QC <- read.csv(here("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#4-7_QC.csv"), check.names = F)
batch3.QC <- read.csv(here("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#8-13_QC.csv"), check.names = F)
batch4.QC <- read.csv(here("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#14-21_QC.csv"), check.names = F)
batch5.QC <- read.csv(here("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#22-35_QC.csv"), check.names = F)

#Loading metabolite files
batch1.data <- read.csv(here("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#1-3_Data.csv"))
batch2.data <- read.csv(here("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#4-7_Data.csv"))
batch3.data <- read.csv(here("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#8-13_Data.csv"))
batch4.data <- read.csv(here("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#14-21_Data.csv"))
batch5.data <- read.csv(here("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/Batch#22-35_Data.csv"))

#Sample Info file for renaming
sampleinfo <- read.csv(here("data/metabolome/Original Raw Data Files/Data File Used in CleanUp/sampleinfo.csv"))

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

#Each Individual Batch
#Reorder each batch
#Reorder column names based on sample ID
batch1.data <- batch1.data[ , c(1, 1 + mixedorder(names(batch1.data[ , -1])))]
batch2.data <- batch2.data[ , c(1, 1 + mixedorder(names(batch2.data[ , -1])))]
batch3.data <- batch3.data[ , c(1, 1 + mixedorder(names(batch3.data[ , -1])))]
batch4.data <- batch4.data[ , c(1, 1 + mixedorder(names(batch4.data[ , -1])))]
batch5.data <- batch5.data[ , c(1, 1 + mixedorder(names(batch5.data[ , -1])))]

#Reorder rows based off alphabet
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


#Merging all batches together
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


# Set-Up
dat <- log10(dat)
dat <- as.data.frame(t(dat))
panel <- colnames(dat) # names of all metabolites on panel
dat$run_order <- samps$new_run_order
dat$run <- samps$run
dat$group <- samps$group
dat$grade <- samps$graded_hemolysis
dat$prep <- samps$prep_batch
dat$group <- factor(dat$group)

################################################
save.image(file=here('data/metabolome/environment.1.RData')) # save environment
load(here('data/metabolome/environment.1.RData'))


save(file=here('data/metabolome/logged.data'), dat, samps, panel) # save mz data and sample info




rm(list=ls())

################################################################################################
## quality control, normalization, imputation
################################################################################################
load(here('data/metabolome/logged.data'))

# Exclude Hemolyzed Samples (remove score =4)--------------------------------------------
dat <- dat[dat$group == 'Sample', ]
dat <- dat[dat$grade != 4, ]
dat <- dat[, c(1:361)]
samps <- samps[samps$group == 'Sample',]
samps <- samps[samps$graded_hemolysis != 4, ]


# Removing the most- missing metabolites ----------------------------------------------
table(colSums(is.na(dat))) # NA count per metabolite, 108 complete metabolites, 180 completely missing metabolites
par(mfrow=c(1,1))
hist(colSums(is.na(dat)/nrow(dat)), border=0, col='grey40', xlab='missingness (0 to 1)', main='', 30)

K <-print(round(nrow(dat) *0.1)) # 10% of 1274 samples = 127
dat <- dat[ ,colSums(is.na(dat)) <= K] # 10% missingness threshold, no more than 10% missing

dogmzs <- colnames(dat) # these are the dog metabolites that we will work with


# plots to show batch (prep) and LCMS run order effects:

par(mfrow=c(5,4), mar=c(4,4,1,1)+0.5)
for (i in 1:length(dogmzs)) {
  if(!all(is.na(dat[,dogmzs[i]]))) {
    plot(1:1274, dat[, i], pch=20, cex=0.5,
         main = colnames(dat)[i],
         xlab = "", ylab = "Abundance", xlim = c(1,1274), las = 2, col = samps$prep_batch) +
      abline(v = c(74,229,444,726)) }} # vertical lines delineate LCMS runs


######################################################################
######################################################################
# outlier (metabolites, sic) within batch.
# find w/ Grub Test: G = extremeX-mean(x)/var(x)
# replace with NA (for now, will be imputed back)

dat$run_order <- samps$new_run_order
dat$run <- samps$run
dat$prep <- samps$prep_batch

b1 <- dat[dat$prep == 1, ]
b2 <- dat[dat$prep == 2, ]
b3 <- dat[dat$prep == 3, ]

outliers1 <- data.frame(matrix(ncol = 5, nrow = 137))
dd1 <- vector()
dd2 <- vector()
for (i in c(1:49,51:87,89:137)){
  m <- grubbs.test(b1[, i], type = 10)
  outliers1[i, 2] <- m$alternative
  outliers1[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd1 <- rbind(dd1, i)
  }
  m <- grubbs.test(b1[, i], type = 10, opposite = T)
  outliers1[i, 4] <- m$alternative
  outliers1[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd2 <- rbind(dd2, i)
  }
  rownames(outliers1)[i] <- colnames(b1)[i]
  outliers1[i, 1] <- i
}

outliers2 <- data.frame(matrix(ncol = 5, nrow = 137))
dd3 <- vector()
dd4 <- vector()
for (i in c(1:49,51:87,89:137)){
  m <- grubbs.test(b2[, i], type = 10)
  outliers2[i, 2] <- m$alternative
  outliers2[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd3 <- rbind(dd3, i)
  }
  m <- grubbs.test(b2[, i], type = 10, opposite = T)
  outliers2[i, 4] <- m$alternative
  outliers2[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd4 <- rbind(dd4, i)
  }
  rownames(outliers2)[i] <- colnames(b2)[i]
  outliers2[i, 1] <- i
}

outliers3 <- data.frame(matrix(ncol = 5, nrow = 137))
dd5 <- vector()
dd6 <- vector()
for (i in c(1:49,51:87,89:137)){
  m <- grubbs.test(b3[, i], type = 10)
  outliers3[i, 2] <- m$alternative
  outliers3[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd5 <- rbind(dd5, i)
  }
  m <- grubbs.test(b3[, i], type = 10, opposite = T)
  outliers3[i, 4] <- m$alternative
  outliers3[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd6 <- rbind(dd6, i)
  }
  rownames(outliers3)[i] <- colnames(b3)[i]
  outliers3[i, 1] <- i
}


dd1 <- as.data.frame(dd1)
colnames(dd1) <- 'X1'
dd1 <- merge(dd1, outliers1[, c(1,2)], by = 'X1')
dd1 <- dd1 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd1$X1)){
    if (dd1$X1[j] == i){
      if (dd1$`max or min`[j] == 'highest'){
        x <- which.max(b1[,i])
        b1[x, i] <- NA
      }
      else{
        x <- which.min(b1[, i])
        b1[x, i] <- NA
      }
    }
    else{}
  }
}

dd2 <- as.data.frame(dd2)
colnames(dd2) <- 'X1'
dd2 <- merge(dd2, outliers1[, c(1,4)], by = 'X1')
dd2 <- dd2 %>% separate(X4, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd2$X1)){
    if (dd2$X1[j] == i){
      if (dd2$`max or min`[j] == 'highest'){
        x <- which.max(b1[,i])
        b1[x, i] <- NA
      }
      else{
        x <- which.min(b1[, i])
        b1[x, i] <- NA
      }
    }
    else{}
  }
}

dd3 <- as.data.frame(dd3)
colnames(dd3) <- 'X1'
dd3 <- merge(dd3, outliers2[, c(1,2)], by = 'X1')
dd3 <- dd3 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd3$X1)){
    if (dd3$X1[j] == i){
      if (dd3$`max or min`[j] == 'highest'){
        x <- which.max(b2[,i])
        b2[x, i] <- NA
      }
      else{
        x <- which.min(b2[, i])
        b2[x, i] <- NA
      }
    }
    else{}
  }
}

dd5 <- as.data.frame(dd5)
colnames(dd5) <- 'X1'
dd5 <- merge(dd5, outliers3[, c(1,2)], by = 'X1')
dd5 <- dd5 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd5$X1)){
    if (dd5$X1[j] == i){
      if (dd5$`max or min`[j] == 'highest'){
        x <- which.max(b3[,i])
        b3[x, i] <- NA
      }
      else{
        x <- which.min(b3[, i])
        b3[x, i] <- NA
      }
    }
    else{}
  }
}

batch1 <- rbind(b1,b2,b3)


b1 <- dat[dat$prep == 4, ]
b2 <- dat[dat$prep == 5, ]
b3 <- dat[dat$prep == 6, ]
b4 <- dat[dat$prep == 7, ]

outliers1 <- data.frame(matrix(ncol = 5, nrow = 137))
dd1 <- vector()
dd2 <- vector()
for (i in 1:137){
  m <- grubbs.test(b1[, i], type = 10)
  outliers1[i, 2] <- m$alternative
  outliers1[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd1 <- rbind(dd1, i)
  }
  m <- grubbs.test(b1[, i], type = 10, opposite = T)
  outliers1[i, 4] <- m$alternative
  outliers1[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd2 <- rbind(dd2, i)
  }
  rownames(outliers1)[i] <- colnames(b1)[i]
  outliers1[i, 1] <- i
}

outliers2 <- data.frame(matrix(ncol = 5, nrow = 137))
dd3 <- vector()
dd4 <- vector()
for (i in 1:137){
  m <- grubbs.test(b2[, i], type = 10)
  outliers2[i, 2] <- m$alternative
  outliers2[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd3 <- rbind(dd3, i)
  }
  m <- grubbs.test(b2[, i], type = 10, opposite = T)
  outliers2[i, 4] <- m$alternative
  outliers2[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd4 <- rbind(dd4, i)
  }
  rownames(outliers2)[i] <- colnames(b2)[i]
  outliers2[i, 1] <- i
}

outliers3 <- data.frame(matrix(ncol = 5, nrow = 137))
dd5 <- vector()
dd6 <- vector()
for (i in 1:137){
  m <- grubbs.test(b3[, i], type = 10)
  outliers3[i, 2] <- m$alternative
  outliers3[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd5 <- rbind(dd5, i)
  }
  m <- grubbs.test(b3[, i], type = 10, opposite = T)
  outliers3[i, 4] <- m$alternative
  outliers3[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd6 <- rbind(dd6, i)
  }
  rownames(outliers3)[i] <- colnames(b3)[i]
  outliers3[i, 1] <- i
}

outliers4 <- data.frame(matrix(ncol = 5, nrow = 137))
dd7 <- vector()
dd8 <- vector()
for (i in 1:137){
  m <- grubbs.test(b4[, i], type = 10)
  outliers4[i, 2] <- m$alternative
  outliers4[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd7 <- rbind(dd7, i)
  }
  m <- grubbs.test(b4[, i], type = 10, opposite = T)
  outliers4[i, 4] <- m$alternative
  outliers4[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd8 <- rbind(dd8, i)
  }
  rownames(outliers4)[i] <- colnames(b4)[i]
  outliers4[i, 1] <- i
}


dd1 <- as.data.frame(dd1)
colnames(dd1) <- 'X1'
dd1 <- merge(dd1, outliers1[, c(1,2)], by = 'X1')
dd1 <- dd1 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd1$X1)){
    if (dd1$X1[j] == i){
      if (dd1$`max or min`[j] == 'highest'){
        x <- which.max(b1[,i])
        b1[x, i] <- NA
      }
      else{
        x <- which.min(b1[, i])
        b1[x, i] <- NA
      }
    }
    else{}
  }
}

dd3 <- as.data.frame(dd3)
colnames(dd3) <- 'X1'
dd3 <- merge(dd3, outliers2[, c(1,2)], by = 'X1')
dd3 <- dd3 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd3$X1)){
    if (dd3$X1[j] == i){
      if (dd3$`max or min`[j] == 'highest'){
        x <- which.max(b2[,i])
        b2[x, i] <- NA
      }
      else{
        x <- which.min(b2[, i])
        b2[x, i] <- NA
      }
    }
    else{}
  }
}

dd5 <- as.data.frame(dd5)
colnames(dd5) <- 'X1'
dd5 <- merge(dd5, outliers3[, c(1,2)], by = 'X1')
dd5 <- dd5 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd5$X1)){
    if (dd5$X1[j] == i){
      if (dd5$`max or min`[j] == 'highest'){
        x <- which.max(b3[,i])
        b3[x, i] <- NA
      }
      else{
        x <- which.min(b3[, i])
        b3[x, i] <- NA
      }
    }
    else{}
  }
}

dd7 <- as.data.frame(dd7)
colnames(dd7) <- 'X1'
dd7 <- merge(dd7, outliers4[, c(1,2)], by = 'X1')
dd7 <- dd7 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd7$X1)){
    if (dd7$X1[j] == i){
      if (dd7$`max or min`[j] == 'highest'){
        x <- which.max(b4[,i])
        b4[x, i] <- NA
      }
      else{
        x <- which.min(b4[, i])
        b4[x, i] <- NA
      }
    }
    else{}
  }
}

batch2 <- rbind(b1,b2,b3,b4)




b1 <- dat[dat$prep == 8, ]
b2 <- dat[dat$prep == 9, ]
b3 <- dat[dat$prep == 10, ]
b4 <- dat[dat$prep == 11, ]
b5 <- dat[dat$prep == 12, ]
b6 <- dat[dat$prep == 13, ]

outliers1 <- data.frame(matrix(ncol = 5, nrow = 137))
dd1 <- vector()
dd2 <- vector()
for (i in 1:137){
  m <- grubbs.test(b1[, i], type = 10)
  outliers1[i, 2] <- m$alternative
  outliers1[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd1 <- rbind(dd1, i)
  }
  m <- grubbs.test(b1[, i], type = 10, opposite = T)
  outliers1[i, 4] <- m$alternative
  outliers1[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd2 <- rbind(dd2, i)
  }
  rownames(outliers1)[i] <- colnames(b1)[i]
  outliers1[i, 1] <- i
}

outliers2 <- data.frame(matrix(ncol = 5, nrow = 137))
dd3 <- vector()
dd4 <- vector()
for (i in 1:137){
  m <- grubbs.test(b2[, i], type = 10)
  outliers2[i, 2] <- m$alternative
  outliers2[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd3 <- rbind(dd3, i)
  }
  m <- grubbs.test(b2[, i], type = 10, opposite = T)
  outliers2[i, 4] <- m$alternative
  outliers2[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd4 <- rbind(dd4, i)
  }
  rownames(outliers2)[i] <- colnames(b2)[i]
  outliers2[i, 1] <- i
}

outliers3 <- data.frame(matrix(ncol = 5, nrow = 137))
dd5 <- vector()
dd6 <- vector()
for (i in 1:137){
  m <- grubbs.test(b3[, i], type = 10)
  outliers3[i, 2] <- m$alternative
  outliers3[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd5 <- rbind(dd5, i)
  }
  m <- grubbs.test(b3[, i], type = 10, opposite = T)
  outliers3[i, 4] <- m$alternative
  outliers3[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd6 <- rbind(dd6, i)
  }
  rownames(outliers3)[i] <- colnames(b3)[i]
  outliers3[i, 1] <- i
}

outliers4 <- data.frame(matrix(ncol = 5, nrow = 137))
dd7 <- vector()
dd8 <- vector()
for (i in 1:137){
  m <- grubbs.test(b4[, i], type = 10)
  outliers4[i, 2] <- m$alternative
  outliers4[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd7 <- rbind(dd7, i)
  }
  m <- grubbs.test(b4[, i], type = 10, opposite = T)
  outliers4[i, 4] <- m$alternative
  outliers4[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd8 <- rbind(dd8, i)
  }
  rownames(outliers4)[i] <- colnames(b4)[i]
  outliers4[i, 1] <- i
}

outliers5 <- data.frame(matrix(ncol = 5, nrow = 137))
dd9 <- vector()
dd10 <- vector()
for (i in 1:137){
  m <- grubbs.test(b5[, i], type = 10)
  outliers5[i, 2] <- m$alternative
  outliers5[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd9 <- rbind(dd9, i)
  }
  m <- grubbs.test(b5[, i], type = 10, opposite = T)
  outliers5[i, 4] <- m$alternative
  outliers5[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd10 <- rbind(dd10, i)
  }
  rownames(outliers5)[i] <- colnames(b5)[i]
  outliers5[i, 1] <- i
}

outliers6 <- data.frame(matrix(ncol = 5, nrow = 137))
dd11 <- vector()
dd12 <- vector()
for (i in 1:137){
  m <- grubbs.test(b6[, i], type = 10)
  outliers6[i, 2] <- m$alternative
  outliers6[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd11 <- rbind(dd11, i)
  }
  m <- grubbs.test(b6[, i], type = 10, opposite = T)
  outliers6[i, 4] <- m$alternative
  outliers6[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd12 <- rbind(dd12, i)
  }
  rownames(outliers6)[i] <- colnames(b6)[i]
  outliers6[i, 1] <- i
}



dd1 <- as.data.frame(dd1)
colnames(dd1) <- 'X1'
dd1 <- merge(dd1, outliers1[, c(1,2)], by = 'X1')
dd1 <- dd1 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd1$X1)){
    if (dd1$X1[j] == i){
      if (dd1$`max or min`[j] == 'highest'){
        x <- which.max(b1[,i])
        b1[x, i] <- NA
      }
      else{
        x <- which.min(b1[, i])
        b1[x, i] <- NA
      }
    }
    else{}
  }
}

dd3 <- as.data.frame(dd3)
colnames(dd3) <- 'X1'
dd3 <- merge(dd3, outliers2[, c(1,2)], by = 'X1')
dd3 <- dd3 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd3$X1)){
    if (dd3$X1[j] == i){
      if (dd3$`max or min`[j] == 'highest'){
        x <- which.max(b2[,i])
        b2[x, i] <- NA
      }
      else{
        x <- which.min(b2[, i])
        b2[x, i] <- NA
      }
    }
    else{}
  }
}

dd5 <- as.data.frame(dd5)
colnames(dd5) <- 'X1'
dd5 <- merge(dd5, outliers3[, c(1,2)], by = 'X1')
dd5 <- dd5 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd5$X1)){
    if (dd5$X1[j] == i){
      if (dd5$`max or min`[j] == 'highest'){
        x <- which.max(b3[,i])
        b3[x, i] <- NA
      }
      else{
        x <- which.min(b3[, i])
        b3[x, i] <- NA
      }
    }
    else{}
  }
}

dd7 <- as.data.frame(dd7)
colnames(dd7) <- 'X1'
dd7 <- merge(dd7, outliers4[, c(1,2)], by = 'X1')
dd7 <- dd7 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd7$X1)){
    if (dd7$X1[j] == i){
      if (dd7$`max or min`[j] == 'highest'){
        x <- which.max(b4[,i])
        b4[x, i] <- NA
      }
      else{
        x <- which.min(b4[, i])
        b4[x, i] <- NA
      }
    }
    else{}
  }
}

dd9 <- as.data.frame(dd9)
colnames(dd9) <- 'X1'
dd9 <- merge(dd9, outliers5[, c(1,2)], by = 'X1')
dd9 <- dd9 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd9$X1)){
    if (dd9$X1[j] == i){
      if (dd9$`max or min`[j] == 'highest'){
        x <- which.max(b5[,i])
        b5[x, i] <- NA
      }
      else{
        x <- which.min(b5[, i])
        b5[x, i] <- NA
      }
    }
    else{}
  }
}

dd10 <- as.data.frame(dd10)
colnames(dd10) <- 'X1'
dd10 <- merge(dd10, outliers5[, c(1,4)], by = 'X1')
dd10 <- dd10 %>% separate(X4, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd10$X1)){
    if (dd10$X1[j] == i){
      if (dd10$`max or min`[j] == 'highest'){
        x <- which.max(b5[,i])
        b5[x, i] <- NA
      }
      else{
        x <- which.min(b5[, i])
        b5[x, i] <- NA
      }
    }
    else{}
  }
}

dd11 <- as.data.frame(dd11)
colnames(dd11) <- 'X1'
dd11 <- merge(dd11, outliers6[, c(1,2)], by = 'X1')
dd11 <- dd11 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd11$X1)){
    if (dd11$X1[j] == i){
      if (dd11$`max or min`[j] == 'highest'){
        x <- which.max(b6[,i])
        b6[x, i] <- NA
      }
      else{
        x <- which.min(b6[, i])
        b6[x, i] <- NA
      }
    }
    else{}
  }
}

batch3 <- rbind(b1,b2,b3,b4,b5,b6)




b1 <- dat[dat$prep == 14, ]
b2 <- dat[dat$prep == 15, ]
b3 <- dat[dat$prep == 16, ]
b4 <- dat[dat$prep == 17, ]
b5 <- dat[dat$prep == 18, ]
b6 <- dat[dat$prep == 19, ]
b7 <- dat[dat$prep == 20, ]
b8 <- dat[dat$prep == 21, ]

outliers1 <- data.frame(matrix(ncol = 5, nrow = 137))
dd1 <- vector()
dd2 <- vector()
for (i in 1:137){
  m <- grubbs.test(b1[, i], type = 10)
  outliers1[i, 2] <- m$alternative
  outliers1[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd1 <- rbind(dd1, i)
  }
  m <- grubbs.test(b1[, i], type = 10, opposite = T)
  outliers1[i, 4] <- m$alternative
  outliers1[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd2 <- rbind(dd2, i)
  }
  rownames(outliers1)[i] <- colnames(b1)[i]
  outliers1[i, 1] <- i
}

outliers2 <- data.frame(matrix(ncol = 5, nrow = 137))
dd3 <- vector()
dd4 <- vector()
for (i in 1:137){
  m <- grubbs.test(b2[, i], type = 10)
  outliers2[i, 2] <- m$alternative
  outliers2[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd3 <- rbind(dd3, i)
  }
  m <- grubbs.test(b2[, i], type = 10, opposite = T)
  outliers2[i, 4] <- m$alternative
  outliers2[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd4 <- rbind(dd4, i)
  }
  rownames(outliers2)[i] <- colnames(b2)[i]
  outliers2[i, 1] <- i
}

outliers3 <- data.frame(matrix(ncol = 5, nrow = 137))
dd5 <- vector()
dd6 <- vector()
for (i in 1:137){
  m <- grubbs.test(b3[, i], type = 10)
  outliers3[i, 2] <- m$alternative
  outliers3[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd5 <- rbind(dd5, i)
  }
  m <- grubbs.test(b3[, i], type = 10, opposite = T)
  outliers3[i, 4] <- m$alternative
  outliers3[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd6 <- rbind(dd6, i)
  }
  rownames(outliers3)[i] <- colnames(b3)[i]
  outliers3[i, 1] <- i
}

outliers4 <- data.frame(matrix(ncol = 5, nrow = 137))
dd7 <- vector()
dd8 <- vector()
for (i in 1:137){
  m <- grubbs.test(b4[, i], type = 10)
  outliers4[i, 2] <- m$alternative
  outliers4[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd7 <- rbind(dd7, i)
  }
  m <- grubbs.test(b4[, i], type = 10, opposite = T)
  outliers4[i, 4] <- m$alternative
  outliers4[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd8 <- rbind(dd8, i)
  }
  rownames(outliers4)[i] <- colnames(b4)[i]
  outliers4[i, 1] <- i
}

outliers5 <- data.frame(matrix(ncol = 5, nrow = 137))
dd9 <- vector()
dd10 <- vector()
for (i in 1:137){
  m <- grubbs.test(b5[, i], type = 10)
  outliers5[i, 2] <- m$alternative
  outliers5[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd9 <- rbind(dd9, i)
  }
  m <- grubbs.test(b5[, i], type = 10, opposite = T)
  outliers5[i, 4] <- m$alternative
  outliers5[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd10 <- rbind(dd10, i)
  }
  rownames(outliers5)[i] <- colnames(b5)[i]
  outliers5[i, 1] <- i
}

outliers6 <- data.frame(matrix(ncol = 5, nrow = 137))
dd11 <- vector()
dd12 <- vector()
for (i in 1:137){
  m <- grubbs.test(b6[, i], type = 10)
  outliers6[i, 2] <- m$alternative
  outliers6[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd11 <- rbind(dd11, i)
  }
  m <- grubbs.test(b6[, i], type = 10, opposite = T)
  outliers6[i, 4] <- m$alternative
  outliers6[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd12 <- rbind(dd12, i)
  }
  rownames(outliers6)[i] <- colnames(b6)[i]
  outliers6[i, 1] <- i
}

outliers7 <- data.frame(matrix(ncol = 5, nrow = 137))
dd13 <- vector()
dd14 <- vector()
for (i in 1:137){
  m <- grubbs.test(b7[, i], type = 10)
  outliers7[i, 2] <- m$alternative
  outliers7[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd13 <- rbind(dd13, i)
  }
  m <- grubbs.test(b7[, i], type = 10, opposite = T)
  outliers7[i, 4] <- m$alternative
  outliers7[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd14 <- rbind(dd14, i)
  }
  rownames(outliers7)[i] <- colnames(b7)[i]
  outliers7[i, 1] <- i
}

outliers8 <- data.frame(matrix(ncol = 5, nrow = 137))
dd15 <- vector()
dd16 <- vector()
for (i in 1:137){
  m <- grubbs.test(b8[, i], type = 10)
  outliers8[i, 2] <- m$alternative
  outliers8[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd15 <- rbind(dd15, i)
  }
  m <- grubbs.test(b8[, i], type = 10, opposite = T)
  outliers8[i, 4] <- m$alternative
  outliers8[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd16 <- rbind(dd16, i)
  }
  rownames(outliers8)[i] <- colnames(b8)[i]
  outliers8[i, 1] <- i
}




dd1 <- as.data.frame(dd1)
colnames(dd1) <- 'X1'
dd1 <- merge(dd1, outliers1[, c(1,2)], by = 'X1')
dd1 <- dd1 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd1$X1)){
    if (dd1$X1[j] == i){
      if (dd1$`max or min`[j] == 'highest'){
        x <- which.max(b1[,i])
        b1[x, i] <- NA
      }
      else{
        x <- which.min(b1[, i])
        b1[x, i] <- NA
      }
    }
    else{}
  }
}

dd3 <- as.data.frame(dd3)
colnames(dd3) <- 'X1'
dd3 <- merge(dd3, outliers2[, c(1,2)], by = 'X1')
dd3 <- dd3 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd3$X1)){
    if (dd3$X1[j] == i){
      if (dd3$`max or min`[j] == 'highest'){
        x <- which.max(b2[,i])
        b2[x, i] <- NA
      }
      else{
        x <- which.min(b2[, i])
        b2[x, i] <- NA
      }
    }
    else{}
  }
}

dd5 <- as.data.frame(dd5)
colnames(dd5) <- 'X1'
dd5 <- merge(dd5, outliers3[, c(1,2)], by = 'X1')
dd5 <- dd5 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd5$X1)){
    if (dd5$X1[j] == i){
      if (dd5$`max or min`[j] == 'highest'){
        x <- which.max(b3[,i])
        b3[x, i] <- NA
      }
      else{
        x <- which.min(b3[, i])
        b3[x, i] <- NA
      }
    }
    else{}
  }
}

dd6 <- as.data.frame(dd6)
colnames(dd6) <- 'X1'
dd6 <- merge(dd6, outliers3[, c(1,4)], by = 'X1')
dd6 <- dd6 %>% separate(X4, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd6$X1)){
    if (dd6$X1[j] == i){
      if (dd6$`max or min`[j] == 'highest'){
        x <- which.max(b3[,i])
        b3[x, i] <- NA
      }
      else{
        x <- which.min(b3[, i])
        b3[x, i] <- NA
      }
    }
    else{}
  }
}

dd7 <- as.data.frame(dd7)
colnames(dd7) <- 'X1'
dd7 <- merge(dd7, outliers4[, c(1,2)], by = 'X1')
dd7 <- dd7 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd7$X1)){
    if (dd7$X1[j] == i){
      if (dd7$`max or min`[j] == 'highest'){
        x <- which.max(b4[,i])
        b4[x, i] <- NA
      }
      else{
        x <- which.min(b4[, i])
        b4[x, i] <- NA
      }
    }
    else{}
  }
}

dd9 <- as.data.frame(dd9)
colnames(dd9) <- 'X1'
dd9 <- merge(dd9, outliers5[, c(1,2)], by = 'X1')
dd9 <- dd9 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd9$X1)){
    if (dd9$X1[j] == i){
      if (dd9$`max or min`[j] == 'highest'){
        x <- which.max(b5[,i])
        b5[x, i] <- NA
      }
      else{
        x <- which.min(b5[, i])
        b5[x, i] <- NA
      }
    }
    else{}
  }
}

dd10 <- as.data.frame(dd10)
colnames(dd10) <- 'X1'
dd10 <- merge(dd10, outliers5[, c(1,4)], by = 'X1')
dd10 <- dd10 %>% separate(X4, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd10$X1)){
    if (dd10$X1[j] == i){
      if (dd10$`max or min`[j] == 'highest'){
        x <- which.max(b5[,i])
        b5[x, i] <- NA
      }
      else{
        x <- which.min(b5[, i])
        b5[x, i] <- NA
      }
    }
    else{}
  }
}

dd11 <- as.data.frame(dd11)
colnames(dd11) <- 'X1'
dd11 <- merge(dd11, outliers6[, c(1,2)], by = 'X1')
dd11 <- dd11 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd11$X1)){
    if (dd11$X1[j] == i){
      if (dd11$`max or min`[j] == 'highest'){
        x <- which.max(b6[,i])
        b6[x, i] <- NA
      }
      else{
        x <- which.min(b6[, i])
        b6[x, i] <- NA
      }
    }
    else{}
  }
}

dd13 <- as.data.frame(dd13)
colnames(dd13) <- 'X1'
dd13 <- merge(dd13, outliers7[, c(1,2)], by = 'X1')
dd13 <- dd13 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd13$X1)){
    if (dd13$X1[j] == i){
      if (dd13$`max or min`[j] == 'highest'){
        x <- which.max(b7[,i])
        b7[x, i] <- NA
      }
      else{
        x <- which.min(b7[, i])
        b7[x, i] <- NA
      }
    }
    else{}
  }
}

dd14 <- as.data.frame(dd14)
colnames(dd14) <- 'X1'
dd14 <- merge(dd14, outliers7[, c(1,4)], by = 'X1')
dd14 <- dd14 %>% separate(X4, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd14$X1)){
    if (dd14$X1[j] == i){
      if (dd14$`max or min`[j] == 'highest'){
        x <- which.max(b7[,i])
        b7[x, i] <- NA
      }
      else{
        x <- which.min(b7[, i])
        b7[x, i] <- NA
      }
    }
    else{}
  }
}

dd15 <- as.data.frame(dd15)
colnames(dd15) <- 'X1'
dd15 <- merge(dd15, outliers8[, c(1,2)], by = 'X1')
dd15 <- dd15 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd15$X1)){
    if (dd15$X1[j] == i){
      if (dd15$`max or min`[j] == 'highest'){
        x <- which.max(b8[,i])
        b8[x, i] <- NA
      }
      else{
        x <- which.min(b8[, i])
        b8[x, i] <- NA
      }
    }
    else{}
  }
}

batch4 <- rbind(b1,b2,b3,b4,b5,b6,b7,b8)




b1 <- dat[dat$prep == 22, ]
b2 <- dat[dat$prep == 23, ]
b3 <- dat[dat$prep == 24, ]
b4 <- dat[dat$prep == 25, ]
b5 <- dat[dat$prep == 26, ]
b6 <- dat[dat$prep == 27, ]
b7 <- dat[dat$prep == 28, ]
b8 <- dat[dat$prep == 29, ]
b9 <- dat[dat$prep == 30, ]
b10 <- dat[dat$prep == 31, ]
b11 <- dat[dat$prep == 32, ]
b12 <- dat[dat$prep == 33, ]
b13 <- dat[dat$prep == 34, ]
b14 <- dat[dat$prep == 35, ]

outliers1 <- data.frame(matrix(ncol = 5, nrow = 137))
dd1 <- vector()
dd2 <- vector()
for (i in 1:137){
  m <- grubbs.test(b1[, i], type = 10)
  outliers1[i, 2] <- m$alternative
  outliers1[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd1 <- rbind(dd1, i)
  }
  m <- grubbs.test(b1[, i], type = 10, opposite = T)
  outliers1[i, 4] <- m$alternative
  outliers1[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd2 <- rbind(dd2, i)
  }
  rownames(outliers1)[i] <- colnames(b1)[i]
  outliers1[i, 1] <- i
}

outliers2 <- data.frame(matrix(ncol = 5, nrow = 137))
dd3 <- vector()
dd4 <- vector()
for (i in 1:137){
  m <- grubbs.test(b2[, i], type = 10)
  outliers2[i, 2] <- m$alternative
  outliers2[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd3 <- rbind(dd3, i)
  }
  m <- grubbs.test(b2[, i], type = 10, opposite = T)
  outliers2[i, 4] <- m$alternative
  outliers2[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd4 <- rbind(dd4, i)
  }
  rownames(outliers2)[i] <- colnames(b2)[i]
  outliers2[i, 1] <- i
}

outliers3 <- data.frame(matrix(ncol = 5, nrow = 137))
dd5 <- vector()
dd6 <- vector()
for (i in 1:137){
  m <- grubbs.test(b3[, i], type = 10)
  outliers3[i, 2] <- m$alternative
  outliers3[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd5 <- rbind(dd5, i)
  }
  m <- grubbs.test(b3[, i], type = 10, opposite = T)
  outliers3[i, 4] <- m$alternative
  outliers3[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd6 <- rbind(dd6, i)
  }
  rownames(outliers3)[i] <- colnames(b3)[i]
  outliers3[i, 1] <- i
}

outliers4 <- data.frame(matrix(ncol = 5, nrow = 137))
dd7 <- vector()
dd8 <- vector()
for (i in 1:137){
  m <- grubbs.test(b4[, i], type = 10)
  outliers4[i, 2] <- m$alternative
  outliers4[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd7 <- rbind(dd7, i)
  }
  m <- grubbs.test(b4[, i], type = 10, opposite = T)
  outliers4[i, 4] <- m$alternative
  outliers4[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd8 <- rbind(dd8, i)
  }
  rownames(outliers4)[i] <- colnames(b4)[i]
  outliers4[i, 1] <- i
}

outliers5 <- data.frame(matrix(ncol = 5, nrow = 137))
dd9 <- vector()
dd10 <- vector()
for (i in 1:137){
  m <- grubbs.test(b5[, i], type = 10)
  outliers5[i, 2] <- m$alternative
  outliers5[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd9 <- rbind(dd9, i)
  }
  m <- grubbs.test(b5[, i], type = 10, opposite = T)
  outliers5[i, 4] <- m$alternative
  outliers5[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd10 <- rbind(dd10, i)
  }
  rownames(outliers5)[i] <- colnames(b5)[i]
  outliers5[i, 1] <- i
}

outliers6 <- data.frame(matrix(ncol = 5, nrow = 137))
dd11 <- vector()
dd12 <- vector()
for (i in 1:137){
  m <- grubbs.test(b6[, i], type = 10)
  outliers6[i, 2] <- m$alternative
  outliers6[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd11 <- rbind(dd11, i)
  }
  m <- grubbs.test(b6[, i], type = 10, opposite = T)
  outliers6[i, 4] <- m$alternative
  outliers6[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd12 <- rbind(dd12, i)
  }
  rownames(outliers6)[i] <- colnames(b6)[i]
  outliers6[i, 1] <- i
}

outliers7 <- data.frame(matrix(ncol = 5, nrow = 137))
dd13 <- vector()
dd14 <- vector()
for (i in 1:137){
  m <- grubbs.test(b7[, i], type = 10)
  outliers7[i, 2] <- m$alternative
  outliers7[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd13 <- rbind(dd13, i)
  }
  m <- grubbs.test(b7[, i], type = 10, opposite = T)
  outliers7[i, 4] <- m$alternative
  outliers7[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd14 <- rbind(dd14, i)
  }
  rownames(outliers7)[i] <- colnames(b7)[i]
  outliers7[i, 1] <- i
}

outliers8 <- data.frame(matrix(ncol = 5, nrow = 137))
dd15 <- vector()
dd16 <- vector()
for (i in 1:137){
  m <- grubbs.test(b8[, i], type = 10)
  outliers8[i, 2] <- m$alternative
  outliers8[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd15 <- rbind(dd15, i)
  }
  m <- grubbs.test(b8[, i], type = 10, opposite = T)
  outliers8[i, 4] <- m$alternative
  outliers8[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd16 <- rbind(dd16, i)
  }
  rownames(outliers8)[i] <- colnames(b8)[i]
  outliers8[i, 1] <- i
}

outliers9 <- data.frame(matrix(ncol = 5, nrow = 137))
dd17 <- vector()
dd18 <- vector()
for (i in 1:137){
  m <- grubbs.test(b9[, i], type = 10)
  outliers9[i, 2] <- m$alternative
  outliers9[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd17 <- rbind(dd17, i)
  }
  m <- grubbs.test(b9[, i], type = 10, opposite = T)
  outliers9[i, 4] <- m$alternative
  outliers9[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd18 <- rbind(dd18, i)
  }
  rownames(outliers9)[i] <- colnames(b9)[i]
  outliers9[i, 1] <- i
}

outliers10 <- data.frame(matrix(ncol = 5, nrow = 137))
dd19 <- vector()
dd20 <- vector()
for (i in 1:137){
  m <- grubbs.test(b10[, i], type = 10)
  outliers10[i, 2] <- m$alternative
  outliers10[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd19 <- rbind(dd19, i)
  }
  m <- grubbs.test(b10[, i], type = 10, opposite = T)
  outliers10[i, 4] <- m$alternative
  outliers10[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd20 <- rbind(dd20, i)
  }
  rownames(outliers10)[i] <- colnames(b10)[i]
  outliers10[i, 1] <- i
}

outliers11 <- data.frame(matrix(ncol = 5, nrow = 137))
dd21 <- vector()
dd22 <- vector()
for (i in 1:137){
  m <- grubbs.test(b11[, i], type = 10)
  outliers11[i, 2] <- m$alternative
  outliers11[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd21 <- rbind(dd21, i)
  }
  m <- grubbs.test(b11[, i], type = 10, opposite = T)
  outliers11[i, 4] <- m$alternative
  outliers11[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd22 <- rbind(dd22, i)
  }
  rownames(outliers11)[i] <- colnames(b11)[i]
  outliers11[i, 1] <- i
}

outliers12 <- data.frame(matrix(ncol = 5, nrow = 137))
dd23 <- vector()
dd24 <- vector()
for (i in 1:137){
  m <- grubbs.test(b12[, i], type = 10)
  outliers12[i, 2] <- m$alternative
  outliers12[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd23 <- rbind(dd23, i)
  }
  m <- grubbs.test(b12[, i], type = 10, opposite = T)
  outliers12[i, 4] <- m$alternative
  outliers12[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd24 <- rbind(dd24, i)
  }
  rownames(outliers12)[i] <- colnames(b12)[i]
  outliers12[i, 1] <- i
}

outliers13 <- data.frame(matrix(ncol = 5, nrow = 137))
dd25 <- vector()
dd26 <- vector()
for (i in 1:137){
  m <- grubbs.test(b13[, i], type = 10)
  outliers13[i, 2] <- m$alternative
  outliers13[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd25 <- rbind(dd25, i)
  }
  m <- grubbs.test(b13[, i], type = 10, opposite = T)
  outliers13[i, 4] <- m$alternative
  outliers13[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd26 <- rbind(dd26, i)
  }
  rownames(outliers13)[i] <- colnames(b13)[i]
  outliers13[i, 1] <- i
}

outliers14 <- data.frame(matrix(ncol = 5, nrow = 137))
dd27 <- vector()
dd28 <- vector()
for (i in 1:137){
  m <- grubbs.test(b14[, i], type = 10)
  outliers14[i, 2] <- m$alternative
  outliers14[i, 3] <- m$p.value
  if (m$p.value <= 0.05){
    dd27 <- rbind(dd27, i)
  }
  m <- grubbs.test(b14[, i], type = 10, opposite = T)
  outliers14[i, 4] <- m$alternative
  outliers14[i, 5] <- m$p.value
  if (m$p.value <= 0.05){
    dd28 <- rbind(dd28, i)
  }
  rownames(outliers14)[i] <- colnames(b14)[i]
  outliers14[i, 1] <- i
}




dd1 <- as.data.frame(dd1)
colnames(dd1) <- 'X1'
dd1 <- merge(dd1, outliers1[, c(1,2)], by = 'X1')
dd1 <- dd1 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd1$X1)){
    if (dd1$X1[j] == i){
      if (dd1$`max or min`[j] == 'highest'){
        x <- which.max(b1[,i])
        b1[x, i] <- NA
      }
      else{
        x <- which.min(b1[, i])
        b1[x, i] <- NA
      }
    }
    else{}
  }
}

dd3 <- as.data.frame(dd3)
colnames(dd3) <- 'X1'
dd3 <- merge(dd3, outliers2[, c(1,2)], by = 'X1')
dd3 <- dd3 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd3$X1)){
    if (dd3$X1[j] == i){
      if (dd3$`max or min`[j] == 'highest'){
        x <- which.max(b2[,i])
        b2[x, i] <- NA
      }
      else{
        x <- which.min(b2[, i])
        b2[x, i] <- NA
      }
    }
    else{}
  }
}

dd5 <- as.data.frame(dd5)
colnames(dd5) <- 'X1'
dd5 <- merge(dd5, outliers3[, c(1,2)], by = 'X1')
dd5 <- dd5 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd5$X1)){
    if (dd5$X1[j] == i){
      if (dd5$`max or min`[j] == 'highest'){
        x <- which.max(b3[,i])
        b3[x, i] <- NA
      }
      else{
        x <- which.min(b3[, i])
        b3[x, i] <- NA
      }
    }
    else{}
  }
}

dd7 <- as.data.frame(dd7)
colnames(dd7) <- 'X1'
dd7 <- merge(dd7, outliers4[, c(1,2)], by = 'X1')
dd7 <- dd7 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd7$X1)){
    if (dd7$X1[j] == i){
      if (dd7$`max or min`[j] == 'highest'){
        x <- which.max(b4[,i])
        b4[x, i] <- NA
      }
      else{
        x <- which.min(b4[, i])
        b4[x, i] <- NA
      }
    }
    else{}
  }
}

dd8 <- as.data.frame(dd8)
colnames(dd8) <- 'X1'
dd8 <- merge(dd8, outliers4[, c(1,4)], by = 'X1')
dd8 <- dd8 %>% separate(X4, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd8$X1)){
    if (dd8$X1[j] == i){
      if (dd8$`max or min`[j] == 'highest'){
        x <- which.max(b4[,i])
        b4[x, i] <- NA
      }
      else{
        x <- which.min(b4[, i])
        b4[x, i] <- NA
      }
    }
    else{}
  }
}

dd9 <- as.data.frame(dd9)
colnames(dd9) <- 'X1'
dd9 <- merge(dd9, outliers5[, c(1,2)], by = 'X1')
dd9 <- dd9 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd9$X1)){
    if (dd9$X1[j] == i){
      if (dd9$`max or min`[j] == 'highest'){
        x <- which.max(b5[,i])
        b5[x, i] <- NA
      }
      else{
        x <- which.min(b5[, i])
        b5[x, i] <- NA
      }
    }
    else{}
  }
}

dd11 <- as.data.frame(dd11)
colnames(dd11) <- 'X1'
dd11 <- merge(dd11, outliers6[, c(1,2)], by = 'X1')
dd11 <- dd11 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd11$X1)){
    if (dd11$X1[j] == i){
      if (dd11$`max or min`[j] == 'highest'){
        x <- which.max(b6[,i])
        b6[x, i] <- NA
      }
      else{
        x <- which.min(b6[, i])
        b6[x, i] <- NA
      }
    }
    else{}
  }
}

dd12 <- as.data.frame(dd12)
colnames(dd12) <- 'X1'
dd12 <- merge(dd12, outliers6[, c(1,4)], by = 'X1')
dd12 <- dd12 %>% separate(X4, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd12$X1)){
    if (dd12$X1[j] == i){
      if (dd12$`max or min`[j] == 'highest'){
        x <- which.max(b6[,i])
        b6[x, i] <- NA
      }
      else{
        x <- which.min(b6[, i])
        b6[x, i] <- NA
      }
    }
    else{}
  }
}

dd13 <- as.data.frame(dd13)
colnames(dd13) <- 'X1'
dd13 <- merge(dd13, outliers7[, c(1,2)], by = 'X1')
dd13 <- dd13 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd13$X1)){
    if (dd13$X1[j] == i){
      if (dd13$`max or min`[j] == 'highest'){
        x <- which.max(b7[,i])
        b7[x, i] <- NA
      }
      else{
        x <- which.min(b7[, i])
        b7[x, i] <- NA
      }
    }
    else{}
  }
}

dd15 <- as.data.frame(dd15)
colnames(dd15) <- 'X1'
dd15 <- merge(dd15, outliers8[, c(1,2)], by = 'X1')
dd15 <- dd15 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd15$X1)){
    if (dd15$X1[j] == i){
      if (dd15$`max or min`[j] == 'highest'){
        x <- which.max(b8[,i])
        b8[x, i] <- NA
      }
      else{
        x <- which.min(b8[, i])
        b8[x, i] <- NA
      }
    }
    else{}
  }
}

dd17 <- as.data.frame(dd17)
colnames(dd17) <- 'X1'
dd17 <- merge(dd17, outliers9[, c(1,2)], by = 'X1')
dd17 <- dd17 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd17$X1)){
    if (dd17$X1[j] == i){
      if (dd17$`max or min`[j] == 'highest'){
        x <- which.max(b9[,i])
        b9[x, i] <- NA
      }
      else{
        x <- which.min(b9[, i])
        b9[x, i] <- NA
      }
    }
    else{}
  }
}

dd18 <- as.data.frame(dd18)
colnames(dd18) <- 'X1'
dd18 <- merge(dd18, outliers9[, c(1,4)], by = 'X1')
dd18 <- dd18 %>% separate(X4, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd18$X1)){
    if (dd18$X1[j] == i){
      if (dd18$`max or min`[j] == 'highest'){
        x <- which.max(b9[,i])
        b9[x, i] <- NA
      }
      else{
        x <- which.min(b9[, i])
        b9[x, i] <- NA
      }
    }
    else{}
  }
}

dd19 <- as.data.frame(dd19)
colnames(dd19) <- 'X1'
dd19 <- merge(dd19, outliers10[, c(1,2)], by = 'X1')
dd19 <- dd19 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd19$X1)){
    if (dd19$X1[j] == i){
      if (dd19$`max or min`[j] == 'highest'){
        x <- which.max(b10[,i])
        b10[x, i] <- NA
      }
      else{
        x <- which.min(b10[, i])
        b10[x, i] <- NA
      }
    }
    else{}
  }
}

dd21 <- as.data.frame(dd21)
colnames(dd21) <- 'X1'
dd21 <- merge(dd21, outliers11[, c(1,2)], by = 'X1')
dd21 <- dd21 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd21$X1)){
    if (dd21$X1[j] == i){
      if (dd21$`max or min`[j] == 'highest'){
        x <- which.max(b11[,i])
        b11[x, i] <- NA
      }
      else{
        x <- which.min(b11[, i])
        b11[x, i] <- NA
      }
    }
    else{}
  }
}

dd22 <- as.data.frame(dd22)
colnames(dd22) <- 'X1'
dd22 <- merge(dd22, outliers11[, c(1,4)], by = 'X1')
dd22 <- dd22 %>% separate(X4, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd22$X1)){
    if (dd22$X1[j] == i){
      if (dd22$`max or min`[j] == 'highest'){
        x <- which.max(b11[,i])
        b11[x, i] <- NA
      }
      else{
        x <- which.min(b11[, i])
        b11[x, i] <- NA
      }
    }
    else{}
  }
}

dd23 <- as.data.frame(dd23)
colnames(dd23) <- 'X1'
dd23 <- merge(dd23, outliers12[, c(1,2)], by = 'X1')
dd23 <- dd23 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd23$X1)){
    if (dd23$X1[j] == i){
      if (dd23$`max or min`[j] == 'highest'){
        x <- which.max(b12[,i])
        b12[x, i] <- NA
      }
      else{
        x <- which.min(b12[, i])
        b12[x, i] <- NA
      }
    }
    else{}
  }
}

dd24 <- as.data.frame(dd24)
colnames(dd24) <- 'X1'
dd24 <- merge(dd24, outliers12[, c(1,4)], by = 'X1')
dd24 <- dd24 %>% separate(X4, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd24$X1)){
    if (dd24$X1[j] == i){
      if (dd24$`max or min`[j] == 'highest'){
        x <- which.max(b12[,i])
        b12[x, i] <- NA
      }
      else{
        x <- which.min(b12[, i])
        b12[x, i] <- NA
      }
    }
    else{}
  }
}

dd25 <- as.data.frame(dd25)
colnames(dd25) <- 'X1'
dd25 <- merge(dd25, outliers13[, c(1,2)], by = 'X1')
dd25 <- dd25 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd25$X1)){
    if (dd25$X1[j] == i){
      if (dd25$`max or min`[j] == 'highest'){
        x <- which.max(b13[,i])
        b13[x, i] <- NA
      }
      else{
        x <- which.min(b13[, i])
        b13[x, i] <- NA
      }
    }
    else{}
  }
}

dd26 <- as.data.frame(dd26)
colnames(dd26) <- 'X1'
dd26 <- merge(dd26, outliers13[, c(1,4)], by = 'X1')
dd26 <- dd26 %>% separate(X4, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd26$X1)){
    if (dd26$X1[j] == i){
      if (dd26$`max or min`[j] == 'highest'){
        x <- which.max(b13[,i])
        b13[x, i] <- NA
      }
      else{
        x <- which.min(b13[, i])
        b13[x, i] <- NA
      }
    }
    else{}
  }
}

dd27 <- as.data.frame(dd27)
colnames(dd27) <- 'X1'
dd27 <- merge(dd27, outliers14[, c(1,2)], by = 'X1')
dd27 <- dd27 %>% separate(X2, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd27$X1)){
    if (dd27$X1[j] == i){
      if (dd27$`max or min`[j] == 'highest'){
        x <- which.max(b14[,i])
        b14[x, i] <- NA
      }
      else{
        x <- which.min(b14[, i])
        b14[x, i] <- NA
      }
    }
    else{}
  }
}

dd28 <- as.data.frame(dd28)
colnames(dd28) <- 'X1'
dd28 <- merge(dd28, outliers14[, c(1,4)], by = 'X1')
dd28 <- dd28 %>% separate(X4, "max or min", " ")
for (i in 1:137){
  for (j in 1:length(dd28$X1)){
    if (dd28$X1[j] == i){
      if (dd28$`max or min`[j] == 'highest'){
        x <- which.max(b14[,i])
        b14[x, i] <- NA
      }
      else{
        x <- which.min(b14[, i])
        b14[x, i] <- NA
      }
    }
    else{}
  }
}

batch5 <- rbind(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14)


save.image(file=here('data/metabolome/environment.2.RData')) # save environment
######################################################################
######################################################################

load(here('data/metabolome/environment.2.RData'))

#Redefining Samples
sample1 <- samps[samps$run == 1, ]
sample2 <- samps[samps$run == 2, ]
sample3 <- samps[samps$run == 3, ]
sample4 <- samps[samps$run == 4, ]
sample5 <- samps[samps$run == 5, ]

batch1$prep <- as.factor(batch1$prep)
batch2$prep <- as.factor(batch2$prep)
batch3$prep <- as.factor(batch3$prep)
batch4$prep <- as.factor(batch4$prep)
batch5$prep <- as.factor(batch5$prep)

batch1$new_SampleID <- sample1$new_SampleID
batch2$new_SampleID <- sample2$new_SampleID
batch3$new_SampleID <- sample3$new_SampleID
batch4$new_SampleID <- sample4$new_SampleID
batch5$new_SampleID <- sample5$new_SampleID

#Normalizing by Prep Batch
batch1[, c(1:137,140)] <- batch1[, c(1:137,140)] %>% 
  group_by(prep) %>%
  mutate_if(is.numeric, scale)
batch1 <- as.data.frame(batch1)
rownames(batch1) <- batch1$new_SampleID

batch2[, c(1:137,140)] <- batch2[, c(1:137,140)] %>% 
  group_by(prep) %>%
  mutate_if(is.numeric, scale)
batch2 <- as.data.frame(batch2)
rownames(batch2) <- batch2$new_SampleID

batch3[, c(1:137,140)] <- batch3[, c(1:137,140)] %>% 
  group_by(prep) %>%
  mutate_if(is.numeric, scale)
batch3 <- as.data.frame(batch3)
rownames(batch3) <- batch3$new_SampleID

batch4[, c(1:137,140)] <- batch4[, c(1:137,140)] %>% 
  group_by(prep) %>%
  mutate_if(is.numeric, scale)
batch4 <- as.data.frame(batch4)
rownames(batch4) <- batch4$new_SampleID

batch5[, c(1:137,140)] <- batch5[, c(1:137,140)] %>% 
  group_by(prep) %>%
  mutate_if(is.numeric, scale)
batch5 <- as.data.frame(batch5)
rownames(batch5) <- batch5$new_SampleID


postOut <- rbind(batch1, batch2, batch3, batch4, batch5)

save(postOut, dat, samps, dogmzs, file='post_outliers')




rm(list=ls())
############################################################################
load('post_outliers')

dat$prep <- samps$prep_batch

par(mfrow=c(5,2), mar=c(4,4,1,1)+0.5)
for (i in 1:137) {
  if(!all(is.na(dat[, i]))) {
    plot(1:1274, dat[, i],
         main = colnames(dat)[i],
         xlab = "", ylab = "Abundance", xlim = c(1,1274), las = 2, col = dat$prep) +
      abline(v = c(74,229,444,726))
  }
  if(!all(is.na(postOut[, i]))) {
    plot(1:1274, postOut[, i],
         main = colnames(postOut)[i],
         xlab = "", ylab = "Abundance", xlim = c(1,1274), las = 2, col = postOut$prep) +
      abline(v = c(74,229,444,726))
  }
} 


## batchwise run order adjustment
runAdj <- postOut

for (i in c(1:49,51:87,89:137)){
  for (j in 1:35){
    m <- lm(postOut[postOut$prep == j, i] ~ run_order, postOut[postOut$prep == j, ],
            na.action = na.exclude)
    runAdj[runAdj$prep == j, i] <- resid(m)
  }
}

for (i in c(50,88)){
  for (j in 4:35){
    m <- lm(postOut[postOut$prep == j, i] ~ run_order, postOut[postOut$prep == j, ],
            na.action = na.exclude)
    runAdj[runAdj$prep == j, i] <- resid(m)
  }
}


par(mfrow=c(5,3), mar=c(4,4,1,1)+0.5)
for (i in 1:137) {
  if(!all(is.na(dat[, i]))) {
    plot(1:1274, dat[, i],
         main = colnames(dat)[i],
         xlab = "", ylab = "Abundance", xlim = c(1,1274), las = 2, col = dat$prep) +
      abline(v = c(74,229,444,726))
  }
  if(!all(is.na(postOut[, i]))) {
    plot(1:1274, postOut[, i],
         main = colnames(postOut)[i],
         xlab = "", ylab = "Abundance", xlim = c(1,1274), las = 2, col = postOut$prep) +
      abline(v = c(74,229,444,726))
  }
  if(!all(is.na(runAdj[, i]))) {
    plot(1:1274, runAdj[, i],
         main = colnames(runAdj)[i],
         xlab = "", ylab = "Abundance", xlim = c(1,1274), las = 2, col = runAdj$prep) +
      abline(v = c(74,229,444,726))
  }
} 


# Imputation --------------------------------------------------------------
runAdj <- as.matrix(runAdj[, dogmzs])
runAdj.impute <- (impute.knn(runAdj))$data
runAdj.impute <- as.data.frame(t(runAdj.impute))
dat.colnames <- samps$elab_id[match(names(runAdj.impute), samps$new_SampleID)]
names(runAdj.impute)[!is.na(dat.colnames)] <- dat.colnames[!is.na(dat.colnames)]
runAdj.impute <- runAdj.impute
runAdj.impute <- runAdj.impute[, mixedorder(colnames(runAdj.impute))]

dat <- as.data.frame(t(scale(runAdj.impute)))
dat$elab_id <- rownames(dat)
dat <- merge(samps, dat)

dat$dog_id <- as.character(dat$dog_id)
dat$elab_id <- as.character(dat$elab_id)
dat$run_order <- as.numeric(dat$run_order)

saveRDS(dat, here("data/metabolome/ProcessedData/TechnicalEffectsRemoved.RDS"))
###############################
rm(list=ls())


      