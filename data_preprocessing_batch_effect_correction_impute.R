library(readxl)
library(stringr)
library(ggplot2)
library(missForest)
library(vsn)
library(sva)

setwd("C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/Dataset 1-from Marco for ARDS")

# read raw data
data <- read_excel("011420-and-080919-and-052521-data-combined-Ji-BAL-EV-quant_proteins - Copy.xlsx", )
gene <- data$`Gene Name`
accession <- data$Accession
data.info <- data[,1:10]
data <- data[,-c(1:10)]

# count missing for control and patient
missing.gene.control <- rowSums(is.na(data[,1:49]))
missing.gene.patient <- rowSums(is.na(data[,50:106]))

# keep genes expressed more than 80% in both groups
rate <- 0.2
index.remove.control <- which(missing.gene.control > 49*rate)
index.remove.patient <- which(missing.gene.patient > 57*rate)
index.remove <- union(index.remove.control, index.remove.patient)
index.keep <- setdiff(1:dim(data)[1], index.remove)
data <- data[index.keep,]
saveRDS(index.keep, 'index.keep.rds')

# save kept gene
gene.keep <- gene[index.keep]
saveRDS(gene.keep, 'gene.keep.rds')

# save kept accession
accession.keep <- accession[index.keep]
saveRDS(accession.keep, 'accession.keep.rds')

# mean expression per sample
sample.mean <- colMeans(data, na.rm = T)

# order samples by batch
order <- read_excel("order.xlsx", col_names = T)
sample <- colnames(data)
sample.short <- sapply(str_split(sample, " - "),tail,1)
time <- match(sample.short, order$sample)

# save sample order
saveRDS(time, 'time.rds')

##############################################################
# vsn normalization
##############################################################
set.seed(2023)
data.norm <- justvsn(as.matrix(data))
dim(data.norm)

# combine raw info
data.norm.info <- cbind(data.info[index.keep,], data.norm)
write.csv(data.norm.info, 'data.norm.raw.csv', row.names = F, na = '')

##############################################################
# Random Forest imputation
##############################################################
set.seed(2023)
doParallel::registerDoParallel(cores = 12) # set based on number of CPU cores
doRNG::registerDoRNG(seed = 123)
system.time({
  result.impute <- missForest(t(data.norm), verbose = T, parallelize = 'variables')
})
data.norm.impute <- result.impute$ximp; dim(data.norm.impute)
saveRDS(data.norm.impute, 'data.norm.impute.80.rds')

# combine raw info
data.norm.impute <- readRDS('data.norm.impute.80.rds'); dim(data.norm.impute)
data.norm.impute <- t(data.norm.impute); dim(data.norm.impute)
data.norm.impute.info <- cbind(data.info[index.keep,], data.norm.impute)
write.csv(data.norm.impute.info, 'data.norm.impute.raw.csv', row.names = F, na = '')

data.complete <- t(readRDS('data.norm.impute.80.rds')); dim(data.complete)

##############################################################
# Batch effect correction
##############################################################
#data.norm.impute <- readRDS('data.norm.impute.80.rds')
dim(data.norm.impute)
batch <- c(rep(2,3), rep(3,29), rep(4,17), rep(1,3), rep(2,6), rep(3,31), rep(4,17))

set.seed(2023)
data.complete <- ComBat(dat=data.norm.impute, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
dim(data.complete)
saveRDS(data.complete, 'data.complete.80.rds')

#data.complete <- readRDS('data.complete.80.rds')
data.complete.info <- cbind(data.info[index.keep,], data.complete)
write.csv(data.complete.info, 'data.complete.raw.csv', row.names = F, na = '')

##############################################################
# PCA Dimension Reduction
##############################################################
data.complete <- readRDS('data.complete.80.rds'); dim(data.complete)
pr.out <- prcomp(t(data.complete), center = T, scale. = T)
data.pca <- pr.out$x; dim(data.pca)
PVE.matrix <- summary(pr.out)$importance; PVE.matrix
PVE <- PVE.matrix[2,]
plot(PVE, xlab='Principle Components', ylab='Proportion of Variance Explained', cex.lab=1.5, pch = 19)

# plot first two PCs for all
color <- c(rep('red',3), rep('green',9), rep('blue',60), rep('orange',34))
plot(data.pca[order(time),2]~data.pca[order(time),1], col = color, pch = 19, main='First Two PCs after Batech Effect Correction (All)',
     xlab='PC1 (19.3%)', ylab='PC2 (6.8%)')

# plot for patients
color <- c(rep('red',3), rep('green',6), rep('blue',31), rep('orange',17))
plot(data.pca[50:106,2]~data.pca[50:106,1], col = color, pch = 19, main='First Two PCs after Batech Effect Correction (Patients)',
     xlab='PC1 (19.3%)', ylab='PC2 (6.8%)')

# plot for control
color <- c(rep('green',3), rep('blue',29), rep('orange',17))
plot(data.pca[1:49,2]~data.pca[1:49,1], col = color, pch = 19, main='First Two PCs after Batech Effect Correction (Control)',
     xlab='PC1 (19.3%)', ylab='PC2 (6.8%)')

##############################################################
# Sample mean across batch
##############################################################
# mean expression per sample
sample.mean <- colMeans(data.complete, na.rm = T)

# plot mean abundance across sample time (median result is very similar)
plot(sample.mean[order(time)], pch = 19, ylab = 'Mean Abundance', xlab = 'Order', main = 'All (After Batch Effect Correction)'); abline(v=c(3.5,12.5,72.5))

# plot patients across sample time
plot(sample.mean[50:106], pch = 19, ylab = 'Mean Abundance', xlab = 'Order', main = 'Patients (After Batch Effect Correction)'); abline(v=c(3.5,9.5,40.5))

# plot control across sample time
plot(sample.mean[1:49], pch = 19, ylab = 'Mean Abundance', xlab = 'Order', main = 'Control (After Batch Effect Correction)'); abline(v=c(3.5,32.5))

##############################################################
# correlations within batch 
##############################################################
order$date <- as.character(order$date)
date <- unique(order$date)

# sample per batch
sample.batch <- list()
for(d in date){
  sample.batch <- append(sample.batch, list(order$sample[order$date == d]))
}

# within sample cor
cor.within <- c()
for(batch in sample.batch){
  data.batch <- data.complete[,which(sample.short%in%batch)]
  cor.batch <- unique(as.vector(cor(data.batch, use = 'complete.obs')))
  cor.batch <- cor.batch[cor.batch!=1]
  cor.within <- c(cor.within, cor.batch)
}

##############################################################
# correlations between batches
##############################################################
comb <- combn(1:length(sample.batch), 2)

# between batch samples
sample.between <- data.frame()
for(i in 1:ncol(comb)){
  sample.between <- rbind(sample.between,expand.grid(sample.batch[[comb[,i][1]]], sample.batch[[comb[,i][2]]],stringsAsFactors = F))
}

# between batch cor
cor.between <- c()
for(i in 1:nrow(sample.between)){
  sample.1 <- sample.between[i,]$Var1
  sample.2 <- sample.between[i,]$Var2
  index.1 <- which(sample.1==sample.short)
  index.2 <- which(sample.2==sample.short)
  cor.between <- c(cor.between, as.numeric(cor(data.complete[,index.1], data.complete[,index.2], use = 'complete.obs')))
}

boxplot(cor.within, cor.between, names=c('Within Batch','Between Batches'), ylab='Correlation', main='Sample Pearson Correlation (Batch Effect Correction)')

##############################################################
# sample distribution
##############################################################
# boxplot for all
color <- c(rep('red',3), rep('green',9), rep('blue',60), rep('orange',34))
boxplot(data.complete[,order(time)], outline = F, col = color, whisklty = 1, xaxt = "n",
        main="Abundance Distribution of All Samples (After Batch Effect Correction)", ylab="Abundance", xlab="Order")

# boxplot for patients
color <- c(rep('red',3), rep('green',6), rep('blue',31), rep('orange',17))
boxplot(data.complete[,50:106], outline = F,  col = color, whisklty = 1, xaxt = "n",
        main="Abundance Distribution of Patients ((After Batch Effect Correction))", ylab="Abundance", xlab="Order")

# boxplot for control
color <- c(rep('green',3), rep('blue',29), rep('orange',17))
boxplot(data.complete[,1:49], outline = F,  col = color, whisklty = 1, xaxt = "n",
        main="Abundance Distribution of Control (After Batch Effect Correction)", ylab="Abundance", xlab="Order")








