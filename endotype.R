library(readxl)
library(stringr)
library(ggplot2)
library(missForest)
library(vsn)
library(sva)

setwd("C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/Dataset 1-from Marco for ARDS")

# read clustering result
#cluster <- readRDS('kmeans.imp.outcome.rds')
cluster <- readRDS('kmeans.biomarker.Aspir.Pneumo.Sepsis.rds')
#cluster <- readRDS('kmeans.biomarker.death.shock.rds')
#cluster <- readRDS('kmeans.imp.Etiology.rds')
#cluster <- readRDS('kmeans.imp.all.rds')
#cluster <- readRDS('kmeans.imp.severity.rds')
#cluster <- readRDS('kmeans.imp.baseline.2.rds')
ptid <- read_excel("ptid.xlsx")

# assign cluster to ptid
ptid.cluster <- data.frame()
for(i in 1:length(cluster)){
  ptid.cluster <- rbind(ptid.cluster, cbind(cluster[[i]],rep(i, length(cluster[[i]]))))
}

colnames(ptid.cluster) <- c('index', 'cluster')
ptid.cluster <- ptid.cluster[order(ptid.cluster$index),]
ptid <- cbind(ptid, ptid.cluster$cluster)
colnames(ptid)[4] <- 'cluster'
ptid <- ptid[,c(2,4)]
table(ptid$cluster)

# creat ptid+endotype table
# ptid[ptid$cluster==1, 'cluster'] <- 'E1'
# ptid[ptid$cluster==2, 'cluster'] <- 'E2'
# write.csv(ptid, 'ptid_E1_E2.csv', row.names = F)
# 
# ptid[ptid$cluster==1, 'cluster'] <- 'E3'
# ptid[ptid$cluster==2, 'cluster'] <- 'E4'
# write.csv(ptid, 'ptid_E3_E4.csv', row.names = F)

# ptid[ptid$cluster==1, 'cluster'] <- 'E8'
# ptid[ptid$cluster==2, 'cluster'] <- 'E9'
# table(ptid$cluster)
# write.csv(ptid, 'ptid_E8_E9.csv', row.names = F)

#ptid[ptid$cluster==1, 'cluster'] <- 'E2\''
#ptid[ptid$cluster==2, 'cluster'] <- 'E1\''
#write.csv(ptid, 'ptid_E1_E2_biomarker.csv', row.names = F)

#######################################################
# baseline
#######################################################
# assign cluster to clinical data
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 1)
clinical[clinical=='.'] <- NA
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)
clinical$gender <- ifelse(clinical$gender==2, 1, 0)
clinical$Race <- ifelse(clinical$Race=="3", 1, 0)

# continuous variables
var.cont <- c('age', 'height', 'weight', 'pbw', 'bmi', 'packyr', 'quads', 'intubdt', 'pao2screen', 'fio2screen','qualpfdt','critdt',
              "hrate","sysbp","diabp","cvp","map",'temp')
p.wil <- c()
p.t <- c()
mean.1 <- c(); mean.2 <- c()
median.1 <- c(); median.2 <- c()

for(v in var.cont){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]) 
  var.2 <- as.numeric(clinical[clinical$cluster==2,v]) 
  
  w <- wilcox.test(var.1, var.2, alternative = "two.sided", na.action='na.omit')
  t <- t.test(var.1, var.2, alternative = "two.sided", na.action='na.omit')
  
  p.wil <- c(p.wil, w$p.value)
  p.t <- c(p.t, t$p.value)
  
  mean.1 <- c(mean.1, mean(var.1, na.rm=T))
  mean.2 <- c(mean.2, mean(var.2, na.rm=T))
  
  median.1 <- c(median.1, median(var.1, na.rm=T))
  median.2 <- c(median.2, median(var.2, na.rm=T))
}

result.endotype.cont <- data.frame(var=var.cont, p.wil, p.t, mean.1, mean.2, median.1, median.2)


# categorical variables
var.cat <- c('gender','ethnic','Race','vaso','smoker','cursmoker', 'death90', 'intfeed')

p.f <- c(); p.chi <- c()
freq.1 <- c(); freq.2 <- c()

for(v in var.cat){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v])
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  tab <- clinical[clinical$cluster!=3,c(v,'cluster')]
  t <- table(tab[,v], tab[,'cluster'])
  
  f <- fisher.test(t)
  chi <- chisq.test(t)
  
  p.f <- c(p.f, f$p.value)
  p.chi <- c(p.chi, chi$p.value)
  
  freq.1 <- c(freq.1, mean(var.1, na.rm=T))
  freq.2 <- c(freq.2, mean(var.2, na.rm=T))
}

result.endotype.cat <- data.frame(var=var.cat, p.f, p.chi, freq.1, freq.2)


#######################################################
# Etiology
#######################################################
# read clinical data
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 2)
clinical[clinical=='.'] <- NA
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)

# reorder to match protein data
clinical <- clinical[match(ptid$ptid, clinical$ptid),]

# categorical variable
var.cat <- c("Trauma", "Sepsis", "Transf", "Aspir", "Pneumo")

# trauma to binary
table(clinical$Trauma)
clinical$Trauma <- ifelse(clinical$Trauma==0, 0, 1)
table(clinical$Trauma)

#  transf to binary
table(clinical$Transf)
clinical$Transf <- ifelse(clinical$Transf==0, 0, 1)
table(clinical$Transf)

#  Aspir to binary
table(clinical$Aspir)
clinical$Aspir <- ifelse(clinical$Aspir==0, 0, 1)
table(clinical$Aspir)

#  Pneumo to binary
table(clinical$Pneumo)
clinical$Pneumo <- ifelse(clinical$Pneumo==0, 0, 1)
table(clinical$Pneumo)

#  Sepsis to binary
table(clinical$Sepsis)
clinical$Sepsis <- ifelse(clinical$Sepsis==0, 0, 1)
table(clinical$Sepsis)

p.f <- c()
p.chi <- c()
freq.1 <- c(); freq.2 <- c()

for(v in var.cat){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v])
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  tab <- clinical[clinical$cluster!=3,c(v,'cluster')]
  t <- table(tab[,v], tab[,'cluster'])
  
  f <- fisher.test(t)
  chi <- chisq.test(t)
  
  p.f <- c(p.f, f$p.value)
  p.chi <- c(p.chi, chi$p.value)
  
  freq.1 <- c(freq.1, mean(var.1, na.rm=T))
  freq.2 <- c(freq.2, mean(var.2, na.rm=T))
}

result.endotype.cat <- data.frame(var=var.cat, p.f, p.chi, freq.1, freq.2)

#######################################################
# severity
#######################################################
# read clinical data
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 5)
clinical[clinical=='.mild'] <- NA
clinical[clinical=='.'] <- NA

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# reorder to match protein data
clinical <- clinical[match(ptid$ptid, clinical$ptid),]

# assign cluster 
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)

# continuous variable
var.cont <- c('pf0', 'pf1', 'pf2', 'pf3', 'pf4','pf5','pf6','pf7')

p.wil <- c()
p.t <- c()

mean.1 <- c(); mean.2 <- c()
median.1 <- c(); median.2 <- c()

for(v in var.cont){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]) 
  var.2 <- as.numeric(clinical[clinical$cluster==2,v]) 
  
  w <- wilcox.test(var.1, var.2, alternative = "two.sided", na.action='na.omit')
  t <- t.test(var.1, var.2, alternative = "two.sided", na.action='na.omit')
  
  p.wil <- c(p.wil, w$p.value)
  p.t <- c(p.t, t$p.value)
  
  mean.1 <- c(mean.1, mean(var.1, na.rm=T))
  mean.2 <- c(mean.2, mean(var.2, na.rm=T))

  median.1 <- c(median.1, median(var.1, na.rm=T))
  median.2 <- c(median.2, median(var.2, na.rm=T))
}

result.endotype.cont <- data.frame(var=var.cont, p.wil, p.t, mean.1, mean.2, median.1, median.2)
#write.csv(result.endotype.cont, 'endotype.pf.csv', row.names = F)

# median
#m <- result.endotype.cont[,6:7]
#par(mfrow=c(1,1))
#matplot(m, type = c("b"),pch=1,col = 1:2)

# mean
# m <- result.endotype.cont[,4:5]
# par(mfrow=c(1,1))
# matplot(m, type = c("b"), pch=1, col = 1:2, ylab = '', xaxt = "n", lwd = 2, ylim = c(100, max(m)))
# axis(1, at=1:8, labels=var.cont)
# legend("topleft", legend = c("Endotype 1", "Endotype 2"), col=1:2, pch=1)

#######################################################
# outcome by pf
#######################################################
# read clinical data
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 13)
clinical[clinical=='.'] <- NA

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# reorder to match protein data
clinical <- clinical[match(ptid$ptid, clinical$ptid),]

# assign cluster 
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)

# continuous variable
var.cont <- c('days2dth', 'apache', "vfd","icufd","cardio28","cns28","coag28","renal28","hepatic28", "orgfree28")

p.wil <- c(); p.t <- c()
mean.1 <- c(); mean.2 <- c()
median.1 <- c(); median.2 <- c()

for(v in var.cont){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]);var.1
  var.2 <- as.numeric(clinical[clinical$cluster==2,v]);var.2
  
  w <- wilcox.test(var.1, var.2, alternative = "two.sided", na.action='na.omit'); w
  t <- t.test(var.1, var.2, alternative = "two.sided", na.action='na.omit'); t

  p.wil <- c(p.wil, w$p.value)
  p.t <- c(p.t, t$p.value)
  
  mean.1 <- c(mean.1, mean(var.1, na.rm=T))
  mean.2 <- c(mean.2, mean(var.2, na.rm=T))
  
  median.1 <- c(median.1, median(var.1, na.rm=T))
  median.2 <- c(median.2, median(var.2, na.rm=T))
}

result.endotype.cont <- data.frame(var=var.cont, p.wil, p.t, mean.1, mean.2, median.1, median.2)


# categorical variable
var.cat <- c('shock','death60')
clinical$shock <- ifelse(clinical$shock=='Yes', 1, 0)
table(clinical$shock)

p.f <- c()
p.chi <- c()
freq.1 <- c(); freq.2 <- c()

for(v in var.cat){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v])
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  tab <- clinical[clinical$cluster!=3,c(v,'cluster')]
  t <- table(tab[,v], tab[,'cluster'])
  
  f <- fisher.test(t)
  chi <- chisq.test(t)
  
  p.f <- c(p.f, f$p.value)
  p.chi <- c(p.chi, chi$p.value)
  
  freq.1 <- c(freq.1, mean(var.1, na.rm=T))
  freq.2 <- c(freq.2, mean(var.2, na.rm=T))
}

result.endotype.cat <- data.frame(var=var.cat, p.f, p.chi, freq.1, freq.2)
