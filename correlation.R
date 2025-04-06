library(readxl)
library(stringr)
library(ggplot2)
library(missForest)
library(vsn)
library(sva)

setwd("C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/Dataset 1-from Marco for ARDS")

data.DE <- readRDS('data.DE.rds'); dim(data.DE)
ptid <- read_excel("ptid.xlsx")

# keep ARDS
data.ARDS <- data.DE[,-c(1, 3, 4, 7:55)]; dim(data.ARDS)
max(as.numeric(apply(data.ARDS[,-c(1,2,3)],2,max)))
min(as.numeric(apply(data.ARDS[,-c(1,2,3)],2,min)))

#######################################################
# baseline
#######################################################
# read clinical data
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 1)
clinical[clinical=='.'] <- NA

# reorder to match protein data
clinical <- clinical[match(ptid$ptid, clinical$ptid),]

# continuous variable
var.cont <- c('age', 'height', 'weight', 'pbw', 'bmi', 'quads', 'pao2screen', 'fio2screen',
          "hrate","sysbp","diabp","cvp","map",'temp')
#var.cont <- c('quads', 'pao2screen', 'fio2screen')
result.cor.cont <- data.frame()

for(v in var.cont){
  r.pearson <- c()
  r.spearman <- c()
  p.pearson <- c()
  p.spearman <- c()

  for(i in 1:dim(data.ARDS)[1]){
    
    var.p <- as.numeric(data.ARDS[i,-c(1:3)])
    var.c <- as.numeric(clinical[[`v`]])
    
    # pearson correlation
    t.pearson <- cor.test(var.p, var.c, method = 'pearson')
    p.pearson <- c(p.pearson, t.pearson$p.value)
    r.pearson <- c(r.pearson, t.pearson$estimate)
    
    # spearman correlation
    t.spearman <- cor.test(var.p, var.c, method = 'spearman')
    p.spearman <- c(p.spearman, t.spearman$p.value)
    r.spearman <- c(r.spearman, t.spearman$estimate)
  }
  
  # FDR control
  p.values.adjust.pearson <- p.adjust(p.pearson, method = "BH")
  p.values.adjust.spearman <- p.adjust(p.spearman, method = "BH")
  
  # combine
  d <- cbind(type=rep('baseline', dim(data.ARDS)[1]), var=rep(v, dim(data.ARDS)[1]), data.ARDS[,1:3], 
             r.pearson, p.values.adjust.pearson, p.pearson,
             r.spearman, p.values.adjust.spearman, p.spearman)
  result.cor.cont <- rbind(result.cor.cont, d)
}
write.csv(result.cor.cont, 'result.cor.cont.baseline.csv', row.names = F)


# categorical variable
var.cat <- c('gender','ethnic','Race','vaso', 'death90', 'intfeed')
result.cor.cat <- data.frame()

for(v in var.cat){
  p <- c()
  fc <- c()
  for(i in 1:dim(data.ARDS)[1]){
    var.p <- as.numeric(data.ARDS[i,-c(1:3)])
    var.c <- as.numeric(clinical[[`v`]])
    l <- split(var.p, var.c)
    w <- wilcox.test(l[[1]], l[[2]])
    p <- c(p, w$p.value)
    fc <- c(fc, (mean(l[[1]]) - mean(l[[2]])))
  }
  
  # FDR control
  p.values.adjust <- p.adjust(p, method = "BH")
  
  d <- cbind(type=rep('baseline', dim(data.ARDS)[1]),var=rep(v, dim(data.ARDS)[1]), 
             data.ARDS[,1:3], fc, p.values.adjust.vil=p.values.adjust, p)
  result.cor.cat <- rbind(result.cor.cat, d)
}
write.csv(result.cor.cat, 'result.cor.cat.baseline.csv', row.names = F)

#######################################################
# Etiology
#######################################################
# read clinical data
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 2)
clinical[clinical=='.'] <- NA

# reorder to match protein data
clinical <- clinical[match(ptid$ptid, clinical$ptid),]

# categorical variable
var.cat <- c("Trauma", "Sepsis", "Transf", "Aspir", "Pneumo")
result.cor.cat <- data.frame()

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

for(v in var.cat){
  p <- c()
  fc <- c()
  for(i in 1:dim(data.ARDS)[1]){
    var.p <- as.numeric(data.ARDS[i,-c(1:3)])
    var.c <- as.numeric(clinical[[`v`]])
    l <- split(var.p, var.c)
    w <- wilcox.test(l[[1]], l[[2]])
    p <- c(p, w$p.value)
    fc <- c(fc, (mean(l[[1]]) - mean(l[[2]])))
  }
  
  # FDR control
  p.values.adjust <- p.adjust(p, method = "BH")
  
  d <- cbind(type=rep('Etiology', dim(data.ARDS)[1]),var=rep(v, dim(data.ARDS)[1]), 
             data.ARDS[,1:3], fc, p.values.adjust.vil=p.values.adjust, p)
  result.cor.cat <- rbind(result.cor.cat, d)
}
write.csv(result.cor.cat, 'result.cor.cat.Etiology.csv', row.names = F)

#######################################################
# lab test
#######################################################
# read clinical data
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 3)
clinical[clinical=='.'] <- NA

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# reorder to match protein data
clinical <- clinical[match(ptid$ptid, clinical$ptid),]

# continuous variable
var.cont <- c('hgb_0', 'sodium_0', 'potas_0' ,'gluc_0', 'bicarb_0', 'phos_0', 'mg_0', 
              'protein_0', 'album_0', 'insulinrt_1', 'insulinsq_1')
result.cor.cont <- data.frame()

for(v in var.cont){
  r.pearson <- c()
  r.spearman <- c()
  p.pearson <- c()
  p.spearman <- c()
  
  for(i in 1:dim(data.ARDS)[1]){
    
    var.p <- as.numeric(data.ARDS[i,-c(1:3)])
    var.c <- as.numeric(clinical[[`v`]])
    
    # pearson correlation
    t.pearson <- cor.test(var.p, var.c, method = 'pearson')
    p.pearson <- c(p.pearson, t.pearson$p.value)
    r.pearson <- c(r.pearson, t.pearson$estimate)
    
    # spearman correlation
    t.spearman <- cor.test(var.p, var.c, method = 'spearman')
    p.spearman <- c(p.spearman, t.spearman$p.value)
    r.spearman <- c(r.spearman, t.spearman$estimate)
  }
  
  # FDR control
  p.values.adjust.pearson <- p.adjust(p.pearson, method = "BH")
  p.values.adjust.spearman <- p.adjust(p.spearman, method = "BH")
  
  # combine
  d <- cbind(type=rep('lab.test', dim(data.ARDS)[1]), var=rep(v, dim(data.ARDS)[1]), data.ARDS[,1:3], 
             r.pearson, p.values.adjust.pearson, p.pearson,
             r.spearman, p.values.adjust.spearman, p.spearman)
  result.cor.cont <- rbind(result.cor.cont, d)
}
write.csv(result.cor.cont, 'result.cor.cont.labtest.csv', row.names = F)

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

# continuous variable
var.cont <- c('pf0', 'pf1', 'pf2','pf3', 'pf4')
result.cor.cont <- data.frame()

for(v in var.cont){
  r.pearson <- c()
  r.spearman <- c()
  p.pearson <- c()
  p.spearman <- c()
  
  for(i in 1:dim(data.ARDS)[1]){
    
    var.p <- as.numeric(data.ARDS[i,-c(1:3)])
    var.c <- as.numeric(clinical[[`v`]])
    
    # pearson correlation
    t.pearson <- cor.test(var.p, var.c, method = 'pearson')
    p.pearson <- c(p.pearson, t.pearson$p.value)
    r.pearson <- c(r.pearson, t.pearson$estimate)
    
    # spearman correlation
    t.spearman <- cor.test(var.p, var.c, method = 'spearman')
    p.spearman <- c(p.spearman, t.spearman$p.value)
    r.spearman <- c(r.spearman, t.spearman$estimate)
  }
  
  # FDR control
  p.values.adjust.pearson <- p.adjust(p.pearson, method = "BH")
  p.values.adjust.spearman <- p.adjust(p.spearman, method = "BH")
  
  # combine
  d <- cbind(type=rep('severity', dim(data.ARDS)[1]), var=rep(v, dim(data.ARDS)[1]), data.ARDS[,1:3], 
             r.pearson, p.values.adjust.pearson, p.pearson,
             r.spearman, p.values.adjust.spearman, p.spearman)
  result.cor.cont <- rbind(result.cor.cont, d)
}
write.csv(result.cor.cont, 'result.cor.cont.severity.csv', row.names = F)


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

# continuous variable
var.cont <- c('days2dth', 'apache', "vfd","icufd","cardio28","cns28","coag28","renal28","hepatic28", "orgfree28")
result.cor.cont <- data.frame()

for(v in var.cont){
  r.pearson <- c()
  r.spearman <- c()
  p.pearson <- c()
  p.spearman <- c()
  
  for(i in 1:dim(data.ARDS)[1]){
    
    var.p <- as.numeric(data.ARDS[i,-c(1:3)])
    var.c <- as.numeric(clinical[[`v`]])
    
    # pearson correlation
    t.pearson <- cor.test(var.p, var.c, method = 'pearson')
    p.pearson <- c(p.pearson, t.pearson$p.value)
    r.pearson <- c(r.pearson, t.pearson$estimate)
    
    # spearman correlation
    t.spearman <- cor.test(var.p, var.c, method = 'spearman')
    p.spearman <- c(p.spearman, t.spearman$p.value)
    r.spearman <- c(r.spearman, t.spearman$estimate)
  }
  
  # FDR control
  p.values.adjust.pearson <- p.adjust(p.pearson, method = "BH")
  p.values.adjust.spearman <- p.adjust(p.spearman, method = "BH")
  
  # combine
  d <- cbind(type=rep('outcome', dim(data.ARDS)[1]), var=rep(v, dim(data.ARDS)[1]), data.ARDS[,1:3], 
             r.pearson, p.values.adjust.pearson, p.pearson,
             r.spearman, p.values.adjust.spearman, p.spearman)
  result.cor.cont <- rbind(result.cor.cont, d)
}
write.csv(result.cor.cont, 'result.cor.cont.outcomepf.csv', row.names = F)


# categorical variable
var.cat <- c('shock','death60')
clinical$shock <- ifelse(clinical$shock=='Yes', 1, 0)
result.cor.cat <- data.frame()

for(v in var.cat){
  p <- c()
  fc <- c()
  for(i in 1:dim(data.ARDS)[1]){
    var.p <- as.numeric(data.ARDS[i,-c(1:3)])
    var.c <- as.numeric(clinical[[`v`]])
    l <- split(var.p, var.c)
    w <- wilcox.test(l[[1]], l[[2]])
    p <- c(p, w$p.value)
    fc <- c(fc, (mean(l[[1]]) - mean(l[[2]])))
  }
  
  # FDR control
  p.values.adjust <- p.adjust(p, method = "BH")
  
  d <- cbind(type=rep('outcome', dim(data.ARDS)[1]),var=rep(v, dim(data.ARDS)[1]), 
             data.ARDS[,1:3], fc, p.values.adjust.vil=p.values.adjust, p)
  result.cor.cat <- rbind(result.cor.cat, d)
}
write.csv(result.cor.cat, 'result.cor.cat.outcomepf.csv', row.names = F)


#######################################################
# vent
#######################################################
# read clinical data
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 11)
clinical[clinical=='.'] <- NA

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# reorder to match protein data
clinical <- clinical[match(ptid$ptid, clinical$ptid),]

# continuous variable
var.cont <- c(colnames(clinical)[4:16], "fio2abg_1", 'pao2abg_1', 'paco2abg_1', 'phabg_1', 'spo2abg_1')
result.cor.cont <- data.frame()

for(v in var.cont){
  r.pearson <- c()
  r.spearman <- c()
  p.pearson <- c()
  p.spearman <- c()
  
  for(i in 1:dim(data.ARDS)[1]){
    
    var.p <- as.numeric(data.ARDS[i,-c(1:3)])
    var.c <- as.numeric(clinical[[`v`]])
    
    # pearson correlation
    t.pearson <- cor.test(var.p, var.c, method = 'pearson')
    p.pearson <- c(p.pearson, t.pearson$p.value)
    r.pearson <- c(r.pearson, t.pearson$estimate)
    
    # spearman correlation
    t.spearman <- cor.test(var.p, var.c, method = 'spearman')
    p.spearman <- c(p.spearman, t.spearman$p.value)
    r.spearman <- c(r.spearman, t.spearman$estimate)
  }
  
  # FDR control
  p.values.adjust.pearson <- p.adjust(p.pearson, method = "BH")
  p.values.adjust.spearman <- p.adjust(p.spearman, method = "BH")
  
  # combine
  d <- cbind(type=rep('vent', dim(data.ARDS)[1]), var=rep(v, dim(data.ARDS)[1]), data.ARDS[,1:3], 
             r.pearson, p.values.adjust.pearson, p.pearson,
             r.spearman, p.values.adjust.spearman, p.spearman)
  result.cor.cont <- rbind(result.cor.cont, d)
}
write.csv(result.cor.cont, 'result.cor.cont.vent.csv', row.names = F)

