library(readxl)
library(stringr)
library(ggplot2)
library(missForest)
library(vsn)
library(sva)
library(NbClust) 
library(ComplexHeatmap)

setwd("C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/Dataset 1-from Marco for ARDS")

cluster <- readRDS('kmeans.rds')
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
#write.csv(ptid, 'ptid.cluster.csv', row.names = F)
ptid <- ptid[,c(2,4)]

# assign cluster to clinical data
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 1)
clinical[clinical=='.'] <- NA
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)
#colnames(clinical)

#######################################################
# baseline
#######################################################
# continuous variables
var.cont <- c('age', 'height', 'weight', 'pbw', 'bmi', 'packyr', 'quads', 'intubdt', 'pao2screen', 'fio2screen','qualpfdt','critdt',
              "hrate","sysbp","diabp","cvp","map",'temp')
for(v in var.cont){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]) 
  var.2 <- as.numeric(clinical[clinical$cluster==2,v]) 
  w1 <- wilcox.test(var.1, var.2, alternative = "two.sided", na.action='na.omit')
  pvalues <- c(w1$p.value)
  
  if(sum(pvalues < 0.1) > 0){
    print('============================='); print(v); print('=============================')
    print(mean(var.1, na.rm=T)); print(mean(var.2, na.rm=T))
    print(w1)
  }
}

# categorical variables
var.cat <- c('gender','ethnic','Race','vaso','smoker','cursmoker', 'death90', 'intfeed')
for(v in var.cat){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v])
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  tab <- clinical[clinical$cluster!=3,c(v,'cluster')]
  t <- table(tab[,v], tab[,'cluster']);t
  #f1 <- fisher.test(t)
  f1 <- chisq.test(t)
  pvalues <- c(f1$p.value)
  
  if(sum(pvalues < 0.1) > 0){
    print('============================='); print(v); print('=============================')
    print(mean(var.1, na.rm=T)); print(mean(var.2, na.rm=T))
    print(f1)
  }
}

#######################################################
# Etiology
#######################################################
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 2)
clinical[clinical=='.'] <- NA
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)
colnames(clinical)

# continuous variables
var.cont <- c('hasddt', 'icudt')
for(v in var.cont){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]) 
  var.2 <- as.numeric(clinical[clinical$cluster==2,v]) 
  w1 <- wilcox.test(var.1, var.2, alternative = "two.sided", na.action='na.omit')
  pvalues <- c(w1$p.value)
  
  if(sum(pvalues < 0.1) > 0){
    print('============================='); print(v); print('=============================')
    print(mean(var.1, na.rm=T)); print(mean(var.2, na.rm=T))
    print(w1)
  }
}

# categorical variables
#clinical[clinical$Trauma==2,'Trauma'] <- 1
var.cat <- c("Trauma", "Sepsite", "Sepsis", "Transf", "Aspir", "Pneumo", "Otherlung", "admtype", "admitfrom")
for(v in var.cat){
  #v <- var.cat[2]
  var.1 <- as.numeric(clinical[clinical$cluster==1,v])
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  tab <- clinical[clinical$cluster!=3,c(v,'cluster')]
  t <- table(tab[,v], tab[,'cluster']); print(t)
  f1 <- fisher.test(t)
  #f1 <- chisq.test(t)
  pvalues <- c(f1$p.value)
  
  if(sum(pvalues < 0.1) > 0){
    print('============================='); print(v); print('=============================')
    print(mean(var.1, na.rm=T)); print(mean(var.2, na.rm=T))
    print(f1)
  }
}

#######################################################
# lab test
#######################################################
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 3)
clinical[clinical=='.'] <- NA
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)
#colnames(clinical)

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# continuous variables
var.cont <- colnames(clinical)[2:(length(colnames(clinical))-1)]
for(v in var.cont){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]) 
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  if(sum(is.na(var.1))<length(var.1) & sum(is.na(var.2))<length(var.2)){
    w1 <- wilcox.test(var.1, var.2, alternative = "two.sided", na.action='na.omit')
    pvalues <- c(w1$p.value)
    if(sum(pvalues < 0.05) > 0){
      print('============================='); print(v); print('=============================')
      print(mean(var.1, na.rm=T)); print(mean(var.2, na.rm=T))
      print(w1)
    }
  }
}

#######################################################
# cardio
#######################################################
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 4)
clinical[clinical=='.'] <- NA
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)
#colnames(clinical)

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# continuous variables
var.cont <- colnames(clinical)[2:(length(colnames(clinical))-1)]
for(v in var.cont){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]) 
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  if(sum(is.na(var.1))<length(var.1) & sum(is.na(var.2))<length(var.2)){
    w1 <- wilcox.test(var.1, var.2, alternative = "two.sided", na.action='na.omit')
    pvalues <- c(w1$p.value)
    if(sum(pvalues < 0.05) > 0){
      print('============================='); print(v); print('=============================')
      print(mean(var.1, na.rm=T)); print(mean(var.2, na.rm=T))
      print(w1)
    }
  }
}

#######################################################
# pf ratio
#######################################################
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 5)
clinical[clinical=='.'] <- NA
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)
#colnames(clinical)

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# continuous variables
var.cont <- colnames(clinical)[2:(length(colnames(clinical))-2)]
for(v in var.cont){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]) 
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  if(sum(is.na(var.1))<length(var.1) & sum(is.na(var.2))<length(var.2)){
    w1 <- wilcox.test(var.1, var.2, alternative = "two.sided", na.action='na.omit')
    pvalues <- c(w1$p.value)
    if(sum(pvalues < 0.05) > 0){
      print('============================='); print(v); print('=============================')
      print(mean(var.1, na.rm=T)); print(mean(var.2, na.rm=T))
      print(w1)
    }
  }
}

#######################################################
# systbp
#######################################################
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 6)
clinical[clinical=='.'] <- NA
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)
#colnames(clinical)

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# continuous variables
var.cont <- colnames(clinical)[2:(length(colnames(clinical))-1)]
for(v in var.cont){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]) 
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  if(sum(is.na(var.1))<length(var.1) & sum(is.na(var.2))<length(var.2)){
    w1 <- wilcox.test(var.1, var.2, alternative = "two.sided", na.action='na.omit')
    pvalues <- c(w1$p.value)
    if(sum(pvalues < 0.05) > 0){
      print('============================='); print(v); print('=============================')
      print(mean(var.1, na.rm=T)); print(mean(var.2, na.rm=T))
      print(w1)
    }
  }
}

#######################################################
# plate
#######################################################
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 7)
clinical[clinical=='.'] <- NA
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)
#colnames(clinical)

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# continuous variables
var.cont <- colnames(clinical)[2:(length(colnames(clinical))-1)]
for(v in var.cont){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]) 
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  if(sum(is.na(var.1))<length(var.1) & sum(is.na(var.2))<length(var.2)){
    w1 <- wilcox.test(var.1, var.2, alternative = "two.sided", na.action='na.omit')
    pvalues <- c(w1$p.value)
    if(sum(pvalues < 0.05) > 0){
      print('============================='); print(v); print('=============================')
      print(mean(var.1, na.rm=T)); print(mean(var.2, na.rm=T))
      print(w1)
    }
  }
}

#######################################################
# creat
#######################################################
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 8)
clinical[clinical=='.'] <- NA
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)
#colnames(clinical)

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# continuous variables
var.cont <- colnames(clinical)[2:(length(colnames(clinical))-1)]
for(v in var.cont){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]) 
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  if(sum(is.na(var.1))<length(var.1) & sum(is.na(var.2))<length(var.2)){
    w1 <- wilcox.test(var.1, var.2, alternative = "two.sided", na.action='na.omit')
    pvalues <- c(w1$p.value)
    if(sum(pvalues < 0.05) > 0){
      print('============================='); print(v); print('=============================')
      print(mean(var.1, na.rm=T)); print(mean(var.2, na.rm=T))
      print(w1)
    }
  }
}

#######################################################
# bili
#######################################################
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 9)
clinical[clinical=='.'] <- NA
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)
#colnames(clinical)

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# continuous variables
var.cont <- colnames(clinical)[2:(length(colnames(clinical))-1)]
for(v in var.cont){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]) 
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  if(sum(is.na(var.1))<length(var.1) & sum(is.na(var.2))<length(var.2)){
    w1 <- wilcox.test(var.1, var.2, alternative = "two.sided", na.action='na.omit')
    pvalues <- c(w1$p.value)
    if(sum(pvalues < 0.05) > 0){
      print('============================='); print(v); print('=============================')
      print(mean(var.1, na.rm=T)); print(mean(var.2, na.rm=T))
      print(w1)
    }
  }
}

#######################################################
# vaso
#######################################################
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 10)
clinical[clinical=='.'] <- NA
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)
#colnames(clinical)

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# categorical variables
var.cat <- colnames(clinical)[2:(length(colnames(clinical))-1)]
for(v in var.cat){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v])
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  tab <- clinical[clinical$cluster!=3,c(v,'cluster')]
  t <- table(tab[,v], tab[,'cluster']);t
  if(length(t)>=4){
    f1 <- fisher.test(t)
    pvalues <- c(f1$p.value)
    
    if(sum(pvalues < 0.1) > 0){
      print('============================='); print(v); print('=============================')
      print(mean(var.1, na.rm=T)); print(mean(var.2, na.rm=T))
      print(f1)
    }
  }
}

#######################################################
# vent
#######################################################
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 11)
clinical[clinical=='.'] <- NA
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)
#colnames(clinical)

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# continuous variables
var.cont <- colnames(clinical)[c(4:16, 25:(length(colnames(clinical))-1))]
for(v in var.cont){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]) 
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  if(sum(is.na(var.1))<length(var.1) & sum(is.na(var.2))<length(var.2)){
    w1 <- wilcox.test(var.1, var.2, alternative = "two.sided", na.action='na.omit')
    pvalues <- c(w1$p.value)
    if(sum(pvalues < 0.05) > 0){
      print('============================='); print(v); print('=============================')
      print(mean(var.1, na.rm=T)); print(mean(var.2, na.rm=T))
      print(w1)
    }
  }
}

# categorical variables
var.cat <- colnames(clinical)[c(2,17:24)]
for(v in var.cat){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v])
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  tab <- clinical[clinical$cluster!=3,c(v,'cluster')]
  t <- table(tab[,v], tab[,'cluster']);t
  if(length(t)>=4){
    f1 <- fisher.test(t)
    pvalues <- c(f1$p.value)
    
    if(sum(pvalues < 0.1) > 0){
      print('============================='); print(v); print('=============================')
      print(t)
      print(f1)
    }
  }
}

#######################################################
# funnel
#######################################################
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 12)
clinical[clinical=='.'] <- NA
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)
#colnames(clinical)

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# continuous variables
var.cont <- c("homedt","deathdt","othstatdt","hospdcdt","uabdt","uabdt1","uabdt2",      
              "dischargedt1", "dischargedt2", "readmitdt1","deathdate","discharge")
for(v in var.cont){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]) 
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  if(sum(is.na(var.1))<length(var.1) & sum(is.na(var.2))<length(var.2)){
    w1 <- wilcox.test(var.1, var.2, alternative = "two.sided", na.action='na.omit')
    pvalues <- c(w1$p.value)
    if(sum(pvalues < 0.05) > 0){
      print('============================='); print(v); print('=============================')
      print(mean(var.1, na.rm=T)); print(mean(var.2, na.rm=T))
      print(w1)
    }
  }
}

# categorical variables
var.cat <- 'status'
for(v in var.cat){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v])
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  tab <- clinical[clinical$cluster!=3,c(v,'cluster')]
  t <- table(tab[,v], tab[,'cluster']);t
  if(length(t)>=4){
    f1 <- fisher.test(t)
    pvalues <- c(f1$p.value)
    
    if(sum(pvalues < 0.1) > 0){
      print('============================='); print(v); print('=============================')
      print(t)
      print(f1)
    }
  }
}

#######################################################
# outcomes by pf
#######################################################
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 13)
clinical[clinical=='.'] <- NA
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)
#colnames(clinical)

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# continuous variables
var.cont <- c("flbal_0","flbal_1","flbal_2","flbal_3","flbal_4","flbal_5","flbal_6","flbal_7","flbal_8","days2dth",
              "apache","vfd","icufd","cardio28","cns28","coag28","renal28","hepatic28", "orgfree28")
for(v in var.cont){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]) 
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  if(sum(is.na(var.1))<length(var.1) & sum(is.na(var.2))<length(var.2)){
    w1 <- wilcox.test(var.1, var.2, alternative = "two.sided", na.action='na.omit')
    pvalues <- c(w1$p.value)
    if(sum(pvalues < 0.05) > 0){
      print('============================='); print(v); print('=============================')
      print(mean(var.1, na.rm=T)); print(mean(var.2, na.rm=T))
      print(w1)
    }
  }
}

# categorical variables
clinical$shock <- ifelse(clinical$shock=='Yes',1,0)
var.cat <- c('shock','death60','death90')
for(v in var.cat){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v])
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  tab <- clinical[clinical$cluster!=3,c(v,'cluster')]
  t <- table(tab[,v], tab[,'cluster']);t
  if(length(t)>=4){
    f1 <- fisher.test(t)
    pvalues <- c(f1$p.value)
    
    if(sum(pvalues < 0.1) > 0){
      print('============================='); print(v); print('=============================')
      print(t)
      print(f1)
    }
  }
}

#######################################################
# outcomes by sepsis
#######################################################
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 14)
clinical[clinical=='.'] <- NA
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)
#colnames(clinical)

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# categorical variables
var.cat <- 'Sepsis'
for(v in var.cat){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v])
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  tab <- clinical[clinical$cluster!=3,c(v,'cluster')]
  t <- table(tab[,v], tab[,'cluster']);t
  if(length(t)>=4){
    f1 <- fisher.test(t)
    pvalues <- c(f1$p.value)
    
    if(sum(pvalues < 0.1) > 0){
      print('============================='); print(v); print('=============================')
      print(t)
      print(f1)
    }
  }
}

#######################################################
# apachep_fluids
#######################################################
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 16)
clinical[clinical=='.'] <- NA
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)
#colnames(clinical)

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

col.name <- colnames(clinical)
names(col.name) <- 1:length(col.name)

# continuous variables
var.cont <- col.name[c(2:11,16:45,48:53,56:61,64:69,72:77,80:85,88:93,96:101,104:107)]
for(v in var.cont){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]) 
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  if(sum(is.na(var.1))<length(var.1) & sum(is.na(var.2))<length(var.2)){
    w1 <- wilcox.test(var.1, var.2, alternative = "two.sided", na.action='na.omit')
    pvalues <- c(w1$p.value)
    if(sum(pvalues < 0.05) > 0 & is.na(pvalues) == F){
      print('============================='); print(v); print('=============================')
      print(mean(var.1, na.rm=T)); print(mean(var.2, na.rm=T))
      print(w1)
    }
  }
}

# categorical variables
var.cat <- setdiff(col.name[-c(1,108)], var.cont)
var.cat <- c("ventl","venth","notedenpt_1","notedenpt_2","notedenpt_3", 
              "notedenpt_4", "notedenpt_5", "notedenpt_6","notedenpt_7", "notedenpt_8")
for(v in var.cat){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v])
  var.2 <- as.numeric(clinical[clinical$cluster==2,v])
  
  tab <- clinical[clinical$cluster!=3,c(v,'cluster')]
  t <- table(tab[,v], tab[,'cluster']);t
  if(length(t)>=4){
    f1 <- fisher.test(t)
    pvalues <- c(f1$p.value)
    
    if(sum(pvalues < 0.1) > 0){
      print('============================='); print(v); print('=============================')
      print(t)
      print(f1)
    }
  }
}
