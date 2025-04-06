library(readxl)
library(stringr)
library(ggplot2)
library(missForest)
library(vsn)
library(sva)

setwd("C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/Dataset 1-from Marco for ARDS")

files <- list.files("C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/Dataset 1-from Marco for ARDS")

#######################################################
# combine continuous variables
#######################################################
files.cont <- files[which(grepl("result.cor.cont", files))]

# Create an empty list to store the data frames
result.cont.list <- list()

# Loop through the file names and read each file
for (file in files.cont) {
  data <- read.csv(file)
  result.cont.list[[file]] <- data
}

# combine
result.cont <- do.call(rbind, result.cont.list)
row.names(result.cont) <- NULL
table(result.cont$type)

# nominal significance weak
result.nominal <- result.cont[(abs(result.cont$r.pearson)>0.4 & result.cont$p.pearson<0.05) | 
                                (abs(result.cont$r.spearman)>0.4 & result.cont$p.spearman<0.05),]
sort(table(result.nominal$var))
sort(table(result.nominal$type))
sort(table(result.nominal$gene.keep))
length(unique(result.nominal$gene.keep))

# remove redundance; change column name
#result.nominal <- result.nominal[,-c(4,5)]
#colnames(result.nominal)[6] <- 'p.value.raw.pearson'
#colnames(result.nominal)[9] <- 'p.value.raw.spearman'

# choose important variable
var.imp <- c('apache','vfd','icufd','cardio28','cns28','coag28','renal28','hepatic28','orgfree28','days2dth',
             'pf0','pf1','pf2','pf3','pf4','age','height','weight','bmi')

# var.imp <- c('orgfree28','pf1','vfd','days2dth', 'icufd', 'cardio28',
#              'renal28','coag28','hepatic28', 'cns28', 
#              'pf0','pf2','apache', 'age', 'height', 'weight', 'bmi',
#              'pf3','pf4')
result.imp <- result.nominal[result.nominal$var%in%var.imp,]

# merge pf
result.imp[result.imp$var%in%c('pf0','pf1','pf2', 'pf3','pf4'),'var'] <- 'pf'

# remove duplicate protein in pf
result.imp.pf <- result.imp[result.imp$var=='pf',]
table(result.imp.pf$gene.keep); length(table(result.imp.pf$gene.keep))
index.dup <- which(duplicated(result.imp.pf$gene.keep) | duplicated(result.imp.pf$gene.keep, fromLast = TRUE))
result.imp.pf.dup <- result.imp.pf[index.dup,]  # duplicated proteins

# keep proteins with larger ave between pearson and spearman
index.keep <- c()
for(g in unique(result.imp.pf.dup$gene.keep)){
  #g <- unique(result.imp.pf.dup$gene.keep)[1]
  d <- result.imp.pf.dup[result.imp.pf.dup$gene.keep==g,]
  r <- apply(d, 1, function(x){
    mean(c(abs(as.numeric(x[6])), abs(as.numeric(x[9]))))
  })
  index.keep <- c(index.keep, as.numeric(names(which.max(r))))
}
index.remove <- setdiff(rownames(result.imp.pf.dup), index.keep)
result.imp <- result.imp[!rownames(result.imp)%in%index.remove,]

sort(table(result.imp$var))
sort(table(result.imp$type))
sort(table(result.imp$gene.keep))
length(unique(result.imp$gene.keep))

# save important variables
write.csv(result.imp, 'cor.cont.imp.csv', row.names = F)

# other variables
result.noimp <- result.nominal[!result.nominal$var%in%var.imp,]; dim(result.noimp); dim(result.imp); dim(result.nominal)
write.csv(result.noimp, 'cor.cont.noimp.csv', row.names = F)

var.cont <- data.frame(sort(table(result.imp$var), decreasing = T))
var.cont$Type <- rep('continuous', dim(var.cont)[1])
#length(unique(var.cont$gene.keep))
#length(unique(result.cont$var))
write.csv(var.cont, 'var.cont.csv', row.names = F)

#######################################################
# combine categorical variables
#######################################################
files <- list.files("C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/Dataset 1-from Marco for ARDS")
files.cat <- files[which(grepl("result.cor.cat", files))]

# Create an empty list to store the data frames
result.cat.list <- list()

# Loop through the file names and read each file
for (file in files.cat) {
  data <- read.csv(file)
  result.cat.list[[file]] <- data
}

# combine
result.cat <- do.call(rbind, result.cat.list)
table(result.cat$type)

# remove redundance; change column name
result.cat <- result.cat[,-c(4,5)]
colnames(result.cat)[6] <- 'p.value.raw.wil'

# significance
result.cat.sig <- result.cat[result.cat$p.values.adjust.vil<0.05 | 
                               (abs(result.cat$fc)>1 & result.cat$p.value.raw.wil<0.05),]
sort(table(result.cat.sig$var))
sort(table(result.cat.sig$gene.keep))
length(unique(result.cat.sig$gene.keep))

# change logFC to FC
result.cat.sig$fc <- ifelse(result.cat.sig$fc>0, 2^result.cat.sig$fc, -1/2^result.cat.sig$fc)

# keep important variables
var.imp <- c('Sepsis', 'shock', 'death60', 'death90', 'Trauma', 'Sepsite', 'Transf', 'Aspir', 'Pneumo', 'gender', 'Race')
result.cat.imp <- result.cat.sig[result.cat.sig$var%in%var.imp,]; dim(result.cat.imp)

sort(table(result.cat.imp$var))
sort(table(result.cat.imp$gene.keep))
length(unique(result.cat.imp$gene.keep))
sort(table(result.cat.imp$type))

# save important variables
write.csv(result.cat.imp, 'cor.cat.imp.csv',row.names = F)

# other variables
result.cat.noimp <- result.cat.sig[!result.cat.sig$var%in%var.imp,]; dim(result.cat.noimp)

sort(table(result.cat.noimp$var))
sort(table(result.cat.noimp$gene.keep))
length(unique(result.cat.noimp$gene.keep))

# save others
write.csv(result.cat.noimp, 'cor.cat.noimp.csv',row.names = F)

# check genes and variables
sort(table(c(result.imp$gene.keep, result.cat.sig$gene.keep)))
length(unique(c(result.imp$gene.keep, result.cat.sig$gene.keep)))
length(unique(c(result.imp$var, result.cat.sig$var)))

var.cat <- data.frame(sort(table(result.cat.imp$var), decreasing = T))
var.cat$Type <- rep('categorical', dim(var.cat)[1])
write.csv(var.cat, 'var.cat.csv', row.names = F)

# check unique proteins
gene.all <- unique(c(result.imp$gene.keep, result.cat.imp$gene.keep))
gene.all <- c(result.imp$gene.keep, result.cat.imp$gene.keep)
gene.all <- c(result.imp$gene.keep[result.imp$type!='baseline'], result.cat.imp$gene.keep[result.cat.imp$type!='baseline'])
gene.all <- unique(c(result.imp$gene.keep[result.imp$type!='baseline'], result.cat.imp$gene.keep[result.cat.imp$type!='baseline']))
gene.all <- unique(result.imp$gene.keep)
gene.all <- unique(result.cat.imp$gene.keep)


# #######################################################
# # summary
# #######################################################
# var.summary <- rbind(var.cont, var.cat)
# colnames(var.summary)[1:2] <- c('Variable', 'Correlated Protein')
# write.csv(var.summary, 'var.summary.csv', row.names = F)
