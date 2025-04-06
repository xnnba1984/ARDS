library(readxl)
library(stringr)
library(ggplot2)
library(missForest)
library(vsn)
library(sva)

setwd("C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/Dataset 1-from Marco for ARDS")

data.complete <- readRDS('data.complete.80.rds'); dim(data.complete)
index.keep <- readRDS('index.keep.rds')
gene.keep <- readRDS('gene.keep.rds')
accession.keep <- readRDS('accession.keep.rds')

# wilcoxon rank test
set.seed(2023)
p.values <- apply(data.complete, 1, function(x){
  wilcox.test(x[1:49], x[50:106])$p.value
} )

# FDR control
p.values.adjust <- p.adjust(p.values , method = "BH")
index.sign <- which(p.values.adjust<0.05); length(index.sign)

# log2 fold change
diff.log2 <- apply(data.complete, 1, function(x){
  mean(x[50:106]) - mean(x[1:49])
})

# save test result for visualization
DE.summary <- data.frame(gene=gene.keep, adjusted_pvalue=p.values.adjust, log2_fold_change=diff.log2)
saveRDS(DE.summary, 'DE.summary.rds')

index.up <- which(diff.log2 > 1); length(index.up)
index.down <- which(diff.log2 < -1); length(index.down)

# significant and large change
index.sign.up <- intersect(index.sign, index.up); length(index.sign.up)
index.sign.down <- intersect(index.sign, index.down); length(index.sign.down)

# combine result
result.up <- data.frame(index_original=index.keep[index.sign.up], index_complete=index.sign.up,
                        adjusted_pvalue=p.values.adjust[index.sign.up], log2_fold=diff.log2[index.sign.up])
result.down <- data.frame(index_original=index.keep[index.sign.down], index_complete=index.sign.down,
                          adjusted_pvalue=p.values.adjust[index.sign.down], log2_fold=diff.log2[index.sign.down])

# save DE matrix
data.up <- cbind(result.up, data.complete[result.up$index_complete,])
data.down <- cbind(result.down, data.complete[result.down$index_complete,])
data.DE <- rbind(data.up, data.down)
data.DE$gene.keep <- gene.keep[data.DE$index_complete]
data.DE$accession.keep <- accession.keep[data.DE$index_complete]

#data.up <- data.frame(data.complete[index.sign.up,]); dim(data.up)
#data.down <- data.frame(data.complete[index.sign.down,]); dim(data.down)
#data.DE <- cbind(gene.keep, accession.keep, p.values.adjust, diff.log2, data.frame(data.complete))
#data.DE <- data.DE[c(index.sign.up, index.sign.down),]; dim(data.DE)
#data.DE <- data.frame(data.DE)
#summary(data.DE)


# remove duplicated proteins
data.dup <- data.DE[duplicated(data.DE$gene.keep)|duplicated(data.DE$gene.keep, fromLast=T), 
                    c('accession.keep', 'gene.keep', 'log2_fold', 'adjusted_pvalue')]; dim(data.dup)
gene.dup <- unique(data.dup$gene.keep)
index.dup.drop <- c()

for(g in gene.dup){
  d <- data.dup[data.dup$gene.keep==g,]; d
  i <- which.max(abs(d$log2_fold))
  index.dup.drop <- c(index.dup.drop, setdiff(rownames(d), rownames(d[i,])))
}
data.DE <- data.DE[!rownames(data.DE)%in%index.dup.drop,]; dim(data.DE)

# change column order
data.DE <- cbind(data.DE[,c('accession.keep','gene.keep')], 
                         subset(data.DE, select=-c(accession.keep,gene.keep)))

saveRDS(data.DE, 'data.DE.rds')

# save up and down regulated proteins for visualization
data.up <- data.DE[data.DE$log2_fold>0,]
data.down <- data.DE[data.DE$log2_fold<0,]

saveRDS(data.up, 'data.up.rds')
saveRDS(data.down, 'data.down.rds')

##############################################################
# match raw data
##############################################################
# read raw data
data <- read_excel("011420-and-080919-and-052521-data-combined-Ji-BAL-EV-quant_proteins - Copy.xlsx", )
data$index_original <- as.numeric(rownames(data))

# merge
result.down.raw <- merge(data, result.down, by = 'index_original')
result.up.raw <- merge(data, result.up, by = 'index_original')

# remove NA genes
result.down.raw <- result.down.raw[complete.cases(result.down.raw$`Gene Name`),]
result.up.raw <- result.up.raw[complete.cases(result.up.raw$`Gene Name`),]

# check accession duplicate
result.down.raw[duplicated(result.down.raw$Accession), c('Accession','Gene Name','adjusted_pvalue','log2_fold')]
result.up.raw[duplicated(result.up.raw$Accession), c('Accession','Gene Name','adjusted_pvalue','log2_fold')]

# check gene duplicate
result.down.raw[duplicated(result.down.raw$`Gene Name`)|duplicated(result.down.raw$`Gene Name`, fromLast=T), 
                c('Accession','Gene Name','adjusted_pvalue','log2_fold')]
result.up.raw[duplicated(result.up.raw$`Gene Name`)|duplicated(result.up.raw$`Gene Name`, fromLast=T), 
                c('Accession','Gene Name','adjusted_pvalue','log2_fold')]

# change column order
result.down.raw <- cbind(result.down.raw[,c('adjusted_pvalue','log2_fold')], 
                         subset(result.down.raw, select=-c(adjusted_pvalue,log2_fold)))
result.up.raw <- cbind(result.up.raw[,c('adjusted_pvalue','log2_fold')], 
                         subset(result.up.raw, select=-c(adjusted_pvalue,log2_fold)))

# save result
write.csv(result.up.raw, 'DP_up.csv', row.names = F, na = '')
write.csv(result.down.raw, 'DP_down.csv', row.names = F, na = '')


##############################################################################################################
# death DEP
##############################################################################################################
data.complete <- readRDS('data.complete.80.rds'); dim(data.complete)

# only patients
data.complete <- data.complete[,50:106]; dim(data.complete)

index.keep <- readRDS('index.keep.rds')
gene.keep <- readRDS('gene.keep.rds')
accession.keep <- readRDS('accession.keep.rds')

# read ptid
ptid <- read_excel("ptid.xlsx")

# read clinical data
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 1)
clinical[clinical=='.'] <- NA

# reorder to match protein data
clinical <- clinical[match(ptid$ptid, clinical$ptid),]
table(clinical$death90)

# death vs suvival index
index.death <- which(clinical$death90==0)
index.survial <- which(clinical$death90==1)

# wilcoxon rank test
set.seed(2023)
p.values <- apply(data.complete, 1, function(x){
  wilcox.test(x[index.survial], x[index.death])$p.value
} )
sort(p.values, decreasing = F)

# FDR control
#p.values.adjust <- p.adjust(p.values , method = "BH")
index.sign <- which(p.values<0.05); length(index.sign)

# log2 fold change
diff.log2 <- apply(data.complete, 1, function(x){
  mean(x[index.death]) - mean(x[index.survial])
})

index.up <- which(diff.log2 > 1); length(index.up)
index.down <- which(diff.log2 < -1); length(index.down)

# significant and large change
index.sign.up <- intersect(index.sign, index.up); length(index.sign.up)
index.sign.down <- intersect(index.sign, index.down); length(index.sign.down)

# combine result
result.up <- data.frame(index_original=index.keep[index.sign.up], index_complete=index.sign.up,
                        adjusted_pvalue=p.values[index.sign.up], log2_fold=diff.log2[index.sign.up])
result.down <- data.frame(index_original=index.keep[index.sign.down], index_complete=index.sign.down,
                          adjusted_pvalue=p.values[index.sign.down], log2_fold=diff.log2[index.sign.down])

# save DE matrix
data.up <- cbind(result.up, data.complete[result.up$index_complete,])
data.down <- cbind(result.down, data.complete[result.down$index_complete,])
data.DE <- rbind(data.up, data.down)
data.DE$gene.keep <- gene.keep[data.DE$index_complete]
data.DE$accession.keep <- accession.keep[data.DE$index_complete]

# remove duplicated proteins
data.dup <- data.DE[duplicated(data.DE$gene.keep)|duplicated(data.DE$gene.keep, fromLast=T), 
                    c('accession.keep', 'gene.keep', 'log2_fold', 'adjusted_pvalue')]; dim(data.dup)
gene.dup <- unique(data.dup$gene.keep)
index.dup.drop <- c()

for(g in gene.dup){
  d <- data.dup[data.dup$gene.keep==g,]; d
  i <- which.max(abs(d$log2_fold))
  index.dup.drop <- c(index.dup.drop, setdiff(rownames(d), rownames(d[i,])))
}
data.DE <- data.DE[!rownames(data.DE)%in%index.dup.drop,]; dim(data.DE)

# change column order
data.DE <- cbind(data.DE[,c('accession.keep','gene.keep')], 
                 subset(data.DE, select=-c(accession.keep,gene.keep)))

saveRDS(data.DE, 'data.DE.death.rds')
