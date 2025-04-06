library(readxl)
library(stringr)
library(ggplot2)
library(vsn)
library(sva)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(readr)
library(tidyr)
library(dplyr)
library(ggtree)
library(forcats)
library(utils)
library(ggsignif)
library(ggpubr)
library(ggprism)
library(rstatix)
library(scales)
library(circlize)
library(venn)
library(ggsankey)
library(ggalluvial)

setwd("C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/Dataset 1-from Marco for ARDS")

############################################################################################################
############################################################################################################
# biomarker death sepsis
############################################################################################################
############################################################################################################
library(cvms)
library(tibble)
library(caret)
library(PRROC)
library(Rtsne)

# read result
#prob.matrix <- read.csv('biomarker_death_combine_all_prob.csv') # death
prob.matrix <- read.csv('biomarker_shock_combine_all_prob.csv') # E34
prob.pred <- prob.matrix[,1]

# transform to binary prediction
y.pred <- ifelse(prob.pred>0.5, 1, 0); table(y.pred)

# read ptid
ptid <- read_excel("ptid.xlsx")

# read clinical data
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 1)
clinical[clinical=='.'] <- NA

# reorder to match protein data
clinical <- clinical[match(ptid$ptid, clinical$ptid),]
#y.truth <- clinical$death90
clinical$shock <- ifelse(clinical$shock=='Yes', 1, 0); y.truth <- clinical$shock
table(y.truth)

############################################################################################################
# AUC curve
############################################################################################################
# death
output.auroc <- roc.curve(scores.class0 = prob.pred[y.truth==1], 
                          scores.class1 = prob.pred[y.truth==0], curve = TRUE)

pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/AUC_death.pdf")
plot(output.auroc, legend=F, color=F, xlab='1-Specificity', main='Death', 
     auc.main=F, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, lwd=2, xaxs="i", yaxs="i"); lines(x = c(0,1), y = c(0,1), lty=2); text('AUC = 0.926', x=0.6, y=0.2, cex=1.2)
dev.off()

# shock
output.auroc <- roc.curve(scores.class0 = prob.pred[y.truth==1], 
                          scores.class1 = prob.pred[y.truth==0], curve = TRUE)

pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/AUC_shock.pdf")
plot(output.auroc, legend=F, color=F, xlab='1-Specificity', main='Shock', 
     auc.main=F, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, lwd=2, xaxs="i", yaxs="i"); lines(x = c(0,1), y = c(0,1), lty=2); text('AUC = 0.898', x=0.6, y=0.2, cex=1.2)
dev.off()


############################################################################################################
############################################################################################################
# endotype etiology biomarker
############################################################################################################
############################################################################################################
# read clustering result
cluster <- readRDS('kmeans.biomarker.Aspir.Pneumo.Sepsis.rds')
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

# read clinical data
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 1)
clinical[clinical=='.'] <- NA

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# reorder to match protein data
clinical <- clinical[match(ptid$ptid, clinical$ptid),]

# assign cluster 
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)

# relevant continuous variable
var.cont <- c('cvp', 'pao2screen', 'pf3')
result <- data.frame()

# create visualization data
for(v in var.cont){
  variable <- rep(v, dim(clinical)[1])
  var.1 <- as.numeric(clinical[clinical$cluster==2,v]);var.1
  var.2 <- as.numeric(clinical[clinical$cluster==1,v]);var.2
  endotype <- c(rep('E4\'', length(var.1)), rep('E3\'', length(var.2)))
  d <- data.frame(variable, endotype, value=c(var.1, var.2))
  result <- rbind(result, d)
}
result$endotype[result$endotype=='E3\''] <- 'E1\''
result$endotype[result$endotype=='E4\''] <- 'E2\''

############################################################################################################
# cvp
############################################################################################################
result.list <- list()
v <- 'cvp'
r <- result[result$variable==v,]

p <- ggplot(r, aes(x=endotype, y=value, fill=endotype)) + 
  geom_boxplot(width=.35, outlier.shape = NA) + 
  theme_classic() +
  geom_jitter(shape=16, size=1, position=position_jitter(width = 0.1, height = 0, seed = 2023), aes(alpha=0.9)) +
  geom_signif(comparisons = list(c("E1\'", "E2\'")), 
              map_signif_level=c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05), 
              test = 'wilcox.test', textsize = 10, y_position = 39) +
  ylab('mmHg') +
  theme(axis.title.x=element_blank(), legend.position = "none", text = element_text(size = 19), 
        axis.title.y = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(1,0,0,1), 'lines'), axis.text = element_text(color = "black")) +
  scale_y_continuous(limits = c(0, 45)) +
  scale_fill_manual(values=c("#7fbc41", "#bf812d")) +
  ggtitle(v); result.list <- append(result.list, list(p)); p


############################################################################################################
# pao2screen
############################################################################################################
v <- 'pao2screen'
r <- result[result$variable==v,]
p <- ggplot(r, aes(x=endotype, y=value, fill=endotype)) + 
  geom_boxplot(width=.35, outlier.shape = NA) + 
  theme_classic() +
  geom_jitter(shape=16, size=1, position=position_jitter(width = 0.1, height = 0, seed = 2023), aes(alpha=0.9)) +
  geom_signif(comparisons = list(c("E1\'", "E2\'")), 
              map_signif_level=c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05), 
              test = 'wilcox.test', textsize = 10, y_position = 305) +
  ylab('mmHg') +
  theme(axis.title.x=element_blank(), legend.position = "none", text = element_text(size = 19), 
        axis.title.y = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(1,0,0,1), 'lines'), axis.text = element_text(color = "black")) +
  scale_y_continuous(limits = c(0, 350)) +
  scale_fill_manual(values=c("#7fbc41", "#bf812d")) +
  ggtitle(v); result.list <- append(result.list, list(p)); p


############################################################################################################
# pf3
############################################################################################################
v <- 'pf3'
r <- result[result$variable==v,]
p <- ggplot(r, aes(x=endotype, y=value, fill=endotype)) + 
  geom_boxplot(width=.35, outlier.shape = NA) + 
  theme_classic() +
  geom_jitter(shape=16, size=1, position=position_jitter(width = 0.1, height = 0, seed = 2023), aes(alpha=0.9)) +
  geom_signif(comparisons = list(c("E1\'", "E2\'")), 
              map_signif_level=c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.1), 
              test = 'wilcox.test', textsize = 10, y_position = 380) +
  ylab('mmHg') +
  theme(axis.title.x=element_blank(), legend.position = "none", text = element_text(size = 19), 
        axis.title.y = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(1,0,0,1), 'lines'), axis.text = element_text(color = "black")) +
  scale_y_continuous(limits = c(0, 430)) +
  scale_fill_manual(values=c("#7fbc41", "#bf812d")) +
  ggtitle(v); result.list <- append(result.list, list(p)); p


############################################################################################################
# barplot categorical
############################################################################################################
variable <- c(rep('Pneumo', 4), rep('Sepsis', 4))
outcome <- c(c('Pneumo','Pneumo','No Pneumo','No Pneumo'), c('Sepsis','Sepsis','No Sepsis','No Sepsis'))
endotype <- rep(c('E1\'','E2\''), 4)
value <- c(91.67, 42.42, 8.33, 57.58, 25.00, 69.70, 75.00, 30.3)
result.cat <- data.frame(variable,outcome,endotype,value);result.cat

# Pneumo
v <- 'Pneumo'
r <- result.cat[result.cat$variable==v,]
r$outcome <- factor(r$outcome, levels = c('Pneumo','No Pneumo'))
p <- ggplot(r, aes(y=value, x=endotype)) + 
  geom_bar(aes(fill=outcome), position="fill", stat="identity", width=.35) +
  scale_y_continuous(limits = c(0, 1.25), labels=scales::percent_format(suffix = ""), 
                     oob = rescale_none, breaks = seq(0, 1, by=0.25), expand = c(0, 0)) +
  theme_classic() +
  geom_signif(aes(xmin = 'E1\'', xmax = 'E2\'', annotations = '***', y_position = 1.05),
              textsize = 10, tip_length = .0003) + 
  scale_fill_manual(values=c("#cc0000", "#56B4E9")) + 
  ggtitle(paste(v, '(%)')) + 
  theme(axis.title.x=element_blank(), text = element_text(size = 19), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(1,0,0,1), 'lines'), axis.text = element_text(color = "black"),
        legend.title=element_blank(),legend.position="right",
        legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=10),
        legend.margin = margin(t = 0, r = 10, b = 0, l = -5))
result.list <- append(result.list, list(p)); p

# Sepsis
v <- 'Sepsis'
r <- result.cat[result.cat$variable==v,]
r$outcome <- factor(r$outcome, levels = c('Sepsis','No Sepsis'))
p <- ggplot(r, aes(y=value, x=endotype)) + 
  geom_bar(aes(fill=outcome), position="fill", stat="identity", width=.35) +
  scale_y_continuous(limits = c(0, 1.25), labels=scales::percent_format(suffix = ""), 
                     oob = rescale_none, breaks = seq(0, 1, by=0.25), expand = c(0, 0)) +
  theme_classic() +
  geom_signif(aes(xmin = 'E1\'', xmax = 'E2\'', annotations = '**', y_position = 1.05),
              textsize = 10, tip_length = .0003) + 
  scale_fill_manual(values=c("#cc0000", "#56B4E9")) + 
  ggtitle(paste(v, '(%)')) + 
  theme(axis.title.x=element_blank(), text = element_text(size = 19), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(1,0,0,1), 'lines'), axis.text = element_text(color = "black"),
        legend.title=element_blank(),legend.position="right",
        legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=10),
        legend.margin = margin(t = 0, r = 10, b = 0, l = -5))
result.list <- append(result.list, list(p)); p

# plot all in pannels
pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/relevance_etiology_biomarrker.pdf", width = 10, height = 6.5)
ggarrange(plotlist=result.list, ncol = 3, nrow = 2)
dev.off()

############################################################################################################
############################################################################################################
# endotype outcome biomarker
############################################################################################################
############################################################################################################
# read clustering result
cluster <- readRDS('kmeans.biomarker.death.shock.rds')
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

# read clinical data
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 1)
clinical[clinical=='.'] <- NA

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# reorder to match protein data
clinical <- clinical[match(ptid$ptid, clinical$ptid),]

# assign cluster 
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)

# relevant continuous variable
var.cont <- c('days2dth', 'temp', 'diabp', 'hrate','age')
result <- data.frame()

# create visualization data
for(v in var.cont){
  variable <- rep(v, dim(clinical)[1])
  var.1 <- as.numeric(clinical[clinical$cluster==2,v]);var.1
  var.2 <- as.numeric(clinical[clinical$cluster==1,v]);var.2
  endotype <- c(rep('E1\'', length(var.1)), rep('E2\'', length(var.2)))
  d <- data.frame(variable, endotype, value=c(var.1, var.2))
  result <- rbind(result, d)
}
result$endotype[result$endotype=='E1\''] <- 'E3\''
result$endotype[result$endotype=='E2\''] <- 'E4\''

############################################################################################################
# day2death
############################################################################################################
result.list <- list()
v <- 'days2dth'
r <- result[result$variable==v,]

p <- ggplot(r, aes(x=endotype, y=value, fill=endotype)) + 
  geom_boxplot(width=.35, outlier.shape = NA) + 
  theme_classic() +
  geom_jitter(shape=16, size=1, position=position_jitter(width = 0.1, height = 0, seed = 2023), aes(alpha=0.9)) +
  geom_signif(comparisons = list(c("E3\'", "E4\'")), 
              map_signif_level=c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05), 
              test = 'wilcox.test', textsize = 10, y_position = 95) +
  ylab('Day') +
  theme(axis.title.x=element_blank(), legend.position = "none", text = element_text(size = 19), 
        axis.title.y = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(1,0,0,1), 'lines'), axis.text = element_text(color = "black")) +
  scale_y_continuous(limits = c(10, 110), breaks = seq(10, 90, 20), guide = guide_prism_minor(), minor_breaks = seq(10, 90, 10)) +
  scale_fill_brewer(palette="Set2") + 
  ggtitle(v); result.list <- append(result.list, list(p)); p


############################################################################################################
# temp
############################################################################################################
v <- 'temp'
r <- result[result$variable==v,]
p <- ggplot(r, aes(x=endotype, y=value, fill=endotype)) + 
  geom_boxplot(width=.35, outlier.shape = NA) + 
  theme_classic() +
  geom_jitter(shape=16, size=1, position=position_jitter(width = 0.1, height = 0, seed = 2023), aes(alpha=0.9)) +
  geom_signif(comparisons = list(c("E3\'", "E4\'")), 
              map_signif_level=c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05), 
              test = 'wilcox.test', textsize = 10, y_position = 39.5) +
  ylab("Celsius") +
  theme(axis.title.x=element_blank(), legend.position = "none", text = element_text(size = 19), 
        axis.title.y = element_text(size=15, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(1,0,0,1), 'lines'), axis.text = element_text(color = "black")) +
  scale_y_continuous(limits = c(32, 40.5), guide = guide_prism_minor()) +
  scale_fill_brewer(palette="Set2") + 
  ggtitle(v); result.list <- append(result.list, list(p)); p


############################################################################################################
# diabp
############################################################################################################
v <- 'diabp'
r <- result[result$variable==v,]
p <- ggplot(r, aes(x=endotype, y=value, fill=endotype)) + 
  geom_boxplot(width=.35, outlier.shape = NA) + 
  theme_classic() +
  geom_jitter(shape=16, size=1, position=position_jitter(width = 0.1, height = 0, seed = 2023), aes(alpha=0.9)) +
  geom_signif(comparisons = list(c("E1\'", "E2\'")), 
              map_signif_level=c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05), 
              test = 'wilcox.test', textsize = 10, y_position = 92) +
  ylab("mmHg") +
  theme(axis.title.x=element_blank(), legend.position = "none", text = element_text(size = 19), 
        axis.title.y = element_text(size=15, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(2,1,2,1), 'lines'), axis.text = element_text(color = "black")) +
  scale_y_continuous(limits = c(33, 102), guide = guide_prism_minor()) +
  scale_fill_brewer(palette="Set2") + 
  ggtitle(v); result.list <- append(result.list, list(p)); p


############################################################################################################
# hrate
############################################################################################################
v <- 'hrate'
r <- result[result$variable==v,]
p <- ggplot(r, aes(x=endotype, y=value, fill=endotype)) + 
  geom_boxplot(width=.35, outlier.shape = NA) + 
  theme_classic() +
  geom_jitter(shape=16, size=1, position=position_jitter(width = 0.1, height = 0, seed = 2023), aes(alpha=0.9)) +
  geom_signif(comparisons = list(c("E1\'", "E2\'")), 
              map_signif_level=c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05), 
              test = 'wilcox.test', textsize = 10, y_position = 139) +
  ylab("mmHg") +
  theme(axis.title.x=element_blank(), legend.position = "none", text = element_text(size = 19), 
        axis.title.y = element_text(size=15, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(2,1,2,1), 'lines'), axis.text = element_text(color = "black")) +
  scale_y_continuous(limits = c(60, 150), guide = guide_prism_minor()) +
  scale_fill_brewer(palette="Set2") + 
  ggtitle(v); result.list <- append(result.list, list(p)); p


############################################################################################################  
# age
############################################################################################################
v <- 'age'
r <- result[result$variable==v,]
p <- ggplot(r, aes(x=endotype, y=value, fill=endotype)) + 
  geom_boxplot(width=.35, outlier.shape = NA) + 
  theme_classic() +
  geom_jitter(shape=16, size=1, position=position_jitter(width = 0.1, height = 0, seed = 2023), aes(alpha=0.9)) +
  geom_signif(comparisons = list(c("E3\'", "E4\'")), 
              map_signif_level=c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05), 
              test = 'wilcox.test', textsize = 10, y_position = 91) +
  ylab("Year") +
  theme(axis.title.x=element_blank(), legend.position = "none", text = element_text(size = 19), 
        axis.title.y = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(1,0,0,1), 'lines'), axis.text = element_text(color = "black")) +
  scale_y_continuous(limits = c(10, 105), breaks = seq(10, 90, 20), guide = guide_prism_minor(), minor_breaks = seq(10, 90, 10)) +
  scale_fill_brewer(palette="Set2") + 
  ggtitle(v); result.list <- append(result.list, list(p)); p


############################################################################################################
# barplot categorical
############################################################################################################
variable <- c(rep('death60', 4), rep('death90', 4), rep('shock', 4))
outcome <- c(rep(c('Death','Death','Survival','Survival'), 2), c('Shock','Shock','No shock','No shock'))
endotype <- rep(c('E3\'','E4\''), 6)
value <- c(29.63, 6.67, 70.37, 93.33, 37.04, 6.67, 62.96, 93.33, 70.37, 33.33, 29.63, 66.67)
result.cat <- data.frame(variable,outcome,endotype,value);result.cat

# death60
v <- 'death60'
r <- result.cat[result.cat$variable==v,]
p <- ggplot(r, aes(y=value, x=endotype)) + 
  geom_bar(aes(fill=outcome), position="fill", stat="identity", width=.35) +
  scale_y_continuous(limits = c(0, 1.25), labels=scales::percent_format(suffix = ""), 
                     oob = rescale_none, breaks = seq(0, 1, by=0.25), expand = c(0, 0)) +
  theme_classic() +
  geom_signif(aes(xmin = 'E3\'', xmax = 'E4\'', annotations = '*', y_position = 1.05),
              textsize = 10, tip_length = .0003) + 
  scale_fill_manual(values=c("#cc0000", "#56B4E9")) + 
  ggtitle(paste(v, '(%)')) + 
  theme(axis.title.x=element_blank(), text = element_text(size = 19), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(1,0,0,1), 'lines'), axis.text = element_text(color = "black"),
        legend.title=element_blank(),legend.position="right",
        legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=10),
        legend.margin = margin(t = 0, r = 10, b = 0, l = -5))
result.list <- append(result.list, list(p)); p

# death90
v <- 'death90'
r <- result.cat[result.cat$variable==v,]
p <- ggplot(r, aes(y=value, x=endotype)) + 
  geom_bar(aes(fill=outcome), position="fill", stat="identity", width=.35) +
  scale_y_continuous(limits = c(0, 1.25), labels=scales::percent_format(suffix = ""), 
                     oob = rescale_none, breaks = seq(0, 1, by=0.25), expand = c(0, 0)) +
  theme_classic() +
  geom_signif(aes(xmin = 'E3\'', xmax = 'E4\'', annotations = '**', y_position = 1.05),
              textsize = 10, tip_length = .0003) + 
  scale_fill_manual(values=c("#cc0000", "#56B4E9")) + 
  ggtitle(paste(v, '(%)')) + 
  theme(axis.title.x=element_blank(), text = element_text(size = 19), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(1,0,0,1), 'lines'), axis.text = element_text(color = "black"),
        legend.title=element_blank(),legend.position="right",
        legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=10),
        legend.margin = margin(t = 0, r = 10, b = 0, l = -5))
result.list <- append(result.list, list(p));p

# shock
v <- 'shock'
r <- result.cat[result.cat$variable==v,]
r$outcome <- factor(r$outcome, levels = c('Shock','No shock'))
p <- ggplot(r, aes(y=value, x=endotype)) + 
  geom_bar(aes(fill=outcome), position="fill", stat="identity", width=.35) +
  scale_y_continuous(limits = c(0, 1.25), labels=scales::percent_format(suffix = ""), 
                     oob = rescale_none, breaks = seq(0, 1, by=0.25), expand = c(0, 0)) +
  theme_classic() +
  geom_signif(aes(xmin = 'E3\'', xmax = 'E4\'', annotations = '**', y_position = 1.05),
              textsize = 10, tip_length = .0003) + 
  scale_fill_manual(values=c("#cc0000", "#56B4E9")) + 
  ggtitle(paste(v, '(%)')) + 
  theme(axis.title.x=element_blank(), text = element_text(size = 19), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(1,0,0,1), 'lines'), axis.text = element_text(color = "black"),
        legend.title=element_blank(),legend.position="right",
        legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=10),
        legend.margin = margin(t = 0, r = 10, b = 0, l = -5))
result.list <- append(result.list, list(p));p

# plot all in pannels
pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/relevance_outcome_biomarrker.pdf", width = 14, height = 7)
ggarrange(plotlist=result.list, ncol = 4, nrow = 2)
dev.off()


############################################################################################################
############################################################################################################
# variable contribution to endotype
############################################################################################################
############################################################################################################
############################################################################################################
# E12
############################################################################################################
result <- read.csv('var_contr_E12.csv')
result$protein <- factor(result$protein, levels=unique(result$protein))

# visualize
pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/var_contr_E12.pdf", width = 14, height = 6)
ggplot(result, aes(x=protein, y=value, group=endotype)) +
  geom_line(aes(color=endotype), size=1) +
  geom_point(size = 2, aes(color=endotype)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values=c('#66c2a4','#fc8d59'), labels = c("E1", "E2")) +
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size=15), legend.title=element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        legend.key = element_blank(), panel.background = element_blank()) + 
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5), expand = c(0,0)) +
  labs(x="Protein and Clinical Variables", y='Standardized Variable Value')
dev.off()

############################################################################################################
# E34
############################################################################################################
result <- read.csv('var_contr_E34.csv')
result$protein <- factor(result$protein, levels=unique(result$protein))

# visualize
pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/var_contr_E34.pdf", width = 14, height = 6)
ggplot(result, aes(x=protein, y=value, group=endotype)) +
  geom_line(aes(color=endotype), size=1) +
  geom_point(size = 2, aes(color=endotype)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values=c('#7fbc41','#bf812d'), labels = c("E3", "E4")) +
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size=15), legend.title=element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        legend.key = element_blank(), panel.background = element_blank()) + 
  scale_y_continuous(limits = c(-1.5, 1.5), breaks = seq(-1.5, 1.5, 0.5), expand = c(0,0)) +
  labs(x="Protein and Clinical Variables", y='Standardized Variable Value')
dev.off()


############################################################################################################
# E56
############################################################################################################
result <- read.csv('var_contr_E56.csv')
result$protein <- factor(result$protein, levels=unique(result$protein))

# visualize
pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/var_contr_E56.pdf", width = 6, height = 6)
ggplot(result, aes(x=protein, y=value, group=endotype)) +
  geom_line(aes(color=endotype), size=1) +
  geom_point(size = 2, aes(color=endotype)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values=c('#fee391','#fa9fb5'), labels = c("E5", "E6")) +
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size=15), legend.title=element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        legend.key = element_blank(), panel.background = element_blank()) + 
  scale_y_continuous(limits = c(-1.2, 1.2), expand = c(0,0)) +
  labs(x="Protein and Clinical Variables", y='Standardized Variable Value')
dev.off()


############################################################################################################
# E67
############################################################################################################
result <- read.csv('var_contr_E67.csv')
result$protein <- factor(result$protein, levels=unique(result$protein))

# visualize
pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/var_contr_E67.pdf", width = 7, height = 6)
ggplot(result, aes(x=protein, y=value, group=endotype)) +
  geom_line(aes(color=endotype), size=1) +
  geom_point(size = 2, aes(color=endotype)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values=c('#fa9fb5','grey'), labels = c("E6", "E7")) +
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size=15), legend.title=element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        legend.key = element_blank(), panel.background = element_blank()) + 
  scale_y_continuous(limits = c(-1.2, 1.2), expand = c(0,0)) +
  labs(x="Protein and Clinical Variables", y='Standardized Variable Value')
dev.off()

############################################################################################################
# E57
############################################################################################################
result <- read.csv('var_contr_E57.csv')
result$protein <- factor(result$protein, levels=unique(result$protein))

# visualize
pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/var_contr_E57.pdf", width = 7, height = 6)
ggplot(result, aes(x=protein, y=value, group=endotype)) +
  geom_line(aes(color=endotype), size=1) +
  geom_point(size = 2, aes(color=endotype)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values=c('#fee391','grey'), labels = c("E5", "E7")) +
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size=15), legend.title=element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        legend.key = element_blank(), panel.background = element_blank()) + 
  scale_y_continuous(limits = c(-1.2, 1.2), expand = c(0,0)) +
  labs(x="Protein and Clinical Variables", y='Standardized Variable Value')
dev.off()

############################################################################################################
# E89
############################################################################################################
result <- read.csv('var_contr_E89.csv')
result$protein <- factor(result$protein, levels=unique(result$protein))

# visualize
pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/var_contr_E89.pdf", width = 14, height = 6)
ggplot(result, aes(x=protein, y=value, group=endotype)) +
  geom_line(aes(color=endotype), size=1) +
  geom_point(size = 2, aes(color=endotype)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values=c("#c0dcc0", "#cc6677"), labels = c("E8", "E9")) +
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 45, hjust=1),
        text = element_text(size=15), legend.title=element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank(),
        legend.key = element_blank(), panel.background = element_blank()) + 
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5), expand = c(0,0)) +
  labs(x="Protein and Clinical Variables", y='Standardized Variable Value')
dev.off()

############################################################################################################
############################################################################################################
# survival analysis
############################################################################################################
############################################################################################################
library("survminer")
library("survival")

############################################################################################################
# survival curve
############################################################################################################
# E12
result <- read.csv('survival_E12.csv')

# construct KM curve
fit.surv <- survfit(Surv(time, status) ~ endotype, data = result)

# plot 
p <- ggsurvplot(fit.surv, conf.int = TRUE, censor=F, palette = c("#66c2a4", "#fc8d59"), size=1.5, font.tickslab=20,
                pval.size=10, legend.labs = c("E1", "E2"), legend = "right", legend.title = '', pval = T)
p <- p + xlab('Days') + ylab('Survival Probability'); p

# save to pdf
ggsave(file = "C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/surv_E12.pdf", units = 'in', width = 6, height = 4)


# E12
result <- read.csv('survival_E12_biomarker.csv')

# construct KM curve
fit.surv <- survfit(Surv(time, status) ~ endotype, data = result)

# plot 
p <- ggsurvplot(fit.surv, conf.int = TRUE, censor=F, palette = c("#66c2a4", "#fc8d59"), size=1.5,font.tickslab=20,
                pval.size=10, legend.labs = c("E3\'", "E4\'"), legend = "right", legend.title = '', pval = T)
p <- p + xlab('Days') + ylab('Survival Probability'); p

# save to pdf
ggsave(file = "C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/surv_E12_biomarker.pdf", units = 'in', width = 7, height = 4)


# E34
result <- read.csv('survival_E34.csv')

# construct KM curve
fit.surv <- survfit(Surv(time, status) ~ endotype, data = result)

# plot 
p <- ggsurvplot(fit.surv, conf.int = TRUE, censor=F, palette = c('#7fbc41','#bf812d'), size=0.75,
                legend.labs = c("E3", "E4"), legend = "right", legend.title = '', pval = T)
p <- p + xlab('Days') + ylab('Survival Probability'); p

# save to pdf
ggsave(file = "C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/surv_E34.pdf", units = 'in', width = 7, height = 4)


# E567
result <- read.csv('survival_E567.csv')

# construct KM curve
fit.surv <- survfit(Surv(time, status) ~ endotype, data = result)

# plot 
p <- ggsurvplot(fit.surv, conf.int = TRUE, censor=F, palette = c('#fee391','#fa9fb5', 'grey'), size=0.75,
                legend.labs = c("E5", "E6", 'E7'), legend = "right", legend.title = '', pval = T)
p <- p + xlab('Days') + ylab('Survival Probability'); p

# save to pdf
ggsave(file = "C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/surv_E567.pdf", units = 'in', width = 7, height = 4)


# E89
result <- read.csv('survival_E89.csv')

# construct KM curve
fit.surv <- survfit(Surv(time, status) ~ endotype, data = result)

# plot 
p <- ggsurvplot(fit.surv, conf.int = TRUE, censor=F, palette = c("#c0dcc0", "#cc6677"), size=0.75,
                legend.labs = c("E8", "E9"), legend = "right", legend.title = '', pval = T)
p <- p + xlab('Days') + ylab('Survival Probability'); p

# save to pdf
ggsave(file = "C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/surv_E89.pdf", units = 'in', width = 7, height = 4)

############################################################################################################
# discharge home curve
############################################################################################################
# E12
result <- read.csv('discharge_E12.csv')

# construct KM curve
fit.surv <- survfit(Surv(time, status) ~ endotype, data = result)

# plot 
p <- ggsurvplot(fit.surv, conf.int = TRUE, censor=F, palette = c("#66c2a4", "#fc8d59"), size=1.5, font.tickslab=20,
                pval.size=10, legend.labs = c("E1", "E2"), legend = "right", legend.title = '', pval = T, fun = 'event',
                pval.coord = c(0, 0.9))
p <- p + xlab('Days') + ylab('Discharge Home Probability'); p

# save to pdf
ggsave(file = "C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/hosp_E12.pdf", units = 'in', width = 6, height = 4)


# E34
result <- read.csv('discharge_E34.csv')

# construct KM curve
fit.surv <- survfit(Surv(time, status) ~ endotype, data = result)

# plot 
p <- ggsurvplot(fit.surv, conf.int = TRUE, censor=F, palette = c('#7fbc41','#bf812d'), size=0.75,
                legend.labs = c("E3", "E4"), legend = "right", legend.title = '', pval = T, fun = 'event',
                pval.coord = c(0, 0.9))
p <- p + xlab('Days') + ylab('Discharge Home Probability'); p

# save to pdf
ggsave(file = "C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/hosp_E34.pdf", units = 'in', width = 7, height = 4)


# E567
result <- read.csv('discharge_E567.csv')

# construct KM curve
fit.surv <- survfit(Surv(time, status) ~ endotype, data = result)

# plot 
p <- ggsurvplot(fit.surv, conf.int = TRUE, censor=F, palette = c('#fee391','#fa9fb5', 'grey'), size=0.75,
                legend.labs = c("E5", "E6", 'E7'), legend = "right", legend.title = '', pval = T, fun = 'event',
                pval.coord = c(0, 0.9))
p <- p + xlab('Days') + ylab('Discharge Home Probability'); p

# save to pdf
ggsave(file = "C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/hosp_E567.pdf", units = 'in', width = 7, height = 4)


# E89
result <- read.csv('discharge_E89.csv')

# construct KM curve
fit.surv <- survfit(Surv(time, status) ~ endotype, data = result)

# plot 
p <- ggsurvplot(fit.surv, conf.int = TRUE, censor=F, palette = c("#c0dcc0", "#cc6677"), size=0.75,
                legend.labs = c("E8", "E9"), legend = "right", legend.title = '', pval = T, fun = 'event',
                pval.coord = c(0, 0.9))
p <- p + xlab('Days') + ylab('Discharge Home Probability'); p

# save to pdf
ggsave(file = "C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/hosp_E89.pdf", units = 'in', width = 7, height = 4)

############################################################################################################
# discharge icu curve
############################################################################################################
# E12
result <- read.csv('icu_E12.csv')

# construct KM curve
fit.surv <- survfit(Surv(time, status) ~ endotype, data = result)

# plot 
p <- ggsurvplot(fit.surv, conf.int = TRUE, censor=F, palette = c("#66c2a4", "#fc8d59"), size=0.75,
                legend.labs = c("E1", "E2"), legend = "right", legend.title = '', pval = T, fun = 'event',
                pval.coord = c(0, 0.9))
p <- p + xlab('Days') + ylab('Discharge ICU Probability'); p

# save to pdf
ggsave(file = "C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/icu_E12.pdf", units = 'in', width = 7, height = 4)


# E34
result <- read.csv('icu_E34.csv')

# construct KM curve
fit.surv <- survfit(Surv(time, status) ~ endotype, data = result)

# plot 
p <- ggsurvplot(fit.surv, conf.int = TRUE, censor=F, palette = c('#7fbc41','#bf812d'), size=0.75,
                legend.labs = c("E3", "E4"), legend = "right", legend.title = '', pval = T, fun = 'event',
                pval.coord = c(0, 0.9))
p <- p + xlab('Days') + ylab('Discharge ICU Probability'); p

# save to pdf
ggsave(file = "C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/icu_E34.pdf", units = 'in', width = 7, height = 4)


# E567
result <- read.csv('icu_E567.csv')

# construct KM curve
fit.surv <- survfit(Surv(time, status) ~ endotype, data = result)

# plot 
p <- ggsurvplot(fit.surv, conf.int = TRUE, censor=F, palette = c('#fee391','#fa9fb5', 'grey'), size=0.75,
                legend.labs = c("E5", "E6", 'E7'), legend = "right", legend.title = '', pval = T, fun = 'event',
                pval.coord = c(0, 0.9))
p <- p + xlab('Days') + ylab('Discharge ICU Probability'); p

# save to pdf
ggsave(file = "C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/icu_E567.pdf", units = 'in', width = 7, height = 4)


# E89
result <- read.csv('icu_E89.csv')

# construct KM curve
fit.surv <- survfit(Surv(time, status) ~ endotype, data = result)

# plot 
p <- ggsurvplot(fit.surv, conf.int = TRUE, censor=F, palette = c("#c0dcc0", "#cc6677"), size=0.75,
                legend.labs = c("E8", "E9"), legend = "right", legend.title = '', pval = T, fun = 'event',
                pval.coord = c(0, 0.9))
p <- p + xlab('Days') + ylab('Discharge ICU Probability'); p

# save to pdf
ggsave(file = "C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/icu_E89.pdf", units = 'in', width = 7, height = 4)

# log ratio test
#survdiff(Surv(time, status) ~ endotype, data=result)


############################################################################################################
############################################################################################################
# biomarker E12 E34
############################################################################################################
############################################################################################################
library(cvms)
library(tibble)
library(caret)
library(PRROC)
library(Rtsne)

# read result
#prob.matrix <- read.csv('biomarker_endotype_E12_combine_rf_prob.csv') # E12
prob.matrix <- read.csv('biomarker_endotype_E34_combine_rf_prob.csv') # E34
prob.pred <- prob.matrix[,19]


# transform to binary prediction
y.pred <- ifelse(prob.pred>=0.5, 1, 0); table(y.pred)

# read patient clusters
#ptid <- read.csv('ptid_E1_E2.csv') # E12
ptid <- read.csv('ptid_E3_E4.csv') # E34

# transform endotype to class
#ptid$cluster <- ifelse(ptid$cluster=='E2', 1, 0); table(ptid$cluster) # E12
ptid$cluster <- ifelse(ptid$cluster=='E4', 1, 0); table(ptid$cluster) # E34
y.truth <- ptid$cluster

############################################################################################################
# confusion matrix
############################################################################################################
cm <- table(y.truth, y.pred); cm
cm[1,] <- cm[1,]/sum(cm[1,]); cm[2,] <- cm[2,]/sum(cm[2,])
cm <- round(cm, 3); cm 
cm <- as.data.frame(cm); cm

#pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/cm_E12.pdf")
pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/cm_E34.pdf")
ggplot(cm, aes(y.truth, y.pred, fill= Freq)) +
  geom_tile() + geom_text(aes(label=sprintf("%0.3f", round(Freq, digits = 3))), size=15) +
  scale_fill_gradient(low="whitesmoke", high="#009194") +
  labs(x = "Truth",y = "Prediction") +
  scale_x_discrete(labels=c("E1","E2"), expand = c(0, 0)) + # E12
  scale_y_discrete(labels=c("E1","E2"), expand = c(0, 0)) +
  #scale_x_discrete(labels=c("E3","E4"), expand = c(0, 0)) + # E34
  #scale_y_discrete(labels=c("E3","E4"), expand = c(0, 0)) +
  theme(legend.position="none", panel.background = element_blank(), 
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
        text = element_text(size=30),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
dev.off()

############################################################################################################
# AUC curve
############################################################################################################
# E12
output.auroc <- roc.curve(scores.class0 = prob.pred[y.truth==1], 
                          scores.class1 = prob.pred[y.truth==0], curve = TRUE)
pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/AUC_E12.pdf")
plot(output.auroc, legend=F, color=F, xlab='1-Specificity', main='E1 vs E2', 
     auc.main=F, cex.lab=2, cex.axis=2, cex.main=2, lwd=5); lines(x = c(0,1), y = c(0,1), lty=2, lwd=5); text('AUC = 0.947', x=0.6, y=0.2, cex=2)
dev.off()

# E34
output.auroc <- roc.curve(scores.class0 = prob.pred[y.truth==1], 
                          scores.class1 = prob.pred[y.truth==0], curve = TRUE)

pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/AUC_E34.pdf")
plot(output.auroc, legend=F, color=F, xlab='1-Specificity', main='E3 vs E4', 
     auc.main=F, cex.lab=2, cex.axis=2, cex.main=2, lwd=5); lines(x = c(0,1), y = c(0,1), lty=2, lwd=5); text('AUC = 0.998', x=0.6, y=0.2, cex=2)
dev.off()

############################################################################################################
# PCA plot
############################################################################################################
# read top proteins
result <- read.csv('biomarker_endotype_E12_combine_rf.csv') # E12
#result <- read.csv('biomarker_endotype_E34_combine_rf.csv') # E34
protein.top <- result$protein[1]
protein.top <- unlist(strsplit(protein.top, ","))

# read normalized DE data 
data.DE <- readRDS('data.DE.rds'); dim(data.DE)
#which(data.DE$gene.keep%in%protein.top)
data.DE <- data.DE[data.DE$gene.keep%in%protein.top,]; dim(data.DE)
data.DE <- data.DE[56:112]; dim(data.DE)

# PCA
pr.out <- prcomp(t(data.DE), center = T, scale. = T)
data.pca <- pr.out$x; dim(data.pca)
PVE.matrix <- summary(pr.out)$importance; PVE.matrix
PVE <- PVE.matrix[2,]
plot(PVE, xlab='Principle Components', ylab='Proportion of Variance Explained', cex.lab=1.5, pch = 19)

# plot E12
data.pca <- data.frame(data.pca[,1:2])
data.pca$label <- as.factor(ptid$cluster); table(data.pca$label)

pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/PCA_E12.pdf", height = 6)
ggplot(data.pca,aes(x=PC1,y=PC2, color=label)) + 
  geom_point(size=5) + 
  scale_color_manual(values=c('#66c2a4','#fc8d59'), labels = c("E1", "E2")) +
  xlim(-6.5, 6) + ylim(-4.5, 2.5) +
  theme(panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
        axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), legend.title=element_blank(),
        legend.text=element_text(size=25), legend.key = element_rect(fill = "transparent"), 
        axis.text=element_text(size=20)) +
  xlab('PC1 (40.6%)') + ylab('PC2 (10.1%)')
dev.off()

# plot E34
data.pca <- data.frame(data.pca[,1:2])
data.pca$label <- as.factor(ptid$cluster); table(data.pca$label)

pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/pca_E34.pdf", height = 6)
ggplot(data.pca,aes(x=PC1,y=PC2, color=label)) + 
  geom_point(size=5) + 
  scale_color_manual(values=c('#7fbc41','#bf812d'), labels = c("E1", "E2")) +
  xlim(-3.5, 6) + ylim(-1.6, 2.2) +
  theme(panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
        axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), legend.title=element_blank(),
        legend.text=element_text(size=25), legend.key = element_rect(fill = "transparent"),
        axis.text=element_text(size=20)) +
  xlab('PC1 (66.1%)') + ylab('PC2 (10.4%)')
dev.off()

############################################################################################################
############################################################################################################
# biomarker E567
############################################################################################################
############################################################################################################
library(cvms)
library(tibble)
library(caret)
library(PRROC)
library(Rtsne)

# read result
prob.matrix <- readRDS('biomarker_endotype_E567_combine_rf_prob.rds')
prob.pred <- prob.matrix[[17]]

# final endotype prediction
y.pred <- apply(prob.pred, 1, which.max)

# read patient clusters
ptid <- read.csv('ptid_E5_E6_E7.csv')

# transform endotype to class
mapping <- c('E5' = 1, 'E6' = 2, 'E7' = 3)
y.truth <- unname(mapping[ptid$cluster])
mean(y.pred==y.truth)

############################################################################################################
# confusion matrix
############################################################################################################
cm <- table(y.truth, y.pred); cm
cm[1,] <- cm[1,]/sum(cm[1,]); cm[2,] <- cm[2,]/sum(cm[2,]); cm[3,] <- cm[3,]/sum(cm[3,])
cm <- round(cm, 3); cm 
cm <- as.data.frame(cm); cm

svg(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/cm_E567.svg")
ggplot(cm, aes(y.truth, y.pred, fill= Freq)) +
  geom_tile() + geom_text(aes(label=sprintf("%0.3f", round(Freq, digits = 3))), size=7) +
  scale_fill_gradient(low="whitesmoke", high="#009194") +
  labs(x = "Truth",y = "Prediction") +
  scale_x_discrete(labels=c("E5","E6","E7"), expand = c(0, 0)) + 
  scale_y_discrete(labels=c("E5","E6","E7"), expand = c(0, 0)) +
  theme(legend.position="none", panel.background = element_blank(), 
        axis.ticks.x=element_blank(), axis.ticks.y=element_blank(),
        text = element_text(size=22),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
dev.off()

############################################################################################################
# AUC curve
############################################################################################################
# E5 vs others
y.truth.5 <- ifelse(y.truth==1,0,1)
prob.pred.5 <- prob.pred[,2] + prob.pred[,3]
auroc.5 <- roc.curve(scores.class0 = prob.pred.5[y.truth.5==1], 
                     scores.class1 = prob.pred.5[y.truth.5==0], curve = TRUE)

# E6 vs others
y.truth.6 <- ifelse(y.truth==2,0,1)
prob.pred.6 <- prob.pred[,1] + prob.pred[,3]
auroc.6 <- roc.curve(scores.class0 = prob.pred.6[y.truth.6==1], 
                     scores.class1 = prob.pred.6[y.truth.6==0], curve = TRUE)

# E7 vs others
y.truth.7 <- ifelse(y.truth==3,0,1)
prob.pred.7 <- prob.pred[,1] + prob.pred[,2]
auroc.7 <- roc.curve(scores.class0 = prob.pred.7[y.truth.7==1], 
                     scores.class1 = prob.pred.7[y.truth.7==0], curve = TRUE)

svg(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/AUC_E567.svg")
plot(auroc.5, legend=F, color='#fee391', xlab='1-Specificity', main='E5 vs E6 vs E7', 
     auc.main=F, cex.lab=1.5, cex.axis=1.5, lwd=2, cex.main=1.5); lines(x = c(0,1), y = c(0,1), lty=2)
plot(auroc.6, legend=F, add=T, lwd=2, color='#fa9fb5')
plot(auroc.7, legend=F, add=T, lwd=2, color='grey')
lines(x = c(0.35,0.45), y = c(0.3,0.3), col='#fee391', lwd=2); text('AUC (E5 vs Others) = 0.994', x=0.75, y=0.3, cex=1.2)
lines(x = c(0.35,0.45), y = c(0.2,0.2), col='#fa9fb5', lwd=2); text('AUC (E6 vs Others) = 0.907', x=0.75, y=0.2, cex=1.2)
lines(x = c(0.35,0.45), y = c(0.1,0.1), col='grey', lwd=2); text('AUC (E7 vs Others) = 0.930', x=0.75, y=0.1, cex=1.2)
text('Accuracy (E5 vs E6 vs E7) = 0.877', x=0.7, y=0, cex=1.2)
dev.off()

############################################################################################################
# PCA plot
############################################################################################################
# read top proteins
result <- read.csv('biomarker_endotype_E567_combine_rf.csv')
protein.top <- result$protein[1]
protein.top <- unlist(strsplit(protein.top, ","))

# read normalized DE data 
data.DE <- readRDS('data.DE.rds'); dim(data.DE)
data.DE <- data.DE[data.DE$gene.keep%in%protein.top,]; dim(data.DE)
data.DE <- data.DE[56:112]

# PCA
pr.out <- prcomp(t(data.DE), center = T, scale. = T)
data.pca <- pr.out$x; dim(data.pca)
PVE.matrix <- summary(pr.out)$importance; PVE.matrix
PVE <- PVE.matrix[2,]
plot(PVE, xlab='Principle Components', ylab='Proportion of Variance Explained', cex.lab=1.5, pch = 19)

# plot E567
data.pca <- data.frame(data.pca[,1:2])
data.pca$label <- as.factor(ptid$cluster)

svg(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/pca_E567.svg")
ggplot(data.pca,aes(x=PC1,y=PC2, color=label)) + 
  geom_point(size=3) + 
  scale_color_manual(values=c('#fee391','#fa9fb5', 'grey'), labels = c("E5", "E6", 'E7')) +
  #xlim(-6.5, 6) + ylim(-4.5, 2.5) +
  theme(panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA, linewidth = .5),
        axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.title=element_blank(),
        legend.text=element_text(size=17), legend.key = element_rect(fill = "transparent"),
        axis.text=element_text(size=15)) +
  xlab('PC1 (59.4%)') + ylab('PC2 (7.6%)')
dev.off()

############################################################################################################
############################################################################################################
# sankey plot
############################################################################################################
############################################################################################################
# read patient clusters
ptid.12 <- read.csv('ptid_E1_E2.csv')
ptid.34 <- read.csv('ptid_E3_E4.csv')
ptid.567 <- read.csv('ptid_E5_E6_E7.csv')
#ptid.89 <- read.csv('ptid_E8_E9.csv')
#ptid.56[ptid.56$cluster=='Not valid endotype', 'cluster'] <- 'Invalid'
ptid <- merge(ptid.12, ptid.34, by='ptid')
ptid <- merge(ptid, ptid.567, by='ptid')
#ptid <- merge(ptid, ptid.89, by='ptid')
ptid$freq <- 1

# read clinical data
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 1)
clinical[clinical=='.'] <- NA
ptid <- merge(ptid, clinical[,c('ptid', 'death90', 'cvp', 'gender')], by = 'ptid')

# outcome
ptid[ptid$death90==0, 'death90'] <- 'Survival'
ptid[ptid$death90=='1', 'death90'] <- 'Death'
ptid$death90 <- as.factor(ptid$death90)

# cvp
ptid$cvp <- as.numeric(ptid$cvp)
ptid$cvp <- as.character(cut(ptid$cvp, breaks = c(0, 12, 200), labels = c('Normal','Abnormal'), include.lowest = F))
ptid[is.na(ptid$cvp), 'cvp'] <- 'Missing'
ptid$cvp <- factor(ptid$cvp, levels = c('Normal','Abnormal','Missing'))

# gender
ptid[ptid$gender==2, 'gender'] <- 'Female'
ptid[ptid$gender=='1', 'gender'] <- 'Male'
ptid$gender <- as.factor(ptid$gender)

# read etiology
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 2)
clinical[clinical=='.'] <- NA
clinical$Pneumo <- ifelse(clinical$Pneumo==0, 0, 1)
ptid <- merge(ptid, clinical[,c('ptid', 'Pneumo')], by = 'ptid')
ptid[ptid$Pneumo==0, 'Pneumo'] <- 'No\n Pneumonia'
ptid[ptid$Pneumo=='1', 'Pneumo'] <- 'Pneumonia'
ptid$Pneumo <- as.factor(ptid$Pneumo)

# read severity
clinical.severity <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 17)
ptid <- merge(ptid, clinical.severity, by = 'ptid')
ptid$severity <- factor(ptid$severity, labels = c('Mild', 'Moderate', 'Severe'))

# reorder by severity
order.row <- ptid[order(ptid$severity, decreasing = F), 'ptid']
ptid$ptid <- factor(ptid$ptid, levels = order.row)

colnames(ptid)[2:4] <- c('endotype.12', 'endotype.34', 'endotype.567')
summary(ptid)

# E1 and E2
pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/Sankey_E1_E2.pdf", height = 9)
p <- ggplot(ptid, aes(y = freq, axis1=ptid, axis2 = endotype.12, axis3 = death90)) +
  geom_alluvium(width = 1/3, aes(fill = endotype.12)) +
  geom_stratum(width = 1/3, alpha = .0001) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Patient", "Endotype", 'Outcome'), expand = c(0, 0)) +
  scale_fill_manual(values=c('#66c2a4','#fc8d59')) + 
  theme(panel.background = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.title.y=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_text(size=13),
        legend.text=element_text(size=11), legend.title=element_blank(), plot.margin = unit(c(0, 0, 0, 1), "inches")) +
  scale_y_continuous(expand=c(0.01,0)); p
  #guides(fill=guide_legend(title="Endotype"))
dev.off()


# E3 and E4
pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/Sankey_E3_E4.pdf", height = 9, width = 8)
ggplot(ptid, aes(y = freq, axis1=ptid, axis2 = endotype.34, axis3 = Pneumo)) +
  geom_alluvium(width = 1/3, aes(fill = endotype.34)) +
  geom_stratum(width = 1/3, alpha = .0001) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Patient", "Endotype", 'Etiology'), expand = c(0, 0)) +
  scale_fill_manual(values=c('#7fbc41','#bf812d')) + 
  theme(panel.background = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.title.y=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_text(size=13),
        plot.margin = unit(c(0, 0, 0, 1), "inches"), legend.title=element_blank()) +
  scale_y_continuous(expand=c(0.01,0)) 
  #guides(fill=guide_legend(title="Endotype"))
dev.off()

# E5, E6, E7
pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/Sankey_E5_E6_E7.pdf", height = 9)
ggplot(ptid, aes(y = freq, axis1=ptid, axis2 = endotype.567, axis3 = cvp)) +
  geom_alluvium(width = 1/3, aes(fill = endotype.567)) +
  geom_stratum(width = 1/3, alpha = .0001) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Patient", "Endotype", 'CVP'), expand = c(0, 0)) +
  scale_fill_manual(values=c("#fee391", "#fa9fb5", 'grey')) + 
  theme(panel.background = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.title.y=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_text(size=13),
        plot.margin = unit(c(0, 0, 0, 1), "inches"), legend.title=element_blank()) +
  scale_y_continuous(expand=c(0.01,0)) 
  #guides(fill=guide_legend(title="Endotype"))
dev.off()

# E8 and E9
#svg(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/Sankey_E1_E2.svg", height = 9)
ggplot(ptid, aes(y = freq, axis1=ptid, axis2 = endotype.89, axis3 = gender)) +
  geom_alluvium(width = 1/3, aes(fill = endotype.89)) +
  geom_stratum(width = 1/3, alpha = .0001) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Patient", "Endotype", 'Gender'), expand = c(0, 0)) +
  scale_fill_manual(values=c("#c0dcc0", "#cc6677")) + 
  theme(panel.background = element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.title.y=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_text(size=13),
        legend.text=element_text(size=11), legend.title=element_blank(), plot.margin = unit(c(0, 0, 0, 1), "inches")) +
  scale_y_continuous(expand=c(0.01,0)) 
#guides(fill=guide_legend(title="Endotype"))
#dev.off()

############################################################################################################
############################################################################################################
# overlap plot
############################################################################################################
############################################################################################################
cluster.outcome <- readRDS('kmeans.imp.outcome.rds'); cluster.outcome
E1 <- cluster.outcome[[1]]
E2 <- cluster.outcome[[2]]

cluster.etiology <- readRDS('kmeans.imp.etiology.rds'); cluster.etiology
E3 <- cluster.etiology[[1]]
E4 <- cluster.etiology[[2]]

cluster.severity <- readRDS('kmeans.imp.severity.rds'); cluster.severity
E5 <- cluster.severity[[2]]
E6 <- cluster.severity[[3]]

# creat data frame for endotype assignment
data.venn <- list(E1=E1, E2=E2, E3=E3, E4=E4, E5=E5, E6=E6)
table.venn <- t(plyr::ldply(data.venn, rbind))
colnames(table.venn) <- table.venn[1,]; table.venn <- table.venn[-1,]

# match ptid
ptid <- read_excel("ptid.xlsx")
list.venn <- list(E1=ptid$ptid[data.venn$E1], E2=ptid$ptid[data.venn$E2], E3=ptid$ptid[data.venn$E3],
                  E4=ptid$ptid[data.venn$E4], E5=ptid$ptid[data.venn$E5], E6=ptid$ptid[data.venn$E6])
table.ptid <- data.frame(t(plyr::ldply(list.venn, rbind)))
colnames(table.ptid) <- table.ptid[1,]; table.ptid <- table.ptid[-1,]

#write.csv(table.ptid, 'endotype_ptid.csv', na = "", row.names = F)


# find overlap for each ptid
overlap.list <- list()
for(p in ptid$ptid){
  #p <- ptid$ptid[1]
  m <- unlist(lapply(table.ptid, function(x){
    p%in%x
  }))
  l <- list(names(which(m==T)))
  names(l) <- p
  overlap.list <- append(overlap.list, l)
}
overlap.table <- data.frame(t(plyr::ldply(overlap.list, rbind)))
colnames(overlap.table) <- overlap.table[1,]; overlap.table <- overlap.table[-1,]


overlap.table <- t(overlap.table)
heatmap.table <- data.frame()
heatmap.table <- data.frame(matrix(0, ncol = 6, nrow = 57))
colnames(heatmap.table) <- c('E1', 'E2', 'E3', 'E4', 'E5', 'E6')
rownames(heatmap.table) <- rownames(overlap.table)
split.column <- c('E1', 'E2', 'E3', 'E4', 'E5', 'E6')

for(i in 1:dim(heatmap.table)[1]){
  #i <- 3
  set.index <- unlist(sapply(overlap.table[i,], function(x){
    which(x==colnames(heatmap.table))
  }))
  heatmap.table[i, set.index] <- length(set.index) - 1
}

top = HeatmapAnnotation(
  type=anno_block(gp = gpar(fill = c('#66c2a4','#fc8d59','#7fbc41','#bf812d',"#fee391", "#fa9fb5")), 
                  labels = c('E1', 'E2', 'E3', 'E4', 'E5', 'E6'), labels_gp = gpar(fontface = "bold"))
)

hm <- Heatmap(as.matrix(heatmap.table), cluster_rows = T, cluster_columns = F, col=structure(c('white','yellow','red'), names = c("0", "1", "2")), 
              rect_gp = gpar(col = "black", lwd = 1), column_split = split.column, top_annotation = top, use_raster=F,
              column_title = NULL, show_column_names = F, row_names_side = "left", show_row_dend = F, column_gap = unit(0, "mm"),
              heatmap_legend_param = list(
                title = "", at = c(0, 1, 2), border = "black",  
                grid_height = unit(.7, "cm"), grid_width = unit(.7, "cm"),
                labels = c("No Overlapping", "Two-Endotype Overlapping", "Three-Endotype Overlapping"),
                labels_gp = gpar(fontsize = 11)
                )); hm <- draw(hm)

#write.csv(overlap.table, 'ptid_shared.csv', na = "", row.names = F)

############################################################################################################
############################################################################################################
# endotype heatmap
############################################################################################################
############################################################################################################
# read clean data
data.DE <- readRDS('data.DE.rds'); dim(data.DE)

# read significant proteins
cor.cont.imp <- read.csv('cor.cont.imp.csv')
cor.cat.imp <- read.csv('cor.cat.imp.csv')

# death90 is outcome
cor.cat.imp[cor.cat.imp$var=='death90', 'type'] <- 'outcome'

############################################################################################################
# outcome
############################################################################################################
gene.imp <- unique(c(cor.cont.imp$gene.keep[cor.cont.imp$type=='outcome'], cor.cat.imp$gene.keep[cor.cat.imp$type=='outcome']))
data.cor.imp <- data.DE[data.DE$gene.keep%in%gene.imp,]

# row name is protein
rownames(data.cor.imp) <- data.cor.imp$gene.keep

# up-regulated proteins
protein <- unlist(readRDS('protein.outcome.rds'))
data.protein <- data.cor.imp[data.cor.imp$gene.keep%in%protein,56:112]; dim(data.protein)

# row split
protein.1 <- readRDS('protein.outcome.rds')[[1]]
protein.2 <- readRDS('protein.outcome.rds')[[2]]
split.row <- rownames(data.protein)%in%protein.1
split.row[split.row==T] <- 'E1'
split.row[split.row=='FALSE'] <- 'E2'
table(split.row)

# col name is patient id
ptid <- read_excel("ptid.xlsx")
colnames(data.protein) <- ptid$ptid

# standardize
data.protein.norm <- t(scale(t(data.protein))); dim(data.protein)
sum(data.protein.norm[1,]); sd(data.protein.norm[1,])

cluster <- readRDS('kmeans.imp.outcome.rds'); cluster
split.column <- rep('E1', length(unlist(cluster)))
split.column[cluster[[2]]] <- 'E2'

#data.protein.norm <- data.protein.norm / max(abs(data.protein.norm))
data.protein.norm[data.protein.norm > 1] <- 1
data.protein.norm[data.protein.norm < -1] <- -1

top = HeatmapAnnotation(
  type=anno_block(gp = gpar(fill = c('#66c2a4', '#fc8d59')), labels = c('E1', 'E2'), labels_gp = gpar(fontface = "bold"))
)

left = rowAnnotation(
  type=anno_block(gp = gpar(fill = c('#66c2a4', '#fc8d59')), labels = c('E1 Up-regulated', 'E2 Up-regulated'), 
                  labels_gp = gpar(fontface = "bold"))
)

ht_opt$TITLE_PADDING = unit(c(15, 15), "points")
col_fun = colorRamp2(c(-1, 0, 1), c("red", "white", "yellow"))
column_lines <- columnAnnotation(
  line = c(2.5), 
  lwd = 2
)

set.seed(2023)
svg(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/heatmap_E1_E2.svg", width = 9, height = 5)
hm1 <- Heatmap(data.protein.norm, column_split = split.column, cluster_rows = T, cluster_columns = T, use_raster=F, 
        show_column_dend = F, show_row_dend = F,  column_title = NULL,
        clustering_distance_columns = "euclidean", name = "Abundance \n (Z-Score)", top_annotation = top,
        show_heatmap_legend = T, column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        width = ncol(data.protein.norm)*unit(3, "mm"), height = nrow(data.protein.norm)*unit(3, "mm"),
        show_row_names = T, show_column_names = T, column_gap = unit(2, "mm"), border = T,
        rect_gp = gpar(col = "white", lwd = 1), row_names_gp = grid::gpar(fontsize = 9),
        column_names_gp = grid::gpar(fontsize = 9), row_gap = unit(2, "mm"),
        row_title = NULL, col = col_fun, border_gp=gpar(lwd = 2),
        heatmap_legend_param = list(at = c(-1, -0.5, 0, 0.5, 1), labels = c("<-1", '-0.5', "0", '0.5', ">1"),
                                    labels_gp = gpar(fontsize = 12, fontface = "bold"), title_gp = gpar(fontsize = 12, fontface = "bold"))); hm1 <- draw(hm1)
dev.off()

############################################################################################################
# etiology
############################################################################################################
gene.imp <- cor.cat.imp$gene.keep[cor.cat.imp$type=='Etiology']
data.cor.imp <- data.DE[data.DE$gene.keep%in%gene.imp,]

# row name is protein
rownames(data.cor.imp) <- data.cor.imp$gene.keep

# up-regulated proteins
protein <- unlist(readRDS('protein.etiology.rds'))
data.protein <- data.cor.imp[data.cor.imp$gene.keep%in%protein,56:112]; dim(data.protein)

# row split
protein.1 <- readRDS('protein.etiology.rds')[[1]]
protein.2 <- readRDS('protein.etiology.rds')[[2]]
split.row <- rownames(data.protein)%in%protein.1
split.row[split.row==T] <- 'E3'
split.row[split.row=='FALSE'] <- 'E4'
table(split.row)

# col name is patient id
ptid <- read_excel("ptid.xlsx")
colnames(data.protein) <- ptid$ptid

# standardize
data.protein.norm <- t(scale(t(data.protein))); dim(data.protein)
sum(data.protein.norm[1,]); sd(data.protein.norm[1,])

cluster <- readRDS('kmeans.imp.etiology.rds'); cluster
split.column <- rep('E3', length(unlist(cluster)))
split.column[cluster[[2]]] <- 'E4'

data.protein.norm[data.protein.norm > 1] <- 1
data.protein.norm[data.protein.norm < -1] <- -1

top = HeatmapAnnotation(
  type=anno_block(gp = gpar(fill = c("#7fbc41", "#bf812d")), labels = c('E3', 'E4'), labels_gp = gpar(fontface = "bold"))
)

left = rowAnnotation(
  type=anno_block(gp = gpar(fill = c('#7fbc41', '#bf812d')), labels = c('E3 Up-regulated', 'E4 Up-regulated'), 
                  labels_gp = gpar(fontface = "bold"))
)

ht_opt$TITLE_PADDING = unit(c(15, 15), "points")
col_fun = colorRamp2(c(-1, 0, 1), c("red", "white", "yellow"))

set.seed(2023)
pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/heatmap_E3_E4.pdf", width = 9, height = 5.5)
hm2 <- Heatmap(data.protein.norm, column_split = split.column, cluster_rows = T, cluster_columns = T, use_raster=F, 
               show_column_dend = F, show_row_dend = F,  column_title = NULL,
               clustering_distance_columns = "euclidean", name = "Abundance \n (Z-Score)", top_annotation = top,
               show_heatmap_legend = T, column_title_gp = gpar(fontsize = 15, fontface = "bold"),
               width = ncol(data.protein.norm)*unit(3, "mm"), height = nrow(data.protein.norm)*unit(3, "mm"),
               show_row_names = T, show_column_names = T, column_gap = unit(2, "mm"), border = T,
               rect_gp = gpar(col = "white", lwd = 1), row_names_gp = grid::gpar(fontsize = 9),
               column_names_gp = grid::gpar(fontsize = 9), col = col_fun, row_title = NULL, border_gp=gpar(lwd = 2),
               row_gap = unit(2, "mm"),
               heatmap_legend_param = list(at = c(-1, -0.5, 0, 0.5, 1), labels = c("<-1", '-0.5', "0", '0.5', ">1"),
                                           labels_gp = gpar(fontsize = 12, fontface = "bold"), title_gp = gpar(fontsize = 12, fontface = "bold"))); hm2 <- draw(hm2)
dev.off()

############################################################################################################
# severity
############################################################################################################
gene.imp <- cor.cont.imp$gene.keep[cor.cont.imp$type=='severity']
data.cor.imp <- data.DE[data.DE$gene.keep%in%gene.imp,]

# row name is protein
rownames(data.cor.imp) <- data.cor.imp$gene.keep

# up-regulated proteins
protein <- unique(unlist(readRDS('protein.severity.rds')))
data.protein <- data.cor.imp[data.cor.imp$gene.keep%in%protein,56:112]; dim(data.protein)

# row split
protein.1 <- readRDS('protein.severity.rds')[[1]]
protein.2 <- readRDS('protein.severity.rds')[[2]]
protein.3 <- readRDS('protein.severity.rds')[[3]]
split.row <- rep('E5', dim(data.protein)[1])
split.row[which(rownames(data.protein)%in%protein.1)] <- 'E5'
split.row[which(rownames(data.protein)%in%protein.2)] <- 'E6'
split.row[which(rownames(data.protein)%in%protein.3)] <- 'E7'
# split.row <- rownames(data.protein)%in%protein.1
# split.row[split.row==T] <- 'E5'
# split.row[split.row=='FALSE'] <- 'E6'
table(split.row)

# col name is patient id
ptid <- read_excel("ptid.xlsx")
colnames(data.protein) <- ptid$ptid

# standardize
data.protein.norm <- t(scale(t(data.protein))); dim(data.protein)
sum(data.protein.norm[1,]); sd(data.protein.norm[1,])

cluster <- readRDS('kmeans.imp.severity.rds'); cluster
split.column <- rep('NE', length(unlist(cluster)))
split.column <- rep('E7', length(unlist(cluster)))
split.column[cluster[[2]]] <- 'E5'
split.column[cluster[[3]]] <- 'E6'
table(split.column)

#data.protein.norm <- data.protein.norm[,split.column!='NE']; dim(data.protein.norm)
#split.column <- split.column[split.column!='NE']; table(split.column)
#split.column <- factor(split.column, levels = c('E5', 'E6'))
split.column <- factor(split.column, levels = c('E5', 'E6', 'E7'))

data.protein.norm[data.protein.norm > 1] <- 1
data.protein.norm[data.protein.norm < -1] <- -1

#data.protein.norm <- data.protein.norm[which(rownames(data.protein.norm)!='GALM'),]

top = HeatmapAnnotation(
  type=anno_block(gp = gpar(fill = c("#fee391", "#fa9fb5", 'grey')), labels = c('E5', 'E6', 'E7'), labels_gp = gpar(fontface = "bold"))
)

left = rowAnnotation(
  type=anno_block(gp = gpar(fill = c('#fee391', '#fa9fb5', 'grey')), labels = c('E5 Up-\n regulated', 'E6 Up- \n regulated', 'E7 Up- \n regulated'), 
                  labels_gp = gpar(fontsize=9, fontface = "bold"))
)

ht_opt$TITLE_PADDING = unit(c(15, 15), "points")
col_fun = colorRamp2(c(-1, 0, 1), c("red", "white", "yellow"))

set.seed(2023)
svg(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/heatmap_E5_E6.svg", width = 7, height = 3.5)
hm3 <- Heatmap(data.protein.norm, column_split = split.column, cluster_rows = F, cluster_columns = T, use_raster=F, 
               show_column_dend = F, show_row_dend = F,  column_title = NULL,
               clustering_distance_columns = "euclidean", name = "Abundance \n (Z-Score)", top_annotation = top,
               show_heatmap_legend = T, column_title_gp = gpar(fontsize = 15, fontface = "bold"),
               width = ncol(data.protein.norm)*unit(3, "mm"), height = nrow(data.protein.norm)*unit(3, "mm"),
               show_row_names = T, show_column_names = T, column_gap = unit(2, "mm"), border = T,
               rect_gp = gpar(col = "white", lwd = 1), row_names_gp = grid::gpar(fontsize = 9),
               column_names_gp = grid::gpar(fontsize = 9), col = col_fun,row_title = NULL, border_gp=gpar(lwd = 2),
               row_gap = unit(2, "mm"), cluster_column_slices = F, 
               row_order = c("PIK3CD","SF3A2","AIF1" ,"PSMB8","ME1","CRNN","MARCKS","C6orf132", "KTN1","GALM","IGHV3-72", "IGHV1-18", "IGKV2-24","PCLO", 'NAP1L1', 'ASPSCR1'),
               heatmap_legend_param = list(at = c(-1, -0.5, 0, 0.5, 1), labels = c("<-1", '-0.5', "0", '0.5', ">1"),
                                           labels_gp = gpar(fontsize = 12, fontface = "bold"), title_gp = gpar(fontsize = 12, fontface = "bold"))); hm3 <- draw(hm3)
dev.off()

############################################################################################################
# baseline
############################################################################################################
gene.imp <- unique(c(cor.cont.imp$gene.keep[cor.cont.imp$type=='baseline'], cor.cat.imp$gene.keep[cor.cat.imp$type=='baseline']))
data.cor.imp <- data.DE[data.DE$gene.keep%in%gene.imp,]

# row name is protein
rownames(data.cor.imp) <- data.cor.imp$gene.keep

# up-regulated proteins
protein <- unlist(readRDS('protein.baseline.rds'))
data.protein <- data.cor.imp[data.cor.imp$gene.keep%in%protein,56:112]; dim(data.protein)

# row split
protein.1 <- readRDS('protein.baseline.rds')[[1]]
protein.2 <- readRDS('protein.baseline.rds')[[2]]
split.row <- rownames(data.protein)%in%protein.1
split.row[split.row==T] <- 'E8'
split.row[split.row=='FALSE'] <- 'E9'
table(split.row)

# col name is patient id
ptid <- read_excel("ptid.xlsx")
colnames(data.protein) <- ptid$ptid

# standardize
data.protein.norm <- t(scale(t(data.protein))); dim(data.protein)
sum(data.protein.norm[1,]); sd(data.protein.norm[1,])

cluster <- readRDS('kmeans.imp.baseline.2.rds'); cluster
split.column <- rep('E8', length(unlist(cluster)))
split.column[cluster[[2]]] <- 'E9'
table(split.column)

#data.protein.norm <- data.protein.norm / max(abs(data.protein.norm))
data.protein.norm[data.protein.norm > 1] <- 1
data.protein.norm[data.protein.norm < -1] <- -1

top = HeatmapAnnotation(
  type=anno_block(gp = gpar(fill = c("#c0dcc0", "#cc6677")), labels = c('E8', 'E9'), labels_gp = gpar(fontface = "bold"))
)

left = rowAnnotation(
  type=anno_block(gp = gpar(fill = c("#c0dcc0", "#cc6677")), labels = c('E8 Up-regulated', 'E9 Up-regulated'), 
                  labels_gp = gpar(fontface = "bold"))
)

ht_opt$TITLE_PADDING = unit(c(15, 15), "points")
col_fun = colorRamp2(c(-1, 0, 1), c("red", "white", "yellow"))
column_lines <- columnAnnotation(
  line = c(2.5), 
  lwd = 2
)

set.seed(2023)
#svg(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/heatmap_E1_E2.svg", width = 9, height = 5)
hm4 <- Heatmap(data.protein.norm, column_split = split.column, cluster_rows = T, cluster_columns = T, use_raster=F, 
               show_column_dend = F, show_row_dend = F,  column_title = NULL,
               clustering_distance_columns = "euclidean", name = "Abundance \n (Z-Score)", top_annotation = top,
               show_heatmap_legend = T, column_title_gp = gpar(fontsize = 15, fontface = "bold"),
               width = ncol(data.protein.norm)*unit(3, "mm"), height = nrow(data.protein.norm)*unit(3, "mm"),
               show_row_names = T, show_column_names = T, column_gap = unit(2, "mm"), border = T,
               rect_gp = gpar(col = "white", lwd = 1), row_names_gp = grid::gpar(fontsize = 9),
               column_names_gp = grid::gpar(fontsize = 9), row_gap = unit(2, "mm"),
               row_title = NULL, col = col_fun, border_gp=gpar(lwd = 2),
               heatmap_legend_param = list(at = c(-1, -0.5, 0, 0.5, 1), labels = c("<-1", '-0.5', "0", '0.5', ">1"),
                                           labels_gp = gpar(fontsize = 12, fontface = "bold"), title_gp = gpar(fontsize = 12, fontface = "bold"))); hm4 <- draw(hm4)
#dev.off()





hm3@ht_list[["Abundance 
 (Z-Score)"]]@row_order

# galm
cluster[[2]]
p <- data.protein.norm[which(rownames(data.protein.norm)=='GALM'),]
p5 <- p[which(split.column=='E5')]; mean(p5)
p6 <- p[which(split.column=='E6')]; mean(p6)

p <- data.protein[which(rownames(data.protein.norm)=='GALM'),]
p5 <- unlist(p[cluster[[2]]]); mean(p5)
p6 <- unlist(p[cluster[[3]]]); mean(p6)



# plot all in pannels
grid.newpage()
pushViewport(viewport(x = 0, y=.7, width = 0.5, just = c("left")))
draw(hm1, newpage = FALSE)
upViewport()

pushViewport(viewport(x = 0, y = .5, width = 0.5, height = 0.5, just = c("left", "top")))
draw(hm3, newpage = FALSE)
popViewport()

pushViewport(viewport(x = .5, y = .882, width = 0.5, height = 0.5, just = c("left", "top")))
draw(hm2, newpage = FALSE)
popViewport()


############################################################################################################
############################################################################################################
# endotype Etiology
############################################################################################################
############################################################################################################
# read clustering result
cluster <- readRDS('kmeans.imp.etiology.rds')
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

# read clinical data
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 5, na = c('.','.mild'))

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# reorder to match protein data
clinical <- clinical[match(ptid$ptid, clinical$ptid),]

# assign cluster 
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)

# continuous variable
var.cont <- c('pf2', 'pf3', 'pf5')

# create visualization data
result <- data.frame()
for(v in var.cont){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]);var.1
  var.2 <- as.numeric(clinical[clinical$cluster==2,v]);var.2
  endotype <- c(rep('E3', length(var.1)), rep('E4', length(var.2)))
  variable <- rep(v, length(var.1)+length(var.2))
  d <- data.frame(variable, endotype, value=c(var.1, var.2))
  result <- rbind(result, d)
}

# pf ratio
result.list <- list()

#pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/pf_etiology.pdf")
p <- ggplot(result, aes(x=variable, y=value, fill=endotype)) + 
  geom_boxplot(width=.5, outlier.shape = NA) + 
  theme_classic() + ylab('pf ratio') + xlab("Day") +
  geom_point(size=.6, position=position_jitterdodge(jitter.width = 0.05, jitter.height = 0, seed = 2023),alpha=0.5) +
  scale_fill_manual(values=c("#7fbc41", "#bf812d")) + 
  theme(text = element_text(size = 25), 
        axis.title.y = element_blank(), plot.margin = unit(c(2,1,2,1), 'lines'),
        axis.text = element_text(color = "black"), legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18, margin = margin(0,0,10,0)), 
        legend.margin = margin(t = 0, r = 10, b = 0, l = -5)) +
  geom_signif(y_position = c(400,400,400), xmin = c(0.8, 1.8, 2.8), xmax = c(1.2, 2.2, 3.2), 
              annotation = c("**", "*", '*'), textsize = c(10), tip_length = .01) +
  scale_y_continuous(limits = c(0, 450), breaks = seq(0, 450, 100), guide = guide_prism_minor(), minor_breaks = seq(0, 400, 50),
                     expand = c(0,0)) +
  guides(fill = guide_legend(override.aes = list(linetype = 0, shape=NA))) +
  ggtitle("P/F ratio (mmHg)") +
  scale_x_discrete(labels=c("pf5" = "5", "pf2" = "2", "pf3" = "3")); p
#dev.off()

result.list <- append(result.list, list(p))

############################################################################################################
# Pneumo
############################################################################################################
variable <- rep('Pneumo', 4)
outcome <- rep(c('Pneumo', 'No Pneumo'), 2)
endotype <- c('E3','E3','E4','E4')
value <- c(91.66, 8.34, 55.55, 44.45)
result.cat <- data.frame(variable,outcome,endotype,value)
result.cat$outcome <- factor(result.cat$outcome, levels=c('Pneumo','No Pneumo'))

v <- 'Pneumo'
r <- result.cat[result.cat$variable==v,]
p <- ggplot(r, aes(y=value, x=endotype)) + 
  geom_bar(aes(fill=outcome), position="fill", stat="identity", width=.35) +
  scale_y_continuous(limits = c(0, 1.25), labels=scales::percent_format(suffix = ""), 
                     oob = rescale_none, breaks = seq(0, 1, by=0.25), expand = c(0, 0)) +
  theme_classic() +
  geom_signif(aes(xmin = "E3", xmax = "E4", annotations = '*', y_position = 1.05),
              textsize = 10, tip_length = .0003) + 
  scale_fill_manual(values=c("#cc0000", "#56B4E9")) + 
  ggtitle(paste(v, '(%)')) + 
  theme(axis.title.x=element_blank(), text = element_text(size = 25), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(2,0,2,1), 'lines'), axis.text = element_text(color = "black"),
        legend.title=element_blank(),legend.position="right",
        legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=10),
        legend.margin = margin(t = 0, r = 10, b = 0, l = -5)); p
result.list <- append(result.list, list(p))

# plot all in pannels
pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/relevance_etiology.pdf", height = 5, width = 9)
ggarrange(plotlist=result.list)
dev.off()

############################################################################################################
############################################################################################################
# endotype severity
############################################################################################################
############################################################################################################
# read clustering result
cluster <- readRDS('kmeans.imp.severity.rds')
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

# read clinical data
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 5, na = c('.','.mild'))

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# reorder to match protein data
clinical <- clinical[match(ptid$ptid, clinical$ptid),]

# assign cluster 
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)

# continuous variable
var.cont <- c('pf0', 'pf1', 'pf2', 'pf3', 'pf4', 'pf6')

# create visualization data
result <- data.frame()
for(v in var.cont){
  var.1 <- as.numeric(clinical[clinical$cluster==2,v]);var.1
  var.2 <- as.numeric(clinical[clinical$cluster==3,v]);var.2
  var.3 <- as.numeric(clinical[clinical$cluster==1,v]);var.3
  endotype <- c(rep('E5', length(var.1)), rep('E6', length(var.2)), rep('E7', length(var.3)))
  variable <- rep(v, length(var.1)+length(var.2)+length(var.3))
  d <- data.frame(variable, endotype, value=c(var.1, var.2, var.3))
  result <- rbind(result, d)
}

############################################################################################################
# pf ratio
############################################################################################################
result.list <- list()
p <- ggplot(result, aes(x=variable, y=value, fill=endotype)) + 
  geom_boxplot(width=.5, outlier.shape = NA) + 
  theme_classic() + ylab('pf ratio') + xlab("Day") +
  geom_point(size=.6, position=position_jitterdodge(jitter.width = 0.05, jitter.height = 0, seed = 2023),alpha=0.5) +
  scale_fill_manual(values=c("#fee391", "#fa9fb5", 'grey')) + 
  theme(text = element_text(size = 19), 
        axis.title.y = element_blank(), plot.margin = unit(c(2,1,2,1), 'lines'),
        axis.text = element_text(color = "black"), legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18, margin = margin(0,0,10,0)), 
        legend.margin = margin(t = 0, r = 10, b = 0, l = -5)) +
  geom_signif(y_position = c(500,500,500,500,500,500,500,500), 
              xmin = c(0.8,1.8,2.8,3.01,3.8,4.01,5,6), 
              xmax = c(1,2,2.99,3.2,3.99,4.2,5.2,6.2), 
              annotation = c('**','*','*','*','*','*','**','*'), 
              textsize = c(10), 
              tip_length = .01) +
  scale_y_continuous(limits = c(0, 600), breaks = seq(0, 500, 100), guide = guide_prism_minor(), minor_breaks = seq(0, 500, 50),
                     expand = c(0,0)) +
  guides(fill = guide_legend(override.aes = list(linetype = 0, shape=NA))) +
  ggtitle("P/F ratio (mmHg)") +
  scale_x_discrete(labels=c("pf0" = "0", "pf1" = "1", "pf2" = "2", "pf3" = "3", "pf4" = "4", "pf6" = "6")); p

result.list <- append(result.list, list(p))

############################################################################################################
# cvp
############################################################################################################
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 1)
clinical[clinical=='.'] <- NA
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)

var.cont <- c('cvp')
result <- data.frame()
for(v in var.cont){
  var.1 <- as.numeric(clinical[clinical$cluster==2,v]);var.1
  var.2 <- as.numeric(clinical[clinical$cluster==3,v]);var.2
  var.3 <- as.numeric(clinical[clinical$cluster==1,v]);var.3
  endotype <- c(rep('E5', length(var.1)), rep('E6', length(var.2)), rep('E7', length(var.3)))
  variable <- rep(v, length(var.1)+length(var.2)+length(var.3))
  d <- data.frame(variable, endotype, value=c(var.1, var.2, var.3))
  result <- rbind(result, d)
}

v <- 'cvp'
r <- result[result$variable==v,]
p <- ggplot(r, aes(x=endotype, y=value, fill=endotype)) + 
  geom_boxplot(outlier.shape = NA, width=0.35) + 
  theme_classic() +
  geom_jitter(shape=16, size=1, position=position_jitter(width = 0.1, height = 0, seed = 2023), aes(alpha=0.9)) +
  geom_signif(comparisons = list(c("E5","E6", 'E7')), 
              map_signif_level=c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05), 
              test = 'wilcox.test', textsize = 10, y_position = 40) +
  ylab("mmHg") +
  theme(axis.title.x=element_blank(), legend.position = "none", text = element_text(size = 19), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(2,1,2,1), 'lines'), axis.text = element_text(color = "black")) +
  scale_y_continuous(limits = c(0, 45), breaks = seq(0, 45, 10), guide = guide_prism_minor(), minor_breaks = seq(0, 40, 5)) +
  scale_fill_manual(values=c("#fee391", "#fa9fb5", 'grey')) + 
  ggtitle('CVP (mmHg)'); p

result.list <- append(result.list, list(p))

# plot all in pannels
svg(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/relevance_severity.svg", height = 5, width = 9)
ggarrange(plotlist=result.list)
dev.off()

############################################################################################################
############################################################################################################
# endotype outcome
############################################################################################################
############################################################################################################
# read clustering result
cluster <- readRDS('kmeans.imp.outcome.rds')
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

# relevant continuous variable
var.cont <- c('days2dth',"hepatic28", 'icufd', 'vfd', 'coag28', 'cns28', 'cardio28', 'renal28', 'apache')
result <- data.frame()

# create visualization data
for(v in var.cont){
  variable <- rep(v, dim(clinical)[1])
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]);var.1
  var.2 <- as.numeric(clinical[clinical$cluster==2,v]);var.2
  endotype <- c(rep('E1', length(var.1)), rep('E2', length(var.2)))
  d <- data.frame(variable, endotype, value=c(var.1, var.2))
  result <- rbind(result, d)
}
result$endotype[result$endotype=='E1'] <- 'E3'
result$endotype[result$endotype=='E2'] <- 'E4'

############################################################################################################
# day2death
############################################################################################################
result.list <- list()
v <- 'days2dth'
r <- result[result$variable==v,]

p <- ggplot(r, aes(x=endotype, y=value, fill=endotype)) + 
  geom_boxplot(width=.35, outlier.shape = NA) + 
  theme_classic() +
  geom_jitter(shape=16, size=1, position=position_jitter(width = 0.1, height = 0, seed = 2023), aes(alpha=0.9)) +
  geom_signif(comparisons = list(c("E3", "E4")), 
              map_signif_level=c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05), 
              test = 'wilcox.test', textsize = 10, y_position = 90) +
  ylab('Day') +
  theme(axis.title.x=element_blank(), legend.position = "none", text = element_text(size = 19), 
        axis.title.y = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(1,0,0,1), 'lines'), axis.text = element_text(color = "black")) +
  scale_y_continuous(limits = c(10, 110), breaks = seq(10, 90, 20), guide = guide_prism_minor(), minor_breaks = seq(10, 90, 10)) +
  scale_fill_brewer(palette="Set2") + 
  ggtitle(v); p; result.list <- append(result.list, list(p))


############################################################################################################
# XXX free days
############################################################################################################
var.cont <- c("hepatic28", 'icufd', 'vfd', 'coag28', 'cns28', 'cardio28', 'renal28')

for(i in 1:length(var.cont)){
  v <- var.cont[i]
  r <- result[result$variable==v,]
  p <- ggplot(r, aes(x=endotype, y=value, fill=endotype)) + 
    geom_boxplot(width=.35, outlier.shape = NA) + 
    theme_classic() +
    geom_jitter(shape=16, size=1, position=position_jitter(width = 0.1, height = 0, seed = 2023), aes(alpha=0.9)) +
    geom_signif(comparisons = list(c("E3", "E4")), 
                map_signif_level=c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05), 
                test = 'wilcox.test', textsize = 10, y_position = 30) +
    ylab('Day') +
    theme(axis.title.x=element_blank(), legend.position = "none", text = element_text(size = 19), 
          axis.title.y = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          plot.title = element_text(hjust = 0.5, size = 18),
          plot.margin = unit(c(1,0,0,1), 'lines'), axis.text = element_text(color = "black")) +
    scale_y_continuous(limits = c(0, 37), breaks = seq(0, 30, 10), guide = guide_prism_minor(), minor_breaks = seq(0, 30, 5)) +
    scale_fill_brewer(palette="Set2") + 
    ggtitle(v); p
  result.list <- append(result.list, list(p))
}

############################################################################################################
# apache
############################################################################################################
v <- 'apache'
r <- result[result$variable==v,]

p <- ggplot(r, aes(x=endotype, y=value, fill=endotype)) + 
  geom_boxplot(width=.35, outlier.shape = NA) + 
  theme_classic() +
  geom_jitter(shape=16, size=1, position=position_jitter(width = 0.1, height = 0, seed = 2023), aes(alpha=0.9)) +
  geom_signif(comparisons = list(c("E3", "E4")), 
              map_signif_level=c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05), 
              test = 'wilcox.test', textsize = 10, y_position = 150) +
  ylab("Score") +
  theme(axis.title.x=element_blank(), legend.position = "none", text = element_text(size = 19), 
        axis.title.y = element_text(size=15, margin = margin(t = 0, r = 5, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(1,0,0,1), 'lines'), axis.text = element_text(color = "black")) +
  scale_y_continuous(limits = c(40, 170), breaks = seq(40, 160, 40), guide = guide_prism_minor(), minor_breaks = seq(40, 160, 20)) +
  scale_fill_brewer(palette="Set2") + 
  ggtitle(v); result.list <- append(result.list, list(p))

############################################################################################################  
# age
############################################################################################################
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 1)
clinical[clinical=='.'] <- NA
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)

var.cont <- c('age')
result <- data.frame()

for(v in var.cont){
  variable <- rep(v, dim(clinical)[1])
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]);var.1
  var.2 <- as.numeric(clinical[clinical$cluster==2,v]);var.2
  endotype <- c(rep('E1', length(var.1)), rep('E2', length(var.2)))
  d <- data.frame(variable, endotype, value=c(var.1, var.2))
  result <- rbind(result, d)
}
result$endotype[result$endotype=='E1'] <- 'E3'
result$endotype[result$endotype=='E2'] <- 'E4'

v <- 'age'
r <- result[result$variable==v,]

p <- ggplot(r, aes(x=endotype, y=value, fill=endotype)) + 
  geom_boxplot(width=.35, outlier.shape = NA) + 
  theme_classic() +
  geom_jitter(shape=16, size=1, position=position_jitter(width = 0.1, height = 0, seed = 2023), aes(alpha=0.9)) +
  geom_signif(comparisons = list(c("E3", "E4")), 
              map_signif_level=c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05), 
              test = 'wilcox.test', textsize = 10, y_position = 90) +
  ylab("Year") +
  theme(axis.title.x=element_blank(), legend.position = "none", text = element_text(size = 19), 
        axis.title.y = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(1,0,0,1), 'lines'), axis.text = element_text(color = "black")) +
  scale_y_continuous(limits = c(10, 105), breaks = seq(10, 90, 20), guide = guide_prism_minor(), minor_breaks = seq(10, 90, 10)) +
  scale_fill_brewer(palette="Set2") + 
  ggtitle(v); p; result.list <- append(result.list, list(p))


############################################################################################################
# barplot categorical
############################################################################################################
variable <- c(rep('death60', 4), rep('death90', 4))
outcome <- rep(c('Death','Death','Survival','Survival'), 2)
endotype <- rep(c('E3','E4'), 4)
value <- c(47.36, 2.63, 52.64, 97.37, 47.63, 7.89, 52.37, 92.11)
result.cat <- data.frame(variable,outcome,endotype,value)

# death60
v <- 'death60'
r <- result.cat[result.cat$variable==v,]
p <- ggplot(r, aes(y=value, x=endotype)) + 
  geom_bar(aes(fill=outcome), position="fill", stat="identity", width=.35) +
  scale_y_continuous(limits = c(0, 1.25), labels=scales::percent_format(suffix = ""), 
                     oob = rescale_none, breaks = seq(0, 1, by=0.25), expand = c(0, 0)) +
  theme_classic() +
  geom_signif(aes(xmin = "E3", xmax = "E4", annotations = '***', y_position = 1.05),
              textsize = 10, tip_length = .0003) + 
  scale_fill_manual(values=c("#cc0000", "#56B4E9")) + 
  ggtitle(paste(v)) + 
  ylab("Outcome (%)") +
  theme(axis.title.x=element_blank(), text = element_text(size = 19), 
        #axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(1,0,0,1), 'lines'), axis.text = element_text(color = "black"),
        legend.title=element_blank(),legend.position="none",
        legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=10),
        legend.margin = margin(t = 0, r = 10, b = 0, l = -5)); p
result.list <- append(result.list, list(p))

# death90
v <- 'death90'
r <- result.cat[result.cat$variable==v,]
p <- ggplot(r, aes(y=value, x=endotype)) + 
  geom_bar(aes(fill=outcome), position="fill", stat="identity", width=.35) +
  scale_y_continuous(limits = c(0, 1.25), labels=scales::percent_format(suffix = ""), 
                     oob = rescale_none, breaks = seq(0, 1, by=0.25), expand = c(0, 0)) +
  theme_classic() +
  geom_signif(aes(xmin = "E3", xmax = "E4", annotations = '**', y_position = 1.05),
              textsize = 10, tip_length = .0003) + 
  scale_fill_manual(values=c("#cc0000", "#56B4E9")) + 
  ggtitle(paste(v)) + 
  ylab("Outcome (%)") +
  theme(axis.title.x=element_blank(), text = element_text(size = 19), 
        #axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(1,0,0,1), 'lines'), axis.text = element_text(color = "black"),
        legend.title=element_blank(),legend.position="none",
        legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=10),
        legend.margin = margin(t = 0, r = 10, b = 0, l = -5))
result.list <- append(result.list, list(p))

# plot all in pannels
pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/relevance_outcome.pdf", width = 11, height = 9)
ggarrange(plotlist=result.list)
dev.off()


############################################################################################################
############################################################################################################
# endotype baseline
############################################################################################################
############################################################################################################
# read clustering result
cluster <- readRDS('kmeans.imp.baseline.2.rds')
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
result.list <- list()


############################################################################################################  
# pbw height packyr temp pao2screen
############################################################################################################
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 1)
clinical[clinical=='.'] <- NA
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# reorder to match protein data
clinical <- clinical[match(ptid$ptid, clinical$ptid),]


var.cont <- c('pbw', 'height', 'packyr','temp','pao2screen')
result <- data.frame()

for(v in var.cont){
  variable <- rep(v, dim(clinical)[1])
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]);var.1
  var.2 <- as.numeric(clinical[clinical$cluster==2,v]);var.2
  endotype <- c(rep('E8', length(var.1)), rep('E9', length(var.2)))
  d <- data.frame(variable, endotype, value=c(var.1, var.2))
  result <- rbind(result, d)
}

# pbw
r <- result[result$variable=='pbw',]
p <- ggplot(r, aes(x=endotype, y=value, fill=endotype)) + 
  geom_boxplot(width=.35, outlier.shape = NA) + 
  theme_classic() +
  geom_jitter(shape=16, size=1, position=position_jitter(width = 0.1, height = 0, seed = 2023), aes(alpha=0.9)) +
  geom_signif(comparisons = list(c("E8", "E9")), 
              map_signif_level=c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05, "~"=0.1), 
              test = 'wilcox.test', textsize = 10, y_position = 95) +
  ylab('Kg') +
  theme(axis.title.x=element_blank(), legend.position = "none", text = element_text(size = 19), 
        axis.title.y = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(2,1,2,1), 'lines'), axis.text = element_text(color = "black")) +
  scale_y_continuous(limits = c(40, 105), breaks = seq(40, 105, 10), guide = guide_prism_minor(), minor_breaks = seq(40, 100, 5)) +
  scale_fill_manual(values=c("#c0dcc0", "#cc6677")) + 
  ggtitle('pbw'); p
result.list <- append(result.list, list(p))


# height
r <- result[result$variable=='height',]
p <- ggplot(r, aes(x=endotype, y=value, fill=endotype)) + 
  geom_boxplot(width=.35, outlier.shape = NA) + 
  theme_classic() +
  geom_jitter(shape=16, size=1, position=position_jitter(width = 0.1, height = 0, seed = 2023), aes(alpha=0.9)) +
  geom_signif(comparisons = list(c("E8", "E9")), 
              map_signif_level=c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05, "~"=0.1), 
              test = 'wilcox.test', textsize = 10, y_position = 200) +
  ylab('Cm') +
  theme(axis.title.x=element_blank(), legend.position = "none", text = element_text(size = 19), 
        axis.title.y = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(2,1,2,1), 'lines'), axis.text = element_text(color = "black")) +
  scale_y_continuous(limits = c(150, 209), breaks = seq(150, 209, 10), guide = guide_prism_minor(), minor_breaks = seq(150, 200, 5)) +
  scale_fill_manual(values=c("#c0dcc0", "#cc6677")) + 
  ggtitle('height'); p
result.list <- append(result.list, list(p))


# packyr
r <- result[result$variable=='packyr',]
p <- ggplot(r, aes(x=endotype, y=value, fill=endotype)) + 
  geom_boxplot(width=.35, outlier.shape = NA) + 
  theme_classic() +
  geom_jitter(shape=16, size=1, position=position_jitter(width = 0.1, height = 0, seed = 2023), aes(alpha=0.9)) +
  geom_signif(comparisons = list(c("E8", "E9")), 
              map_signif_level=c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05, "*"=1), 
              test = 'wilcox.test', textsize = 10, y_position = 125) +
  ylab('Packages') +
  theme(axis.title.x=element_blank(), legend.position = "none", text = element_text(size = 19), 
        axis.title.y = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(2,1,2,1), 'lines'), axis.text = element_text(color = "black")) +
  scale_y_continuous(limits = c(10, 145), breaks = seq(10, 130, 20), guide = guide_prism_minor(), minor_breaks = seq(10, 120, 10)) +
  scale_fill_manual(values=c("#c0dcc0", "#cc6677")) + 
  ggtitle('packyr'); p
result.list <- append(result.list, list(p))


# temp
r <- result[result$variable=='temp',]
p <- ggplot(r, aes(x=endotype, y=value, fill=endotype)) + 
  geom_boxplot(width=.35, outlier.shape = NA) + 
  theme_classic() +
  geom_jitter(shape=16, size=1, position=position_jitter(width = 0.1, height = 0, seed = 2023), aes(alpha=0.9)) +
  geom_signif(comparisons = list(c("E8", "E9")), 
              map_signif_level=c("****"=0.0001, "***"=0.001, "**"=0.01, "~"=0.05, "~"=0.1), 
              test = 'wilcox.test', textsize = 10, y_position = 39.4) +
  ylab('Celsius') +
  theme(axis.title.x=element_blank(), legend.position = "none", text = element_text(size = 19), 
        axis.title.y = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(2,1,2,1), 'lines'), axis.text = element_text(color = "black")) +
  scale_y_continuous(limits = c(35, 40), breaks = seq(35, 40.5, 1), guide = guide_prism_minor(), minor_breaks = seq(35, 40, .5)) +
  scale_fill_manual(values=c("#c0dcc0", "#cc6677")) + 
  ggtitle('temp'); p
result.list <- append(result.list, list(p))


# pao2screen
r <- result[result$variable=='pao2screen',]
p <- ggplot(r, aes(x=endotype, y=value, fill=endotype)) + 
  geom_boxplot(width=.35, outlier.shape = NA) + 
  theme_classic() +
  geom_jitter(shape=16, size=1, position=position_jitter(width = 0.1, height = 0, seed = 2023), aes(alpha=0.9)) +
  geom_signif(comparisons = list(c("E8", "E9")), 
              map_signif_level=c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05, "~"=0.1), 
              test = 'wilcox.test', textsize = 10, y_position = 300) +
  ylab('mm Hg') +
  theme(axis.title.x=element_blank(), legend.position = "none", text = element_text(size = 19), 
        axis.title.y = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(2,1,2,1), 'lines'), axis.text = element_text(color = "black")) +
  scale_y_continuous(limits = c(45, 340), breaks = seq(45, 300, 50), guide = guide_prism_minor(), minor_breaks = seq(45, 300, 25)) +
  scale_fill_manual(values=c("#c0dcc0", "#cc6677")) + 
  ggtitle('pao2screen'); p
result.list <- append(result.list, list(p))


############################################################################################################
# pf7
############################################################################################################
# read clinical data
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 5, na = c('.','.mild'))

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# reorder to match protein data
clinical <- clinical[match(ptid$ptid, clinical$ptid),]

# assign cluster 
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)

# continuous variable
var.cont <- c('pf7')

# create visualization data
result <- data.frame()
for(v in var.cont){
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]);var.1
  var.2 <- as.numeric(clinical[clinical$cluster==2,v]);var.2
  endotype <- c(rep('E8', length(var.1)), rep('E9', length(var.2)))
  variable <- rep(v, length(var.1)+length(var.2))
  d <- data.frame(variable, endotype, value=c(var.1, var.2))
  result <- rbind(result, d)
}

# plot pf7
p <- ggplot(result, aes(x=endotype, y=value, fill=endotype)) + 
  geom_boxplot(width=.35, outlier.shape = NA) + 
  theme_classic() +
  geom_jitter(shape=16, size=1, position=position_jitter(width = 0.1, height = 0, seed = 2023), aes(alpha=0.9)) +
  geom_signif(comparisons = list(c("E8", "E9")), 
              map_signif_level=c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05, "~"=0.1), 
              test = 'wilcox.test', textsize = 10, y_position = 400) +
  ylab('mm Hg') +
  theme(axis.title.x=element_blank(), legend.position = "none", text = element_text(size = 19), 
        axis.title.y = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(2,1,2,1), 'lines'), axis.text = element_text(color = "black")) +
  scale_y_continuous(limits = c(90, 460), breaks = seq(90, 400, 50), guide = guide_prism_minor(), minor_breaks = seq(90, 400, 25)) +
  scale_fill_manual(values=c("#c0dcc0", "#cc6677")) + 
  ggtitle('pf7'); p
result.list <- append(result.list, list(p))


###############################################################################################
# renal28 coag28 vfd
###############################################################################################
clinical <- read_excel('Metadata of Cases_Ji summarized - copy.xlsx', sheet = 13)
clinical[clinical=='.'] <- NA

# remove columns with all NAs
clinical <- clinical[,colSums(is.na(clinical)) < nrow(clinical)]

# reorder to match protein data
clinical <- clinical[match(ptid$ptid, clinical$ptid),]

# assign cluster 
clinical <- merge(clinical, ptid, by = 'ptid')
table(clinical$cluster)

# relevant continuous variable
var.cont <- c('renal28', 'coag28', 'vfd')
result <- data.frame()

# create visualization data
for(v in var.cont){
  variable <- rep(v, dim(clinical)[1])
  var.1 <- as.numeric(clinical[clinical$cluster==1,v]);var.1
  var.2 <- as.numeric(clinical[clinical$cluster==2,v]);var.2
  endotype <- c(rep('E8', length(var.1)), rep('E9', length(var.2)))
  d <- data.frame(variable, endotype, value=c(var.1, var.2))
  result <- rbind(result, d)
}

# plot
for(i in 1:length(var.cont)){
  v <- var.cont[i]
  r <- result[result$variable==v,]
  p <- ggplot(r, aes(x=endotype, y=value, fill=endotype)) + 
    geom_boxplot(width=.35, outlier.shape = NA) + 
    theme_classic() +
    geom_jitter(shape=16, size=1, position=position_jitter(width = 0.1, height = 0, seed = 2023), aes(alpha=0.9)) +
    geom_signif(comparisons = list(c("E8", "E9")), 
                map_signif_level=c("****"=0.0001, "***"=0.001, "**"=0.01, "*"=0.05, "~"=0.1), 
                test = 'wilcox.test', textsize = 10, y_position = 30) +
    ylab('Day') +
    theme(axis.title.x=element_blank(), legend.position = "none", text = element_text(size = 19), 
          axis.title.y = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          plot.title = element_text(hjust = 0.5, size = 18),
          plot.margin = unit(c(2,1,2,1), 'lines'), axis.text = element_text(color = "black")) +
    scale_y_continuous(limits = c(0, 35), breaks = seq(0, 30, 10), guide = guide_prism_minor(), minor_breaks = seq(0, 30, 5)) +
    scale_fill_manual(values=c("#c0dcc0", "#cc6677")) + 
    ggtitle(v); p
  result.list <- append(result.list, list(p))
}
############################################################################################################
# barplot race gender
############################################################################################################
variable <- c(rep('Race', 4), rep('Gender', 4))
outcome <- c('Black','Black','White','White', 'Female','Female','Male','Male')
endotype <- rep(c('E8','E9'), 4)
value <- c(41, 8.5, 59, 91.5, 65, 29, 35, 71)
result.cat <- data.frame(variable,outcome,endotype,value)

# Race
r <- result.cat[result.cat$variable=='Race',]
p <- ggplot(r, aes(y=value, x=endotype)) + 
  geom_bar(aes(fill=outcome), position="fill", stat="identity", width=.35) +
  scale_y_continuous(limits = c(0, 1.25), labels=scales::percent_format(suffix = ""), 
                     oob = rescale_none, breaks = seq(0, 1, by=0.25), expand = c(0, 0)) +
  theme_classic() +
  geom_signif(aes(xmin = "E8", xmax = "E9", annotations = '*', y_position = 1.05),
              textsize = 10, tip_length = .0003) + 
  scale_fill_manual(values=c("#cc0000", "#56B4E9")) + 
  ggtitle('Race (%)') + 
  theme(axis.title.x=element_blank(), text = element_text(size = 19), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(2,0,2,1), 'lines'), axis.text = element_text(color = "black"),
        legend.title=element_blank(),legend.position="right",
        legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=10),
        legend.margin = margin(t = 0, r = 10, b = 0, l = -5)); p

result.list <- append(result.list, list(p))


# Gender
r <- result.cat[result.cat$variable=='Gender',]
p <- ggplot(r, aes(y=value, x=endotype)) + 
  geom_bar(aes(fill=outcome), position="fill", stat="identity", width=.35) +
  scale_y_continuous(limits = c(0, 1.25), labels=scales::percent_format(suffix = ""), 
                     oob = rescale_none, breaks = seq(0, 1, by=0.25), expand = c(0, 0)) +
  theme_classic() +
  geom_signif(aes(xmin = "E8", xmax = "E9", annotations = '*', y_position = 1.05),
              textsize = 10, tip_length = .0003) + 
  scale_fill_manual(values=c("#cc0000", "#56B4E9")) + 
  ggtitle('Gender (%)') + 
  theme(axis.title.x=element_blank(), text = element_text(size = 19), 
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = unit(c(2,0,2,1), 'lines'), axis.text = element_text(color = "black"),
        legend.title=element_blank(),legend.position="right",
        legend.key.size = unit(0.4, 'cm'), legend.text = element_text(size=10),
        legend.margin = margin(t = 0, r = 10, b = 0, l = -5)); p

result.list <- append(result.list, list(p))

# plot all in pannels
#svg(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/relevance_outcome.svg", width = 11, height = 10)
ggarrange(plotlist=result.list)
#dev.off()

############################################################################################################
############################################################################################################
# heat map of DE genes
############################################################################################################
############################################################################################################
data.up <- readRDS('data.up.rds'); dim(data.up)
data.down <- readRDS('data.down.rds'); dim(data.down)

# standardization 
data.DE <- rbind(data.up, data.down)[,-c(1:6)]; dim(data.DE)
data.DE.znorm <- t(scale(t(data.DE))); dim(data.DE.znorm)
mean(data.DE.znorm[1,])
sd(data.DE.znorm[1,])

top = HeatmapAnnotation(
  type=anno_block(gp = gpar(fill = c('green', 'deeppink')), labels = c('Control', 'ARDS'), labels_gp = gpar(fontface = "bold"))
)

left = rowAnnotation(
  type=anno_block(gp = gpar(fill = c('darkorange', 'purple')), labels = c('Up-Regulated', 'Down-Regulated'), labels_gp = gpar(fontface = "bold"))
)

split.column <- c(rep(1, 49), rep(2, 57))
split.row <- c(rep(1, 426), rep(2, 435))

pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/DE.pdf", height = 8)
ht.DE <- Heatmap(data.DE.znorm, name = "Abundance \n (Z-Score)", km = 1, top_annotation = top, left_annotation = left,
        show_row_names = FALSE, show_column_names = FALSE, 
        show_column_dend = F, show_row_dend = F, cluster_columns = F, cluster_rows = F, 
        column_split = split.column, column_title = NULL,
        row_split = split.row, row_title = NULL, 
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 12, fontface = "bold"), title_gp = gpar(fontsize = 12, fontface = "bold"))); ht.DE
dev.off()

############################################################################################################
# volcano plot
############################################################################################################
DE.summary <- readRDS('DE.summary.rds')

pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/vocano.pdf", width = 8)
EnhancedVolcano(DE.summary, lab = NA, x = 'log2_fold_change', y = 'adjusted_pvalue',
                FCcutoff = 1, pCutoff = 0.05, pointSize = 2, legendPosition = 'none',
                xlim = c(-6,6), ylim = c(0,16.5), subtitle = NULL, title = NULL,
                caption = NULL,  cutoffLineWidth = 0.8, gridlines.major = F,
                gridlines.minor = F) 
  #annotate("text", x=5, y=2, label='FDR = 0.05', size=5, fontface = "bold") +
  #annotate("text", x=-4, y=16.5, label='Down-Regulated (2-Fold)', size=5, fontface = "bold") +
  #annotate("text", x=4, y=16.5, label='Up-Regulated (2-Fold)', size=5, fontface = "bold") 
dev.off()  

############################################################################################################
# Correlation
############################################################################################################
# continuous
result.cont <- read.csv('cor.cont.imp.csv')
sort(table(result.cont$var), decreasing = T)

# create 3 types of cor
result.cont$sig.type <- 'Both'
result.cont[abs(result.cont$r.pearson)>0.4 & abs(result.cont$r.spearman)<0.4,]$sig.type <- 'Pearson'
result.cont[abs(result.cont$r.pearson)<0.4 & abs(result.cont$r.spearman)>0.4,]$sig.type <- 'Spearman'
table(result.cont$sig.type)

# assign cor for unique
result.cont$cor <- 0
result.cont$cor[result.cont$sig.type=='Pearson'] <- result.cont[result.cont$sig.type=='Pearson', 'r.pearson']
result.cont$cor[result.cont$sig.type=='Spearman'] <- result.cont[result.cont$sig.type=='Spearman', 'r.spearman']

# assign cor for double
result.cont$cor[result.cont$sig.type=='Both' & result.cont$r.pearson>0] <- 
  apply(result.cont[result.cont$sig.type=='Both' & result.cont$r.pearson>0, c('r.pearson', 'r.spearman')], 1, max)

result.cont$cor[result.cont$sig.type=='Both' & result.cont$r.pearson<0] <- 
  apply(result.cont[result.cont$sig.type=='Both' & result.cont$r.pearson<0, c('r.pearson', 'r.spearman')], 1, min)

# add sign
result.cont$sign <- 'Positive'
result.cont[result.cont$cor<0,]$sign <- 'Negative'

# creat data frame with all combinations
d <- expand.grid(unique(result.cont$var), unique(result.cont$gene.keep))
colnames(d) <- c('var', 'gene.keep')
result.cont <- merge(d, result.cont, by = c('var', 'gene.keep'), all = T)

# change legend order
result.cont$sign <- factor(result.cont$sign, levels = c("Positive", "Negative"))
result.cont$sig.type <- factor(result.cont$sig.type, levels = c("Pearson", "Spearman", "Both"))

# order variable by significant protein
result.cont$var <- factor(result.cont$var, levels=names(sort(table(result.cont[complete.cases(result.cont),]$var), decreasing = T)))

# visualize
svg(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/cor_cont.svg", height = 10, width = 6)
ggplot(result.cont, aes(x = var, y = fct_reorder(gene.keep, cor))) +
  geom_tile(color= "grey", fill= "white", size = 0.1) +
  geom_point(aes(shape=sig.type, size=abs(cor), color = sign), alpha=0.8) + 
  theme(axis.ticks = element_blank(), legend.text=element_text(size=11), legend.title=element_text(size=11)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(shape="Significance", color="Correlation Sign", size='Absolute Correlation') + 
  scale_size_continuous(breaks=rep(c(0.4, 0.45, 0.5, 0.55),3), labels = c(rep("",8),'0.40','0.45','0.50','0.55'),
                        limits = c(0.4, 0.55), range = c(1,4.5)) +
  scale_color_manual(values = c("Negative" = "#1B9E77", "Positive"="#D95F02")) +
  guides(shape = guide_legend(override.aes = list(size = 4), order=1), 
         color = guide_legend(override.aes = list(size = 4), order=2),
         size  = guide_legend(order = 3)) +
  theme(legend.key=element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.background = element_blank()) +
  guides(size=guide_legend(ncol=3, override.aes = list(shape=c(rep(16,4),rep(17,4),rep(15,4))))) +
  scale_shape_manual(values=c(16,17,15), na.translate = FALSE)
dev.off()
    
summary(result.cont$r.pearson)
summary(abs(result.cont$r.pearson))
summary(abs(result.cont$r.spearman))

# categorical
result.cat <- read.csv('cor.cat.imp.csv')

# add sign
result.cat$sign <- 'Positive'
result.cat[result.cat$fc<0,]$sign <- 'Negative'

# Capitalize 
result.cat$var <- str_to_title(result.cat$var)

# change legend order
result.cat$sign <- factor(result.cat$sign, levels = c("Positive", "Negative"))

# order variable by significant protein
result.cat$var <- factor(result.cat$var, levels=names(sort(table(result.cat[complete.cases(result.cat),]$var), decreasing = T)))

# creat data frame with all combinations
d <- expand.grid(unique(result.cat$var), unique(result.cat$gene.keep))
colnames(d) <- c('var', 'gene.keep')
result.cat <- merge(d, result.cat, by = c('var', 'gene.keep'), all = T)

# visualize
svg(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/cor_cat.svg", height = 11, width = 5)
ggplot(result.cat, aes(x = var, y = fct_reorder(gene.keep, fc))) +
  geom_tile(color= "grey", fill= "white", size = 0.1) +
  geom_point(aes(size=abs(fc), color = sign), alpha=0.8) +
  theme(axis.ticks = element_blank(), legend.text=element_text(size=11), legend.title=element_text(size=11)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(color="Fold Change Sign", size='Absolute Fold Change') +
  scale_size_continuous(breaks=c(2,4,7,10),labels = c('2.0','4.0','7.0','10.0'),
                        range = c(1,4), limits = c(2, 10)) + 
  scale_color_manual(values = c("Negative" = "#1B9E77", "Positive"="#D95F02")) +
  guides(color = guide_legend(override.aes = list(size = 3.5), order=1),
         size = guide_legend(order = 2)) +
  theme(legend.key=element_blank(), axis.text.y = element_text(size = 7),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.background = element_blank()) 
dev.off()

summary(abs(result.cat$fc))


############################################################################################################
# pca plot on DEP
############################################################################################################
data.DE <- readRDS('data.DE.rds'); dim(data.DE)
data.DE <- data.DE[,-c(1:6)]; dim(data.DE)

pr.out <- prcomp(t(data.DE), center = T, scale. = T)
data.pca <- pr.out$x; dim(data.pca)
PVE.matrix <- summary(pr.out)$importance; PVE.matrix
PVE <- PVE.matrix[2,]
plot(PVE, xlab='Principle Components', ylab='Proportion of Variance Explained', cex.lab=1.5, pch = 19)

color <- c(rep('green',49), rep('red',57))

# pc1 vs pc2
pdf(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/PCA.pdf", width = 10)
plot(data.pca[,2]~data.pca[,1], col = color, pch = 19, cex.lab=1.3, cex.axis=1.3, cex.sub=1.3, cex=1.5,
     xlab='PC1 (38.7%)', ylab='PC2 (6.0%)', xlim=c(-30,30))
legend("topleft", legend = c("Control", "ARDS"), col=c('green','red'), pch=16, cex=1.3, pt.cex = 1)
dev.off()

# pc2 vs pc3
plot(data.pca[,2]~data.pca[,3], col = color, pch = 19, main='PC2 vs. PC3 (All)',
     xlab='PC2 (6.0%)', ylab='PC3 (3.5%)')
legend("topleft", legend = c("Control", "ARDS"), col=c('green','red'), pch=16)


############################################################################################################
############################################################################################################
# batch effect examination
############################################################################################################
############################################################################################################
# order samples by batch
order <- read_excel("order.xlsx", col_names = T)
sample <- colnames(data)
sample.short <- sapply(str_split(sample, " - "),tail,1)
time <- match(sample.short, order$sample)

############################################################################################################
# mean abundance before correction
############################################################################################################
# read data 
data <- readRDS('data.norm.impute.80.rds'); dim(data)
data <- t(data); dim(data)

# mean expression per sample
sample.mean <- colMeans(data, na.rm = T)
summary(sample.mean)

# plot mean abundance across sample time (median result is very similar)
svg(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/Batch_mean_before.svg", height = 5, width = 11)
plot(sample.mean[order(time)], pch = 19, ylab = 'Mean Abundance', xlab = 'Order', 
     main = 'Before Batch Effect Correction',cex.lab=1.5, cex.main=1.5, cex.axis=1.3,
     ylim=c(25.36, 25.57)); abline(v=c(3.5,12.5,72.5), col=2, lty=2, lwd=2)
dev.off()

##############################################################
# sample distribution boxplot before correction
##############################################################
# boxplot for all
svg(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/Batch_boxplot_before.svg", height = 5, width = 11)
color <- c(rep('red',3), rep('green',9), rep('blue',60), rep('orange',34))
boxplot(data[,order(time)], outline = F, col = color, whisklty = 1, xaxt = "n", ylim=c(17, 34),
        cex.lab=1.5, cex.main=1.5, cex.axis=1.3,
        main="Before Batch Effect Correction", ylab="Abundance", xlab="Order")
dev.off()

##############################################################
# correlation before correction
##############################################################
# correlations within batch 
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
  data.batch <- data[,which(sample.short%in%batch)]
  cor.batch <- unique(as.vector(cor(data.batch, use = 'complete.obs')))
  cor.batch <- cor.batch[cor.batch!=1]
  cor.within <- c(cor.within, cor.batch)
}

# correlations between batches
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
  cor.between <- c(cor.between, as.numeric(cor(data[,index.1], data[,index.2], use = 'complete.obs')))
}

svg(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/Batch_cor_before.svg", width = 6, height = 6)
boxplot(cor.within, cor.between, names=c('Within Batch','Between Batches'), ylab='Pearson Correlation', 
        main='Before Batch Effect Correction', cex.lab=1.5, cex.main=1.5, cex.axis=1.3, ylim=c(0.5,1))
dev.off()

##############################################################
# PCA Dimension Reduction before correction
##############################################################
pr.out <- prcomp(t(data), center = T, scale. = T)
data.pca <- pr.out$x; dim(data.pca)
PVE.matrix <- summary(pr.out)$importance; PVE.matrix
PVE <- PVE.matrix[2,]
plot(PVE, xlab='Principle Components', ylab='Proportion of Variance Explained', cex.lab=1.5, pch = 19)

# plot first two PCs
svg(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/Batch_PCA_before.svg", width = 6, height = 6)
color <- c(rep('red',3), rep('green',9), rep('blue',60), rep('orange',34))
plot(data.pca[order(time),2]~data.pca[order(time),1], col = color, pch = 19, main='Before Batech Effect Correction',
     xlab='PC1 (18.5%)', ylab='PC2 (10.7%)', cex.lab=1.5, cex.main=1.5, cex.axis=1.3, ylim=c(-21,50), xlim=c(-31, 34))
dev.off()

############################################################################################################
# mean abundance after correction
############################################################################################################
# read data after correction
data <- readRDS('data.complete.80.rds'); dim(data)

# mean expression per sample
sample.mean <- colMeans(data, na.rm = T)
summary(sample.mean)

# plot mean abundance across sample time (median result is very similar)
svg(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/Batch_mean_after.svg", height = 5, width = 11)
plot(sample.mean[order(time)], pch = 19, ylab = 'Mean Abundance', xlab = 'Order', 
     main = 'After Batch Effect Correction',cex.lab=1.5, cex.main=1.5, cex.axis=1.3,
     ylim=c(25.36, 25.57)); abline(v=c(3.5,12.5,72.5), col=2, lty=2, lwd=2)
dev.off()

##############################################################
# sample distribution boxplot after correction
##############################################################
# boxplot for all
svg(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/Batch_boxplot_after.svg", height = 5, width = 11)
color <- c(rep('red',3), rep('green',9), rep('blue',60), rep('orange',34))
boxplot(data[,order(time)], outline = F, col = color, whisklty = 1, xaxt = "n", ylim=c(17, 34),
        cex.lab=1.5, cex.main=1.5, cex.axis=1.3,
        main="After Batch Effect Correction", ylab="Abundance", xlab="Order")
dev.off()

##############################################################
# correlation after correction
##############################################################
# correlations within batch 
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
  data.batch <- data[,which(sample.short%in%batch)]
  cor.batch <- unique(as.vector(cor(data.batch, use = 'complete.obs')))
  cor.batch <- cor.batch[cor.batch!=1]
  cor.within <- c(cor.within, cor.batch)
}

# correlations between batches
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
  cor.between <- c(cor.between, as.numeric(cor(data[,index.1], data[,index.2], use = 'complete.obs')))
}

svg(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/Batch_cor_after.svg", width = 6, height = 6)
boxplot(cor.within, cor.between, names=c('Within Batch','Between Batches'), ylab='Pearson Correlation', 
        main='After Batch Effect Correction', cex.lab=1.5, cex.main=1.5, cex.axis=1.3, ylim=c(0.5, 1))
dev.off()

##############################################################
# PCA Dimension Reduction after correction
##############################################################
pr.out <- prcomp(t(data), center = T, scale. = T)
data.pca <- pr.out$x; dim(data.pca)
PVE.matrix <- summary(pr.out)$importance; PVE.matrix
PVE <- PVE.matrix[2,]
plot(PVE, xlab='Principle Components', ylab='Proportion of Variance Explained', cex.lab=1.5, pch = 19)

# plot first two PCs
svg(file="C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/figures/Batch_PCA_after.svg", width = 6, height = 6)
color <- c(rep('red',3), rep('green',9), rep('blue',60), rep('orange',34))
plot(data.pca[order(time),2]~data.pca[order(time),1], col = color, pch = 19, main='After Batech Effect Correction',
     xlab='PC1 (20.7%)', ylab='PC2 (7.8%)', cex.lab=1.5, cex.main=1.5, cex.axis=1.3, ylim=c(-21,50), xlim=c(-31, 34))
dev.off()


