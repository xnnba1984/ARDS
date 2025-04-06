library(readxl)
library(stringr)
library(ggplot2)
library(missForest)
library(vsn)
library(sva)
library(NbClust) 
library(ComplexHeatmap)

setwd("C:/Users/mxi1/OneDrive - Loyola University Chicago/Dr Ji/Dataset 1-from Marco for ARDS")

# read clean data
data.DE <- readRDS('data.DE.rds'); dim(data.DE)

# read significant proteins
cor.cont.imp <- read.csv('cor.cont.imp.csv')
cor.cat.imp <- read.csv('cor.cat.imp.csv')

# death90 is outcome
cor.cat.imp[cor.cat.imp$var=='death90', 'type'] <- 'outcome'

# data based on importain proteins
#data.cor.imp <- data.DE[data.DE$gene.keep%in%cor.cont.imp$gene.keep | data.DE$gene.keep%in%cor.cat.imp$gene.keep,]
gene.imp <- unique(cor.cont.imp$gene.keep[cor.cont.imp$type=='severity'])
gene.imp <- unique(c(cor.cont.imp$gene.keep[cor.cont.imp$type=='outcome'], cor.cat.imp$gene.keep[cor.cat.imp$type=='outcome']))
gene.imp <- unique(cor.cat.imp$gene.keep[cor.cat.imp$type=='Etiology'])
gene.imp <- unique(c(cor.cont.imp$gene.keep[cor.cont.imp$type=='baseline'], cor.cat.imp$gene.keep[cor.cat.imp$type=='baseline']))

# read death biomarker
biomarker.death <- read.csv('biomarker_death_combine_all.csv')
gene.death <- unlist(strsplit(biomarker.death$protein[1], ","))

# read shock biomarker
biomarker.shock <- read.csv('biomarker_shock_combine_all.csv')
gene.shock <- unlist(strsplit(biomarker.shock$protein[2], ","))

# combine
gene.imp <- unique(c(gene.death, gene.shock))

# read Aspir biomarker
biomarker.Aspir <- read.csv('biomarker_Aspir_combine_all.csv')
gene.Aspir <- unlist(strsplit(biomarker.Aspir$protein[1], ","))

# read Pneumo biomarker
biomarker.Pneumo <- read.csv('biomarker_Pneumo_combine_all.csv')
gene.Pneumo <- unlist(strsplit(biomarker.Pneumo$protein[1], ","))

# read Sepsis biomarker
biomarker.Sepsis <- read.csv('biomarker_Sepsis_combine_all.csv')
gene.Sepsis <- unlist(strsplit(biomarker.Sepsis$protein[1], ","))

# combine
gene.imp <- unique(c(gene.Aspir, gene.Pneumo, gene.Sepsis))

data.cor.imp <- data.DE[data.DE$gene.keep%in%gene.imp,]
data.cor.imp.ARDS <- data.cor.imp[,56:112]; dim(data.cor.imp.ARDS)

# standardize
data.cor.imp.ARDS.norm <- t(scale(t(data.cor.imp.ARDS))); dim(data.cor.imp.ARDS.norm)
sum(data.cor.imp.ARDS.norm[1,]); sd(data.cor.imp.ARDS.norm[1,])

########################################################################################################
# find number of clusters k
########################################################################################################
index.all <- c("kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew", "friedman", "rubin", "cindex", "db", "silhouette", 
            "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "gamma", "gplus", "tau", "dunn", 
            "hubert", "sdindex", "dindex", "sdbw")

#index.all <- c('ch','duda','pseudot2','beale','gap','mcclain', 'cindex', 'dunn','marriot','trcovw','tracew','rubin','ball','ratkowsky')
val.list <- list()
set.seed(2023)
system.time({
  for(i in index.all){
    tryCatch({
      #i <- index[2]
      #print('====================='); print(i); print('=====================')
      nb <- NbClust(t(data.cor.imp.ARDS.norm), distance = "euclidean", min.nc = 2, max.nc = 4, method = "kmeans", index = i)
      #print(nb$All.index)
      #print(table(nb$Best.partition))
      l <- length(table(nb$Best.partition))
      if(l==2){
        print('====================='); print(i); print('=====================')
        print(l)
      }
      val.list[`i`] <- list(nb)
    }, error=function(e){print(e)})
  }
})
#saveRDS(val.list, 'val.kmeans.imp.rds')

val.result <- unlist(lapply(val.list, function(x){
  length(table(x$Best.partition))
})); table(val.result)

########################################################################################################
# clustering and heatmap
########################################################################################################
set.seed(2023)
hm <- Heatmap(data.cor.imp.ARDS.norm, show_column_names = FALSE, cluster_rows = T, cluster_columns = T, show_row_names = FALSE,
               use_raster=F, show_column_dend = T, show_row_dend = F, column_km = 2, column_km_repeats = 100,
               clustering_distance_columns = "euclidean", name = "Abundance \n (Z-Score)"); hm <- draw(hm)

cluster <- lapply(column_order(hm), sort); cluster
unlist(lapply(column_order(hm), length))

# save clustering result
#saveRDS(cluster, 'kmeans.imp.outcome.rds')
#saveRDS(cluster, 'kmeans.imp.baseline.2.rds')
#saveRDS(cluster, 'kmeans.biomarker.outcome.rds')
#saveRDS(cluster, 'kmeans.biomarker.death.shock.rds')
saveRDS(cluster, 'kmeans.biomarker.Aspir.Pneumo.Sepsis.rds')



