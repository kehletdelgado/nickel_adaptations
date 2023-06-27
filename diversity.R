setwd("~/Documents/meso_stats/")

library(vegan)
library("ape")
library(ggplot2)
library(ggtree)
library(readxl)
library(readr)
library(data.table)
library(tidyverse)
library(tidygenomes)
library("dplyr")


#read in tree
tree <- read.tree(file = "pruned_302_phylophlan_faa.tre")

#read in metadata
df = read.csv("Delgado_Meso_Mastersheet_by_Strain_2023-02-17.csv", header = T)
#read in filter list
h = read.csv("keep.csv", header = T)
#filter list to strains we are keeping
filtered_df <- filter(df, Unique.ID %in% h$Genome)
dim(filtered_df)

###Read distance matrix from EMBOSS
matty <- read_phylip_distmat("fas_acc.trim.distmat", skip = 7, include_diagonal = F)
matty2 <- pairs2matrix(matty, object_1 = sequence_1, object_2 = sequence_2,
             measure = distance, diag = TRUE)
matt3 <- as.data.frame(matty2)
dim(matt3)
matt3 <- data.frame(names = row.names(matt3), matt3)
names(matt3)[names(matt3) == 'names'] <- 'Genome'

#join meta data and distance 
meta_distance <- inner_join(filtered_df, matt3,  by= c("Unique.ID"='Genome'), 
                            keep = FALSE,
                            na_matches = c( "never"))
dim(meta_distance)

dist <- meta_distance %>% select(all_of(filtered_df$Unique.ID)) %>% as.dist()

set.seed(1991)

#PERMANOVA with adonis2 
adonis2(dist ~  Reserve * PlantSpecies.Strain * SerpentineType  , data = meta_distance, na.action=na.exclude)


#Test for differences in dispersion
#For reserve
disp <- betadisper(dist, meta_distance$Reserve)
anova(disp)
#For serpentine type
disp <- betadisper(dist, meta_distance$SerpentineType) 
anova(disp) 
#for host
disp <- betadisper(dist, meta_distance$PlantSpecies.Strain) 
anova(disp)

