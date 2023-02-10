setwd("~/Documents/meso_stats/")

library("ape")
library(ggplot2)
library(ggtree)
library(readxl)
library(dplyr)


#Load tree
tree <- read.tree(file = "RAxML_bipartitions.protein_fas2.accurate_rename.tre")

#Root tree
tree_r <- root(tree, "S_meliloti")

#read tip names to 
h = read.csv("keep.csv", header = T)

#prune
pruned <- keep.tip(tree_r, h$Genome)

#base ggtree 
z<- ggtree(pruned, ladderize = TRUE, ) +  geom_rootedge(rootedge = 0.1)+
  geom_tiplab(size = 1.2, align = TRUE, linesize = 0.107) + 
  geom_treescale(width=0.1, x=0, y=288, offset=1) 

#rotate it how I want 
z1 <- rotate(z, 303) %>% rotate(450) %>% rotate(304)

#read in mastersheet
D1 = read.csv("Delgado_Meso_Mastersheet_by_Strain_2023-01-30.csv", header = T)


####For tree + MIC

mic <- D1 %>%
  select(Unique.ID,   Ni.MIC_mean
  )
mic2 <- mic
mic$Unique.ID <- NULL
mic$Ni.MIC_mean <- as.numeric(mic$Ni.MIC_mean)
row.names(mic) <- mic2$Unique.ID

z2 <- gheatmap(z1, mic, offset=.04, width=.1,
               colnames = FALSE)+ 
  scale_colour_gradient( low = "darkolivegreen1",
                         high = "darkolivegreen",name = "MIC",
                         aesthetics = "fill" , na.value = "white")
z2




####For tree + OD

od <- D1 %>%
  select(Unique.ID,   OD72hr.ni_mean
  )
od2 <- od
od$Unique.ID <- NULL
od$OD72hr.ni_mean <- as.numeric(od$OD72hr.ni_mean)
row.names(od) <- od2$Unique.ID

z2 <- gheatmap(z1, od, offset=.04, width=.1,
               colnames = FALSE) + 
  scale_colour_gradient( low = "lightblue1",
                         high = "blue3",name = "OD",
                         aesthetics = "fill" , na.value = "white")
z2



####For tree plus Ni ppm

ni <- D1 %>%
  select(Unique.ID,   Ni.ppm
  )
ni2 <- ni
ni$Unique.ID <- NULL
ni$Ni.ppm <- as.numeric(ni$Ni.ppm)
row.names(ni) <- ni2$Unique.ID

z2 <- gheatmap(z1, ni, offset=.04, width=.1,
               colnames = FALSE) + 
  scale_colour_gradient( low = "darkseagreen1",
                         high = "darkgreen",name = "Soil Ni",
                         aesthetics = "fill" , na.value = "white")
z2




####For tree plus Soil Type

palette <-  c("navajowhite","forestgreen")

soil <- D1 %>%
  select(Unique.ID,   SerpentineType
  )
soil2 <- soil
soil$Unique.ID <- NULL
row.names(soil) <- soil2$Unique.ID

z2 <- gheatmap(z1, soil, offset=.04, width=.1,
               colnames = FALSE) + 
  scale_fill_manual(values = palette, na.value = "white", name = "Soil Type") 

z2
