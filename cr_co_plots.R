

setwd("~/Documents/meso_stats/")
library(readxl)
library(ggplot2)


D1 = read.csv("Delgado_Meso_Mastersheet_by_Strain_2023-02-17.csv", header = T)




ggplot(D1, aes(x=Ni.MIC_mean, y=Co.MIC_mean)) + 
  geom_jitter(width = 0.1, height = 0.1, color = "blue3")+
  theme(panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text= element_text(size=20),
        axis.title = element_text(size=20)) +
  labs(x = "Nickel MIC (mM)", 
       y = "Cobalt MIC (mM)")


ggplot(D1, aes(x=Ni.MIC_mean, y=Cr.MIC_mean)) +
  geom_jitter(width = 0.1, height = 0.1)  +
  theme(panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text= element_text(size=20),
        axis.title = element_text(size=20)) +
  labs(x = "Nickel MIC (mM)", 
       y = "Chromium MIC (mM)")




