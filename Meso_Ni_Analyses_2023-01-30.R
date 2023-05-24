#### Cami's analysis of soil chemistry and Meso traits for Hanna Delgado #####

#### LOAD LIBRARIES ###########################

library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
library(lme4)
library(DHARMa)
library(emmeans)
library(multcompView)
library(tidyr)
library(stringr)
library(plyr)


#### CAMI'S CODE ######################################################

# Note for Cami's code: run entire "Import Data" section before running any code in "Figures" or "Analyses" sections

# function to remove NA's from select columns
complete_fun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

CI <- function(x) {mean(x)+sd(x)/sqrt(length(x))}
CI2 <- function(x) {mean(x)-sd(x)/sqrt(length(x))}



#function for the standard error of a binary variable
binomial.SE = function(x) {sqrt((mean(x)*(1-mean(x)))/length(x))}

#function for logistic curve
logistic = function(x) {exp(x)/(1+exp(x))}


# IMPORT DATA #########################################################

# *** D1: PCA of soil data ################################################

# what soil variables discriminate most among field sites? 

setwd("C:/Users/camil/Washington State University (email.wsu.edu)/Plant Microbe Lab - Camille Wendlandt/Meso Nickel Tolerance Project/MASTERSHEET/")
D1 = read.csv("Meso_Mastersheet_by_Soil_2023-01-30.csv", header = T)


# keep just soil chemistry columns and identifying/grouping info
D1 <- select(D1, c("Reserve", "Site", "Subsite",
                   "SAMPLEID", "Ni.ppm", "OM", "ENR", "CEC",
                   "NO3_N", "P1", "HCO3_P", "PH", "BUFFER_PH", "S", "H",
                   "H_PCT", "K", "K_PCT", "MG", "MG_PCT", "CA", "CA_PCT",
                   "NA.", "NA_PCT"))

# further exclude redundant columns
D1 <- select(D1, -c("ENR", "BUFFER_PH", "H", "H_PCT", "K_PCT", "MG_PCT",
                    "CA_PCT", "NA_PCT"))

# choose one phosphorus measurement for each sample, based on soil pH
# this site recommends HCO3_P method if pH is 7.4 or greater; otherwise, use P1
# https://blog-crop-news.extension.umn.edu/2021/02/what-is-best-soil-test-option-for.html
D1$PH %>% unique %>% sort # soil pH ranges 5.6 to 8.0
D1$Phos <- ifelse(D1$PH<7.4, D1$P1, 
                  ifelse(D1$PH>=7.4, D1$HCO3_P, NA))

# exclude redundant phosphorus columns
D1 <- select(D1, -c("P1", "HCO3_P"))
names(D1) # 11 soil chemistry measurements
# [1] "Reserve"      "Site"         "Subsite.Soil" "SAMPLEID"     "Ni.ppm"       "OM"          
# [7] "CEC"          "NO3_N"        "PH"           "S"            "K"            "MG"          
# [13] "CA"           "NA."          "Phos"    


# how many reserves & sites contributed soil samples to the study?
D1 %>% nrow # 114
D1$Reserve %>% unique %>% length # 18 reserves
D1$Site %>% unique %>% length # 55 sites

# check 11 soil variables for outliers before doing PCA
# dotchart(D1$Ni) # no outliers
# dotchart(D1$OM) # one high value > 14; poss outlier
# dotchart(D1$CEC) # no outliers
# dotchart(D1$NO3_N) # no outliers
# dotchart(D1$PH) # no outliers
# dotchart(D1$S) # two high values > 25; poss outliers
# dotchart(D1$K) # no outliers
# dotchart(D1$MG) # no outliers
# dotchart(D1$CA) # no outliers
# dotchart(D1$NA.) # some high values, one > 80; poss outliers
# dotchart(D1$Phos) # no outliers

# look at potential outliers
D1[which(D1$OM>14),]$SAMPLEID # C0020
D1[which(D1$S>25),]$SAMPLEID # A0016 C0003
D1[which(D1$NA.>80),]$SAMPLEID # A0038

# for now, retain all samples
D1 %>% dim # 114 15

# for PCA, also drop any soils (rows) having NA for any variable
D1 <- na.omit(D1)
D1 %>% dim # 114 15
D1 %>% names
rownames(D1) <- D1$SAMPLEID

# do PCA, using just soil chemistry data columns
pca.soil = princomp(D1[,5:15], cor = TRUE)
summary(pca.soil)

# screeplot
screeplot(pca.soil, type = "l", npcs = 11, main = "Screeplot of the PCs")
# first 2 PCs most important, followed by next 2

# plot first 2 PCs
plot(pca.soil$scores[,1],pca.soil$scores[,2], xlab="PC1", ylab = "PC2")

# look at biplot
factoextra::fviz_pca(pca.soil, label = "var", repel = TRUE)

# save PC1 and PC2 scores
pcscores <- pca.soil$scores
pcscores <- as.data.frame(pcscores)
pcscores$SAMPLEID <- rownames(pcscores)
rownames(pcscores) <- NULL
pcscores <- select(pcscores, c("SAMPLEID", "Comp.1", "Comp.2"))
pcscores <- dplyr::rename(pcscores, PC1 = Comp.1, PC2 = Comp.2)


# look at loadings
loadings(pca.soil)
l = loadings(pca.soil)
l <- as.matrix(l[,])
l <- as.data.frame(l)
l$Variable <- rownames(l)
rownames(l) <- NULL
head(l)


# PC1 (OM, K)
l$Abs.Comp.1 <- abs(l$Comp.1) # absolute value of loadings
l$Abs.Comp.1 %>% sort(decreasing = TRUE) -> top # save in order from highest down
l[which(l$Abs.Comp.1==top[1]),c("Variable", "Comp.1")] # top loading = OM
l[which(l$Abs.Comp.1==top[2]),c("Variable", "Comp.1")] # 2nd top loading = K

# PC2 (Ca, Ni)
l$Abs.Comp.2 <- abs(l$Comp.2) # absolute value of loadings
l$Abs.Comp.2 %>% sort(decreasing = TRUE) -> top # save in order from highest down
l[which(l$Abs.Comp.2==top[1]),c("Variable", "Comp.2")] # top loading = CA
l[which(l$Abs.Comp.2==top[2]),c("Variable", "Comp.2")] # 2nd top loading = Ni.ppm

# PC3 (pH, Ni)
l$Abs.Comp.3 <- abs(l$Comp.3) # absolute value of loadings
l$Abs.Comp.3 %>% sort(decreasing = TRUE) -> top # save in order from highest down
l[which(l$Abs.Comp.3==top[1]),c("Variable", "Comp.3")] # top loading = PH
l[which(l$Abs.Comp.3==top[2]),c("Variable", "Comp.3")] # 2nd top loading = Ni.ppm

# PC4 (Phos, MG)
l$Abs.Comp.4 <- abs(l$Comp.4) # absolute value of loadings
l$Abs.Comp.4 %>% sort(decreasing = TRUE) -> top # save in order from highest down
l[which(l$Abs.Comp.4==top[1]),c("Variable", "Comp.4")] # top loading = Phos
l[which(l$Abs.Comp.4==top[2]),c("Variable", "Comp.4")] # 2nd top loading = MG



# *** D2: MIC data by soil traits #########################################

# is MIC predicted by soil chemistry at a site?
# as predictors, try Ni.ppm, Ca:Mg ratio, and PC2
# note, need pcscores dataframe from PCA section to incorporate PC1 and PC2 into data

setwd("C:/Users/camil/Washington State University (email.wsu.edu)/Plant Microbe Lab - Camille Wendlandt/Meso Nickel Tolerance Project/MASTERSHEET/")
D2 = read.csv("Delgado_Meso_Mastersheet_by_Replicate_2023-01-30.csv", header = T)

D2 <- complete_fun(D2, "Ni.MIC")
D2 <- complete_fun(D2, "Ni.ppm")
D2 %>% nrow() # 2034
D2$Unique.ID %>% unique %>% length # 668

# get familiar with structure of Nickel MIC data
D2$Ni.MIC %>% unique %>% sort # 1 2 3 4 5 6
D2$Unique.ID %>% unique %>% length # 671 strains

# merge in PCscores
D2 <- merge(D2, pcscores, by = "SAMPLEID", all.x = T)
D2 %>% dim # 2034 118
D2[which(!is.na(D2$PC1)),]$Unique.ID %>% unique %>% length # 457 strains have a soil SAMPLEID, and so got matched with PC scores
D2[which(!is.na(D2$SAMPLEID)),]$SAMPLEID %>% unique %>% length # 109 unique SAMPLEIDs




# *** D3: nre genotype data by soil traits ##################################

# does soil nickel level affect the probability of a strain having the Nre locus?
# note, need pcscores dataframe from PCA section to incorporate PC1 and PC2 into data

setwd("C:/Users/camil/Washington State University (email.wsu.edu)/Plant Microbe Lab - Camille Wendlandt/Meso Nickel Tolerance Project/MASTERSHEET/")
D3 = read.csv("Delgado_Meso_Mastersheet_by_Replicate_2023-01-30.csv", header = T)
dim(D3) # 2308 116

# keep just first record (row) for each Unique.ID (get rid of replicate info)
D3 <- D3[!duplicated(D3$Unique.ID),]
D3 %>% dim # 715 116
D3$Unique.ID %>% duplicated %>% any # FALSE; good

# add PC scores to data
D3 <- merge(D3, pcscores, by = "SAMPLEID", all.x = T)
D3 %>% dim # 715 118
D3[which(!is.na(D3$PC1)),] %>% nrow # 469 strains have a soil SAMPLEID, and so got matched with PC scores
D3[which(!is.na(D3$SAMPLEID)),]$SAMPLEID %>% unique %>% length # 109 unique SAMPLEIDs; good

# note: all strains that have nreA have nreX, and all that have nreX have nreA
# do all strains with nreA have nreX?
D3[which(D3$nreX=="yes"),]$nreA %>% unique # "yes"
D3[which(D3$nreX=="no"),]$nreA %>% unique # "no"

# change nreA genotypes to binary data
D3$nreA %>% unique # no yes NA
D3[which(D3$nreA=="yes"),]$nreA = 1
D3[which(D3$nreA=="no"),]$nreA = 0
D3$nreA <- as.numeric(D3$nreA)
D3$nreA %>% unique # 0 1 NA

# change nreY genotypes to binary data
D3$nreY %>% unique # no yes NA
D3[which(D3$nreY=="yes"),]$nreY = 1
D3[which(D3$nreY=="no"),]$nreY = 0
D3$nreY <- as.numeric(D3$nreY)
D3$nreY %>% unique # 0 1 NA



# FIGURES ##########################################################


# *** Fig 1a: PCA biplot #### 

D1$Serp <- ifelse(D1$CA/D1$MG<1,"Serp", "Non-serp")
D1 %>% nrow # 114

# note: to make this plot you need to install ggbiplot; use code below to install
# install.packages("devtools")
# library(devtools)
# install_github("vqv/ggbiplot")
library(ggbiplot)

g <- ggbiplot(pca.soil, 
              obs.scale = 1, 
              var.scale = 1, 
              # groups = D1$Serp, ellipse = TRUE, # this line of code isn't working
              var.axes = FALSE,
              varname.adjust = 2,
              varname.size = 6) +
  xlab("PC1") +
  ylab("PC2") +
  geom_point(aes(shape=D1$Serp, color = D1$Ni.ppm), size = 3) +
  scale_shape_discrete(name  ="Soil type",
                       breaks=c("Non-serp", "Serp"),
                       labels=c("Non-serp", "Serp")) +
  scale_colour_gradient(name = "Nickel (ppm)",
                        low = "darkolivegreen1",
                        high = "darkolivegreen") +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size = 15, color = "black"),
        axis.title = element_text(size = 20),
        legend.position = "right",
        legend.title = element_text(size = 20),
        legend.text=element_text(size = 15),
        legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.5, 'cm'))
print(g)

# if plot does not print, run line of code below and try again
# dev.off()

#setwd("C:/Users/camil/Washington State University (email.wsu.edu)/Plant Microbe Lab - Camille Wendlandt/Meso Nickel Tolerance Project/Figures/")
#ggsave("Fig1a.pdf",  g,  width = 5, height = 3, dpi = 1200)

# need to unload ggbiplot or other plots may not print
detach(package:ggbiplot, unload==TRUE)



# *** Fig 1c: Nickel vs serpentine status ####
D1$Serp <- ifelse(D1$CA/D1$MG<1,"Serp", "Non-serp")
D1 %>% nrow # 114
D1$SAMPLEID %>% unique %>% length # 114

temp <- D1 %>%
  group_by(Serp) %>%
  dplyr::summarize(CIup = CI(Ni.ppm), CIlow = CI2(Ni.ppm), Ni.ppm = mean(Ni.ppm, na.rm = TRUE))
head(temp)

barcolors <- c("darkolivegreen1", "darkolivegreen")

p1 <- ggplot(D1, aes(x = Serp, y = Ni.ppm)) + 
  xlab("") +
  ylab("Soil nickel (ppm)") +
  geom_bar(data = temp, stat = "identity", color = "black", fill = barcolors) +
  geom_errorbar(data = temp, aes(ymin=CIlow, ymax=CIup), width = 0.2,
                position = position_dodge(0.6), lwd = 0.8) +
  geom_jitter(pch = 16, size = 1.5, color = "gray20") +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.text.x =element_text(size = 20, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 20))
print(p1)

#setwd("C:/Users/camil/Washington State University (email.wsu.edu)/Plant Microbe Lab - Camille Wendlandt/Meso Nickel Tolerance Project/Figures/")
#ggsave("Fig1c.pdf",  p1,  width = 5, height = 3, dpi = 1200)


# *** Fig S1a: MIC vs Ni ####

usedata <- complete_fun(D2, "Ni.MIC")
usedata <- complete_fun(usedata, "Ni.ppm")
usedata %>% nrow() # 2034
usedata$Unique.ID %>% unique %>% length # 668

# summarize raw data by strain, getting a confidence interval for each MIC value
strainmeans <- usedata %>%
  group_by(Unique.ID, Ni.ppm) %>%
  dplyr::summarize(Mean = mean(Ni.MIC), CIup = CI(Ni.MIC), CIlow = CI2(Ni.MIC))

head(strainmeans)
strainmeans %>% nrow # 668 strains

p <- ggplot(strainmeans, aes(x = Ni.ppm, y = Mean))

p1 <- p + 
  geom_point(col = "gray60") +
  geom_smooth(method = "lm") +
  theme(panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text= element_text(size=20),
        axis.title = element_text(size=20)) +
  labs(x = "Soil nickel (ppm)", 
       y = "Nickel MIC (mM)")
# p1 <- arrangeGrob(p1, top = textGrob("A", x = unit(0.02, "npc"),
#                                      y = unit(0.95, "npc"), just=c("left","top"),
#                                      gp=gpar(col="black", fontsize=20)))

print(p1)

setwd("C:/Users/camil/Washington State University (email.wsu.edu)/Plant Microbe Lab - Camille Wendlandt/Meso Nickel Tolerance Project/Figures/")
ggsave(paste0("FigS1a_", Sys.Date(), ".pdf"),  p1,  width = 3, height = 3, dpi = 1200)


# *** Fig S1b: MIC vs Ca:Mg ####

usedata <- complete_fun(D2, "Ni.MIC")
usedata <- complete_fun(usedata, "CA.MG.RATIO")
usedata %>% nrow() # 2034
usedata$Unique.ID %>% unique %>% length # 668 strains

# summarize raw data by strain, getting a confidence interval for each MIC value
strainmeans <- usedata %>%
  group_by(Unique.ID, log(CA.MG.RATIO)) %>%
  dplyr::summarize(Mean = mean(Ni.MIC), CIup = CI(Ni.MIC), CIlow = CI2(Ni.MIC))

head(strainmeans)
strainmeans %>% nrow # 668 strains

p <- ggplot(strainmeans, aes(x = strainmeans$`log(CA.MG.RATIO)`, y = Mean))

p1 <- p + 
  geom_point(col = "gray60") +
  geom_smooth(method = "lm") +
  theme(panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text= element_text(size=20),
        axis.title = element_text(size=20)) +
  labs(x = "log(Ca:Mg ratio)", 
       y = "Nickel MIC (mM)")
# p1 <- arrangeGrob(p1, top = textGrob("B", x = unit(0.02, "npc"),
#                                      y = unit(0.95, "npc"), just=c("left","top"),
#                                      gp=gpar(col="black", fontsize=20)))

print(p1)

setwd("C:/Users/camil/Washington State University (email.wsu.edu)/Plant Microbe Lab - Camille Wendlandt/Meso Nickel Tolerance Project/Figures/")
ggsave(paste0("FigS1b_", Sys.Date(), ".pdf"),  p1,  width = 3, height = 3, dpi = 1200)



# *** Fig S1c: MIC vs PC2 ####

usedata <- complete_fun(D2, "Ni.MIC")
usedata <- complete_fun(usedata, "PC2")
usedata %>% nrow() # 1404
usedata$Unique.ID %>% unique %>% length # 457

# summarize raw data by strain, getting a confidence interval for each MIC value
strainmeans <- usedata %>%
  group_by(Unique.ID, PC2) %>%
  dplyr::summarize(Mean = mean(Ni.MIC), CIup = CI(Ni.MIC), CIlow = CI2(Ni.MIC))

head(strainmeans)
strainmeans %>% nrow # 457 strains

p <- ggplot(strainmeans, aes(x = PC2, y = Mean))

p1 <- p + 
  geom_point(col = "gray60") +
  geom_smooth(method = "lm") +
  theme(panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        axis.text= element_text(size=20),
        axis.title = element_text(size=20)) +
  labs(x = "Soil PC2", 
       y = "Nickel MIC (mM)")
# p1 <- arrangeGrob(p1, top = textGrob("C", x = unit(0.02, "npc"),
#                                      y = unit(0.95, "npc"), just=c("left","top"),
#                                      gp=gpar(col="black", fontsize=20)))

print(p1)

setwd("C:/Users/camil/Washington State University (email.wsu.edu)/Plant Microbe Lab - Camille Wendlandt/Meso Nickel Tolerance Project/Figures/")
ggsave(paste0("FigS1c_", Sys.Date(), ".pdf"),  p1,  width = 3, height = 3, dpi = 1200)


# *** Fig S10a: nreAX vs Ni ####

usedata <- complete_fun(D3, "nreA")
usedata %>% nrow # 302
usedata <- complete_fun(usedata, "Ni.ppm")
usedata %>% nrow # 302 strains total
usedata[which(usedata$nreA==1),] %>% nrow # 208 strains have nreA

# model
M1 <- glmer(nreA ~ Ni.ppm + PlantSpecies.Strain + (1|Site) + (1|Reserve), data = usedata, family = binomial)

# save model terms
intercept <- summary(M1)$coefficients[1]
slope <- summary(M1)$coefficients[2]

# show p-value for Nickel term
drop1(M1, test = "Chisq")
# npar    AIC    LRT   Pr(Chi)    
# <none>                   232.06                     
# Ni.ppm                 1 278.26 48.198 3.853e-12 ***
# PlantSpecies.Strain    1 230.08  0.019    0.8901    


# calculate mean + SE probability of having nre for various nickel bins
breaks = with(usedata, seq(min(Ni.ppm), max(Ni.ppm), length = 9)) #define 8 bins, using a sequence of 9 breaks
cut.Ni = cut(usedata$Ni.ppm, breaks = breaks) #make a binning factor from Ni.ppm
means = with(usedata, tapply(nreA, cut.Ni, mean)) #calculate the proportion of presence by bin
ses = with(usedata, tapply(nreA, cut.Ni, binomial.SE)) #calculate the standard error for each bin

# set where pdf image will be saved
setwd("C:/Users/camil/Washington State University (email.wsu.edu)/Plant Microbe Lab - Camille Wendlandt/Meso Nickel Tolerance Project/Figures/")
pdf(file = paste0("FigS10a_", Sys.Date(), ".pdf"), width = 6, height = 5)


#plot the bin means
plot(means ~ breaks[1:8], pch = 19, ylim = c(0,1), xlim =
       range(usedata$Ni.ppm), xlab = 'Soil Ni (ppm)', ylab = 'Presence of nreAX in genome', col = 'red')

#plot the bin SEs
segments(breaks[1:8], means+ses, breaks[1:8], means-ses)

#plot the raw data
with(usedata, points(nreA ~ Ni.ppm))

#plot the fitted curve from model
curve(logistic(intercept+slope*x), add = T, col = 'blue', lwd = 3)

# add text annotation of p-value
text(x = 60, y = 0.6, labels="P < 0.0001")

# save pdf image
dev.off()




# *** Fig S10b: nreY vs Ni ####

usedata <- complete_fun(D3, "nreY")
usedata %>% nrow # 302
usedata <- complete_fun(usedata, "Ni.ppm")
usedata %>% nrow # 302 strains total
usedata[which(usedata$nreY==1),] %>% nrow # 164 strains have nreY

# model
M1 <- glmer(nreY ~ Ni.ppm + PlantSpecies.Strain + (1|Site) + (1|Reserve), data = usedata, family = binomial)

# save model terms
intercept <- summary(M1)$coefficients[1]
slope <- summary(M1)$coefficients[2]

# show p-value for Nickel term
drop1(M1, test = "Chisq")
# npar    AIC    LRT   Pr(Chi)    
# <none>                   257.60                     
# Ni.ppm                 1 292.38 36.780 1.322e-09 ***
# PlantSpecies.Strain    1 255.77  0.171    0.6793 

# calculate mean + SE probability of having nre for various nickel bins
breaks = with(usedata, seq(min(Ni.ppm), max(Ni.ppm), length = 9)) #define 8 bins, using a sequence of 9 breaks
cut.Ni = cut(usedata$Ni.ppm, breaks = breaks) #make a binning factor from Ni.ppm
means = with(usedata, tapply(nreY, cut.Ni, mean)) #calculate the proportion of presence by bin
ses = with(usedata, tapply(nreY, cut.Ni, binomial.SE)) #calculate the standard error for each bin

# set where pdf image will be saved
setwd("C:/Users/camil/Washington State University (email.wsu.edu)/Plant Microbe Lab - Camille Wendlandt/Meso Nickel Tolerance Project/Figures/")
pdf(file = paste0("FigS10b_", Sys.Date(), ".pdf"), width = 6, height = 5)


#plot the bin means
plot(means ~ breaks[1:8], pch = 19, ylim = c(0,1), xlim =
       range(usedata$Ni.ppm), xlab = 'Soil Ni (ppm)', ylab = 'Presence of nreY in genome', col = 'red')

#plot the bin SEs
segments(breaks[1:8], means+ses, breaks[1:8], means-ses)

#plot the raw data
with(usedata, points(nreY ~ Ni.ppm))

#plot the fitted curve from model
intercept <- summary(M1)$coefficients[1]
slope <- summary(M1)$coefficients[2]
curve(logistic(intercept+slope*x), add = T, col = 'blue', lwd = 3)

# add text annotation of p-value
text(x = 60, y = 0.6, labels="P < 0.0001")

# save pdf image
dev.off()





# ANALYSES #############################################################


### *** MIC ~ Ni.ppm (Table S1) ####

# set usedata
usedata <- complete_fun(D2, "Ni.MIC")
usedata <- complete_fun(usedata, "Ni.ppm")
usedata %>% nrow() # 2034
usedata$Unique.ID %>% unique %>% length # 668


M1 <- lmer(Ni.MIC ~ Ni.ppm + PlantSpecies.Strain + 
              (1|Reserve) + (1|Site) + (1|Unique.ID), 
            data = usedata)

# check singularity
model = M1
isSingular(model, tol = 1e-05) # FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) # e-05, converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# ok

# look at parameter estimates
summary(M1)

# test main effects
drop1(M1, test = "Chisq")
# npar    AIC    LRT   Pr(Chi)    
# <none>                   4854.0                     
# Ni.ppm                 1 4885.7 33.679 6.499e-09 ***
# PlantSpecies.Strain    1 4852.1  0.090    0.7638 

# test significance of Site
M1.noSite <- update(M1, .~. -(1|Site))
anova(M1, M1.noSite)
#           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# M1.noSite    6 4884.4 4918.1 -2436.2   4872.4                         
# M1           7 4854.0 4893.3 -2420.0   4840.0 32.359  1  1.282e-08 ***

# test significance of Reserve 
M1.noRes <- update(M1, .~. -(1|Reserve))
anova(M1, M1.noRes)
#          npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)
# M1.noRes    6 4852.1 4885.8 -2420.1   4840.1                    
# M1          7 4854.0 4893.3 -2420.0   4840.0 0.106  1     0.7448

# test significance of Strain
M1.noStrain <- update(M1, .~. -(1|Unique.ID))
anova(M1, M1.noStrain)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# M1.noStrain    6 5962.2 5995.9 -2975.1   5950.2                         
# M1             7 4854.0 4893.3 -2420.0   4840.0 1110.2  1  < 2.2e-16 ***


### *** MIC ~ Ca:Mg (Table S1) ####

# set usedata
usedata <- complete_fun(D2, "Ni.MIC")
usedata <- complete_fun(usedata, "CA.MG.RATIO")
usedata %>% nrow() # 2034
usedata$Unique.ID %>% unique %>% length # 668 strains


M1 <- lmer(Ni.MIC ~ log(CA.MG.RATIO) + PlantSpecies.Strain + 
             (1|Reserve) + (1|Site) + (1|Unique.ID), 
           data = usedata)

# check singularity
model = M1
isSingular(model, tol = 1e-05) # FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) # converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# ok, tried poisson, but poisson does not improve

# look at parameter estimates
summary(M1)

# test main effects
drop1(M1, test = "Chisq")
#                     npar    AIC    LRT   Pr(Chi)    
# <none>                   4821.5                     
# log(CA.MG.RATIO)       1 4885.7 66.224 4.025e-16 ***
# PlantSpecies.Strain    1 4819.5  0.090    0.7642 

# test significance of Site
M1.noSite <- update(M1, .~. -(1|Site))
anova(M1, M1.noSite)
#           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)   
# M1.noSite    6 4829.0 4862.7 -2408.5   4817.0                        
# M1           7 4821.5 4860.8 -2403.7   4807.5 9.5021  1   0.002052 **

# test significance of Reserve
M1.noRes <- update(M1, .~. -(1|Reserve))
anova(M1, M1.noRes)
#          npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)
# M1.noRes    6 4819.5 4853.2 -2403.7   4807.5                    
# M1          7 4821.5 4860.8 -2403.7   4807.5     0  1          1

# test significance of Strain
M1.noStrain <- update(M1, .~. -(1|Unique.ID))
anova(M1, M1.noStrain)
#             npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# M1.noStrain    6 5915.7 5949.4 -2951.8   5903.7                         
# M1             7 4821.5 4860.8 -2403.7   4807.5 1096.2  1  < 2.2e-16 ***
  

### *** MIC ~ PC2 (Table S1) ####

# set usedata
usedata <- complete_fun(D2, "Ni.MIC")
usedata <- complete_fun(usedata, "PC2")
usedata %>% nrow() # 1404
usedata$Unique.ID %>% unique %>% length # 457


M1 <- lmer(Ni.MIC ~ PC2 + PlantSpecies.Strain + 
             (1|Reserve) + (1|Site) + (1|Unique.ID), 
           data = usedata)

# check singularity
model = M1
isSingular(model, tol = 1e-05) # FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad)) # e-05, converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# ok

# look at parameter estimates
summary(M1)

# test main effects
drop1(M1, test = "Chisq")
#                     npar    AIC     LRT   Pr(Chi)    
# <none>                   3338.8                      
# PC2                    1 3362.3 25.4692 4.495e-07 ***
# PlantSpecies.Strain    1 3337.0  0.2687    0.6042    
 

# test significance of Site
M1.noSite <- update(M1, .~. -(1|Site))
anova(M1, M1.noSite)
#           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# M1.noSite    6 3355.5 3387.0 -1671.8   3343.5                         
# M1           7 3338.8 3375.5 -1662.4   3324.8 18.717  1  1.516e-05 ***

# test significance of Reserve
M1.noRes <- update(M1, .~. -(1|Reserve))
anova(M1, M1.noRes)
#          npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)
# M1.noRes    6 3336.8 3368.3 -1662.4   3324.8                    
# M1          7 3338.8 3375.5 -1662.4   3324.8     0  1     0.9946

# test significance of Strain
M1.noStrain <- update(M1, .~. -(1|Unique.ID))
anova(M1, M1.noStrain)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# M1.noStrain    6 4132.2 4163.7 -2060.1   4120.2                         
# M1             7 3338.8 3375.5 -1662.4   3324.8 795.47  1  < 2.2e-16 ***


### *** nreAX ~ Ni.ppm (Table S7) ####

usedata <- complete_fun(D3, "nreA")
usedata %>% nrow # 302
usedata <- complete_fun(usedata, "Ni.ppm")
usedata %>% nrow # 302 strains total
usedata[which(usedata$nreA==1),] %>% nrow # 208 strains have nreA

M1 <- glmer(nreA ~ Ni.ppm + PlantSpecies.Strain + (1|Site) + (1|Reserve), data = usedata, family = binomial)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# e-05, converged

# look at residuals
model = M1
simulationOutput <- simulateResiduals(fittedModel = model)
plot(simulationOutput)
# ok residuals

summary(M1)


# log-odds of having nreA increases by beta for every 1 ppm increase in soil nickel
# so the odds of having nreA are multipled by exp(beta) for every 1 ppm increase in soil nickel

exp(0.13816) # 1.148159

# test main effects
options(digits = 10)
drop1(M1, test = "Chisq")
#                     npar       AIC       LRT    Pr(Chi)    
# <none>                   232.06312                         
# Ni.ppm                 1 278.26081 48.197687 3.8534e-12 ***
# PlantSpecies.Strain    1 230.08222  0.019100    0.89008  

# test significance of Site
M1.noSite <- update(M1, .~. -(1|Site))
anova(M1, M1.noSite)
#           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)   
# M1.noSite    4 240.81 255.65 -116.41   232.81                        
# M1           5 232.06 250.62 -111.03   222.06 10.749  1   0.001044 **

# test significance of Reserve
M1.noRes <- update(M1, .~. -(1|Reserve))
anova(M1, M1.noRes)
# npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)
# M1.noRes    4 230.06 244.91 -111.03   222.06                    
# M1          5 232.06 250.62 -111.03   222.06     0  1     0.9978



# *** nreAX ~ Ca:Mg (Table S7) #### 

usedata <- complete_fun(D3, "nreA")
usedata %>% nrow # 302
usedata <- complete_fun(usedata, "CA.MG.RATIO")
usedata %>% nrow # 302 strains total
usedata[which(usedata$nreA==1),] %>% nrow # 208 strains have nreA

M1 <- glmer(nreA ~ log(CA.MG.RATIO) + PlantSpecies.Strain + (1|Site) + (1|Reserve), data = usedata, family = binomial)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# e-05, converged

# look at residuals
model = M1
simulationOutput <- simulateResiduals(fittedModel = model)
plot(simulationOutput)
# great residuals

summary(M1)

# log-odds of having nreA increases by beta for every 1 ppm increase in soil nickel
# so the odds of having nreA are multipled by exp(beta) for every 1 unit increase in PC2

# test main effects
drop1(M1, test = "Chisq")
#                     npar    AIC    LRT   Pr(Chi)    
# <none>                   225.73                     
# log(CA.MG.RATIO)       1 278.26 54.534 1.528e-13 ***
# PlantSpecies.Strain    1 224.31  0.580    0.4463  

# test significance of Site
M1.noSite <- update(M1, .~. -(1|Site))
anova(M1, M1.noSite)
#           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# M1.noSite    4 240.87 255.71 -116.43   232.87                         
# M1           5 225.73 244.28 -107.86   215.73 17.141  1   3.47e-05 ***

# test significance of Reserve
M1.noRes <- update(M1, .~. -(1|Reserve))
anova(M1, M1.noRes)
# npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)
# M1.noRes    4 223.73 238.57 -107.86   215.73                    
# M1          5 225.73 244.28 -107.86   215.73     0  1     0.9983




# *** nreAX ~ PC2 (Table S7) #### 

usedata <- complete_fun(D3, "nreA")
usedata %>% nrow # 302
usedata <- complete_fun(usedata, "PC2")
usedata %>% nrow # 221 strains total
usedata[which(usedata$nreA==1),] %>% nrow # 148 strains have nreA


M1 <- glmer(nreA ~ PC2 + PlantSpecies.Strain + (1|Site) + (1|Reserve), data = usedata, family = binomial)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# e-06, converged

# look at residuals
model = M1
simulationOutput <- simulateResiduals(fittedModel = model)
plot(simulationOutput)
# ok residuals

summary(M1)

# log-odds of having nreA increases by beta for every 1 ppm increase in soil nickel
# so the odds of having nreA are multipled by exp(beta) for every 1 unit increase in PC2

# test main effects
drop1(M1, test = "Chisq")
# npar    AIC     LRT   Pr(Chi)    
# <none>                   205.10                      
# PC2                    1 224.35 21.2510 4.029e-06 ***
# PlantSpecies.Strain    1 203.12  0.0189    0.8906  

# test significance of Site
M1.noSite <- update(M1, .~. -(1|Site))
anova(M1, M1.noSite)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# M1.noSite    4 215.56 229.16 -103.78   207.56                         
# M1           5 205.10 222.09  -97.55   195.10 12.463  1  0.0004152 ***

# test significance of Reserve
M1.noRes <- update(M1, .~. -(1|Reserve))
anova(M1, M1.noRes)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# M1.noRes    4 204.09 217.68 -98.043   196.09                     
# M1          5 205.10 222.09 -97.550   195.10 0.9853  1     0.3209




### *** nreY ~ Ni.ppm (Table S8) ####

usedata <- complete_fun(D3, "nreY")
usedata %>% nrow # 302
usedata <- complete_fun(usedata, "Ni.ppm")
usedata %>% nrow # 302 strains total
usedata[which(usedata$nreY==1),] %>% nrow # 164 strains have nreY

M1 <- glmer(nreY ~ Ni.ppm + PlantSpecies.Strain + (1|Site) + (1|Reserve), data = usedata, family = binomial)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# e-06, converged

# look at residuals
model = M1
simulationOutput <- simulateResiduals(fittedModel = model)
plot(simulationOutput)
# ok residuals

summary(M1)



# log-odds of having nreY increases by beta for every 1 ppm increase in soil nickel
# so the odds of having nreY are multipled by exp(beta) for every 1 ppm increase in soil nickel

exp(0.1044) # 1.110044

# test main effects
drop1(M1, test = "Chisq")
                    # npar    AIC    LRT   Pr(Chi)    
# <none>                   257.60                     
# Ni.ppm                 1 292.38 36.780 1.322e-09 ***
# PlantSpecies.Strain    1 255.77  0.171    0.6793 

# test significance of Site
M1.noSite <- update(M1, .~. -(1|Site))
anova(M1, M1.noSite)
#           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# M1.noSite    4 278.77 293.61 -135.38   270.77                         
# M1           5 257.60 276.15 -123.80   247.60 23.175  1  1.479e-06 ***

# test significance of Reserve
M1.noRes <- update(M1, .~. -(1|Reserve))
anova(M1, M1.noRes)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# M1.noRes    4 256.26 271.11 -124.13   248.26                     
# M1          5 257.60 276.15 -123.80   247.60 0.6678  1     0.4138


### *** nreY ~ Ca:Mg (Table S8) ####

usedata <- complete_fun(D3, "nreY")
usedata %>% nrow # 302
usedata <- complete_fun(usedata, "CA.MG.RATIO")
usedata %>% nrow # 302 strains total
usedata[which(usedata$nreY==1),] %>% nrow # 164 strains have nreY

M1 <- glmer(nreY ~ log(CA.MG.RATIO) + PlantSpecies.Strain + (1|Site) + (1|Reserve), data = usedata, family = binomial)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# e-06, converged

# look at residuals
model = M1
simulationOutput <- simulateResiduals(fittedModel = model)
plot(simulationOutput)
# good residuals

summary(M1)

# log-odds of having nreY increases by beta for every 1 ppm increase in soil nickel
# so the odds of having nreY are multipled by exp(beta) for every 1 ppm increase in soil nickel

# test main effects
drop1(M1, test = "Chisq")
#                     npar       AIC       LRT    Pr(Chi)    
# <none>                   232.11037                         
# log(CA.MG.RATIO)       1 292.37617 62.265798 3.0009e-15 ***
# PlantSpecies.Strain    1 230.46039  0.350022     0.5541    

# test significance of Site
M1.noSite <- update(M1, .~. -(1|Site))
anova(M1, M1.noSite)
# npar       AIC       BIC     logLik  deviance   Chisq Df Pr(>Chisq)    
# M1.noSite    4 250.69307 265.53478 -121.34653 242.69307                          
# M1           5 232.11037 250.66250 -111.05518 222.11037 20.5827  1  5.711e-06 ***

# test significance of Reserve
M1.noRes <- update(M1, .~. -(1|Reserve))
anova(M1, M1.noRes)
# npar       AIC       BIC     logLik  deviance   Chisq Df Pr(>Chisq)
# M1.noRes    4 230.12203 244.96374 -111.06102 222.12203                      
# M1          5 232.11037 250.66250 -111.05518 222.11037 0.01167  1    0.91399


### *** nreY ~ PC2 (Table S8) ####

usedata <- complete_fun(D3, "nreY")
usedata %>% nrow # 302
usedata <- complete_fun(usedata, "PC2")
usedata %>% nrow # 221 strains total
usedata[which(usedata$nreY==1),] %>% nrow # 114 strains have nreY

M1 <- glmer(nreY ~ PC2 + PlantSpecies.Strain + (1|Site) + (1|Reserve), data = usedata, family = binomial)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# e-05, converged

# look at residuals
model = M1
simulationOutput <- simulateResiduals(fittedModel = model)
plot(simulationOutput)
# ok residuals

summary(M1)

# log-odds of having nreY increases by beta for every 1 ppm increase in soil nickel
# so the odds of having nreY are multipled by exp(beta) for every 1 ppm increase in soil nickel

# test main effects
drop1(M1, test = "Chisq")
#                     npar    AIC     LRT   Pr(Chi)    
# <none>                   199.76                      
# PC2                    1 219.20 21.4399 3.651e-06 ***
# PlantSpecies.Strain    1 197.79  0.0309    0.8604   

# test significance of Site
M1.noSite <- update(M1, .~. -(1|Site))
anova(M1, M1.noSite)
#           npar    AIC    BIC   logLik deviance  Chisq Df Pr(>Chisq)    
# M1.noSite    4 214.08 227.67 -103.040   206.08                         
# M1           5 199.76 216.75  -94.879   189.76 16.322  1  5.344e-05 ***

# test significance of Reserve
M1.noRes <- update(M1, .~. -(1|Reserve))
anova(M1, M1.noRes)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)  
# M1.noRes    4 201.19 214.79 -96.597   193.19                       
# M1          5 199.76 216.75 -94.879   189.76 3.4375  1    0.06373 .




