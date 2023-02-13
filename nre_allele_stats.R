#Delgado et al. 2023 analyses 
#Angeliqua Montoya
#2/1/2023 Final check (ie fix nreAXY allele types analysis and run all code one more time)
#use most figs I made yesterday but just new ones for fig 3B,C,D after saving data 


setwd("C:/Users/13605/Desktop/Lab/MIC.UV.analysis_Spring 2021/Delgado 2023_final figs and data")

#load libraries 
library(lme4)
library(DHARMa)
library(ggplot2)
library(dplyr)
library(emmeans)
library(tidyr)
library(stringr)
library(gridExtra)
library(plyr)


#Meso Nickel Mastersheet
gegd = read.csv("Delgado_Meso_Mastersheet_by_Replicate_2023-01-30.csv")

#change relevant columns to factors
gegd$Unique.ID = as.factor(gegd$Unique.ID)
gegd$nreA = as.factor(gegd$nreA)
gegd$nreX = as.factor(gegd$nreX)
gegd$nreY = as.factor(gegd$nreY)
gegd$nreA_alleles = as.factor(gegd$nreA_alleles)
gegd$nreX_alleles = as.factor(gegd$nreX_alleles)
gegd$nreY_alleles = as.factor(gegd$nreY_alleles)
gegd$Clade = as.factor(gegd$Clade)
gegd$SerpentineType = as.factor(gegd$SerpentineType)

#change soil type from "n" "s" to Serp, Non-serp
gegd$SerpentineType = ifelse(gegd$SerpentineType =="N", "Non-serp","Serp")

#### Growth tradeoff's analysis: Soil type (Serp/Non-serp) and nreA, X, or Y presence/absence ####

#### data prep 
head(gegd)

#Now change data format so growth data (OD) is long (not OD in Ni and OD in no Ni in separate columns)
# growth in nickel
g.ni = subset(gegd, select = c(Unique.ID, OD72hr.ni))
head(g.ni)
colnames(g.ni)[2]= "OD"
# add column for nickel mM
g.ni$ni.mm = "1"
str(g.ni)

# growth in no nickel
g.no = subset(gegd, select = c(Unique.ID, OD72hr.no))
head(g.no)
colnames(g.no)[2]= "OD"
#add column for nickel mM
g.no$ni.mm = "0"
str(g.no)

#before merging, change Unique.ID to character
g.no$Unique.ID = as.character(g.no$Unique.ID)
g.ni$Unique.ID = as.character(g.ni$Unique.ID)

#combine ODlong growth data
OD.long=rbind(g.ni, g.no, by="Unique.ID")
head(OD.long)
str(OD.long)

#change Unique.ID back to factor
OD.long$Unique.ID = as.factor(OD.long$Unique.ID)

#need to subset gegd again (keep one row per Unique.ID)
# remove duplicate Unique.IDs
gegd3= gegd[!duplicated(gegd$Unique.ID),]
head(gegd3)
nrow(gegd3) #715 Unique.IDs 

g.long=merge(gegd3,OD.long, by="Unique.ID")
nrow(g.long) #4616 rows (each Unique.ID listed at least 6 times - 3 growth in ni, and 3 growth in no ni, but some Unique.IDs assayed extra times as controls or other have no growth data)
length(unique(g.long$Unique.ID)) #715 but some Unique.IDs missing growth data 

#how many Unique.IDs have growth data?
g.longc= g.long[complete.cases(g.long$OD),]
length(unique(g.longc$Unique.ID)) #296 Unique.IDs have growth data
nrow(g.longc) #expect at least 1776 rows, but 1818 (probably some Unique.IDs have extra rows- tested multiple times)
head(g.longc)

str(g.longc)
#change relevant columns 
g.longc$ni.mm = as.factor(g.longc$ni.mm)
g.longc$OD = as.numeric(g.longc$OD)


#### Make subsets for different clades 
# Note: g.longc = all clades 

#only clade 1
g.longc1 = g.longc[ which(g.longc$Clade=='1'),]

#only clade 2
g.longc2 = g.longc[ which(g.longc$Clade=='2'),]

#other clades (not clade 1 or 2)
g.longco = g.longc[ which(g.longc$Clade!='1' & g.longc$Clade!='2'),]

#drop levels
g.longc = droplevels(g.longc)

#### ****end of clade changes

length(unique(g.longc$Unique.ID)) 
#296 Unique.IDs in all clades
length(unique(g.longc1$Unique.ID))
#86 strains in clade 1
length(unique(g.longc2$Unique.ID))
#143 strains in clade 2
length(unique(g.longco$Unique.ID))
#67 strains in all other clades


# Calculate Unique.ID growth means (all clades)
ODs = aggregate(OD ~ Unique.ID*ni.mm, FUN=mean, data=g.longc)
ODs$se = aggregate(OD ~ Unique.ID*ni.mm, FUN=function(x) sd(x)/sqrt(length(x)), data=g.longc)[,3]
head(ODs)
#merge with gegd3 (Unique.IDs only listed once but has soil data)
Unique.ID.mean= merge(gegd3, ODs, by="Unique.ID", all.X=F, all.y=T)

head(Unique.ID.mean)
length(unique(Unique.ID.mean$Unique.ID)) #296 strains

#repeat for other clade subsets

#clade 1
ODs1 = aggregate(OD ~ Unique.ID*ni.mm, FUN=mean, data=g.longc1)
ODs1$se = aggregate(OD ~ Unique.ID*ni.mm, FUN=function(x) sd(x)/sqrt(length(x)), data=g.longc1)[,3]
head(ODs1)
#merge with gegd3 (Unique.IDs only listed once but has soil data)
Unique.ID.mean1= merge(gegd3, ODs1, by="Unique.ID", all.X=F, all.y=T)

head(Unique.ID.mean1)
length(unique(Unique.ID.mean1$Unique.ID)) #86 strains

#clade 2
ODs2 = aggregate(OD ~ Unique.ID*ni.mm, FUN=mean, data=g.longc2)
ODs2$se = aggregate(OD ~ Unique.ID*ni.mm, FUN=function(x) sd(x)/sqrt(length(x)), data=g.longc2)[,3]
head(ODs2)
#merge with gegd3 (Unique.IDs only listed once but has soil data)
Unique.ID.mean2= merge(gegd3, ODs2, by="Unique.ID", all.X=F, all.y=T)

head(Unique.ID.mean2)
length(unique(Unique.ID.mean2$Unique.ID)) #143 strains

#other clades 
ODso = aggregate(OD ~ Unique.ID*ni.mm, FUN=mean, data=g.longco)
ODso$se = aggregate(OD ~ Unique.ID*ni.mm, FUN=function(x) sd(x)/sqrt(length(x)), data=g.longco)[,3]
head(ODso)
#merge with gegd3 (Unique.IDs only listed once but has soil data)
Unique.ID.meano= merge(gegd3, ODso, by="Unique.ID", all.X=F, all.y=T)

head(Unique.ID.meano)
length(unique(Unique.ID.meano$Unique.ID)) #67 strains


#### Done prepping data

#### Soil type tradeoff's analysis ####

#### Serp/Non-serp model across all clades
M0=lm(OD ~ SerpentineType*ni.mm, na.action=na.omit, data=Unique.ID.mean) 
summary(M0)

#Residuals:
#     Min       1Q   Median       3Q      Max 
#-0.55082 -0.16802  0.01111  0.16151  0.64803 
#Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                0.69095    0.02386  28.953  < 2e-16 ***
#SerpentineTypeSerp        -0.02638    0.02875  -0.918    0.359    
#ni.mm1                    -0.22643    0.03375  -6.709 4.61e-11 ***
#SerpentineTypeSerp:ni.mm1  0.22845    0.04065   5.619 2.96e-08 ***
#Residual standard error: 0.2289 on 588 degrees of freedom
#Multiple R-squared:  0.09774,	Adjusted R-squared:  0.09314 
#F-statistic: 21.23 on 3 and 588 DF,  p-value: 4.535e-13


# check residuals
simulationOutput = simulateResiduals(fittedModel = M0)
plot(simulationOutput)
# overall good fit

# test fixed effects
drop1(M0, test = "F")

#Single term deletions
#Model:
#OD ~ SerpentineType * ni.mm
#                     Df Sum of Sq    RSS     AIC F value    Pr(>F)    
#<none>                            30.808 -1741.8                      
#SerpentineType:ni.mm  1    1.6545 32.463 -1712.8  31.577 2.962e-08 ****
  
# test which groups differ 
M0 %>% emmeans(~ SerpentineType*ni.mm) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp #holm adjustment for multiple comparisons

# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]

#                     contrast   estimate         SE  df   t.ratio      p.value
#2 (Non-serp 0) - (Non-serp 1)  0.2264283 0.03374938 588  6.709108 1.844285e-10
#4       Serp 0 - (Non-serp 1)  0.2000482 0.02874628 588  6.959098 4.584328e-11
#6       (Non-serp 1) - Serp 1 -0.2020661 0.02874628 588 -7.029297 3.465289e-11

# Results summary: Strains from serpentine have greater growth in the presence of nickel
# No tradeoff's detected: strains from serp or non-serp soil do not differ in growth in the absence of nickel

#### **Figure 1B ####
s.ns = ggplot(data = emmip(M0, SerpentineType ~ ni.mm, CIs=T, plotit=F), aes(x = ni.mm, y = yvar, group = SerpentineType, color=SerpentineType)) +
    geom_point(aes(shape=SerpentineType), size=3, position=position_dodge(0.05)) +
    scale_shape_manual(values=c(16, 17))+
    geom_line(size=1, position=position_dodge(0.05)) +
    scale_color_manual(values=c("darkolivegreen1", "darkolivegreen"))+
    geom_errorbar(aes(ymin = LCL, ymax = UCL), width=0.15, position=position_dodge(0.05))+
    ylab(expression(Fitness~(OD[600]))) + xlab("Nickel (mM)")+
    ylim(0.3, 0.9) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5), 
        axis.ticks.length = unit(0.1, "cm"),text = element_text(size=15, colour = "black"), axis.title = element_text(size = 20), 
        axis.text.x = element_text(size=15, vjust=0.95, colour = "black"), 
        axis.text.y = element_text(size=15,hjust=0.95, colour = "black"),
        legend.title = element_blank(), legend.text=element_text(size=20), legend.position="top")

s.ns 
#ggsave("Fig1b_updated_1.31.23.pdf",  s.ns,  width = 3.5, height = 3, dpi = 1200) 




#### nre tradeoff's analysis: Supplementary Figure S12 A-L####

#### ***nreA ####
MA=lm(OD ~ nreA*ni.mm, na.action=na.omit, data=Unique.ID.mean) 

summary(MA)
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-0.56442 -0.15663  0.00806  0.15648  0.65149 

#Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     0.666267   0.023498  28.355  < 2e-16 ***
#nreAyes         0.009389   0.028235   0.333     0.74    
#ni.mm1         -0.234592   0.033231  -7.060 4.73e-12 ***
#nreAyes:ni.mm1  0.239120   0.039931   5.988 3.69e-09 ***
#Residual standard error: 0.2242 on 588 degrees of freedom
#Multiple R-squared:  0.1348,	Adjusted R-squared:  0.1304 
#F-statistic: 30.53 on 3 and 588 DF,  p-value: < 2.2e-16


# check residuals
simulationOutput = simulateResiduals(fittedModel = MA)
plot(simulationOutput)
# overall good fit, some heterogeneity of variance

# test fixed effects
drop1(MA, test = "F")
#Single term deletions
#Model:
#OD ~ nreA * ni.mm
#           Df Sum of Sq    RSS     AIC F value    Pr(>F)    
#<none>                  29.544 -1766.6                      
#nreA:ni.mm  1    1.8018 31.345 -1733.6  35.861 3.691e-09 ***

# test which groups differ  
MA %>% emmeans(~ nreA*ni.mm) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp 

temp
# contrast      estimate     SE  df t.ratio p.value
# no 0 - yes 0  -0.00939 0.0282 588  -0.333  1.0000
# no 0 - no 1    0.23459 0.0332 588   7.060  <.0001
# no 0 - yes 1  -0.01392 0.0282 588  -0.493  1.0000
# yes 0 - no 1   0.24398 0.0282 588   8.641  <.0001
# yes 0 - yes 1 -0.00453 0.0221 588  -0.204  1.0000
# no 1 - yes 1  -0.24851 0.0282 588  -8.801  <.0001

# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]

#      contrast   estimate         SE  df   t.ratio      p.value
#2  no 0 - no 1  0.2345923 0.03323052 588  7.059545 1.890865e-11
#4 yes 0 - no 1  0.2439812 0.02823521 588  8.641027 2.633538e-16
#6 no 1 - yes 1 -0.2485088 0.02823521 588 -8.801381 9.017015e-17


# nreA reaction norm plot for supplementary figure S12
nreA.p = ggplot(data = emmip(MA, nreA ~ ni.mm, CIs=T, plotit=F), aes(x = ni.mm, y = yvar, group = nreA)) +
    geom_point() +
    geom_line(aes(linetype = nreA), size=1) +
    scale_linetype_manual(values=c("dashed", "solid"), labels=c("absent","present")) +
    geom_errorbar(aes(ymin = yvar - SE, ymax = yvar + SE), width=0.1)+
    ylab(expression(Fitness~(OD[600]))) + xlab("Nickel (mM)")+
    ylim(0.1, 0.95) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5), 
        axis.ticks.length = unit(0.1, "cm"),text = element_text(size=20), 
        axis.text.x = element_text(size=20, vjust=0.95), 
        axis.text.y = element_text(size=20,hjust=0.95),
        legend.key.size =  unit(0.4, "in"))

nreA.p
#ggsave("Fig_S12_A_updated_1.31.23.pdf",  nreA.p,  width = 5, height = 4, dpi = 1200)

# clade 1 only 
MA1=lm(OD ~ nreA*ni.mm, na.action=na.omit, data=Unique.ID.mean1) 

nreA.p1 = ggplot(data = emmip(MA1, nreA ~ ni.mm, CIs=T, plotit=F), aes(x = ni.mm, y = yvar, group = nreA)) +
    geom_point() +
    geom_line(aes(linetype = nreA), size=1) +
    scale_linetype_manual(values=c("dashed", "solid"), labels=c("absent","present")) +
    geom_errorbar(aes(ymin = yvar - SE, ymax = yvar + SE), width=0.1)+
    ylab(expression(Fitness~(OD[600]))) + xlab("Nickel (mM)")+
    ylim(0.1, 0.95) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5), 
        axis.ticks.length = unit(0.1, "cm"),text = element_text(size=20), 
        axis.text.x = element_text(size=20, vjust=0.95), 
        axis.text.y = element_text(size=20,hjust=0.95),
        legend.key.size =  unit(0.4, "in"))

nreA.p1

#ggsave("Fig_S12_B_updated_1.31.23.pdf",  nreA.p1,  width = 5, height = 4, dpi = 1200)

# clade 2
MA2=lm(OD ~ nreA*ni.mm, na.action=na.omit, data=Unique.ID.mean2) 

nreA.p2 = ggplot(data = emmip(MA2, nreA ~ ni.mm, CIs=T, plotit=F), aes(x = ni.mm, y = yvar, group = nreA)) +
    geom_point() +
    geom_line(aes(linetype = nreA), size=1) +
    scale_linetype_manual(values=c("dashed", "solid"), labels=c("absent","present")) +
    geom_errorbar(aes(ymin = yvar - SE, ymax = yvar + SE), width=0.1)+
    ylab(expression(Fitness~(OD[600]))) + xlab("Nickel (mM)")+
    ylim(0.1, 0.95) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5), 
        axis.ticks.length = unit(0.1, "cm"),text = element_text(size=20), 
        axis.text.x = element_text(size=20, vjust=0.95), 
        axis.text.y = element_text(size=20,hjust=0.95),
        legend.key.size =  unit(0.4, "in"))

nreA.p2

#ggsave("Fig_S12_C_updated_1.31.23.pdf",  nreA.p2,  width = 5, height = 4, dpi = 1200)

# other clades
MAo=lm(OD ~ nreA*ni.mm, na.action=na.omit, data=Unique.ID.meano) 
#plot
nreA.po = ggplot(data = emmip(MAo, nreA ~ ni.mm, CIs=T, plotit=F), aes(x = ni.mm, y = yvar, group = nreA)) +
    geom_point() +
    geom_line(aes(linetype = nreA), size=1) +
    scale_linetype_manual(values=c("dashed", "solid"), labels=c("absent","present")) +
    geom_errorbar(aes(ymin = yvar - SE, ymax = yvar + SE), width=0.1)+
    ylab(expression(Fitness~(OD[600]))) + xlab("Nickel (mM)")+
    ylim(0.1, 0.95) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5), 
        axis.ticks.length = unit(0.1, "cm"),text = element_text(size=20), 
        axis.text.x = element_text(size=20, vjust=0.95), 
        axis.text.y = element_text(size=20,hjust=0.95),
        legend.key.size =  unit(0.4, "in"))

nreA.po

#ggsave("Fig_S12_D_updated_1.31.23.pdf",  nreA.po,  width = 5, height = 4, dpi = 1200)



#### ***nreX ####
#model 
MX=lm(OD ~ nreX*ni.mm, na.action=na.omit, data=Unique.ID.mean) 

summary(MX)
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-0.56442 -0.15663  0.00806  0.15648  0.65149 

#Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     0.666267   0.023498  28.355  < 2e-16 ***
#nreXyes         0.009389   0.028235   0.333     0.74    
#ni.mm1         -0.234592   0.033231  -7.060 4.73e-12 ***
#nreXyes:ni.mm1  0.239120   0.039931   5.988 3.69e-09 ***

#Residual standard error: 0.2242 on 588 degrees of freedom
#Multiple R-squared:  0.1348,	Adjusted R-squared:  0.1304 
#F-statistic: 30.53 on 3 and 588 DF,  p-value: < 2.2e-16

# check residuals
simulationOutput = simulateResiduals(fittedModel = MX)
plot(simulationOutput)
# overall good fit

# test fixed effects
drop1(MX, test = "F")
#Single term deletions
#Model:
#OD ~ nreX * ni.mm
#           Df Sum of Sq    RSS     AIC F value    Pr(>F)    
#<none>                  29.544 -1766.6                      
#nreX:ni.mm  1    1.8018 31.345 -1733.6  35.861 3.691e-09 ***

# test which groups differ 
MX %>% emmeans(~ nreX*ni.mm) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp 

temp
# contrast      estimate     SE  df t.ratio p.value
# no 0 - yes 0  -0.00939 0.0282 588  -0.333  1.0000
# no 0 - no 1    0.23459 0.0332 588   7.060  <.0001
# no 0 - yes 1  -0.01392 0.0282 588  -0.493  1.0000
# yes 0 - no 1   0.24398 0.0282 588   8.641  <.0001
# yes 0 - yes 1 -0.00453 0.0221 588  -0.204  1.0000
# no 1 - yes 1  -0.24851 0.0282 588  -8.801  <.0001


# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]

#      contrast   estimate         SE  df   t.ratio      p.value
#2  no 0 - no 1  0.2345923 0.03323052 588  7.059545 1.890865e-11
#4 yes 0 - no 1  0.2439812 0.02823521 588  8.641027 2.633538e-16
#6 no 1 - yes 1 -0.2485088 0.02823521 588 -8.801381 9.017015e-17

# nreX reaction norm plot for supplementary figure S12
nreX.p = ggplot(data = emmip(MX, nreX ~ ni.mm, CIs=T, plotit=F), aes(x = ni.mm, y = yvar, group = nreX)) +
    geom_point() +
    geom_line(aes(linetype = nreX), size=1) +
    scale_linetype_manual(values=c("dashed", "solid"), labels=c("absent","present")) +
    geom_errorbar(aes(ymin = yvar - SE, ymax = yvar + SE), width=0.1)+
    ylab(expression(Fitness~(OD[600]))) + xlab("Nickel (mM)")+
    ylim(0.1, 0.95) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5), 
        axis.ticks.length = unit(0.1, "cm"),text = element_text(size=20), 
        axis.text.x = element_text(size=20, vjust=0.95), 
        axis.text.y = element_text(size=20,hjust=0.95),
        legend.key.size =  unit(0.4, "in"))


nreX.p

#ggsave("Fig_S12_E_updated_1.31.23.pdf",  nreX.p,  width = 5, height = 4, dpi = 1200)

# clade 1 only 
MX1=lm(OD ~ nreX*ni.mm, na.action=na.omit, data=Unique.ID.mean1) 

nreX.p1 = ggplot(data = emmip(MX1, nreX ~ ni.mm, CIs=T, plotit=F), aes(x = ni.mm, y = yvar, group = nreX)) +
    geom_point() +
    geom_line(aes(linetype = nreX), size=1) +
    scale_linetype_manual(values=c("dashed", "solid"), labels=c("absent","present")) +
    geom_errorbar(aes(ymin = yvar - SE, ymax = yvar + SE), width=0.1)+
    ylab(expression(Fitness~(OD[600]))) + xlab("Nickel (mM)")+
    ylim(0.1, 0.95) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5), 
        axis.ticks.length = unit(0.1, "cm"),text = element_text(size=20), 
        axis.text.x = element_text(size=20, vjust=0.95), 
        axis.text.y = element_text(size=20,hjust=0.95),
        legend.key.size =  unit(0.4, "in"))

nreX.p1

#ggsave("Fig_S12_F_updated_1.31.23.pdf",  nreX.p1,  width = 5, height = 4, dpi = 1200)

# clade 2
MX2=lm(OD ~ nreX*ni.mm, na.action=na.omit, data=Unique.ID.mean2) 

nreX.p2 = ggplot(data = emmip(MX2, nreX ~ ni.mm, CIs=T, plotit=F), aes(x = ni.mm, y = yvar, group = nreX)) +
    geom_point() +
    geom_line(aes(linetype = nreX), size=1) +
    scale_linetype_manual(values=c("dashed", "solid"), labels=c("absent","present")) +
    geom_errorbar(aes(ymin = yvar - SE, ymax = yvar + SE), width=0.1)+
    ylab(expression(Fitness~(OD[600]))) + xlab("Nickel (mM)")+
    ylim(0.1, 0.95) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5), 
        axis.ticks.length = unit(0.1, "cm"),text = element_text(size=20), 
        axis.text.x = element_text(size=20, vjust=0.95), 
        axis.text.y = element_text(size=20,hjust=0.95),
        legend.key.size =  unit(0.4, "in"))

nreX.p2

#ggsave("Fig_S12_G_updated_1.31.23.pdf",  nreX.p2,  width = 5, height = 4, dpi = 1200)

# other clades
MXo=lm(OD ~ nreX*ni.mm, na.action=na.omit, data=Unique.ID.meano) 

nreX.po = ggplot(data = emmip(MXo, nreX ~ ni.mm, CIs=T, plotit=F), aes(x = ni.mm, y = yvar, group = nreX)) +
    geom_point() +
    geom_line(aes(linetype = nreX), size=1) +
    scale_linetype_manual(values=c("dashed", "solid"), labels=c("absent","present")) +
    geom_errorbar(aes(ymin = yvar - SE, ymax = yvar + SE), width=0.1)+
    ylab(expression(Fitness~(OD[600]))) + xlab("Nickel (mM)")+
    ylim(0.1, 0.95) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5), 
        axis.ticks.length = unit(0.1, "cm"),text = element_text(size=20), 
        axis.text.x = element_text(size=20, vjust=0.95), 
        axis.text.y = element_text(size=20,hjust=0.95),
        legend.key.size =  unit(0.4, "in"))

nreX.po

#ggsave("Fig_S12_H_updated_1.31.23.pdf",  nreX.po,  width = 5, height = 4, dpi = 1200)

#### ***nreY ####
#model 
MY=lm(OD ~ nreY*ni.mm, na.action=na.omit, data=Unique.ID.mean) 

summary(MY)
#Residuals:
#     Min       1Q   Median       3Q      Max 
#-0.57536 -0.17212  0.01924  0.16149  0.64093 

#Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     0.674102   0.019743  34.145  < 2e-16 ***
#nreYyes        -0.002435   0.026687  -0.091    0.927    
#ni.mm1         -0.175915   0.027920  -6.301 5.81e-10 ***
#nreYyes:ni.mm1  0.195377   0.037740   5.177 3.10e-07 ***

#Residual standard error: 0.2285 on 588 degrees of freedom
#Multiple R-squared:  0.1006,	Adjusted R-squared:  0.09601 
#F-statistic: 21.92 on 3 and 588 DF,  p-value: 1.81e-13

# check residuals
simulationOutput = simulateResiduals(fittedModel = MY)
plot(simulationOutput)
# overall good fit, some heterogeneity of variance

# test fixed effects
drop1(MY, test = "F")
#Single term deletions
#Model:
#OD ~ nreY * ni.mm
#           Df Sum of Sq    RSS     AIC F value    Pr(>F)    
#<none>                  30.711 -1743.7                      
#nreY:ni.mm  1    1.3997 32.110 -1719.3    26.8 3.103e-07 ***

# test which groups differ from one 
MY %>% emmeans(~ nreY*ni.mm) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp 

temp
# contrast      estimate     SE  df t.ratio p.value
# no 0 - yes 0   0.00243 0.0267 588   0.091  1.0000
# no 0 - no 1    0.17591 0.0279 588   6.301  <.0001
# no 0 - yes 1  -0.01703 0.0267 588  -0.638  1.0000
# yes 0 - no 1   0.17348 0.0267 588   6.501  <.0001
# yes 0 - yes 1 -0.01946 0.0254 588  -0.766  1.0000
# no 1 - yes 1  -0.19294 0.0267 588  -7.230  <.0001


# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]

#      contrast   estimate         SE  df   t.ratio      p.value
#2  no 0 - no 1  0.1759147 0.02792021 588  6.300621 2.325822e-09
#4 yes 0 - no 1  0.1734802 0.02668653 588  6.500665 8.541652e-10
#6 no 1 - yes 1 -0.1929421 0.02668653 588 -7.229944 9.059346e-12

# nreY reaction norm plot for supplementary figure S12 
nreY.p = ggplot(data = emmip(MY, nreY ~ ni.mm, CIs=T, plotit=F), aes(x = ni.mm, y = yvar, group = nreY)) +
    geom_point() +
    geom_line(aes(linetype = nreY), size=1) +
    scale_linetype_manual(values=c("dashed", "solid"), labels=c("absent","present")) +
    geom_errorbar(aes(ymin = yvar - SE, ymax = yvar + SE), width=0.1)+
    ylab(expression(Fitness~(OD[600]))) + xlab("Nickel (mM)")+
    ylim(0.1, 0.95) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5), 
        axis.ticks.length = unit(0.1, "cm"),text = element_text(size=20), 
        axis.text.x = element_text(size=20, vjust=0.95), 
        axis.text.y = element_text(size=20,hjust=0.95),
        legend.key.size =  unit(0.4, "in"))
nreY.p

#ggsave("Fig_S12_I_updated_1.31.23.pdf",  nreY.p,  width = 5, height = 4, dpi = 1200)


# clade 1 only 
MY1=lm(OD ~ nreY*ni.mm, na.action=na.omit, data=Unique.ID.mean1) 

nreY.p1 = ggplot(data = emmip(MY1, nreY ~ ni.mm, CIs=T, plotit=F), aes(x = ni.mm, y = yvar, group = nreY)) +
    geom_point() +
    geom_line(aes(linetype = nreY), size=1) +
    scale_linetype_manual(values=c("dashed", "solid"), labels=c("absent","present")) +
    geom_errorbar(aes(ymin = yvar - SE, ymax = yvar + SE), width=0.1)+
    ylab(expression(Fitness~(OD[600]))) + xlab("Nickel (mM)")+
    ylim(0.1, 0.95) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5), 
        axis.ticks.length = unit(0.1, "cm"),text = element_text(size=20), 
        axis.text.x = element_text(size=20, vjust=0.95), 
        axis.text.y = element_text(size=20,hjust=0.95),
        legend.key.size =  unit(0.4, "in"))

nreY.p1

#ggsave("Fig_S12_J_updated_1.31.23.pdf",  nreY.p1,  width = 5, height = 4, dpi = 1200)


# clade 2
MY2=lm(OD ~ nreY*ni.mm, na.action=na.omit, data=Unique.ID.mean2) 

nreY.p2 = ggplot(data = emmip(MY2, nreY ~ ni.mm, CIs=T, plotit=F), aes(x = ni.mm, y = yvar, group = nreY)) +
    geom_point() +
    geom_line(aes(linetype = nreY), size=1) +
    scale_linetype_manual(values=c("dashed", "solid"), labels=c("absent","present")) +
    geom_errorbar(aes(ymin = yvar - SE, ymax = yvar + SE), width=0.1)+
    ylab(expression(Fitness~(OD[600]))) + xlab("Nickel (mM)")+
    ylim(0.1, 0.95) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5), 
        axis.ticks.length = unit(0.1, "cm"),text = element_text(size=20), 
        axis.text.x = element_text(size=20, vjust=0.95), 
        axis.text.y = element_text(size=20,hjust=0.95),
        legend.key.size =  unit(0.4, "in"))

nreY.p2

#ggsave("Fig_S12_K_updated_1.31.23.pdf",  nreY.p2,  width = 5, height = 4, dpi = 1200)


# other clades
MYo=lm(OD ~ nreY*ni.mm, na.action=na.omit, data=Unique.ID.meano) 

nreY.po = ggplot(data = emmip(MYo, nreY ~ ni.mm, CIs=T, plotit=F), aes(x = ni.mm, y = yvar, group = nreY)) +
    geom_point() +
    geom_line(aes(linetype = nreY), size=1) +
    scale_linetype_manual(values=c("dashed", "solid"), labels=c("absent","present")) +
    geom_errorbar(aes(ymin = yvar - SE, ymax = yvar + SE), width=0.1)+
    ylab(expression(Fitness~(OD[600]))) + xlab("Nickel (mM)")+
    ylim(0.1, 0.95) +
    theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5), 
        axis.ticks.length = unit(0.1, "cm"),text = element_text(size=20), 
        axis.text.x = element_text(size=20, vjust=0.95), 
        axis.text.y = element_text(size=20,hjust=0.95),
        legend.key.size =  unit(0.4, "in"))

nreY.po

#ggsave("Fig_S12_L_updated_1.31.23.pdf",  nreY.po,  width = 5, height = 4, dpi = 1200)



#### C089B tradeoff's analysis: wildtype & nreXY deletion ####

to=read.csv("C089B_wt_del_tradeoffs_data_KJ.csv")
# seperate metadata from data
meta.to = select(to, "Meta.Data")
to=select(to,-"Meta.Data")

str(to)
head(to)

#change to factors
to$Strain=as.factor(to$Strain)
to$Ni.mM=as.factor(to$Ni.mM)
to$Strain.ID=as.factor(to$Strain.ID)
to$Time.hr=as.factor(to$Time.hr)
to$nreXY=as.factor(to$nreXY)

# only keep 0 and 1mM nickel
to2=to[which(to$Ni.mM=="0" | to$Ni.mM=="1"),]

to2=droplevels(to2)

# subset by time
# 48 hours
to48= to2[which(to2$Time.hr=="48"),]
to48=droplevels(to48)
head(to48)
str(to48) #64 obs. of  8 variables

# Model 
M48=lm(OD ~ Strain*Ni.mM, na.action=na.omit, data=to48) 
summary(M48)

#Residuals:
#      Min        1Q    Median        3Q       Max 
#-0.043625 -0.008594  0.000625  0.007547  0.036812 
#
#Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
#(Intercept)      0.127063   0.004067  31.241  < 2e-16 ***
#Strain.IDWT        -0.005000   0.005752  -0.869    0.388    
#Ni.mM1          -0.062438   0.005752 -10.855 8.73e-16 ***
#Strain.IDWT:Ni.mM1  0.058562   0.008134   7.199 1.14e-09 ***

#Residual standard error: 0.01627 on 60 degrees of freedom
#Multiple R-squared:  0.7195,	Adjusted R-squared:  0.7055 
#F-statistic: 51.31 on 3 and 60 DF,  p-value: < 2.2e-16

# check residuals
simulationOutput = simulateResiduals(fittedModel = M48)
plot(simulationOutput)

# test fixed effects
drop1(M48, test = "F")
#Single term deletions
#Model:
#OD ~ Strain.ID * Ni.mM
#             Df Sum of Sq      RSS     AIC F value    Pr(>F)    
#<none>                    0.015880 -523.30                      
#Strain.ID:Ni.mM  1  0.013718 0.029599 -485.45  51.831 1.143e-09 ***

M48 %>% emmeans(~ Strain*Ni.mM) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp 

# contrast      estimate      SE df t.ratio p.value
# Del 0 - WT 0   0.00500 0.00575 60   0.869  0.7763
# Del 0 - Del 1  0.06244 0.00575 60  10.855  <.0001
# Del 0 - WT 1   0.00887 0.00575 60   1.543  0.3843
# WT 0 - Del 1   0.05744 0.00575 60   9.986  <.0001
# WT 0 - WT 1    0.00387 0.00575 60   0.674  0.7763
# Del 1 - WT 1  -0.05356 0.00575 60  -9.312  <.0001

# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]

#       contrast   estimate          SE df   t.ratio      p.value
#2 Del 0 - Del 1  0.0624375 0.005751868 60 10.855169 5.237467e-15
#4  WT 0 - Del 1  0.0574375 0.005751868 60  9.985886 1.127514e-13
#6  Del 1 - WT 1 -0.0535625 0.005751868 60 -9.312192 1.182256e-12

#### **Figure 2E ####
# Growth at 48 hours emmip plot
to48.plot = ggplot(data = emmip(M48, Strain ~ Ni.mM, CIs=T, plotit=F), aes(x = Ni.mM, y = yvar, group = Strain)) +
    geom_point() +
    geom_line(aes(linetype = Strain), size=1) +
    scale_linetype_manual(values=c("dashed", "solid"), labels = c((expression(Delta*"nreXY")), "Wild-type")) +
    geom_errorbar(aes(ymin = yvar - SE, ymax = yvar + SE), width=0.1)+
    ylab(expression(Fitness~(OD[600])))+ xlab("Nickel (mM)")+
    theme(panel.background = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),
        axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5), 
        axis.ticks.length = unit(0.1, "cm"),text = element_text(size=15, colour = "black"), axis.title = element_text(size = 20), 
        axis.text.x = element_text(size=15, vjust=0.95, colour = "black"), 
        axis.text.y = element_text(size=15,hjust=0.95, colour = "black"),
        legend.title = element_blank(), legend.text=element_text(size=20), legend.position="top", legend.key.size =  unit(0.5, "in"))

to48.plot

#ggsave("Tradeoff_C089Bwt.nreXYdel.pdf", to48.plot, width = 5, height = 5, dpi = 1200)


#### nre alleles analysis ####
# Make separate data frames for nreA, X or Y
# Note: each time this code is ran it will generate new data sets for nreA, X, or Y which may cause the figures to not appear as they do in the manuscript (because we randomly drop to one allele type if a strain has more than one copy)
# but we have a copy of the data subset used for figures in the manuscript that will need to be read in to regenerate figures exactly -- see below 

#make data subset (only those with "Genodata")
gegd2 = gegd[which(gegd$GenoData == "has data"),]
nrow(gegd2) #1048
length(unique(gegd2$Unique.ID)) #302

#### ***nreA alleles ####
# if a Unique.ID has more than 1 nreA allele copy, make new rows for those Unique.IDs
A.gegd = gegd2 %>%
    mutate(nreA_alleles = strsplit(as.character(nreA_alleles), ",")) %>%
    tidyr::unnest(nreA_alleles)

A.gegd = as.data.frame(A.gegd)
nrow(A.gegd) #1111

# change all nreA_alleles=NA (Unique.ID lacks nreA) to allele=0
A.gegd$nreA_alleles = as.character(A.gegd$nreA_alleles)
A.gegd$nreA_alleles[is.na(A.gegd$nreA_alleles)]= "0"
# change back to factor
A.gegd$nreA_alleles = as.factor(A.gegd$nreA_alleles)
head(A.gegd)

#before dropping to one allele, make sure I keep rows where there is mic data 
A.gegdc=A.gegd[complete.cases(A.gegd$Ni.MIC),]
nrow(A.gegdc) #1003
length(unique(A.gegdc$Unique.ID)) #296 Unique.IDs

#arbitrarily keep one nreA allele for each Unique.ID 
nreA.data = ddply(A.gegdc, .(Unique.ID), function(x) x[sample(nrow(x),1),])
nrow(nreA.data) #296

str(nreA.data)
nreA.data=droplevels(nreA.data)
levels(nreA.data$nreA_alleles)

#data for mixed model prep
#now use one nreA allele type per Unique.ID but use multiple reps of MIC in mixed model
#subset nreA.data 
nreA.data2= subset(A.gegdc, select = c(Unique.ID, Ni.MIC))

head(nreA.data2)
nrow(nreA.data2) #1003 rows

#now remove Ni MIC col from nreA.data, then merge nreA.data and nreA.data2 (each Unique.ID will have 1 nreA variant and 3 reps of MIC)
head(nreA.data)
nreA.data3 = subset(nreA.data, select = -Ni.MIC)

nreA=merge(nreA.data3,nreA.data2, by="Unique.ID")
nrow(nreA) #1003  
length(unique(nreA$Unique.ID)) #296

# Saved data set here 
#write.csv(x = nreA, file = paste("nreA_data_dropped_extra_alleles.csv"), row.names = F)
#If regenerating figure from manuscript, use commented line below in code to read in our version of nreA
#nreA = read.csv("nreA_data_dropped_extra_alleles.csv")
nreA$Unique.ID = as.factor(nreA$Unique.ID)
nreA$nreA_alleles = as.factor(nreA$nreA_alleles)

#model nreA alleles
M.nreA = lme4::lmer(Ni.MIC ~ nreA_alleles + (1|Unique.ID), na.action=na.omit, data = nreA)

# for mixed model
# check singularity
model = M.nreA
lme4::isSingular(model, tol = 1e-05)
# FALSE

relgrad = with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# e-08, convergence was fine

# check residuals
simulationOutput = DHARMa::simulateResiduals(fittedModel = M.nreA)
plot(simulationOutput)
# significant deviation for KS test, dispersion and outliers good
# some heterogeneity of variance


# test fixed effects
drop1(M.nreA, test = "Chisq")
#Model:
#Ni.MIC ~ nreA_alleles + (1 | Unique.ID)
#             npar    AIC    LRT   Pr(Chi)    
#<none>            2407.8                     
#nreA_alleles    6 2528.5 132.69 < 2.2e-16 ***


#emmeans 
M.nreA %>% emmeans::emmeans(~ nreA_alleles) %>%
  emmeans::contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp #holm adjustment for multiple comparisons


# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]
#   contrast   estimate        SE       df    t.ratio      p.value
#1     0 - 1 -2.1098417 0.4116100 281.9161  -5.125827 1.046803e-05
#2     0 - 2 -1.6060418 0.1538194 290.2725 -10.441085 1.414427e-20
#3     0 - 3 -1.3423313 0.1281627 291.4025 -10.473654 1.126722e-20
#5     0 - 5 -0.8809763 0.2487926 278.5212  -3.541007 7.940112e-03
#11    1 - 6  1.5486009 0.4526814 280.8302   3.420951 1.119398e-02
#15    2 - 6  1.0448010 0.2432244 281.4251   4.295626 4.322491e-04
#18    3 - 6  0.7810906 0.2278670 280.5612   3.427836 1.119398e-02


#make a plot with the emmeans estimates and SE's 
nrea.emmeans = as.data.frame(M.nreA %>% emmeans::emmeans(~ nreA_alleles))

# recall note above: figure may not appear as it does in the manuscript if an allele type was dropped while re-running code
#if you want exact same figure you must read in the nreA file
#emmeans plot
nreAplot = ggplot(nrea.emmeans, aes(x=reorder(nreA_alleles,+emmean), y=emmean, fill=nreA_alleles)) +
  geom_bar(stat='identity', color="black") + 
  geom_errorbar(aes(ymin=emmean-SE, ymax = emmean+SE), width=0.15) + 
  xlab('nreA allele') + ylab('MIC') +
  scale_y_continuous(breaks=seq(0, 5.8, 1)) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),text = element_text(size=30), axis.text.x = element_text(size=28, vjust=0.5), 
        axis.text.y = element_text(size=28,hjust=0.95), plot.title = element_text(hjust = 0.5), legend.position = "none") +
  scale_fill_manual(values=c("1" = "#4A1486",
                             "2" = "#6A51A3" ,
                             "3"  ="#807DBA" ,
                             "4" = "#9E9AC8",
                             "5" ="#BCBDDC" ,
                             "6" = "#DADAEB",
                             "0"= "white")) +
  annotate("text", x = 1, y = 3.7, label = "a", size = 5) +
  annotate("text", x = 2, y = 4.2, label = "ab", size = 5) +
  annotate("text", x = 3, y = 4.6, label = "bc", size = 5) +
  annotate("text", x = 4, y = 5.2, label = "abc", size = 5) +
  annotate("text", x = 5, y = 5.2, label = "c", size = 5) +
  annotate("text", x = 6, y = 5.3, label = "c", size = 5) +
  annotate("text", x = 7, y = 5.7, label = "c", size = 5) 

#### **Figure 3B ####
nreAplot
#ggsave("Fig3b_updated_1.31.2023.pdf",  nreAplot,  width = 3.5, height = 3.5, dpi = 1200)



#### ***nreX alleles ####

# if Unique.ID has more than 1 type of nreX allele, change so that there are new rows for those Unique.IDs
X.gegd=gegd2 %>% 
    mutate(nreX_alleles = strsplit(as.character(nreX_alleles), ",")) %>% 
    tidyr::unnest(nreX_alleles)

X.gegd = as.data.frame(X.gegd)
nrow(X.gegd) #1096

# change all nreX_alleles=NA (Unique.ID lacks nreX) to allele=0
X.gegd$nreX_alleles = as.character(X.gegd$nreX_alleles)
X.gegd$nreX_alleles[is.na(X.gegd$nreX_alleles)]= "0"
#change back to factor
X.gegd$nreX_alleles = as.factor(X.gegd$nreX_alleles)
head(X.gegd)

#before dropping to one allele, make sure I keep rows where there is mic data 
X.gegdc=X.gegd[complete.cases(X.gegd$Ni.MIC),]
nrow(X.gegdc) #988
length(unique(X.gegdc$Unique.ID)) #296 Unique.IDs


#arbitrarily keep one nreX allele for each Unique.ID 
nreX.data = ddply(X.gegdc, .(Unique.ID), function(x) x[sample(nrow(x),1),])

nrow(nreX.data) #296
nreX.data = droplevels(nreX.data)
levels(nreX.data$nreX_alleles)

#data for mixed model prep
#now use one nreX allele type per Unique.ID but use multiple reps of MIC in mixed model
#subset nreX.data 
nreX.data2= subset(X.gegdc, select = c(Unique.ID, Ni.MIC))

head(nreX.data2)
nrow(nreX.data2) #988 rows

#now remove Ni MIC col from nreX.data, then merge nreX.data and nreX.data2 (each Unique.ID will have 1 nreX variant and 3 reps of MIC)
head(nreX.data)
nreX.data3 = subset(nreX.data, select = -Ni.MIC)

nreX=merge(nreX.data3,nreX.data2, by="Unique.ID")
nrow(nreX) #988  
length(unique(nreX$Unique.ID)) #296

# Saved data set here 
#write.csv(x = nreX, file = paste("nreX_data_dropped_extra_alleles.csv"), row.names = F)
#If regenerating figure from manuscript, use commented line below in code to read in our version of nreX
#nreX = read.csv("nreX_data_dropped_extra_alleles.csv")
nreX$Unique.ID = as.factor(nreX$Unique.ID)
nreX$nreX_alleles = as.factor(nreX$nreX_alleles)

#model nreX alleles
M.nreX = lme4::lmer(Ni.MIC ~ nreX_alleles + (1|Unique.ID), na.action=na.omit, data = nreX)

# for mixed model
# check singularity
model = M.nreX
lme4::isSingular(model, tol = 1e-05)
# FALSE

relgrad = with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# e-09, convergence was fine

# check residuals
simulationOutput = DHARMa::simulateResiduals(fittedModel = M.nreX)
plot(simulationOutput)
# significant deviation for KS test, dispersion and outliers good
# some heterogeneity of variance

# test fixed effects
drop1(M.nreX, test = "Chisq")
#Model:
#Ni.MIC ~ nreX_alleles + (1 | Unique.ID)
#             npar    AIC   LRT   Pr(Chi)    
#<none>            2402.9                     
#nreX_alleles    3 2511.4 114.56 < 2.2e-16 ***


#emmeans 
M.nreX %>% emmeans::emmeans(~ nreX_alleles) %>%
  emmeans::contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp #holm adjustment for multiple comparisons


# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]
#  contrast   estimate        SE       df    t.ratio      p.value
#1    0 - 1 -1.4068554 0.1199335 293.7846 -11.730293 1.527278e-25
#3    0 - 3 -0.9268288 0.1871349 288.1381  -4.952730 6.243741e-06
#5    1 - 3  0.4800266 0.1751869 286.5186   2.740083 2.611172e-02



#make a plot with the emmeans estimates and SE's
nrex.emmeans = as.data.frame(M.nreX %>% emmeans::emmeans(~ nreX_alleles))

# recall note above: figure may not appear as it does in the manuscript if an allele type was dropped while re-running code
#if you want exact same figure you must read in the nreX file
#emmeans plot
nreXplot = ggplot(nrex.emmeans, aes(x=reorder(nreX_alleles,+emmean), y=emmean, fill=nreX_alleles)) +
  geom_bar(stat='identity', color="black") + 
  geom_errorbar(aes(ymin=emmean-SE, ymax = emmean+SE), width=0.15) + 
  xlab('nreX allele') + ylab('MIC') +
  scale_y_continuous(breaks=seq(0, 5.8, 1)) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),text = element_text(size=30), axis.text.x = element_text(size=28, vjust=0.5), 
        axis.text.y = element_text(size=28,hjust=0.95), plot.title = element_text(hjust = 0.5), legend.position = "none") +
  scale_fill_manual(values=c("1" = "#8C510A",
                             "2" = "#BF812D" ,
                             "3"  ="#DFC27D" ,
                             "0"= "white")) +
  annotate("text", x = 1, y = 3.6, label = "a", size = 5) +
  annotate("text", x = 2, y = 4.5, label = "b", size = 5) +
  annotate("text", x = 3, y = 5.2, label = "abc", size = 5) +
  annotate("text", x = 4, y = 5.2, label = "c", size = 5)

#### **Figure 3C ####
nreXplot
#ggsave("Fig3c_updated_1.31.2023.pdf",  nreXplot,  width = 2.8, height = 3.5, dpi = 1200) 





#### ***nreY alleles ####
# if Unique.ID has more than 1 type of nreY allele, change so that there are new rows for those Unique.IDs
Y.gegd=gegd2 %>% 
    mutate(nreY_alleles = strsplit(as.character(nreY_alleles), ",")) %>% 
    tidyr::unnest(nreY_alleles)

Y.gegd = as.data.frame(Y.gegd)
nrow(Y.gegd) #1072

# change all nreY_alleles=NA (Unique.ID lacks nreY) to allele=0
Y.gegd$nreY_alleles = as.character(Y.gegd$nreY_alleles)
Y.gegd$nreY_alleles[is.na(Y.gegd$nreY_alleles)]= "0"
#change back to factor
Y.gegd$nreY_alleles = as.factor(Y.gegd$nreY_alleles)
head(Y.gegd)

#only 2 Unique.IDs have nreY alleles 7 and 8 
#drop nreY_alleles 7 and 8
Y.gegd = Y.gegd[!(Y.gegd$nreY_allele=="7" | Y.gegd$nreY_allele=="8"),]


#before dropping to one allele, make sure I keep rows where there is mic data 
Y.gegdc=Y.gegd[complete.cases(Y.gegd$Ni.MIC),]
nrow(Y.gegdc) #952
length(unique(Y.gegdc$Unique.ID)) #292 Unique.IDs


#arbitrarily keep one nreY allele for each Unique.ID 
nreY.data = ddply(Y.gegdc, .(Unique.ID), function(x) x[sample(nrow(x),1),])
 
nrow(nreY.data) #292
nreY.data = droplevels(nreY.data)
levels(nreY.data$nreY_alleles)


#data for mixed model prep
#now use one nreY allele type per Unique.ID but use multiple reps of MIC in mixed model
#subset nreY.data 
nreY.data2= subset(Y.gegdc, select = c(Unique.ID, Ni.MIC))

head(nreY.data2)
nrow(nreY.data2) #952 rows

#now remove Ni MIC col from nreY.data, then merge nreY.data and nreY.data2 (each Unique.ID will have 1 nreY variant and 3 reps of MIC)
head(nreY.data)
nreY.data3 = subset(nreY.data, select = -Ni.MIC)

nreY=merge(nreY.data3,nreY.data2, by="Unique.ID")
nrow(nreY) #952  
length(unique(nreY$Unique.ID)) #292

# Saved data set here 
#write.csv(x = nreY, file = paste("nreY_data_dropped_extra_alleles.csv"), row.names = F)
#If regenerating figure from manuscript, use commented line below in code to read in our version of nreY
#nreY = read.csv("nreY_data_dropped_extra_alleles.csv")
nreY$Unique.ID = as.factor(nreY$Unique.ID)
nreY$nreY_alleles = as.factor(nreY$nreY_alleles)

#model nreY alleles 
M.nreY = lme4::lmer(Ni.MIC ~ nreY_alleles + (1|Unique.ID), na.action=na.omit, data = nreY)

# for mixed model
# check singularity
model = M.nreY
lme4::isSingular(model, tol = 1e-05)
# FALSE

relgrad = with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# e-07, convergence was fine

# check residuals
simulationOutput = DHARMa::simulateResiduals(fittedModel = M.nreY)
plot(simulationOutput)
# significant deviation for KS test
# some heterogeneity of variance, outlier

# test fixed effects
drop1(M.nreY, test = "Chisq")
#Model:
#Ni.MIC ~ nreY_alleles + (1 | Unique.ID)
#             npar    AIC    LRT   Pr(Chi)    
#<none>            2370.4                     
#nreY_alleles    6 2459.2 100.78 < 2.2e-16 ***


#emmeans 
M.nreY %>% emmeans::emmeans(~ nreY_alleles) %>%
  emmeans::contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp #holm adjustment for multiple comparisons


# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]
#   contrast  estimate        SE       df   t.ratio      p.value
#1     0 - 1 -1.5023073 0.1889833 284.0571 -7.949420 9.339561e-13
#2     0 - 2 -1.2565746 0.1970659 284.1694 -6.376418 1.389923e-08
#3     0 - 3 -1.0808409 0.1412334 285.0998 -7.652869 6.106190e-12
#11    1 - 6  1.3457234 0.2861590 283.1935  4.702712 7.229101e-05
#15    2 - 6  1.0999907 0.2915601 283.2762  3.772775 3.342058e-03
#18    3 - 6  0.9242569 0.2571365 283.3005  3.594421 6.133349e-03



#make a plot with the emmeans estimates and SE's
nrey.emmeans = as.data.frame(M.nreY %>% emmeans::emmeans(~ nreY_alleles))

# recall note above: figure may not appear as it does in the manuscript if an allele type was dropped while re-running code
#if you want exact same figure you must read in the nreY file
#emmeans plot
nreYplot = ggplot(nrey.emmeans, aes(x=reorder(nreY_alleles,+emmean), y=emmean, fill=nreY_alleles)) +
  geom_bar(stat='identity', color="black") + 
  geom_errorbar(aes(ymin=emmean-SE, ymax = emmean+SE), width=0.15) + 
  xlab('nreY allele') + ylab('MIC') +
  scale_y_continuous(breaks=seq(0, 5.8, 1)) +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank(),axis.line = element_line(colour = "black"),axis.ticks = element_line(colour = "black", size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),text = element_text(size=30), axis.text.x = element_text(size=28, vjust=0.5), 
        axis.text.y = element_text(size=28,hjust=0.95), plot.title = element_text(hjust = 0.5), legend.position = "none") +
scale_fill_manual(values=c("1" = "#08306B",
                             "2" = "#08519C" ,
                             "3"  ="#2171B5" ,
                             "4" = "#4292C6",
                             "5" ="#6BAED6" ,
                             "6" = "#9ECAE1",
                             "0"= "white")) +
  annotate("text", x = 1, y = 3.9, label = "a", size = 5) +
  annotate("text", x = 2, y = 4.2, label = "ab", size = 5) +
  annotate("text", x = 3, y = 4.6, label = "abc", size = 5) +
  annotate("text", x = 4, y = 5.2, label = "abc", size = 5) +
  annotate("text", x = 5, y = 5.2, label = "c", size = 5) +
  annotate("text", x = 6, y = 5.3, label = "c", size = 5) +
  annotate("text", x = 7, y = 5.5, label = "c", size = 5) 

#### **Figure 3D ####
nreYplot

#ggsave("Fig3d_updated_1.31.2023.pdf", nreYplot, width = 3.5, height = 3.5, dpi = 1200) 



# *****nreAXY allele iteration models ####
# Does randomly dropping alleles change the results? --overall no
# run this code to generate a new 100 iterations or use csv files:
# nreA_table_randomized.csv
# nreX_table_randomized.csv
# nreY_table_randomized.csv

# ******nreA iterations ####
#A.gegdC is the data set for nreA
#each Unique.ID listed multiple times (may have only 1 allele type or a few Unique.IDs have more then 1) and this data only has rows with MIC data

nreAtable = data.frame()
# Run iterations on 
for (i in 1:100) {
  
  # create data set with only one observation per Unique.ID 
  A.gegdc.rand = ddply(A.gegdc, .(Unique.ID), function(x) x[sample(nrow(x),1),])
  
  #now remove Ni MIC col from nreA.data, then merge nreA.data and nreA.data2 (each Unique.ID will have 1 nreA variant and 3 reps of MIC)
  A.gegdc.rand=subset(A.gegdc.rand, select =-Ni.MIC)

  nreA.rand=merge(A.gegdc.rand,nreA.data2, by="Unique.ID", all.x=T, all.y=T)

  # run model on randomized dataset
  MA.rand = lmer(Ni.MIC ~ nreA_alleles + (1|Unique.ID), na.action=na.omit, data = nreA.rand)

  #emmeans 
  temp = MA.rand %>% emmeans(~ nreA_alleles) %>%
    contrast(., method = "pairwise") %>%
    summary(by = NULL, adjust = "holm")
  
  temp2=as.data.frame(temp)
  nreAtable = rbind(temp2, nreAtable)

}

head(nreAtable)
nrow(nreAtable) 

#change dash in contrasts to underscore
nreAtable$contrast = gsub('-', '_', nreAtable$contrast)

#write.csv(x = nreAtable, file = paste("nreA_table_randomized.csv"), row.names = F)
str(nreAtable)
nreAtable$contrast = as.factor(nreAtable$contrast)

# for manuscript significant results, check results of iterations
#first list percent of time the contrast was present out of 100 iterations
#then list percent of time the contrast was significant out of number of times it was present
#contrast 0_1
sum(nreAtable$contrast == "0 _ 1") #100% present
length(intersect(which(nreAtable['contrast']== '0 _ 1'), which(nreAtable['p.value']<=0.05))) #100% sig
#contrast 0_2
sum(nreAtable$contrast == "0 _ 2") #100% present
length(intersect(which(nreAtable['contrast']== '0 _ 2'), which(nreAtable['p.value']<=0.05))) #100% sig
#contrast 0_3
sum(nreAtable$contrast == "0 _ 3") #100% present
length(intersect(which(nreAtable['contrast']== '0 _ 3'), which(nreAtable['p.value']<=0.05))) #100% sig
#contrast 0_5
sum(nreAtable$contrast == "0 _ 5") #100% present
length(intersect(which(nreAtable['contrast']== '0 _ 5'), which(nreAtable['p.value']<=0.05))) #95% sig
#contrast 1_6
sum(nreAtable$contrast == "1 _ 6") #100% present
length(intersect(which(nreAtable['contrast']== '1 _ 6'), which(nreAtable['p.value']<=0.05))) #66% sig
#contrast 2_6
sum(nreAtable$contrast == "2 _ 6") #100% present
length(intersect(which(nreAtable['contrast']== '2 _ 6'), which(nreAtable['p.value']<=0.05))) #100% sig
#contrast 3_6
sum(nreAtable$contrast == "3 _ 6") #100% present
length(intersect(which(nreAtable['contrast']== '3 _ 6'), which(nreAtable['p.value']<=0.05))) #56% sig



# ******nreX iterations ####

nreXtable = data.frame()
# Run iterations on 
for (i in 1:100) {
  
  # create data set with only one observation per Unique.ID
  X.gegdc.rand = ddply(X.gegdc, .(Unique.ID), function(x) x[sample(nrow(x),1),])
  
  #now remove Ni MIC col from nreX.data, then merge nreX.data and nreX.data2 (each Unique.ID will have 1 nreX variant and 3 reps of MIC)
  X.gegdc.rand=subset(X.gegdc.rand, select =-Ni.MIC)

  nreX.rand=merge(X.gegdc.rand,nreX.data2, by="Unique.ID", all.x=T, all.y=T)

  # run model on randomized dataset
  MX.rand = lmer(Ni.MIC ~ nreX_alleles + (1|Unique.ID), na.action=na.omit, data = nreX.rand)

  #emmeans 
  temp = MX.rand %>% emmeans(~ nreX_alleles) %>%
    contrast(., method = "pairwise") %>%
    summary(by = NULL, adjust = "holm")
  
  temp2=as.data.frame(temp)
  nreXtable = rbind(temp2, nreXtable)

}

head(nreXtable)
nrow(nreXtable)

#change dash in contrasts to underscore
nreXtable$contrast = gsub('-', '_', nreXtable$contrast)

#write.csv(x = nreXtable, file = paste("nreX_table_randomized.csv"), row.names = F)

# for manuscript significant results, check results of iterations
#first list percent of time the contrast was present out of 100 iterations
#then list percent of time the contrast was significant out of number of times it was present
#contrast 0_1
sum(nreXtable$contrast == "0 _ 1") #100% present
length(intersect(which(nreXtable['contrast']== '0 _ 1'), which(nreXtable['p.value']<=0.05))) #100% sig
#contrast 0_3
sum(nreXtable$contrast == "0 _ 3") #100% present
length(intersect(which(nreXtable['contrast']== '0 _ 3'), which(nreXtable['p.value']<=0.05))) #100% sig
#contrast 1_3
sum(nreXtable$contrast == "1 _ 3") #100% present
length(intersect(which(nreXtable['contrast']== '1 _ 3'), which(nreXtable['p.value']<=0.05))) #87% sig



# ******nreY iterations ####
nreYtable = data.frame()
# Run iterations on 
for (i in 1:100) {
  
  # create data set with only one observation per Unique.ID
  Y.gegdc.rand = ddply(Y.gegdc, .(Unique.ID), function(x) x[sample(nrow(x),1),])
  
  #now remove Ni MIC col from nreY.data, then merge nreY.data and nreY.data2 (each Unique.ID will have 1 nreY variant and 3 reps of MIC)
  Y.gegdc.rand=subset(Y.gegdc.rand, select =-Ni.MIC)

  nreY.rand=merge(Y.gegdc.rand,nreY.data2, by="Unique.ID", all.x=T, all.y=T)

  # run model on randomized dataset
  MY.rand = lmer(Ni.MIC ~ nreY_alleles + (1|Unique.ID), na.action=na.omit, data = nreY.rand)

  #emmeans 
  temp = MY.rand %>% emmeans(~ nreY_alleles) %>%
    contrast(., method = "pairwise") %>%
    summary(by = NULL, adjust = "holm")
  
  temp2=as.data.frame(temp)
  nreYtable = rbind(temp2, nreYtable)

}

head(nreYtable)
nrow(nreYtable)

#change dash in contrasts to underscore
nreYtable$contrast = gsub('-', '_', nreYtable$contrast)

#write.csv(x = nreYtable, file = paste("nreY_table_randomized.csv"), row.names = F)

# for manuscript significant results, check results of iterations
#first list percent of time the contrast was present out of 100 iterations
#then list percent of time the contrast was significant out of number of times it was present
#contrast 0_1
sum(nreYtable$contrast == "0 _ 1") #100% present
length(intersect(which(nreYtable['contrast']== '0 _ 1'), which(nreYtable['p.value']<=0.05))) #100% sig
#contrast 0_2
sum(nreYtable$contrast == "0 _ 2") #100% present
length(intersect(which(nreYtable['contrast']== '0 _ 2'), which(nreYtable['p.value']<=0.05))) #100% sig
#contrast 0_3
sum(nreYtable$contrast == "0 _ 3") #100% present
length(intersect(which(nreYtable['contrast']== '0 _ 3'), which(nreYtable['p.value']<=0.05))) #100% sig
#contrast 1_6
sum(nreYtable$contrast == "1 _ 6") #100% present
length(intersect(which(nreYtable['contrast']== '1 _ 6'), which(nreYtable['p.value']<=0.05))) #100% sig
#contrast 2_6
sum(nreYtable$contrast == "2 _ 6") #100% present
length(intersect(which(nreYtable['contrast']== '2 _ 6'), which(nreYtable['p.value']<=0.05))) #99% sig
#contrast 3_6
sum(nreYtable$contrast == "3 _ 6") #100% present
length(intersect(which(nreYtable['contrast']== '3 _ 6'), which(nreYtable['p.value']<=0.05))) #73% sig


