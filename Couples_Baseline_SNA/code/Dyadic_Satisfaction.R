library(lme4)
library(lavaan)
library(foreign)
library(psych)
library(plyr)
library(ggplot2)
library(afex)
library(reshape2)
setwd("/Users/michael/Tresors/PD_SNA/Couples_Baseline_SNA")
source(file.path(getMainDir(), "Miscellaneous", "Global_Functions.R"))
source("code/APIM_Functions.R")
load(file="data/SNA_Processed_27Jul2017.RData")
load(file="data/selfreports/couples_baseline_clinical_27Jul2017.RData")

#look at relationship satisfaction as a function of BPD etc.
couples_baseline_clin <- subset(couples_baseline_clin, PTNUM %in% unique(sna$PTNUM))

corstarsl(couples_baseline_clin[,c("DASCon", "DASSat", "DASCoh", "DASAffExp", "DASTotal", "CTS_PhysAssV", "CTS_PhysAssP", "CTS_PsychAggV", "CTS_PsychAggP", "CTS_Perp", "CTS_Victim")])

corstarsl(couples_baseline_clin[,c("DASCon", "DASSat", "DASCoh", "DASAffExp", "DASTotal", "IIP_agency", "IIP_communion", "IIP_elevation", "IIP_pd1", "IIP_pd2", "IIP_pd3")])

predstocenter <- c("PAIBORtot", "IIP_agency", "IIP_communion", "IIP_elevation", "IIP_pd1", "IIP_pd2", "IIP_pd3", "IIP_c1", "IIP_c2", "ECRanx", "ECRavoid", "DASTotal", 
    "CTS_Victim", "CTS_Perp", "CTS_PsychAggV", "CTS_PsychAggP", "CTS_PhysAssV", "CTS_PhysAssP", "bordl_sidp", "narci_sidp", "antso_sidp", "avoid_sidp", "OPD_sidp", "depen_sidp", "nobpd", "nonarc")

couples_baseline_clin <- f_centerPredictors(couples_baseline_clin, predstocenter, addsuffix=".c")
couples_clin_wide <- f_centerPredictors(couples_clin_wide, as.vector(outer(predstocenter, c("0", "1"), paste, sep="_")), addsuffix=".c")


#free <- indistinguishable <- " %WITHIN%\n "
#for (p in 1:length(predictors)) {
#  free <- paste0(free, DV, partners[1], " ON ", predictors[p], partners[1], " ", predictors[p], partners[2], ";\n ",
#      DV, partners[2], " ON ", predictors[p], partners[1], " ", predictors[p], partners[2], ";\n")
#  indistinguishable <- paste0(indistinguishable, DV, partners[1], " ON ", predictors[p], partners[1], " (a", p, ")\n ", predictors[p], partners[2], " (p", p, ");\n ",
#      DV, partners[2], " ON ", predictors[p], partners[1], " (p", p, ")\n ", predictors[p], partners[2], " (a", p, ");\n")
#}
#free <- paste0(free, " %BETWEEN%\n")
#indistinguishable <- paste0(indistinguishable, " %BETWEEN%\n")





xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor="IIP_pd1", additional=".c", printall=FALSE) #interpersonal sensitivity
xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor="IIP_pd2", additional=".c", printall=FALSE) #interpersonal ambivalence
xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor="IIP_pd3", additional=".c") #aggression

#not much on c scales
xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor="IIP_c1", additional=".c") #need for social approval
xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor="IIP_c2", additional=".c") #lack of sociability


xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor=c("IIP_pd1", "IIP_pd2", "IIP_pd3"), additional=".c", printall=FALSE) #all PD scales

#conclusion: aggression, but not sensitivity or ambivalence was most strongly associated with relationship dissatisfaction,
#both in self and other.


#switch to circle

xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor="IIP_elevation", additional=".c", printall=FALSE)
xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor="IIP_agency", additional=".c", printall=TRUE)
xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor="IIP_communion", additional=".c", printall=FALSE)

xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor=c("IIP_communion", "IIP_agency", "IIP_elevation"), additional=".c", printall=FALSE)

xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor="narci_sidp", additional=".c")
xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor=c("narci_sidp", "nonarc"), additional=".c", printall=FALSE) #holds up
anova(xx$indistinguishable, xx$afree)
xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor="bordl_sidp", additional=".c", printall=FALSE)
xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor=c("bordl_sidp", "nobpd"), additional=".c", printall=TRUE)
xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor="bordlCount", printall=FALSE)
xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor=c("nobpdCount", "bordlCount"), printall=FALSE)
xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor=c("pdtot"), printall=FALSE)
xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor=c("allpdCount"), printall=FALSE)
xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor="antso_sidp", additional=".c")
#xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor="Victim", additional=".c")
#xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor="Perp", additional=".c", printall=TRUE) #patient perp contributes to greater patient and partner dissatisfaction (but converse does not hold
#xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor="ECRanx", additional=".c")
#xx <- runAPIM(couples_clin_wide, DV="DASTotal", predictor="ECRavoid", additional=".c")

#wide
#aidan <- read.csv("data/Couples Individual Level Data Crossed Mplus.csv")
#aidan <- aidan[,c("")]

#aidan <- read.csv("data/Couples Individual Level Data.csv")
#aidan_melt <- melt(aidan[,sapply(aidan, is.numeric)], id.vars=c("ptnum", "dyadid"))
#aidan_wide <- dcast(aidan_melt, ptnum ~ variable + dyadid)
#
##this matches Aidan's Mplus output
#xx <- runAPIM(aidan_wide, DV="dastotal", predictor="DOM1", additional="", printall=TRUE)

#range(aidan_wide$dastotal_0, na.rm=T)
#range(aidan_wide$dastotal_1, na.rm=T)

#mycomp <- subset(couples_clin_wide, PTNUM %in% aidan_wide$ptnum)

#combo <- merge(mycomp[,c("PTNUM", "DASTotal_0", "DASTotal_1", "iip_agency_0", "iip_agency_1")], 
#    aidan_wide[,c("ptnum", "dastotal_0", "dastotal_1", "DOM1_0", "DOM1_1")], by.x="PTNUM", by.y="ptnum", all.y=TRUE)
#
#with(combo, all.equal(DASTotal_0, dastotal_0))
#with(combo, DASTotal_0 - dastotal_0)
#with(combo, DASTotal_1 - dastotal_1)
#with(combo, cor(iip_agency_0, DOM1_0, use="pairwise.complete.obs"))
#with(combo, cor(iip_agency_1, DOM1_1, use="pairwise.complete.obs"))

#library(MplusAutomation)
#test <- mplusObject(MODEL = "DASTotal_0.c ON iip_agency_0.c iip_agency_1.c; DASTotal_1.c ON iip_agency_0.c iip_agency_1.c;", rdata = couples_clin_wide,
#    usevariables=c("DASTotal_0.c", "DASTotal_1.c", "iip_agency_0.c", "iip_agency_1.c", "PTNUM"))
#
## estimate the model in Mplus and read results back into R
#res <- mplusModeler(test, "couples.dat", modelout = "model1.inp", run = 1L)


#mean diffs in agency or communion?
with(couples_clin_wide, t.test(iip_agency_0, iip_agency_1))
with(couples_clin_wide, t.test(iip_communion_0, iip_communion_1))
with(couples_clin_wide, t.test(iip_elevation_0, iip_elevation_1))

#descriptives
library(effsize)

vars <- c("DASTotal", "iip_elevation", "iip_communion", "iip_agency", "allpdCount", "bordlCount", "narciCount", "nobpdCount", "nonarcCount")
for (v in vars) {
  cat("\n\n-----\n", v, "\n\n----\n")
  print(cohen.d(couples_baseline_clin[[v]], factor(couples_baseline_clin$DyadID), pooled=TRUE, na.rm=TRUE))
  print(tapply(couples_baseline_clin[[v]], couples_baseline_clin$DyadID, mean, na.rm=TRUE))
  print(tapply(couples_baseline_clin[[v]], couples_baseline_clin$DyadID, sd, na.rm=TRUE))
}

#combine attachment style with PDs
#with all these predictors, even though VIF is acceptable, they really knock each other out
#what about the possibility that attachment style effects are exacerbated by BPD?
combined <- '
DASTotal_0.c ~ a1*IIP_PD3_0.c + p1*IIP_PD3_1.c + 
   a2*ECRavoid_0.c + p2*ECRavoid_1.c + 
   a3*bordl_sidp_0.c + p3*bordl_sidp_1.c +
   a4*ECRanx_0.c + a4*ECRanx_1.c +
   a5*narci_sidp_0.c + p5*narci_sidp_1.c + p_age_0 + p_age_1

DASTotal_1.c ~ a1*IIP_PD3_1.c + p1*IIP_PD3_0.c + 
   a2*ECRavoid_1.c + p2*ECRavoid_0.c +
   a3*bordl_sidp_1.c + p3*bordl_sidp_0.c +
   a4*ECRanx_1.c + p4*ECRanx_0.c +
	 a5*narci_sidp_1.c + p5*narci_sidp_0.c + p_age_1 + p_age_0
'

combined <- '
    DASTotal_0.c ~ iip_agency_0.c + iip_agency_1.c     
    DASTotal_1.c ~ iip_agency_1.c + iip_agency_0.c
    '

m <- sem(combined, couples_clin_wide, missing="ML", estimator="MLR")
summary(m, fit.measures=TRUE)
#modindices(m)

couples_clin_wide$dummy <- rnorm(nrow(couples_clin_wide))
dummy_lm <- lm(dummy ~ IIP_PD3_0.c + IIP_PD3_1.c + ECRavoid_0.c + ECRavoid_1.c + ECRanx_0.c + ECRanx_1.c + bordl_sidp_0.c + bordl_sidp_1.c + narci_sidp_0.c + narci_sidp_1.c, couples_clin_wide)
car::vif(dummy_lm)

cor(couples_clin_wide[,c("IIP_PD3_0.c", "IIP_PD3_1.c", "ECRavoid_0.c", "ECRavoid_1.c", "ECRanx_0.c", "ECRanx_1.c", "bordl_sidp_0.c", "bordl_sidp_1.c", "narci_sidp_0.c", "narci_sidp_1.c")], use="pairwise.complete.obs")

mixed(DASTotal ~ IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + Victim.c + narci_sidp.c + bordl_sidp.c + antso_sidp.c + (1 | PTNUM), couples_baseline_clin, method="PB")



##ANALYSES OF DAS SUBSCALES
mixed(DASTotal ~ IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + Victim.c + narci_sidp.c + bordl_sidp.c + antso_sidp.c + (1 | PTNUM), couples_baseline_clin, method="PB")

summary(m <- lmer(DASCon ~ IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | PTNUM), couples_baseline_clin))
summary(m <- lmer(DASCon ~ IIP_PD1.c*IIP_PD2.c*IIP_PD3.c + (1 | PTNUM), couples_baseline_clin))
summary(m <- lmer(DASCon ~ PAIBORtot.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | PTNUM), couples_baseline_clin))
summary(m <- lmer(DASCon ~ PAIBORtot.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + Victim.c + Perp.c + (1 | PTNUM), couples_baseline_clin))

summary(m <- lmer(DASCon ~ PAIBORtot.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + PsychAggV.c + PsychAggP.c + (1 | PTNUM), couples_baseline_clin))
summary(m <- lmer(DASCon ~ PAIBORtot.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + PhysAssV.c + PhysAssP.c + (1 | PTNUM), couples_baseline_clin))

#narcissistic Sx, problems with aggression (IIP PD3), and being the victim of psychological and physical aggression/assault predict lower consensus
summary(m <- lmer(DASCon ~ narci_sidp.c + bordl_sidp.c + antso_sidp.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | PTNUM), couples_baseline_clin))

summary(m <- lmer(DASCon ~ narci_sidp.c + bordl_sidp.c + antso_sidp.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + Victim.c + Perp.c + (1 | PTNUM), couples_baseline_clin))
mixed(DASCon ~ narci_sidp.c + bordl_sidp.c + antso_sidp.c + IIP_PD3.c + (1 | PTNUM), couples_baseline_clin)

summary(m <- lmer(DASSat ~ IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | PTNUM), couples_baseline_clin))
summary(m <- lmer(DASSat ~ PAIBORtot.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | PTNUM), couples_baseline_clin)) #knock each couples_baseline_clin out
summary(m <- lmer(DASSat ~ narci_sidp.c + bordl_sidp.c + antso_sidp.c + avoid_sidp.c + (1 | PTNUM), couples_baseline_clin))
summary(m <- lmer(DASSat ~ narci_sidp.c + bordl_sidp.c + antso_sidp.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | PTNUM), couples_baseline_clin))

summary(m <- lmer(DASCoh ~ IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | PTNUM), couples_baseline_clin))
summary(m <- lmer(DASCoh ~ PAIBORtot.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | PTNUM), couples_baseline_clin)) #knock each couples_baseline_clin out
summary(m <- lmer(DASCoh ~ PAIBORtot.c + (1 | PTNUM), couples_baseline_clin))
summary(m <- lmer(DASCoh ~ PAIBORtot.c + bordl_sidp.c + OPD_sidp.c +  (1 | PTNUM), couples_baseline_clin))
summary(m <- lmer(DASCoh ~ PAIBORtot.c + OPD_sidp.c +  (1 | PTNUM), couples_baseline_clin))
summary(m <- lmer(DASCoh ~ bordl_sidp.c + OPD_sidp.c + (1 | PTNUM), couples_baseline_clin))
summary(m <- lmer(DASCoh ~ narci_sidp.c + bordl_sidp.c + antso_sidp.c + avoid_sidp.c + depen_sidp.c + (1 | PTNUM), couples_baseline_clin))
summary(m <- lmer(DASCoh ~ narci_sidp.c + bordl_sidp.c + antso_sidp.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | PTNUM), couples_baseline_clin))

#applies to total, too... effect of perp on greater relationship satisfaction goes away when victim removed (likely collinearity problem given .61 correlation between victim and perp)
summary(m <- lmer(DASTotal ~ narci_sidp.c + bordl_sidp.c + antso_sidp.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + Victim.c + Perp.c + (1 | PTNUM), couples_baseline_clin))

summary(m <- lmer(DASTotal ~ narci_sidp.c + bordl_sidp.c + antso_sidp.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + Perp.c + (1 | PTNUM), couples_baseline_clin))
summary(m <- lmer(DASTotal ~ narci_sidp.c + bordl_sidp.c + antso_sidp.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + Victim.c + (1 | PTNUM), couples_baseline_clin)) #this holds up
summary(m <- lmer(DASTotal ~ bordl_sidp.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + Victim.c + (1 | PTNUM), couples_baseline_clin)) #BPD Sx alone not informative




mixed(DASTotal ~ IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + Victim.c + narci_sidp.c + bordl_sidp.c + antso_sidp.c + (1 | PTNUM), couples_baseline_clin, method="PB")

