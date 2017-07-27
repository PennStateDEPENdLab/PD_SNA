#Associations between vertex measures of ego-centered graphs and PDs
#For network-wide metrics of ego networks, see Ego_Graph_Analysis.R
library(plyr)
library(lme4)
library(igraph)
library(afex)
library(ggplot2)
library(lattice)
setwd("/Users/michael/Tresors/PD_SNA/Couples_Baseline_SNA")
source("code/APIM_Functions.R")
source(file.path(getMainDir(), "Miscellaneous", "Global_Functions.R"))
load(file="data/SNA_Processed_6Oct2015.RData")
load(file="data/ego_graph_measures_6Oct2015.RData")
load(file="data/selfreports/couples_baseline_clinical_9Oct2015.RData")

mvertex <- merge(vertex_agg, couples_baseline_clin, by=c("UsrID", "PTNUM"), all=FALSE)

#merge back in alter ratings to vertex data.frame to allow for tests of happy with degree, for example
#note that these are embedded in the graphlist, so we have some work to do
vertexAttribs <- do.call(rbind, lapply(graphlist, function(x) {
      df <- x$actors
      df$UsrID <- x$UsrID
      df$PTNUM <- x$PTNUM
      df$alter_name <- df$name #for name match in merge
      df
    }))

mvertex <- merge(mvertex, vertexAttribs[,c("UsrID", "PTNUM", "alter_name", "CutOff", "Close", "Happy", "Angry", "attach_total")],
    by=c("UsrID", "PTNUM", "alter_name"))

predstocenter <- c("PAIBORtot", "iip_agency", "iip_communion", "iip_elevation", "IIP_PD1", "IIP_PD2", "IIP_PD3", "ECRanx", "ECRavoid", "DASTotal", 
    "Victim", "Perp", "PsychAggV", "PsychAggP", "PhysAssV", "PhysAssP", "bordl_sidp", "narci_sidp", "antso_sidp", "avoid_sidp", "OPD_sidp", "depen_sidp",
    "nonarc", "nobpd", "pdtot", "Happy", "Angry")

#mvertex$UsrID <- factor(mvertex$UsrID)
#mvertex$PTNUM <- factor(mvertex$PTNUM)
mvertex <- f_centerPredictors(mvertex, predstocenter, addsuffix=".c")
mvertex <- ddply(mvertex, .(UsrID), function(subdf) {
      subdf$alter_num <- 1:nrow(subdf)
      subdf
    })

mvertex_melt <- melt(mvertex[,sapply(mvertex, is.numeric)], id.vars=c("PTNUM", "DyadID", "alter_num"))
mvertex_wide <- dcast(mvertex_melt, PTNUM + alter_num ~ variable + DyadID) #drops UsrID because not mentioned
rm(mvertex_melt)

#duped from Couples_SNA.R... should clean up redundancy later
mvertex_wide <- plyr::rename(mvertex_wide, c(iip_communion_0="aff_0", iip_communion_1="aff_1",
        iip_agency_0="dom_0", iip_agency_1="dom_1",
        iip_elevation_0="iipe_0", iip_elevation_1="iipe_1",
        szoid_sidp_0="sz_0", szoid_sidp_1="sz_1",
        stypl_sidp_0="styp_0", stypl_sidp_1="styp_1",
        parnd_sidp_0="par_0", parnd_sidp_1="par_1",
        bordl_sidp_0="bpd_0", bordl_sidp_1="bpd_1",
        histr_sidp_0="his_0", histr_sidp_1="his_1",
        narci_sidp_0="nar_0", narci_sidp_1="nar_1",
        antso_sidp_0="ant_0", antso_sidp_1="ant_1",
        avoid_sidp_0="avd_0", avoid_sidp_1="avd_1",
        depen_sidp_0="depen_0", depen_sidp_1="depen_1",
        obcmp_sidp_0="obcmp_0", obcmp_sidp_1="obcmp_1",
        IIP_PD1_0="pd1_0", IIP_PD1_1="pd1_1",
        IIP_PD2_0="pd2_0", IIP_PD2_1="pd2_1",
        IIP_PD3_0="pd3_0", IIP_PD3_1="pd3_1",
        PAIBORtot_0="pbor_0", PAIBORtot_1="pbor_1",
        strength_2_0="str2_0", strength_2_1="str2_1",
        strength_3_0="str3_0", strength_3_1="str3_1"))

#vertex analyses
# graph thresholding: _2 removes all links for "have never met" and "not close at all"
#                     _3 removes all links for "have never met", "not close at all", and "slightly close"
#
# measures:
#   degree_2: number of connections at 2 threshold
#   degree_3: number of connections at 3 threshold
#   strength_2: summed strengths (closeness) at 2 threshold
#   strength_3: summed strengths (closeness) at 3 threshold
#   evcent_weighted_2: weighted eigenvector centrality without 1s and 2s
#   evcent_binary_2: binary evcent at 2 thresh
#   evcent_weighted_3: weighted evcent at 3 thresh
#   evcent_binary_3: binary evcent at 3 thresh
#   closeness_weighted_2: closeness centrality (1/sum of distances to all other nodes)
#   closeness_binary_2:
#   closeness_weighted_3:
#   closeness_binary_3:
#   betweenness_weighted_2:
#   betweenness_binary_2:
#   betweenness_weighted_3:
#   betweenness_binary_3:
#   eccentricity_2: maximum of the shortest distance to all other vertices 
#   eccentricity_3
#   locclust_binary_2
#   locclust_weighted_2
#   locclust_binary_3
#   locclust_weighted_3

# TODO: look at betweenness centrality for patient and partner specifically (how much traffic flows through them)

# DEGREE CENTRALITY

#summary(lmer(degree_2 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))
#summary(lmer(degree_2 ~ PAIBORtot.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | UsrID) + (1 | PTNUM), mvertex))
#summary(lmer(degree_2 ~ IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | UsrID) + (1 | PTNUM), mvertex))
mvertex$degree_2_sqrt <- sqrt(mvertex$degree_2) #positive skew
summary(lmer(degree_2 ~ iip_elevation.c + iip_communion.c + iip_agency.c + DyadID + (1 | UsrID) + (1 | PTNUM), mvertex)) #higher communion associated with greater degree
(mixed(degree_2 ~ iip_elevation.c + iip_communion.c + iip_agency.c + DyadID + (1 | UsrID) + (1 | PTNUM), mvertex))
(mixed(degree_2_sqrt ~ iip_elevation.c + iip_communion.c + iip_agency.c + DyadID + (1 | UsrID) + (1 | PTNUM), mvertex)) #holds with transform
(mixed(degree_2 ~ pdtot + DyadID + (1 | UsrID) + (1 | PTNUM), mvertex)) #non-bpd sig .01
summary(lmer(degree_2 ~ bordl_sidp.c + DyadID + (1 | UsrID) + (1 | PTNUM), mvertex)) #null effect
summary(lmer(degree_2 ~ narci_sidp.c + DyadID + (1 | UsrID) + (1 | PTNUM), mvertex)) #null
summary(lmer(degree_2 ~ narci_sidp.c + nonarc.c + DyadID + (1 | UsrID) + (1 | PTNUM), mvertex)) #null
summary(glmer(degree_2 ~ bordl_sidp.c + nobpd.c + DyadID + (1 | UsrID) + (1 | PTNUM), mvertex, family=poisson)) #non-bpd sig
(mixed(degree_2 ~ bordl_sidp.c + nobpd.c + DyadID + (1 | UsrID) + (1 | PTNUM), mvertex)) #non-bpd sig .01
(mixed(degree_2_sqrt ~ bordl_sidp.c + nobpd.c + DyadID + (1 | UsrID) + (1 | PTNUM), mvertex)) #remains after sqrt transform

#should look at Clifton analysis: are more central people in the network also rated more positively?
summary(glmer(degree_2 ~ bordl_sidp.c + nobpd.c + DyadID + (1 | UsrID) + (1 | PTNUM), mvertex, family=poisson)) #non-bpd sig
summary(glmer(degree_2 ~ Happy.c + DyadID + (1 | UsrID) + (1 | PTNUM), mvertex, family=poisson)) #strong coupling
summary(glmer(degree_2 ~ Angry.c*bordl_sidp.c + DyadID + (1 | UsrID) + (1 | PTNUM), mvertex, family=poisson)) #why is this positive?

#this looks like a replication -- take into APIM
mvertex$DyadID_fac <- factor(mvertex$DyadID, levels=c(0,1), labels=c("partner", "patient"))
summary(m <- glmer(degree_2 ~ Happy.c*bordl_sidp.c + DyadID_fac + (1 | UsrID) + (1 | PTNUM), mvertex, family=poisson)) 

#holds when controlling for OPD
summary(m <- glmer(degree_2 ~ Happy.c*bordl_sidp.c + Happy.c*nobpd.c + DyadID_fac + (1 | UsrID) + (1 | PTNUM), mvertex, family=poisson)) 

summary(m <- glmer(degree_2 ~ Happy.c*bordl_sidp.c + Happy.c*nobpd.c + DyadID_fac + (1 + Happy.c | UsrID) + (1 | PTNUM), mvertex, family=poisson))

summary(m <- glmer(degree_2 ~ Happy.c*bordl_sidp.c + Happy.c*nobpd.c + DyadID_fac + (1 + Happy.c | UsrID) + (1 | PTNUM), mvertex, family=poisson))

#without centering (for ease of plot axes)
summary(m <- glmer(degree_2 ~ Happy*bordl_sidp + DyadID_fac + (1 | UsrID) + (1 | PTNUM), mvertex, family=poisson)) #why is this positive?
cm <- lmerCellMeans(m, divide="bordl_sidp")
#average patient and partner estimates (since this is just offset)
cm <- aggregate(cm, list(cm$bordl_sidp, cm$Happy), function(x) { if (is.numeric(x)) mean(x) else head(x, n=1) })
ggplot(cm, aes(x=Happy, y=degree_2, color=bordl_sidp)) + geom_line()

#look at anger
#huh, anger is associated with greater centrality
#this may represent the tendency to be mad at a couple of people who know each other
summary(m <- glmer(degree_2 ~ Angry.c*bordl_sidp.c + DyadID_fac + (1 | UsrID) + (1 | PTNUM), mvertex, family=poisson))
summary(m <- glmer(degree_2 ~ Angry.c*bordl_sidp.c + Angry.c*nobpd.c + DyadID_fac + (1 | UsrID) + (1 | PTNUM), mvertex, family=poisson)) 

summary(m <- glmer(degree_2 ~ Angry.c*narci_sidp.c + DyadID_fac + (1 | UsrID) + (1 | PTNUM), mvertex, family=poisson))
summary(m <- glmer(degree_2 ~ Angry.c*narci_sidp.c + Angry.c*nonarc.c + DyadID_fac + (1 | UsrID) + (1 | PTNUM), mvertex, family=poisson))

summary(m <- glmer(degree_2 ~ Angry*narci_sidp + DyadID_fac + (1 | UsrID) + (1 | PTNUM), mvertex, family=poisson)) #why is this positive?
cm <- lmerCellMeans(m, divide="narci_sidp")
#average patient and partner estimates (since this is just offset)
cm <- aggregate(cm, list(cm$narci_sidp, cm$Angry), function(x) { if (is.numeric(x)) mean(x) else head(x, n=1) })
ggplot(cm, aes(x=Angry, y=degree_2, color=narci_sidp)) + geom_line()

summary(m <- glmer(degree_3 ~ Angry*narci_sidp + DyadID_fac + (1 | UsrID) + (1 | PTNUM), mvertex, family=poisson)) #why is this positive?
cm <- lmerCellMeans(m, divide="narci_sidp")
#average patient and partner estimates (since this is just offset)
cm <- aggregate(cm, list(cm$narci_sidp, cm$Angry), function(x) { if (is.numeric(x)) mean(x) else head(x, n=1) })
ggplot(cm, aes(x=Angry, y=degree_3, color=narci_sidp)) + geom_line()

summary(m <- glmer(strength_2 ~ Angry.c*narci_sidp.c + Angry.c*nonarc.c + DyadID_fac + (1 | UsrID) + (1 | PTNUM), mvertex, family=poisson))

#both Narc and BPD Sx associated with weaker anger - degree coupling 
#random slopes for anger -> strength coupling yields MUCH better AIC
#summary(m <- glmer(strength_2 ~ Angry.c*narci_sidp.c + Angry.c*bordl_sidp.c + DyadID_fac + (1 + Angry.c | UsrID) + (1 | PTNUM), mvertex, family=poisson))
#summary(m <- glmer(strength_2 ~ Angry.c*narci_sidp.c + Angry.c*bordl_sidp.c + DyadID_fac + (1 | UsrID) + (1 | PTNUM), mvertex, family=poisson))
#this makes it fall apart
summary(m <- glmer(strength_2 ~ Angry.c*narci_sidp.c + Angry.c*nonarc.c + DyadID_fac + (1 +Angry.c 	| UsrID) + (1 | PTNUM), mvertex, family=poisson))
summary(m <- glmer(degree_3 ~ Angry.c*narci_sidp.c + Angry.c*nonarc.c + DyadID_fac + (1 +Angry.c 	| UsrID) + (1 | PTNUM), mvertex, family=poisson))


#cutoffs
summary(m <- glmer(degree_2 ~ CutOff*bordl_sidp.c + DyadID_fac + (1 | UsrID) + (1 | PTNUM), mvertex, family=poisson))
summary(m <- glmer(degree_2 ~ CutOff*bordl_sidp.c + CutOff*nobpd.c + DyadID_fac + (1 | UsrID) + (1 | PTNUM), mvertex, family=poisson))
summary(m <- glmer(degree_3 ~ CutOff*bordl_sidp.c + DyadID_fac + (1 | UsrID) + (1 | PTNUM), mvertex, family=poisson))
summary(m <- glmer(degree_3 ~ CutOff*bordl_sidp.c + CutOff*nobpd.c + DyadID_fac + (1 | UsrID) + (1 | PTNUM), mvertex, family=poisson))

mvertex_wide$bpd0c <- mvertex_wide$bpd_0 - mean(mvertex_wide$bpd_0, na.rm=TRUE)
mvertex_wide$bpd1c <- mvertex_wide$bpd_1 - mean(mvertex_wide$bpd_1, na.rm=TRUE)
mvertex_wide$hap0c <- mvertex_wide$Happy_0 - mean(mvertex_wide$Happy_0, na.rm=TRUE)
mvertex_wide$hap1c <- mvertex_wide$Happy_1 - mean(mvertex_wide$Happy_1, na.rm=TRUE)
mvertex_wide$ang0c <- mvertex_wide$Angry_0 - mean(mvertex_wide$Angry_0, na.rm=TRUE)
mvertex_wide$ang1c <- mvertex_wide$Angry_1 - mean(mvertex_wide$Angry_1, na.rm=TRUE)

#here's the output for Mplus MLAPIM
#the challenge is that we're talking about random slopes at within happy -> degree for actor and partner
#so it's a very slow model
prepareMplusData(mvertex_wide, filename="hap_bpd_mlapim.dat", keepCols=c("PTNUM", "alter_num",
        "degree_2_0", "degree_2_1", "bpd0c", "bpd1c", "hap0c", "hap1c", "ang0c", "ang1c"),
    inpfile="hap_bpd_mlapim.inp")
#output <- mplusMLAPIM(mvertex_wide, "degree_2", "dom") #null



#MLAPIM
#degree probably best modeled by Poisson since it is a count distribution (not overly inflated)
output <- mplusMLAPIM(mvertex_wide, "degree_2", "dom") #null
output <- mplusMLAPIM(mvertex_wide, "degree_2", "aff") #wow, all negative partner effect!
output <- mplusMLAPIM(mvertex_wide, "degree_2", "aff", count=TRUE) #wow, all partner effect!
output <- mplusMLAPIM(mvertex_wide, "degree_2", c("aff", "dom", "iipe")) #altogether: aff holds up
output <- mplusMLAPIM(mvertex_wide, "degree_2", c("aff", "dom", "iipe"), count=TRUE) #altogether: aff holds up

output <- mplusMLAPIM(mvertex_wide, "degree_2", c("aff", "dom"), count=TRUE) #altogether: aff holds up
output$free$parameters
output <- mplusMLAPIM(mvertex_wide, "degree_2", c("bpd"), count=TRUE)
output <- mplusMLAPIM(mvertex_wide, "degree_2", c("bpd", "nobpd"), count=TRUE) #same direction as MLM, but marginal
output <- mplusMLAPIM(mvertex_wide, "degree_2", c("nobpd"), count=TRUE) #ugh, flip
output <- mplusMLAPIM(mvertex_wide, "degree_2", c("nobpd")) #ugh, flip
output <- mplusMLAPIM(mvertex_wide, "degree_2", c("pdtot"), count=TRUE) #weird 
output <- mplusMLAPIM(mvertex_wide, "degree_2", c("nar", "nonarc"), count=TRUE) #null

#check effects at threshold=3
output <- mplusMLAPIM(mvertex_wide, "degree_3", "aff", count=TRUE) #still significant negative partner
output <- mplusMLAPIM(mvertex_wide, "degree_3", "dom", count=TRUE) #weak asymmetric partner effect
output <- mplusMLAPIM(mvertex_wide, "degree_3", "iipe", count=TRUE) #weak asymmetric partner effect
output <- mplusMLAPIM(mvertex_wide, "degree_3", c("aff", "dom", "iipe"), count=TRUE) #holds up

#summary(lmer(degree_3 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))
#summary(lmer(degree_3 ~ PAIBORtot.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | UsrID) + (1 | PTNUM), mvertex))
#summary(lmer(degree_3 ~ iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))
#summary(lmer(degree_3 ~ PAIBORtot.c + ECRanx.c + ECRavoid.c + (1 | UsrID) + (1 | PTNUM), mvertex))


#strength (weighted degree)
#summary(lmer(strength_2 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + DyadID + (1 | UsrID) + (1 | PTNUM), mvertex))
#summary(lmer(strength_3 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + DyadID + (1 | UsrID) + (1 | PTNUM), mvertex))

histogram(~strength_2 | DyadID, mvertex)
histogram(~sqrt(strength_2) | DyadID, mvertex) #much more normal...

#probably also use a Poisson here, but the mean is much higher
output <- mplusMLAPIM(mvertex_wide, "str2", "aff", count=TRUE) #consistent with degree
output <- mplusMLAPIM(mvertex_wide, "str2", "dom", count=TRUE) #null
output <- mplusMLAPIM(mvertex_wide, "str2", "iipe", count=TRUE) #null
output <- mplusMLAPIM(mvertex_wide, "str2", c("iipe", "aff", "dom"), count=FALSE) #washed out a bit
output <- mplusMLAPIM(mvertex_wide, "str2", c("aff"), count=FALSE)

output <- mplusMLAPIM(mvertex_wide, "str2", c("bpd", "nobpd"), count=TRUE)
output <- mplusMLAPIM(mvertex_wide, "str2", c("bpd"), count=TRUE)
output <- mplusMLAPIM(mvertex_wide, "str2", c("pdtot"), count=FALSE)
output <- mplusMLAPIM(mvertex_wide, "str2", c("nar", "nonarc"), count=TRUE)



#eigenvector centrality
summary(lmer(evcent_weighted_2 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))
summary(lmer(evcent_weighted_2 ~ PAIBORtot.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | UsrID) + (1 | PTNUM), mvertex))

summary(lmer(evcent_weighted_3 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))
summary(lmer(evcent_weighted_3 ~ PAIBORtot.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | UsrID) + (1 | PTNUM), mvertex))

summary(lmer(evcent_binary_2 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))
summary(lmer(evcent_binary_2 ~ PAIBORtot.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | UsrID) + (1 | PTNUM), mvertex))

summary(lmer(evcent_binary_3 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))
summary(lmer(evcent_binary_3 ~ PAIBORtot.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | UsrID) + (1 | PTNUM), mvertex))

#closeness centrality: both communion and agency increase closeness
summary(lmer(closeness_weighted_2 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))
summary(lmer(closeness_weighted_3 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))

mixed(closeness_weighted_3 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex, method="KR")

summary(lmer(closeness_binary_2 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))
summary(lmer(closeness_binary_3 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))

#betweenness centrality
summary(lmer(betweenness_binary_2 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))

#agency goes with > betweenness; BPD with weaker betweenness at thresh=3 
summary(lmer(betweenness_binary_3 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))

summary(lmer(betweenness_weighted_2 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))
summary(lmer(betweenness_weighted_3 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))

summary(lmer(eccentricity_2 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))

#longest distance between two nodes is weaker in BPD and greater in agency? Would have predicted > in BPD
summary(lmer(eccentricity_3 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))
summary(lmer(eccentricity_3 ~ PAIBORtot.c + (1 | UsrID) + (1 | PTNUM), mvertex))

#weaker local clustering on average in BPD (thresh = 2) -- matches transitivity above
summary(lmer(locclust_binary_2 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))
summary(lmer(locclust_binary_3 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))
summary(lmer(locclust_weighted_2 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))
summary(lmer(locclust_weighted_3 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))

