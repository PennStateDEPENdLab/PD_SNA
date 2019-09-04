library(lme4)
library(lavaan)
library(igraph)
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
corstarsl(couples_baseline_clin[,c("IIP_agency", "IIP_communion", "IIP_pd1", "IIP_pd2", "IIP_pd3")])
predstocenter <- c("PAIBORtot", "IIP_agency", "IIP_communion", "IIP_pd1", "IIP_pd2", "IIP_pd3", "ECRanx", "ECRavoid", "DASTotal", 
    "CTS_Victim", "CTS_Perp", "CTS_PsychAggV", "CTS_PsychAggP", "CTS_PhysAssV", "CTS_PhysAssP", "bordl_sidp", "narci_sidp", 
    "antso_sidp", "avoid_sidp", "OPD_sidp", "depen_sidp", "parnd_sidp", "histr_sidp", "nobpd")

couples_baseline_clin <- f_centerPredictors(couples_baseline_clin, predstocenter, addsuffix=".c")
couples_clin_wide <- f_centerPredictors(couples_clin_wide, as.vector(outer(predstocenter, c("0", "1"), paste, sep="_")), addsuffix=".c")

#the merge in Import_Couples_SNA.R doesn't include SIDP etc.
sna_merge <- merge(sna, couples_baseline_clin, by=c("UsrID", "PTNUM"))
sna_merge <- f_centerPredictors(sna_merge, predstocenter, addsuffix=".c")
sna_merge <- ddply(sna_merge, .(UsrID), function(subdf) {
      subdf$alter_num <- 1:nrow(subdf)
      subdf
    })

#dplyr version
library(dplyr)
#sna_merge2 <- sna_merge %>% group_by(UsrID) %>% mutate(alter_num2=row_number())
sna_merge2 <- sna_merge %>% group_by(UsrID) %>% mutate(alter_num2=1:length(UsrID)) #OR use 1:length of a variable


sna_merge$YrsKnown[which(sna_merge$YrsKnown == -1)] <- 0 #less than one year shows up as -1 
sna_merge$PTNUM <- as.numeric(as.character(sna_merge$PTNUM))
sna_merge$UsrID <- as.numeric(as.character(sna_merge$UsrID))
sna_merge$pastRom <- as.numeric(sna_merge$Romance=="past") #only about 5%... not a lot to work with
sna_merge$pdtot <- sna_merge$bordl_sidp + sna_merge$nobpd
sna_merge$nonarc <- with(sna_merge, szoid_sidp + stypl_sidp + parnd_sidp + histr_sidp + bordl_sidp +
        antso_sidp + obcmp_sidp + depen_sidp + avoid_sidp) #non-narc Sx
sna_merge$attach_figure <- as.numeric(sna_merge$attach_figure) #was stored as logical upstream

vertex_melt <- melt(sna_merge[,sapply(sna_merge, is.numeric)], id.vars=c("PTNUM", "DyadID", "alter_num"))
vertex_wide <- dcast(vertex_melt, PTNUM + alter_num ~ variable + DyadID)

##Conduct ML SEM analyses of vertex attributes
##Needs to be multilevel to allow for an APIM approach since alters are nested within individual

##At the moment, Mplus seems like the only way to get this done
#shorten some key variable names
vertex_wide <- plyr::rename(vertex_wide, c(IIP_communion_0="aff_0", IIP_communion_1="aff_1",
        IIP_agency_0="dom_0", IIP_agency_1="dom_1",
        IIP_elevation_0="IIPe_0", IIP_elevation_1="IIPe_1",
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
        IIP_pd1_0="pd1_0", IIP_pd1_1="pd1_1",
        IIP_pd2_0="pd2_0", IIP_pd2_1="pd2_1",
        IIP_pd3_0="pd3_0", IIP_pd3_1="pd3_1",
        PAIBORtot_0="pbor_0", PAIBORtot_1="pbor_1",
        attach_figure_0="afig_0", attach_figure_1="afig_1",
        attach_total_0="arat_0", attach_total_1="arat_1",
        FaceContact_0="face_0", FaceContact_1="face_1"
        ))


##Cutoff analyses
output <- mplusMLAPIM(vertex_wide, "CutOff", "dom", categorical=TRUE)
output <- mplusMLAPIM(vertex_wide, "CutOff", "aff", categorical=TRUE)
output <- mplusMLAPIM(vertex_wide, "CutOff", c("aff", "dom", "IIPe"), categorical=TRUE) #dom and IIPe associated with more cutoffs
output <- mplusMLAPIM(vertex_wide, "CutOff", "pd1", categorical=TRUE) #interpersonal sensitivity; all actor effects
output <- mplusMLAPIM(vertex_wide, "CutOff", "pd2", categorical=TRUE) #ambivalence; all actor effects
output <- mplusMLAPIM(vertex_wide, "CutOff", "pd3", categorical=TRUE) #aggression; all actor effects
output <- mplusMLAPIM(vertex_wide, "CutOff", c("pd1", "pd2", "pd3"), categorical=TRUE) #combined

output <- mplusMLAPIM(vertex_wide, "CutOff", "par", categorical=TRUE) #both actor effects, indistinguishable
output <- mplusMLAPIM(vertex_wide, "CutOff", "sz", categorical=TRUE) #both actor effects, mostly indistinguishable (slight partner effect in free)
output <- mplusMLAPIM(vertex_wide, "CutOff", "styp", categorical=TRUE) #both actor and partner effects.
output <- mplusMLAPIM(vertex_wide, "CutOff", "bpd", categorical=TRUE) #both actor effects, indistinguishable
output <- mplusMLAPIM(vertex_wide, "CutOff", "ant", categorical=TRUE) #Mostly actor effects, free model slightly better
output <- mplusMLAPIM(vertex_wide, "CutOff", "nar", categorical=TRUE) #both actor effects, indistinguishable
output <- mplusMLAPIM(vertex_wide, "CutOff", "his", categorical=TRUE) #both actor effects, indistinguishable
output <- mplusMLAPIM(vertex_wide, "CutOff", "avd", categorical=TRUE) #NULL effects
output <- mplusMLAPIM(vertex_wide, "CutOff", "depen", categorical=TRUE) #NULL effects
output <- mplusMLAPIM(vertex_wide, "CutOff", "obcmp", categorical=TRUE) #NULL effects
output <- mplusMLAPIM(vertex_wide, "CutOff", "OPD_sidp", categorical=TRUE) #both actor effects, indistinguishable
output <- mplusMLAPIM(vertex_wide, "CutOff", c("pd1", "pd3"), categorical=TRUE) #aggression; mostly actor effects
output <- mplusMLAPIM(vertex_wide, "CutOff", c("bpd", "nar", "ant"), categorical=TRUE) #aggression; mostly actor effects
output <- mplusMLAPIM(vertex_wide, "CutOff", "OPD_sidp", categorical=TRUE) #both actor and partner effects here.
output <- mplusMLAPIM(vertex_wide, "CutOff", "bpd", categorical=TRUE)
output <- mplusMLAPIM(vertex_wide, "CutOff", "pdtot", categorical=TRUE) #all PDs
output <- mplusMLAPIM(vertex_wide, "CutOff", c("bpd", "nobpd"), categorical=TRUE) #this may be the most relevant approach given the sampling scheme
output <- mplusMLAPIM(vertex_wide, "CutOff", c("nar", "nonarc"), categorical=TRUE) #this may be the most relevant approach given the sampling scheme

#Anger analyses (treat as Poisson)
#TODO: double check that we exclude the patient and partner from ratings
output <- mplusMLAPIM(vertex_wide, "Angry", "OPD_sidp", count=TRUE) #NULL
output <- mplusMLAPIM(vertex_wide, "Angry", "bpd", count=TRUE) #NULL  
output <- mplusMLAPIM(vertex_wide, "Angry", "pdtot", count=TRUE) #NULL
output <- mplusMLAPIM(vertex_wide, "Angry", c("bpd", "OPD_sidp"), count=TRUE) #this may be the most relevant approach given the sampling scheme
output <- mplusMLAPIM(vertex_wide, "Angry", "pbor", count=TRUE) #NULL (self-report)
output <- mplusMLAPIM(vertex_wide, "Angry", "nar", count=TRUE) #only marginal here -- see ZIP model output created manually
output <- mplusMLAPIM(vertex_wide, "Angry", c("nar", "nonarc"), count=TRUE) #include nonarc

output_poi <- mplusMLAPIM(vertex_wide, "Angry", c("dom", "aff", "IIPe"), count=TRUE)
output_zip <- mplusMLAPIM(vertex_wide, "Angry", c("dom", "aff", "IIPe"), zip=TRUE)
output_poi$indistinguishable$parameters$unstandardized
output_zip$indistinguishable$parameters$unstandardized
output_poi$indistinguishable$summaries$AICC
output_zip$indistinguishable$summaries$AICC

output <- mplusMLAPIM(vertex_wide, "Angry", c("dom"), count=TRUE)
output <- mplusMLAPIM(vertex_wide, "Angry", c("dom"), zip=TRUE)
output <- mplusMLAPIM(vertex_wide, "Angry", c("aff"), count=TRUE) #no dice
output <- mplusMLAPIM(vertex_wide, "Angry", c("IIPe"), count=TRUE) #actor and partner effects, some evidence of distinguishability of dyads

#MLM approach
#anger as a function of BPD Sx
#NB: this approach is flawed because PAI/BPD is confounded with DyadID -- there is a mean shift in Angry between partners and patients!
#summary(m <- lmer(Angry ~ PAIBORtot.c + (1 | UsrID), data=sna_merge, REML=TRUE))

#angry looks more like a ZIP or Poisson-distributed variable
#need to account for patient x partner effects here to accurately represent anger levels as a function of BPD
summary(mnest <- glmer(Angry ~ PAIBORtot.c + DyadID + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=poisson)) # NS .11

#treated as Gaussian, the PAI effect is p = .04
#(mnest <- mixed(Angry ~ PAIBORtot.c + DyadID + (1 | UsrID) + (1 | PTNUM), data=sna_merge, REML=TRUE))

#does not hold up for interview-based Sx
#(mnest <- mixed(Angry ~ bordl_sidp.c + DyadID + (1 | UsrID) + (1 | PTNUM), data=sna_merge, REML=TRUE))
#(mnest <- mixed(Angry ~ narci_sidp.c + DyadID + (1 | UsrID) + (1 | PTNUM), data=sna_merge, REML=TRUE))

#okay, look at angry as a function of other PD Sx like Narc
summary(mnest <- glmer(Angry ~ bordl_sidp.c + DyadID + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=poisson))
summary(mnest <- glmer(Angry ~ narci_sidp.c + DyadID + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=poisson)) #sig at .04
summary(mnest2 <- glmer(Angry ~ narci_sidp.c * DyadID + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=poisson))
anova(mnest, mnest2) #not a better fit

summary(mnest <- glmer(Angry ~ OPD_sidp.c + DyadID + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=poisson)) #NS
summary(mnest <- glmer(Angry ~ antso_sidp.c + DyadID + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=poisson)) #NS
summary(mnest <- glmer(Angry ~ histr_sidp.c + DyadID + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=poisson)) #NS
summary(mnest <- glmer(Angry ~ avoid_sidp.c + DyadID + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=poisson)) #NS
summary(mnest <- glmer(Angry ~ depen_sidp.c + DyadID + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=poisson)) #NS
summary(mnest <- glmer(Angry ~ parnd_sidp.c + DyadID + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=poisson)) #NS

#conclusion: both Poisson MLM and ZIP APIM indicate that Narci is the main driver of anger ratings

#BPD and OPD correlate at .73...
#summary(mnest <- glmer(Angry ~ bordl_sidp.c + DyadID + OPD_sidp.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=poisson))
#summary(mnest2 <- glmer(Angry ~ bordl_sidp.c*DyadID + OPD_sidp.c*DyadID + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=poisson))
#anova(mnest, mnest2) #interaction model not better

#Happy analyses
output <- mplusMLAPIM(vertex_wide, "Happy", "OPD_sidp") #lower actor happiness in OPD
output <- mplusMLAPIM(vertex_wide, "Happy", "bpd") #lower actor happiness in BPD
output <- mplusMLAPIM(vertex_wide, "Happy", "pdtot") #lower actor happiness in all PDs
output <- mplusMLAPIM(vertex_wide, "Happy", c("bpd", "nobpd")) #When they duke it out, OPD wins
output <- mplusMLAPIM(vertex_wide, "Happy", c("pbor", "OPD_sidp")) #OPD beats self-reported BPD, too
output <- mplusMLAPIM(vertex_wide, "Happy", c("nar"))
output <- mplusMLAPIM(vertex_wide, "Happy", c("nar", "nonarc")) #Partner effect remains after controlling for lower happiness nonarc actor effect

output <- mplusMLAPIM(vertex_wide, "Happy", c("IIPe", "aff", "dom"))

output <- mplusMLAPIM(vertex_wide, "Happy", c("IIPe", "aff", "dom"))


summary(mnest <- lmer(Happy ~ narci_sidp.c * DyadID + (1 | UsrID) + (1 | PTNUM), data=sna_merge, REML=TRUE))
##older stuff


#Romance analyses
output <- mplusMLAPIM(vertex_wide, "pastRom", c("pdtot"), categorical=TRUE) #more past partners in PDs
output <- mplusMLAPIM(vertex_wide, "pastRom", c("OPD_sidp"), categorical=TRUE) #more past partners in PDs
output <- mplusMLAPIM(vertex_wide, "pastRom", c("bpd"), categorical=TRUE) #holds for BPD proper
output <- mplusMLAPIM(vertex_wide, "pastRom", c("bpd", "OPD_sidp"), categorical=TRUE) #basically washes out when controlling for OPD
output <- mplusMLAPIM(vertex_wide, "pastRom", c("nar", "nonarc"), categorical=TRUE) #kind of a mess here... no consistency

output <- mplusMLAPIM(vertex_wide, "pastRom", c("IIPe", "aff", "dom"), categorical=TRUE) #kind of a mess here... no consistency

summary(mnest <- glmer(pastRom ~ narci_sidp.c * DyadID + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=binomial))
summary(mnest <- glmer(pastRom ~ bordl_sidp.c * DyadID + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=binomial)) #holds up in MLM
summary(mnest <- glmer(pastRom ~ bordl_sidp.c + OPD_sidp.c + DyadID + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=binomial)) #washes out here, too, goes to marginal

#Relationship duration
histogram(vertex_wide$YrsKnown_0)
histogram(vertex_wide$YrsKnown_1)

#positive skew that is amenable to sqrt transformation... or could just control for participant age
t.test(vertex_wide$YrsKnown_0, vertex_wide$YrsKnown_1) #patients have shorter relationships.... sig effect
sd(vertex_wide$YrsKnown_0)
sd(vertex_wide$YrsKnown_1, na.rm=TRUE) #patients have shorter relationships.... big effect

#BPD Sx associated with shorter durations, but not specific... so PDs associated with shorter relationships
summary(mnest <- lmer(YrsKnown ~ bordl_sidp.c + DyadID + p_age + (1 | UsrID) + (1 | PTNUM), data=sna_merge)) #sig
summary(mnest <- lmer(YrsKnown ~ pdtot + DyadID + p_age + (1 | UsrID) + (1 | PTNUM), data=sna_merge)) #sig
(mnest <- mixed(YrsKnown ~ pdtot + DyadID + p_age + (1 | UsrID) + (1 | PTNUM), data=sna_merge)) #sig
summary(mnest <- lmer(YrsKnown ~ narci_sidp.c + DyadID + p_age + (1 | UsrID) + (1 | PTNUM), data=sna_merge)) #NS
summary(mnest <- lmer(YrsKnown ~ narci_sidp.c + nonarc + DyadID + p_age + (1 | UsrID) + (1 | PTNUM), data=sna_merge)) #NS
summary(mnest <- lmer(YrsKnown ~ bordl_sidp.c + OPD_sidp.c + DyadID + p_age + (1 | UsrID) + (1 | PTNUM), data=sna_merge)) #knocked out by OPD

#attachment figures
#look at logistic MLSEM for attachment figures
output <- mplusMLAPIM(vertex_wide, "afig", c("nar", "nonarc"), categorical=TRUE) 
output <- mplusMLAPIM(vertex_wide, "afig", c("bpd", "nobpd"), categorical=TRUE)
output <- mplusMLAPIM(vertex_wide, "afig", c("bpd"), categorical=TRUE) 
output <- mplusMLAPIM(vertex_wide, "afig", c("avd"), categorical=TRUE) 
output <- mplusMLAPIM(vertex_wide, "afig", c("depen"), categorical=TRUE) 

#attachment ratings
output <- mplusMLAPIM(vertex_wide, "arat", c("nar", "nonarc"))
output <- mplusMLAPIM(vertex_wide, "arat", c("bpd", "OPD_sidp"))
output <- mplusMLAPIM(vertex_wide, "arat", c("dom", "aff"))
output <- mplusMLAPIM(vertex_wide, "arat", c("pd1"))
output <- mplusMLAPIM(vertex_wide, "arat", c("pd2"))
output <- mplusMLAPIM(vertex_wide, "arat", c("pd3")) #interesting: partner effect where more problems with anger linked with higher partner attach ratings (buffer?)

#face-to-face contact
output <- mplusMLAPIM(vertex_wide, "face", c("pdtot")) #not indistinguable! only patient actor effect
output <- mplusMLAPIM(vertex_wide, "face", c("bpd", "OPD_sidp"))
output <- mplusMLAPIM(vertex_wide, "face", c("pbor", "OPD_sidp"))
output <- mplusMLAPIM(vertex_wide, "face", c("nar", "nonarc")) #same here -- nonarc sx associated with less contact in patients only
output <- mplusMLAPIM(vertex_wide, "face", c("avd"))
output <- mplusMLAPIM(vertex_wide, "face", c("depen"))
output <- mplusMLAPIM(vertex_wide, "face", c("obcmp"))

#closeness
output <- mplusMLAPIM(vertex_wide, "Close", c("pdtot"))
output <- mplusMLAPIM(vertex_wide, "Close", c("bpd"))
output <- mplusMLAPIM(vertex_wide, "Close", c("bpd", "OPD_sidp")) #opd Sx associated with lower closeness ratings
output <- mplusMLAPIM(vertex_wide, "Close", c("nar", "nonarc")) #
output <- mplusMLAPIM(vertex_wide, "Close", c("aff", "dom")) #



#nested random effect: UsrID nested within PTNUM (Dyad)
xtabs(~ UsrID + PTNUM, sna_merge, drop = TRUE, sparse = TRUE)

#see here for specification of nested factors: need to simply have random intercepts as a function of both UsrID and PTNUM (similar to crossed)
#http://lme4.r-forge.r-project.org/book/Ch2.pdf












#model attachment figures
mnest <- glmer(attach_figure ~ ECRanx.c*PAIBORtot.c + ECRavoid.c*PAIBORtot.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=binomial, nAGQ=1)
summary(mnest)

#dimensional ratings of attachment
mnest <- lmer(attach_total ~ ECRanx.c*PAIBORtot.c + ECRavoid.c*PAIBORtot.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, REML=TRUE)
summary(mnest)

mnest <- lmer(attach_total ~ PAIBORtot.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, REML=TRUE)
summary(mnest)


mnest <- lmer(attach_total ~ PAIBORtot.c*Cutoff*ECRanx.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, REML=TRUE)
summary(mnest)

cm <- lmerCellMeans(mnest, divide=c("ECRanx.c"))

ggplot(cm, aes(x=PAIBORtot.c, y=attach_total, color=ECRanx.c)) + geom_line() + facet_wrap(~Cutoff)





#interesting BPD x Attachment Anxiety effect (washes out when Avoid included): More BPD Sx and greater Attachment anxiety goes with higher attachment ratings, p = .02 
library(afex)
mnest <- lmer(attach_total ~ ECRanx.c*PAIBORtot.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, REML=TRUE)
mnest <- mixed(attach_total ~ ECRanx.c*PAIBORtot.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, REML=TRUE)
summary(mnest)

mnest <- lmer(attach_total ~ ECRavoid.c*PAIBORtot.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, REML=TRUE)
mnest <- mixed(attach_total ~ ECRavoid.c*PAIBORtot.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, REML=TRUE)
summary(mnest)

sna_merge$Cutoff_num <- sapply(sna_merge$Cutoff, function(x) { if (is.na(x)) x else if (x=="yes") 1 else 0 })

mcutoff <- glmer(CutOff ~ ECRanx.c*PAIBORtot.c + ECRavoid.c*PAIBORtot.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=binomial, nAGQ=1)
summary(mcutoff)

mcutoff <- glmer(CutOff ~ IIP_agency.c*PAIBORtot.c + IIP_agency.c*PAIBORtot.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=binomial, nAGQ=1)
summary(mcutoff)

mcutoff <- glmer(CutOff ~ IIP_agency.c + IIP_communion.c + ECRanx.c + ECRavoid.c + PAIBORtot.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=binomial, nAGQ=1)
summary(mcutoff)

mcutoff <- glmer(CutOff ~ IIP_agency.c + IIP_communion.c + bordl_sidp.c + narci_sidp.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=binomial, nAGQ=1)
summary(mcutoff)

mcutoff <- glmer(CutOff ~ IIP_agency.c + IIP_communion.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=binomial, nAGQ=1)
summary(mcutoff)

prepareMplusData(sna_merge, filename="cutoff_test.dat", keepCols=c("UsrID", "PTNUM", "CutOff", "IIP_agency.c", "IIP_communion.c",
        "bordl_sidp.c", "narci_sidp.c"), inpfile="cutoff_test.inp")

prepareMplusData(vertex_wide, filename="cutoff_test_wide.dat", keepCols=c("PTNUM", "CutOff_0", "CutOff_1", "IIP_agency.c_0", 
        "IIP_agency.c_1"), inpfile="cutoff_test_wide.inp")
    
    
    
    mcutoff <- glmer(CutOff ~ IIP_communion.c*PAIBORtot.c + IIP_communion.c*PAIBORtot.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=binomial, nAGQ=1)
summary(mcutoff)



#not much action for attachment anxiety in cutoff
mcutoff_noanx <- glmer(Cutoff_num ~ ECRavoid.c*PAIBORtot.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=binomial, nAGQ=1)
summary(mcutoff_noanx)

mcutoff_noanx <- glmer(Cutoff_num ~ PAIBORtot.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=binomial, nAGQ=1)
summary(mcutoff_noanx)


anova(mcutoff_noanx, mcutoff)

cm <- lmerCellMeans(mcutoff_noanx, divide=c("ECRavoid.c", "PAIBORtot.c"))

#convert to probability
cm$Cutoff_prob <- exp(cm$Cutoff_num)/(1+exp(cm$Cutoff_num))
cm$sehigh <- exp(cm$Cutoff_num + cm$se)/(1+exp(cm$Cutoff_num + cm$se))
cm$selow <- exp(cm$Cutoff_num - cm$se)/(1+exp(cm$Cutoff_num - cm$se))

pdf("figures/Cutoff by ECR Avoid and PAI-BOR.pdf", width=9, height=6)
ggplot(cm, aes(x=ECRavoid.c, y=Cutoff_prob, ymin=selow, ymax=sehigh, color=PAIBORtot.c)) + geom_pointrange(position=position_dodge(width=0.5), size=1.5) + 
    ggtitle("Probability of Alter Cutoff as a function of ECR Avoidance and PAI-BOR Sx")
dev.off()

mcutoff_close <- glmer(Cutoff_num ~ Close*PAIBORtot.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=binomial, nAGQ=1)

#yes, strong 3-way interaction
#interpretation: Those with low attachment avoidance and high BPD Sx tended to nominate more folks in their network who were distant and cutoff.
mcutoff_close <- glmer(Cutoff_num ~ Close*ECRavoid.c*PAIBORtot.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=binomial, nAGQ=1)
summary(mcutoff_close)


cm <- lmerCellMeans(mcutoff_close, divide=c("ECRavoid.c", "PAIBORtot.c"), n.cont=10)


mcutoff_close <- glmer(Cutoff_num ~ Close*ECRavoid.c + Close*ECRanx.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=binomial, nAGQ=1)


#
pdf("figures/Cutoff by ECR Avoid and PAI-BOR and Close.pdf", width=9, height=6)
ggplot(cm, aes(x=Close, y=Cutoff_prob, ymin=selow, ymax=sehigh, color=PAIBORtot.c)) + facet_wrap(~ECRavoid.c) + geom_pointrange(position=position_dodge(width=0.5), size=1.5) + 
    ggtitle("Probability of Alter Cutoff as a function of ECR Avoidance and PAI-BOR Sx")
dev.off()









####
graphlist[[2]]$actors$attach_total <- sapply(graphlist[[2]]$actors$attach_total, function(x) {if (is.na(x)) 10 else x })



centralization.degree(g2) #8000.0 matches
V(g2)


betweenness(g2, normalize=FALSE) #8000.0 matches

#plotting color in terms of self, partner, all others



V(g)$partner <- sapply(V(g)$partner, function(x) { 
      if(is.na(x)) { 
        "green" #self 
      } else if (x=="yes") {
        "pink"
      } else {
        "gray90"
      }
    })

g <- delete.edges(g, which(E(g)$weight == 1)) #remove "they have never met"

#plot(g, layout=layout.fruchterman.reingold, edge.width=E(g)$weight, vertex.size=V(g)$attach_total)
plot(g, layout=layout.fruchterman.reingold, edge.width=E(g)$weight, vertex.size=V(g)$attach_total,
    vertex.color=V(g)$partner)

#binarize edges at 1/2 versus 3/4/5
g2 <- delete.edges(g, which(E(g)$weight <= 2))
g2 <- remove.edge.attribute(g2, "weight")

pdf("figures/example_graph.pdf", width=11, height=8.5)
plot(g2, vertex.size=V(g2)$attach_total,
    vertex.color=V(g2)$partner)
dev.off()

## The opposite operation
get.data.frame(g2, what="vertices")
get.data.frame(g2, what="edges")

rglplot(g)


centralization.degree(g2)



edge.betweenness(g, e=E(g), directed = FALSE, weights = NULL)


pdf("figures/all_couple_graphs_6Oct2015.pdf", width=25, height=15)
for (i in 1:length(couple_graphs)) {
  #deidentify by removing last initials
  thiscouple <- couple_graphs[[i]]
  thiscouple$actors$name <- sub("(.*)_[A-z]+(\\.[0-9]+)", "\\1\\2", thiscouple$actors$name, perl=TRUE)
  thiscouple$edges_avgratings$a1 <- sub("(.*)_[A-z]+(\\.[0-9]+)", "\\1\\2", thiscouple$edges_avgratings$a1, perl=TRUE)
  thiscouple$edges_avgratings$a2 <- sub("(.*)_[A-z]+(\\.[0-9]+)", "\\1\\2", thiscouple$edges_avgratings$a2, perl=TRUE)
  g <- graph.data.frame(thiscouple$edges_avgratings, directed=FALSE, vertices=thiscouple$actors)
  #print(g, e=TRUE, v=TRUE)
  
  g <- delete.edges(g, which(E(g)$weight <= 2)) #remove "they have never met" and "not close at all"
  
  V(g)$role_color <- sapply(V(g)$role, function(x) { 
        if(x == "patient") { 
          "slateblue3" 
        } else if (x=="partner") {
          "indianred2"
        } else if (x=="patient_only") {
          "skyblue1"
        } else if (x=="partner_only") {
          "lightsalmon1"
        } else if (x=="shared") {
          "violet"
        } else { stop("unable to match") }
      })
  
  #g2 <- delete.edges(g, which(E(g)$weight <= 3))
  #g2 <- remove.edge.attribute(g2, "weight")

  #try to adjust vertex size
  co <- layout.fruchterman.reingold(g)
  
  plot(0, type="n", ann=FALSE, axes=FALSE, xlim=extendrange(co[,1]), 
      ylim=extendrange(co[,2]))
  plot(g, layout=co, rescale=FALSE, add=TRUE,
      vertex.shape="rectangle",
      vertex.size=(strwidth(V(g)$name) + strwidth("oo")) * 190,
      vertex.size2=strheight("I") * 2 * 230,
      edge.width=E(g)$weight*2, vertex.label.color="black",
      vertex.color=V(g)$role_color, vertex.label.family="sans", vertex.label.cex=2, main=names(couple_graphs)[i]
  )
#  plot(g, layout=layout.fruchterman.reingold, edge.width=E(g)$weight*2, vertex.label.color="black", #vertex.size=V(g)$attach_total,
#      vertex.color=V(g)$role_color, vertex.label.family="sans", vertex.label.cex=1.3, vertex.size=13, main=names(couple_graphs)[i],
#      vertex.shape="rectangle")
  
  #to control font size of tite: title("This is my first igraph",cex.main=3,col.main="green")
}

dev.off()






g <- graph.data.frame(couple_graphs[[1]]$edges, directed=FALSE, vertices=couple_graphs[[1]]$actors)
print(g, e=TRUE, v=TRUE)

g <- delete.edges(g, which(E(g)$weight <= 2)) #remove "they have never met" and "not close at all"

V(g)$role_color <- sapply(V(g)$role, function(x) { 
      if(x == "patient") { 
        "slateblue3" 
      } else if (x=="partner") {
        "indianred2"
      } else if (x=="patient_only") {
        "skyblue1"
      } else if (x=="partner_only") {
        "lightsalmon1"
      } else if (x=="shared") {
        "violet"
      } else { stop("unable to match") }
    })

g2 <- delete.edges(g, which(E(g)$weight <= 3))
g2 <- remove.edge.attribute(g2, "weight")

pdf("figures/8000_couple_SNA.pdf", width=25, height=25)
plot(g, layout=layout.fruchterman.reingold, edge.width=E(g)$weight, vertex.label.color="black", #vertex.size=V(g)$attach_total,
    vertex.color=V(g)$role_color, vertex.label.family="sans", vertex.label.cex=1.3, vertex.size=13)
dev.off()


#binarize edges at 1/2 versus 3/4/5
g2 <- delete.edges(g, which(E(g)$weight <= 2))
g2 <- remove.edge.attribute(g2, "weight")

pdf("figures/example_graph.pdf", width=11, height=8.5)
plot(g2, vertex.size=V(g2)$attach_total,
    vertex.color=V(g2)$partner)
dev.off()


#look at whether attachment items form a single factor
cor(sna$Close, sna$attach_total, use="pairwise.complete.obs")
msyn <- 'attach=~prxsk + sepds + sfhvn + secbs
    sfhvn ~~ secbs'
mcfa <- cfa(msyn, sna, missing="ML", estimator="MLR")

summary(mcfa, fit.measures=TRUE)
standardizedSolution(mcfa, type="std.all")
modificationIndices(mcfa)
