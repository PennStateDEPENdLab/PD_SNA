library(lme4)
library(igraph)
library(afex)
library(ggplot2)
setwd("/Users/michael/Tresors/PD_SNA/Couples_Baseline_SNA")
source("code/APIM_Functions.R")
source(file.path(getMainDir(), "Miscellaneous", "Global_Functions.R"))
load(file="data/SNA_Processed_27Jul2017.RData")
load(file="data/ego_graph_measures_27Jul2017.RData")
load(file="data/selfreports/couples_baseline_clinical_27Jul2017.RData")


#this script is for running analyses of ego-centered network metrics as a function of PD variables
cor(couples_baseline_clin[,c("IIP_communion", "IIP_agency", "IIP_elevation")], use="pairwise.complete.obs")
#analyses of alter/vertex attributes
mgraph <- merge(graph_agg, couples_baseline_clin, by=c("UsrID", "PTNUM"), all=FALSE)

predstocenter <- c("PAIBORtot", "IIP_agency", "IIP_communion", "IIP_elevation", "IIP_pd1", "IIP_pd2", "IIP_pd3", "ECRanx", "ECRavoid", "DASTotal", 
    "Victim", "Perp", "PsychAggV", "PsychAggP", "PhysAssV", "PhysAssP", "bordl_sidp", "narci_sidp", "antso_sidp", "avoid_sidp", "OPD_sidp", "depen_sidp",
    "nonarc", "nobpd", "pdtot")

#couples_baseline_clin <- f_centerPredictors(couples_baseline_clin, predstocenter, addsuffix=".c")
#couples_clin_wide <- f_centerPredictors(couples_clin_wide, as.vector(outer(predstocenter, c("0", "1"), paste, sep="_")), addsuffix=".c")

mgraph <- f_centerPredictors(mgraph, predstocenter, addsuffix=".c")

mgraph$nm3sqrt <- sqrt(mgraph$num_maximal_3cliques_3) #long tail here
mgraph$nm2sqrt <- sqrt(mgraph$num_maximal_3cliques_2) #long tail here

#analysis of graph-wide metrics

mgraph$DyadID <- as.numeric(substr(as.character(mgraph$UsrID), 5, 5))
mgraph_melt <- melt(mgraph[,sapply(mgraph, is.numeric)], id.vars=c("PTNUM", "DyadID", "UsrID"))
mgraph_wide <- dcast(mgraph_melt, PTNUM ~ variable + DyadID) #drops UsrID because not mentioned

graph_wide_short <- plyr::rename(mgraph_wide, c(IIP_communion_0="aff_0", IIP_communion_1="aff_1",
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
        number_optcommunity_weighted_2_0="ncom2_0", number_optcommunity_weighted_2_1="ncom2_1",
        largest_clique_size_2_0="lcliq2_0", largest_clique_size_2_1="lcliq2_1",
        largest_clique_size_3_0="lcliq3_0", largest_clique_size_3_1="lcliq3_1",
        number_optcommunity_weighted_3_0="ncom3_0", number_optcommunity_weighted_3_1="ncom3_1",
        avgsize_optcommunity_weighted_2_0="asize2_0", avgsize_optcommunity_weighted_2_1="asize2_1",
        avgsize_optcommunity_weighted_3_0="asize3_0", avgsize_optcommunity_weighted_3_1="asize3_1",
        modularity_optcommunity_weighted_2_0="modw2_0", modularity_optcommunity_weighted_2_1="modw2_1",
        modularity_optcommunity_weighted_3_0="modw3_0", modularity_optcommunity_weighted_3_1="modw3_1",
        modularity_optcommunity_binary_2_0="modb2_0", modularity_optcommunity_binary_2_1="modb2_1",
        modularity_optcommunity_binary_3_0="modb3_0", modularity_optcommunity_binary_3_1="modb3_1",
        yearsknown_assortativity_2_0="yka2_0", yearsknown_assortativity_2_1="yka2_1",
        romance_assortativity_2_0="ra2_0", romance_assortativity_2_1="ra2_1"
        ))

#mgraph$UsrID <- factor(mgraph$UsrID)
#mgraph$PTNUM <- factor(mgraph$PTNUM)

#now we're ready for APIM land...
options(error=NULL)


#MEASURES OF GRAPH CENTRALIZATION (PROPORTIONAL TO AVERAGE CENTRALITY OF NODES)
#basically normally distributed
#xx <- runAPIM(mgraph_wide, DV="centralization_degree_2", predictor="PAIBORtot", additional="", printall=FALSE)
#xx <- runAPIM(mgraph_wide, DV="centralization_degree_3", predictor="PAIBORtot", additional="", printall=FALSE)

xx <- runAPIM(mgraph_wide, DV="centralization_degree_2", predictor="bordl_sidp.c", additional="", printall=FALSE)
xx <- runAPIM(mgraph_wide, DV="centralization_degree_3", predictor="bordl_sidp.c", additional="", printall=FALSE)

xx <- runAPIM(mgraph_wide, DV="centralization_degree_2", predictors=c("bordl_sidp.c", "nobpd.c"), additional="", printall=FALSE)
xx <- runAPIM(mgraph_wide, DV="centralization_degree_3", predictors=c("bordl_sidp.c", "nobpd.c"), additional="", printall=FALSE)

xx <- runAPIM(mgraph_wide, DV="centralization_degree_2", predictors=c("narci_sidp.c", "nonarc.c"), additional="", printall=FALSE)

xx <- runAPIM(mgraph_wide, DV="centralization_degree_2", predictors=c("pdtot"), printall=FALSE)
xx <- runAPIM(mgraph_wide, DV="centralization_degree_3", predictors=c("pdtot"), printall=FALSE)

#weird that I get different p-values despite constraints. Try in Mplus
xx <- runAPIM(mgraph_wide, DV="centralization_degree_3", predictors=c("pdtot"), printall=FALSE)
tmp <- plyr::rename(mgraph_wide, c(centralization_degree_3_0="cd3_0", centralization_degree_3_1="cd3_1",
        centralization_degree_2_0="cd2_0", centralization_degree_2_1="cd2_1",
        centralization_eigenvector_2_0="ce2_0", centralization_eigenvector_2_1="ce2_1",
        centralization_eigenvector_3_0="ce3_0", centralization_eigenvector_3_1="ce3_1",
        bordl_sidp.c_0="bpdc_0", bordl_sidp.c_1="bpdc_1",
        nobpd.c_0="nobpdc_0", nobpd.c_1="nobpdc_1"))
output <- mplusAPIM(tmp, "cd3", "pdtot")
output <- mplusAPIM(tmp, "cd2", "pdtot")

#okay, this is a match
output <- mplusAPIM(tmp, "cd2", c("bpdc", "nobpdc"))
output <- mplusAPIM(tmp, "cd3", c("bpdc", "nobpdc"))

output <- mplusAPIM(tmp, "cd2", c("nobpdc")) 
output <- mplusAPIM(tmp, "cd3", c("nobpdc"))#weird reversal of actor/partner

xx <- runAPIM(mgraph_wide, DV="centralization_degree_2", predictors=c("narci_sidp.c", "nonarc.c"), printall=TRUE) #null
xx <- runAPIM(mgraph_wide, DV="centralization_degree_3", predictors=c("narci_sidp.c", "nonarc.c"), printall=FALSE) #null

histogram(mgraph_wide$centralization_degree_2_0)
xx <- runAPIM(mgraph_wide, DV="centralization_degree_2", predictors=c("IIP_agency.c", "IIP_communion.c"), additional="", printall=FALSE)
xx <- runAPIM(mgraph_wide, DV="centralization_degree_3", predictors=c("IIP_agency.c", "IIP_communion.c"), additional="", printall=FALSE) #some agency effects here

#indeed, the parameter estimate here is a very close match to APIM with indistinguishable dyads 
summary(m <- lmer(centralization_degree_2 ~ pdtot.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(centralization_degree_3 ~ pdtot.c + DyadID + (1 | PTNUM), mgraph))

summary(m <- lmer(centralization_degree_2 ~ bordl_sidp.c + nobpd.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(centralization_degree_3 ~ bordl_sidp.c + nobpd.c + DyadID + (1 | PTNUM), mgraph))

#eigenvector centralization
xx <- runAPIM(mgraph_wide, DV="centralization_eigenvector_2", predictors=c("pdtot"))
xx <- runAPIM(mgraph_wide, DV="centralization_eigenvector_3", predictors=c("pdtot"))

xx <- runAPIM(mgraph_wide, DV="centralization_eigenvector_2", predictors=c("nobpd.c"))

xx <- runAPIM(mgraph_wide, DV="centralization_eigenvector_2", predictors=c("bordl_sidp.c", "nobpd.c"))
xx <- runAPIM(mgraph_wide, DV="centralization_eigenvector_3", predictors=c("bordl_sidp.c", "nobpd.c"))

xx <- runAPIM(mgraph_wide, DV="centralization_eigenvector_2", predictors=c("narci_sidp.c", "nonarc.c"))
xx <- runAPIM(mgraph_wide, DV="centralization_eigenvector_3", predictors=c("narci_sidp.c", "nonarc.c"))

output <- mplusAPIM(tmp, "ce2", "pdtot")
output <- mplusAPIM(tmp, "ce3", "pdtot")

output <- mplusAPIM(tmp, "ce2", c("bpdc", "nobpdc"))
output <- mplusAPIM(tmp, "ce3", c("bpdc", "nobpdc"))

#effects here, but hard to interpret
xx <- runAPIM(mgraph_wide, DV="centralization_eigenvector_2", predictors=c("IIP_agency.c", "IIP_communion.c"))
xx <- runAPIM(mgraph_wide, DV="centralization_eigenvector_3", predictors=c("IIP_agency.c", "IIP_communion.c"))

xx <- runAPIM(mgraph_wide, DV="centralization_eigenvector_2", predictors=c("IIP_pd1.c", "IIP_pd2.c", "IIP_pd3.c"))
xx <- runAPIM(mgraph_wide, DV="centralization_eigenvector_2", predictors=c("IIP_pd2.c"))
xx <- runAPIM(mgraph_wide, DV="centralization_eigenvector_3", predictors=c("IIP_pd2.c"), printall=TRUE)


summary(m <- lmer(centralization_eigenvector_2 ~ PAIBORtot.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(centralization_eigenvector_2 ~ PAIBORtot.c*IIP_agency.c*IIP_communion.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(centralization_eigenvector_2 ~ IIP_pd1.c + (1 | PTNUM), mgraph))
summary(m <- lmer(centralization_eigenvector_2 ~ IIP_pd2.c + (1 | PTNUM), mgraph))
summary(m <- lmer(centralization_eigenvector_2 ~ IIP_pd3.c + (1 | PTNUM), mgraph))

mixed(centralization_eigenvector_2 ~ PAIBORtot.c*IIP_agency.c*IIP_communion.c + (1 | PTNUM), data=mgraph, method="KR")

mixed(centralization_eigenvector_2 ~ PAIBORtot.c*IIP_pd1.c*IIP_pd2.c*IIP_pd3.c + (1 | PTNUM), data=mgraph, method="KR")
mixed(centralization_eigenvector_2 ~ PAIBORtot.c*IIP_pd2.c + (1 | PTNUM), data=mgraph, method="KR")
summary(m)

cm <- lmerCellMeans(m, divide=c("IIP_agency.c", "IIP_communion.c"), n.cont=10)
ggplot(cm, aes(x=PAIBORtot.c, y=centralization_eigenvector_2, color=IIP_agency.c, shape=IIP_communion.c)) + geom_point() + geom_line()

m <- lmer(centralization_eigenvector_3 ~ PAIBORtot + (1 | PTNUM), mgraph)
summary(m)



##MEASURES OF ASSORTATIVITY


#angry assortativity
#summary(m <- lmer(angry_assortativity_2 ~ PAIBORtot.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(angry_assortativity_2 ~ IIP_agency.c + IIP_communion.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(angry_assortativity_2 ~ IIP_pd1.c + IIP_pd2.c + IIP_pd3.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(angry_assortativity_2 ~ narci_sidp.c + nonarc.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(angry_assortativity_2 ~ bordl_sidp.c + nobpd.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(angry_assortativity_2 ~ nobpd.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(angry_assortativity_2 ~ pdtot.c + DyadID + (1 | PTNUM), mgraph))
(m <- mixed(angry_assortativity_2 ~ nobpd.c + DyadID + (1 | PTNUM), mgraph)) #p = .04
(m <- mixed(angry_assortativity_2 ~ bordl_sidp.c + nobpd.c + DyadID + (1 | PTNUM), mgraph)) #p = .009

xx <- runAPIM(mgraph_wide, DV="angry_assortativity_2", predictors=c("nobpd.c"))
xx <- runAPIM(mgraph_wide, DV="angry_assortativity_2", predictors=c("bordl_sidp.c", "nobpd.c"))
xx <- runAPIM(mgraph_wide, DV="angry_assortativity_2", predictors=c("narci_sidp.c", "nonarc.c"))
xx <- runAPIM(mgraph_wide, DV="angry_assortativity_3", predictors=c("bordl_sidp.c", "nobpd.c"))
xx <- runAPIM(mgraph_wide, DV="angry_assortativity_2", predictors=c("pdtot"))
xx <- runAPIM(mgraph_wide, DV="angry_assortativity_2", predictors=c("IIP_agency.c", "IIP_communion.c")) #null

#summary(m <- lmer(angry_assortativity_2 ~ PAIBORtot.c + ECRanx.c + ECRavoid.c + (1 | PTNUM), mgraph)) #no attachment stuff for now
#summary(m <- lmer(angry_assortativity_2 ~ ECRanx.c + (1 | PTNUM), mgraph))
#summary(m <- lmer(angry_assortativity_2 ~ ECRavoid.c + (1 | PTNUM), mgraph))

#happy assortativity: no effects
#summary(m <- lmer(happy_assortativity_2 ~ PAIBORtot.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(happy_assortativity_2 ~ IIP_agency.c + IIP_communion.c + DyadID + (1 | PTNUM), mgraph))
#summary(m <- lmer(happy_assortativity_2 ~ IIP_pd1.c + IIP_pd2.c + IIP_pd3.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(happy_assortativity_2 ~ narci_sidp.c + nonarc.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(happy_assortativity_2 ~ bordl_sidp.c + nobpd.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(happy_assortativity_2 ~ nobpd.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(happy_assortativity_2 ~ pdtot.c + DyadID + (1 | PTNUM), mgraph))
(m <- mixed(happy_assortativity_2 ~ nobpd.c + DyadID + (1 | PTNUM), mgraph)) #p = .04
(m <- mixed(happy_assortativity_2 ~ bordl_sidp.c + nobpd.c + DyadID + (1 | PTNUM), mgraph)) #p = .009

xx <- runAPIM(mgraph_wide, DV="happy_assortativity_2", predictors=c("nobpd.c"))
xx <- runAPIM(mgraph_wide, DV="happy_assortativity_2", predictors=c("bordl_sidp.c", "nobpd.c"))
xx <- runAPIM(mgraph_wide, DV="happy_assortativity_2", predictors=c("pdtot"))
xx <- runAPIM(mgraph_wide, DV="happy_assortativity_2", predictors=c("IIP_agency.c", "IIP_communion.c"))

#attach assortativity: no effects
#summary(m <- lmer(attach_assortativity_2 ~ PAIBORtot.c + DyadID +  (1 | PTNUM), mgraph))
summary(m <- lmer(attach_assortativity_2 ~ IIP_agency.c + IIP_communion.c + DyadID + (1 | PTNUM), mgraph))
#summary(m <- lmer(attach_assortativity_2 ~ IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(attach_assortativity_2 ~ pdtot + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(attach_assortativity_2 ~ bordl_sidp.c + nobpd.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(attach_assortativity_2 ~ narci_sidp.c + nonarc.c + DyadID + (1 | PTNUM), mgraph))

#summary(m <- lmer(attach_assortativity_2 ~ PAIBORtot.c*ECRanx.c*ECRavoid.c + (1 | PTNUM), mgraph))
#mixed(attach_assortativity_2 ~ PAIBORtot.c*ECRanx.c*ECRavoid.c + (1 | PTNUM), mgraph)
#summary(m <- lmer(attach_assortativity_2 ~ PAIBORtot.c*ECRanx.c + (1 | PTNUM), mgraph))

#cm <- lmerCellMeans(m, divide=c("ECRanx.c"))
#ggplot(cm, aes(x=PAIBORtot.c, y=attach_assortativity_2, color=ECRanx.c)) + geom_line()

#summary(m <- lmer(attach_assortativity_2 ~ ECRanx.c + (1 | PTNUM), mgraph))
#summary(m <- lmer(attach_assortativity_2 ~ ECRavoid.c + (1 | PTNUM), mgraph))
#hist(mgraph$attach_assortativity_2)

#years known assortativity (people you've known longer tend to have similar connection strengths)
summary(m <- lmer(yearsknown_assortativity_2 ~ PAIBORtot.c + DyadID + (1 | PTNUM), mgraph))
#summary(m <- lmer(yearsknown_assortativity_2 ~ ECRanx.c + (1 | PTNUM), mgraph))
#summary(m <- lmer(yearsknown_assortativity_2 ~ ECRavoid.c + (1 | PTNUM), mgraph))

#so being more domineering and overly nurturant

summary(m <- lmer(yearsknown_assortativity_2 ~ IIP_agency.c + IIP_communion.c + DyadID + (1 | PTNUM), mgraph))
library(sjPlot)

sjp.lmer(m, type = "fe.slope", vars = c("IIP_agency.c", "IIP_communion.c"))

#is there an outlier?
mgraph$IIP_communion_wins <- winsor(mgraph$IIP_communion, trim=0.01) #trim a few outliers
mgraph$IIP_agency_wins <- winsor(mgraph$IIP_agency, trim=0.01) #trim a few outliers
summary(m <- lmer(yearsknown_assortativity_2 ~ IIP_communion_wins + IIP_agency_wins + DyadID + (1 | PTNUM), mgraph)) #holds up
sjp.lmer(m, type = "fe.slope", vars = c("IIP_agency_wins", "IIP_communion_wins"))

#PD scales don't pan out here
summary(m <- lmer(yearsknown_assortativity_2 ~ IIP_pd1.c + IIP_pd2.c + IIP_pd3.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(yearsknown_assortativity_2 ~ IIP_pd1.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(yearsknown_assortativity_2 ~ IIP_pd2.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(yearsknown_assortativity_2 ~ IIP_pd3.c + DyadID + (1 | PTNUM), mgraph))


summary(m <- lmer(yearsknown_assortativity_2 ~ IIP_agency.c * IIP_communion.c * DyadID + (1 | PTNUM), mgraph)) #no interactions
#summary(m <- lmer(yearsknown_assortativity_2 ~ IIP_PD1.c + DyadID + (1 | PTNUM), mgraph))
#summary(m <- lmer(yearsknown_assortativity_2 ~ IIP_PD2.c + DyadID + (1 | PTNUM), mgraph))
#summary(m <- lmer(yearsknown_assortativity_2 ~ IIP_PD3.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(yearsknown_assortativity_2 ~ bordl_sidp.c + nobpd.c + DyadID + (1 | PTNUM), mgraph))
(m <- mixed(yearsknown_assortativity_2 ~ bordl_sidp.c + nobpd.c + DyadID + (1 | PTNUM), mgraph)) #p=.01 for BPD
summary(m <- lmer(yearsknown_assortativity_2 ~ narci_sidp.c + nonarc.c + DyadID + (1 | PTNUM), mgraph))

#try agency + communion with BPD
(m <- mixed(yearsknown_assortativity_2 ~ bordl_sidp.c + nobpd.c + IIP_communion.c + IIP_agency.c + DyadID + (1 | PTNUM), mgraph)) #p=.01 for BPD
summary(m <- lmer(yearsknown_assortativity_2 ~ bordl_sidp.c + nobpd.c + IIP_communion.c + IIP_agency.c + DyadID + (1 | PTNUM), mgraph)) #p=.01 for BPD

#looks like marginal BPD negative effect
#sig negative communion actor effect
#sig negative agency actor effect
#sig negative agency partner effect 
xx <- runAPIM(mgraph_wide, DV="yearsknown_assortativity_2", predictors=c("bordl_sidp.c", "nobpd.c", "IIP_communion.c", "IIP_agency.c"))
output <- mplusAPIM(graph_wide_short, "yka2", c("bpd", "nobpd"), count=FALSE)
output <- mplusAPIM(graph_wide_short, "yka2", c("pdtot"), count=FALSE)
mean(mgraph$yearsknown_assortativity_2)
sd(mgraph$yearsknown_assortativity_2)
output <- mplusAPIM(graph_wide_short, "yka2", c("aff", "dom"), count=FALSE) #no effect for IIPe
output <- mplusAPIM(graph_wide_short, "yka2", c("aff", "dom", "IIPe"), count=FALSE) #no effect for IIPe

hist(mgraph$yearsknown_assortativity_2)

#so, those low in BPD and IIP scores tend to have a positive association between years known and alter connections. This is disrupted in PD 

#cutoff assortativity
#the distribution is too crazy, likely because some people have no cutoffs (leading to NaN) and others have many (leading to 1.0)

#emotion assortativity (also very peaked and with weird skew. skip for now
histogram(mgraph$emotion_assortativity_2)
summary(m <- lmer(emotion_assortativity_2 ~ IIP_agency.c + IIP_communion.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(emotion_assortativity_2 ~ pdtot.c + DyadID + (1 | PTNUM), mgraph))

#past romance assortativity.
histogram(mgraph$romance_assortativity_2)

#those with higher communion (more nurturant) style tend to organize their networks more around past romantic partners

summary(m <- lmer(romance_assortativity_2 ~ IIP_agency.c + IIP_communion.c + IIP_elevation.c + DyadID + (1 | PTNUM), mgraph))
sjp.lmer(m, type = "fe.slope", vars = c("IIP_communion.c"))
summary(m <- lmer(romance_assortativity_2 ~ pdtot.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(romance_assortativity_2 ~ bordl_sidp.c + nobpd.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(romance_assortativity_2 ~ narci_sidp.c + nonarc.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(romance_assortativity_2 ~ pdtot.c + DyadID + (1 | PTNUM), mgraph))


summary(m <- lmer(romance_assortativity_2 ~ IIP_agency.c + IIP_communion.c + IIP_elevation.c + DyadID + (1 | PTNUM), mgraph))
output <- mplusAPIM(graph_wide_short, "ra2", c("dom", "aff"), count=FALSE)


#largest clique size
summary(m <- lmer(largest_clique_size_2 ~ bordl_sidp.c + (1 | PTNUM), mgraph))
summary(m <- lmer(largest_clique_size_3 ~ bordl_sidp.c + (1 | PTNUM), mgraph))

summary(m <- lmer(largest_clique_size_2 ~ bordl_sidp.c + nobpd.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(largest_clique_size_3 ~ bordl_sidp.c + nobpd.c + DyadID + (1 | PTNUM), mgraph))

summary(m <- lmer(largest_clique_size_2 ~ narci_sidp.c + nonarc.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(largest_clique_size_3 ~ narci_sidp.c + nonarc.c + DyadID + (1 | PTNUM), mgraph))

summary(m <- lmer(largest_clique_size_2 ~ pdtot.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(largest_clique_size_3 ~ pdtot.c + DyadID + (1 | PTNUM), mgraph))

#marginal
xx <- runAPIM(mgraph_wide, DV="largest_clique_size_2", predictors=c("pdtot.c"))
xx <- runAPIM(mgraph_wide, DV="largest_clique_size_3", predictors=c("pdtot.c"))

#non-BPD PD Sx associated with smaller largest cliques. But these effects become marginal if BPD Sx included in PD sum
xx <- runAPIM(mgraph_wide, DV="largest_clique_size_2", predictors=c("bordl_sidp.c", "nobpd.c"))
xx <- runAPIM(mgraph_wide, DV="largest_clique_size_3", predictors=c("bordl_sidp.c", "nobpd.c"))




output <- mplusAPIM(graph_wide_short, "lcliq2", c("bpd", "nobpd"), count=FALSE)
output <- mplusAPIM(graph_wide_short, "lcliq2", c("bpd", "nobpd"), count=FALSE)

xx <- runAPIM(mgraph_wide, DV="largest_clique_size_2", predictors=c("nobpd.c"))
xx <- runAPIM(mgraph_wide, DV="largest_clique_size_3", predictors=c("nobpd.c"))

xx <- runAPIM(mgraph_wide, DV="largest_clique_size_2", predictors=c("narci_sidp.c", "nonarc.c"))
xx <- runAPIM(mgraph_wide, DV="largest_clique_size_3", predictors=c("narci_sidp.c", "nonarc.c")) #non-narc Sx same as non-BPD above

#large partner effect here: more problems with partner's affiliation associated with own smaller clique size
histogram(sqrt(mgraph$largest_clique_size_2))
histogram(mgraph$largest_clique_size_3)
summary(m <- lmer(sqrt(largest_clique_size_2) ~ IIP_agency.c + IIP_communion.c + DyadID + (1 | PTNUM), mgraph))
summary(m <- lmer(sqrt(largest_clique_size_3) ~ IIP_agency.c + IIP_communion.c + DyadID + (1 | PTNUM), mgraph))

summary(m <- glmer(largest_clique_size_2 ~ IIP_agency.c + IIP_communion.c + DyadID + (1 | PTNUM), mgraph, family=poisson))
summary(m <- glmer(largest_clique_size_3 ~ IIP_agency.c + IIP_communion.c + DyadID + (1 | PTNUM), mgraph, family=poisson))


xx <- runAPIM(mgraph_wide, DV="largest_clique_size_2", predictors=c("IIP_agency.c", "IIP_communion.c"))
xx <- runAPIM(mgraph_wide, DV="largest_clique_size_3", predictors=c("IIP_agency.c", "IIP_communion.c")) #same direction, but marginal
xx <- runAPIM(mgraph_wide, DV="largest_clique_size_2", predictors=c("IIP_elevation.c", "IIP_agency.c", "IIP_communion.c")) #more 
xx <- runAPIM(mgraph_wide, DV="largest_clique_size_3", predictors=c("IIP_elevation.c", "IIP_agency.c", "IIP_communion.c")) #same direction, but marginal

#summary(m <- lmer(largest_clique_size_2 ~ PAIBORtot.c + DyadID + (1 | PTNUM), mgraph))
#summary(m <- lmer(largest_clique_size_3 ~ PAIBORtot.c + DyadID + (1 | PTNUM), mgraph))

#largest clique is smaller with BPD Sx. NB: does not control for patient/partner, so invalid!
#mixed(largest_clique_size_2 ~ PAIBORtot.c + (1 | PTNUM), mgraph, method="KR")
#mixed(largest_clique_size_3 ~ PAIBORtot.c + (1 | PTNUM), mgraph, method="KR")

#mixed(largest_clique_size_2 ~ PAIBORtot.c*ECRanx.c + (1 | PTNUM), mgraph, method="KR")
#mixed(largest_clique_size_3 ~ PAIBORtot.c*ECRanx.c + (1 | PTNUM), mgraph, method="KR")

#mixed(largest_clique_size_3 ~ PAIBORtot.c*ECRanx.c + (1 | PTNUM), mgraph, method="KR")
#mixed(largest_clique_size_3 ~ PAIBORtot.c*ECRavoid.c + (1 | PTNUM), mgraph, method="KR")

summary(m <- lmer(largest_clique_size_2 ~ IIP_agency.c + IIP_communion.c + IIP_elevation.c + DyadID + (1 | PTNUM), mgraph))
#summary(m <- lmer(largest_clique_size_2 ~ IIP_PD1.c + (1 | PTNUM), mgraph))
#summary(m <- lmer(largest_clique_size_2 ~ IIP_PD2.c + (1 | PTNUM), mgraph))
#summary(m <- lmer(largest_clique_size_2 ~ IIP_PD3.c + (1 | PTNUM), mgraph))
#summary(m <- lmer(largest_clique_size_2 ~ IIP_PD1.c*IIP_PD2.c*IIP_PD3.c + (1 | PTNUM), mgraph))

#number of maximal cliques of size 3 or larger
#more maximal cliques in those with elevated affiliation problems (enmeshment)

#summary(lmer(num_maximal_3cliques_2 ~ PAIBORtot.c + (1 | PTNUM), mgraph))
#summary(lmer(num_maximal_3cliques_2 ~ ECRanx.c +  (1 | PTNUM), mgraph)) #IIP_communion.c +
histogram(mgraph$num_maximal_3cliques_2) #skewed
histogram(mgraph$num_maximal_3cliques_3)
summary(lmer(num_maximal_3cliques_2 ~ IIP_communion.c + DyadID + (1 | PTNUM), mgraph)) #large effect
summary(lmer(sqrt(num_maximal_3cliques_2) ~ IIP_communion.c + DyadID + (1 | PTNUM), mgraph)) #large effect

summary(lmer(num_maximal_3cliques_2 ~ IIP_agency.c + DyadID + (1 | PTNUM), mgraph))
summary(lmer(sqrt(num_maximal_3cliques_2) ~ IIP_agency.c + DyadID + (1 | PTNUM), mgraph))

summary(lmer(num_maximal_3cliques_2 ~ IIP_elevation.c + DyadID + (1 | PTNUM), mgraph))

#summary(m <- lmer(num_maximal_3cliques_2 ~ IIP_agency.c*IIP_communion.c*PAIBORtot.c + (1 | PTNUM), mgraph))
#mixed(num_maximal_3cliques_2 ~ IIP_agency.c*IIP_communion.c*PAIBORtot.c + (1 | PTNUM), mgraph, method="KR")
#cm <- lmerCellMeans(m, divide=c("IIP_agency.c", "PAIBORtot.c"))
#ggplot(cm, aes(x=IIP_communion.c, y=num_maximal_3cliques_2, color=IIP_agency.c, shape=PAIBORtot.c)) + geom_line() + geom_point(size=3)

#holds at 3 threshold
histogram(mgraph$num_maximal_3cliques_3)
histogram(mgraph$num_maximal_3cliques_2)
summary(lmer(num_maximal_3cliques_3 ~ IIP_communion.c + DyadID + (1 | PTNUM), mgraph))
mixed(num_maximal_3cliques_3 ~ IIP_communion.c + DyadID + (1 | PTNUM), mgraph, method="KR")
mixed(nm3sqrt ~ IIP_communion.c + DyadID + (1 | PTNUM), mgraph, method="KR")
mixed(nm2sqrt ~ IIP_communion.c + DyadID + (1 | PTNUM), mgraph, method="KR")

xx <- runAPIM(mgraph_wide, DV="nm2sqrt", predictors=c("IIP_agency.c", "IIP_communion.c")) #all on the actor side
xx <- runAPIM(mgraph_wide, DV="nm2sqrt", predictors=c("IIP_agency.c", "IIP_communion.c")) #all on the actor side
#looks like pfree may be better
summary(xx$pfree, fit.measures=TRUE) #yes, this is better, but really, the actor effect is all that matters anyhow (all others NS)
xx <- runAPIM(mgraph_wide, DV="num_maximal_3cliques_3", predictors=c("IIP_agency.c", "IIP_communion.c")) #same direction, but marginal
summary(xx$pfree, fit.measures=TRUE) #yes, this is better, but really, the actor effect is all that matters anyhow (all others NS)

#summary(lmer(num_maximal_3cliques_3 ~ PAIBORtot.c + (1 | PTNUM), mgraph))
#summary(lmer(num_maximal_3cliques_3 ~ IIP_communion.c*PAIBORtot.c + (1 | PTNUM), mgraph))

#size of maximal cliques
histogram(mgraph$avgsize_maximal_3cliques_2)
#histogram(sqrt(mgraph$avgsize_maximal_3cliques_2))
summary(lmer(avgsize_maximal_3cliques_2 ~ IIP_communion.c + (1 | PTNUM), mgraph))
summary(lmer(avgsize_maximal_3cliques_3 ~ IIP_communion.c + (1 | PTNUM), mgraph))
summary(lmer(avgsize_maximal_3cliques_2 ~ IIP_agency.c + (1 | PTNUM), mgraph))

summary(lmer(avgsize_maximal_3cliques_2 ~ PAIBORtot.c + DyadID + (1 | PTNUM), mgraph))
summary(lmer(avgsize_maximal_3cliques_3 ~ PAIBORtot.c + DyadID + (1 | PTNUM), mgraph)) #holds at 3 threshold

#size of maximal cliques also smaller in those with high attachment anxiety (only at higher threshold)
#summary(lmer(avgsize_maximal_3cliques_2 ~ ECRanx.c + (1 | PTNUM), mgraph))
#summary(lmer(avgsize_maximal_3cliques_3 ~ ECRanx.c + (1 | PTNUM), mgraph))
#summary(lmer(avgsize_maximal_3cliques_2 ~ ECRanx.c*PAIBORtot.c + (1 | PTNUM), mgraph))
#summary(lmer(avgsize_maximal_3cliques_3 ~ ECRanx.c*PAIBORtot.c + (1 | PTNUM), mgraph)) #BPD diff explained/goes away when controlling for ECR anx
#summary(lmer(avgsize_maximal_3cliques_3 ~ ECRanx.c + PAIBORtot.c + (1 | PTNUM), mgraph)) #BPD diff explained/goes away when controlling for ECR anx
#summary(lmer(avgsize_maximal_3cliques_2 ~ PAIBORtot.c*IIP_communion.c + (1 | PTNUM), mgraph))
#summary(lmer(avgsize_maximal_3cliques_3 ~ PAIBORtot.c*IIP_communion.c + (1 | PTNUM), mgraph)) #not explained by communion

#higher number of communities in BPD
#summary(lmer(number_optcommunity_weighted_2 ~ PAIBORtot.c + DyadID + (1 | PTNUM), mgraph))
#summary(lmer(number_optcommunity_weighted_3 ~ PAIBORtot.c + DyadID + (1 | PTNUM), mgraph))

histogram(~number_optcommunity_weighted_2 | DyadID, mgraph)

#xx <- runAPIM(mgraph_wide, DV="number_optcommunity_weighted_2", predictor="PAIBORtot", additional="", printall=FALSE)
#xx <- runAPIM(mgraph_wide, DV="number_optcommunity_weighted_3", predictor="PAIBORtot", additional="", printall=FALSE)

xx <- runAPIM(mgraph_wide, DV="number_optcommunity_weighted_2", predictor=c("bordl_sidp", "nobpd"), printall=FALSE)
xx <- runAPIM(mgraph_wide, DV="number_optcommunity_weighted_3", predictor=c("bordl_sidp", "nobpd"), printall=FALSE)
xx <- runAPIM(mgraph_wide, DV="number_optcommunity_weighted_2", predictor=c("pdtot"), printall=FALSE)
xx <- runAPIM(mgraph_wide, DV="number_optcommunity_weighted_3", predictor=c("pdtot"), printall=FALSE)
summary(xx$pfree)
xx <- runAPIM(mgraph_wide, DV="number_optcommunity_weighted_2", predictor=c("narci_sidp", "nonarc"), additional="", printall=FALSE)

xx <- runAPIM(mgraph_wide, DV="number_optcommunity_weighted_2", predictor=c("IIP_agency"), additional="", printall=FALSE)
xx <- runAPIM(mgraph_wide, DV="number_optcommunity_weighted_2", predictor=c("IIP_agency"), additional="", printall=FALSE)

#model as Poisson
output <- mplusAPIM(graph_wide_short, "ncom2", "pdtot", count=TRUE)
output <- mplusAPIM(graph_wide_short, "ncom3", "pdtot", count=TRUE)

output <- mplusAPIM(graph_wide_short, "ncom2", c("bpd", "nobpd"), count=TRUE)
output <- mplusAPIM(graph_wide_short, "ncom2", c("nar", "nonarc"), count=TRUE)
output <- mplusAPIM(graph_wide_short, "ncom2", c("dom", "aff"), count=TRUE)

#some discrepancy between thresholds here... one is actor, one is partner
with(mgraph_wide, cor(number_optcommunity_weighted_3_0, number_optcommunity_weighted_3_1))
with(mgraph_wide, cor(number_optcommunity_weighted_2_0, number_optcommunity_weighted_2_1))

#should correspond with smaller sizes... YES!
output <- mplusAPIM(graph_wide_short, "asize2", c("bpd", "nobpd"))
output <- mplusAPIM(graph_wide_short, "asize2", c("nar", "nonarc"))
output <- mplusAPIM(graph_wide_short, "asize2", c("pdtot"))
output <- mplusAPIM(graph_wide_short, "asize2", c("dom", "aff"))


output <- mplusAPIM(graph_wide_short, "ncom2", c("dom", "aff"), count=TRUE)

#summary(lmer(avgsize_optcommunity_weighted_2 ~ PAIBORtot.c + DyadID + (1 | PTNUM), mgraph))
#summary(lmer(avgsize_optcommunity_weighted_3 ~ PAIBORtot.c + DyadID + (1 | PTNUM), mgraph))

xx <- runAPIM(mgraph_wide, DV="avgsize_optcommunity_weighted_2", predictor="PAIBORtot", additional="", printall=FALSE)
xx <- runAPIM(mgraph_wide, DV="avgsize_optcommunity_weighted_3", predictor="PAIBORtot", additional="", printall=FALSE)

#weakly positive greater modularity in BPD (not significant)
#modularity has a nice normal distribution
summary(lmer(modularity_optcommunity_weighted_2 ~ PAIBORtot.c + (1 | PTNUM), mgraph))
summary(lmer(modularity_optcommunity_weighted_3 ~ PAIBORtot.c + (1 | PTNUM), mgraph))
summary(lmer(modularity_optcommunity_binary_2 ~ PAIBORtot.c + (1 | PTNUM), mgraph))
summary(lmer(modularity_optcommunity_binary_3 ~ PAIBORtot.c + (1 | PTNUM), mgraph))

xx <- runAPIM(mgraph_wide, DV="modularity_optcommunity_weighted_2", predictor="pdtot", printall=FALSE)

xx <- runAPIM(mgraph_wide, DV="modularity_optcommunity_weighted_2", predictor="bordl_sidp", printall=FALSE)
xx <- runAPIM(mgraph_wide, DV="modularity_optcommunity_weighted_3", predictor="bordl_sidp", printall=FALSE)

output <- mplusAPIM(graph_wide_short, "modw2", "pdtot")
output <- mplusAPIM(graph_wide_short, "modw2", c("bpd", "nobpd"))
output <- mplusAPIM(graph_wide_short, "modb2", c("bpd", "nobpd")) #lots of effect here
output <- mplusAPIM(graph_wide_short, "modw3", c("bpd", "nobpd")) #lots of effect here

#BPD Sx suppress modularity; OPD symptoms enhance.


output <- mplusAPIM(graph_wide_short, "modb2", c("nar", "nonarc"))

output <- mplusAPIM(graph_wide_short, "modw2", c("aff", "dom", "IIPe")) #both dom and aff suppress modularity
output <- mplusAPIM(graph_wide_short, "modw2", c("aff", "dom")) #both dom and aff suppress modularity

#density differences?
summary(lmer(density_2 ~ PAIBORtot.c + DyadID + (1 | PTNUM), mgraph))
summary(lmer(density_3 ~ PAIBORtot.c + DyadID + (1 | PTNUM), mgraph))

xx <- runAPIM(mgraph_wide, DV="density_2", predictor="bordl_sidp", additional="", printall=TRUE)
summary(lmer(density_2 ~ PAIBORtot.c + (1 | PTNUM), mgraph))

#not really explained by ECR
#summary(lmer(density_2 ~ PAIBORtot.c*ECRanx.c*ECRavoid.c + (1 | PTNUM), mgraph))
#summary(lmer(density_3 ~ PAIBORtot.c*ECRanx.c*ECRavoid.c + (1 | PTNUM), mgraph))
#summary(lmer(density_2 ~ PAIBORtot.c*IIP_communion.c*IIP_agency.c + (1 | PTNUM), mgraph))
#summary(lmer(density_3 ~ PAIBORtot.c*IIP_communion.c + PAIBORtot.c*IIP_agency.c + (1 | PTNUM), mgraph))

#transitivity is weaker in BPD
summary(lmer(transitivity_2 ~ PAIBORtot.c + DyadID + (1 | PTNUM), mgraph))
summary(lmer(transitivity_3 ~ PAIBORtot.c + DyadID + (1 | PTNUM), mgraph))

xx <- runAPIM(mgraph_wide, DV="transitivity_2", predictor=c("bordl_sidp", "nobpd"), additional="", printall=FALSE)
xx <- runAPIM(mgraph_wide, DV="transitivity_2", predictor=c("pdtot"), additional="", printall=FALSE)
xx <- runAPIM(mgraph_wide, DV="transitivity_2", predictor=c("IIP_agency", "IIP_communion", "IIP_elevation"), additional="", printall=FALSE)

summary(lmer(transitivity_2 ~ bordl_sidp.c + DyadID + (1 | PTNUM), mgraph))
summary(lmer(transitivity_3 ~ bordl_sidp.c + DyadID + (1 | PTNUM), mgraph))

summary(lmer(transitivity_2 ~ PAIBORtot.c*IIP_communion.c + PAIBORtot.c*IIP_agency.c + (1 | PTNUM), mgraph))
summary(lmer(transitivity_3 ~ PAIBORtot.c*IIP_communion.c + PAIBORtot.c*IIP_agency.c + (1 | PTNUM), mgraph))

#average path length not related to BPD
summary(lmer(avgpath_2 ~ PAIBORtot.c + (1 | PTNUM), mgraph))
summary(lmer(avgpath_3 ~ PAIBORtot.c + (1 | PTNUM), mgraph))

#but is related to communion
#so more communality means there are fewer hops to connect any 2 people
summary(lmer(avgpath_2 ~ IIP_communion.c*IIP_agency.c + (1 | PTNUM), mgraph))
summary(lmer(avgpath_2 ~ IIP_communion.c + (1 | PTNUM), mgraph))
summary(lmer(avgpath_2 ~ IIP_agency.c + (1 | PTNUM), mgraph))

xx <- runAPIM(mgraph_wide, DV="density_2", predictor="IIP_agency", additional="", printall=TRUE)
summary(lmer(density_2 ~ PAIBORtot.c + (1 | PTNUM), mgraph))


#check small worldness
#have some infinite values for this...
#set to highest observed
mgraph$smallworld_2[which(mgraph$smallworld_2 == Inf)] <- max(mgraph$smallworld_2[!mgraph$smallworld_2==Inf])
mgraph$smallworld_3[which(mgraph$smallworld_3 == Inf)] <- max(mgraph$smallworld_3[!mgraph$smallworld_3==Inf])
summary(lmer(smallworld_2 ~ IIP_communion.c*IIP_agency.c + (1 | PTNUM), mgraph))
summary(lmer(smallworld_2 ~ PAIBORtot.c + (1 | PTNUM), mgraph)) #spurious?

summary(lmer(smallworld_3 ~ IIP_communion.c*IIP_agency.c + (1 | PTNUM), mgraph))
summary(lmer(smallworld_3 ~ PAIBORtot.c + (1 | PTNUM), mgraph))
