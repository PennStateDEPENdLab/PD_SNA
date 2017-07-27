library(lme4)
library(igraph)
library(afex)
library(ggplot2)
library(reshape2)
library(lattice)
setwd("/Users/michael/Tresors/PD_SNA/Couples_Baseline_SNA")
source(file.path(getMainDir(), "Miscellaneous", "Global_Functions.R"))
load(file="data/SNA_Processed_6Oct2015.RData")
load(file="data/couple_graph_measures_6Oct2015.RData") #has dyad members in graph
patientBw <- with(couple_vertex_metrics, betweenness_weighted_2[role=="patient"])
partnerBw <- with(couple_vertex_metrics, betweenness_weighted_2[role=="partner"])

#load(file="data/couple_graph_measures_6Oct2015.RData") #has dyad members in graph
load(file="data/selfreports/couples_baseline_clinical_9Oct2015.RData")
load("data/couple_graph_measures_nodyad_15Oct2015.RData")

#get vertex betweenness for each partner and get relationship duration
coupleLength <- sapply(couple_graphs, function(p) {
  mean(p$actors$patient.YrsKnown[p$actors$role=="partner"], p$actors$partner.YrsKnown[p$actors$role=="patient"]) 
})
ptnums <- sapply(couple_graphs, "[[", "PTNUM")
df <- data.frame(PTNUM=ptnums, reldur=coupleLength, pbw_1=patientBw, pbw_0=partnerBw)
couples_baseline_clin <- merge(couples_baseline_clin, df, by="PTNUM")
couples_clin_wide <- merge(couples_clin_wide, df, by="PTNUM")
mean(df$reldur)
sd(df$reldur)

#remove rule-out IDs (per email from Emily Landis 7Oct2015)
ro_ids <- c(8002, 8012, 8021, 8025, 8037, 8076, 8077, 8079, 8089, 8090, 8040, 8047, 8048, 8097, 8119, 8125)
couples_baseline_clin <- subset(couples_baseline_clin, UsrID < 90000 & !PTNUM %in% ro_ids)
#sink("private/clinical_missing_report_7Oct2015.txt")
#missingDataReport(couples_baseline_clin, idVars="UsrID")
#sink()

mgraph <- merge(couple_graph_metrics, couples_clin_wide, by=c("PTNUM"), all.x=TRUE)
mvertex <- merge(couple_vertex_metrics, couples_clin_wide, by=c("PTNUM"), all.x=TRUE)


mgraph$PTNUM <- factor(mgraph$PTNUM)
mvertex$PTNUM <- factor(mvertex$PTNUM)
predstocenter <- c("PAIBORtot", "iip_agency", "iip_communion", "iip_elevation", "IIP_PD1", "IIP_PD2", "IIP_PD3", "ECRanx", "ECRavoid", "DASTotal", 
    "Victim", "Perp", "bordl_sidp", "narci_sidp", "antso_sidp", "avoid_sidp", "OPD_sidp", "nonarc", "pdtot", "nobpd", "reldur")
mgraph <- f_centerPredictors(mgraph, apply(expand.grid(predstocenter, c("0", "1")), 1, paste, collapse="_"), addsuffix=".c")
mvertex <- f_centerPredictors(mvertex, apply(expand.grid(predstocenter, c("0", "1")), 1, paste, collapse="_"), addsuffix=".c")

##Betweenness of couple
histogram(mgraph$pbw_0)
histogram(mgraph$pbw_1)
t.test(mgraph$pbw_0, mgraph$pbw_1) #NS

tmp <- plyr::rename(mgraph, c(iip_agency_0.c="domc_0", iip_agency_1.c="domc_1", iip_communion_0.c="affc_0", iip_communion_1.c="affc_1",
        iip_elevation_0.c="iipe_0", iip_elevation_1.c="iipe_1", bordl_sidp_0.c="bpdc_0", bordl_sidp_1.c="bpdc_1"))

#summary(glm(pbw ~ p_age_0 + iip_communion_0.c + iip_communion_1.c, x, family=poisson))
output <- mplusAPIM(tmp, "pbw", c("domc", "affc", "iipe"))
output <- mplusAPIM(tmp, "pbw", c("domc", "affc"))
output <- mplusAPIM(tmp, "pbw", c("domc"))
output <- mplusAPIM(tmp, "pbw", c("affc"), exogCorr=TRUE)

output <- mplusAPIM(tmp, "pbw", c("pdtot"), exogCorr=TRUE)
output <- mplusAPIM(tmp, "pbw", c("bpdc", "nobpd"), exogCorr=TRUE) #sig neg partner Bw for non-bpd Sx

###NEW: DYADIC METRICS WITHOUT COUPLE INCLUDED

#just number of shared alters
x <- ddply(mvertex, .(PTNUM), function(subdf) {
      c(nshared=sum(subdf$role=="shared"))
    })

histogram(x$nshared)
x <- merge(x, couples_clin_wide, by="PTNUM")
x <- f_centerPredictors(x, apply(expand.grid(predstocenter, c("0", "1")), 1, paste, collapse="_"), addsuffix=".c")

summary(glm(nshared ~ p_age_0 + iip_communion_0.c + iip_communion_1.c, x, family=poisson))
summary(glm(nshared ~ p_age_0 + p_age_1 + iip_agency_0.c + iip_agency_1.c, x, family=poisson))
summary(m <- glm(nshared ~ iip_agency_0 + iip_agency_1, x, family=poisson)) #p_age_0 + 
newdf <- with(mgraph, expand.grid(iip_agency_0=round(seq(min(iip_agency_0, na.rm=T), max(iip_agency_0, na.rm=T), length.out=10), 2), 
        iip_agency_1=round(seq(min(iip_agency_1, na.rm=T), max(iip_agency_1, na.rm=T), length.out=10), 2)))

ggplot(data=transform(newdf, yp=exp(predict(m, newdf))), 
    aes(y=yp, x=iip_agency_1, color=factor(iip_agency_0))) + stat_smooth(method=lm)

summary(m <- glm(nshared ~ p_age_0 + p_age_1 + iip_elevation_0 + iip_elevation_1, x, family=poisson)) #p_age_0 
summary(m <- glm(nshared ~ iip_elevation_0 * iip_elevation_1, x, family=poisson)) #p_age_0
newdf <- with(mgraph, expand.grid(iip_elevation_0=round(seq(min(iip_elevation_0, na.rm=T), max(iip_elevation_0, na.rm=T), length.out=10), 2), 
        iip_elevation_1=round(seq(min(iip_elevation_1, na.rm=T), max(iip_elevation_1, na.rm=T), length.out=10), 2)))

ggplot(data=transform(newdf, yp=exp(predict(m, newdf))), 
    aes(y=yp, x=iip_elevation_1, color=factor(iip_elevation_0))) + stat_smooth(method=lm)

summary(glm(nshared ~ p_age_0 + pdtot_0.c + pdtot_1.c, x, family=poisson))
summary(glm(nshared ~ p_age_0 + narci_sidp_0.c + narci_sidp_1.c, x, family=poisson))
summary(glm(nshared ~ p_age_0 + bordl_sidp_0.c + bordl_sidp_1.c, x, family=poisson))
summary(glm(nshared ~ p_age_0 + bordl_sidp_0.c + bordl_sidp_1.c + nobpd_0.c + nobpd_1.c, x, family=poisson))
summary(glm(nshared ~ p_age_0 + p_age_1 + narci_sidp_0.c + narci_sidp_1.c + nonarc_0.c + nonarc_1.c, x, family=poisson))
summary(glm(nshared ~ p_age_0 + p_age_1 + pdtot_0.c + pdtot_1.c, x, family=poisson))
summary(lm(nshared ~ p_age_0 + p_age_1 + pdtot_0.c + pdtot_1.c, x))
summary(lm(nshared ~ p_age_0 + p_age_1 + narci_sidp_0.c + narci_sidp_1.c + nonarc_0.c + nonarc_1.c, x))
cor.test(x$nshared, x$narci_sidp_0)
cor.test(x$p_age_0, x$p_age_1)
cor(x$nshared, x$narci_sidp_1, use="pairwise.complete.obs")



#graph-wide analysis
cor.test(~ PAIBORtot_0.c + PAIBORtot_1.c, mgraph) #weak correlation
cor.test(~ iip_agency_0.c + iip_agency_1.c, mgraph) #weak correlation
cor.test(~ iip_communion_0.c + iip_communion_1.c, mgraph) #weak correlation
cor.test(~ ECRanx_0.c + ECRanx_1.c, mgraph) #weak correlation
cor.test(~ ECRavoid_0.c + ECRavoid_1.c, mgraph) #weak correlation
cor(mgraph[,c("IIP_PD1_0.c", "IIP_PD1_1.c", "IIP_PD2_0.c", "IIP_PD2_1.c", "IIP_PD3_0.c", "IIP_PD3_1.c")], use="pairwise.complete.obs")
cor(mgraph[,c("antso_sidp_0.c", "antso_sidp_1.c", "bordl_sidp_0.c", "bordl_sidp_1.c", "narci_sidp_0.c", "narci_sidp_1.c")], use="pairwise.complete.obs")

#partner BPD Sx contribute to weaker centralization

#most interested in number of modules, module size, modularity, and transitivity
#2 threshold makes sense because this represents at least some plausible connection among alters. 3 is strict
histogram(mgraph$avgsize_optcommunity_weighted_2)
histogram(sqrt(mgraph$avgsize_optcommunity_weighted_2))
moments::skewness(mgraph$avgsize_optcommunity_weighted_2)
moments::skewness(sqrt(mgraph$avgsize_optcommunity_weighted_2))
#num communities looks roughly poisson at thresh 2
summary(m <- lm(sqrt(avgsize_optcommunity_weighted_2) ~ iip_agency_0.c * iip_agency_1.c + iip_communion_0.c * iip_communion_1.c, mgraph)) #null
summary(m <- lm(sqrt(avgsize_optcommunity_weighted_2) ~ iip_elevation_0.c * iip_elevation_1.c, mgraph))
summary(m <- lm(sqrt(avgsize_optcommunity_weighted_2) ~ pdtot_0.c * pdtot_1.c, mgraph)) #null
summary(m <- lm(sqrt(avgsize_optcommunity_weighted_2) ~ bordl_sidp_0.c * bordl_sidp_1.c, mgraph))
summary(m <- lm(sqrt(avgsize_optcommunity_weighted_2) ~ bordl_sidp_0.c * bordl_sidp_1.c + nobpd_0.c * nobpd_1.c, mgraph))

summary(m <- lm(sqrt(avgsize_optcommunity_weighted_2) ~ bordl_sidp_0 * bordl_sidp_1, mgraph))
newdf <- with(mgraph, expand.grid(bordl_sidp_0=round(seq(min(bordl_sidp_0, na.rm=T), max(bordl_sidp_0, na.rm=T), length.out=10), 2), 
        bordl_sidp_1=round(seq(min(bordl_sidp_1, na.rm=T), max(bordl_sidp_1, na.rm=T), length.out=10), 2)))

library(ggplot2)
#weirdness in the asymmetry
p <- ggplot(data=transform(newdf, yp=predict(m, newdf)^2), 
    aes(y=yp, x=bordl_sidp_1, color=factor(bordl_sidp_0))) + stat_smooth(method=lm)
p


summary(m <- lm(number_optcommunity_weighted_2 ~ iip_agency_0.c * iip_agency_1.c + iip_communion_0.c * iip_communion_1.c, mgraph)) #null
summary(m <- glm(number_optcommunity_weighted_2 ~ iip_agency_0.c * iip_agency_1.c + iip_communion_0.c * iip_communion_1.c, mgraph, family=poisson)) #null
summary(m <- glm(number_optcommunity_weighted_3 ~ iip_agency_0.c * iip_agency_1.c + iip_communion_0.c * iip_communion_1.c, mgraph, family=poisson)) #null
summary(m <- glm(number_optcommunity_weighted_2 ~ iip_elevation_0.c * iip_elevation_1.c, mgraph, family=poisson)) #trend toward interaction
summary(m <- glm(number_optcommunity_weighted_3 ~ iip_elevation_0.c * iip_elevation_1.c, mgraph, family=poisson)) #interaction

summary(m <- glm(number_optcommunity_weighted_3 ~ iip_elevation_0 * iip_elevation_1, mgraph, family=poisson)) #interaction
newdf <- with(mgraph, expand.grid(iip_elevation_0=round(seq(min(iip_elevation_0, na.rm=T), max(iip_elevation_0, na.rm=T), length.out=10), 2), 
        iip_elevation_1=round(seq(min(iip_elevation_1, na.rm=T), max(iip_elevation_1, na.rm=T), length.out=10), 2)))

library(ggplot2)
ggplot(data=transform(newdf, yp=exp(predict(m, newdf))), 
    aes(y=yp, x=iip_elevation_1, color=factor(iip_elevation_0))) + stat_smooth(method=lm) + theme_bw()

summary(m <- lm(number_optcommunity_weighted_2 ~ iip_elevation_0.c * iip_elevation_1.c, mgraph)) #null
summary(m <- lm(number_optcommunity_weighted_3 ~ iip_agency_0.c * iip_agency_1.c +  iip_communion_0.c * iip_communion_1.c + iip_elevation_0.c * iip_elevation_1.c, mgraph)) #null
summary(m <- glm(number_optcommunity_weighted_2 ~ iip_agency_0.c * iip_agency_1.c +  iip_communion_0.c * iip_communion_1.c + iip_elevation_0.c * iip_elevation_1.c, mgraph, family=poisson)) #null
summary(m <- glm(number_optcommunity_weighted_3 ~ iip_agency_0.c * iip_agency_1.c +  iip_communion_0.c * iip_communion_1.c + iip_elevation_0.c * iip_elevation_1.c, mgraph, family=poisson)) #null

histogram(mgraph$number_optcommunity_weighted_2)
histogram(mgraph$number_optcommunity_weighted_3)

summary(m <- glm(number_optcommunity_weighted_2 ~ bordl_sidp_0.c * bordl_sidp_1.c + nobpd_0.c * nobpd_1.c, mgraph, family=poisson)) 
summary(m <- glm(number_optcommunity_weighted_2 ~ bordl_sidp_0.c * bordl_sidp_1.c, mgraph, family=poisson))
#with(mgraph, interaction.plot(bordl_sidp_1, bordl_sidp_0, response=number_optcommunity_weighted_2))

summary(m <- glm(number_optcommunity_weighted_2 ~ bordl_sidp_0 * bordl_sidp_1, mgraph, family=poisson))
newdf <- with(mgraph, expand.grid(bordl_sidp_0=round(seq(min(bordl_sidp_0, na.rm=T), max(bordl_sidp_0, na.rm=T), length.out=10), 2), 
        bordl_sidp_1=round(seq(min(bordl_sidp_1, na.rm=T), max(bordl_sidp_1, na.rm=T), length.out=10), 2)))

library(ggplot2)
#weirdness in the asymmetry
p <- ggplot(data=transform(newdf, yp=exp(predict(m, newdf))), 
    aes(y=yp, x=bordl_sidp_1, color=factor(bordl_sidp_0))) + stat_smooth(method=lm)

p

summary(m <- lm(number_optcommunity_weighted_2 ~ bordl_sidp_0 + bordl_sidp_1, mgraph))
summary(m <- lm(largest_clique_size_2 ~ bordl_sidp_0.c * bordl_sidp_1.c + nobpd_0.c * nobpd_1.c , mgraph))
summary(m <- lm(largest_clique_size_2 ~ pdtot_0.c * pdtot_1.c, mgraph))


summary(m <- lm(transitivity_2 ~ bordl_sidp_0.c + bordl_sidp_1.c + nobpd_0.c + nobpd_1.c, mgraph))
summary(m <- lm(transitivity_2 ~ pdtot_0.c + pdtot_1.c, mgraph))
summary(m <- lm(transitivity_2 ~ narci_sidp_0.c + narci_sidp_1.c + nonarc_0.c + nonarc_1.c, mgraph))
#summary(m <- lm(transitivity_2 ~ avoid_sidp_0.c + avoid_sidp_1.c + pdtot_0.c + pdtot_1.c, mgraph))

#non-narc Sx in partner associated with weaker transitivity



#m <- lm(centralization_degree_2 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph)
#summary(m)

#partner Sx contributes to weaker centralization, patient to greater.
#m <- lm(centralization_degree_3 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph)
#summary(m)

#summary(m <- lm(centralization_degree_2 ~ ECRanx_0.c*ECRanx_1.c, mgraph))
#summary(m <- lm(centralization_degree_3 ~ ECRanx_0.c*ECRanx_1.c, mgraph))

summary(m <- lm(centralization_degree_2 ~ iip_agency_0.c + iip_agency_1.c + iip_communion_0.c + iip_communion_1.c, mgraph))
summary(m <- lm(centralization_degree_3 ~ iip_agency_0.c + iip_agency_1.c + iip_communion_0.c + iip_communion_1.c, mgraph))
summary(m <- lm(centralization_degree_2 ~ nobpd_0.c + nobpd_1.c + bordl_sidp_0.c + bordl_sidp_1.c, mgraph))
summary(m <- lm(centralization_degree_2 ~ pdtot_0 + pdtot_1, mgraph))

#use lavaan to check patient/partner constraints
m <- 'centralization_degree_2 ~ c1*nobpd_0.c + c1*nobpd_1.c + c2*bordl_sidp_0.c + c2*bordl_sidp_1.c'
m <- sem(m, mgraph, estimator="MLR", missing="ML")

summary(m <- lm(centralization_degree_3 ~ iip_agency_0.c + iip_agency_1.c + iip_communion_0.c + iip_communion_1.c, mgraph))




summary(m <- lm(centralization_eigenvector_2 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))
summary(m <- lm(centralization_eigenvector_3 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))
summary(m <- lm(centralization_eigenvector_2 ~ iip_agency_0.c*iip_agency_1.c + iip_communion_0.c*iip_communion_1.c, mgraph))
summary(m <- lm(centralization_eigenvector_3 ~ iip_agency_0.c*iip_agency_1.c + iip_communion_0.c*iip_communion_1.c, mgraph))
summary(m <- lm(centralization_eigenvector_3 ~ ECRanx_0.c*ECRanx_1.c + ECRavoid_0.c*ECRavoid_1.c, mgraph))
summary(m <- lm(centralization_eigenvector_3 ~ ECRanx_0.c*ECRanx_1.c, mgraph)) #significant interaction such that both partners' attachment anxiety contribute to greater centralization

#assortativity of nomination (patient, partner, or both)
#something messed up with this at the moment: assortativity all zero
summary(m <- lm(nomination_assortativity_2 ~ ECRanx_0.c*ECRanx_1.c, mgraph))
summary(m <- lm(nomination_assortativity_3 ~ ECRanx_0.c*ECRanx_1.c, mgraph))


summary(m <- lm(largest_clique_size_2 ~ ECRanx_0.c*ECRanx_1.c, mgraph))
summary(m <- lm(largest_clique_size_3 ~ ECRanx_0.c*ECRanx_1.c, mgraph))

summary(m <- lm(largest_clique_size_2 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))
summary(m <- lm(largest_clique_size_3 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))

summary(m <- lm(largest_clique_size_2 ~ PAIBORtot_0.c*PAIBORtot_1.c + ECRanx_0.c*ECRanx_1.c + ECRavoid_0.c*ECRavoid_1.c, mgraph))
summary(m <- lm(largest_clique_size_3 ~ PAIBORtot_0.c*PAIBORtot_1.c + ECRanx_0.c*ECRanx_1.c + ECRavoid_0.c*ECRavoid_1.c, mgraph))

summary(m <- lm(largest_clique_size_2 ~ PAIBORtot_0.c*PAIBORtot_1.c + iip_communion_0.c*iip_communion_1.c, mgraph))
summary(m <- lm(largest_clique_size_3 ~ PAIBORtot_0.c*PAIBORtot_1.c + iip_communion_0.c*iip_communion_1.c, mgraph))
summary(m <- lm(largest_clique_size_3 ~ PAIBORtot_0.c*PAIBORtot_1.c + iip_communion_0.c*iip_communion_1.c + 
            PAIBORtot_0.c*iip_communion_0.c + PAIBORtot_0.c*iip_communion_1.c +
            PAIBORtot_0.c*iip_agency_0.c + PAIBORtot_0.c*iip_agency_1.c, mgraph))

summary(m2 <- lm(largest_clique_size_3 ~ PAIBORtot_0.c*PAIBORtot_1.c*iip_communion_0.c*iip_communion_1.c, mgraph))


summary(m <- lm(num_maximal_3cliques_2 ~ PAIBORtot_0.c*PAIBORtot_1.c*iip_communion_0.c*iip_communion_1.c, mgraph))
summary(m <- lm(num_maximal_3cliques_3 ~ PAIBORtot_0.c*PAIBORtot_1.c*iip_communion_0.c*iip_communion_1.c, mgraph))

summary(m <- lm(num_maximal_3cliques_2 ~ PAIBORtot_0.c*PAIBORtot_1.c*iip_agency_0.c*iip_agency_1.c, mgraph))
summary(m <- lm(num_maximal_3cliques_3 ~ PAIBORtot_0.c*PAIBORtot_1.c*iip_agency_0.c*iip_agency_1.c, mgraph))

summary(m <- lm(num_maximal_3cliques_2 ~ iip_agency_0.c*iip_agency_1.c, mgraph))
summary(m <- lm(num_maximal_3cliques_3 ~ iip_agency_0.c*iip_agency_1.c, mgraph))

summary(m <- lm(num_maximal_3cliques_2 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))
summary(m <- lm(num_maximal_3cliques_3 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))

summary(m <- lm(avgsize_maximal_3cliques_2 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))
summary(m <- lm(avgsize_maximal_3cliques_3 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))

summary(m <- lm(avgsize_maximal_3cliques_2 ~ iip_communion_0.c*iip_communion_1.c, mgraph))
summary(m <- lm(avgsize_maximal_3cliques_3 ~ iip_communion_0.c*iip_communion_1.c, mgraph))
summary(m <- lm(avgsize_maximal_3cliques_2 ~ iip_agency_0.c*iip_agency_1.c, mgraph))
summary(m <- lm(avgsize_maximal_3cliques_3 ~ iip_agency_0.c*iip_agency_1.c, mgraph))

#partner agency tends to go with larger couple number max cliques and larger average size of cliques

summary(m <- lm(number_optcommunity_weighted_2 ~ bordl_sidp_0.c*bordl_sidp_1.c, mgraph))
summary(m <- lm(number_optcommunity_weighted_3 ~ bordl_sidp_0.c*bordl_sidp_1.c, mgraph))
summary(m <- lm(number_optcommunity_weighted_2 ~ narci_sidp_0.c*narci_sidp_1.c, mgraph))
summary(m <- lm(number_optcommunity_weighted_3 ~ narci_sidp_0.c*narci_sidp_1.c, mgraph))
summary(m <- lm(number_optcommunity_weighted_2 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))
summary(m <- lm(number_optcommunity_weighted_3 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))
summary(m <- lm(number_optcommunity_weighted_2 ~ iip_communion_0.c*iip_communion_1.c, mgraph))
summary(m <- lm(number_optcommunity_weighted_3 ~ iip_communion_0.c*iip_communion_1.c, mgraph))
summary(m <- lm(number_optcommunity_weighted_2 ~ iip_agency_0.c*iip_agency_1.c, mgraph))
summary(m <- lm(number_optcommunity_weighted_3 ~ iip_agency_0.c*iip_agency_1.c, mgraph))
summary(m <- lm(number_optcommunity_weighted_2 ~ DASTotal_0.c*DASTotal_1.c, mgraph))
summary(m <- lm(number_optcommunity_weighted_3 ~ DASTotal_0.c*DASTotal_1.c, mgraph))

summary(m <- lm(avgsize_optcommunity_weighted_2 ~ bordl_sidp_0.c*bordl_sidp_1.c, mgraph))
summary(m <- lm(avgsize_optcommunity_weighted_3 ~ bordl_sidp_0.c*bordl_sidp_1.c, mgraph))
summary(m <- lm(avgsize_optcommunity_weighted_2 ~ narci_sidp_0.c*narci_sidp_1.c, mgraph)) #**
summary(m <- lm(avgsize_optcommunity_weighted_3 ~ narci_sidp_0.c*narci_sidp_1.c, mgraph)) #*
summary(m <- lm(avgsize_optcommunity_weighted_2 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))
summary(m <- lm(avgsize_optcommunity_weighted_3 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))
summary(m <- lm(avgsize_optcommunity_weighted_2 ~ iip_communion_0.c*iip_communion_1.c, mgraph))
summary(m <- lm(avgsize_optcommunity_weighted_3 ~ iip_communion_0.c*iip_communion_1.c, mgraph))
summary(m <- lm(avgsize_optcommunity_weighted_2 ~ iip_agency_0.c*iip_agency_1.c, mgraph))
summary(m <- lm(avgsize_optcommunity_weighted_3 ~ iip_agency_0.c*iip_agency_1.c, mgraph))
summary(m <- lm(avgsize_optcommunity_weighted_2 ~ DASTotal_0.c*DASTotal_1.c, mgraph))
summary(m <- lm(avgsize_optcommunity_weighted_3 ~ DASTotal_0.c*DASTotal_1.c, mgraph))

summary(m <- lm(modularity_optcommunity_weighted_2 ~ bordl_sidp_0.c*bordl_sidp_1.c, mgraph))
summary(m <- lm(modularity_optcommunity_weighted_3 ~ bordl_sidp_0.c*bordl_sidp_1.c, mgraph))
summary(m <- lm(modularity_optcommunity_weighted_2 ~ narci_sidp_0.c*narci_sidp_1.c, mgraph))
summary(m <- lm(modularity_optcommunity_weighted_3 ~ narci_sidp_0.c*narci_sidp_1.c, mgraph))
summary(m <- lm(modularity_optcommunity_weighted_2 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))
summary(m <- lm(modularity_optcommunity_weighted_3 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))
summary(m <- lm(modularity_optcommunity_weighted_2 ~ iip_communion_0.c*iip_communion_1.c, mgraph))
summary(m <- lm(modularity_optcommunity_weighted_3 ~ iip_communion_0.c*iip_communion_1.c, mgraph))
summary(m <- lm(modularity_optcommunity_weighted_2 ~ iip_agency_0.c*iip_agency_1.c, mgraph))
summary(m <- lm(modularity_optcommunity_weighted_3 ~ iip_agency_0.c*iip_agency_1.c, mgraph))

#density
#summary(m <- lm(density_2 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))
#summary(m <- lm(density_3 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))
summary(m <- lm(density_2 ~ iip_communion_0.c + iip_communion_1.c, mgraph))
summary(m <- lm(density_3 ~ iip_communion_0.c + iip_communion_1.c, mgraph))
summary(m <- lm(density_2 ~ iip_agency_0.c + iip_agency_1.c, mgraph))
summary(m <- lm(density_3 ~ iip_agency_0.c + iip_agency_1.c, mgraph))
summary(m <- lm(density_2 ~ iip_elevation_0.c + iip_elevation_1.c, mgraph))
summary(m <- lm(density_3 ~ iip_elevation_0.c + iip_elevation_1.c + reldur, mgraph))

newdf <- expand.grid(iip_elevation_0.c=seq(min(mgraph$iip_elevation_0.c, na.rm=T), max(mgraph$iip_elevation_0.c, na.rm=T)),
    iip_elevation_1.c=seq(min(mgraph$iip_elevation_1.c, na.rm=T), max(mgraph$iip_elevation_1.c, na.rm=T)))

library(ggplot2)
p <- ggplot(data=transform(newdf, yp=predict(m, newdf)), 
    aes(y=yp, x=iip_elevation_0.c, color=factor(iip_elevation_1.c))) + stat_smooth(method=lm)


summary(m <- lm(density_2 ~ iip_agency_0.c + iip_agency_1.c + iip_elevation_0.c + iip_elevation_1.c, mgraph))
summary(m <- lm(density_3 ~ iip_agency_0.c + iip_agency_1.c + iip_elevation_0.c + iip_elevation_1.c, mgraph))

summary(m <- lm(density_2 ~ pdtot_0.c + pdtot_1.c, mgraph))
summary(m <- lm(density_3 ~ pdtot_0.c + pdtot_1.c, mgraph))

summary(m <- lm(density_3 ~ iip_agency_0.c + iip_agency_1.c + iip_elevation_0.c + iip_elevation_1.c, mgraph))

#summary(m <- lm(density_2 ~ ECRanx_0.c*ECRanx_1.c, mgraph))
#summary(m <- lm(density_3 ~ ECRavoid_0.c*ECRavoid_1.c, mgraph))

#transitivity
summary(m <- lm(transitivity_2 ~ bordl_sidp_0.c + bordl_sidp_1.c + nobpd_0.c + nobpd_1.c, mgraph))
summary(m <- lm(transitivity_2 ~ pdtot_0.c * pdtot_1.c, mgraph))
#summary(m <- lm(transitivity_3 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))


summary(m <- lm(transitivity_2 ~ iip_communion_0.c + iip_communion_1.c, mgraph))
summary(m <- lm(transitivity_3 ~ iip_communion_0.c + iip_communion_1.c, mgraph))
summary(m <- lm(transitivity_2 ~ iip_agency_0.c * iip_agency_1.c, mgraph))
summary(m <- lm(transitivity_3 ~ iip_agency_0.c + iip_agency_1.c, mgraph))
#summary(m <- lm(transitivity_2 ~ ECRanx_0.c*ECRanx_1.c, mgraph))
#summary(m <- lm(transitivity_3 ~ ECRavoid_0.c*ECRavoid_1.c, mgraph))

#edge betweenness between partners (sqrt transform to normalize)
summary(m <- lm(sqrt(edge_betweenness_binary_2) ~ bordl_sidp_0.c*bordl_sidp_1.c, mgraph))
summary(m <- lm(sqrt(edge_betweenness_weighted_2) ~ bordl_sidp_0.c*bordl_sidp_1.c, mgraph))

newdf <- expand.grid(bordl_sidp_0.c=seq(min(mgraph$bordl_sidp_0.c, na.rm=T), max(mgraph$bordl_sidp_0.c, na.rm=T)),
    bordl_sidp_1.c=seq(min(mgraph$bordl_sidp_1.c, na.rm=T), max(mgraph$bordl_sidp_1.c, na.rm=T)))

library(ggplot2)
p <- ggplot(data=transform(newdf, yp=predict(m, newdf)), 
    aes(y=yp, x=bordl_sidp_0.c, color=factor(bordl_sidp_1.c))) + stat_smooth(method=lm)

interaction.plot(mgraph$bordl_sidp_0.c, mgraph$bordl_sidp_1.c, sqrt(mgraph$edge_betweenness_binary_2))





summary(m <- lm(sqrt(edge_betweenness_binary_2) ~ narci_sidp_0.c*narci_sidp_1.c, mgraph))
summary(m <- lm(sqrt(edge_betweenness_weighted_2) ~ narci_sidp_0.c*narci_sidp_1.c, mgraph))
summary(m <- lm(sqrt(edge_betweenness_binary_2) ~ antso_sidp_0.c*antso_sidp_1.c, mgraph))
summary(m <- lm(sqrt(edge_betweenness_weighted_2) ~ antso_sidp_0.c*antso_sidp_1.c, mgraph))

summary(m <- lm(sqrt(edge_betweenness_binary_2) ~ avoid_sidp_0.c*avoid_sidp_1.c, mgraph))
summary(m <- lm(sqrt(edge_betweenness_weighted_2) ~ avoid_sidp_0.c*avoid_sidp_1.c, mgraph))

summary(m <- lm(sqrt(edge_betweenness_binary_2) ~ OPD_sidp_0.c*OPD_sidp_1.c, mgraph))
summary(m <- lm(sqrt(edge_betweenness_weighted_2) ~ OPD_sidp_0.c*OPD_sidp_1.c, mgraph))


summary(m <- lm(sqrt(edge_betweenness_binary_2) ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))
summary(m <- lm(sqrt(edge_betweenness_weighted_2) ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))
summary(m <- lm(sqrt(edge_betweenness_binary_2) ~ iip_communion_0.c*iip_communion_1.c, mgraph))
summary(m <- lm(sqrt(edge_betweenness_weighted_2) ~ iip_communion_0.c*iip_communion_1.c, mgraph))
summary(m <- lm(sqrt(edge_betweenness_binary_2) ~ iip_agency_0.c*iip_agency_1.c, mgraph))
summary(m <- lm(sqrt(edge_betweenness_weighted_2) ~ iip_agency_0.c*iip_agency_1.c, mgraph))
summary(m <- lm(sqrt(edge_betweenness_binary_2) ~ ECRanx_0.c*ECRanx_1.c, mgraph))
summary(m <- lm(sqrt(edge_betweenness_weighted_2) ~ ECRavoid_0.c*ECRavoid_1.c, mgraph))

summary(m <- lm(sqrt(edge_betweenness_binary_2) ~ ECRanx_0.c*ECRanx_1.c + PAIBORtot_0.c*PAIBORtot_1.c, mgraph))



#vertex measures
summary(m <- lmer(degree_2 ~ PAIBORtot_0.c*PAIBORtot_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(degree_3 ~ PAIBORtot_0.c*PAIBORtot_1.c + (1 | PTNUM), mvertex))

summary(m <- lmer(degree_2 ~ PAIBORtot_0.c*PAIBORtot_1.c + iip_communion_0.c*iip_communion_1.c + iip_agency_0.c*iip_agency_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(degree_3 ~ PAIBORtot_0.c*PAIBORtot_1.c + iip_communion_0.c*iip_communion_1.c + iip_agency_0.c*iip_agency_1.c + (1 | PTNUM), mvertex))

#partner, but not patient, BPD Sx associated with lower degree. Partner agency goes with higher degree, and both partners having more BPD Sx goes with lower degree
mixed(degree_2 ~ PAIBORtot_0.c*PAIBORtot_1.c + iip_communion_0.c*iip_communion_1.c + iip_agency_0.c*iip_agency_1.c + (1 | PTNUM), mvertex)
mixed(degree_3 ~ PAIBORtot_0.c*PAIBORtot_1.c + iip_communion_0.c*iip_communion_1.c + iip_agency_0.c*iip_agency_1.c + (1 | PTNUM), mvertex)

#clinician-rated BPD, but not OPD, Sx in partner associated with weaker degree in dyadic network
summary(m <- lmer(degree_2 ~ bordl_sidp_0.c*bordl_sidp_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(degree_3 ~ bordl_sidp_0.c*bordl_sidp_1.c + (1 | PTNUM), mvertex))

summary(m <- lmer(degree_2 ~ iip_communion_0.c*iip_communion_1.c + iip_agency_0.c*iip_agency_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(degree_3 ~ iip_communion_0.c*iip_communion_1.c + iip_agency_0.c*iip_agency_1.c + (1 | PTNUM), mvertex))

#greater narcissism in both partners associated with weaker degree
summary(m <- lmer(degree_2 ~ narci_sidp_0.c + narci_sidp_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(degree_2 ~ narci_sidp_0.c*role + narci_sidp_1.c*role + (1 | PTNUM), mvertex))
car::Anova(m)
cm <- lmerCellMeans(m, divide=c("narci_sidp_0.c", "narci_sidp_1.c"), divide.prefix=FALSE)
ggplot(cm, aes(x=narci_sidp_0.c, y=degree_2, color=narci_sidp_1.c, ymin=degree_2-se, ymax=degree_2+se)) + 
    geom_pointrange(size=2) + theme_bw() +
    facet_wrap(~role, scales="free")


summary(m <- lmer(degree_2 ~ pdtot_0.c*role + pdtot_1.c*role + (1 | PTNUM), mvertex, REML=FALSE))
cm <- lmerCellMeans(m, divide=c("pdtot_0.c", "pdtot_1.c"), divide.prefix=FALSE)
ggplot(cm, aes(x=pdtot_0.c, y=degree_2, color=pdtot_1.c, ymin=degree_2-se, ymax=degree_2+se)) + 
    geom_pointrange(size=2) + theme_bw() +
    facet_wrap(~role)

#random slope for role definitely fits better, but washes out some PD effects...
summary(m <- lmer(degree_2 ~ pdtot_0.c*role + pdtot_1.c*role + (1 + role | PTNUM), mvertex, REML=FALSE))
car::Anova(m)

summary(m <- lmer(degree_2 ~ bordl_sidp_0.c*role + bordl_sidp_1.c*role + nobpd_0.c*role + nobpd_1.c*role + (1 + role | PTNUM), mvertex, REML=FALSE))
car::Anova(m)

summary(m <- lmer(degree_2 ~ narci_sidp_0.c*role + narci_sidp_1.c*role + nonarc_0.c*role + nonarc_1.c*role + (1 + role | PTNUM), mvertex, REML=FALSE))
summary(m <- lmer(degree_3 ~ narci_sidp_0.c*role + narci_sidp_1.c*role + nonarc_0.c*role + nonarc_1.c*role + (1 + role | PTNUM), mvertex, REML=FALSE))
car::Anova(m)
summary(m <- lmer(degree_2 ~ iip_communion_0.c + bordl_sidp_1.c + nobpd_0.c*role + nobpd_1.c*role + (1 + role | PTNUM), mvertex, REML=FALSE))

#holds
summary(m <- lmer(strength_3 ~ narci_sidp_0.c*role + narci_sidp_1.c*role + nonarc_0.c*role + nonarc_1.c*role + reldur + (1 + role | PTNUM), mvertex, REML=FALSE))
car::Anova(m)


summary(m <- lmer(strength_3 ~ iip_communion_0.c*role + iip_communion_1.c*role + iip_agency_0.c*role + iip_agency_1.c*role + reldur + (1 + role | PTNUM), mvertex, REML=FALSE))
car::Anova(m)



#huh, as partner has higher iip communion, his own alters tend to have higher degree compard to shared or patient only
#as patient tends to be more dominant, his own and shared alters tend to have higher degree
summary(m <- lmer(degree_2 ~ iip_agency_0.c*role + iip_agency_1.c*role + iip_communion_0.c*role + iip_communion_1.c*role + (1 + role | PTNUM), mvertex, REML=FALSE))
car::Anova(m)
summary(m <- lmer(strength_3 ~ iip_agency_0.c*role + iip_agency_1.c*role + iip_communion_0.c*role + iip_communion_1.c*role + (1 + role | PTNUM), mvertex, REML=FALSE))

#horrible skew
summary(m <- lmer(log(betweenness_weighted_2 + 1) ~ narci_sidp_0.c*role + narci_sidp_1.c*role + nonarc_0.c*role + nonarc_1.c*role + (1 + role | PTNUM), mvertex, REML=FALSE))


#way to look at how central  patient and partner alters are in the network
#so here, partner's alters are more peripheral as patient's have higher PD Sx (because eccentricity higher in partners)
summary(m <- lmer(eccentricity_2 ~ pdtot_0.c*role + pdtot_1.c*role + (1 + role | PTNUM), mvertex, REML=FALSE))
car::Anova(m)

summary(m <- lmer(eccentricity_3 ~ pdtot_0.c*role + pdtot_1.c*role + (1 + role | PTNUM), mvertex, REML=FALSE))
car::Anova(m)


summary(m <- lmer(eccentricity_2 ~ narci_sidp_0.c*role + narci_sidp_1.c*role + nonarc_0.c*role + nonarc_1.c*role + (1 + role | PTNUM), mvertex, REML=FALSE))
summary(m <- lmer(eccentricity_3 ~ narci_sidp_0.c*role + narci_sidp_1.c*role + nonarc_0.c*role + nonarc_1.c*role + (1 + role | PTNUM), mvertex, REML=FALSE))



summary(m <- lmer(degree_2 ~ narci_sidp_0 + narci_sidp_1 + (1 | PTNUM), mvertex))


cm <- lmerCellMeans(m, divide="narci_sidp_0", n.cont=8, divide.prefix=FALSE)
png("narci interact placeholder.png", width=6, height=5, units="in", res=150)
ggplot(cm, aes(x=narci_sidp_1, y=degree_2, color=narci_sidp_0, ymin=degree_2 - se, ymax=degree_2 + se)) + geom_line() + geom_errorbar() +
    ylab("Average degree centrality") + xlab("Patient Narcissistic Features") + scale_color_brewer("Partner\nNarcisstic Sx", palette="Dark2") +
    theme_bw(base_size=24)
dev.off()

summary(m <- lmer(degree_2 ~ narci_sidp_0.c*narci_sidp_1.c + nonarc_0.c*nonarc_1.c + (1 | PTNUM), mvertex))

summary(m <- lmer(degree_2 ~ bordl_sidp_0.c*bordl_sidp_1.c + nobpd_0.c*nobpd_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(degree_2 ~ pdtot_0.c*pdtot_1.c + (1 | PTNUM), mvertex))

summary(m <- lmer(degree_3 ~ narci_sidp_0.c*narci_sidp_1.c + (1 | PTNUM), mvertex))

cm <- lmerCellMeans(m, divide=c("narci_sidp_0.c", "narci_sidp_1.c"))

pdf("figures/narci_sidp degree interaction.pdf", width=9, height=9)
ggplot(cm, aes(x=narci_sidp_0.c, y=degree_2, color=narci_sidp_1.c, ymin=degree_2 - se, ymax=degree_2 + se)) + geom_pointrange(position=position_dodge(width=0.2))
dev.off()

mixed(degree_2 ~ narci_sidp_0.c*narci_sidp_1.c + (1 | PTNUM), mvertex)
mixed(degree_3 ~ narci_sidp_0.c*narci_sidp_1.c + (1 | PTNUM), mvertex)

summary(m <- lmer(degree_2 ~ narci_sidp_0.c*narci_sidp_1.c + iip_communion_0.c*iip_communion_1.c + iip_agency_0.c*iip_agency_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(degree_3 ~ narci_sidp_0.c*narci_sidp_1.c + iip_communion_0.c*iip_communion_1.c + iip_agency_0.c*iip_agency_1.c + (1 | PTNUM), mvertex))

summary(m <- lmer(degree_2 ~ OPD_sidp_0.c*OPD_sidp_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(degree_3 ~ OPD_sidp_0.c*OPD_sidp_1.c + (1 | PTNUM), mvertex))

summary(m <- lmer(degree_2 ~ OPD_sidp_0.c*OPD_sidp_1.c + iip_communion_0.c*iip_communion_1.c + iip_agency_0.c*iip_agency_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(degree_3 ~ OPD_sidp_0.c*OPD_sidp_1.c + iip_communion_0.c*iip_communion_1.c + iip_agency_0.c*iip_agency_1.c + (1 | PTNUM), mvertex))

summary(m <- lmer(degree_2 ~ IIP_PD1_0.c*IIP_PD1_1.c + IIP_PD2_0.c*IIP_PD2_1.c + IIP_PD3_0.c*IIP_PD3_1.c + (1 | PTNUM), mvertex))

#strength
summary(m <- lmer(strength_2 ~ narci_sidp_0.c*narci_sidp_1.c + bordl_sidp_0.c*bordl_sidp_1.c + (1 | PTNUM), mvertex)) #2 makes more sense for strength since it's continuous
summary(m <- lmer(strength_3 ~ narci_sidp_0.c*narci_sidp_1.c + bordl_sidp_0.c*bordl_sidp_1.c + (1 | PTNUM), mvertex))


summary(m <- lmer(strength_2 ~ narci_sidp_0.c*narci_sidp_1.c + (1 | PTNUM), mvertex)) #2 makes more sense for strength since it's continuous
cm <- lmerCellMeans(m, divide=c("narci_sidp_0.c", "narci_sidp_1.c"))

pdf("figures/narci_sidp strength interaction.pdf", width=9, height=9)
ggplot(cm, aes(x=narci_sidp_0.c, y=strength_2, color=narci_sidp_1.c, ymin=strength_2 - se, ymax=strength_2 + se)) + geom_pointrange(position=position_dodge(width=0.2))
dev.off()

#partner agency still associated with greater strength
summary(m <- lmer(strength_2 ~ iip_communion_0.c*iip_communion_1.c + iip_agency_0.c*iip_agency_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(strength_3 ~ iip_communion_0.c*iip_communion_1.c + iip_agency_0.c*iip_agency_1.c + (1 | PTNUM), mvertex))


#eigenvector 
summary(m <- lmer(evcent_binary_2 ~ bordl_sidp_0.c*bordl_sidp_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(evcent_binary_3 ~ bordl_sidp_0.c*bordl_sidp_1.c + (1 | PTNUM), mvertex))

summary(m <- lmer(evcent_binary_2 ~ PAIBORtot_0.c*PAIBORtot_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(evcent_binary_3 ~ PAIBORtot_0.c*PAIBORtot_1.c + (1 | PTNUM), mvertex))

summary(m <- lmer(evcent_weighted_2 ~ PAIBORtot_0.c*PAIBORtot_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(evcent_weighted_3 ~ PAIBORtot_0.c*PAIBORtot_1.c + (1 | PTNUM), mvertex))

summary(m <- lmer(evcent_weighted_2 ~ narci_sidp_0.c*narci_sidp_1.c + bordl_sidp_0.c*bordl_sidp_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(evcent_weighted_3 ~ narci_sidp_0.c*narci_sidp_1.c + bordl_sidp_0.c*bordl_sidp_1.c + (1 | PTNUM), mvertex))

summary(m <- lmer(evcent_binary_2 ~ narci_sidp_0.c*narci_sidp_1.c + bordl_sidp_0.c*bordl_sidp_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(evcent_binary_3 ~ narci_sidp_0.c*narci_sidp_1.c + bordl_sidp_0.c*bordl_sidp_1.c + (1 | PTNUM), mvertex))

summary(m <- lmer(evcent_weighted_2 ~ iip_communion_0.c*iip_communion_1.c + iip_agency_0.c*iip_agency_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(evcent_weighted_3 ~ iip_communion_0.c*iip_communion_1.c + iip_agency_0.c*iip_agency_1.c + (1 | PTNUM), mvertex))

#closeness
summary(m <- lmer(closeness_binary_2 ~ bordl_sidp_0.c*bordl_sidp_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(closeness_binary_3 ~ bordl_sidp_0.c*bordl_sidp_1.c + (1 | PTNUM), mvertex))

summary(m <- lmer(closeness_binary_2 ~ PAIBORtot_0.c*PAIBORtot_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(closeness_binary_3 ~ PAIBORtot_0.c*PAIBORtot_1.c + (1 | PTNUM), mvertex))

summary(m <- lmer(closeness_weighted_2 ~ PAIBORtot_0.c*PAIBORtot_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(closeness_weighted_3 ~ PAIBORtot_0.c*PAIBORtot_1.c + (1 | PTNUM), mvertex))

summary(m <- lmer(closeness_weighted_2 ~ narci_sidp_0.c*narci_sidp_1.c + bordl_sidp_0.c*bordl_sidp_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(closeness_weighted_3 ~ narci_sidp_0.c*narci_sidp_1.c + bordl_sidp_0.c*bordl_sidp_1.c + (1 | PTNUM), mvertex))

summary(m <- lmer(closeness_binary_2 ~ narci_sidp_0.c*narci_sidp_1.c + bordl_sidp_0.c*bordl_sidp_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(closeness_binary_3 ~ narci_sidp_0.c*narci_sidp_1.c + bordl_sidp_0.c*bordl_sidp_1.c + (1 | PTNUM), mvertex))

summary(m <- lmer(closeness_weighted_2 ~ iip_communion_0.c*iip_communion_1.c + iip_agency_0.c*iip_agency_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(closeness_weighted_3 ~ iip_communion_0.c*iip_communion_1.c + iip_agency_0.c*iip_agency_1.c + (1 | PTNUM), mvertex))


