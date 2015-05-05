library(lme4)
library(igraph)
library(afex)
library(ggplot2)
library(reshape2)
setwd(file.path(getMainDir(), "BPD_Couples", "Baseline_SNA"))
source(file.path(getMainDir(), "Miscellaneous", "Global_Functions.R"))
load(file="SNA_Processed_6Apr2015.RData")
load(file="couple_graph_measures.RData")

other_melt <- melt(other[,c("PTNUM", "DyadID", "age", "DASCon", "DASSat", "DASCoh", "DASAffExp", "DASTotal", "ECRanx", "ECRavoid", 
            "PAIBORnegr", "PAIBORtot", "IIP_PD1", "IIP_PD2", "IIP_PD3", "iip_agency", "iip_communion", "Victim", "Perp", "bordl_sidp", "narci_sidp", "antso_sidp", "avoid_sidp", "OPD_sidp")], id.vars=c("PTNUM", "DyadID"))

other_wide <- dcast(other_melt, PTNUM ~ variable + DyadID)

mgraph <- merge(couple_graph_metrics, other_wide, by=c("PTNUM"), all.x=TRUE)
mvertex <- merge(couple_vertex_metrics, other_wide, by=c("PTNUM"), all.x=TRUE)

mgraph$PTNUM <- factor(mgraph$PTNUM)
mvertex$PTNUM <- factor(mvertex$PTNUM)
predstocenter <- c("PAIBORtot", "iip_agency", "iip_communion", "IIP_PD1", "IIP_PD2", "IIP_PD3", "ECRanx", "ECRavoid", "DASTotal", "Victim", "Perp", "bordl_sidp", "narci_sidp", "antso_sidp", "avoid_sidp", "OPD_sidp")
mgraph <- f_centerPredictors(mgraph, apply(expand.grid(predstocenter, c("0", "1")), 1, paste, collapse="_"), addsuffix=".c")
mvertex <- f_centerPredictors(mvertex, apply(expand.grid(predstocenter, c("0", "1")), 1, paste, collapse="_"), addsuffix=".c")

#graph-wide analysis
cor.test(~ PAIBORtot_0.c + PAIBORtot_1.c, mgraph) #weak correlation
cor.test(~ iip_agency_0.c + iip_agency_1.c, mgraph) #weak correlation
cor.test(~ iip_communion_0.c + iip_communion_1.c, mgraph) #weak correlation
cor.test(~ ECRanx_0.c + ECRanx_1.c, mgraph) #weak correlation
cor.test(~ ECRavoid_0.c + ECRavoid_1.c, mgraph) #weak correlation
cor(mgraph[,c("IIP_PD1_0.c", "IIP_PD1_1.c", "IIP_PD2_0.c", "IIP_PD2_1.c", "IIP_PD3_0.c", "IIP_PD3_1.c")], use="pairwise.complete.obs")
cor(mgraph[,c("antso_sidp_0.c", "antso_sidp_1.c", "bordl_sidp_0.c", "bordl_sidp_1.c", "narci_sidp_0.c", "narci_sidp_1.c")], use="pairwise.complete.obs")

#partner BPD Sx contribute to weaker centralization
m <- lm(centralization_degree_2 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph)
summary(m)

#partner Sx contributes to weaker centralization, patient to greater.
m <- lm(centralization_degree_3 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph)
summary(m)

summary(m <- lm(centralization_degree_2 ~ ECRanx_0.c*ECRanx_1.c, mgraph))
summary(m <- lm(centralization_degree_3 ~ ECRanx_0.c*ECRanx_1.c, mgraph))

summary(m <- lm(centralization_eigenvector_2 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))
summary(m <- lm(centralization_eigenvector_3 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))
summary(m <- lm(centralization_eigenvector_2 ~ iip_agency_0.c*iip_agency_1.c + iip_communion_0.c*iip_communion_1.c, mgraph))
summary(m <- lm(centralization_eigenvector_3 ~ iip_agency_0.c*iip_agency_1.c + iip_communion_0.c*iip_communion_1.c, mgraph))
summary(m <- lm(centralization_eigenvector_3 ~ ECRanx_0.c*ECRanx_1.c + ECRavoid_0.c*ECRavoid_1.c, mgraph))
summary(m <- lm(centralization_eigenvector_3 ~ ECRanx_0.c*ECRanx_1.c, mgraph)) #significant interaction such that both partners' attachment anxiety contribute to greater centralization

#assortativity of nomination (patient, partner, or both)
#something messed up with this at the moment: assortativity all zero
#summary(m <- lm(nomination_assortativity_2 ~ ECRanx_0.c*ECRanx_1.c, mgraph)) #significant interaction such that both partners' attachment anxiety contribute to greater centralization

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
summary(m <- lm(density_2 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))
summary(m <- lm(density_3 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))
summary(m <- lm(density_2 ~ iip_communion_0.c*iip_communion_1.c, mgraph))
summary(m <- lm(density_3 ~ iip_communion_0.c*iip_communion_1.c, mgraph))
summary(m <- lm(density_2 ~ iip_agency_0.c*iip_agency_1.c, mgraph))
summary(m <- lm(density_3 ~ iip_agency_0.c*iip_agency_1.c, mgraph))
summary(m <- lm(density_2 ~ ECRanx_0.c*ECRanx_1.c, mgraph))
summary(m <- lm(density_3 ~ ECRavoid_0.c*ECRavoid_1.c, mgraph))

#transitivity
summary(m <- lm(transitivity_2 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))
summary(m <- lm(transitivity_3 ~ PAIBORtot_0.c*PAIBORtot_1.c, mgraph))
summary(m <- lm(transitivity_2 ~ iip_communion_0.c*iip_communion_1.c, mgraph))
summary(m <- lm(transitivity_3 ~ iip_communion_0.c*iip_communion_1.c, mgraph))
summary(m <- lm(transitivity_2 ~ iip_agency_0.c*iip_agency_1.c, mgraph))
summary(m <- lm(transitivity_3 ~ iip_agency_0.c*iip_agency_1.c, mgraph))
summary(m <- lm(transitivity_2 ~ ECRanx_0.c*ECRanx_1.c, mgraph))
summary(m <- lm(transitivity_3 ~ ECRavoid_0.c*ECRavoid_1.c, mgraph))

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

summary(m <- lmer(degree_2 ~ bordl_sidp_0.c*bordl_sidp_1.c + iip_communion_0.c*iip_communion_1.c + iip_agency_0.c*iip_agency_1.c + (1 | PTNUM), mvertex))
summary(m <- lmer(degree_3 ~ bordl_sidp_0.c*bordl_sidp_1.c + iip_communion_0.c*iip_communion_1.c + iip_agency_0.c*iip_agency_1.c + (1 | PTNUM), mvertex))

#greater narcissism in both partners associated with weaker degree
summary(m <- lmer(degree_2 ~ narci_sidp_0.c*narci_sidp_1.c + (1 | PTNUM), mvertex))
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


