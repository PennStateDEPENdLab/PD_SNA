library(lme4)
library(igraph)
library(afex)
library(ggplot2)
setwd(file.path(getMainDir(), "BPD_Couples", "Baseline_SNA"))
source(file.path(getMainDir(), "Miscellaneous", "Global_Functions.R"))
load(file="SNA_Processed_6Apr2015.RData")
load(file="ego_graph_measures.RData")

mgraph <- merge(graph_agg, other, by=c("UsrID", "PTNUM"), all.x=TRUE)
mvertex <- merge(vertex_agg, other, by=c("UsrID", "PTNUM"), all.x=TRUE)

mgraph$UsrID <- factor(mgraph$UsrID)
mgraph$PTNUM <- factor(mgraph$PTNUM)
mgraph$PAIBORtot.c <- mgraph$PAIBORtot - mean(mgraph$PAIBORtot, na.rm=TRUE)
mgraph$iip_agency.c <- mgraph$iip_agency - mean(mgraph$iip_agency, na.rm=TRUE)
mgraph$iip_communion.c <- mgraph$iip_communion - mean(mgraph$iip_communion, na.rm=TRUE)
mgraph$IIP_PD1.c <- mgraph$IIP_PD1 - mean(mgraph$IIP_PD1, na.rm=TRUE)
mgraph$IIP_PD2.c <- mgraph$IIP_PD2 - mean(mgraph$IIP_PD2, na.rm=TRUE)
mgraph$IIP_PD3.c <- mgraph$IIP_PD3 - mean(mgraph$IIP_PD3, na.rm=TRUE)
mgraph$ECRanx.c <- mgraph$ECRanx - mean(mgraph$ECRanx, na.rm=TRUE)
mgraph$ECRavoid.c <- mgraph$ECRavoid - mean(mgraph$ECRavoid, na.rm=TRUE)


mvertex$UsrID <- factor(mvertex$UsrID)
mvertex$PTNUM <- factor(mvertex$PTNUM)
mvertex$PAIBORtot.c <- mvertex$PAIBORtot - mean(mvertex$PAIBORtot, na.rm=TRUE)
mvertex$iip_agency.c <- mvertex$iip_agency - mean(mvertex$iip_agency, na.rm=TRUE)
mvertex$iip_communion.c <- mvertex$iip_communion - mean(mvertex$iip_communion, na.rm=TRUE)
mvertex$IIP_PD1.c <- mvertex$IIP_PD1 - mean(mvertex$IIP_PD1, na.rm=TRUE)
mvertex$IIP_PD2.c <- mvertex$IIP_PD2 - mean(mvertex$IIP_PD2, na.rm=TRUE)
mvertex$IIP_PD3.c <- mvertex$IIP_PD3 - mean(mvertex$IIP_PD3, na.rm=TRUE)
mvertex$ECRanx.c <- mvertex$ECRanx - mean(mvertex$ECRanx, na.rm=TRUE)
mvertex$ECRavoid.c <- mvertex$ECRavoid - mean(mvertex$ECRavoid, na.rm=TRUE)


#analysis of graph-wide metrics

#graph centralization
m <- lmer(centralization_degree_2 ~ PAIBORtot + (1 | PTNUM), mgraph)
summary(m)

m <- lmer(centralization_degree_3 ~ PAIBORtot + (1 | PTNUM), mgraph)
summary(m)

summary(m <- lmer(centralization_eigenvector_2 ~ PAIBORtot + (1 | PTNUM), mgraph))
summary(m <- lmer(centralization_eigenvector_2 ~ PAIBORtot.c*iip_agency.c*iip_communion.c + (1 | PTNUM), mgraph))
summary(m <- lmer(centralization_eigenvector_2 ~ IIP_PD1.c + (1 | PTNUM), mgraph))
summary(m <- lmer(centralization_eigenvector_2 ~ IIP_PD2.c + (1 | PTNUM), mgraph))
summary(m <- lmer(centralization_eigenvector_2 ~ IIP_PD3.c + (1 | PTNUM), mgraph))

mixed(centralization_eigenvector_2 ~ PAIBORtot.c*iip_agency.c*iip_communion.c + (1 | PTNUM), data=mgraph, method="KR")

mixed(centralization_eigenvector_2 ~ PAIBORtot.c*IIP_PD1.c*IIP_PD2.c*IIP_PD3.c + (1 | PTNUM), data=mgraph, method="KR")
mixed(centralization_eigenvector_2 ~ PAIBORtot.c*IIP_PD2.c + (1 | PTNUM), data=mgraph, method="KR")
summary(m)

cm <- lmerCellMeans(m, divide=c("iip_agency.c", "iip_communion.c"), n.cont=10)
ggplot(cm, aes(x=PAIBORtot.c, y=centralization_eigenvector_2, color=iip_agency.c, shape=iip_communion.c)) + geom_point() + geom_line()

m <- lmer(centralization_eigenvector_3 ~ PAIBORtot + (1 | PTNUM), mgraph)
summary(m)


#angry assortativity
summary(m <- lmer(angry_assortativity ~ PAIBORtot.c + (1 | PTNUM), mgraph))
summary(m <- lmer(angry_assortativity ~ iip_agency.c*iip_communion.c + (1 | PTNUM), mgraph))
summary(m <- lmer(angry_assortativity ~ IIP_PD1.c*IIP_PD2.c*IIP_PD3.c + (1 | PTNUM), mgraph))
summary(m <- lmer(angry_assortativity ~ PAIBORtot.c*ECRanx.c*ECRavoid.c + (1 | PTNUM), mgraph))
summary(m <- lmer(angry_assortativity ~ ECRanx.c + (1 | PTNUM), mgraph))
summary(m <- lmer(angry_assortativity ~ ECRavoid.c + (1 | PTNUM), mgraph))

#attach assortativity
summary(m <- lmer(attach_assortativity ~ PAIBORtot.c + (1 | PTNUM), mgraph))
summary(m <- lmer(attach_assortativity ~ iip_agency.c*iip_communion.c + (1 | PTNUM), mgraph))
summary(m <- lmer(attach_assortativity ~ IIP_PD1.c*IIP_PD2.c*IIP_PD3.c + (1 | PTNUM), mgraph))
summary(m <- lmer(attach_assortativity ~ PAIBORtot.c*ECRanx.c*ECRavoid.c + (1 | PTNUM), mgraph))
mixed(attach_assortativity ~ PAIBORtot.c*ECRanx.c*ECRavoid.c + (1 | PTNUM), mgraph)
summary(m <- lmer(attach_assortativity ~ PAIBORtot.c*ECRanx.c + (1 | PTNUM), mgraph))

cm <- lmerCellMeans(m, divide=c("ECRanx.c"))
ggplot(cm, aes(x=PAIBORtot.c, y=attach_assortativity, color=ECRanx.c)) + geom_line()

summary(m <- lmer(attach_assortativity ~ ECRanx.c + (1 | PTNUM), mgraph))
summary(m <- lmer(attach_assortativity ~ ECRavoid.c + (1 | PTNUM), mgraph))
hist(mgraph$attach_assortativity)

summary(m <- lmer(yearsknown_assortativity ~ PAIBORtot.c + (1 | PTNUM), mgraph))
summary(m <- lmer(yearsknown_assortativity ~ ECRanx.c + (1 | PTNUM), mgraph))
summary(m <- lmer(yearsknown_assortativity ~ ECRavoid.c + (1 | PTNUM), mgraph))
hist(mgraph$yearsknown_assortativity)

summary(m <- lmer(largest_clique_size_2 ~ PAIBORtot.c + (1 | PTNUM), mgraph))
summary(m <- lmer(largest_clique_size_3 ~ PAIBORtot.c + (1 | PTNUM), mgraph))

#largest clique is smaller with BPD Sx
mixed(largest_clique_size_2 ~ PAIBORtot.c + (1 | PTNUM), mgraph, method="KR")
mixed(largest_clique_size_3 ~ PAIBORtot.c + (1 | PTNUM), mgraph, method="KR")

mixed(largest_clique_size_2 ~ PAIBORtot.c*ECRanx.c + (1 | PTNUM), mgraph, method="KR")
mixed(largest_clique_size_3 ~ PAIBORtot.c*ECRanx.c + (1 | PTNUM), mgraph, method="KR")

mixed(largest_clique_size_3 ~ PAIBORtot.c*ECRanx.c + (1 | PTNUM), mgraph, method="KR")
mixed(largest_clique_size_3 ~ PAIBORtot.c*ECRavoid.c + (1 | PTNUM), mgraph, method="KR")

summary(m <- lmer(largest_clique_size_2 ~ iip_agency.c*iip_communion.c + (1 | PTNUM), mgraph))
summary(m <- lmer(largest_clique_size_2 ~ IIP_PD1.c + (1 | PTNUM), mgraph))
summary(m <- lmer(largest_clique_size_2 ~ IIP_PD2.c + (1 | PTNUM), mgraph))
summary(m <- lmer(largest_clique_size_2 ~ IIP_PD3.c + (1 | PTNUM), mgraph))
summary(m <- lmer(largest_clique_size_2 ~ IIP_PD1.c*IIP_PD2.c*IIP_PD3.c + (1 | PTNUM), mgraph))

summary(m <- lmer(largest_clique_size_2 ~ iip_agency.c*iip_communion.c + (1 | PTNUM), mgraph))
summary(m <- lmer(largest_clique_size_2 ~ IIP_PD1.c + (1 | PTNUM), mgraph))
summary(m <- lmer(largest_clique_size_2 ~ IIP_PD2.c + (1 | PTNUM), mgraph))
summary(m <- lmer(largest_clique_size_2 ~ IIP_PD3.c + (1 | PTNUM), mgraph))
summary(m <- lmer(largest_clique_size_2 ~ IIP_PD1.c*IIP_PD2.c*IIP_PD3.c + (1 | PTNUM), mgraph))

summary(lmer(num_maximal_3cliques_2 ~ PAIBORtot.c + (1 | PTNUM), mgraph))
summary(lmer(num_maximal_3cliques_2 ~ ECRanx.c +  (1 | PTNUM), mgraph)) #iip_communion.c + 
summary(lmer(num_maximal_3cliques_2 ~ iip_communion.c + (1 | PTNUM), mgraph))
summary(lmer(num_maximal_3cliques_2 ~ iip_agency.c + (1 | PTNUM), mgraph))

summary(lmer(num_maximal_3cliques_2 ~ iip_agency.c + (1 | PTNUM), mgraph))

summary(m <- lmer(num_maximal_3cliques_2 ~ iip_agency.c*iip_communion.c*PAIBORtot.c + (1 | PTNUM), mgraph))
mixed(num_maximal_3cliques_2 ~ iip_agency.c*iip_communion.c*PAIBORtot.c + (1 | PTNUM), mgraph, method="KR")
cm <- lmerCellMeans(m, divide=c("iip_agency.c", "PAIBORtot.c"))
ggplot(cm, aes(x=iip_communion.c, y=num_maximal_3cliques_2, color=iip_agency.c, shape=PAIBORtot.c)) + geom_line() + geom_point(size=3)


#holds at 3 threshold
summary(lmer(num_maximal_3cliques_3 ~ iip_communion.c + (1 | PTNUM), mgraph))
mixed(num_maximal_3cliques_3 ~ iip_communion.c + (1 | PTNUM), mgraph, method="KR")
summary(lmer(num_maximal_3cliques_3 ~ PAIBORtot.c + (1 | PTNUM), mgraph))

summary(lmer(num_maximal_3cliques_3 ~ iip_communion.c*PAIBORtot.c + (1 | PTNUM), mgraph))

#size of maximal cliques
summary(lmer(avgsize_maximal_3cliques_2 ~ iip_communion.c + (1 | PTNUM), mgraph))
summary(lmer(avgsize_maximal_3cliques_3 ~ iip_communion.c + (1 | PTNUM), mgraph))
summary(lmer(avgsize_maximal_3cliques_2 ~ iip_agency.c + (1 | PTNUM), mgraph))

#size of maximal cliques also smaller with BPD Sx
summary(lmer(avgsize_maximal_3cliques_2 ~ PAIBORtot.c + (1 | PTNUM), mgraph))
summary(lmer(avgsize_maximal_3cliques_3 ~ PAIBORtot.c + (1 | PTNUM), mgraph)) #holds at 3 threshold

#size of maximal cliques also smaller in those with high attachment anxiety (only at higher threshold)
summary(lmer(avgsize_maximal_3cliques_2 ~ ECRanx.c + (1 | PTNUM), mgraph))
summary(lmer(avgsize_maximal_3cliques_3 ~ ECRanx.c + (1 | PTNUM), mgraph))
summary(lmer(avgsize_maximal_3cliques_2 ~ ECRanx.c*PAIBORtot.c + (1 | PTNUM), mgraph))
summary(lmer(avgsize_maximal_3cliques_3 ~ ECRanx.c*PAIBORtot.c + (1 | PTNUM), mgraph)) #BPD diff explained/goes away when controlling for ECR anx
summary(lmer(avgsize_maximal_3cliques_3 ~ ECRanx.c + PAIBORtot.c + (1 | PTNUM), mgraph)) #BPD diff explained/goes away when controlling for ECR anx

summary(lmer(avgsize_maximal_3cliques_2 ~ PAIBORtot.c*iip_communion.c + (1 | PTNUM), mgraph))
summary(lmer(avgsize_maximal_3cliques_3 ~ PAIBORtot.c*iip_communion.c + (1 | PTNUM), mgraph)) #not explained by communion

#higher number of communities in BPD
summary(lmer(number_optcommunity_weighted_2 ~ PAIBORtot.c + (1 | PTNUM), mgraph))
summary(lmer(number_optcommunity_weighted_3 ~ PAIBORtot.c + (1 | PTNUM), mgraph))

#should correspond with smaller sizes... YES!
summary(lmer(avgsize_optcommunity_weighted_2 ~ PAIBORtot.c + (1 | PTNUM), mgraph))
summary(lmer(avgsize_optcommunity_weighted_3 ~ PAIBORtot.c + (1 | PTNUM), mgraph))


#weakly positive greater modularity in BPD (not significant)
summary(lmer(modularity_optcommunity_weighted_2 ~ PAIBORtot.c + (1 | PTNUM), mgraph))
summary(lmer(modularity_optcommunity_weighted_3 ~ PAIBORtot.c + (1 | PTNUM), mgraph))
summary(lmer(modularity_optcommunity_binary_2 ~ PAIBORtot.c + (1 | PTNUM), mgraph))
summary(lmer(modularity_optcommunity_binary_3 ~ PAIBORtot.c + (1 | PTNUM), mgraph))

#density differences?
summary(lmer(density_2 ~ PAIBORtot.c + (1 | PTNUM), mgraph))
summary(lmer(density_3 ~ PAIBORtot.c + (1 | PTNUM), mgraph))

#not really explained by ECR
summary(lmer(density_2 ~ PAIBORtot.c*ECRanx.c*ECRavoid.c + (1 | PTNUM), mgraph))
summary(lmer(density_3 ~ PAIBORtot.c*ECRanx.c*ECRavoid.c + (1 | PTNUM), mgraph))
summary(lmer(density_2 ~ PAIBORtot.c*iip_communion.c*iip_agency.c + (1 | PTNUM), mgraph))
summary(lmer(density_3 ~ PAIBORtot.c*iip_communion.c + PAIBORtot.c*iip_agency.c + (1 | PTNUM), mgraph))

#transitivity is weaker in BPD
summary(lmer(transitivity_2 ~ PAIBORtot.c + (1 | PTNUM), mgraph))
summary(lmer(transitivity_3 ~ PAIBORtot.c + (1 | PTNUM), mgraph))

summary(lmer(transitivity_2 ~ PAIBORtot.c*iip_communion.c + PAIBORtot.c*iip_agency.c + (1 | PTNUM), mgraph))
summary(lmer(transitivity_3 ~ PAIBORtot.c*iip_communion.c + PAIBORtot.c*iip_agency.c + (1 | PTNUM), mgraph))

#average path length not related to BPD
summary(lmer(avgpath_2 ~ PAIBORtot.c + (1 | PTNUM), mgraph))
summary(lmer(avgpath_3 ~ PAIBORtot.c + (1 | PTNUM), mgraph))

#but is related to communion
#so more communality means there are fewer hops to connect any 2 people
summary(lmer(avgpath_2 ~ iip_communion.c*iip_agency.c + (1 | PTNUM), mgraph))
summary(lmer(avgpath_2 ~ iip_communion.c + (1 | PTNUM), mgraph))
summary(lmer(avgpath_2 ~ iip_agency.c + (1 | PTNUM), mgraph))

#check small worldness
#have some infinite values for this...
#set to highest observed
mgraph$smallworld_2[which(mgraph$smallworld_2 == Inf)] <- max(mgraph$smallworld_2[!mgraph$smallworld_2==Inf])
mgraph$smallworld_3[which(mgraph$smallworld_3 == Inf)] <- max(mgraph$smallworld_3[!mgraph$smallworld_3==Inf])
summary(lmer(smallworld_2 ~ iip_communion.c*iip_agency.c + (1 | PTNUM), mgraph))
summary(lmer(smallworld_2 ~ PAIBORtot.c + (1 | PTNUM), mgraph)) #spurious?

summary(lmer(smallworld_3 ~ iip_communion.c*iip_agency.c + (1 | PTNUM), mgraph))
summary(lmer(smallworld_3 ~ PAIBORtot.c + (1 | PTNUM), mgraph))


#vertex analyses

#weaker degree in BPD, somewhat stronger with communion
summary(lmer(degree_2 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))
summary(lmer(degree_2 ~ PAIBORtot.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | UsrID) + (1 | PTNUM), mvertex))
summary(lmer(degree_2 ~ IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | UsrID) + (1 | PTNUM), mvertex))
summary(lmer(degree_2 ~ iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))

summary(lmer(degree_3 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))
summary(lmer(degree_3 ~ PAIBORtot.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | UsrID) + (1 | PTNUM), mvertex))
summary(lmer(degree_3 ~ iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))

summary(lmer(degree_3 ~ PAIBORtot.c + ECRanx.c + ECRavoid.c + (1 | UsrID) + (1 | PTNUM), mvertex))

#strength (weighted degree)
#basically same pattern (BPD weaker)
summary(lmer(strength_2 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))
summary(lmer(strength_3 ~ PAIBORtot.c + iip_communion.c + iip_agency.c + (1 | UsrID) + (1 | PTNUM), mvertex))

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

