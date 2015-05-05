library(lme4)
library(lavaan)
library(igraph)
library(foreign)
library(psych)
library(plyr)
library(ggplot2)
library(afex)
setwd(file.path(getMainDir(), "BPD_Couples", "Baseline_SNA"))
source(file.path(getMainDir(), "Miscellaneous", "Global_Functions.R"))
load(file="SNA_Processed_6Apr2015.RData")

#look at relationship satisfaction as a function of BPD etc.
other <- subset(other, PTNUM %in% unique(sna$PTNUM))

corstarsl(other[,c("DASCon", "DASSat", "DASCoh", "DASAffExp", "DASTotal", "PhysAssV", "PhysAssP", "PsychAggV", "PsychAggP", "Perp", "Victim")])
predstocenter <- c("PAIBORtot", "iip_agency", "iip_communion", "IIP_PD1", "IIP_PD2", "IIP_PD3", "ECRanx", "ECRavoid", "DASTotal", 
    "Victim", "Perp", "PsychAggV", "PsychAggP", "PhysAssV", "PhysAssP", "bordl_sidp", "narci_sidp", "antso_sidp", "avoid_sidp", "OPD_sidp", "depen_sidp")
other <- f_centerPredictors(other, predstocenter, addsuffix=".c")
summary(m <- lmer(DASCon ~ IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | PTNUM), other))
summary(m <- lmer(DASCon ~ IIP_PD1.c*IIP_PD2.c*IIP_PD3.c + (1 | PTNUM), other))
summary(m <- lmer(DASCon ~ PAIBORtot.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | PTNUM), other))
summary(m <- lmer(DASCon ~ PAIBORtot.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + Victim.c + Perp.c + (1 | PTNUM), other))

summary(m <- lmer(DASCon ~ PAIBORtot.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + PsychAggV.c + PsychAggP.c + (1 | PTNUM), other))
summary(m <- lmer(DASCon ~ PAIBORtot.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + PhysAssV.c + PhysAssP.c + (1 | PTNUM), other))

#narcissistic Sx, problems with aggression (IIP PD3), and being the victim of psychological and physical aggression/assault predict lower consensus
summary(m <- lmer(DASCon ~ narci_sidp.c + bordl_sidp.c + antso_sidp.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | PTNUM), other))

summary(m <- lmer(DASCon ~ narci_sidp.c + bordl_sidp.c + antso_sidp.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + Victim.c + Perp.c + (1 | PTNUM), other))

summary(m <- lmer(DASSat ~ IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | PTNUM), other))
summary(m <- lmer(DASSat ~ PAIBORtot.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | PTNUM), other)) #knock each other out
summary(m <- lmer(DASSat ~ narci_sidp.c + bordl_sidp.c + antso_sidp.c + avoid_sidp.c + (1 | PTNUM), other))
summary(m <- lmer(DASSat ~ narci_sidp.c + bordl_sidp.c + antso_sidp.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | PTNUM), other))

summary(m <- lmer(DASCoh ~ IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | PTNUM), other))
summary(m <- lmer(DASCoh ~ PAIBORtot.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | PTNUM), other)) #knock each other out
summary(m <- lmer(DASCoh ~ PAIBORtot.c + (1 | PTNUM), other))
summary(m <- lmer(DASCoh ~ PAIBORtot.c + bordl_sidp.c + OPD_sidp.c +  (1 | PTNUM), other))
summary(m <- lmer(DASCoh ~ PAIBORtot.c + OPD_sidp.c +  (1 | PTNUM), other))
summary(m <- lmer(DASCoh ~ bordl_sidp.c + OPD_sidp.c + (1 | PTNUM), other))
summary(m <- lmer(DASCoh ~ narci_sidp.c + bordl_sidp.c + antso_sidp.c + avoid_sidp.c + depen_sidp.c + (1 | PTNUM), other))
summary(m <- lmer(DASCoh ~ narci_sidp.c + bordl_sidp.c + antso_sidp.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + (1 | PTNUM), other))


#applies to total, too... effect of perp on greater relationship satisfaction goes away when victim removed (likely collinearity problem given .61 correlation between victim and perp)
summary(m <- lmer(DASTotal ~ narci_sidp.c + bordl_sidp.c + antso_sidp.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + Victim.c + Perp.c + (1 | PTNUM), other))

summary(m <- lmer(DASTotal ~ narci_sidp.c + bordl_sidp.c + antso_sidp.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + Perp.c + (1 | PTNUM), other))
summary(m <- lmer(DASTotal ~ narci_sidp.c + bordl_sidp.c + antso_sidp.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + Victim.c + (1 | PTNUM), other)) #this holds up
summary(m <- lmer(DASTotal ~ bordl_sidp.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + Victim.c + (1 | PTNUM), other)) #BPD Sx alone not informative

mixed(DASTotal ~ narci_sidp.c + bordl_sidp.c + antso_sidp.c + IIP_PD1.c + IIP_PD2.c + IIP_PD3.c + Victim.c + (1 | PTNUM), other, method="PB")


#look at whether attachment items form a single factor
cor(sna$Close_num, sna$attach_total, use="pairwise.complete.obs")
msyn <- 'attach=~prxsk_num + sepds_num + sfhvn_num + secbs_num
    sfhvn_num ~~ secbs_num'
mcfa <- cfa(msyn, sna, missing="ML", estimator="MLR")

summary(mcfa, fit.measures=TRUE)
standardizedSolution(mcfa, type="std.all")
modificationIndices(mcfa)


#anger as a function of BPD Sx
m <- lmer(Angry_num ~ PAIBORtot.c + (1 | UsrID), data=sna_merge, REML=TRUE)
summary(m)

#nested random effect: UsrID nested within PTNUM (Dyad)
xtabs(~ UsrID + PTNUM, sna_merge, drop = TRUE, sparse = TRUE)

#see here for specification of nested factors: need to simply have random intercepts as a function of both UsrID and PTNUM (similar to crossed)
#http://lme4.r-forge.r-project.org/book/Ch2.pdf

#Yes, anger toward members of the social network is more common in BPD
mnest <- lmer(Angry_num ~ PAIBORtot.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, REML=TRUE)
summary(mnest)

mnest <- lmer(Angry_num ~ bordl_sidp.c + OPD_sidp.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, REML=TRUE)
summary(mnest)


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

mcutoff <- glmer(Cutoff_num ~ ECRanx.c*PAIBORtot.c + ECRavoid.c*PAIBORtot.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=binomial, nAGQ=1)
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

mcutoff_close <- glmer(Cutoff_num ~ Close_num*PAIBORtot.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=binomial, nAGQ=1)

#yes, strong 3-way interaction
#interpretation: Those with low attachment avoidance and high BPD Sx tended to nominate more folks in their network who were distant and cutoff.
mcutoff_close <- glmer(Cutoff_num ~ Close_num*ECRavoid.c*PAIBORtot.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=binomial, nAGQ=1)
summary(mcutoff_close)


cm <- lmerCellMeans(mcutoff_close, divide=c("ECRavoid.c", "PAIBORtot.c"), n.cont=10)


mcutoff_close <- glmer(Cutoff_num ~ Close_num*ECRavoid.c + Close_num*ECRanx.c + (1 | UsrID) + (1 | PTNUM), data=sna_merge, family=binomial, nAGQ=1)


#
pdf("figures/Cutoff by ECR Avoid and PAI-BOR and Close.pdf", width=9, height=6)
ggplot(cm, aes(x=Close_num, y=Cutoff_prob, ymin=selow, ymax=sehigh, color=PAIBORtot.c)) + facet_wrap(~ECRavoid.c) + geom_pointrange(position=position_dodge(width=0.5), size=1.5) + 
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


pdf("figures/all_couple_graphs.pdf", width=25, height=25)
for (i in 1:length(couple_graphs)) {
  g <- graph.data.frame(couple_graphs[[i]]$edges_avgratings, directed=FALSE, vertices=couple_graphs[[i]]$actors)
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

    plot(g, layout=layout.fruchterman.reingold, edge.width=E(g)$weight*2, vertex.label.color="black", #vertex.size=V(g)$attach_total,
      vertex.color=V(g)$role_color, vertex.label.family="sans", vertex.label.cex=1.3, vertex.size=13, main=names(couple_graphs)[i])
  
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
