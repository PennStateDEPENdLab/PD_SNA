#at the moment, self-reports and clinical data are spread across different datasets that are incremental, not complete.
library(foreign)
library(plyr)
library(reshape2)
library(dplyr) #load after plyr
setwd("~/Tresors/PD_SNA/Couples_Baseline_SNA")
source(file.path(getMainDir(), "Miscellaneous", "Global_Functions.R"))

#merge together self-reports etc.
#biq <- read.spss("data/selfreports/BIQ 3.31.15.sav", to.data.frame=TRUE)
#biq <- subset(biq, select=c(UsrID, PTNUM, DyadID, age, sex, race_n))
#biq_append <- read.spss("data/selfreports/BIQ_24Sep2015.sav", to.data.frame=TRUE) #NB: this is just an updated copy, not a full file
#biq_append <- plyr::rename(biq_append, c(gender="sex"))
#biq_append <- subset(biq_append, select=c(UsrID, PTNUM, DyadID, age, sex, race_n))
#biq_all <- rbind(biq, biq_append)
#biq_all <- biq_all[!duplicated(biq_all$UsrID, fromLast = FALSE), ] #drop duplicated UsrIDs, preferring earlier row (i.e., older data)

#Feb2017 (should be final sample
biq <- read.spss("data/selfreports/INTAKE_20170223/BIQ.sav", to.data.frame=TRUE)
biq <- biq %>% filter(mth==0) %>% select(UsrID, PTNUM, DyadID, age, sex, race)

#cts <- read.spss("data/selfreports/CTS intake scores 3.4.15.sav", to.data.frame=TRUE)
#cts_append <- read.spss("data/selfreports/CTS_24Sep2015.sav", to.data.frame=TRUE) #NB: this is also just updates to prior file
#cts_append <- plyr::rename(cts_append, c(dyad="DyadID"))
#cts_append2 <- read.spss("data/selfreports/CTS_incremental_8Oct2015.sav", to.data.frame=TRUE)
#cts_append2 <- plyr::rename(cts_append2, c(dyadID="DyadID"))
#cts_append2 <- subset(cts_append2, mth==0) #only baseline
#cts_append2 <- cts_append2[,names(cts_append)] #only keep scores
#cts_append3 <- read.spss("data/selfreports/CTS_incremental_9Oct2015.sav", to.data.frame=TRUE)
#cts_append3 <- plyr::rename(cts_append3, c(dyadID="DyadID"))
#cts_append3 <- subset(cts_append3, mth==0) #only baseline
#cts_append3 <- cts_append3[,names(cts_append)] #only keep scores
#cts_all <- rbind(cts, cts_append, cts_append2, cts_append3) 
#cts_all <- cts_all[!duplicated(cts_all$UsrID, fromLast = FALSE), ] #drop duplicated UsrIDs, preferring earlier row (i.e., older data)

#raw data
cts <- read.spss("data/selfreports/INTAKE_20170223/CTS.sav", to.data.frame=TRUE)

#VARIABLE LABELS 
#PsychAggV 'CTS psychological aggression-victim',
#PsychAggP 'CTS psychological aggression-perpetrator',
#PhysAssV 'CTS physical assault-victim',
#PhysAssP 'CTS physical assault-perpetrator',
#PsychAggVM 'CTS psychological aggression-victim (minor)',
#PsychAggPM 'CTS psychological aggression-perpetrator (minor)',
#PsychAggVS 'CTS psychological aggression-victim (severe)',
#PsychAggPS 'CTS psychological aggression-perpetrator (severe)',
#PhysAssVM 'CTS physical assault-victim (minor)',
#PhysAssPM 'CTS physical assault-perpetrator (minor)',
#PhysAssVS 'CTS physical assault-victim (severe)',
#PhysAssPS 'CTS physical assault-perpetrator (severe)',
#Victim 'CTS victim total',
#Perp 'CTS perpetrator total'.

#recode CTS to represent counts (uneven anchor spacing)
#DEFINITIONS:  v=victim, p=perpetrator, m=minor, s=severe
cts <- cts %>% mutate_at(vars(CTS1:CTS40), funs(recode(. , `0`=0, `1`=1, `2`=2, `3`=4, `4`=8, `5`=15, `6`=25, `7`=98))) %>%
    rename(DyadID=dyadID) %>% #capitalize D for match
    mutate(
        CTS_PsychAggV = CTS1 + CTS11 + CTS15 + CTS19 + CTS27 + CTS33 + CTS35 + CTS37,
        CTS_PsychAggP = CTS2 + CTS12 + CTS16 + CTS20 + CTS28 + CTS34 + CTS36 + CTS38,
        CTS_PhysAssV = CTS3 + CTS5 + CTS7 + CTS9 + CTS13 + CTS17 + CTS21 + CTS23 + CTS25 + CTS29 + CTS31 + CTS39,
        CTS_PhysAssP = CTS4 + CTS6 + CTS8 + CTS10 + CTS14 + CTS18 + CTS22 + CTS24 + CTS26 + CTS30 + CTS32 + CTS40,
        CTS_PsychAggVM = CTS1 + CTS19 + CTS27 + CTS35,
        CTS_PsychAggPM = CTS2 + CTS20 + CTS28 + CTS36,
        CTS_PsychAggVS = CTS11 + CTS15 + CTS33 + CTS37,
        CTS_PsychAggPS = CTS12 + CTS16 + CTS34 + CTS38,
        CTS_PhysAssVM = CTS3 + CTS5 + CTS7 + CTS25 + CTS29,
        CTS_PhysAssPM = CTS4 + CTS6 + CTS8 + CTS26 + CTS30,
        CTS_PhysAssVS = CTS9 + CTS13 + CTS17 + CTS21 + CTS23 + CTS31 + CTS39,
        CTS_PhysAssPS = CTS10 + CTS14 + CTS18 + CTS22 + CTS24 + CTS32 + CTS40,
        CTS_Victim = CTS_PsychAggV + CTS_PhysAssV,
        CTS_Perp = CTS_PsychAggP + CTS_PhysAssP
)

cts <- cts %>% filter(mth==0 & PTNUM >= 8000) %>% select(-one_of(c("ADate", "RecNum", "mth", paste0("CTS", 1:40))))


#das <- read.spss("data/selfreports/DAS intake scores 3.4.15.sav", to.data.frame=TRUE)
#das_append <- read.spss("data/selfreports/DAS_24Sep2015.sav", to.data.frame=TRUE)
#das_append <- plyr::rename(das_append, c(Dyad="DyadID"))
#das_append$DASrelationship <- factor(das_append$DASrelationship, levels=c(0,1), labels=c("no", "yes"))
#das_append2 <- read.spss("data/selfreports/DAS_incremental_8Oct2015.sav", to.data.frame=TRUE) #LOOKS LIKE THERE ARE TWO ENTRIES FOR ONE PERSON. MTH ISSUE
#das_append2 <- subset(das_append2, mth==0)
#das_append2$DASrelationship <- factor(das_append2$DASrelationship, levels=c(0,1), labels=c("no", "yes"))
#das_append2 <- das_append2[,names(das_append)] #only keep scores
#das_append3 <- read.spss("data/selfreports/DAS_incremental_9Oct2015.sav", to.data.frame=TRUE) #LOOKS LIKE THERE ARE TWO ENTRIES FOR ONE PERSON. MTH ISSUE
#das_append3 <- subset(das_append3, mth==0)
#das_append3$DASrelationship <- factor(das_append3$DASrelationship, levels=c(0,1), labels=c("no", "yes"))
#das_append3 <- das_append3[,names(das_append)] #only keep scores
#das_all <- rbind(das, das_append, das_append2, das_append3)
#das_all <- das_all[!duplicated(das_all$UsrID, fromLast = FALSE), ] #drop duplicated UsrIDs, preferring earlier row (i.e., older data)
#
##has 8888 as missing code
#lapply(das_all, function(x) { if (is.numeric(x)) { range(x, na.rm=TRUE) } } )
#
#recodeMissing <- function(vec, targetval=99) {
#  vec[vec==targetval] <- NA
#  return(vec)
#}
#
#das_all <- colwise(recodeMissing, targetval=8888)(das_all)

das <- read.spss("data/selfreports/INTAKE_20170223/DAS.sav", to.data.frame=TRUE)

reverseItems <- paste("DAS", c(1:15, 18, 19, 32), sep="") 
das[,reverseItems] <- lapply(das[,reverseItems], function(x) { 6 - x }) #reverse score and subtract 1 (becomes 0-5)

reverseItems <- paste("DAS", c(23, 24), sep="") 
das[,reverseItems] <- lapply(das[,reverseItems], function(x) { 5 - x }) #reverse score and subtract 1 (becomes 0-4)

subtractItems <- paste("DAS", c(16, 17, 20:22, 25:28, 31), sep="")
das$DAS20[das$DAS20 == 7] <- NA #a '7' in the original coding indicates they are not married/live together 
das[,subtractItems] <- lapply(das[,subtractItems], function(x) { x - 1 })

flipItems <- paste("DAS", c(29, 30), sep="")
das[,flipItems] <- lapply(das[,flipItems], function(x) { ifelse(x==1, 0, 1) })

#if DASrelationship is 0, then the person is not in a relationship with their partner (or at least the one at enrollment)
#invalidate any data
das <- das %>% mutate_at(vars(one_of(paste0("DAS", 1:32))), funs(ifelse(DASrelationship==0, NA, .)))

das <- das %>% mutate(
    DASCon=DAS1 + DAS2 + DAS3 + DAS5 + DAS7 + DAS8 + DAS9 + DAS10 + DAS11 + DAS12 + DAS13 + DAS14 + DAS15, #consensus
    DASSat=DAS16 + DAS17 + DAS18 + DAS19 + DAS20 + DAS21 + DAS22 + DAS23 + DAS31 + DAS32, #satisfaction
    DASCoh=DAS24 + DAS25 + DAS26 + DAS27 + DAS28, #cohesion
    DASAffExp=DAS4 + DAS6 + DAS29 + DAS30, #affectional expression
    DASTotal=DAS1 + DAS2 + DAS3 + DAS4 + DAS5 + DAS6 + DAS7 + DAS8 + DAS9 + DAS10 + DAS11 + DAS12 + 
        DAS13 + DAS14 + DAS15 + DAS16 + DAS17 + DAS18 + DAS19 + DAS20 + DAS21 + DAS22 + DAS23 + DAS24 +
        DAS25 + DAS26 + DAS27 + DAS28 + DAS29 + DAS30 + DAS31 + DAS32
)

das <- das %>% filter(mth==0 & PTNUM >= 8000) %>% select(-one_of(c("ADate", "RecNum", "mth", paste0("DAS", 1:32))))

#ecr <- read.spss("data/selfreports/ECR intake scores 3.4.15.sav", to.data.frame=TRUE)
#ecr_append <- read.spss("data/selfreports/ECR_24Sep2015.sav", to.data.frame=TRUE)
#ecr_append <- plyr::rename(ecr_append, c(Dyad="DyadID"))
#ecr_append2 <- read.spss("data/selfreports/ECR_incremental_8Oct2015.sav", to.data.frame=TRUE)
#ecr_append2 <- subset(ecr_append2, mth==0)
#ecr_append2 <- ecr_append2[,names(ecr_append)] #only keep scores
#ecr_append3 <- read.spss("data/selfreports/ECR_incremental_9Oct2015.sav", to.data.frame=TRUE)
#ecr_append3 <- subset(ecr_append3, mth==0)
#ecr_append3 <- ecr_append3[,names(ecr_append)] #only keep scores
#ecr_all <- rbind(ecr, ecr_append, ecr_append2, ecr_append3)
#ecr_all <- ecr_all[!duplicated(ecr_all$UsrID, fromLast = FALSE), ] #drop duplicated UsrIDs, preferring earlier row (i.e., older data)

ecr <- read.spss("data/selfreports/INTAKE_20170223/ECR.sav", to.data.frame=TRUE)
xx <- score_ecr(ecr, item_prefix="ECR", max_impute = 0.5, keep_reverse_codes = TRUE, drop_items = FALSE)

reverseItems <- paste("ECR", c(2, 22, 3, 5, 11,15, 17, 19, 25, 27, 29, 31, 33, 35), sep="")
ecr[,reverseItems] <- lapply(ecr[,reverseItems], function(x) { 8 - x }) #1-7 scoring
ecr <- ecr %>% rowwise() %>% mutate(ECRanx=mean(c(ECR2, ECR4, ECR6, ECR8, ECR10, ECR12, ECR14,  ECR16, ECR18, ECR20, ECR22, ECR24, ECR26, ECR28, ECR30, ECR32, ECR34, ECR36)),
    ECRavoid=mean(c(ECR1, ECR3, ECR5, ECR7, ECR9, ECR11, ECR13, ECR15, ECR17, ECR19, ECR21, ECR23, ECR25, ECR27, ECR29, ECR31, ECR33, ECR35)))

ecr <- ecr %>% filter(mth==0 & PTNUM >= 8000) %>% select(-one_of(c("ADate", "RecNum", "mth", paste0("ECR", 1:36))))

#pai <- read.spss("data/selfreports/PAI intake scores 3.4.15.sav", to.data.frame=TRUE)
#pai_append <- read.spss("data/selfreports/PAIBOR_24Sep2015.sav", to.data.frame=TRUE)
#pai_append <- plyr::rename(pai_append, c(Dyad="DyadID"))
#pai_append2 <- read.spss("data/selfreports/PAIBOR_incremental_8Oct2015.sav", to.data.frame=TRUE)
#pai_append2 <- subset(pai_append2, mth==0)
#pai_append2 <- pai_append2[,names(pai_append)]
#pai_append3 <- read.spss("data/selfreports/PAIBOR_incremental_9Oct2015.sav", to.data.frame=TRUE)
#pai_append3 <- subset(pai_append3, mth==0)
#pai_append3 <- pai_append3[,names(pai_append)]
#pai_all <- rbind(pai, pai_append, pai_append2, pai_append3)
#pai_all <- pai_all[!duplicated(pai_all$UsrID, fromLast = FALSE), ] #drop duplicated UsrIDs, preferring earlier row (i.e., older data)

pai <- read.spss("data/selfreports/INTAKE_20170223/PAIBOR.sav", to.data.frame=TRUE)
reverseItems <- paste("PAIBOR", c(7, 12, 14, 19, 20, 24), sep="")
pai[,reverseItems] <- lapply(pai[,reverseItems], function(x) { 3 - x }) #0-3 scoring
pai$PAIBORaffe <- with(pai, PAIBOR1 + PAIBOR4 + PAIBOR7 + PAIBOR10 + PAIBOR14 + PAIBOR18)
pai$PAIBORiden <- with(pai, PAIBOR2 + PAIBOR5 + PAIBOR8 + PAIBOR11 + PAIBOR15 + PAIBOR19)
pai$PAIBORnegr <- with(pai, PAIBOR3 + PAIBOR6 + PAIBOR9 + PAIBOR12 + PAIBOR16 + PAIBOR20)
pai$PAIBORharm <- with(pai, PAIBOR13 + PAIBOR17 + PAIBOR21 + PAIBOR22 + PAIBOR23 + PAIBOR24)
pai$PAIBORtot <- with(pai, PAIBORaffe + PAIBORiden + PAIBORnegr + PAIBORharm)
pai <- pai %>% filter(mth==0 & PTNUM >= 8000) %>% select(-one_of(c("ADate", "RecNum", "mth", paste0("PAIBOR", 1:24))))

##iip <- read.spss("data/IIP intake scores.sav", to.data.frame=TRUE)
#iip <- read.spss("data/selfreports/IIP90 intake.sav", to.data.frame=TRUE)
#iip <- subset(iip, select=c(UsrID, PTNUM, DyadID, IIP_PD1, IIP_PD2, IIP_PD3, IIP_C1, IIP_C2, iip_pa,
#        iip_bc, iip_de, iip_fg, iip_hi, iip_jk, iip_lm, iip_no, iip_bpd, PD_on_IIP))
#iip_append <- read.spss("data/selfreports/IIP_24Sep2015.sav", to.data.frame=TRUE)
#iip_append <- plyr::rename(iip_append, c(Dyad="DyadID"))
#iip_append$IIPSumValue <- NULL
#iip_append2 <- read.spss("data/selfreports/IIP90_incremental_8Oct2015.sav", to.data.frame=TRUE)
#iip_append2 <- subset(iip_append2, mth==0)
#iip_append2 <- iip_append2[,names(iip_append)]
#iip_all <- rbind(iip, iip_append, iip_append2)
#iip_all <- iip_all[!duplicated(iip_all$UsrID, fromLast = FALSE), ] #drop duplicated UsrIDs, preferring earlier row (i.e., older data)

iip <- read.spss("data/selfreports/INTAKE_20170223/IIP90.sav", to.data.frame=TRUE)
iip <- iip %>% rowwise() %>% mutate(
    IIP_pa = mean(c(IIP21,IIP40,IIP57,IIP58,IIP65,IIP68,IIP76,IIP80)),
    IIP_bc = mean(c(IIP1,IIP26,IIP28,IIP38,IIP41,IIP50,IIP73,IIP88)),
    IIP_de = mean(c(IIP11,IIP18,IIP20,IIP24,IIP27,IIP31,IIP46,IIP82)),
    IIP_fg = mean(c(IIP3,IIP7,IIP17,IIP22,IIP43,IIP45,IIP71,IIP85)),
    IIP_hi = mean(c(IIP5,IIP6,IIP8,IIP9,IIP12,IIP15,IIP23,IIP49)),
    IIP_jk = mean(c(IIP2,IIP10,IIP29,IIP44,IIP48,IIP54,IIP69,IIP83)),
    IIP_lm = mean(c(IIP25,IIP37,IIP47,IIP59,IIP64,IIP67,IIP70,IIP87)),
    IIP_no = mean(c(IIP4,IIP30,IIP39,IIP52,IIP56,IIP61,IIP62,IIP78)),
    IIP_bpd = mean(c(IIP51,IIP53,IIP55,IIP66,IIP77,IIP80,IIP89,IIP90)), #clifton bpd scale
    IIP_pd1 = mean(c(IIP1,IIP35,IIP36,IIP42,IIP51,IIP55,IIP60,IIP78,IIP79,IIP81,IIP86)), #interpersonal sensitivity
    IIP_pd2 = mean(c(IIP13,IIP14,IIP26,IIP28,IIP32,IIP34,IIP38,IIP40,IIP41,IIP84)), #interpersonal ambivalence
    IIP_pd3 = mean(c(IIP50,IIP53,IIP58,IIP63,IIP77,IIP80,IIP88)), #aggression
    IIP_pd = mean(c(IIP_pd1, IIP_pd2, IIP_pd3)), #overall pd; pd_on_IIP = IIP_pd > 1.1
    IIP_havePD = as.numeric(IIP_pd > 1.1),
    IIP_c1 = mean(c(IIP2,IIP9,IIP16,IIP48,IIP59,IIP66,IIP72,IIP74,IIP75)), #need for social approval
    IIP_c2 = mean(c(IIP3,IIP7,IIP17,IIP19,IIP22,IIP33,IIP43,IIP49,IIP71,IIP85)), #lack of sociability
    IIP_c = mean(c(IIP_c1, IIP_c2)),
    IIP_agency = .25*(IIP_pa - IIP_hi + .707*(IIP_bc + IIP_no - IIP_fg - IIP_jk)),
    IIP_communion = .25*(IIP_lm - IIP_de + .707*(IIP_no + IIP_jk - IIP_bc - IIP_fg)),
    IIP_elevation = (IIP_pa + IIP_bc + IIP_de + IIP_fg + IIP_hi + IIP_jk + IIP_lm + IIP_no)/8
)

iip <- iip %>% filter(mth==0 & PTNUM >= 8000) %>% select(-one_of(c("ADate", "RecNum", "mth", paste0("IIP", 1:90))))

#sidp <- read.spss("data/selfreports/SIDP4 intake consensus scored 4.2.15.sav", to.data.frame=TRUE)
#sidp_append <- read.spss("data/selfreports/SIDP_24Sep2015.sav", to.data.frame=TRUE)
#sidp_append$NumMiss <- NULL
#sidp_all <- rbind(sidp, sidp_append)
#sidp_all <- sidp_all[!duplicated(sidp_all$UsrID, fromLast = FALSE), ] #drop duplicated UsrIDs, preferring earlier row (i.e., older data)
#new dataset has all individuals at baseline
#sidp <- read.spss("data/selfreports/SIDP_scored.10.8.15.sav", to.data.frame=TRUE)
#sidp <- subset(sidp, mth==0 & raterID==0 & ratercode==6) #raterID 0 is case conference; ratercode 6 is case conference; mth 0 is intake
#sidp_all <- sidp[, c("UsrID", "PTNUM", "DyadID", grep(".*_sidp$", names(sidp), value=TRUE), grep(".*Count$", names(sidp), value=TRUE))]

#this is the (final) update from Nate; the above is missing high ids (8140-8155ish) 
sidp <- read.spss("data/selfreports/INTAKE_20170223/SIDP.sav", to.data.frame=TRUE)

#szoid5 and sztypl8 are stored in same column (assuming it's an identical criterion)
#but this blows up the logic of the auto scoring below
sidp <- sidp %>% rename(szoid5=szoid5stypl8) %>% mutate(stypl8=szoid5, UsrID=as.numeric(paste0(PTNUM, DyadID))) %>% 
    select(-initials, -sectionQpatient, -sectionQinterviewer)

sidp <- filter(sidp, mth==0 & raterID==0 & ratercode==6) #raterID 0 is case conference; ratercode 6 is case conference; mth 0 is intake
score_pd <- function(df, prefix, dropitems=TRUE) {
  require(dplyr)
  for (p in prefix) {
    tosum <- df %>% select(matches(paste0(p, "\\d+")))
    df[[paste0(p, "_sidp")]] <- apply(tosum, 1, function(r) { sum(r)})
    df <- df %>% select(-matches(paste0(p, "\\d+")))
  }
  return(df)
}
#sort(names(sidp))
sidp <- score_pd(sidp, c("antso", "avoid", "bordl", "depen", "deprs", "histr", "narci", "negtv", "obcmp", "parnd", "stypl", "szoid"))
#sidp %>% filter(mth==0) %>% select(starts_with("histr"))




#counts are present/absent ratings, whereas dimensional scores account for 0, 1, 2 coding
#use dimensional (_sidp) for now.


couples_baseline_clin = Reduce(function(...) { merge(..., by=c("UsrID", "PTNUM", "DyadID"), all=TRUE) },
    list(biq, cts, das, ecr, pai, iip, sidp))

#rename age and sex variables to be clear that these pertain to participant, not alter
couples_baseline_clin <- plyr::rename(couples_baseline_clin, c(age="p_age", sex="p_sex"))

#recode 9999 as missing

#outlier check
#lapply(couples_baseline_clin, function(x) { if (is.numeric(x)) { range(x, na.rm=TRUE) } } )

couples_baseline_clin <- subset(couples_baseline_clin, PTNUM < 9000) #pilot IDs

#list of rule-out and break-up IDs from Emily Landis
ro_ids <- c(8002, 8012, 8021, 8025, 8037, 8054, 8068, 8076, 8077, 8079, 8089, 8090, 8040, 8047, 8048, 8097, 8119, 8125)
message("Removing rule-out/incomplete IDs: ", paste(ro_ids, collapse=", "))
couples_baseline_clin <- subset(couples_baseline_clin, !PTNUM %in% ro_ids)

#PD_on_IIP is coded 0/1 for some data and no/yes for others
couples_baseline_clin$IIP_havePD <- sapply(couples_baseline_clin$IIP_havePD, function(x) { 
      if (x %in% c("1", "yes")) { 1 } else { 0 }
    })

couples_baseline_clin$p_female <- as.numeric(couples_baseline_clin$p_sex == "female") #define as 0/1
couples_baseline_clin$DASrelationship <- as.numeric(couples_baseline_clin$DASrelationship == "yes") #recode as 0/1

#total of the DSM-IV PDs (excluding depressive and negativistic)
couples_baseline_clin$pdtot <- with(couples_baseline_clin, szoid_sidp + stypl_sidp + parnd_sidp + histr_sidp + bordl_sidp +
        narci_sidp + antso_sidp + obcmp_sidp + depen_sidp + avoid_sidp)

#total PD features excluding BPD (note that this differs from OPD because OPD includes DSM-III PDs)
couples_baseline_clin$nobpd <- with(couples_baseline_clin, szoid_sidp + stypl_sidp + parnd_sidp + histr_sidp +
        narci_sidp + antso_sidp + obcmp_sidp + depen_sidp + avoid_sidp)

#couples_baseline_clin$nobpdCount <- with(couples_baseline_clin, szoidCount + styplCount + parndCount + histrCount +
#        narciCount + antsoCount + obcmpCount + depenCount + avoidCount)

couples_baseline_clin$nonarc <- with(couples_baseline_clin, szoid_sidp + stypl_sidp + parnd_sidp + histr_sidp + bordl_sidp +
        antso_sidp + obcmp_sidp + depen_sidp + avoid_sidp) #non-narc Sx

#couples_baseline_clin$nonarcCount <- with(couples_baseline_clin, szoidCount + styplCount + parndCount + histrCount + bordlCount +
#        antsoCount + obcmpCount + depenCount + avoidCount) #non-narc Sx

#couples_baseline_clin$allpdCount <- with(couples_baseline_clin, szoidCount + styplCount + parndCount + 
#        bordlCount + narciCount + antsoCount + histrCount + avoidCount + depenCount + obcmpCount)

#save a wide version of the clinical data with _0 and _1 variable name suffix denoting the partner and patient, respectively
couples_melt <- melt(couples_baseline_clin[,sapply(couples_baseline_clin, is.numeric)], id.vars=c("PTNUM", "DyadID"))
couples_clin_wide <- dcast(couples_melt, PTNUM ~ variable + DyadID)

#save(couples_baseline_clin, couples_clin_wide, file="data/selfreports/couples_baseline_clinical_9Oct2015.RData")
#write.csv(x=couples_baseline_clin, file="data/selfreports/couples_baseline_clinical_9Oct2015.csv", row.names=FALSE)

save(couples_baseline_clin, couples_clin_wide, file="data/selfreports/couples_baseline_clinical_27Jul2017.RData")
write.csv(x=couples_baseline_clin, file="data/selfreports/couples_baseline_clinical_27Jul2017.csv", row.names=FALSE)


#sink("private/clinical_missing_report_7Oct2015.txt")
missingDataReport(couples_baseline_clin, idVars="UsrID")
#sink()
