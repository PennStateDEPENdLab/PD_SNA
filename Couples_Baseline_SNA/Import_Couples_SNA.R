library(lme4)
library(igraph)
library(foreign)
library(psych)
library(lavaan)
library(plyr)
setwd("/Users/michael/Tresors/PD_SNA/Couples_Baseline_SNA")
sna <- read.spss("data/SNA alter summary 8000-8094 3.24.15.sav", to.data.frame=TRUE)

#8030 DyadID 0 (partner) missing
sna <- subset(sna, PTNUM != "8030")

#merge together self-reports etc.
biq <- read.spss("data/BIQ 3.31.15.sav", to.data.frame=TRUE)
cts <- read.spss("data/CTS intake scores 3.4.15.sav", to.data.frame=TRUE)
das <- read.spss("data/DAS intake scores 3.4.15.sav", to.data.frame=TRUE)
ecr <- read.spss("data/ECR intake scores 3.4.15.sav", to.data.frame=TRUE)
pai <- read.spss("data/PAI intake scores 3.4.15.sav", to.data.frame=TRUE)
#iip <- read.spss("data/IIP intake scores.sav", to.data.frame=TRUE)
iip <- read.spss("data/IIP90 intake.sav", to.data.frame=TRUE)
iip$iip_agency <- with(iip, .25*(iip_pa - iip_hi + .707*(iip_bc + iip_no - iip_fg - iip_jk)))
iip$iip_communion <- with(iip, .25*(iip_lm - iip_de + .707*(iip_no + iip_jk - iip_bc - iip_fg)))
iip <- subset(iip, select=c(UsrID, PTNUM, DyadID, IIP_PD1, IIP_PD2, IIP_PD3, iip_agency, iip_communion))

sidp <- read.spss("data/SIDP4 intake consensus scored 4.2.15.sav", to.data.frame=TRUE)
sidp <- sidp[, c("UsrID", "PTNUM", "DyadID", grep(".*_sidp$", names(sidp), value=TRUE), grep(".*Count$", names(sidp), value=TRUE))]

#counts are present/absent ratings, whereas dimensional scores account for 0, 1, 2 coding
#use dimensional (_sidp) for now.

#has 8888 as missing code
lapply(das, function(x) { if (is.numeric(x)) { range(x, na.rm=TRUE) } } )

recodeMissing <- function(vec, targetval=99) {
  vec[vec==targetval] <- NA
  return(vec)
}

das <- colwise(recodeMissing, targetval=8888)(das)
            
other = Reduce(function(...) { merge(..., by=c("UsrID", "PTNUM", "DyadID"), all=TRUE) },
    list(biq, cts, das, ecr, pai, iip, sidp))

#recode 9999 as missing
table(sna$prxsk)
table(sna$prxsk_num)

sna$prxsk_num <- as.numeric(sna$prxsk) - 1 #subtract 1 to keep 0-6 scaling from printed form
sna$sepds_num <- as.numeric(sna$sepds) - 1 
sna$sfhvn_num <- as.numeric(sna$sfhvn) - 1
sna$secbs_num <- as.numeric(sna$secbs) - 1
sna$Angry_num <- 6 - as.numeric(sna$Angry) #Reverse code anger so that higher scores indicate more anger
sna$Happy_num <- as.numeric(sna$Happy) #Happy appears to be coded in the right direction, corr with attach_total is .53 
sna$UsrID <- factor(sna$UsrID)
sna$PTNUM <- factor(sna$PTNUM)

#need to code egocentric ties numerically and align with alter tie csv files
#code such that larger weights indicate greater closeness
#need to recode as 2-5, since 1 (have never met) only occurs in alter tie csv files
#thus, reverse code such that 5=Very close, 4=Moderately close, 3=Slightly close, 2=Not at all close
sna$Close_num <- as.numeric(sna$Close) + 1

#trim white space from alter
sna$Alter <- as.character(gsub("(^\\s*|\\s*$)", "", sna$Alter, perl=TRUE))

round(prop.table(table(sna$prxsk_num)), 3)
round(prop.table(table(sna$sepds_num)), 3)
round(prop.table(table(sna$sfhvn_num)), 3)
round(prop.table(table(sna$secbs_num)), 3)

#define attachment figure as 5 (agree somewhat) or 5 (agree strongly) on each of these 4 items
#note that as.numeric above 
sna$attach_figure <- sna$prxsk_num >= 5 & sna$sepds_num >= 5 & sna$sfhvn_num >= 5 & sna$secbs_num >= 5
#sna$attach_figure <- sna$prxsk_num >= 6 & sna$sepds_num >= 6 & sna$sfhvn_num >= 6 & sna$secbs_num >= 6

#this results in about a 10% rate of attachment figures, which is reasonable
table(sna$attach_figure)
prop.table(table(sna$attach_figure))

table(other$PAIBORtot >= 30)
#hist(other$PAIBORtot)

#good reliability for attachment items (0.89)
alpha(sna[,c("prxsk_num", "sepds_num", "sfhvn_num", "secbs_num")])

sna$attach_total <- apply(sna[,c("prxsk_num", "sepds_num", "sfhvn_num", "secbs_num")], 1, sum)

#aggregate alter statistics for each participant
sna_agg <- ddply(sna, .(UsrID), function(subdf) {
      past_rom <- sum(subdf$Romance=="yes, in the past")
      cur_rom <- sum(subdf$Romance=="yes, currently")
      all_rom <- sum(subdf$Romance %in% c("yes, currently", "yes, in the past"))
      tot_cutoff <- sum(subdf$Cutoff=="yes")
      avg_duration <- mean(subdf$YrsKnwn, na.rm=TRUE)
      attach_total <- sum(subdf$attach_figure, na.rm=FALSE)
      data.frame(UsrID=subdf$UsrID[1], PTNUM=subdf$PTNUM[1], Dyad=subdf$Dyad[1],
          past_rom, cur_rom, all_rom, tot_cutoff, avg_duration, attach_total)
    })

sna_merge <- merge(sna, other[,c("UsrID", "sex", "age", "ECRanx", "ECRavoid", "PAIBORaffe", "PAIBORiden", "PAIBORnegr", "PAIBORtot", "IIP_PD1", "IIP_PD2", "IIP_PD3", "iip_agency", "iip_communion")], all.x=TRUE, by="UsrID")
setdiff(unique(sna$PTNUM), unique(other$PTNUM)) #has SNA but no other
setdiff(unique(other$PTNUM), unique(sna$PTNUM)) #has other but no SNA

#add row number to be able to track down alter attributes later (when populating vertices)
sna_merge$rownum <- 1:nrow(sna_merge)

options(error=recover)
#try to line up couples alters

##STEP 1: TRY TO IDENTIFY ALTERS IN COMMON ON THE BASIS OF NAMES AND AGES
couple_split <- split(sna_merge, sna_merge$PTNUM)
couple_match <- lapply(couple_split, function(subdf) {
      
      patient_data <- subset(subdf, Dyad=="patient")
      partner_data <- subset(subdf, Dyad=="partner")
      
      #use adist (string distance) to look for matches on alter
      distmat <- matrix(NA_real_, nrow=length(patient_data$Alter), ncol=length(partner_data$Alter))
      for (i in 1:length(patient_data$Alter)) {
        distmat[i,] <- sapply(partner_data$Alter, function(partner_name) {
              adist(patient_data$Alter[i], partner_name, ignore.case=TRUE, partial=FALSE)
            })        
      }
      
      patient_bestmatch <- apply(distmat, 1, which.min) #for each patient alter (1-25), what is the position of the nearest partner alter (1-25)
      patient_disttomatch <- apply(distmat, 1, min) #how close is the potential match
      
      #don't really need to go in both directions since it is symmetric
      #partner_match <- apply(distmat, 2, which.min)
      #partner_disttomatch <- apply(distmat, 2, min)
      
      within_distance <- patient_disttomatch <= 4 #4 distance seems about right to leave out horrible matches
      #use additional criteria below to refine matches (i.e., above is overly inclusive)
      
      #an alter can only be matched once, and duplicate matches are indicative of a false positive.
      #for any duplicates, only retain the alter with the lowest distance.
      #if two (or more) names have the same minimum distance, keep them in (weed out using other metrics below)
      
      #use a logical and of within_distance and is_exclusive_best to generate alter match      
      is_exclusive_best <- sapply(1:length(patient_bestmatch), function(x) {
            dupes <- which(patient_bestmatch == patient_bestmatch[x])
            if (length(dupes) > 1L) {
              if (patient_disttomatch[x] == min(patient_disttomatch[dupes])) {
                TRUE #best (or equal to another best)
              } else { FALSE } #not among best
            } else { TRUE } #no dupes
          })
      
      #if (subdf$PTNUM[1L] == "8006") { browser() }
      
      alter_match <- which(within_distance & is_exclusive_best)
      
      #also verify that age is within 5 years
      age_diff <- abs(patient_data$Age[alter_match] - partner_data$Age[patient_bestmatch[alter_match]]) <= 5
      #if age is missing, an NA will be generated for age_diff... in this case, default to a match (inclusive)
      age_diff[which(is.na(age_diff))] <- TRUE
      alter_match <- alter_match[age_diff]
      
      #filter on last initial match
      last_initial_patient <- tolower(sub("^\\s*\\w+.*[\\s_]+([A-z]+)[_\\s]*$", "\\1", patient_data$Alter[alter_match], perl=TRUE))
      last_initial_partner <- tolower(sub("^\\s*\\w+.*[\\s_]+([A-z]+)[_\\s]*$", "\\1", partner_data$Alter[patient_bestmatch[alter_match]], perl=TRUE))
      
      #check to make sure this works
      #if(any(grepl("_", patient_data$Alter[name_match], fixed=TRUE))) { browser()}
      
      #interestingly, some of the later couples do not have spaces to separate last names, only underscores...
      #even more confusing: some folks end the name with an underscore (e.g., Merrit_Gi_). Filter out
      
      if (length(last_initial_patient) > 0L) {
        stopifnot(length(last_initial_patient) == length(last_initial_partner))
        last_initial_match <- sapply(1:length(last_initial_patient), function(x) {
              #truncate to same length (so, use two-character last names if available, but fall back to one as needed)
              minlen <- min(c(nchar(last_initial_patient[x]), nchar(last_initial_partner[x])))
              substr(last_initial_patient[x], 1, minlen) == substr(last_initial_partner[x], 1, minlen)
            })
        
        alter_match <- alter_match[last_initial_match]
      }
      
      #filter on whether first letter of first name matches (could lead to problems with nick names?)
      first_initial_patient <- tolower(sub("^\\s*([A-z])\\w+.*[\\s_]+[A-z]+[_\\s]*$", "\\1", patient_data$Alter[alter_match], perl=TRUE))
      first_initial_partner <- tolower(sub("^\\s*([A-z])\\w+.*[\\s_]+[A-z]+[_\\s]*$", "\\1", partner_data$Alter[patient_bestmatch[alter_match]], perl=TRUE))
      
      if (length(first_initial_patient) > 0L) {
        stopifnot(length(first_initial_patient) == length(first_initial_partner))
        first_initial_match <- first_initial_patient == first_initial_partner
        
        alter_match <- alter_match[first_initial_match]
      }
      
      if (length(alter_match) > 0L) {
        common <- data.frame(patient_data[alter_match,c("Alter", "Age", "rownum")], partner_data[patient_bestmatch[alter_match], c("Alter", "Age", "rownum")])
        names(common) <- c("patient_name", "patient_age", "patient_rownum", "partner_name", "partner_age", "partner_rownum")
        patient_only <- patient_data[-1*alter_match,c("Alter", "Age", "rownum")]
        partner_only <- partner_data[-1*patient_bestmatch[alter_match], c("Alter", "Age", "rownum")]
      } else {
        common <- list()
        patient_only <- patient_data[,c("Alter", "Age", "rownum")]
        partner_only <- partner_data[,c("Alter", "Age", "rownum")]
      }
      
      idlist <- rbind(
          common,
          data.frame(patient_name=patient_only$Alter, patient_age=patient_only$Age, patient_rownum=patient_only$rownum,
              partner_name=NA_character_, partner_age=NA_real_, partner_rownum=NA_integer_),
          data.frame(patient_name=NA_character_, patient_age=NA_real_, patient_rownum=NA_integer_, 
              partner_name=partner_only$Alter, partner_age=partner_only$Age, partner_rownum=partner_only$rownum)
      )
      
      idlist$alter_id <- 1:nrow(idlist)
      #generate shared name for adjacency list
      idlist$shared_name <- apply(idlist[,c("patient_name", "partner_name")], 1, function(r) { 
            if (is.na(r["patient_name"])) {
              r["partner_name"]
            } else if (is.na(r["partner_name"])) {
              r["patient_name"]
            } else {
              r["patient_name"] #if both are present (common), use the patient name
            }
          })
      
      #there is a chance of an identical name, but different alters (e.g., shawn d for 8000 -- probably a jr)
      #to avoid confusion downstream, make shared name name.alter_id
      idlist$shared_name <- paste(idlist$shared_name, idlist$alter_id, sep=".")
      
      return(list(common=common, patient_only=patient_only, partner_only=partner_only, idlist=idlist))
    })

#8030 only has 25 alters

#output file of matches for manual checks
sink(file="couple_altermatch_31Mar2015.txt", append=FALSE)
print(couple_match)
sink()

#manual fixes for name matching
#8006:
#patient name: jess ra
#partner name: Jessica Rad
#Jean Ray is a closer string match to Jess Ra, so leading to false negative
#Changed patient jess ra to Jessica Rad (row 292)

#8010:
#patient name: J R L
#partner name: JR Le
#not matching using heuristics above. Manually fixed dataset (row 477) to be JR Le for patient

#8023:
#patient name: jennifer p
#partner name: Jen P
#string distance a bit too far due to nick name. Manually changed dataset (row 1018) to jennifer p for partner

#8031:
# patient name: Mary Jo Da (age 50)
# partner name: Mary Jo Mc (age 50)
# This is the same person according to RA checks. Change patient name to Mary Jo Mc (row 1327)

#8032:
# patient name: Celeste Pe (age 31)
# partner name: Celeste Pi (age 31)
# verified by RAs to be the same person. Change patient last name to Pi (row 1403)

#8041:
# patient name: Kaitlin Ge (age NA)
# partner name: Caitlyn Ge (age 30)
# change patient name to Caitlyn Ge (age 30) (row 1670)
# patient name: BJ Ge (age 32)
# partner name: BJ III Ge (age 30)
# change patient name to BJ III Ge (row 1660)

#8043:
# patient name: Edey Se (age NA)
# partner name: Edith Se (age 51)
# hard to tell based on this info, but Edey seems like the obvious nickname for edith
# HOLDING OFF FOR NOW

#8046:
# patient name: Josefina As (age 85, family)
# partner name: Josephina Pa (age 93, partner's family)
# very likely match based on relationship and uncommon name. Change partner name to Josefina As and age to 85 (row 1832).

#8058: 
#has a problem with jr and sr. 
#partner names: "Carlos_Cr" and "Carlos_Jr_Cr"
#patient names: "Carlos_Sr_Cr" and "Carlos_Jr_Cr"
#leads to Jr matching better than Sr. Manually fixed in dataset (row 2185) to be Carlos_Sr_Cr in both cases

#8061:
# patient name: Lo_WI (age 35)
# partner name: Lo_Wi (age 29)
# change both ages to 32 (midpoint)

#8063:
# patient name: Mary_Ann_Va (age 70)
# partner name: Mary-Anne_VB (age 70)
# change partner name to Mary_Ann_Va

# patient name: Carson_Va (age 8)
# partner name: Carson_VB (age 7)
# change partner name to Carson_Va (row 2443)

# patient name: Jeffrey_Va (age 31)
# partner name: Jeff_VB (age 32)
# change partner name to Jeffrey_Va (row 2438)

# patient name: Dan_Go (age 58)
# partner name: Daniel_Go (age 58)
# change patient name to Daniel_Go (row 2460)

# patient name: Landen_Va (age 1)
# partner name: Landen_VB (age 2)
# change partner name to Landen_Va (row 2444)

# patient name: Tim_Va (age 33)
# partner name: Timmy_VB (age 33)
# change partner name to Tim_Va

# patient name: Kaeden_Va (age 10)
# partner name: Kaeden_VB (age 10)
# change partner name to Kaeden_Va (row 2442)

#8064:
# patient name: Gerry_Ch (age 67)
# partner name: Geraldine_Ch (age 63)
# change patient name to Geraldine_Ch (row 2502)

#8066:
# patient name: Ruth_Ra (age 71)
# partner name: Ruth_Ra (age 78)
# change patient and partner ages to 74.5 (midpoint)

#8067:
# patient name: Don_BA (age 61)
# partner name: Donald_Ba (age 60)
# change patient name to Donald_Ba (row 2652)

# patient name: Dar_BA (age 64)
# partner name: Darlene_Ba (age 55)
# change patient name to Darlene_Ba and both ages to 60

# patient name: Tony_PI (age 50)
# partner name: Anthony_Pi (age 50)
# change patient name to Anthony_Pi (row 2657)

#8069:
# patient name: Debbie_Fo (age 49)
# partner name: Debbie_Fo (age 57)
# change both ages to 53 (midpoint)

# patient name: Ken_Fo (age 47)
# partner name: Ken_Fo (age 60)
# change both ages to 53

# patient name: Zack_Co (age 36)
# partner name: Zack_Go (age 38)
# RA notes indicate relationship match. Change partner name to Zack_Co (row 2686)

#8070:
# patient name: Mark_Go (age 24)
# partner name: Mark_Ge (age 22)
# RAs indicate relationship match. Change partner name to Mark_Go 

#8072:
# patient name: Kathy_Ta (age 57)
# partner name: Kathryn_Ta (age 49)
# RAs indicate relationship match. Change patient name to Kathryn_Ta and both ages to 53.

# patient name: Sean_Ke (age 3)
# partner name: Shawn_De (age 3)
# change partner name to Sean_Ke (row 2845)

# patient name: Carson_Ke (age 2)
# partner name: Carson_De (age 2)
# change partner name to Carson_Ke.

#8074:
# patient name: Fazilla_Sh (age 28)
# partner name: Fazilla_Ja (age 28)
# change partner name to Fazilla_Sh (row 2888)

#8082:
# patient name: Stephanie_Ge (age 56)
# partner name: Stephanie_Ge (age 50)
# change ages to 53

# patient name: Mark_Jr_Mi (age 37)
# partner name: Marky_Mi (age 37)
# change partner name to Mark_Jr_Mi

# patient name: Mark_Sr_Mi (age 60)
# partner name: Mark_Mi (age 60)
# change partner name to Mark_Sr_Mi

#8084:
# patient name: Deb_Ca (age 58)
# partner name: Deb_Ca (age 51)
# change both ages to 54

#8088:
# patient name: Annette_Pa (age 43)
# partner name: Annette_Pa (age 49)
# change both ages to 46

#8091:
# patient name: Mike_Re (age 55)
# partner name: Michael_Re (age 54)
# change patient name to Michael_Re (row 3305)

# patient name: Dave_Al (age 54)
# partner name: David_Al (age 52)
# change patient name to David_Al (3314)



##STEP 2: IMPORT CLOSENESS RATINGS AMONG ALTERS FOR EACH PARTICIPANT TO BUILD EGO/DYAD NETWORK
##  graphlist will contain an element for each participant with actors (vertices) and edges as a list (to match igraph input).
##  Alter names are matched with the couple-matching naming scheme above so that dyadic networks can be built by concatenating
##  the edges for both partners.

#import closeness ratings for all alters to build network
weighted_matrices <- list.files("data/ties", pattern="weighted_matrix.csv", recursive=TRUE, full.names=TRUE)

#filter out follow-up matrices
weighted_matrices <- weighted_matrices[!grepl("6mo", weighted_matrices, ignore.case=TRUE)]

#watch out for > 2 matched files per ID
PTNUMs <- sub("^.*/ties/(\\d+)/.*$", "\\1", weighted_matrices, perl=TRUE)
which(table(PTNUMs) > 2)
which(table(PTNUMs) < 2)

graphlist <- list()
for (f in 1:length(weighted_matrices)) {
  alterties <- read.csv(weighted_matrices[f], header=TRUE)[,-1] #first column contains names of alters (redundant with cols since symmetric)
  #recode 1-5 where higher values indicate greater closeness: 5=very close, 4=moderately close, 3=slightly close, 2=not close at all, 1=have never met
  alterties <- as.data.frame(lapply(alterties, function(x) { 6 - x })) #note: this will produce values of 6 along the diagonal which was previously 0, but we don't use diagonal (self-closeness) anyway
  
  PTNUM <- sub("^.*/ties/(\\d+)/.*$", "\\1", weighted_matrices[f], perl=TRUE)
  alternames <- names(alterties)
  #try to match with either patient or partner from sna alter summary above
  #note: for some early couples (e.g., 8003), the full names (not initials) are present, but these were cleaned for the SPSS file...
  #I could try to build the alter summary myself from the alter_summary csv files, but would lose information about patient partner
  #for now, use summed min adist to match
  if (!is.null(couple_match[[PTNUM]])) { #there are more egonet files in the ties folder than in dataset
    idlist <- couple_match[[PTNUM]][["idlist"]]
    adist_to_patient <- sapply(alternames, function(tiename) {
          #return minimum distance (best match) for this name from the ties file
          dists <- sapply(idlist[,"patient_name"], function(summaryname) {
                adist(summaryname, tiename, ignore.case=TRUE, partial=FALSE) })
          
          c(mindist=min(dists, na.rm=TRUE), whichmin=which.min(dists))
        })
    
    adist_to_partner <- sapply(alternames, function(tiename) {
          #return minimum distance (best match) for this name from the ties file
          dists <- sapply(idlist[,"partner_name"], function(summaryname) {
                adist(summaryname, tiename, ignore.case=TRUE, partial=FALSE) })
          
          c(mindist=min(dists, na.rm=TRUE), whichmin=which.min(dists))
        })
    
    #first row is the minimum distance from each alter in ties csv file to alter in idlist
    #second row is the position within the idlist of the closest match to the alter in the tie csv file 
    #note that the ordering of alters in the SPSS file should roughly correspond to the order of the ties file (but partner is shifted to last position, I think)
    
    #another frustration is that I'll need to calculate the distance between each alter in the ties csv file and the best alter in the idlist
    if (sum(adist_to_partner[1,]) > sum(adist_to_patient[1,])) {
      patient <- 1
      closest_alter <- adist_to_patient
    } else {
      patient <- 0
      closest_alter <- adist_to_partner
    }
    
    #this is the list of shared names (matching alter summary) in the order of the ties matrix
    #use this, not the names in the csv file, to create the adjacency list
    sharednames <- idlist[closest_alter[2,],"shared_name"]
    sharedids <- idlist[closest_alter[2,],"alter_id"]
    
    #each name should only match once -- otherwise, the adist above is failing and we are getting multiple matches for an alter
    #this was happening for some of the same files as the ego alter summaries above.
    #in addition, some of the files seemed like the alters were edited, but then the weighted matrices did not align
    #made manual edits to 8001 anthony and 8023 james to make alter alignment match
    #8026/christina_po_weighted_matrix.csv: had two entries for Liz_S, but alter summary identifies first as Jonathan Sw. Edited weighted matrix to match.
    #8043/MS80431_weighted_matrix.csv: Ruth_Mu was only listed as "_". Manually replaced
    
    stopifnot(length(sharednames) == length(unique(sharednames)))
    
    #TODO: add romance, cutoff, happy, and emotion (categorical) to list of attributes to test for assortativity
    
    #generate vertex list with attach_total, attach_figure, YrsKnwn, Angry_num, partner
    if (patient==1) {
      actors <- data.frame(name=idlist$shared_name[!is.na(idlist$patient_rownum)],
          alter_id=idlist$alter_id[!is.na(idlist$patient_rownum)],
          sna_merge[na.omit(idlist$patient_rownum), c("Alter", "partner", "attach_total", "YrsKnwn", "Angry_num", "Happy_num", "Close_num", "Romance")], stringsAsFactors=FALSE) #put in original name for clarity
    } else {
      actors <- data.frame(name=idlist$shared_name[!is.na(idlist$partner_rownum)],
          alter_id=idlist$alter_id[!is.na(idlist$partner_rownum)],
          sna_merge[na.omit(idlist$partner_rownum), c("Alter", "partner", "attach_total", "YrsKnwn", "Angry_num", "Happy_num", "Close_num", "Romance")], stringsAsFactors=FALSE) #put in original name for clarity
    }
    
    #go over matrix and generate adjacency list by mapping alter names to the idlist and using common alter_ids
    #25 choose 2 should give a total of 300 edges.
    elist <- c()
    for (i in 1:nrow(alterties)) {
      for (j in 1:ncol(alterties)) {
        if (i > j) {
          #only iterate over lower triangle of matrix
          #name the closeness rating "weight" to be recognized by igraph as a weighted edgelist
          elist <- rbind(elist, 
              data.frame(a1=sharednames[i], a2=sharednames[j], weight=alterties[i,j], stringsAsFactors=FALSE)
          )
        }
      }
    }
    
    #also need to add connections for ego
    #looks like these were left out from egonet calculations of closeness, degree, etc.
    egocon <- adply(actors, 1, function(rdf) {
          data.frame(a1="SELF", a2=rdf$name, weight=rdf$Close_num)
        })[,c("a1", "a2", "weight")] #only keep edge fields to match above (plyr returns original fields, too)
    elist_withego <- rbind(elist, egocon)
    
    actors_withego  <- rbind.fill(actors, data.frame(name="SELF", alter_id=0)) #have to add self to vertices
    
    graphlist[[paste(PTNUM, patient, sep=".")]] <- list(PTNUM=as.integer(PTNUM), UsrID=as.integer(paste0(PTNUM, patient)), actors=actors, edges=elist, actors_withego=actors_withego, edges_withego=elist_withego)   
    
  } else {
    warning("No summary data available for file: ", weighted_matrices[f])
  }
}


##STEP 3: BUILD DYADIC NETWORKS BY MERGING ACTORS FOR EACH PARTNER AND COMBINING EDGES
## Need to be careful about handling edges in common across partners so that two edges don't exist for one possible alter-alter connection
library(parallel)
cl_fork <- makeForkCluster(nnodes=8)

##now match up couples' networks by splitting graphlist names and identifying couples
idsplit <- do.call(rbind, strsplit(names(graphlist), ".", fixed=TRUE))

#verify that we have two networks per couple
table(idsplit[,1])

couple_graphs <- list()

for (id in unique(idsplit[,1])) {
  matches <- which(idsplit[,1] == id)
  stopifnot(length(matches) == 2)
  
  #to blend we need to identify the partner within the ego network and add in ego-centered ratings to the edge list
  #also need to extend attributes of actors to differentiate patient versus partner ratings of closeness, years known, etc.
  
  patient <- graphlist[[ matches[which(idsplit[matches,2] == "1")] ]]
  partner <- graphlist[[ matches[which(idsplit[matches,2] == "0")] ]]
  
  patient_actors <- patient$actors
  partner_actors <- partner$actors
  
  notid <- ! names(patient_actors) %in% c("name", "alter_id")
  names(patient_actors)[notid] <- paste("patient", names(patient_actors)[notid], sep=".")
  
  notid <- ! names(partner_actors) %in% c("name", "alter_id")
  names(partner_actors)[notid] <- paste("partner", names(partner_actors)[notid], sep=".")
  
  combined_actors <- merge(patient_actors, partner_actors, by=c("name", "alter_id"), all=TRUE)
  combined_actors <- combined_actors[order(combined_actors$alter_id),] #order by alter id for clarity
  combined_actors$role <- factor(apply(combined_actors[,c("patient.partner", "partner.partner", "patient.Alter", "partner.Alter")], 1, function(row) {
            if (!is.na(row["partner.partner"]) && row["partner.partner"] == "yes") {
              "patient" #the partner's partner is the patient
            } else if (!is.na(row["patient.partner"]) && row["patient.partner"] == "yes") {
              "partner" #the patient's partner is the partner
            } else if (is.na(row["patient.Alter"])) {
              "partner_only" #nominated only by partner
            } else if (is.na(row["partner.Alter"])) {
              "patient_only"
            } else if (! (is.na(row["parter.Alter"]) && is.na(row["patient.Alter"]))) {
              "shared"
            } else {
              stop("Unable to match")
            }
          }))
  
  #drop patient.partner, partner.partner, patient.Alter, and partner.Alter since these are replaced by role
  combined_actors <- combined_actors[, !names(combined_actors) %in% c("patient.partner", "partner.partner", "patient.Alter", "partner.Alter")]
  
  patient_altername <- combined_actors$name[which(combined_actors$role == "patient")]
  partner_altername <- combined_actors$name[which(combined_actors$role == "partner")]
  
  if (length(patient_altername) > 1) stop("More than one patient name matched")
  if (length(partner_altername) > 1) stop("More than one partner name matched")
  
  #add patient closeness ratings to edge list
  patient_egoedges <- adply(patient_actors, 1, function(rdf) {
        data.frame(a1=patient_altername, a2=rdf$name, weight=rdf$patient.Close_num)
      })[,c("a1", "a2", "weight")] #only keep edge fields to match above (plyr returns original fields, too)
  
  partner_egoedges <- adply(partner_actors, 1, function(rdf) {
        data.frame(a1=partner_altername, a2=rdf$name, weight=rdf$partner.Close_num)
      })[,c("a1", "a2", "weight")] #only keep edge fields to match above (plyr returns original fields, too)
  
  #screen for 1) self connections (invalid) and 2) overlapping ratings of alters
  
  #check for matching alter id in patient list based on partner actors
  if ((patient_name <- as.character(with(partner$actors, name[partner=="yes"]))) %in% patient$actors$name) {
    warning("For id: ", id, ", patient has listed self among alters. Removing this self-connection.")
    patient$actors <- subset(patient$actors, name != patient_name)
    patient$edges <- subset(patient$edges, a1 != patient_name & a2 != patient_name)
  }
  
  if ((partner_name <- as.character(with(patient$actors, name[partner=="yes"]))) %in% partner$actors$name) {
    warning("For id: ", id, ", partner has listed self among alters. Removing this self-connection.")
    partner$actors <- subset(partner$actors, name != partner_name)
    #only keep the ego edges, not (redundant) closeness ratings from the weighted matrix.
    partner$edges <- subset(partner$edges, a1 != partner_name & a2 != partner_name)    
  }  
  
  patient_edges <- rbind(patient$edges, patient_egoedges)
  partner_edges <- rbind(partner$edges, partner_egoedges)
  
  clusterExport(cl_fork, c("patient_edges", "partner_edges")) #export edge lists to workers (otherwise may not get correct lists)
  
  #look for edges in common and create variants of edge list that prefer 1) partner, 2) patient, or 3) mean closeness
  #this is slow! (and a lousy way to code it, but don't care at the moment. use parallel to speed up a bit)
  patient_partneredges_match <- parSapply(cl_fork, 1:nrow(patient_edges), function(ri) { #sapply(1:nrow(patient_edges), function(ri) { #
        rpatient <- patient_edges[ri,c("a1", "a2")] 
        
        partner_pos <- sapply(1:nrow(partner_edges), function(rj) {
              rpartner <- partner_edges[rj,c("a1", "a2")]
              
              if ((rpatient["a1"] == rpartner["a1"] && rpatient["a2"] == rpartner["a2"]) || 
                  (rpatient["a1"] == rpartner["a2"] && rpatient["a2"] == rpartner["a1"])) {
                print(rpatient)
                print(rpartner)
                return(c(ri, rj))
              } else {
                return(NA)
              }
            })
        
        if (sum(is.na(partner_pos)) < length(partner_pos)) {
          #match found
          if (length(partner_pos) - sum(is.na(partner_pos)) != 1) { warning("Duplicate ties."); browser() }
          return(na.omit(do.call(c, partner_pos)))
        } else {
          return(NULL)
        }
      })
  
  #make into match x c(patient_row, partner_row) matrix
  patient_partneredges_match <- do.call(rbind, patient_partneredges_match)
  
  #if there is at least one overlapping edge, we have work to do
  if (!is.null(patient_partneredges_match[1L])) {
    patient_edges_removeshared <- patient_edges[-1*patient_partneredges_match[,1],]
    partner_edges_removeshared <- partner_edges[-1*patient_partneredges_match[,2],]
    
    patient_edges_shared <- patient_edges[patient_partneredges_match[,1],]
    partner_edges_shared <- partner_edges[patient_partneredges_match[,2],]
    
    shared_weights <- cbind(patient_edges_shared$weight, partner_edges_shared$weight) 
    
    shared_edge_average <- data.frame(patient_edges_shared[,c("a1", "a2")], weight=apply(shared_weights, 1, mean, na.rm=TRUE))
    
    combined_edges_prefpatient <- rbind(patient_edges, partner_edges_removeshared)
    combined_edges_prefpartner <- rbind(partner_edges, patient_edges_removeshared)
    combined_edges_avgratings <- rbind(patient_edges_removeshared, partner_edges_removeshared, shared_edge_average)
    
  } else {
    combined_edges_prefpatient <- combined_edges_prefpartner <- combined_edges_avgratings <- rbind(patient_edges, partner_edges) #no overlap
  }
  
  #deprecate straight combination since this leads to multiple edges for a given connnection
  #combined_edges <- rbind(patient_edges, partner_edges)
  
  couple_graphs[[id]] <- list(PTNUM=as.integer(id),
      actors=combined_actors, 
      edges_prefpatient=combined_edges_prefpatient,
      edges_prefpartner=combined_edges_prefpartner,
      edges_avgratings=combined_edges_avgratings,
      patient_partnerclose=with(patient$actors, Close_num[partner=="yes"]), 
      partner_patientclose=with(partner$actors, Close_num[partner=="yes"])
  )
  
}

stopCluster(cl_fork)

save(couple_graphs, graphlist, couple_match, sna, sna_agg, sna_merge, other, file="data/SNA_Processed_6Apr2015.RData")

##GENERATE EGO-CENTERED NETWORK METRICS
vertex_agg <- c()
graph_agg <- c()

##EGO-CENTERED ANALYSES
for (p in graphlist) {
  #consider whether to include ego-based edges since we may not consider an alter truly connected to ego if the closeness rating is weak 
  #g <- graph.data.frame(p$edges_withego, directed=FALSE, vertices=p$actors_withego)
  g <- graph.data.frame(p$edges, directed=FALSE, vertices=p$actors)
  #print(g, e=TRUE, v=TRUE)
  
  #remove any connections that are not close at all
  gthresh2 <- delete.edges(g, which(E(g)$weight <= 2)) #1=have never met, 2=not close at all
  gthresh3 <- delete.edges(g, which(E(g)$weight <= 3)) #1=have never met, 2=not close at all, 3=slightly close
  gbin2 <- remove.edge.attribute(gthresh2, "weight") #unweighted graph <= 2
  gbin3 <- remove.edge.attribute(gthresh3, "weight") #unweighted graph <= 3
  
  #node/vertex metrics
  
  #degree centrality (degree does not account for weight)
  degree_2 <- degree(gthresh2) #element for each vertex
  centralization_degree_2 <- centralization.degree(gthresh2)$centralization #graph level centrality index according to degree
  
  degree_3 <- degree(gthresh3) #element for each vertex
  centralization_degree_3 <- centralization.degree(gthresh3)$centralization #graph level centrality index
  
  #weighted equivalent is called strength and represents sum of weights for connected vertices
  #should probably re-weight the closeness ratings to be 1, 2, 3 or something
  strength_2 <- graph.strength(gthresh2)
  strength_3 <- graph.strength(gthresh3)
  
  #eigenvector centrality
  evcent_weighted_2 <- evcent(gthresh2)$vector #element for each vertex
  evcent_binary_2 <- evcent(gbin2)$vector #element for each vertex
  centralization_eigenvector_2 <- centralization.evcent(gbin2)$centralization
  
  #eigenvector centrality
  evcent_weighted_3 <- evcent(gthresh3)$vector #element for each vertex
  evcent_binary_3 <- evcent(gbin3)$vector #element for each vertex
  centralization_eigenvector_3 <- centralization.evcent(gbin3)$centralization
  
  #closeness centrality
  closeness_weighted_2 <- closeness(gthresh2, normalized=TRUE)*100 #normalized and multiplied by 100 to match ucinet.
  closeness_binary_2 <- closeness(gbin2, normalized=TRUE)*100 #normalized and multiplied by 100 to match ucinet.
  
  closeness_weighted_3 <- closeness(gthresh3, normalized=TRUE)*100 #normalized and multiplied by 100 to match ucinet.
  closeness_binary_3 <- closeness(gbin3, normalized=TRUE)*100 #normalized and multiplied by 100 to match ucinet.
  
  betweenness_weighted_2 <- betweenness(gthresh2, normalize=TRUE)*100
  betweenness_binary_2 <- betweenness(gbin2, normalize=TRUE)*100
  
  betweenness_weighted_3 <- betweenness(gthresh3, normalize=TRUE)*100
  betweenness_binary_3 <- betweenness(gbin3, normalize=TRUE)*100
  
  #do alters connect with each other because of similarities on vertex attributes? (esp. attachment
  attach_assortativity <- assortativity(gthresh2, V(g)$attach_total)
  angry_assortativity <- assortativity(gthresh2, V(g)$Angry_num)
  yearsknown_assortativity <- assortativity(gthresh2, V(g)$YrsKnwn)
  
  #TODO: add romance, cutoff, happy, and emotion (categorical) to list of attributes to test for assortativity
  #also consider whether to test each attachment item separately here.
  
  #something to consider adding: number of cohesive blocks
  #b <- cohesive.blocks(gthresh3)
  #plot(b, g)
  
  #cliques (which are based on unweighted edges)
  largest_cliques_2 <- largest.cliques(gbin2)
  maximal_3cliques_2 <- maximal.cliques(gbin2, min=3) #maximal cliques of size 3 or larger
  
  largest_cliques_3 <- largest.cliques(gbin3)
  maximal_3cliques_3 <- maximal.cliques(gbin3, min=3) #maximal cliques of size 3 or larger
  
  #community detection using modularity optimization algorithm
  optcommunity_weighted_2 <- optimal.community(gthresh2)
  optcommunity_binary_2 <- optimal.community(gbin2)
  optcommunity_weighted_3 <- optimal.community(gthresh3)
  optcommunity_binary_3 <- optimal.community(gbin3)
  
  #community detection using edge betweenness algorithm
  #omitting for now
#  ebc_binary_2 <- edge.betweenness.community(gbin2)
#  ebc_weighted_2 <- edge.betweenness.community(gthresh2)
#  ebc_binary_3 <- edge.betweenness.community(gbin3)
#  ebc_weighted_3 <- edge.betweenness.community(gthresh3)
  
  #comparing parcellations
  #compare(ebc_binary_2, ebc_binary_3, method="nmi")
  
  #eccentricity (shortest path from each vertex to furthest point away in graph)
  #part of small worldness -- only sensible for binary graphs
  #note: eccentricity is zero for a vertex that is unconnected to any other
  #when the ego edges are left out, this would mean an alter connected only to self would have eccentricity = 0
  eccentricity_2 <- eccentricity(gbin2)
  eccentricity_3 <- eccentricity(gbin3)
  
  #density of graphs (number of edges versus number of possible edges
  density_2 <- graph.density(gbin2, loops=FALSE)
  density_3 <- graph.density(gbin3, loops=FALSE)
  
  #local clustering
  locclust_binary_2 <- transitivity(gbin2, type="local") #will be NaN for a vertex with no connected triples
  locclust_weighted_2 <- transitivity(gthresh2, type="weighted") #will be NaN for a vertex with no connected triples
  locclust_binary_3 <- transitivity(gbin3, type="local") #will be NaN for a vertex with no connected triples
  locclust_weighted_3 <- transitivity(gthresh3, type="weighted") #will be NaN for a vertex with no connected triples
  
  #overall transitivity (global clustering) of graph
  transitivity_2 <- transitivity(gbin2, type="global")
  transitivity_3 <- transitivity(gbin3, type="global")
  
  #average path length of graph
  avgpath_2 <- average.path.length(gbin2, unconnected=TRUE)
  avgpath_3 <- average.path.length(gbin3, unconnected=TRUE)
  
  #figure out small worldness of network according to Humphries 2008
  #generate an E-R random graph with the same number of nodes and edges
  
  g_rand_2 <- erdos.renyi.game(n=vcount(gbin2), p.or.m=ecount(gbin2), type=c("gnm"),directed = FALSE, loops = FALSE)
  transitivity_rand_2 <- transitivity(g_rand_2, type="global")
  avgpath_rand_2 <- average.path.length(g_rand_2, unconnected=TRUE)
  
  g_rand_3 <- erdos.renyi.game(n=vcount(gbin3), p.or.m=ecount(gbin3), type=c("gnm"),directed = FALSE, loops = FALSE)
  transitivity_rand_3 <- transitivity(g_rand_3, type="global")
  avgpath_rand_3 <- average.path.length(g_rand_3, unconnected=TRUE)
  
  #from humphries 2008
  smallworld_2 <- (transitivity_2/transitivity_rand_2)/(avgpath_2/avgpath_rand_2)
  smallworld_3 <- (transitivity_3/transitivity_rand_3)/(avgpath_3/avgpath_rand_3)
  
  #assemble vertex measures and graph measures into distinct data.frames
  
  #add averaged vertex measures to graph measures for descriptives
  
  vertex_measures <- data.frame(
      PTNUM=p$PTNUM,
      UsrID=p$UsrID,
      alter_name=V(g)$name,
      degree_2=degree_2,
      degree_3=degree_3,
      strength_2=strength_2,
      strength_3=strength_3,
      evcent_weighted_2=evcent_weighted_2,
      evcent_binary_2=evcent_binary_2,
      evcent_weighted_3=evcent_weighted_3,
      evcent_binary_3=evcent_binary_3,
      closeness_weighted_2=closeness_weighted_2,
      closeness_binary_2=closeness_binary_2,
      closeness_weighted_3=closeness_weighted_3,
      closeness_binary_3=closeness_binary_3,
      betweenness_weighted_2=betweenness_weighted_2,
      betweenness_binary_2=betweenness_binary_2,
      betweenness_weighted_3=betweenness_weighted_3,
      betweenness_binary_3=betweenness_binary_3,
      eccentricity_2=eccentricity_2,
      eccentricity_3=eccentricity_3,
      locclust_binary_2=locclust_binary_2,
      locclust_weighted_2=locclust_weighted_2,
      locclust_binary_3=locclust_binary_3,
      locclust_weighted_3=locclust_weighted_3      
  )
  
  graph_measures <- data.frame(
      PTNUM=p$PTNUM,
      UsrID=p$UsrID,
      centralization_degree_2=centralization_degree_2,
      centralization_degree_3=centralization_degree_3,
      centralization_eigenvector_2=centralization_eigenvector_2,
      centralization_eigenvector_3=centralization_eigenvector_3,
      attach_assortativity=attach_assortativity,
      angry_assortativity=angry_assortativity,
      yearsknown_assortativity=yearsknown_assortativity,
      largest_clique_size_2=length(largest_cliques_2[[1]]),
      largest_clique_size_3=length(largest_cliques_3[[1]]),
      num_maximal_3cliques_2=length(maximal_3cliques_2),
      avgsize_maximal_3cliques_2=mean(sapply(maximal_3cliques_2, length)),
      num_maximal_3cliques_3=length(maximal_3cliques_3),
      avgsize_maximal_3cliques_3=mean(sapply(maximal_3cliques_3, length)),
      modularity_optcommunity_weighted_2=modularity(optcommunity_weighted_2),
      modularity_optcommunity_binary_2=modularity(optcommunity_binary_2),
      avgsize_optcommunity_weighted_2=mean(sizes(optcommunity_weighted_2)),
      number_optcommunity_weighted_2=length(sizes(optcommunity_weighted_2)),
      modularity_optcommunity_weighted_3=modularity(optcommunity_weighted_3),
      modularity_optcommunity_binary_3=modularity(optcommunity_binary_3),
      avgsize_optcommunity_weighted_3=mean(sizes(optcommunity_weighted_3)),
      number_optcommunity_weighted_3=length(sizes(optcommunity_weighted_3)),
      density_2=density_2,
      density_3=density_3,
      transitivity_2=transitivity_2,
      transitivity_3=transitivity_3,
      avgpath_2=avgpath_2,
      avgpath_3=avgpath_3,
      smallworld_2=smallworld_2,
      smallworld_3=smallworld_3
  )
  
  graph_agg <- rbind(graph_agg, graph_measures)
  vertex_agg <- rbind(vertex_agg, vertex_measures)
  
}

save(graph_agg, vertex_agg, file="data/ego_graph_measures.RData")


couple_vertex_metrics <- c()
couple_graph_metrics <- c()

##DYADIC METRICS
for (p in couple_graphs) {
  #consider whether to include ego-based edges since we may not consider an alter truly connected to ego if the closeness rating is weak 
  #g <- graph.data.frame(p$edges_withego, directed=FALSE, vertices=p$actors_withego)
  
  p$actors$role_cat <- factor(sapply(p$actors$role, function(x) { 
            if (x %in% c("patient_only", "patient")) { "patient" 
            } else if (x %in% c("partner_only", "partner")) { "partner"
            } else "shared"
          }))
  g <- graph.data.frame(p$edges_avgratings, directed=FALSE, vertices=p$actors)
  #print(g, e=TRUE, v=TRUE)
  
  #remove any connections that are not close at all
  gthresh2 <- delete.edges(g, which(E(g)$weight <= 2)) #1=have never met, 2=not close at all
  gthresh3 <- delete.edges(g, which(E(g)$weight <= 3)) #1=have never met, 2=not close at all, 3=slightly close
  gbin2 <- remove.edge.attribute(gthresh2, "weight") #unweighted graph <= 2
  gbin3 <- remove.edge.attribute(gthresh3, "weight") #unweighted graph <= 3
  
  #node/vertex metrics
  
  #degree centrality (degree does not account for weight)
  degree_2 <- degree(gthresh2) #element for each vertex
  centralization_degree_2 <- centralization.degree(gthresh2)$centralization #graph level centrality index according to degree
  
  degree_3 <- degree(gthresh3) #element for each vertex
  centralization_degree_3 <- centralization.degree(gthresh3)$centralization #graph level centrality index
  
  #weighted equivalent is called strength and represents sum of weights for connected vertices
  #should probably re-weight the closeness ratings to be 1, 2, 3 or something
  strength_2 <- graph.strength(gthresh2)
  strength_3 <- graph.strength(gthresh3)
  
  #eigenvector centrality
  evcent_weighted_2 <- evcent(gthresh2)$vector #element for each vertex
  evcent_binary_2 <- evcent(gbin2)$vector #element for each vertex
  centralization_eigenvector_2 <- centralization.evcent(gbin2)$centralization
  
  #eigenvector centrality
  evcent_weighted_3 <- evcent(gthresh3)$vector #element for each vertex
  evcent_binary_3 <- evcent(gbin3)$vector #element for each vertex
  centralization_eigenvector_3 <- centralization.evcent(gbin3)$centralization
  
  #closeness centrality
  closeness_weighted_2 <- closeness(gthresh2, normalized=TRUE)*100 #normalized and multiplied by 100 to match ucinet.
  closeness_binary_2 <- closeness(gbin2, normalized=TRUE)*100 #normalized and multiplied by 100 to match ucinet.
  
  closeness_weighted_3 <- closeness(gthresh3, normalized=TRUE)*100 #normalized and multiplied by 100 to match ucinet.
  closeness_binary_3 <- closeness(gbin3, normalized=TRUE)*100 #normalized and multiplied by 100 to match ucinet.
  
  betweenness_weighted_2 <- betweenness(gthresh2, normalize=TRUE)*100
  betweenness_binary_2 <- betweenness(gbin2, normalize=TRUE)*100
  
  betweenness_weighted_3 <- betweenness(gthresh3, normalize=TRUE)*100
  betweenness_binary_3 <- betweenness(gbin3, normalize=TRUE)*100
  
  couple_edge <- get.edge.ids(gthresh2, c(V(gthresh2)$name[V(gthresh2)$role=="partner"], V(gthresh2)$name[V(gthresh2)$role == "patient"]))
  edge_betweenness_weighted_2 <- edge.betweenness(gthresh2, couple_edge)
  edge_betweenness_binary_2 <- edge.betweenness(gbin2, couple_edge)
  edge_betweenness_weighted_3 <- edge.betweenness(gthresh3, couple_edge)
  edge_betweenness_binary_3 <- edge.betweenness(gbin3, couple_edge)
  
  #do alters connect with each other because of similarities on vertex attributes? (esp. attachment
  #not relevant to couple? Since we have different ratings
  #attach_assortativity <- assortativity(gthresh2, V(g)$attach_total)
  #angry_assortativity <- assortativity(gthresh2, V(g)$Angry_num)
  #yearsknown_assortativity <- assortativity(gthresh2, V(g)$YrsKnwn)
  nomination_assortativity_2 <- assortativity.nominal(gthresh2, V(g)$role_cat)
  nomination_assortativity_3 <- assortativity.nominal(gthresh3, V(g)$role_cat)
  
  #TODO: add romance, cutoff, happy, and emotion (categorical) to list of attributes to test for assortativity
  #also consider whether to test each attachment item separately here.
  
  #something to consider adding: number of cohesive blocks
  #b <- cohesive.blocks(gthresh3)
  #plot(b, g)
  
  #cliques (which are based on unweighted edges)
  largest_cliques_2 <- largest.cliques(gbin2)
  maximal_3cliques_2 <- maximal.cliques(gbin2, min=3) #maximal cliques of size 3 or larger
  
  largest_cliques_3 <- largest.cliques(gbin3)
  maximal_3cliques_3 <- maximal.cliques(gbin3, min=3) #maximal cliques of size 3 or larger
  
  #community detection using modularity optimization algorithm
  optcommunity_weighted_2 <- optimal.community(gthresh2)
  optcommunity_binary_2 <- optimal.community(gbin2)
  optcommunity_weighted_3 <- optimal.community(gthresh3)
  optcommunity_binary_3 <- optimal.community(gbin3)
  
  #community detection using edge betweenness algorithm
  #omitting for now
#  ebc_binary_2 <- edge.betweenness.community(gbin2)
#  ebc_weighted_2 <- edge.betweenness.community(gthresh2)
#  ebc_binary_3 <- edge.betweenness.community(gbin3)
#  ebc_weighted_3 <- edge.betweenness.community(gthresh3)
  
  #comparing parcellations
  #compare(ebc_binary_2, ebc_binary_3, method="nmi")
  
  #eccentricity (shortest path from each vertex to furthest point away in graph)
  #part of small worldness -- only sensible for binary graphs
  #note: eccentricity is zero for a vertex that is unconnected to any other
  #when the ego edges are left out, this would mean an alter connected only to self would have eccentricity = 0
  eccentricity_2 <- eccentricity(gbin2)
  eccentricity_3 <- eccentricity(gbin3)
  
  #density of graphs (number of edges versus number of possible edges
  density_2 <- graph.density(gbin2, loops=FALSE)
  density_3 <- graph.density(gbin3, loops=FALSE)
  
  #local clustering
  locclust_binary_2 <- transitivity(gbin2, type="local") #will be NaN for a vertex with no connected triples
  locclust_weighted_2 <- transitivity(gthresh2, type="weighted") #will be NaN for a vertex with no connected triples
  locclust_binary_3 <- transitivity(gbin3, type="local") #will be NaN for a vertex with no connected triples
  locclust_weighted_3 <- transitivity(gthresh3, type="weighted") #will be NaN for a vertex with no connected triples
  
  #overall transitivity (global clustering) of graph
  transitivity_2 <- transitivity(gbin2, type="global")
  transitivity_3 <- transitivity(gbin3, type="global")
  
  #average path length of graph
  avgpath_2 <- average.path.length(gbin2, unconnected=TRUE)
  avgpath_3 <- average.path.length(gbin3, unconnected=TRUE)
  
  #figure out small worldness of network according to Humphries 2008
  #generate an E-R random graph with the same number of nodes and edges
  
  g_rand_2 <- erdos.renyi.game(n=vcount(gbin2), p.or.m=ecount(gbin2), type=c("gnm"),directed = FALSE, loops = FALSE)
  transitivity_rand_2 <- transitivity(g_rand_2, type="global")
  avgpath_rand_2 <- average.path.length(g_rand_2, unconnected=TRUE)
  
  g_rand_3 <- erdos.renyi.game(n=vcount(gbin3), p.or.m=ecount(gbin3), type=c("gnm"),directed = FALSE, loops = FALSE)
  transitivity_rand_3 <- transitivity(g_rand_3, type="global")
  avgpath_rand_3 <- average.path.length(g_rand_3, unconnected=TRUE)
  
  #from humphries 2008
  smallworld_2 <- (transitivity_2/transitivity_rand_2)/(avgpath_2/avgpath_rand_2)
  smallworld_3 <- (transitivity_3/transitivity_rand_3)/(avgpath_3/avgpath_rand_3)
  
  #assemble vertex measures and graph measures into distinct data.frames
  
  #add averaged vertex measures to graph measures for descriptives
  
  vertex_measures <- data.frame(
      PTNUM=p$PTNUM,
      alter_name=V(g)$name,
      degree_2=degree_2,
      degree_3=degree_3,
      strength_2=strength_2,
      strength_3=strength_3,
      evcent_weighted_2=evcent_weighted_2,
      evcent_binary_2=evcent_binary_2,
      evcent_weighted_3=evcent_weighted_3,
      evcent_binary_3=evcent_binary_3,
      closeness_weighted_2=closeness_weighted_2,
      closeness_binary_2=closeness_binary_2,
      closeness_weighted_3=closeness_weighted_3,
      closeness_binary_3=closeness_binary_3,
      betweenness_weighted_2=betweenness_weighted_2,
      betweenness_binary_2=betweenness_binary_2,
      betweenness_weighted_3=betweenness_weighted_3,
      betweenness_binary_3=betweenness_binary_3,
      eccentricity_2=eccentricity_2,
      eccentricity_3=eccentricity_3,
      locclust_binary_2=locclust_binary_2,
      locclust_weighted_2=locclust_weighted_2,
      locclust_binary_3=locclust_binary_3,
      locclust_weighted_3=locclust_weighted_3      
  )
  
  graph_measures <- data.frame(
      PTNUM=p$PTNUM,
      centralization_degree_2=centralization_degree_2,
      centralization_degree_3=centralization_degree_3,
      centralization_eigenvector_2=centralization_eigenvector_2,
      centralization_eigenvector_3=centralization_eigenvector_3,
      #attach_assortativity=attach_assortativity,
      #angry_assortativity=angry_assortativity,
      #yearsknown_assortativity=yearsknown_assortativity,
      nomination_assortativity_2=nomination_assortativity_2,
      nomination_assortativity_3=nomination_assortativity_3,
      largest_clique_size_2=length(largest_cliques_2[[1]]),
      largest_clique_size_3=length(largest_cliques_3[[1]]),
      num_maximal_3cliques_2=length(maximal_3cliques_2),
      avgsize_maximal_3cliques_2=mean(sapply(maximal_3cliques_2, length)),
      num_maximal_3cliques_3=length(maximal_3cliques_3),
      avgsize_maximal_3cliques_3=mean(sapply(maximal_3cliques_3, length)),
      modularity_optcommunity_weighted_2=modularity(optcommunity_weighted_2),
      modularity_optcommunity_binary_2=modularity(optcommunity_binary_2),
      avgsize_optcommunity_weighted_2=mean(sizes(optcommunity_weighted_2)),
      number_optcommunity_weighted_2=length(sizes(optcommunity_weighted_2)),
      modularity_optcommunity_weighted_3=modularity(optcommunity_weighted_3),
      modularity_optcommunity_binary_3=modularity(optcommunity_binary_3),
      avgsize_optcommunity_weighted_3=mean(sizes(optcommunity_weighted_3)),
      number_optcommunity_weighted_3=length(sizes(optcommunity_weighted_3)),
      density_2=density_2,
      density_3=density_3,
      transitivity_2=transitivity_2,
      transitivity_3=transitivity_3,
      avgpath_2=avgpath_2,
      avgpath_3=avgpath_3,
      smallworld_2=smallworld_2,
      smallworld_3=smallworld_3,
      edge_betweenness_weighted_2=edge_betweenness_weighted_2,
      edge_betweenness_weighted_3=edge_betweenness_weighted_3,
      edge_betweenness_binary_2=edge_betweenness_binary_2,
      edge_betweenness_binary_3=edge_betweenness_binary_3 #NB: The edge betweenness measures are in graph because they are just for the EBC for the couple edge
  )
  
  couple_graph_metrics <- rbind(couple_graph_metrics, graph_measures)
  couple_vertex_metrics <- rbind(couple_vertex_metrics, vertex_measures)
  
}

save(couple_graph_metrics, couple_vertex_metrics, file="data/couple_graph_measures.RData")
