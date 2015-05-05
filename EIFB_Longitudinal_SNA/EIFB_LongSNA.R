library(foreign)
library(plyr)
setwd("/Users/michael/Tresors/PD_SNA/EIFB_Longitudinal_SNA")

newalters <- read.spss("data/SNA new alters 2.12.14.sav", to.data.frame=TRUE)
sna <- read.spss("data/SNA intake&followup 2.12.14.sav", to.data.frame=TRUE)

#order by ptnum, month
newalters <- newalters[order(newalters$PTNUM, newalters$Month),]
sna <- sna[order(sna$PTNUM, sna$Month),]

head(sna)
head(newalters)

table(sna$ChangeDistance, useNA="always")

#attachment ratings were obtained at each follow-up
table(sna$prxsk, useNA="always") #proximity seeking: it is important to me to see or talk with X regularly.
table(sna$sepds, useNA="always") #separation distance: X is a person I do not like to be away from.
table(sna$sfhvn, useNA="always") #safe haven: X is the first person that I would turn to if I had a problem.
table(sna$secbs, useNA="always") #secure base: X is the person that I would actually count on to always be there for me and care for me no matter what.

#trim fixed width character import
sna$ReasonCutoff <- factor(gsub("^\\s+|\\s+$", "", sna$ReasonCutoff, perl=TRUE))
sna$Name <- factor(gsub("^\\s+|\\s+$", "", sna$Name, perl=TRUE))

#so looks like sna really has the follow-up data of interest...
#ChangeDistance:
#  8 = not applicable (likely due to death)
#  9 = missing data
sna[which(sna$ChangeDistance %in% c(8,9)),]

#7084 is missing 3, 6, 9 month data. There are also a ton of cutoffs at 12 months
#Can't reasonably carry over data at 12 months based on the baseline
sna <- subset(sna, PTNUM!=7084)

#Set change distance for alters who died (ChangeDistance==8)
#Currently entered such that code is 8 at the first assessment where alter died, but then set to 0 (no change) thereafter
snaFill <- ddply(sna, .(PTNUM), function(subdf) {
      alterdeath <- subset(subdf, ChangeDistance==8)
      if (nrow(alterdeath) > 0L) {
        #verify that ratings are missing for this and subsequent ratings
        for (i in 1:nrow(alterdeath)) {
          futureratings <- with(subdf, which(AlterID == alterdeath[i, "AlterID"] & Name == alterdeath[i, "Name"] & Month >= alterdeath[i, "Month"]))
          if (length(futureratings) > 0L) { 
            #for tracking down records to verify
            #print(subdf[futureratings,c("PTNUM", "Month", "AlterID", "ChangeDistance", "Name", "prxsk", "sepds", "sfhvn", "secbs")])
            subdf[futureratings, c("ChangeDistance", "prxsk", "sepds", "sfhvn", "secbs")] <- t(replicate(length(futureratings), c(8, rep(NA_integer_, 4)))) 
          }
        }
      }
      subdf
    })

#to verify:
#    PTNUM Month AlterID ChangeDistance    Name prxsk sepds sfhvn secbs
#52   7154     3      19              8 andre m    NA    NA    NA    NA
#81   7154     6      19             -1 andre m     3     3     0     3
#110  7154     9      19              8 andre m    NA    NA    NA    NA
#139  7154    12      19              8 andre m    NA    NA    NA    NA

#    PTNUM Month AlterID ChangeDistance      Name prxsk sepds sfhvn secbs
#74   7110     9       6              8 Michael H     6     6     6     6
#112  7110    12       6              8 Michael H    NA    NA    NA    NA

#    PTNUM Month AlterID ChangeDistance     Name prxsk sepds sfhvn secbs
#37   7163     3       7              8 debbie s    NA    NA    NA    NA
#67   7163     6       7              8 debbie s    NA    NA    NA    NA
#97   7163     9       7              8 debbie s     6     6     6     6
#127  7163    12       7              8 debbie s    NA    NA    NA    NA

#    PTNUM Month AlterID ChangeDistance       Name prxsk sepds sfhvn secbs
#70   7165     6       2              8 clarence b    NA    NA    NA    NA
#114  7165     9       2              8 clarence b     5     5     3     5
#158  7165    12       2              0 clarence b    NA    NA    NA    NA

#    PTNUM Month AlterID ChangeDistance            Name prxsk sepds sfhvn secbs
#40   7166     3      10              8 Angie Armstrong    NA    NA    NA    NA
#70   7166     6      10              8 Angie Armstrong     0     0     0     0
#104  7166     9      10              8 Angie Armstrong    NA    NA    NA    NA
#138  7166    12      10              8 Angie Armstrong    NA    NA    NA    NA

#when change distance zero, fill in four attachment items from prior assessment, if available.
snaFill <- ddply(snaFill, .(PTNUM), function(subdf) {
      #if (subdf$PTNUM[1L]==7164) {browser()}
      rowsToFill <- with(subdf, which(Month > 0 & ChangeDistance==0))
      if (length(rowsToFill) > 0L) {
        for (r in rowsToFill) {
          curMonth <- subdf[r, "Month"]
          if (curMonth == 3) { prevMonth <- 0 
          } else if (curMonth == 6) { prevMonth <- 3
          } else if (curMonth == 9) { prevMonth <- 6
          } else if (curMonth == 12) { prevMonth <- 9
          } else { stop("Cannot determine previous month.")}
          
          candidateMatch <- subdf[which(subdf$Month <= prevMonth & subdf$AlterID == subdf[r, "AlterID"] & subdf$Name == subdf[r, "Name"]), ]
          
          #if (r == 109 && subdf$PTNUM[1L] == 7164) { browser()}
          
          if (nrow(candidateMatch) > 0L) {
            if (nrow(candidateMatch) > 1L) {
              if (length(unique(candidateMatch$Month)) > nrow(candidateMatch)) { stop ("Multiple possible previous records per month")
              } else {
                candidateMatch <- candidateMatch[order(candidateMatch$Month, decreasing=TRUE)[1L],] #just the most recent match
              }
            }
            
            subdf[r, c("prxsk", "sepds", "sfhvn", "secbs")] <- candidateMatch[, c("prxsk", "sepds", "sfhvn", "secbs")]
            
          } else { print("Cannot find record to fill in") } #; browser() }
        }
      }
      return(subdf)
    }
)

table(subset(snaFill, Month > 0)$prxsk, useNA="always")

missprxsk <- snaFill[which(is.na(snaFill$prxsk) & snaFill$Month > 0 & snaFill$ChangeDistance != 8), c("PTNUM", "AlterID", "Month", "ChangeDistance", "prxsk", "sepds", "sfhvn", "secbs")]

write.table(missprxsk, file="Missing_attachment_ratings_5May2015.txt", row.names=FALSE, quote=FALSE)

#get all rows for PTNUM and AlterID where rating is missing
b <- apply(missprxsk, 1, function(r) {
      snaFill[which(snaFill$PTNUM == r["PTNUM"] & snaFill$AlterID == r["AlterID"]),]
    })

missinglookup <- do.call("rbind", b)[, c("PTNUM", "Month", "ChangeDistance", "AlterID", "prxsk", "sepds", "sfhvn", "secbs")]

write.table(missinglookup, file="All_ratings_for_missing_5May2015.txt", row.names=FALSE, quote=FALSE)

str(sna[which(sna$PTNUM %in% missprxsk$PTNUM & sna$AlterID %in% missprxsk$AlterID),])

#7084 is missing records at 3, 6, and 9 -- maybe worth dropping altogether
subset(sna, PTNUM==7084)
