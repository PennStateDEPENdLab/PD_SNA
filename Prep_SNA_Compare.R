library(RSQLite)
library(foreign)
setwd(file.path(getMainDir(), "EIFB_Longitudinal"))

#Open DB connections
#sna1 <- dbConnect(dbDriver("SQLite"), dbname=file.path(getMainDir(), "Miscellaneous", "Shared_Data", "SocialNetwork_19Apr2012.sqlite"))
#sna2 <- dbConnect(dbDriver("SQLite"), dbname=file.path(getMainDir(), "Miscellaneous", "Shared_Data", "SocialNetwork_Double_19Apr2012.sqlite"))
sna1 <- dbConnect(dbDriver("SQLite"), dbname=file.path(getMainDir(), "Miscellaneous", "Shared_Data", "SocialNetwork_27Nov2012.sqlite"))
sna2 <- dbConnect(dbDriver("SQLite"), dbname=file.path(getMainDir(), "Miscellaneous", "Shared_Data", "SocialNetwork_Double_27Nov2012.sqlite"))


#dbGetInfo(sna1)
#dbGetInfo(sna2)
#dbListTables(sna1)
#dbListTables(sna2)

sna1.alter <- dbReadTable(sna1, "[Alter]")
sna1.attachRating <- dbReadTable(sna1, "AttachmentRating")
sna1.attachScreening <- dbReadTable(sna1, "AttachmentScreening")

sna2.alter <- dbReadTable(sna2, "[Alter]")
sna2.attachRating <- dbReadTable(sna2, "AttachmentRating")
sna2.attachScreening <- dbReadTable(sna2, "AttachmentScreening")

#Access and SQLite don't support full outer joins, but this would be the full SQL syntax to merge the tables properly
#sqlstyle <- dbGetQuery(sna1, "SELECT * FROM (AttachmentScreening a_s OUTER JOIN AttachmentRating a_r ON (a_s.AlterPK = a_r.AlterPK AND a_s.Month = a_r.Month)) LEFT OUTER JOIN [Alter] a ON a.AlterPK = a_s.AlterPK")

#Close DB connections
dbDisconnect(sna1)
dbDisconnect(sna2)

#print out contents of each table
#str(sna1.alter)
#str(sna1.attachRating)
#str(sna1.attachScreening)

fieldOrder <- c("PTNUM", "Month", "Alter_ID", "Name", "ChangeDistance", "Died", "ReasonCutoff", "prxsk", "sepds", "sfhvn", "secbs", "Age", "Sex",
    "rshp", "lngth", "lngyr", "lngmo", "roman", "face", "phone", "email", "cutof", "argue", "close", "trust", "crtcl", "advic", "suprt")

#merge screening and rating first -- because there may be screening records that have no rating, so want to fully propagate alter fields below
sna1.merge <- merge(subset(sna1.attachRating, select=-c(DateCreated, AttachmentRatingID)),
    subset(sna1.attachScreening, select=-c(DateCreated, AttachmentScreeningID)),
    by=c("AlterPK", "Month"), all=TRUE)

sna1.merge <- merge(sna1.merge,
    subset(sna1.alter, select=-c(DateCreated, SO)),
    by=c("AlterPK"), all=TRUE)

#order by PTNUM, Month, Alter_ID; and order
sna1.merge <- sna1.merge[with(sna1.merge, order(PTNUM, Month, Alter_ID)), fieldOrder]

#drop alter pk
sna1.merge$AlterPK <- NULL
#sna1.merge$RowNum <- 1:nrow(sna1.merge)

####sna double
sna2.merge <- merge(subset(sna2.attachRating, select=-c(DateCreated, AttachmentRatingID)),
    subset(sna2.attachScreening, select=-c(DateCreated, AttachmentScreeningID)),
    by=c("AlterPK", "Month"), all=TRUE)

sna2.merge <- merge(sna2.merge,
    subset(sna2.alter, select=-c(DateCreated, SO, Month)), #somehow double DB has dummy month in Alter table
    by=c("AlterPK"), all=TRUE)

#order by PTNUM, Month, Alter_ID; and order fields more intuitively
sna2.merge <- sna2.merge[with(sna2.merge, order(PTNUM, Month, Alter_ID)), fieldOrder]

#drop alter pk
sna2.merge$AlterPK <- NULL
#sna2.merge$RowNum <- 1:nrow(sna2.merge)

#write.dta(sna1.merge, "SNASingleMerge_19Apr2012.dta")
#write.dta(sna2.merge, "SNADoubleMerge_19Apr2012.dta")

write.dta(sna1.merge, "SNASingleMerge_27Nov2012.dta")
write.dta(sna2.merge, "SNADoubleMerge_27Nov2012.dta")

