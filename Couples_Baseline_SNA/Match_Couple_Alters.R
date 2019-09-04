##the txt file has detailed start and stop times for all codes
datadir <- "/Users/mnh5174/Box_Sync/DEPENd/Projects/PD_SNA/ties"
setwd(datadir)
batchmode <- TRUE
if (!batchmode) {
  library(tcltk)
  
  cat("Choose the Alter list for participant 1.\n")
  p1file <- "" 
  while(p1file=="") {
    p1file <- tclvalue(tkgetOpenFile(title="Choose the alter_summary.csv file for participant 1", filetypes = "{{CSV files} {.csv}} {{All files} *}"))
  }
  
  cat("Choose the Alter list for participant 2.\n")
  p2file <- "" 
  while(p2file=="") {
    p2file <- tclvalue(tkgetOpenFile(title="Choose the alter_summary.csv file for participant 2", filetypes = "{{CSV files} {.csv}} {{All files} *}"))
  }  
} else {
  subdirs <- list.dirs(path=datadir, recursive=FALSE)
  pfiles <- list()
  for (p in 1:length(subdirs)) {
    pfiles[[p]] <- list.files(path=subdirs[p], pattern=".*alter_summary.csv", full.names = TRUE, recursive = TRUE)  
  }
 
  nfiles <- sapply(pfiles, length)
}



runmatch <- function(p1file, p2file) {
  stopifnot(file.exists(p1file) && file.exists(p2file))
  
  p1alters <- read.csv(p1file, header=TRUE)
  p2alters <- read.csv(p2file, header=TRUE)
  
  #capitalization inconsistencies
  names(p1alters) <- tolower(names(p1alters))
  names(p2alters) <- tolower(names(p2alters))
  
  #use adist (string distance) to look for matches on alter
  distmat <- matrix(NA_real_, nrow=length(p1alters$alter_name), ncol=length(p2alters$alter_name))
  for (i in 1:length(p1alters$alter_name)) {
    distmat[i,] <- sapply(p2alters$alter_name, function(partner_name) {
      adist(p1alters$alter_name[i], partner_name, ignore.case=TRUE, partial=FALSE)
    })
  }
  
  distmat_score <- as.numeric(distmat <= 4) #0/1 for whether the match is close
  
  #patient_bestmatch <- apply(distmat, 1, which.min) #for each patient alter (1-25), what is the position of the nearest partner alter (1-25)
  #patient_disttomatch <- apply(distmat, 1, min) #how close is the potential match
  
  #within_distance <- as.numeric(patient_disttomatch <= 4) #4 distance seems about right to leave out horrible matches
  #use additional criteria below to refine matches (i.e., above is overly inclusive)
  
  #an alter can only be matched once, and duplicate matches are indicative of a false positive.
  #for any duplicates, only retain the alter with the lowest distance.
  #if two (or more) names have the same minimum distance, keep them in (weed out using other metrics below)
  
  #use a logical and of within_distance and is_exclusive_best to generate alter match
  exclusivemat <- matrix(NA_integer_, nrow=length(p1alters$alter_name), ncol=length(p2alters$alter_name))
  for (i in 1:length(p1alters$alter_name)) {
    for (j in 1:length(p2alters$alter_name)) {
      if (distmat[i,j] == min(distmat[i,])) { exclusivemat[i,j] <- 1
      } else exclusivemat[i,j] <- 0
    }
  }
  
  
  #
  ##also verify that age is within 10 years
  #age_diff <- abs(p1alters$age[alter_match] - p2alters$age[patient_bestmatch[alter_match]]) <= 10
  #
  ##if age is missing, an NA will be generated for age_diff... in this case, default to non-match (no point)
  #age_diff[which(is.na(age_diff))] <- FALSE
  #
  #age_diff <- as.numeric(age_diff)
  
  
  lastinitialmat <- matrix(NA_integer_, nrow=length(p1alters$alter_name), ncol=length(p2alters$alter_name))
  
  #filter on last initial match
  last_initial_patient <- tolower(sub("^\\s*\\w+.*[\\s_]+([A-z]+)[_\\s]*$", "\\1", p1alters$alter_name, perl=TRUE))
  last_initial_partner <- tolower(sub("^\\s*\\w+.*[\\s_]+([A-z]+)[_\\s]*$", "\\1", p2alters$alter_name, perl=TRUE))
  
  for (i in 1:length(last_initial_patient)) {
    for (j in 1:length(last_initial_partner)) {
      minlen <- min(c(nchar(last_initial_patient[i]), nchar(last_initial_partner[j])))
      lastinitialmat[i,j] <- substr(last_initial_patient[i], 1, minlen) == substr(last_initial_partner[j], 1, minlen)    
    } 
  }
  
  
  firstinitialmat <- matrix(NA_integer_, nrow=length(p1alters$alter_name), ncol=length(p2alters$alter_name))
  
  #filter on first initial match
  first_initial_patient <- tolower(sub("^\\s*([A-z])\\w+.*[\\s_]+[A-z]+[_\\s]*$", "\\1", p1alters$alter_name, perl=TRUE))
  first_initial_partner <- tolower(sub("^\\s*([A-z])\\w+.*[\\s_]+[A-z]+[_\\s]*$", "\\1", p2alters$alter_name, perl=TRUE))
  
  for (i in 1:length(first_initial_patient)) {
    for (j in 1:length(first_initial_partner)) {
      minlen <- min(c(nchar(first_initial_patient[i]), nchar(first_initial_partner[j])))
      firstinitialmat[i,j] <- substr(first_initial_patient[i], 1, minlen) == substr(first_initial_partner[j], 1, minlen)    
    } 
  }
  
  matchscore <- distmat_score + exclusivemat + lastinitialmat + firstinitialmat
  
  clunkyRbind <- function(mat, vec) {
    if (!is.null(mat) && length(vec) < ncol(mat)) {
      vec <- c(vec, rep(NA_character_, ncol(mat) - length(vec)))
    }
    if (!is.null(mat) && ncol(mat) < length(vec)) {
      mat <- cbind(mat, matrix(NA_character_, nrow=nrow(mat), ncol=length(vec) - ncol(mat)))
    }  
    mat <- rbind(mat, vec)
  }
  
  #write to file
  saveFile <- ""
  while (saveFile=="") {
    saveFile <- tclvalue(tkgetSaveFile(title="Save alter match file as...", defaultextension="txt", filetypes = "{{TXT files} {.txt}} {{All files} *}", initialfile="Alter_match.txt"))  
  }
  
  sink(file=saveFile, split=TRUE)
  
  p1best <- apply(matchscore, 1, function(r) {
    as.character(p2alters$alter_name[which(r==max(r))])      
  })
  
  
  matp1 <- c()
  
  for (i in 1:length(p1best)) {
    matp1 <- clunkyRbind(matp1, p1best[[i]])
  }
  
  #limit to 5 matches
  if (ncol(matp1) > 5) {
    matp1 <- matp1[,1:5]
  }
  
  matp1 <- cbind(as.character(p1alters$alter_name), matp1)
  
  colnames(matp1) <- c("p1_name", paste0("p2_match", 1:(ncol(matp1) - 1)))
  rownames(matp1) <- NULL
  
  cat("Best alter matches with master: ", p1file, "\n  and comparator: ", p2file, "\n\n")
  print(data.frame(matp1), row.names=FALSE)
  
  #now comparison with p2 as master
  p2best <- apply(matchscore, 2, function(r) {
    as.character(p1alters$alter_name[which(r==max(r))])      
  })
  
  matp2 <- c()
  for (i in 1:length(p2best)) {
    matp2 <- clunkyRbind(matp2, p2best[[i]])
  }
  
  #limit to 5 matches
  if (ncol(matp2) > 5) {
    matp2 <- matp2[,1:5]
  }
  
  matp2 <- cbind(as.character(p2alters$alter_name), matp2)
  
  colnames(matp2) <- c("p2_name", paste0("p1_match", 1:(ncol(matp2) - 1)))
  rownames(matp2) <- NULL
  
  cat("\n\nBest alter matches with master: ", p2file, "\n  and comparator: ", p1file, "\n\n")
  print(data.frame(matp2), row.names=FALSE)
  
  sink()
}
