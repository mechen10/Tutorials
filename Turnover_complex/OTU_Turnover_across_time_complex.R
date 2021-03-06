#!/bin/bash

# This script is to find number/percent of OTUs gained and lost in a time series

######################### Data variables ############################

# Set path. Change this path to folder that you want output to be generated
setwd("/Users/melissachen/Documents/Masters/UNIVERSALCODE_git/Tutorials/Turnover/")

# This is the file path to your OTU table 
# MUST BE COUNT DATA; NOT relative abundance.
# If relative abundance, please comment out section indicated below.
OTUFP <- "./OTU_Table_text.txt"

# This is the file path to your Mapping file.
mfFP <- "./MF_nochlpmito_m1000.txt"

# Minimum number of samples an OTU must occur in for it NOT to be filtered out
MinOccur <- 3

# Mapping file header name with treatments
# Can be numeric or character
treatment <- "Type"

# Mapping file header name with Replicates
# Can be numeric or character
replicates <- "Morph"
#* Note, I'm using 'Morph' as replicates because that's what makes sense for my data-- you adjust yours as necessary.

# Mapping file header name with time (or other metric you are plotting by)
# MUST BE NUMERIC
byFactor <- "Time"

# Mapping file header name and factor that designates whether there are samples you don't want to keep
# For example, if you only want to compare rat fecal swabs ('Fecal')
# but also have food swab samples ('Food') in your data set, you would do something like:
# Type:Fecal
# Can be numeric or character
toKeep <- "Description:AM"

# Color list; in same order as treatments. Comma separated.
colList <- "red,blue"

######################### Load data ############################

OTUTable <- read.delim(paste0(OTUFP)
                       , header = TRUE
                       , skip = 1 # First row should be a comment made by QIIME
                       , row.names = 1
                       , check.names = FALSE # Prevents R from changing header names
                       )

MF <- read.delim(paste0(mfFP)
                 , header = TRUE
                 , row.names = 1
                 , stringsAsFactors = FALSE
)

# Make sure all sites are in OTU and MF; and vice versa

OTUTable <- OTUTable[,colnames(OTUTable) %in% rownames(MF)]
MF <- MF[rownames(MF) %in% colnames(OTUTable),]

# Split of 'toKeep' to two two parts; the header name and factor name
toKeep <- unlist(strsplit(toKeep, ":"))

# Make colList
colList <- unlist(strsplit(colList,","))
######################### Set up data ############################

# This data has TWO treatment types across time: Hakai and Port Moody (aka Reed Point)
# It calculates and plots turn over for BOTH sites in a single plot. Adjust as necessary for your own data.

# Filter out all samples that are NOT AM swabs (water, other surfaces, etc)
MF <- MF[grep(paste0(toKeep[2]), MF[,paste0(toKeep[1])]),]
OTUTable <- OTUTable[,match(rownames(MF), colnames(OTUTable))]

# Get list of Replicates
replicatesList <- unique(MF[,replicates])

# Split up two datasets; Hakai and Port Moody
# First, get all treatment names
treatmentList <- levels(factor(MF[,treatment]))

# Add names to colors
names(colList) <- treatmentList

allTreatments <- list()
for (treat in treatmentList) {
  allTreatments[[treat]] <- list()
  allTreatments[[treat]][["overall"]] <- list()
  
  # Filter out all samples for treatment
  MF.temp <- MF[which(MF[,treatment] == treat),]
  allTreatments[[treat]][["overall"]][["MF"]]<- MF.temp
  # Filter OTU table based on site names in MF.H
  OTUTable.temp <- OTUTable[,match(rownames(MF.temp), colnames(OTUTable))]
  allTreatments[[treat]][["overall"]][["OTUTable"]]  <- OTUTable.temp
   
  # Delete things that don't occur more than 'MinOccur' times in the filtered data set
  toDelete <- c()
  for (r in 1:nrow(OTUTable.temp)) {
    nOccurance <- sum(OTUTable.temp[r,] > 0.0)
    if (nOccurance < MinOccur) {
      toDelete <- c(toDelete, r)
    }
  }
  if (length(toDelete) > 0) {
    OTUTable.temp.filt <- OTUTable.temp[-toDelete,]
  }
  
  # Make presence/absence
  # HAKAI
  OTUTable.temp.presabs <- OTUTable.temp.filt
  for (c in 1:ncol(OTUTable.temp.filt)) {
    for (r in 1:nrow(OTUTable.temp.filt)) {
      if (OTUTable.temp.filt[r,c] > 0.0) {
        OTUTable.temp.presabs[r,c] <- 1
      } else {
        OTUTable.temp.presabs[r,c] <- 0
      }
    }
  }
  allTreatments[[treat]][["overall"]][["OTUTable.presabs"]]  <- OTUTable.temp.presabs
  
  # Order factors into sorted levels
  sortedFactors <- sort(as.numeric(unique(MF.temp[,byFactor])))
  MF.temp[,byFactor] <- factor(MF.temp[,byFactor], levels = sortedFactors) 
  
  # Make into list
  timesList.temp <- levels(MF.temp[,byFactor])
  
  # Empty matrix with 4 types of OTUs; total, lost, gained, and retained
  turnOver.temp <- matrix(nrow = 4, ncol = length(timesList.temp))
  rownames(turnOver.temp) <- c("total","lost","gained", "retained")
  colnames(turnOver.temp) <- timesList.temp
  
  allold <- c()
  for (i in 1:length(timesList.temp)) {
    ttemp <- timesList.temp[i]
    namesKeep <- rownames(MF.temp)[which(MF.temp[,byFactor] == ttemp)]
    allonly <- OTUTable.temp.presabs[,match(namesKeep, colnames(OTUTable.temp.presabs))]
    
    # Filter by which ones are actually there
    allonly <- allonly[which(rowSums(allonly) != 0),]
    
    # Get total OTU count for this time
    allcount <- nrow(allonly)
    
    # Get "old" list of OTUs and compare which ones are the same
    alllost <- length(allold)-sum(allold %in% rownames(allonly))
    allgain <- allcount-sum(rownames(allonly) %in% allold)
    allretained <- sum(rownames(allonly) %in% allold)
    
    # Load into matrices
    turnOver.temp["total", paste0(ttemp)] <- allcount
    turnOver.temp["lost", paste0(ttemp)] <- round(alllost/allcount,2)
    turnOver.temp["gained", paste0(ttemp)] <- round(allgain/allcount,2)
    turnOver.temp["retained", paste0(ttemp)] <- round(allretained/allcount,2)
    
    # Finally, set "old" as the new ones
    allold <- names(which(rowSums(allonly) != 0))
  }
  allTreatments[[treat]][["overall"]][["turnOver"]] <- turnOver.temp
  
  
  ## REPLICATES SEPARATELY ##
  # Make list for replicates
  replicateTurnOver.temp <- list()
  for (r in replicatesList) {
    # Empty matrix with 4 types of OTUs; total, lost, gained, and retained
    replicateTurnOver.temp[[r]] <- matrix(nrow = 4, ncol = length(timesList.temp))
    rownames(replicateTurnOver.temp[[r]]) <- c("total","lost","gained", "retained")
    colnames(replicateTurnOver.temp[[r]]) <- timesList.temp
    allold <- c()
    i <- 5
    for (i in 1:length(timesList.temp)) {
      ttemp <- timesList.temp[i]
      namesKeep <- rownames(MF.temp)[which((MF.temp[,byFactor] == ttemp) & (MF.temp[,replicates] == r))]
      allonly <- OTUTable.temp.presabs[,match(namesKeep, colnames(OTUTable.temp.presabs))]
      
      # Filter by which ones are actually there
      allonly <- allonly[which(rowSums(allonly) != 0),]
      
      # Get total OTU count for this time
      allcount <- nrow(allonly)
      
      # Get "old" list of OTUs and compare which ones are the same
      alllost <- length(allold)-sum(allold %in% rownames(allonly))
      allgain <- allcount-sum(rownames(allonly) %in% allold)
      allretained <- sum(rownames(allonly) %in% allold)
      
      # Load into matrices
      replicateTurnOver.temp[[r]]["total", paste0(ttemp)] <- allcount
      replicateTurnOver.temp[[r]]["lost", paste0(ttemp)] <- round(alllost/allcount,2)
      replicateTurnOver.temp[[r]]["gained", paste0(ttemp)] <- round(allgain/allcount,2)
      replicateTurnOver.temp[[r]]["retained", paste0(ttemp)] <- round(allretained/allcount,2)
      
      # Finally, set "old" as the new ones
      allold <- names(which(rowSums(allonly) != 0))
    }
  }
  allTreatments[[treat]][["replicateTurnOver"]] <- replicateTurnOver.temp
}



######### MAKE PLOTTING FILES ############
# Make turnOver files into a single table for each metric
metrics <- c("retained","gained","lost")

# Get number of byFactor types
allbyFactors <- c()
for (i in allTreatments) {
  allbyFactors <- c(allbyFactors, colnames(i[["overall"]][["turnOver"]]))
}
lengthFactors <- length(unique(allbyFactors))

# Get names of max lengths, ordered
orderedFactors <- factor(allbyFactors, levels = sort(as.numeric(unique(allbyFactors))))

# Get number of replicates
nameReps <- list()
for (ntype in names(allTreatments)) {
  nameReps[[ntype]] <- names(allTreatments[[ntype]][["replicateTurnOver"]])
}
nReps <- length(unlist(nameReps))

# Get names of overall + replicates
namesList <- c()
for (ntype in names(allTreatments)) {
  namesList <- c(namesList, paste0(ntype,"_"))
  for (names in names(allTreatments[[ntype]][["replicateTurnOver"]])) {
    namesList <- c(namesList, paste0(ntype,"_",names))
  }
}

combinedMetrics <- list()
for (m in metrics) {
  # Make empty matrix
  combinedMetrics[[m]] <- matrix(nrow = length(allTreatments) + nReps, ncol = lengthFactors)
  rownames(combinedMetrics[[m]]) <- namesList
  colnames(combinedMetrics[[m]]) <- levels(orderedFactors)
  
  # Add values
  for (treat in names(allTreatments)) {
    timeList <- colnames(allTreatments[[treat]][["overall"]][["turnOver"]])
    for (t in timeList) {
      Temp <- allTreatments[[treat]][["overall"]][["turnOver"]][paste0(m), t]
      combinedMetrics[[m]][paste0(treat,"_"),paste0(t)] <- Temp
      for (rep in names(allTreatments[[treat]][["replicateTurnOver"]])) {
        Temp <- allTreatments[[treat]][["replicateTurnOver"]][[rep]][paste0(m), t]
        combinedMetrics[[m]][paste0(treat,"_",rep),paste0(t)] <- Temp
      }
    }
  }
}

# # Save files into tab delimited tables
for (metric in names(combinedMetrics)) {
  write.table(combinedMetrics[[metric]]
              , sep = "\t"
              , file = paste0("turnOver_",metric,".txt")
              , col.names = NA # Adds a blank cell in first column so that headers line up
  )
}

####### PLOTTING ##########
# Plot all metrics
xnum <- ncol(combinedMetrics[[1]])-1
xlabBlank <- seq(xnum)
xlabValues <- colnames(combinedMetrics[[1]])[2:ncol(combinedMetrics[[1]])]

for (metric in names(combinedMetrics)) {
  miny <- min(combinedMetrics[[metric]][,2:ncol(combinedMetrics[[metric]])], na.rm = TRUE)
  maxy <- max(combinedMetrics[[metric]][,2:ncol(combinedMetrics[[metric]])], na.rm = TRUE)

  pdf(paste0(metric,".pdf"), pointsize = 14, width = 7, height = 5)
  par(fig = c(0,0.8,0,1), mar = c(6,4,2,4))
  plot(xlabBlank,NULL
       , xlim = c(1,6)
       , ylim = c(miny,maxy)
       , pch= ""
       , xaxt = "n"
       , xlab = ""
       , ylab = paste0("Percent community ",metric)
  )
  title(xlab = "Time"
        , line = 4)
  axis(side = 1
       , at = xlabBlank
       , labels = xlabValues
       , las = 2)
  for (treat in treatmentList) {
    tempTreat <- combinedMetrics[[metric]][grep(paste0(treat,"_"), rownames(combinedMetrics[[metric]])),]
    tempTreat <- tempTreat[,2:ncol(tempTreat)]
    for (r in 1:nrow(tempTreat)) {
      if (rownames(tempTreat)[r] == paste0(treat,"_")) {
        points(tempTreat[r,] ~ xlabBlank
               , col = colList[treat]
               , pch = 19)
        lines(na.exclude(tempTreat[r,]) ~ xlabBlank[!is.na(tempTreat[r,])] 
              , col = colList[treat]
              , lty = "solid"
              , lwd = 2
        )
      } else {
        points(tempTreat[r,] ~ xlabBlank
               , col = colList[treat]
               , pch = 21)
        lines(na.exclude(tempTreat[r,]) ~ xlabBlank[!is.na(tempTreat[r,])] 
              , col = colList[treat]
              , lty = "dotted"
              , lwd = 2
        )
      }
      
    }
  }
  
  par(fig = c(0.7,1,0,1), mar = c(4,0,4,0), new = TRUE)
  plot(0, NULL
       , pch = ""
       , xaxt = 'n'
       , yaxt = 'n'
       , xlab = ""
       , ylab = ""
       , bty = 'n'
  )
  
  legend( "left"
          , legend = c(treatmentList, "Overall")
          , lty = c(rep("dotted", length(treatmentList)), "solid")
          , col = c(colList, "black")
          , lwd = c(rep(2, length(treatmentList)),2)
          , pch = c(rep(21, length(treatmentList)),19)
  )
  dev.off()
}


