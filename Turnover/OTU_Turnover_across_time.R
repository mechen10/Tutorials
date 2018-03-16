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

######################### Set up data ############################

# This data has TWO treatment types across time: Hakai and Port Moody (aka Reed Point)
# It calculates and plots turn over for BOTH sites in a single plot. Adjust as necessary for your own data.

# Filter out all samples that are NOT AM swabs (water, other surfaces, etc)
MF <- MF[grep("FB|CR|BL", MF$Morph),]
OTUTable <- OTUTable[,match(rownames(MF), colnames(OTUTable))]

# Split up two datasets; Hakai and Port Moody
# HAKAI
# Filter out all Hakai samples
MF.H <- MF[which(MF$Type == "H"),]
# Filter OTU table based on site names in MF.H
OTUTable.H <- OTUTable[,match(rownames(MF.H), colnames(OTUTable))]

# PORTMOODY
# Filter out all Hakai samples
MF.P<- MF[which(MF$Type == "P"),]
# Filter OTU table based on site names in MF.H
OTUTable.P <- OTUTable[,match(rownames(MF.P), colnames(OTUTable))]

# Then, filter tables so that OTUs that are not present in at least 'MinOccur' samples are deleted.
# HAKAI
toDelete <- c()
for (r in 1:nrow(OTUTable.H)) {
  nOccurance <- sum(OTUTable.H[r,] > 0.0)
  if (nOccurance < MinOccur) {
    toDelete <- c(toDelete, r)
  }
}
if (length(toDelete) > 0) {
  OTUTable.H.filt <- OTUTable.H[-toDelete,]
}

# PortMoody
toDelete <- c()
for (r in 1:nrow(OTUTable.P)) {
  nOccurance <- sum(OTUTable.P[r,] > 0.0)
  if (nOccurance < MinOccur) {
    toDelete <- c(toDelete, r)
  }
}
if (length(toDelete) > 0) {
  OTUTable.P.filt <- OTUTable.P[-toDelete,]
}

# Make presence/absence
# HAKAI
OTUTable.H.presabs <- OTUTable.H.filt
for (c in 1:ncol(OTUTable.H.filt)) {
  for (r in 1:nrow(OTUTable.H.filt)) {
    if (OTUTable.H.filt[r,c] > 0.0) {
      OTUTable.H.presabs[r,c] <- 1
    } else {
      OTUTable.H.presabs[r,c] <- 0
    }
  }
}

# PORTMOODY
OTUTable.P.presabs <- OTUTable.P.filt
for (c in 1:ncol(OTUTable.P.filt)) {
  for (r in 1:nrow(OTUTable.P.filt)) {
    if (OTUTable.P.filt[r,c] > 0.0) {
      OTUTable.P.presabs[r,c] <- 1
    } else {
      OTUTable.P.presabs[r,c] <- 0
    }
  }
}

########### Finding presence/absence for turnover #############


#HAKAI

# Order factors into levels
MF.H$Time <- factor(MF.H$Time, levels = c("20","60","360","720","5760")) 

# Make into list
timesList.H <- levels(MF.H$Time)

# Empty matrix with 4 types of OTUs; total, lost, gained, and retained
turnOver.H <- matrix(nrow = 4, ncol = length(timesList.H))
rownames(turnOver.H) <- c("total","lost","gained", "retained")
colnames(turnOver.H) <- timesList.H

allold <- c()

for (i in 1:length(timesList.H)) {
  ttemp <- timesList.H[i]
  allonly <- OTUTable.H.presabs[,grep(paste0("(FB|BL|CR)[-]", ttemp,"[-]"), colnames(OTUTable.H.presabs))]
  
  # Filter by which ones are actually there

  allonly <- allonly[which(rowSums(allonly) != 0),]
  
  # Get total OTU count for this time
  allcount <- nrow(allonly)
  
  # Get "old" list of OTUs and compare which ones are the same
  alllost <- length(allold)-sum(allold %in% rownames(allonly))
  
  allgain <- allcount-sum(rownames(allonly) %in% allold)
  
  # Get "retained" OTUs between each
  
  allretain <- sum(rownames(allonly) %in% allold)
  
  
  # Load into matrices
  turnOver.H["total", paste0(ttemp)] <- allcount
  turnOver.H["lost", paste0(ttemp)] <- round(alllost/allcount,2)
  turnOver.H["gained", paste0(ttemp)] <- round(allgain/allcount,2)
  turnOver.H["retained", paste0(ttemp)] <- round(allretain/allcount,2)
  
  # Finally, set "old" as the new ones
  allold <- names(which(rowSums(allonly) != 0))
}

#HAKAI

# Order factors into levels
MF.P$Time <- factor(MF.P$Time, levels = c("20","60","180","360","720","1440")) 

# Make into list
timesList.P <- levels(MF.P$Time)

# Empty matrix with 4 types of OTUs; total, lost, gained, and retained
turnOver.P <- matrix(nrow = 4, ncol = length(timesList.P))
rownames(turnOver.P) <- c("total","lost","gained", "retained")
colnames(turnOver.P) <- timesList.P

allold <- c()

for (i in 1:length(timesList.P)) {
  ttemp <- timesList.P[i]
  allonly <- OTUTable.P.presabs[,grep(paste0(ttemp,"[.](FB|BL|CR)[.]"), colnames(OTUTable.P.presabs))]
  
  # Filter by which ones are actually there
  allonly <- allonly[which(rowSums(allonly) != 0),]
  
  # Get total OTU count for this time
  allcount <- nrow(allonly)
  
  # Get "old" list of OTUs and compare which ones are the same
  alllost <- length(allold)-sum(allold %in% rownames(allonly))
  
  allgain <- allcount-sum(rownames(allonly) %in% allold)
  
  # Get "retained" OTUs between each
  allretain <- sum(rownames(allonly) %in% allold)
  
  
  # Load into matrices
  
  turnOver.P["total", paste0(ttemp)] <- allcount
  turnOver.P["lost", paste0(ttemp)] <- round(alllost/allcount,2)
  turnOver.P["gained", paste0(ttemp)] <- round(allgain/allcount,2)
  turnOver.P["retained", paste0(ttemp)] <- round(allretain/allcount,2)
  
  # Finally, set "old" as the new ones
  allold <- names(which(rowSums(allonly) != 0))
}

# Print out two turn Over files in case you want to do other things with them
write.table(turnOver.H
            , sep = "\t"
            , file = "turnOver_H.txt"
            , col.names = NA # This adds space to first column so that headers line up
            )
write.table(turnOver.P
            , sep = "\t"
            , file = "turnOver_P.txt"
            , col.names = NA # This adds space to first column so that headers line up
)

# Make the two turnOver files into a single table of retained OTUs
retained <- matrix(nrow = 2, ncol = 6)
# New rownames
rownames(retained) <- c("H","P")
# ALL levels 
colnames(retained) <- c("60","180","360","720","1440","5760")
  for (t in c("60","180","360","720","1440")) {
    TempP <- turnOver.P["retained",t]
    retained["P", t] <- TempP
  }
  for (t in c("60","360","720","5760")) {
    TempH <- turnOver.H["retained",t]
    retained["H", t] <- TempH
  }

# Make the two turnOver files into a single table of lost OTUs
lost <- matrix(nrow = 2, ncol = 6)
# New rownames
rownames(lost) <- c("H","P")
# ALL levels 
colnames(lost) <- c("60","180","360","720","1440","5760")
for (t in c("60","180","360","720","1440")) {
  TempP <- turnOver.P["lost",t]
  lost["P", t] <- TempP
}
for (t in c("60","360","720","5760")) {
  TempH <- turnOver.H["lost",t]
  lost["H", t] <- TempH
}

# Make the two turnOver files into a single table of gained OTUs
gained <- matrix(nrow = 2, ncol = 6)
# New rownames
rownames(gained) <- c("H","P")
# ALL levels 
colnames(gained) <- c("60","180","360","720","1440","5760")
for (t in c("60","180","360","720","1440")) {
  TempP <- turnOver.P["gained",t]
  gained["P", t] <- TempP
}
for (t in c("60","360","720","5760")) {
  TempH <- turnOver.H["gained",t]
  gained["H", t] <- TempH
}



####### PLOTTING ##########
# PLOT RETAINED OTUs (OTUs that are consistent between current time point and previous)
miny <- min(retained, na.rm = TRUE)
maxy <- max(retained, na.rm = TRUE)
pdf("retainedOTUs.pdf", pointsize = 14, width = 7, height = 5)
par(fig = c(0,0.8,0,1), mar = c(6,4,2,4))
plot(c(1,2,3,4,5,6),NULL
     , xlim = c(1,6)
     , ylim = c(miny,maxy)
     , pch= ""
     , xaxt = "n"
     , xlab = ""
     , ylab = "Percent community retained"
)
title(xlab = "Time"
      , line = 4)
axis(side = 1
     , at = c(1,2,3,4,5,6)
     , labels = c("1 hour","3 hours","6 hours","12 hours","1 day","4 days")
     , las = 2)
points(retained[1,] ~ c(1,2,3,4,5,6)
       , col = "dodgerblue"
       , pch = 19)
points(retained[2,] ~ c(1,2,3,4,5,6)
       , col = "dodgerblue"
       , pch = 21)

lines(na.exclude(retained[1,]) ~ c(1,2,3,4,5,6)[!is.na(retained[1,])] 
      , col = "dodgerblue"
      , lty = "solid"
      , lwd = 2
)

lines(na.exclude(retained[2,]) ~ c(1,2,3,4,5,6)[!is.na(retained[2,])] 
      , col = "dodgerblue"
      , lty = "dotted"
      , lwd = 2
)

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
        , legend = c("Reed Point","Hakai")
        , lty = c("dotted","solid")
        , col = "dodgerblue"
        , lwd = c(2,2)
        , pch = c(21,19)
)
dev.off()


# PLOT LOST OTUs (OTUs that are lost between previous time point and current time point)
miny <- min(lost, na.rm = TRUE)
maxy <- max(lost, na.rm = TRUE)
pdf("lostOTUs.pdf", pointsize = 14, width = 7, height = 5)
par(fig = c(0,0.8,0,1), mar = c(6,4,2,4))
plot(c(1,2,3,4,5,6),NULL
     , xlim = c(1,6)
     , ylim = c(miny,maxy)
     , pch= ""
     , xaxt = "n"
     , xlab = ""
     , ylab = "Percent community lost"
)
title(xlab = "Time"
      , line = 4)
axis(side = 1
     , at = c(1,2,3,4,5,6)
     , labels = c("1 hour","3 hours","6 hours","12 hours","1 day","4 days")
     , las = 2)
points(lost[1,] ~ c(1,2,3,4,5,6)
       , col = "dodgerblue"
       , pch = 19)
points(lost[2,] ~ c(1,2,3,4,5,6)
       , col = "dodgerblue"
       , pch = 21)

lines(na.exclude(lost[1,]) ~ c(1,2,3,4,5,6)[!is.na(lost[1,])] 
      , col = "dodgerblue"
      , lty = "solid"
      , lwd = 2
)

lines(na.exclude(lost[2,]) ~ c(1,2,3,4,5,6)[!is.na(lost[2,])] 
      , col = "dodgerblue"
      , lty = "dotted"
      , lwd = 2
)

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
        , legend = c("Reed Point","Hakai")
        , lty = c("dotted","solid")
        , col = "dodgerblue"
        , lwd = c(2,2)
        , pch = c(21,19)
)
dev.off()



# PLOT GAINED OTUs (OTUs that are gained between previous time point and current time point)
miny <- min(gained, na.rm = TRUE)
maxy <- max(gained, na.rm = TRUE) 
pdf("gainedOTUs.pdf", pointsize = 14, width = 7, height = 5)
par(fig = c(0,0.8,0,1), mar = c(6,4,2,4))
plot(c(1,2,3,4,5,6),NULL
     , xlim = c(1,6)
     , ylim = c(miny,maxy)
     , pch= ""
     , xaxt = "n"
     , xlab = ""
     , ylab = "Percent community gained"
)
title(xlab = "Time"
      , line = 4)
axis(side = 1
     , at = c(1,2,3,4,5,6)
     , labels = c("1 hour","3 hours","6 hours","12 hours","1 day","4 days")
     , las = 2)
points(gained[1,] ~ c(1,2,3,4,5,6)
       , col = "salmon"
       , pch = 19)
points(gained[2,] ~ c(1,2,3,4,5,6)
       , col = "salmon"
       , pch = 21)

lines(na.exclude(gained[1,]) ~ c(1,2,3,4,5,6)[!is.na(gained[1,])] 
      , col = "salmon"
      , lty = "solid"
      , lwd = 2
)

lines(na.exclude(gained[2,]) ~ c(1,2,3,4,5,6)[!is.na(gained[2,])] 
      , col = "salmon"
      , lty = "dotted"
      , lwd = 2
)

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
        , legend = c("Reed Point","Hakai")
        , lty = c("dotted","solid")
        , col = "salmon"
        , lwd = c(2,2)
        , pch = c(21,19)
)
dev.off()

