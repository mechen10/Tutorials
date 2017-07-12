#!/bin/bash

# This is a tutorial to make barplots by time and treatment type
# I also go over how to randomly assign colors.

######################### Data variables ############################

# Set path. Change this path to folder that you want output to be generated
setwd("/Users/melissachen/Documents/Masters/UNIVERSALCODE_git/Tutorials/")

# This is the file path to your taxa summaries table (collapsed at some level; I've used Order here).
# Generated using summarize_taxa.py OR summarize_taxa_through_plots.py
OTUFP <- "./OTU_L4.txt"

# This is the file path to your Mapping file.
mfFP <- "./MF_AM.txt"

######################### Load data ############################

OTUTable <- read.delim(paste0(OTUFP)
                       , header = TRUE
                       , skip = 1
                       , row.names = 1)

MF <- read.delim(paste0(mfFP)
                 , header = TRUE
                 , row.names = 1
                 , stringsAsFactors = FALSE
                 )

# Make sure all sites are in OTU and MF; and vice versa

OTUTable <- OTUTable[,colnames(OTUTable) %in% rownames(MF)]
MF <- MF[rownames(MF) %in% colnames(OTUTable),]

######################### Set up data ############################

# In my data, I want to plot time series (Time) and treatment (Morph).
# First and foremost, I get rid of all data that is NOT one of three treatment types (FB,BL,CR)

MF <- MF[grep("FB|BL|CR", MF$Morph),]

# Now, make sure all samples are found in both MF and OTU Table
OTUTable <- OTUTable[,which(colnames(OTUTable) %in% rownames(MF))]
MF <- MF[which(rownames(MF) %in% colnames(OTUTable)),]

# And make sure it's in correct order.
OTUTable <- OTUTable[,match(rownames(MF), colnames(OTUTable))]

# Sort treatments into the order you want
MF$Morph <- factor(MF$Morph
                   , levels = c("FB","BL","CR") # This list needs to include ALL treatments and it is in the order you want them to be
                   )

# Sort time into order you want. Make it a 'factor', not a numerical variable

MF$Time <- factor(MF$Time
                  , levels = c("20","60","360","720","5760") # This list needs to include ALL treatments and it is in the order you want them to be
                  )

# Sort entire MF by Time and then by replicate
MF <- MF[with(MF, order(MF$Time, MF$Rep)),]

# Now, separate OTU table into treatment groups

MF.FB <- MF[grep("^FB$", MF$Morph),] # Get FB from MF
OTUTable.FB <- OTUTable[,c(match(rownames(MF.FB), colnames(OTUTable)))] # Filter OTU table by MF

MF.BL <- MF[grep("^BL$", MF$Morph),] # Get FB from MF
OTUTable.BL <- OTUTable[,c(match(rownames(MF.BL), colnames(OTUTable)))] # Filter OTU table by MF

MF.CR <- MF[grep("^CR$", MF$Morph),] # Get FB from MF
OTUTable.CR <- OTUTable[,c(match(rownames(MF.CR), colnames(OTUTable)))] # Filter OTU table by MF

######################### Colors ############################

# Get all colors, but grep out greys and whites
colors.filtered <- colors()[-grep("white|gr(a|e)y", colors())]

# Sample random colors for each OTU in OTU table
colorsRandom <- sample(colors.filtered, nrow(OTUTable))
# NOTE: There may not be enough colours if you have very diverse samples, or if you use high levels of taxonomy (ie genus)

######################### Set up plot parameters ############################

# We have 3 treatments at 5 time points.

# Find out how many reps in each treatment; to use in next step
counts.FB <- table(MF.FB$Time)
counts.BL <- table(MF.BL$Time)
counts.CR <- table(MF.CR$Time)

# Find maximum number of treatments at each time point
MAXcounts <- apply(rbind(counts.FB, counts.BL, counts.CR), MARGIN = 2, FUN = max)

# Need to adjust number of bars in each treatment so that they all equal MAXcounts
# Adjustments below:

# FB needs to add:
needToAdd.FB <- MAXcounts-counts.FB
# Only needs to add one at the end
OTUTable.FB.edited <- cbind(OTUTable.FB, matrix(ncol = 1, nrow = nrow(OTUTable.FB)))

# BL needs to add 
needToAdd.BL <- MAXcounts-counts.BL
# Add 2 at second time point
OTUTable.BL.edited <- cbind(OTUTable.BL[,1:18], matrix(ncol = 2, nrow = nrow(OTUTable.BL)), OTUTable.BL[,19:ncol(OTUTable.BL)])

# CR needs to add 
needToAdd.CR <- MAXcounts-counts.CR
# Add 1 at second time point; 1 at 3rd time point; 2 at 4th time point; 3 at last time point
OTUTable.CR.edited <- cbind(OTUTable.CR[,1:19]
                            , matrix(ncol = 1, nrow = nrow(OTUTable.CR))
                            , OTUTable.CR[,20:28]
                            , matrix(ncol = 1, nrow = nrow(OTUTable.CR))
                            , OTUTable.CR[,28:34]
                            , matrix(ncol = 2, nrow = nrow(OTUTable.CR))
                            ,OTUTable.CR[,35:ncol(OTUTable.CR)]
                            , matrix(ncol = 2, nrow = nrow(OTUTable.CR)))

# Vector that dictates spacing; spacing BEFORE each bar
spacingVector <- c(rep(0,10),1,rep(0,9),1,rep(0,9),1,rep(0,8),1,rep(0,10))

pdf("BarPlot.pdf" # Use pdf instead of jpeg because you can adjust size without worrying about resolution.
    , pointsize = 14) # Make font size bigger
par(fig = c(0,1,0,1)
    )
plot(0,0 # PLOT EMPTY PLOT (necessary for 'title' and 'axis')
    , xlab = "" # No xlabel
    , ylab = ""# No ylabel 
    , xaxt = "n" # no x axis
    , yaxt = "n" # no y axis
    , bty = "n" # not box outline
    , pch = "" # no symbol to plot
    )
axis(1, # bottom axis
     , at = c(-0.65, -0.35, 0, 0.3, 0.65) # This is determined just by trying diff values until it looks right
     , labels = c("20 min","1 h","6 h","12 h","4 d")
     , tick = FALSE # suppress the tick line 
     , line = 2 # move down 3 lines
     )
title(ylab = "Relative Abudance"
      , line = 2 # Move out 2 lines
      )
par(fig = c(0,1,0.6,1) # c(x1,x2,y1,y2); This sets up fraction of page this plot will take up.
    , mar = c(2.1,5.1,2.1,4.1) # Sets up margin size. Order is bottom, left, top, right.
    , new = TRUE
    , oma = c(2,2,0,0)
    )
barplot(as.matrix(OTUTable.FB.edited) # Must be matrix
        , ylab = "Finely Branched"
        , col = colorsRandom
        , space = spacingVector # Spacing; space is BEFORE bar
        , xaxt = 'n' # Turn off x axis, because I will add in manual one at bottom.
        , xlab = '' # Blank x axis label; will add manula one at bottom
        )
par(fig = c(0,1,0.3,0.7) # c(x1,x2,y1,y2); This sets up fraction of page this plot will take up.
    , mar = c(2.1,5.1,2.1,4.1) # Sets up margin size. Order is bottom, left, top, right.
    , new = TRUE
    , oma = c(2,2,0,0)
)
barplot(as.matrix(OTUTable.BL.edited) 
        , ylab = "Bladed"
        , col = colorsRandom
        , space = spacingVector
        , xaxt = 'n' # Turn off x axis, because I will add in manual one at bottom.
        , xlab = '' # Blank x axis label; will add manula one at bottom
)
par(fig = c(0,1,0,0.4) # c(x1,x2,y1,y2); This sets up fraction of page this plot will take up.
    , mar = c(2.1,5.1,2.1,4.1) # Sets up margin size. Order is bottom, left, top, right.
    , new = TRUE
    , oma = c(2,2,0,0)
    )
barplot(as.matrix(OTUTable.CR.edited)
        , ylab = "Crustose"
        , col = colorsRandom
        , space = spacingVector
        , xaxt = 'n' # Turn off x axis, because I will add in manual one at bottom.
        , xlab = '' # Blank x axis label; will add manula one at bottom
)
dev.off() # Have at end of plot to tell it you're done drawing
