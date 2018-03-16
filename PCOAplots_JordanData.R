#!/bin/bash

# This is to do a pcoa plot with Jordan's Macrocystis paper

################ SET PATHWAYS HERE: ############
# Change pathway to where ever your files are in your computer
# Also, change pathway to folder where you want this to save images

# SET WORKING DIRECTORY:
setwd("/Users/melissachen/Documents/Masters/Other_peoples_stuff/JordanData")

# DM FILEPATH
dmFP <- "./bray_curtis_rarefaction_400_CAZy.txt"

# MF FILEPATH
mfFP <- "./cazy_sample.mapping.txt"

#################################################
############### REQUIRED LIBRARIES #####################
require("ape") # ape has pcoa
require("MASS") # MASS has isoMDS

############### IMPORT DATA #####################

# Import distance matrix
dm <- read.delim(paste0(dmFP) # pastes the filepath set above
                 , header = TRUE # first line should be treated as headers
                 , row.names = 1 # first column should be treated as rownames
                 , stringsAsFactors = FALSE # Prevents R from interpreting numbers as 'factors'-- not crucial, but may prevent issues later
                 )

# Import mappingfile

MF <- read.delim(paste0(mfFP)
                 , header = TRUE # pastes the filepath set above
                 , row.names = 1 # first line should be treated as headers
                 , stringsAsFactors = FALSE # first column should be treated as rownames
                 , na.strings = c("","na","NA","N/A") # tells R that things that should be interpreted as NA
                 )
############### DATA EDITING #####################

# First, make sure all samples found in MF are in dm and vice versa
dm <- dm[rownames(dm) %in% rownames(MF), colnames(dm) %in% rownames(MF)]
MF <- MF[rownames(MF) %in% rownames(dm),]

# Now, make sure they're in correct order

dm <- dm[match(rownames(MF),rownames(dm)),match(rownames(MF),colnames(dm))]

# Check that this is true

all(rownames(dm) == rownames(MF))
all(colnames(dm) == rownames(MF))
# Should equal "TRUE"

############### MAKE COLORS AND SYMBOLS #####################

# Make 'Description' a factor; again, I type out order that I want.
MF$Samples <- factor(MF$Samples, levels = c("Water","Macrocystis-Middle","Macrocystis-Bottom"))

### COLORS ####
samplePch <- c(21 #circle with outline
               , 24 #triangle with outline
               , 24 #trianlge with outline
)
# For more shapes, see ?pch

# For Samples, make them different colors
sampleCol <- c("blue","darkgreen","green")


############### MAKE PCOA #####################

# Make PCOA
dm.pcoa <- pcoa(dm)

# CHOOSE AXIS THAT YOU WANT
axis1 <- 1
axis2 <- 2
toPlot.PCOA <- dm.pcoa$vectors[,c(axis1,axis2)] # Extracts columns axis1 and axis2 (1 and 2, set above)

# GET PERCENT VARIATION OF AXIS
# $values$Relative_eig gives you variation of each axis

perc1 <- round(dm.pcoa$values$Relative_eig[axis1]*100,2)
perc2 <- round(dm.pcoa$values$Relative_eig[axis2]*100,2)
# Note, I *100 and also rounded them to 2 decimal places: round(x,2)


############### MAKE PLOT #####################
pdf("PCOA_plot.pdf", pointsize = 14) # I like to make font big, so I set pointsize to 14
par(mar= c(5.1,5.1,5.1,2.1)) # Use par to adjust margins so it's nice and pretty
plot(toPlot.PCOA
     , bg = sampleCol[factor(MF$Samples)] # Colors by factored levels
     , col = "black" # Outline is black
     , pch = samplePch[factor(MF$Samples)] # I like circles
     , cex = 2 # Makes points bigger
     
     , xlab = paste0("PCO 1 (", perc1 ,"% of variation)") # insert percent var from above
     , ylab = paste0("PCO 2 (", perc2 , "% of variation)")
     
     , main = "Macrocystis and water samples" # Main title
     )
legend("topright" # designate position
       , legend = levels(MF$Samples)
       , pch = samplePch
       , pt.bg = sampleCol
       , col = "black"
       )
dev.off() # Confirms that this is end of plot that it's making


############### MAKE NMDS #####################

# Make PCOA
dm.nmds <- isoMDS(as.dist(dm) # Needs to be distance matrix
                  , k = 2 # two axis
                  )

# EXTRACT AXIS
toPlot.NMDS <- dm.nmds$points # Extracts columns axis1 and axis2 (1 and 2, set above)

# EXTRACT STRESS
stressVal <- round(dm.nmds$stress/100, 2)

# Note, I /100 and also rounded them to 2 decimal places: round(x,2)


############### MAKE PLOT #####################
pdf("NMDS_plot.pdf", pointsize = 14) # I like to make font big, so I set pointsize to 14
par(mar= c(5.1,5.1,5.1,2.1)) # Use par to adjust margins so it's nice and pretty
plot(toPlot.NMDS
     , bg = sampleCol[factor(MF$Samples)] # Colors by factored levels
     , col = "black" # Outline is black
     , pch = samplePch[factor(MF$Samples)] # I like circles
     , cex = 2 # Makes points bigger
     
     , xlab = paste0("NMDS 1") # insert percent var from above
     , ylab = paste0("NMDS 2")
     
     , main = "Macrocystis and water samples" # Main title
     , sub = paste0("Stress: ", stressVal)
)
legend("topright" # designate position
       , legend = levels(MF$Samples)
       , pch = samplePch
       , pt.bg = sampleCol
       , col = "black"
)
dev.off() # Confirms that this is end of plot that it's making
