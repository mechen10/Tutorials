#!/bin/bash

# Random data of taxa summaries plot with 6 samples and 2 species
marine1 <- c(.2,.8)
marine2 <- c(.1,.9)
marine3 <- c(.3,.7)

fresh1 <- c(.9,.1)
fresh2 <- c(.88,.12)
fresh3 <- c(.7,.3)

# Make matrix
taxaTable <-cbind(marine1,marine2,marine2,fresh1,fresh2,fresh3)

# Name OTUs OTU1 and OTU2
rownames(taxaTable) <- c('OTU1','OTU2')

# Set number of colors you want (number of taxa or OTUs)
NUMBEROFCOLORS <- nrow(taxaTable)

# Choose random colors
randomColors <- sample(colors()[-grep("gr(a|e)y|white", colors())], NUMBEROFCOLORS)

# Make barplot where the fresh and marine samples are separated. Make legend on the side, separate from barplot

quartz()
par(fig = c(0,0.7,0,1) # This puts the first plot within the x values of 0-0.7 and the y values of 0-1
    )
barplot(taxaTable
        , col = randomColors[factor(rownames(taxaTable))] # Colors from above; 
        #the format for coloring is colorList[factor(listofIDs)] because this ensures the colors match up to certain IDs.
        # You COULD just do 'col = randomColors' and this will just iterate through all the colors each time, but if you have missing data or uneven samples, this way is better
        , las = 2 # Make axis vertical
        , space = c(0,0,0,1,0,0) # This sets the spacing BEFORE each bar; note there is more space between 3rd and 4th bars
        )
par(fig = c(0.7,1,0,1) # This puts second plot (where legend will go) in remaining space
    , new = TRUE # Tells it to plot on same plot as barplot
    )
plot(0,0 # MAKE A BLANK PLOT
     , pch = '' # No points
     , xaxt = 'n' # Supress x axis
     , yaxt = 'n' # Supress y axis
     , bty = 'n' # Supress box around plot
     , xlab = '' # no xlable
     , ylab = '' # no ylable
     )
legend('center' # puts in center of space
       , legend = rownames(taxaTable) # OTUs
       , pch = 21
       , col = 'black' # outline is black for pch 19
       , pt.bg = randomColors # fills in pch 19 with color
       )
