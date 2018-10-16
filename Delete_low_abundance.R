##### THIS IS NOT AN EXECUTABLE SCRIPT. FOLLOW INSTRUCTIONS MANUALLY #####


# ~~~~~~~~~~~~~~~~~~~ Deleting low abundance OTUs in OTU table ~~~~~~~~~~~~~~~~~~~ #
# Requirements:
# OTU table must be output of 'biom convert' command.  
#  (AKA, columns are samples; rows are OTUs; the first line is '#Constructed from biom file')
# May or may not have taxonomy, but if it has it, the header must be called "taxonomy"

###### LOAD OTU TABLE and set thresholds #########

# Choose OTU table
otu.table <- read.delim(file.choose(), header=TRUE, row.names = 1, stringsAsFactors = FALSE, skip=1)

# get rid of "taxonomy" column, if applicable
taxcol <- which(colnames(otu.table)=="taxonomy")
if (length(taxcol) >0 ) {
    taxa <- data.frame(otu.table[,taxcol], row.names = rownames(otu.table))
    otu.table <- otu.table[,-taxcol]
}

# figure out if table is relative abundance or count
if (max(otu.table) < 1) {
    count <- FALSE
} else {
    count <- TRUE
}

#### Set per sample thresholds ####
## You should set ONE of these, depending on whether your table is relative abundance or threshold.
samp_thresh_count <- 5 # number of OTUs per sample to keep in dataset
samp_thresh_percent <- 0.01 # relative abundance of otu in each sample to keep in dataset
#### Set overall dataset threshold ####
## You should set ONE of these, depending on whether your table is relative abundance or threshold.
overall_thresh_count <- 100 # minimum number of times an OTU appears in OTU table
overall_thresh_percent <- 0.001 # minimum relative abundance an OTU appears in OTU table


##### Convert OTU table to relative abundance, if desired #######
## If your OTU table is counts but you want to convert to relative abundance, use this:
# set want_relabund to TRUE if you want relative abundance and NOT count data
want_relabund <- FALSE
if (want_relabund) {
    otu.table <- t(t(otu.table)/colSums(otu.table))
    count <- FALSE
}


###### (1) CHANGE ALL COUNTS WITH LESS THAN X PER SAMP TO 0 ######
print("DELETING COUNTS WITH LESS THAN 5 PER SAMP")

if ( count ) {
    otu.table[otu.table < samp_thresh_count] <- 0
} else if ( !count) {
    otu.table[otu.table < samp_thresh_percent] <- 0
}


###### (2) REMOVE ALL OTUS WITH LESS THAN X IN WHOLE OTU TABLE ####
# remove all zeros from otu table
print("DELETEING COUNTS WITH LESS THAN 10 per SAMP")


if ( count ) {
    todel <- (which(rowSums(otu.table) < overall_thresh_count))
    if (length(todel) > 0) {
        otu.table <- otu.table[-todel,]
    }
} else if ( !count) {
    todel <- (which(rowSums(otu.table) < overall_thresh_percent))
    if (length(todel) > 0) {
        otu.table <- otu.table[-todel,]
    }
}

####### Double check that the smallest observation is your threshold ########
# Observations per sample minimum
min(otu.table[otu.table >0])
# Observations for otu table minimum
min(rowSums(otu.table))

###### Print resulting OTU table as txt file ######
# Add back taxonomy, if applicable
if ( length(taxcol) > 0 ) {
    # taxa.filt <- taxa[rownames(otu.table),]
    otu.table$taxonomy <- taxa[rownames(otu.table),]
}

## You'll have to manually convert back into biom if desired.
# inserts first column header
otu.temp <- data.frame("#OTU ID"=rownames(otu.table),otu.table)
colnames(otu.temp)[1] <- "#OTU ID"
# write out
write.table(otu.temp, file="OTUTable_lowabundremoved.txt", quote = F, sep="\t", row.names = FALSE )

