library(vegan)

#### FUNCTION ####
pairwise_adonis_MC <- function(dm_data, mf, mf_col, p_corr) {
    # This function creates a list of items: first, it generates pairwise adonis for each pair of treatments to compare.
    # Then, it calculates the adjusted p-value as designated by p_corr
    # It need vegan to work
    
    # Check to make sure rownames and colnames are the same in your distance matrix; sometimes R loads them funny
    if ( any(!(rownames(dm_data) == colnames(dm_data))) ) { # check to see if any row != col name
        print("WARNING: dm not symmetrical. Check to make sure correct file is loaded. This script automatically copies row names and uses as col names")
        colnames(dm_data) <- rownames(dm_data)
    }
    
    # Note: you MUST ensure that all samples in dm are also in mf. This script does not check this for you.
    groups_compare <- levels(factor(mf[,mf_col])) # extract the names of the groups you want to compare
    
    combin_groups <- combn(x=length(groups_compare), m=2) # Make matrix of comparisons to make
    adonis_results <- list()
    # Iterate through all combinations of groups to compare
    for ( comb in 1:ncol(combin_groups) ) {
        grp1 <- combin_groups[,comb][1] # group one
        grp2 <- combin_groups[,comb][2] # group two
        
        grp1_names <- rownames(mf)[which(mf[,mf_col] == groups_compare[grp1])] # list of sample names, group1
        grp2_names <- rownames(mf)[which(mf[,mf_col] == groups_compare[grp2])] # list of sample names, group2
        
        mf_compare <- mf[c(grp1_names,grp2_names),] # make filtered mf with only samples to compare
        dm_compare <- dm_data[c(grp1_names,grp2_names),c(grp1_names,grp2_names)] # make filtered dm with only samples to compare
        adonis_results[[paste0(groups_compare[grp1],".",groups_compare[grp2])]] <- adonis(dist(dm_compare) ~ mf_compare[,mf_col]) # Run adonis
    }
    
    # Name make summary with p-value adjustments
    adj_pvalues <- matrix(ncol=length(adonis_results), nrow=1, dimnames = list(paste0(p_corr),names(adonis_results))) # Empty matrix 
    for ( i in names(adonis_results) ) {
        p <- adonis_results[[i]]$aov.tab$`Pr(>F)`[1] # extract p-value
        adj_pvalues[1,i] <- p.adjust(p, method=p_corr, n = length(adonis_results)) # adjust and save
    }
    
    adonis_results[["Adjusted_pvalues"]] <- adj_pvalues
    
    return(adonis_results) # Return results
}

######### USAGE #############

### ADJUST FOR YOURSELF ###
setwd("Documents/Masters/UNIVERSALCODE_git/Tutorials/pairwise_adonis/") # Set working to your current working that has distance matrix and mapping file 

dm_data <- read.delim("bray_curtis_dm.txt",header=TRUE,row.names = 1) # Load distance matrix

mf <- read.delim("MF_nochlpmito_m800.txt",header=TRUE,row.names=1,na.strings=c("NA","na",""),stringsAsFactors=FALSE)

mf_col <- "SubstrateType" # This is the name of the mf column that you want to do a pairwise anova on

p_corr <- "BH" # The adjustment method to use for p-value adjustments; look in ?p.adjust for more options
###

# Examples:
results_benjaminihochberg <- pairwise_adonis_MC(dm_data=dm_data, mf = mf, mf_col = mf_col, p_corr = p_corr)

results_bonferroni <- pairwise_adonis_MC(dm_data=dm_data, mf = mf, mf_col = mf_col, p_corr = "bonferroni")


###############################

