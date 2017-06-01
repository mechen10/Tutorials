#!/bin/bash

###########################
# INSTRUCTIONS:

# 1. Put the file path of your OTU table in the "INPUT" section
# 2. Save this document
# 3. In terminal, move to your desired folder (where you want the output to be)
# 4. Run this script by dragging it into terminal and pressing enter (Make sure macqiime is working!!)
# 5. Done

###########################
# INPUT: PUT YOUR FILES HERE
# Follow the format provided

# OTUTABLE=/path/to/your/OTUtable.biom
OTUTABLE=

###########################

mkdir FILTERED_OTU_TABLE

# Converts the biom file into a text file
biom convert -i $OTUTABLE -o ./FILTERED_OTU_TABLE/OTUTABLE_original.txt --to-tsv --header-key taxonomy

# Reverse greps out Chloroplasts, Mitochondria, and Eukaryota. 
# Note: the output is OTUs to KEEP, not the discard; done this way bc it is more efficient.
grep -v "Chloroplast" ./FILTERED_OTU_TABLE/OTUTABLE_original.txt | grep -v "Mitochondria" | grep -v "Eukaryota" >> ./FILTERED_OTU_TABLE/nochlpmitoeuk.txt

# Filter them out
filter_otus_from_otu_table.py -i $OTUTABLE -e ./FILTERED_OTU_TABLE/nochlpmitoeuk.txt --negate_ids_to_exclude -o ./FILTERED_OTU_TABLE/OTU_Table_nochlpmitoeuk.biom

rm ./FILTERED_OTU_TABLE/nochlpmitoeuk.txt
rm ./FILTERED_OTU_TABLE/OTUTABLE_original.txt

echo ""
echo 'Done filtering'
echo ""