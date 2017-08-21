# workflow for formatting sample names according to metadata file (Quadrat_Metagenomics_metadata.txt)

# load in directory 
path.to.md <- "/Users/parfreylab/Desktop/lab_member_files/kevin/metagenomics/matt_plate_samples/"
dir <- "/Users/parfreylab/Desktop/lab_member_files/kevin/metagenomics/matt_plate_samples/KEGG_RUN_2/res_kegg_levels/for_plotting/2/"
# load in files to format
files <- c("kegg_level_2_1rm.txt")
# load in metadata file
md <- paste(path.to.md, "Quadrat_Metagenomics_metadata.txt", sep = "")
md <- read.delim(md, sep = "\t", header = TRUE, check.names = FALSE) # check.names removes replacement of separaters using "."
md$`#SampleID` # names of all the samples. already in proper order
sample_order <- c("KEGG Function", as.vector(md$`#SampleID`)) # get the sample names as a vector to reorder columns. add in Protein Family so we know the protein families

########################################################################
# USE IF LOOPING THROUGH MULTIPLE LINES (CMD + SHIFT + C TO UNCOMMENT) #
########################################################################
# For all files in the directory
# for (f in files) {
#   stat <- paste(dir, "res/z_summary/", f, sep = "")
#   stat <- read.delim(stat, sep = "\t", header = TRUE, check.names = FALSE) 
#   # gsub to parse out patterns in the names of the file (columns), replace as necessary
#   # names(stat) is all sample names
#   names(stat) <- gsub("_L001_R1_001_merged_cat_filtered", "", names(stat))
#   names(stat) <- gsub("-", "_", names(stat))
#   names(stat) <- gsub("_S[0-9]+", "", x = names(stat))
#   stat <- stat[,sample_order] # reorder columns based on sample names from metadata file
#   write.table(stat, file = paste(dir, "res/z_summary_formatted/", f, sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t" )
# }

########################
# USE IF ONLY ONE FILE #
########################
stat <- paste(dir, "kegg_level_2.txt", sep = "")
stat <- read.delim(stat, sep = "\t", header = TRUE, check.names = FALSE, na.strings = c("NA"))
names(stat) <- gsub("_L001_R1_001_merged_cat_filtered", "", names(stat))
names(stat) <- gsub("-", "_", names(stat))
names(stat) <- gsub("_S[0-9]+", "", names(stat))
names(stat)
stat <- stat[,sample_order] # reorder columns based on sample names from metadata file
head(stat)
write.table(stat, file = paste(dir, "kegg_level_2_ordered.txt", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t" )


# if needed - deletes blank columns
stat[52] <- NULL
stat[52] <- NULL

stat
names(stat) <- gsub(".txt", "", names(stat))
stat <- stat[,sample_order]
stat$KO <- seq(nrow(stat))
stat
write.table(stat, file = paste(dir, "KEGG_RUN_2/KEGG_stats_4_for_qiime.txt", sep = ""), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t" )

# for later use
names(stat)
names(stat) <- gsub("quadrat", "qd", names(stat))
names(stat) <- gsub("water", "w", names(stat))
names(stat) <- gsub("pre.treatment", "pt", names(stat))

