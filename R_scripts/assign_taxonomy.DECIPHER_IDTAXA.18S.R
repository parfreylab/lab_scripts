#to be done after processing amplicon data with dada2
#last modified: Feb 17th, 2026
#author: Evan Morien

####libraries####
library(tidyverse)
library(data.table)
library(DECIPHER)
library(Biostrings)
library(ensembleTax) #note, please install this from my fork until the devs fix their issues with PR2's updated taxonomy. use the devtools install method here, but replace "dcat4" with "morien" and set the vignettes parameter to FALSE: https://github.com/morien/ensembleTax

#code to install our fork of ensembleTax
#library(devtools)
#devtools::install_github("morien/ensembleTax", build_manual = TRUE, build_vignettes = FALSE, force=TRUE)
#packageVersion("ensembleTax")

####Inputs####
#sequence table from dada2 (see dada2 processing scripts for any amplicon). ASVs not labeled (i.e. raw sequences as ASV names)
seqtab <- fread("sequence_table.18S.merged.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
row.names(seqtab) <- seqtab[,1] #set row names
seqtab <- seqtab[,-1] #remove column with row names in it
seqtab <- as.matrix(seqtab) #cast the object as a matrix
mode(seqtab) <- "numeric"

## creating a DNAStringSet object from the ASVs
dna_18S <- DNAStringSet(getSequences(seqtab))

#there are three choices for reference database (as of Nov. 2025): PR2 (IDTAXA provided), SILVA (Hakai custom db), or MetaZooGene (Hakai custom DB)
#there is a section of this example code for each, read each carefully if you would like to use it

#### loading PR2 DECIPHER train set #this is provided by PR2 developers, see PR2's github page for updated versions ####
trainingset <- readRDS("/mnt/Genomics/Working/databases/amplicon_taxonomyDBs/DECIPHER_IDTAXA/pr2_version_5.1.0_SSU.decipher.trained.rds")

#taxonomic inference #threshold for accepting taxonomic assignments can be adjusted with the 'threshold' parameter. see documentation.
tax_info_18S <- IdTaxa(dna_18S, trainingset, type="extended", strand = "both", processors = 64, threshold = 50) 

## making and writing out standard output files:
# giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs_18S <- colnames(seqtab)

# this creates names (ASV1, ASV2, ..., ASVX) for each ASV
snam <- vector(mode = "character", length = length(dna_18S))
pad <- nchar(length(snam)) #length of zeroes to pad ASV names. we do this so the alphanumeric ASV names are sorted correctly by the idtax2df function
for (i in 1:length(dna_18S)) {
  snam[i] <- str_pad(as.character(i), width=pad, side="left", pad="0")
  snam[i] <- paste0("ASV", snam[i])
}
names(dna_18S) <- snam

tax_tab_18S <- idtax2df(tax_info_18S, 
                        db = "pr2", 
                        ranks = c("rootrank", "domain", "supergroup", "division", "subdivision", "class", "order", "family", "genus", "species"),
                        boot = 50,
                        rubric = dna_18S,
                        return.conf = FALSE)

colnames(tax_tab_18S) <- c("ASV_ID", "sequence", "domain", "supergroup", "division", "subdivision", "class", "order", "family", "genus", "species")
row.names(tax_tab_18S) <- NULL
#we added leading zeroes earlier to allow correct sorting order by the idtax2df function, now remove them so the ASV IDs match the other IDs in the rest of our pipeline
pattern <- "ASV0*"
replacement <- "ASV"
tax_tab_18S$ASV_ID <- gsub(pattern, replacement, tax_tab_18S$ASV_ID)
write.table(tax_tab_18S, "taxonomy_table.18S.DECIPHER.PR2v5.1.txt", sep = "\t", quote = F, row.names = FALSE)

# also create corresponding confidence table. this is useful for interpretation, but assignments can be filtered above by passing a threshold to the IdTaxa() function
ranks <- c("domain", "supergroup", "division", "subdivision", "class", "order", "family", "genus", "species")
confidences <- t(sapply(tax_info_18S, function (x) {
  m<-match(ranks, x$rank)
  confidence<-x$confidence[m]}
))
colnames(confidences) <- ranks
row.names(confidences) <- NULL
confidences <- data.frame("ASV_ID" = sub(">", "", names(dna_18S)), confidences, check.names = FALSE)
pattern <- "ASV0*"
replacement <- "ASV"
confidences$ASV_ID <- gsub(pattern, replacement, confidences$ASV_ID)
write.table(confidences, "taxonomy_table.18S.DECIPHER.PR2v5.1.conf.txt", sep = "\t", quote = F, row.names = FALSE)

#if not also assigning with SILVA138, now skip to section "formatting final output files"

#### slightly different procedure for SILVA v 138 ####
## loading SILVA v138 DECIPHER train set #see file silva_138_reformat_rescript.sh to see how these data were reformatted for this step from the raw silva files
trainingset <- readRDS("/mnt/Genomics/Working/databases/amplicon_taxonomyDBs/DECIPHER_IDTAXA/decipher_SILVA_2024/silva-138.1-ssu-nr99.rescript.corrected.euk.prok-clustered.8000prok.no-accn.trained.24-04-16.rds")

#taxonomic inference #threshold for accepting taxonomic assignments can be adjusted with the 'threshold' parameter. see documentation.
tax_info_18S <- IdTaxa(dna_18S, trainingset, type="extended", strand = "both", processors = 88, threshold = 50) 

## making and writing out standard output files:

# creating table of taxonomy and setting any that are unclassified as "NA". The IdTaxa step collapses a lot of ranks, so we have to also propogate the last known rank along the rows and specify which taxa have been carried down the row so we know it's not the correct taxon (eg, columns listed with the class as the class, order, and family with unclassified genus and no species)
tax_tab_18S <- data.frame() #initialize the df

tax_tab_18S <- idtax2df(tax_info_18S, 
                        db = "pr2", #using "silva"here doesn't work with our custom formulated db, so use "pr2" and make some minor modifications after producing the table
                        ranks = c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"),
                        boot = 50,
                        rubric = dna_18S,
                        return.conf = FALSE)

#modify table to be ready for analysis
tax_tab_18S[,c(2,11)] <- NULL #remove unwanted columns
colnames(tax_tab_18S) <-c("ASV_ID", "Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "genus", "species")
row.names(tax_tab_18S) <- tax_tab_18S[,1]
tax_tab_18S[,1] <- NULL #remove ASV_ID column for now

# Convert entire dataframe to character matrix
mat <- as.matrix(tax_tab_18S)
mat[] <- as.character(mat)

# Preprocess: convert empty strings and "NA" strings to actual NA
mat[mat == ""] <- NA
mat[mat == "NA"] <- NA

# Get column information
n_cols <- ncol(mat)
cols_to_process <- 2:n_cols  # Columns to process (skip first column)

# Process each row
new_mat <- t(apply(mat, 1, function(row_vec) {
  base <- row_vec[1]  # Initialize base with first column value
  run_length <- 0     # Track consecutive replacements
  result <- row_vec[cols_to_process]  # Initialize results
  
  for (j in 2:n_cols) {  # Process columns 2 to end
    current <- row_vec[j]
    trigger <- FALSE
    
    # Check trigger conditions:
    if (is.na(current)) {
      trigger <- TRUE
    } else if (grepl("uncultured|unclassified", current, ignore.case = TRUE)) {
      trigger <- TRUE
    } else if (!is.na(base) && current == base) {
      trigger <- TRUE
    }
    
    if (trigger) {
      run_length <- run_length + 1
      if (is.na(base)) {
        new_val <- NA_character_
      } else {
        # Create collapsed x's (e.g., "_xx" instead of "_x_x")
        new_val <- paste0(base, "_", paste(rep("x", run_length), collapse = ""))
      }
      result[j-1] <- new_val
    } else {
      run_length <- 0
      base <- current  # Update base to current value
      result[j-1] <- current
    }
  }
  return(result)
}))

# Update dataframe columns (skip first column), convert all empty strings to proper NA values
tax_tab_18S[, cols_to_process] <- as.data.frame(new_mat, stringsAsFactors = FALSE)
tax_tab_18S[tax_tab_18S == ""] <- NA
tax_tab_18S[tax_tab_18S == "NA"] <- NA

# creating vector of desired ranks
ranks <- c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "genus", "species")

#we added leading zeroes earlier to allow correct sorting order by the idtax2df function, now remove them so the ASV IDs match the other IDs in the rest of our pipeline
pattern <- "ASV0*"
replacement <- "ASV"
ASV_names <- gsub(pattern, replacement, row.names(tax_tab_18S)) #save for later
row.names(tax_tab_18S) <- ASV_names
#fixing "Incertae_Sedis" rank labels, which aren't informative and should be replaced based on the rank AFTER, rather than the rank PREVIOUS to avoid propogating high ranks (ie. Eukaryota) all the way up to genus rank in the case of ex. Telonema
tax_tab_18S <- tax_tab_18S %>%
  pmap_dfr(~ {
    row_vec <- c(...)
    n <- length(row_vec)
    current <- row_vec[n]
    if (n > 1) {
      for (j in (n-1):1) {
        # Safely check for pattern with isTRUE() to handle NA and other edge cases
        if (isTRUE(str_detect(row_vec[j], regex("incertae_sedis", ignore_case = TRUE)))) {
          row_vec[j] <- current
        }
        current <- row_vec[j]
      }
    }
    setNames(as.list(row_vec), names(tax_tab_18S))
  })
tax_tab_18S <- as.data.frame(tax_tab_18S) #above operations convert taxa table to tibble. change from tibble back to data frame to avoid weird stuff later on
row.names(tax_tab_18S) <- ASV_names #need to reset row names after this step

#add row names back in as a column
df <- cbind(row.names(tax_tab_18S), tax_tab_18S)
colnames(df) <- c("ASV_ID", "Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "genus", "species")
row.names(df) <- NULL # we don't need row names
tax_tab_18S <- df

#additional modifications for species column
#"Clade environmental"
tax_tab_18S[which(grepl("environmental", tax_tab_18S$species) == TRUE),]
species_sub <- as.character(tax_tab_18S[which(grepl("environmental", tax_tab_18S$species) == TRUE),7]) #removing a few odd spots where the species label is "metagenome"
species_sub <- gsub("_x+", "", species_sub) #remove extra Xs for species label
species_sub <- paste0(species_sub, " sp.", sep="")
tax_tab_18S$species[which(grepl("environmental", tax_tab_18S$species) == TRUE)] <- species_sub
#"metagenome"
species_sub <- as.character(tax_tab_18S[which(tax_tab_18S$species == "metagenome"),7]) #removing a few odd spots where the species label is "metagenome"
species_sub <- gsub("_x+", "", species_sub) #remove extra Xs for species label
species_sub <- paste0(species_sub, " sp.", sep="")
tax_tab_18S$species[which(tax_tab_18S$species == "metagenome")] <- species_sub
#View(tax_tab_18S) #checking over manually. looks good
#other various edge cases that should be fixed
tax_tab_18S$species <- gsub("unidentified_", "", tax_tab_18S$species) #remove "unidentified_" prefix (rare, and will be remedied by appending sp. in the following commands)
tax_tab_18S$species <- gsub("cf._", "", tax_tab_18S$species) #remove "cf." labels
tax_tab_18S$species <- gsub("_sp.", " sp.", tax_tab_18S$species) #remove unnecessary underscores
tax_tab_18S$species <- gsub("_x+", " sp.", tax_tab_18S$species) #remove extra Xs and add "sp." for species labels
tax_tab_18S$species <- gsub("sp. sp.", "sp.", tax_tab_18S$species) #remove doubled "sp." suffixes
tax_tab_18S$species <- gsub("_", " ", tax_tab_18S$species) #replace underscores with spaces
tax_tab_18S$species_wc <- stringr::str_count(tax_tab_18S$species, "\\S+")
tax_tab_18S$species[which(tax_tab_18S$species_wc == 1)] <- paste0(tax_tab_18S$species[which(tax_tab_18S$species_wc == 1)], " sp.")
tax_tab_18S$species_wc <- NULL #remove word count column before saving

#save taxa table
write.table(tax_tab_18S, "taxonomy_table.18S.DECIPHER.SILVAv138.txt", sep = "\t", quote = F, row.names = FALSE)

############
####TODO####
############
#  #  #  #  strip extra stuff from species and genus columns. if first word from each column doesn't match, apply first word from species column to genus column. if they do match, make no changes

# also create corresponding confidence table. this is useful for interpretation, but assignments can be filtered above by passing a threshold to the IdTaxa() function
# creating vector of desired ranks
header <- c("ASV_ID", "Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "genus", "species")

conf_tab_18S <- idtax2df(tax_info_18S, 
                        db = "pr2", #using "silva"here doesn't work with our custom formulated db, so use "pr2" and make some minor modifications after producing the table
                        ranks = c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "genus", "species"),
                        boot = 50,
                        rubric = dna_18S,
                        return.conf = TRUE)

conf_tab_18S <- as.data.frame(conf_tab_18S[2])
conf_tab_18S[,c(2,11)] <- NULL
colnames(conf_tab_18S) <- header

#necessary?
pattern <- "ASV0*"
replacement <- "ASV"
confidences$ASV_ID <- gsub(pattern, replacement, confidences$ASV_ID)

write.table(conf_tab_18S, "taxonomy_table.18S.DECIPHER.SILVAv138.conf.txt", sep = "\t", quote = F, row.names = FALSE)

#### metazoogene annotations ####
#### slightly different procedure again for MetaZooGene reference set ####
## loading SILVA v138 DECIPHER train set #see file silva_138_reformat_rescript.sh to see how these data were reformatted for this step from the raw silva files
trainingset <- readRDS("/mnt/Genomics/Working/databases/amplicon_taxonomyDBs/DECIPHER_IDTAXA/decipher_metaZooGene/18S/MZG-18S-concise.trained.rds")

#taxonomic inference #threshold for accepting taxonomic assignments can be adjusted with the 'threshold' parameter. see documentation.
tax_info_18S <- IdTaxa(dna_18S, trainingset, type="extended", strand = "both", processors = 112, threshold = 50) 

## making and writing out standard output files:

# creating table of taxonomy and setting any that are unclassified as "NA". The IdTaxa step collapses a lot of ranks, so we have to also propogate the last known rank along the rows and specify which taxa have been carried down the row so we know it's not the correct taxon (eg, columns listed with the class as the class, order, and family with unclassified genus and no species)
tax_tab_18S <- data.frame() #initialize the df

tax_tab_18S <- idtax2df(tax_info_18S, 
                        db = "pr2", #using "silva"here doesn't work with our custom formulated db, so use "pr2" and make some minor modifications after producing the table
                        ranks = c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"),
                        boot = 50,
                        rubric = dna_18S,
                        return.conf = FALSE)

#modify table to be ready for analysis
tax_tab_18S[,c(2,10,11)] <- NULL #remove unwanted columns
colnames(tax_tab_18S) <-c("ASV_ID", "Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "genus", "species")
row.names(tax_tab_18S) <- tax_tab_18S[,1]
tax_tab_18S[,1] <- NULL #remove ASV_ID column for now

# Convert entire dataframe to character matrix
mat <- as.matrix(tax_tab_18S)
mat[] <- as.character(mat)

# Preprocess: convert empty strings and "NA" strings to actual NA
mat[mat == ""] <- NA
mat[mat == "NA"] <- NA

# Get column information
n_cols <- ncol(mat)
cols_to_process <- 2:n_cols  # Columns to process (skip first column)

# Process each row
new_mat <- t(apply(mat, 1, function(row_vec) {
  base <- row_vec[1]  # Initialize base with first column value
  run_length <- 0     # Track consecutive replacements
  result <- row_vec[cols_to_process]  # Initialize results
  
  for (j in 2:n_cols) {  # Process columns 2 to end
    current <- row_vec[j]
    trigger <- FALSE
    
    # Check trigger conditions:
    if (is.na(current)) {
      trigger <- TRUE
    } else if (grepl("uncultured|unclassified", current, ignore.case = TRUE)) {
      trigger <- TRUE
    } else if (!is.na(base) && current == base) {
      trigger <- TRUE
    }
    
    if (trigger) {
      run_length <- run_length + 1
      if (is.na(base)) {
        new_val <- NA_character_
      } else {
        # Create collapsed x's (e.g., "_xx" instead of "_x_x")
        new_val <- paste0(base, "_", paste(rep("x", run_length), collapse = ""))
      }
      result[j-1] <- new_val
    } else {
      run_length <- 0
      base <- current  # Update base to current value
      result[j-1] <- current
    }
  }
  return(result)
}))

# Update dataframe columns (skip first column), convert all empty strings to proper NA values
tax_tab_18S[, cols_to_process] <- as.data.frame(new_mat, stringsAsFactors = FALSE)
tax_tab_18S[tax_tab_18S == ""] <- NA
tax_tab_18S[tax_tab_18S == "NA"] <- NA


# creating vector of desired ranks
ranks <- c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "genus", "species")

#we added leading zeroes earlier to allow correct sorting order by the idtax2df function, now remove them so the ASV IDs match the other IDs in the rest of our pipeline
pattern <- "ASV0*"
replacement <- "ASV"
ASV_names <- gsub(pattern, replacement, row.names(tax_tab_18S)) #save for later
row.names(tax_tab_18S) <- ASV_names
tax_tab_18S <- as.data.frame(tax_tab_18S) #above operations convert taxa table to tibble. change from tibble back to data frame to avoid weird stuff later on
row.names(tax_tab_18S) <- ASV_names #need to reset row names after this step

#add row names back in as a column
df <- cbind(row.names(tax_tab_18S), tax_tab_18S)
colnames(df) <- c("ASV_ID", "Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "genus", "species")
row.names(df) <- NULL # we don't need row names
tax_tab_18S <- df

#additional modifications for species column
#other various edge cases that should be fixed
tax_tab_18S$species <- gsub("unidentified_", "", tax_tab_18S$species) #remove "unidentified_" prefix (rare, and will be remedied by appending sp. in the following commands)
tax_tab_18S$species <- gsub("cf._", "", tax_tab_18S$species) #remove "cf." labels
tax_tab_18S$species <- gsub("_sp.", " sp.", tax_tab_18S$species) #remove unnecessary underscores
tax_tab_18S$species <- gsub("_x+", " sp.", tax_tab_18S$species) #remove extra Xs and add "sp." for species labels
tax_tab_18S$species <- gsub("sp. sp.", "sp.", tax_tab_18S$species) #remove doubled "sp." suffixes
tax_tab_18S$species <- gsub("_", " ", tax_tab_18S$species) #replace underscores with spaces
tax_tab_18S$species_wc <- stringr::str_count(tax_tab_18S$species, "\\S+")
tax_tab_18S$species[which(tax_tab_18S$species_wc == 1)] <- paste0(tax_tab_18S$species[which(tax_tab_18S$species_wc == 1)], " sp.")
tax_tab_18S$species_wc <- NULL #remove word count column before saving

#save taxa table
write.table(tax_tab_18S, "taxonomy_table.18S.DECIPHER.MetaZooGene-C.txt", sep = "\t", quote = F, row.names = FALSE)

# also create corresponding confidence table. this is useful for interpretation, but assignments can be filtered above by passing a threshold to the IdTaxa() function
# creating vector of desired ranks
header <- c("ASV_ID", "Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "genus", "species")

conf_tab_18S <- idtax2df(tax_info_18S, 
                         db = "pr2", #using "silva"here doesn't work with our custom formulated db, so use "pr2" and make some minor modifications after producing the table
                         ranks = c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "genus", "species"),
                         boot = 50,
                         rubric = dna_18S,
                         return.conf = TRUE)

conf_tab_18S <- as.data.frame(conf_tab_18S[2])
conf_tab_18S[,c(2,10,11)] <- NULL
colnames(conf_tab_18S) <- header

#necessary?
pattern <- "ASV0*"
replacement <- "ASV"
conf_tab_18S$ASV_ID <- gsub(pattern, replacement, conf_tab_18S$ASV_ID)

write.table(conf_tab_18S, "taxonomy_table.18S.DECIPHER.MetaZooGene-C.conf.txt", sep = "\t", quote = F, row.names = FALSE)

#### formatting final output files ####
#### replace the long ASV names (the actual sequences) with human-readable names ####
#save the new names and sequences as a .fasta file in your project working directory, and save a table that shows the mapping of sequences to new ASV names
my_otu_table <- t(as.data.frame(seqtab)) #transposed (OTUs are rows) data frame
ASV.seq <- as.character(unclass(row.names(my_otu_table))) #store sequences in character vector
ASV.num <- paste0("ASV", seq(ASV.seq), sep='') #create new names
write.table(cbind(ASV.num, ASV.seq), "sequence_ASVname_mapping.18S.txt", sep="\t", quote=F, row.names=F, col.names=F)
write.fasta(sequences=as.list(ASV.seq), names=ASV.num, "18S_ASV_sequences.fasta") #save sequences with new names in fasta format

#IMPORTANT: sanity checks
colnames(seqtab) == ASV.seq #only proceed if this tests as true for all elements

#assign new ASV names
colnames(seqtab) <- ASV.num

#re-save sequence and taxonomy tables with updated names
write.table(data.frame("row_names"=rownames(seqtab),seqtab),"sequence_table.seqtab18S.combined.w_ASV_names.txt", row.names=FALSE, quote=F, sep="\t")
