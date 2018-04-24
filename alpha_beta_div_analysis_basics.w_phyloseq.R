#basic procedure for loading in a .biom file, processing data with phyloseq, doing diversity analyses
#NOTE: THIS SCRIPT IS MEANT TO BE CHANGED TO FIT YOUR DATA. PLEASE REMEMBER TO:
#       1. MAKE A COPY OF THIS SCRIPT IF YOU WISH TO MODIFY IT
#       2. CHANGE ALL GENERALIZED PARAMETERS/VARIABLES (EX: "FACTOR_1") TO MATCH YOUR DATA

#### load libraries ####
#you must install these first if you want to load the data in using phyloseq and process with deseq
library(tidyverse)
# Load the reshape2 package for converting between long and wide format data
library(reshape2)
# Load the stringr package for improved filtering of text
library(stringr)
# Load the ape package for reading and modifying phylogenetic trees
library(ape)
# Load the phyloseq package for microbial community analysis
library(phyloseq)
# Load the data.table package for better metadata manipulation
library(data.table)
# Load the viridis package for colour palettes for continuous data
library(viridis)
# Load the qualpalr package for colour palettes for qualitative data
library(qualpalr)
# load the ggplot2 package for visualization of data
library(ggplot2)
#load the vegan library
library(vegan)

#### environment settings ####
#set working directory
setwd("/path/to/working/directory/")

#### load in .biom file, metadata, and (optional) phylogenetic tree ####
# REQUIRED: Read the raw data in .biom format, and store as a phylo_seq object called "rawdata"
rawdata <- import_biom(file.path("/file/path/", "OTU_Table.biom"), # file.path() is used for cross-platform compatibility
                       parallel = TRUE,
                       trim_ws = TRUE) # use multiple processor cores for shorter read times

# REQUIRED: Read the raw metadata in .txt format, and store as a tibble (an improved dataframe) called "rawmetadata"
rawmetadata <- read_delim(file = file.path("/file/path/", "metadata.txt"), # file.path() is used for cross-platform compatibility
                          "\t", # the metadata file must be tab delimited and in .txt format
                          escape_double = FALSE, # the imported text file does not 'escape' quotation marks by wrapping them with more quotation marks
                          trim_ws = TRUE) # remove leading and trailing spaces from character string entries

# IMPORTANT: ROW NAMES OF THE RAW METADATA FILE MUST MATCH SAMPLE NAMES FROM THE OTU TABLE. IF NOT, THE PHYLOSEQ OBJECT WILL NOT BE MERGED PROPERLY AND NON-MATCHING IDS WILL BE THROWN OUT
# RECOMMENDED: preserve a list of samples that are unique to the OTU table and metadata file to make sure there aren't unexpected differences
# NOTE: change MY_SAMPLE_IDs to match the column where sample IDs are stored in the metadata file. in QIIME compatible metadata files, these will be under `#SampleID`
notinmeta <- setdiff(sample_names(rawdata), rawmetadata$MY_SAMPLE_IDs)
# if there are any samples in "notinmeta", remove these samples from the raw data
allsamples <- sample_names(rawdata)
allsamples <- allsamples[!(allsamples %in% notinmeta)]
biomdata <- prune_samples(allsamples, rawdata)
# Find all samples in raw metadata that are not in raw data
notinraw <- setdiff(rawmetadata$MY_SAMPLE_IDs, sample_names(rawdata))
# Remove these samples from the metadata
metadata <- filter(rawmetadata, !(MY_SAMPLE_IDs %in% notinraw))
# Create vector consisting of ALL removed samples; useful for record-keeping
symmdiff <- c(notinmeta, notinraw)
# Format metadata tibble into the phyloseq metadata format
metadata <- sample_data(metadata)
# Add rownames to the phyloseq metadata which correspond to the OTU IDs
rownames(metadata) <- metadata$MY_SAMPLE_IDs

# OPTIONAL: Read in a phylogenetic tree in .tre format, and store as a tree object called "rawtreedata"
# IMPORTANT: TREE TIPS MUST MATCH OTU IDs FROM TAXA TABLE. THIS MEANS IF YOU ARE USING AN EPA PLACEMENT TREE, YOU MUST REMOVE "QUERY___" from each tree tip label BEFORE IMPORTING if you wish to use the tree in a phyloseq object
# file.path() is used for cross-platform compatibility
rawtreedata <- read_tree(file.path("/file/path", "phylo_tree.tre"))

#IMPORTANT NOTE: at this point you should make sure your sample IDs in the data, metadata, and tree data objects match

#### OPTIONAL: drop unwanted levels from metadata now, before converting to phyloseq ####
metadata <- metadata[which(metadata$SOME_VARIABLE != "UNWANTED VALUE"), ] #works with factor-formatted or continuous variables
#### create phyloseq object with completed metadata, otu table, and tree ####
project_data <- merge_phyloseq(biomdata, metadata, rawtreedata)
#filtering steps, if not already done before loading into R
#filter out samples with less than 1000 reads (arbitrary threshold and generally the minimum, choose your own, use sample_counts() to look at distribution of read counts per sample)
project_data <- prune_samples(sample_sums(project_data) >= 1000, project_data) 
# Remove OTUs with less than N total reads. (N = 250 in example) 
project_data <- prune_taxa(taxa_sums(project_data) >= 250, project_data)
# 16S ONLY: Remove mitochondrial and chloroplast OTUs #IMPORTANT: make sure that the filter terms will work with your taxonomy strings, and ranks
project_data <- project_data %>%
  subset_taxa(Rank5 != "__Mitochondria") %>% 
  subset_taxa(Rank3 != "__Chloroplast")

# 18S (and optional for 16S): Remove unwanted clades 
project_data <- project_data %>%
  subset_taxa(Rank5 != "UNWANTED_HOST_FAMILY") %>% 
  subset_taxa(Rank7 != "UNWANTED_CLADE")

# Preserve unassigned taxa
project_data.unassigned <- project_data %>%
  subset_taxa(Rank1 == "Unassigned") #works as long as your OTUs without a taxonomy assignment are labeled as "Unassigned". adjust accordingly.
# Remove unassigned taxa
project_data <- project_data %>%
  subset_taxa(Rank1 != "Unassigned")
# Remove counts that may represent noise, use a threshold (we are using a threshod of 2 reads for most datasets. be sure to choose the correct value for your own data.)
otu <- as.data.frame(otu_table(project_data)) #get OTU table
otu_table(project_data)[otu <= 2] <- 0 #for entries where the raw abundance of an OTU in a sample is less than N (N=2 in the example), set the raw read count to 0

# OPTIONAL: modify Rank labels in taxa table (check the colnames of the tax_table(project_data) object to see if you want to change them)
colnames(tax_table(project_data)) <- c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7")

#### plot rarefaction curves ####
plot(sort(sample_sums(project_data))) #looking at sample read counts
summary(sample_sums(project_data))

#found rarefaction curve function here: https://github.com/mahendra-mariadassou/phyloseq-extended (richness.R) 
#to use this you have to load the "ggrare" function in from richness.R #it will throw errors if you are labeling with groups that have too few categories
source("/path/to/phyloseq-extended/R/richness.R") # load in ggrare function
p <- ggrare(project_data, step = 1000, color = "factor", se = FALSE)
# OPTIONAL: facet_wrap the plot
p + facet_wrap(~FACTOR_1 + FACTOR_2)
#you can use the plot above to judge a rough cutoff for rarefaction. you can also do this with QIIME's alpha rarefaction script

#you can use the "which" command like the example below to see which samples you'll lose for a given cutoff value
which(sample_sums(project_data) < 20000)

#### rarefy data ####
set.seed(24) #you must set a numerical seed like this for reproducibility
project_data.rarefied <- rarefy_even_depth(project_data, sample.size = min(sample_sums(project_data)))

#### Create Colour Palettes ####
# 1. identify the number of colors you need from the factor you want to plot by
numcol <- length(unique(sample_data(project_data)$FACTOR_1)) #EXAMPLE ONLY, adjust per object being plotted
# 2. use a number seed to determine how qualpar samples your colors from its palette
set.seed(13)
# 3. use qualpalr colour palettes for easily distinguishing taxa
newpal <- qualpal(n = numcol, colorspace = "pretty")
# 4. If you want to color only a few things in the plot, try using a more colourblind-friendly palette:
cbPalette <- c("#E69F00", "#56B4E9", "#000000", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#FFFF00", "#999999", "#FF00FF", "#F0E442", "#FFFFFF", "#00FFFF")


#### basic alpha div plot ####
#chao1
pdf("AlphaDiversity.chao1.factor.experiment_name.pdf", #name of file to print. can also include relative or absolute path before filename.
    width = 16, height = 9)# define plot width and height. completely up to user.
p <- plot_richness(project_data.rarefied, x="factor", color = "factor", measures=c("Chao1"))
p + geom_boxplot(outlier.colour = "red", outlier.shape = 13) + 
  facet_grid(~ FACTOR_2) + #use this to divide data into separate plots based on a factor/variable
  theme(strip.background = element_rect(fill="white"), strip.placement = "bottom") +
  theme(strip.text = element_text(colour = 'black')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values=cbPalette) + scale_colour_manual(values=cbPalette) + 
  labs(title="Alpha Diversity, Factor N, Chao1", x="Factor1 ~ Factor2", y="Richness/Alpha Diversity")
dev.off()
#### calculate alpha diversity and add it as a column in the metadata ####
project_data.chao1 = estimate_richness(project_data.rarefied, split = TRUE, measures = c("Chao1")) #estimate richness
sample_data(project_data.rarefied)$chao1 <- project_data.chao1$Chao1 #add to metadata (the rows are in the same order already)
sample_data(project_data.rarefied)$chao1 <- as.numeric(sample_data(project_data.rarefied)$chao1)

#### beta diversity (NMDS, PCoA, etc.) ####
#do ordinations
set.seed(24)
NMDS.bray <- ordinate(
  physeq = project_data, 
  method = "NMDS", 
  distance = "bray"
) # you can choose different methods and distance metrics, see the ordinate help page for details. this function works with "phyloseq" class objects.

#### making beta div plots ####
#we get more plotting control if we don't use the phyloseq plotting functions for ordination plots, and instead add the results of the ordination to our existing metadata
NMDS <- as.data.frame(sample_data(project_data))
bray <- as.data.frame(NMDS.bray$points)
row.names(bray) == row.names(NMDS) #sanity check #tests as true
NMDS$NMDS.bray1 <- bray$MDS1
NMDS$NMDS.bray2 <- bray$MDS2
#OPTIONAL: sort data for better looking easier to read plots
NMDS.sort <- NMDS[order(NMDS$FACTOR_1, NMDS$FACTOR_2),]

#plain NMDS plot colored by "FACTOR_1" and shaped by "FACTOR_2"
pdf("NMDS.expt_name.FACTOR_1_color.FACTOR_2_shape.pdf"
    , width = 16 # Default is 7
    , height = 9 # Change to 10; make it taller
)
p <- ggplot(NMDS, aes(x=NMDS.bray1, y=NMDS.bray2, shape = FACTOR_2, color = FACTOR_1)) # change the first argument to NMDS.sort if the optional command was ran
p + geom_point(size=4) + scale_shape_manual(values=1:nlevels(NMDS.sort$FACTOR_2)) +
  labs(title="NMDS Factor1 & Factor2") + 
  scale_fill_manual(values=cbPalette) + scale_colour_manual(values=cbPalette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

#facet wrapped NMDS plot
pdf("NMDS.expt_name.FACTOR_1~FACTOR_2.pdf"
    , width = 16 # Default is 7
    , height = 8 # Change to 10; make it taller
)
p <- ggplot(NMDS.sort, aes(x=NMDS.bray1, y=NMDS.bray2, color = FACTOR_3, shape = FACTOR_2))
p + facet_wrap(~FACTOR_1) + geom_point(size=4) +
  theme(strip.background = element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0.2, "lines")) +
  labs(title="NMDS Factor2 & Factor3 ~ Factor1") + 
  scale_fill_manual(values=cbPalette) + scale_colour_manual(values=cbPalette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()