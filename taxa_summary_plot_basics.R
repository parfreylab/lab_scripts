#####################################################################################
##General Framework for Making Taxa Summary Plots from a QIIME-formatted .biom File##
#####################################################################################
#recommend using RStudio for this

#### Load R Packages #### #you must install these first if you want to make the plots with ggplot, and load the data in using phyloseq
library(tidyverse)
# Load the reshape2 package for converting between long and wide format data
library(reshape2)
# Load the stringr package for improved filtering of text
library(stringr)
# Load the ape package for reading and modifying phylogenetic trees
library(ape)
# Load the phyloseq package for microbial community analysis
library(phyloseq)
# Load the viridis package for colour palettes for continuous data
library(viridis)
# Load the qualpalr package for colour palettes for qualitative data
library(qualpalr)

#### Environment Settings ####
# Set global ggplot2 theme; this is most common ggplot format for publication
theme_set(theme_bw())
#set working directory
setwd("/path/to/working/directory")

#### Load Required/Optional Files ####
# REQUIRED: Read the raw data in .biom format, and store as a phylo_seq object called "rawdata"
rawdata <- import_biom(file.path("relative/or/absolute/path/to/", "OTU_table.w_taxa.biom"), # file.path() is used for cross-platform compatibility
                       parallel = TRUE) # use multiple processor cores for shorter read times

# REQUIRED: Read the raw metadata in .txt format, and store as a tibble (an improved dataframe) called "rawmetadata"
rawmetadata <- read_delim(file = file.path("relative/or/absolute/path/to/", "metadata.complete.txt"), # file.path() is used for cross-platform compatibility
                          "\t", # the text file is tab delimited
                          escape_double = FALSE, # the imported text file does not 'escape' quotation marks by wrapping them with more quotation marks
                          trim_ws = TRUE, # remove leading and trailing spaces from character string entries
                          )

# OPTIONAL: Read in the raw phylogenetic tree in .tre format, and store as a tree object called "rawtreedata"
# file.path() is used for cross-platform compatibility
rawtreedata <- read_tree(file.path("tree", "16s_makephylo_fasttree.tre"))

#### Initial Data Processing ####
# 1. Find all samples in raw data that are not in raw metadata
notinmeta <- setdiff(sample_names(rawdata), rawmetadata$`#SampleID`)

# 2. if there are any samples in "notinmeta", remove these samples from the raw data
allsamps <- sample_names(rawdata)
allsamps <- allsamps[!(allsamps %in% notinmeta)]
data <- prune_samples(allsamps, rawdata)

# 3a. Find all samples in raw metadata that are not in raw data
notinraw <- setdiff(rawmetadata$`#SampleID`, sample_names(rawdata))

# 3b. Remove these samples from the metadata
metadata <- filter(rawmetadata, !(`#SampleID` %in% notinraw))

# 4. Create vector consisting of ALL removed samples; useful for record-keeping
symmdiff <- c(notinmeta, notinraw)

# 5. Format metadata tibble into the phyloseq metadata format
metadata <- sample_data(metadata)

# 6. Add rownames to the phyloseq metadata which correspond to the OTU IDs
rownames(metadata) <- metadata$`#SampleID`

# 7. merge raw data, metadata, and tree data into single phyloseq object
project_data <- merge_phyloseq(data, metadata, rawtreedata)

# 8. modify Rank labels in taxa table
colnames(tax_table(project_data)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#### Quality and Taxa Filtering ####
# 1. look at minimum, mean, and maximum sample counts, if desired
smin <- min(sample_sums(project_data))
smean <- mean(sample_sums(project_data))
smax <- max(sample_sums(project_data))

# 2. Remove samples with less than N reads. (N = 1000 in example) 
project_data <- prune_samples(sample_sums(project_data) >= 1000, project_data)

# 3. Remove OTUs with less than N total reads. (N = 250 in example) 
project_data <- prune_taxa(taxa_sums(project_data) >= 250, project_data) 

# 4. Remove mitochondrial and chloroplast OTUs
project_data <- project_data %>%
  subset_taxa(Family != "Mitochondria") %>%
  subset_taxa(Class != "Chloroplast")

#### Create Plotting Objects ####
# 1. reshape data based on taxonomic level you are interested in
family_plot <- project_data %>%
  tax_glom(taxrank = "Family") %>%                      # agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>%  # Transform to rel. abundance
  psmelt() %>%                                          # Melt to long format
  arrange(Family)

# 2. make aesthetic changes to taxa rank labels, if desired:
family_plot$Family <- str_replace_all(family_plot$Family, "__", "") # Remove underscores
family_plot$Order <- str_replace_all(family_plot$Order, "__", "") # Remove underscores
family_plot$Class <- str_replace_all(family_plot$Class, "__", "") # Remove underscores
family_plot$Phylum <- str_replace_all(family_plot$Phylum, "__", "") # Remove underscores

# 3. order levels you are interested in the way you would like them plotted. many ways to do this.
family_plot$factor_N <- factor(family_plot$factor_N, levels = c("A", "B", "C", "1", "2", "3", "z", "y", "x", "w")) #method 1
family_plot$factor_N <- factor(family_plot$factor_N, levels = sort(unique(factor(family_plot$factor_N)))) #method 2
mylevels <- c("A", "B", "C", "1", "2", "3", "z", "y", "x", "w") #method 3
family_plot$factor_N <- factor(family_plot$factor_N, levels = mylevels)

#### Create Color Palette(s) ####
# 1. identify the number of colors you need
numcol <- length(unique(rat_family$Family))
# 2. use a number seed to determine how qualpar samples your colors from its palette
set.seed(15)
# 3. use qualpalr colour palettes for easily distinguishing taxa
newpal <- qualpal(n = numcol, colorspace = "pretty")

# 4. If you want to color labels in the plot, try using a more colourblind-friendly palette:
myPalette <- c("#00FFFF", "#56B4E9", "#000000", "#009E73", "#F0E442", "#E69F00", "#D55E00", "#CC79A7", "#999999", "#FF00FF")

#### Make Plots ####
#example plot showing relative abundance of taxa through timepoints (coded as "experiment_day"), with individuals grouped together (coded as "rat_name"), and animals in the same cage grouped together as well (coded as "cage")
ggplot(family_plot, aes(x = experiment_day, y = Abundance, fill = Family)) + 
  facet_wrap(cage~rat_name, ncol=4, strip.position = "top") + #facet_wrap is the function for determining how plots are grouped within a multi-plot space
  theme_bw() +
  theme(strip.background = element_rect(fill="white")) + #these "theme" settings determine how the facet grid looks on the plot
  theme(strip.text = element_text(colour = 'black')) +
  geom_bar(stat = "identity", width = 1.0) + #geom_bar controls the bar plots themselves
  scale_y_continuous(expand = c(0.005,0.005)) + #this controls the y axis scale, for bigger margins above and below, increase the numbers provided
  scale_fill_manual(values = newpal$hex) + #here's where to use the color palate derived above
  ylab("Relative Abundance") + #x and y axis labels
  xlab("Experiment Day") +
  ggtitle("Taxonomic Composition of Rat Gut Bacterial Communities Grouped by Cage ~ Individual") + #plot title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0, "lines")) #another "theme" command. a lot of extra control can be established here. this line ensures that there is no padding in the multi-plot grid

#similar example plot where the number of columns in the multi-plot is specified directly
ggplot(family_plot, aes(x = experiment_day, y = Abundance, fill = Family)) + 
  facet_wrap(~ rat_name, ncol=10) + #here we only group the data by one factor, instead of two like the example above
  theme_bw() +
  theme(strip.background = element_rect(fill="white")) +
  theme(strip.text = element_text(colour = 'black')) +
  geom_bar(stat = "identity", width = 1.0) +
  scale_y_continuous(expand = c(0.005,0.005)) +
  scale_fill_manual(values = newpal$hex) +
  ylab("Relative Abundance") +
  xlab("Experiment Day") +
  ggtitle("Taxonomic Composition of Rat Gut Bacterial Communities Through Time") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0, "lines"))

#in this example the data is grouped in such a way that the relative abundances for each sample adds to the total height of the bar they are grouped into
ggplot(rat_family, aes(x = experiment_day, y = Abundance, fill = Family)) + 
  facet_wrap(~ infected) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white")) +
  theme(strip.text = element_text(colour = 'black')) +
  geom_bar(stat="identity", width = 1.0) +
  scale_y_continuous(expand = c(0.005,0.005)) +
  scale_fill_manual(values = newpal$hex) +
  ylab("Relative Abundance") +
  xlab("Experiment Day") +
  ggtitle("Taxonomic Composition of Rat Gut Bacterial Communities by Hymenolepis Infection Status") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0, "lines"))

#in this example the x axis labels are coloured by an additional colour palette #be careful with this setting, it may not apply the colors correctly. I have had better luck with a string of colors which I define manually, in the order I want them used in the plot.
ggplot(rat_family, aes(x = rat_name, y = Abundance, fill = Family)) + 
  facet_grid(experiment_day~.) +
  theme_bw() +
  theme(strip.background = element_rect(fill="white"), axis.text.x = element_text(colour=myPalette)) +
  theme(strip.text = element_text(colour = 'black')) +
  geom_bar(stat = "identity", width = 1.0) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  scale_fill_manual(values = newpal$hex) +
  ylab("Relative Abundance") +
  xlab("Experiment Day") +
  ggtitle("Taxonomic Composition of Rat Gut Bacterial Communities Grouped by Experiment Day ~ Individual") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0, "lines"))


#### Print Plots to PDF ####
pdf("Taxa_summary.cage~individual.pdf", #name of file to print. can also include relative or absolute path before filename.
    , width = 16 # define plot width and height. completely up to user.
    , height = 9
)
# plot commands for an individual plot or single multi-plot go here
dev.off()