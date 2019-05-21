##########################################################################
##General Framework for Making Taxa Summary Plots From a Phyloseq Object##
##########################################################################
#author: Evan Morien
#last modified: May 17th, 2019

#recommend using RStudio for this
#NOTE: THIS SCRIPT IS MEANT TO BE CHANGED TO FIT YOUR DATA. PLEASE REMEMBER TO:
#       1. MAKE A COPY OF THIS SCRIPT IF YOU WISH TO MODIFY IT
#       2. CHANGE ALL GENERALIZED PARAMETERS/VARIABLES (EX: "FACTOR_1") TO MATCH YOUR DATA

#IMPORTANT NOTE: taxa summary plots should be made with non-rarefied data.

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
#please see the script "load_data_into_phyloseq_object.R" for instructions on loading and filtering your data using phyloseq. Use the phyloseq object you've filtered as input for the next section of this script.

#### Create Plotting Objects ####
# 1. reshape data based on taxonomic level you are interested in
taxonomy_plot_obj <- project_data %>%
  tax_glom(taxrank = "Rank5") %>%                      # agglomerate at family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>%  # Transform to rel. abundance
  psmelt() %>%                                          # Melt to long format
  arrange(Rank5)

# 2. make aesthetic changes to taxa rank labels, if desired:
taxonomy_plot_obj$Rank5 <- str_replace_all(taxonomy_plot_obj$Rank5, "__", "") # Remove underscores if your taxonomy strings have them. you'd do this for every rank in turn.
taxonomy_plot_obj$Rank4 <- str_replace_all(taxonomy_plot_obj$RankN, "STRING", "REPLACEMENT") # Generalized example of string replacement

# 3. order levels you are interested in the way you would like them plotted. many ways to do this.
taxonomy_plot_obj$factor_N <- factor(taxonomy_plot_obj$factor_N, levels = c("A", "B", "C", "1", "2", "3", "z", "y", "x", "w")) #method 1
taxonomy_plot_obj$factor_N <- factor(taxonomy_plot_obj$factor_N, levels = sort(unique(factor(taxonomy_plot_obj$factor_N)))) #method 2
mylevels <- c("A", "B", "C", "1", "2", "3", "z", "y", "x", "w") #method 3
taxonomy_plot_obj$factor_N <- factor(taxonomy_plot_obj$factor_N, levels = mylevels)

#### Create Color Palette(s) ####
# 1. identify the number of colors you need
numcol <- length(unique(taxonomy_plot_obj$Rank5))
# 2. use a number seed to determine how qualpar samples your colors from its palette
set.seed(15)
# 3. use qualpalr colour palettes for easily distinguishing taxa
newpal <- qualpal(n = numcol, colorspace = "pretty")

# 4. If you want make your plot better, try using a more colourblind-friendly palette:
cbPalette <- c("#E69F00", "#56B4E9", "#000000", "#009E73", "#CC79A7", "#D55E00", "#0072B2", "#FFFF00", "#999999", "#FF00FF", "#F0E442", "#FFFFFF", "#00FFFF")

# 5. If you want to make clade-specific colour assignments, and make multiple plots that have the same colours for the same clades, the best way to do this is to create your own palette that assigns colours to clades
#step 1: Identify unique clade names in each of your plots you plan to make (do this with the unique() function, input should be the column of the plotting object containing the clade labels)
#step 2: Use the newpal creation code in the lines just above to make a palette with N colours, where N is the number of unique clades you need colours for.
#step 3: The tedious part. You need to write out the palette manually. When you're done, it will look like the line of code below
myCustomPalette <- c("clade1" = "#333BFF",  "clade2"="#BBB3FF",  "clade3" = "#CC6600" ,  "etc."="#9633FF") 
#this palette can be as large as you like, but the more colours you add, the harder it is to discriminate them visually. I recommend 20 colours MAXIMUM.

#### Make Plots ####
#example plot showing relative abundance of taxa through timepoints (coded as "experiment_day"), with individuals grouped together (coded as "rat_name"), and animals in the same cage grouped together as well (coded as "cage")
ggplot(taxonomy_plot_obj, aes(x = Sample, y = Abundance, fill = Rank5)) + 
  facet_wrap(FACTOR_1~FACTOR_2, ncol=4, strip.position = "top", drop=TRUE, scales="free") + #OPTIONAL LINE: facet_wrap is the function for determining how plots are grouped within a multi-plot space
  theme_bw() +
  theme(strip.background = element_rect(fill="white")) + #these "theme" settings determine how the facet grid looks on the plot
  theme(strip.text = element_text(colour = 'black')) +
  geom_bar(stat = "identity", width = 1.0) + #geom_bar controls the bar plots themselves
  scale_y_continuous(expand = c(0.005,0.005)) + #this controls the y axis scale, for bigger margins above and below, increase the numbers provided
  scale_fill_manual(values = newpal$hex) + #here's where to use the raw colour palate derived above
  scale_fill_manual(values = cbPalette) + #example with colourblind palette
  #NOTE: to preserve colours between plots, you must use the line below, not the ones above
  scale_fill_manual(values = myCustomPalette) + #example with custom palette
  ylab("Relative Abundance") + #x and y axis labels
  xlab("My X Axis Label") +
  ggtitle("My Plot Title") + #plot title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0, "lines")) #another "theme" command. a lot of extra control can be established here. this line ensures that there is no padding in the multi-plot grid

#similar example plot where the number of columns in the multi-plot is specified directly in the facet_wrap command
ggplot(taxonomy_plot_obj, aes(x = Sample, y = Abundance, fill = Rank5)) + 
  facet_wrap(~ FACTOR_1, ncol=5, drop=TRUE, scales="free") + #here we only group the data by one factor, instead of two like the example above
  theme_bw() +
  theme(strip.background = element_rect(fill="white")) +
  theme(strip.text = element_text(colour = 'black')) +
  geom_bar(stat = "identity", width = 1.0) +
  scale_y_continuous(expand = c(0.005,0.005)) +
  scale_fill_manual(values = newpal$hex) +
  ylab("Relative Abundance") +
  xlab("My X Axis Label") +
  ggtitle("My Plot Title") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0, "lines"))

#in this example the data is grouped in such a way that the relative abundances for each sample adds to the total height of the bar they are grouped into
ggplot(taxonomy_plot_obj, aes(x = FACTOR_N, y = Abundance, fill = Rank5)) + 
  facet_wrap(~ FACTOR_X, drop=TRUE, scales="free") +
  theme_bw() +
  theme(strip.background = element_rect(fill="white")) +
  theme(strip.text = element_text(colour = 'black')) +
  geom_bar(stat="identity", width = 1.0) +
  scale_y_continuous(expand = c(0.005,0.005)) +
  scale_fill_manual(values = newpal$hex) +
  ylab("Relative Abundance") +
  xlab("FACTOR N") +
  ggtitle("MY PLOT TITLE") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.spacing = unit(0, "lines"))

#### Print Plots to PDF ####
pdf("taxa_summary_plot_name.pdf", #name of file to print. can also include relative or absolute path before filename.
    , width = 16 # define plot width and height (this is in Inches). completely up to user.
    , height = 9
)
# plot commands for an individual plot or single multi-plot go here
dev.off()
