#### Basic framework for making a tree + heatmap plot with ggtree, starting with a phyloseq object (or a filtered tree) and a set of metadata for the tree ####
#### Libraries ####
library(tidyverse) #for ggplot2, and for using the %>% pipe
library(phyloseq) #if you want to start with a phyloseq object
library(ggtree) #for the plots we're making

#author: evan morien
#revision: Nov 29th, 2017
#more info on different types of plots you can make with ggtree can be found in the documentation section here: https://bioconductor.org/packages/release/bioc/html/ggtree.html
#and specificaly, here: https://bioconductor.org/packages/release/bioc/vignettes/ggtree/inst/doc/advanceTreeAnnotation.html

#NOTE: these figures don't tend to display very accurately in R Studio, so I print them to PDF and adjust by reviewing the PDF on another screen
#some details: a large figure and a small figure will both require adjustments of several parameters, namely...
#1. the size of the plot defined up front (height is most crucial)
#2. the scale_y_continuous function (look this one up, but briefly: adjusting the smaller number will change how big the upper and lower margins are outside the actual tree plot. it's the way to leave room for labels at the bottom/top of the heatmap columns.)
#3. offset: how much spacing to put between the tree tips and the heatmap
#4. width: how wide the heatmap is in relation to the tree
#5. align: i find it hard to read plots with unaligned tip labels in relation to the heatmap.

#NOTE: i used a phyloseq object as my starting point because it was easy to filter it based on the taxonomies I had in one of my metadata columns, there's no need to start with anything but a tree and metadata, though
#where you see the "phy_tree(phyloseq_object)" part of the line where ggtree() is invoked, you can just replace that with an object that stores a single tree that corresponds to your metadata

#### code for a small-ish figure (clade of 20-60 members) ####
dd <- read.table(file="prev.rel_ab.infected_vs_weight_change.wtaxa.plus_deseq_results.rat_colitis.oct_2016.txt", sep="\t", stringsAsFactor=F, header=T, row.names=1) #read in metadata
dd <- dd[grepl("Mollicutes", dd[["Taxonomy"]]), ] #OPTIONAL: filter metadata if necessary
project_data.mollicutes <- prune_taxa(row.names(dd), project_data) #this is a phyloseq function
#get tree from phyloseq object and plot
pdf("taxonomy_deseq_plots.mollicutes.pdf", #name of file to print. can also include relative or absolute path before filename.
    width = 15, # define plot width and height. completely up to user.
    height = 12)
p <- ggtree(phy_tree(project_data.mollicutes), ladderize=TRUE, right=TRUE) + geom_tiplab(size=2, align=TRUE, linetype="dotted", linesize=.05) + 
  theme_tree2() + ggtitle("Mollicutes log2FoldChange Values") + theme(plot.title = element_text(hjust = 0.5, size=20, face="bold"))
pp <- (p + scale_y_continuous(expand=c(0.25, 2))) %>%
  gheatmap(dd[,1:8], offset=0.05, width=0.5, colnames=T, colnames_position="bottom", colnames_angle=45, low="black", high="pink", color="lightgrey", hjust = 1) #note i am plotting a specific subset of my metadata here. adjust as needed.
pp + theme(legend.position="right")
dev.off()

#### code for a really large figure (clade of 500+ members) ####
dd <- read.table(file="prev.rel_ab.infected_vs_weight_change.wtaxa.plus_deseq_results.rat_colitis.oct_2016.txt", sep="\t", stringsAsFactor=F, header=T, row.names=1)
dd <- dd[grepl("Clostridia", dd[["Taxonomy"]]), ]
project_data.clostridiales <- prune_taxa(row.names(dd), project_data) #this is a phyloseq function
#get tree from phyloseq object and plot
pdf("taxonomy_deseq_plots.clostridia.pdf", #name of file to print. can also include relative or absolute path before filename.
    width = 15, # define plot width and height. completely up to user.
    height = 90)
p <- ggtree(phy_tree(project_data.clostridiales), ladderize=TRUE, right=TRUE) + geom_tiplab(size=2, align=TRUE, linetype="dotted", linesize=.05) + 
  theme_tree2() + ggtitle("Clostridia log2FoldChange Values") + theme(plot.title = element_text(hjust = 0.5, size=20, face="bold"))
pp <- (p + scale_y_continuous(expand=c(0.025, 2))) %>%
  gheatmap(dd[,1:8], offset=0.11, width=0.2, colnames=T, colnames_position="bottom", colnames_angle=45, low="black", high="pink", color="lightgrey", hjust = 1)
pp + theme(legend.position="right")
dev.off()