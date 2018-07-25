#basic procedure for calculating relative abundance and prevalence according to a metadata category
#NOTE: THIS SCRIPT IS MEANT TO BE CHANGED TO FIT YOUR DATA. PLEASE REMEMBER TO:
#       1. MAKE A COPY OF THIS SCRIPT INTO YOUR WORKING DIRECTORY
#       2. CHANGE ALL GENERALIZED PARAMETERS/VARIABLES (EX: "FACTOR") TO MATCH YOUR DATA
#       3. THIS CODE SUPPLEMENTS THE OTHER SCRIPTS ON THIS REPO, AND THE "project_data" OBJECT IS ONE OF CLASS PHYLOSEQ CREATED WITH, FOR EXAMPLE, THE ALPHA/BETA DIV R SCRIPT

#AUTHOR:        Evan Morien
#last modified: April 10th, 2018

#### calculate relative abundance for FACTOR ####
ra.FACTOR <- project_data %>%
  merge_samples("FACTOR") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  t() %>%
  otu_table() %>%
  as.data.frame()
colnames(ra.FACTOR) <- paste("relative_abundance", colnames(ra.FACTOR), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables

#### calculate prevalence for FACTOR ####
#prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). threshold for presence 0.001 relative abundance.
prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function, but need to be expressed as a percentage only for samples belonging to a particular FACTOR
  rel <- x/sum(x)
  x[which(rel < 0.001)] <- 0
  x[x >= 2] <- 1
  return(x)
}

allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each FACTOR
  x[x >= 0] <- 1
  return(x)
}

prev_counts.FACTOR <- project_data %>% #this produces prevalence "counts" for each FACTOR, but not percentages
  transform_sample_counts(fun = prevalence) %>%
  merge_samples("FACTOR") %>%
  t() %>%
  otu_table() %>%
  as.data.frame()
colnames(prev_counts.FACTOR) <- paste("prevalence", colnames(prev_counts.FACTOR), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables

prev_possible.FACTOR <- project_data %>% #this produces a maximum possible prevalence count per FACTOR per OTU
  transform_sample_counts(fun = allones) %>%
  merge_samples("FACTOR") %>%
  t() %>%
  otu_table() %>%
  as.data.frame()
colnames(prev_possible.FACTOR) <- paste("prevalence", colnames(prev_possible.FACTOR), sep = ".") #add something to distinguish between relative abundance and prevalence

#dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-FACTOR basis
prev.FACTOR = (prev_counts.FACTOR/prev_possible.FACTOR)*100


#### creating and formatting prevalence/rel. abundance table ####
#merge data frames
ra.prev.FACTOR <- transform(merge(ra.FACTOR,prev.FACTOR,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

#add taxonomy strings
taxa.tomerge <- as.data.frame(tax_table(project_data))
taxa.tomerge <- taxa.tomerge %>% unite("Taxonomy", Rank1, Rank2, Rank3, Rank4, Rank5, Rank6, Rank7, sep="; ") #remember to make this match the column labels of your own taxonomy

#merge data frames
new_prevalence_table <- transform(merge(ra.prev.FACTOR,taxa.tomerge,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

#write otu table to file
write.table(new_prevalence_table, file="rel_ab.prev.FACTOR_collapsed.txt", quote=F, row.names=T, col.names=T, sep="\t")
#modify this with nano to insert tab at beginning of header line

