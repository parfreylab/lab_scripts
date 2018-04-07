#### calculate relative abundance for species ####
ra.species <- project_data %>%
  merge_samples("Rank7") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  t() %>%
  otu_table() %>%
  as.data.frame()
colnames(ra.species) <- paste("relative_abundance", colnames(ra.species), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables

#### calculate prevalence for species ####
#prevalence is the percentage of samples that an OTU shows up in (compared to the total number of samples). threshold for presence 0.001 relative abundance.
prevalence = function(x){ #this only returns prevalence counts per OTU, which can be summed by "merge_samples" function, but need to be expressed as a percentage only for samples belonging to a particular species
  rel <- x/sum(x)
  x[which(rel < 0.001)] <- 0
  x[x >= 2] <- 1
  return(x)
}

allones = function(x){ #this is just a trick i used to make a matrix identical to the prevalence counts matrix, but just counting how many samples belong to each species
  x[x >= 0] <- 1
  return(x)
}

prev_counts.species <- project_data %>% #this produces prevalence "counts" for each species, but not percentages
  transform_sample_counts(fun = prevalence) %>%
  merge_samples("Rank7") %>%
  t() %>%
  otu_table() %>%
  as.data.frame()
colnames(prev_counts.species) <- paste("prevalence", colnames(prev_counts.species), sep = ".") #add something to column headers to distinguish between relative abundance and prevalence tables

prev_possible.species <- project_data %>% #this produces a maximum possible prevalence count per species per OTU
  transform_sample_counts(fun = allones) %>%
  merge_samples("Rank7") %>%
  t() %>%
  otu_table() %>%
  as.data.frame()
colnames(prev_possible.species) <- paste("prevalence", colnames(prev_possible.species), sep = ".") #add something to distinguish between relative abundance and prevalence

#dividing the first matrix by the second will give us prevalence expressed as a percentage on a per-species basis
prev.species = (prev_counts.species/prev_possible.species)*100


#### creating and formatting prevalence/rel. abundance table ####
#merge data frames
ra.prev.species <- transform(merge(ra.species,prev.species,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

#add taxonomy strings
taxa.tomerge <- as.data.frame(tax_table(project_data))
taxa.tomerge <- taxa.tomerge %>% unite("Taxonomy", Rank1, Rank2, Rank3, Rank4, Rank5, Rank6, Rank7, sep="; ")

#merge data frames
new_prevalence_table <- transform(merge(ra.prev.species,taxa.tomerge,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)

#write otu table to file
write.table(new_prevalence_table, file="rel_ab.prev.species_collapsed.txt", quote=F, row.names=T, col.names=T, sep="\t")
#modify this with nano to insert tab at beginning of header line

