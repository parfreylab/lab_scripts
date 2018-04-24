#### PermANOVA and Beta Dispersion Tests ####
#   NOTES:
#   1) if you rarefied your data to do your alpha and beta diversity analysies, you need to do this analysis on the same rarefied data
#   2) this example uses a filtered and rarefied phyloseq object as a starting point

#   NOTE: THIS SCRIPT IS MEANT TO BE CHANGED TO FIT YOUR DATA. PLEASE REMEMBER TO:
#       1. MAKE A COPY OF THIS SCRIPT IF YOU WISH TO MODIFY IT
#       2. CHANGE ALL GENERALIZED PARAMETERS/VARIABLES (EX: "FACTOR_1") TO MATCH YOUR DATA



library(phyloseq)
library(vegan)

#### analysis: PERMANOVA to see which factors contribute significant to the variation observed in the dataset ####
#using the vegan "adonis" function, we can calculate relationship of each factor with variance w/in the data
#need a distance matrix and the sample data (metadata) to start
project_bray.rarefied <- phyloseq::distance(project_data.rarefied, method = "bray")
sample_df <- data.frame(sample_data(project_data.rarefied))

#now the adonis test for things we can test for the whole experiment
res.adonis.rarefied <- adonis(project_bray.rarefied ~ FACTOR1 + FACTOR2 + FACTOR3 + etc..., data=sample_df, method="bray")
res.adonis.rarefied #run for full description of results
summary(res.adonis.rarefied) #run for results summary (this is less informative if I remember correctly)

#### analysis: beta dispersion test ####
beta.FACTOR1 <- betadisper(project_bray.rarefied, sample_df$FACTOR1) #testing differences between group centroids w/in specified factor #see documentation for additional parameters that can be adjusted
permutest(beta.FACTOR1) #do permutation test of above #this is more customizable, see documentation
beta.FACTOR2 <- betadisper(project_bray.rarefied, sample_df$FACTOR2)
permutest(beta.FACTOR2)
beta.FACTOR3 <- betadisper(project_bray.rarefied, sample_df$FACTOR3)
permutest(beta.FACTOR3)
