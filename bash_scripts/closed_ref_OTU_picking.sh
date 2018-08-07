#!/bin/bash

# Input:
# Will be specified below/while running

# A script for processing a NODE_rep.fasta file by Closed Reference OTU picking and taxonomy assignments

# Last edited: August 7, 2018
# authors: Tracy Wang, Evan Morien

# Input:
# paths to NODE_REPRESENTATIVES fasta file from MED, SILVA reference sequence set and taxonomy map, parameters file 
# Parameters file should include (please remove # when copying text over):
#pick_otus:similarity    1.0
#assign_taxonomy:similarity    1.0
#assign_taxonomy:uclust_max_accepts    1
#pick_otus:enable_rev_strand_match    True
# Output:
# A filtered taxa map with the successful taxonomy assignments
# Taxonomy assignments in separate folders. 
# Please view detailed outputs in log file

# Requirements: #links good as of July 2018
# QIIME or macqiime is installed http://qiime.org/

#--------------------------------------------------------------------------------------
### SECTION ONE: Filepaths and parameters ###
# This section prompts the user to specify necessary file paths and parameters
# All files made during this script are specified in the LOG file, as well as all commands that were run.

echo "This script will be conducting closed reference OTU picking and taxonomy assignments using SILVA, default parameters and the filtered fasta file from closed reference OTU picking.
      This script requires the following complete pathways:"
echo "1. NODE REPRESENTATIVES from MED"
echo "2. SILVA reference sequence set and taxonomy map"
echo "3. Parameters file including: 
pick_otus:similarity    1.0
assign_taxonomy:similarity    1.0
assign_taxonomy:uclust_max_accepts    1
pick_otus:enable_rev_strand_match    True
"
echo "Please make sure that all of the above conditions are met before running this script, as well as invoking macqiime"
	echo ""
	read -p "Press enter if you wish to continue. Press CTRL+C to exit script now"
	echo ""

# It also records their specifications in the log so they know exactly what they did
# I'm copying these variables into the LOG file so we know what we used in the past if we need to go back.

# Creating a log file
printf "Closed Reference OTU Picking LOG File - " > closed_reference_otu_picking.log
# appending date to log file
date >> closed_reference_otu_picking.log
echo "" "">> closed_reference_otu_picking.log
echo "Input information: complete pathways to data inputs" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log

echo "Please enter complete pathway to NODE REPRESENTATIES fasta file created by MED"
read path_noderep

echo "Node Representatives File: $path_noderep" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log

echo "Please enter complete pathway to parameters file"
read path_parameters

echo "Parameters File: $path_parameters" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log

echo "Please enter complete pathway to SILVA reference sequence set (unaligned rep set)"
read path_SILVAseq

echo "SILVA Reference Seqs unaligned: $path_SILVAseq" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log


echo "Please enter complete pathway to SILVA taxonomy map"
read path_SILVAmap

echo "SILVA taxonomy map: $path_SILVAmap" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log
#--------------------------------------------------------------------------------------
### SECTION TWO: Closed Reference OTU Picking ###
echo "" >> closed_reference_otu_picking.log
echo "Scripts that were run:" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log
echo "First section: Closed Reference OTU picking" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log

echo "SCRIPT HAS NOW STARTED TO RUN"
echo "Closed reference OTU picking..."
# Run pick_closed_reference_otus.py: allows parallelization and designates number of parallel threads to run, but will hang on certain machines. a known bug. leaving out parallelization.

pick_closed_reference_otus.py -i $path_noderep -r $path_SILVAseq -o closed_ref_OTU_picking/ -t $path_SILVAmap -p $path_parameters 
echo "pick_closed_reference_otus.py -i $path_noderep -r $path_SILVAseq -o closed_ref_OTU_picking/ -t $path_SILVAmap -p $path_parameters" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log

echo "Filtering against list of failures..."
# Filter input fasta files with list of "failures" from pick_closed_reference_otus.py's output
filter_fasta.py -f $path_noderep -o closed_ref_OTU_picking/NODE-REPRESENTATIVES.DOWNSTREAM.failures.fasta -s closed_ref_OTU_picking/uclust_ref_picked_otus/*_failures.txt 
echo "filter_fasta.py -f $path_noderep -o closed_ref_OTU_picking/NODE-REPRESENTATIVES.DOWNSTREAM.failures.fasta -s closed_ref_OTU_picking/uclust_ref_picked_otus/*_failures.txt" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log
#--------------------------------------------------------------------------------------
### SECTION THREE: Assigning taxonomy ###
# Assign taxonomy using SILVA database #
echo "Second section: Assigning taxonomy" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log

# assign_taxonomy: (using --similarity and --uclust_max_accepts parameters)
# Settings:
# 100% similarity, best hit only
# The rest on default run settings
echo "Taxonomy assignment now begins!"
echo "Assigning taxonomy using SILVA database..." 
assign_taxonomy.py -i closed_ref_OTU_picking/NODE-REPRESENTATIVES.DOWNSTREAM.failures.fasta -t $path_SILVAmap -r $path_SILVAseq -m uclust -o ./assign_taxonomy_sim100 --similarity 1.0 --uclust_max_accepts 1 
echo "assign_taxonomy.py -i closed_ref_OTU_picking/NODE-REPRESENTATIVES.DOWNSTREAM.failures.fasta -t $path_SILVAmap -r $path_SILVAseq -m uclust -o ./assign_taxonomy_sim100 --similarity 1.0 --uclust_max_accepts 1" "" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log

# Take the unassigned OTUs from the previous step and run them through a second iteration of taxonomy assignment with less stringent parameters.
echo "Running through second iteration of taxonomy assignment for unassigned OTUs..."
# We keep everything from taxa map except unassigned MED nodes:
grep -v "Unassigned" assign_taxonomy_sim100/NODE-REPRESENTATIVES.DOWNSTREAM.failures_tax_assignments.txt > assign_taxonomy_sim100/sim_100.txt 
echo "grep -v "Unassigned" assign_taxonomy_sim100/NODE-REPRESENTATIVES.DOWNSTREAM.failures_tax_assignments.txt > assign_taxonomy_sim100/sim_100.txt" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log
# We keep OTU IDs only for unassigned nodes:
grep "Unassigned" assign_taxonomy_sim100/NODE-REPRESENTATIVES.DOWNSTREAM.failures_tax_assignments.txt | awk '{print $1}' > assign_taxonomy_sim100/unassigned_100.txt 
echo ""Unassigned" assign_taxonomy_sim100/NODE-REPRESENTATIVES.DOWNSTREAM.failures_tax_assignments.txt | awk '{print $1}' > assign_taxonomy_sim100/unassigned_100.txt" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log
# Use list of unassigned OTU ID's to filter original fasta file:
filter_fasta.py -f closed_ref_OTU_picking/NODE-REPRESENTATIVES.DOWNSTREAM.failures.fasta -o closed_ref_OTU_picking/for_sim99_taxa_assignments.fasta -s assign_taxonomy_sim100/unassigned_100.txt 
echo "filter_fasta.py -f closed_ref_OTU_picking/NODE-REPRESENTATIVES.DOWNSTREAM.failures.fasta -o closed_ref_OTU_picking/for_sim99_taxa_assignments.fasta -s assign_taxonomy_sim100/unassigned_100.txt" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log

# Assign taxonomy using the filtered fasta file #
# Settings:
# 99% similarity, best hit only
echo "Assigning taxonomy using filtered fasta file..."
assign_taxonomy.py -i closed_ref_OTU_picking/for_sim99_taxa_assignments.fasta -t $path_SILVAmap -r $path_SILVAseq  -m uclust -o ./assign_taxonomy_sim99 --similarity 0.99 --uclust_max_accepts 1 
echo "assign_taxonomy.py -i closed_ref_OTU_picking/for_sim99_taxa_assignments.fasta -t $path_SILVAmap -r $path_SILVAseq  -m uclust -o ./assign_taxonomy_sim99 --similarity 0.99 --uclust_max_accepts 1 " >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log

# Take the unassigned OTUs from the previous step, use them to filter the original fasta file, and reassign taxonomy with default settings:
echo "Using unassigned OTUs to filter original fasta file and reassign taxonomy...."
grep -v "Unassigned" assign_taxonomy_sim99/for_sim99_taxa_assignments_tax_assignments.txt > assign_taxonomy_sim99/sim_99.txt 
grep "Unassigned" assign_taxonomy_sim99/for_sim99_taxa_assignments_tax_assignments.txt | awk '{print $1}' > assign_taxonomy_sim99/unassigned_99.txt 
filter_fasta.py -f closed_ref_OTU_picking/for_sim99_taxa_assignments.fasta -o closed_ref_OTU_picking/for_default_taxa_assignments.fasta -s assign_taxonomy_sim99/unassigned_99.txt 
echo "grep -v "Unassigned" assign_taxonomy_sim99/for_sim99_taxa_assignments_tax_assignments.txt > assign_taxonomy_sim99/sim_99.txt" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log
echo "grep "Unassigned" assign_taxonomy_sim99/for_sim99_taxa_assignments_tax_assignments.txt | awk '{print $1}' > assign_taxonomy_sim99/unassigned_99.txt" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log
echo "filter_fasta.py -f closed_ref_OTU_picking/for_sim99_taxa_assignments.fasta -o closed_ref_OTU_picking/for_default_taxa_assignments.fasta -s assign_taxonomy_sim99/unassigned_99.txt " >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log

# Assign taxonomy with default parameters #
# Settings:
# 97% similarity, top 3 hits ‘averaged’ to create a taxonomy string
echo "Assigning taxonomy using default parameters..."
assign_taxonomy.py -i closed_ref_OTU_picking/for_default_taxa_assignments.fasta -t $path_SILVAmap -r $path_SILVAseq -m uclust -o ./assign_taxonomy_defaults 
echo "assign_taxonomy.py -i closed_ref_OTU_picking/for_default_taxa_assignments.fasta -t $path_SILVAmap -r $path_SILVAseq -m uclust -o ./assign_taxonomy_defaults" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log

# Convert the closed-reference OTU table output from .biom to a tab-delimited .txt format #
echo "Converting OTU from .biom to .txt format..."
biom convert -i closed_ref_OTU_picking/otu_table.biom -o closed_ref_OTU_picking/otu_table.txt --to-tsv 
echo "biom convert -i closed_ref_OTU_picking/otu_table.biom -o closed_ref_OTU_picking/otu_table.txt --to-tsv" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log

# Fetch accessions assigned in the closed reference step, then use them to filter the SILVA taxa map:
# The file we get them from has two header lines, skipping them with -n flag
echo "Filtering SILVA taxa map...."
tail -n +3 closed_ref_OTU_picking/otu_table.txt | awk '{print $1}' > closed_ref_OTU_picking/SILVA_accessions.txt 
echo "tail -n +3 closed_ref_OTU_picking/otu_table.txt | awk '{print $1}' > closed_ref_OTU_picking/SILVA_accessions.txt" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log

# Filter taxa map using grep -f (-f stands for file. see grep manual for details... enter "man grep" in terminal window.) 
# NOTE: grep -f is useful for simple tasks but can take a long time if the input list is large (>100 or so entries)
# grep -F is for “fixed strings” and will speed up this grep command considerably.
echo "Filtering SILVA taxa map using grep -f..."
grep -F -f closed_ref_OTU_picking/SILVA_accessions.txt $path_SILVAmap > closed_ref_OTU_picking/inherited_taxonomies.txt 
echo "grep -F -f closed_ref_OTU_picking/SILVA_accessions.txt $path_SILVAmap > closed_ref_OTU_picking/inherited_taxonomies.txt" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log

# Combine the filtered taxa map with the successful taxonomy assignments from the last few steps:
echo "Combine filtered taxa map with taxonomy assignments..."
cat closed_ref_OTU_picking/inherited_taxonomies.txt assign_taxonomy_sim100/sim_100.txt assign_taxonomy_sim99/sim_99.txt assign_taxonomy_defaults/*.txt > closed_ref_OTU_picking/complete_taxonomy_assignments.txt 
echo "cat closed_ref_OTU_picking/inherited_taxonomies.txt assign_taxonomy_sim100/sim_100.txt assign_taxonomy_sim99/sim_99.txt assign_taxonomy_defaults/*.txt > closed_ref_OTU_picking/complete_taxonomy_assignments.txt" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log

# Finished script
echo "The script is now done running."
echo "DONE" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log

# List of output files:
echo "The following are outputs made from script:" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log
echo "closed_ref_OTU_picking folder:" >> closed_reference_otu_picking.log
echo "complete_taxonomy_assignments.txt
inherited_taxonomies.txt
log_20180803171038.txt
NODE-REPRESENTATIVES.DOWNSTREAM.failures.fasta
otu_table.biom
otu_table.txt
SILVA_accessions.txt
uclust_ref_picked_otus (folder)" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log
echo "closed_ref_OTU_picking/uclust_ref_picked_otus folder:" >> closed_reference_otu_picking.log
echo "NODE-REPRESENTATIVES.DOWNSTREAM.no_chimeras_clusters.uc
NODE-REPRESENTATIVES.DOWNSTREAM.no_chimeras_failures.txt
NODE-REPRESENTATIVES.DOWNSTREAM.no_chimeras_otus.log
NODE-REPRESENTATIVES.DOWNSTREAM.no_chimeras_otus.txt" >> closed_reference_otu_picking.log
echo "" "">> closed_reference_otu_picking.log
echo "assign_taxonomy_sim100 folder:" >> closed_reference_otu_picking.log
echo "NODE-REPRESENTATIVES.DOWNSTREAM.failures_tax_assignments.log
NODE-REPRESENTATIVES.DOWNSTREAM.failures_tax_assignments.txt
sim_100.txt
unassigned_100.txt" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log
echo "assign_taxonomy_sim99 folder:" >> closed_reference_otu_picking.log
echo "for_sim99_taxa_assignments_tax_assignments.log
for_sim99_taxa_assignments_tax_assignments.txt
sim_99.txt
unassigned_99.txt" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log
echo "assign_taxonomy_defaults folder" >> closed_reference_otu_picking.log
echo "for_default_taxa_assignments_tax_assignments.log
for_default_taxa_assignments_tax_assignments.txt" >> closed_reference_otu_picking.log
echo "" >> closed_reference_otu_picking.log
echo "Log File is complete" >> closed_reference_otu_picking.log
