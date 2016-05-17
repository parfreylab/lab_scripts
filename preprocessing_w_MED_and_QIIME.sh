#!/bin/bash


# Input:
# Single seq file for all sequences

# a script for processing raw 16s or 18s sequence data in fastq format into an OTU table, phylogenetic tree, and representative sequences
# V 1.0 May 5th 2016
# authors: Melissa Chen, Evan Morien

# Input:
# Folder containing sequence files in .fastq or .fastq.gz format, path of mapping file, paths of database files, parameters for MED and QIIME to process sequence file
# Output:
# filtered OTU table with taxonomy added, representative sequences for each OTU, phylogenetic tree of sequences, various log files, other intermediate files from MED and QIIME

# Requirements: #links good as of May 2016
# QIIME or macqiime is installed http://qiime.org/
# fastx toolkit is installed http://hannonlab.cshl.edu/fastx_toolkit/
# MED is installed http://merenlab.org/2014/08/16/installing-the-oligotyping-pipeline/
# folder with raw fastq files must not contain any non-fastq files

#------------------------
# SECTION ONE: Filepaths and parameters
# This section prompts the user to specify necessary file paths and parameters
# It also records their specifications in the log so they know exactly what they did
# Last thing it does is count the seqs in each file prior to processing.
#------------------------
# Making Project directory: script will prompt user to type in a folder name and press enter. Then, we move inside the folder.
# All files will be created inside this directory. This means I don't need complete filepaths in order to do all commands.
# First, there will be a series of questions to ask. These answers will specify the parameters you want.

macqiime

echo "Here is a list of the info you will need to enter. Collect all this and have it ready before running this script. If you don't have it, you can press ctrl+c to kill the script now."
echo "1. a name for your project (a directory named with the input here will be created in the folder you run the script in)"
echo "2. is the data 16s or 18s"
echo "3. complete pathway to fastq files (just the folder they are in. they don't have to be unzipped if the aren't already)"
echo "4. complete filepath for properly formatted mapping file"
echo "5. complete filepaths for the reference database rep set fasta, map of ref tree, and aligned rep set"
echo "6. the trimming length for fastx toolkit (we ususally use 250)"
echo "7. the alignment gap filtering parameter: a value from 0.0 to 1.0, representing the percentage of gaps that are required to filter out an alignment position. recommend 0.9 for 'safe bet'. please see filter_alignment.py documentation for more details."
echo "8. the alignment entropy filtering paramter. a value from 0.0 to 1.0 representing the N/1.0 percentage of positions to be removed based on entropy. recommend 0.05 for 'safe bet'. please see filter_alignment.py documentation for more details."
echo "9. the tree building method. Valid choices are: clustalw, raxml_v730, muscle, fasttree, clearcut"
echo "10. the minimum number of reads per OTU to include in final OTU table"
echo "11. INPUT PROMPT COMES LATER: minimum substantive abundance (numeric. the -M paramteter for MED. a good baseline is 100, but your dataset may require a lower number. see MED documentation for additional details.)"

echo ""
echo ""
echo "ARE YOU OPERATING IN MACQIIME? If not, you are about to waste a lot of time waiting. Do it now if you haven't."
	echo ""
	read -p "Press enter if you wish to continue. Press CTRL+C to exit script now"
	echo ""

echo "Please enter project name"
read projectname
mkdir "$projectname"
cd "$projectname"
# Now, we are in the project directory.

echo "Is this 16s or 18s? [type 16 or 18]"
read projecttype

# Copying raw data from folder and putting it in $projectname folder.
# I wanted to copy it because then you can't mess up your raw data. Just in case!
# cat is used so I'm not tampering with original file.

echo "Enter complete pathway to raw fastq files."
read fastqlocation

	
# You should end up with a folder called "raw_fastq" inside your project directory.
# Next, we will save other complete file paths for the database as variables to be used later.
# I used the SILVA111 database for this, but other databases can be used as well.
# I'm also copying these variables into the LOG file so we know what we used in the past if we need to go back.

echo "Enter complete pathway to mapping file (txt)"
read mappingfile

echo "Mappingfile: $mappingfile" >> LOG
echo "" >> LOG
	
	
echo "Enter complete pathway to database rep set (fasta)"
read reftreefasta
echo "Enter complete pathway to reference taxonomy (txt)"
read reftreemap
echo "Enter complete pathway to aligned rep set (fasta)"
read reftreealigned
	
	
	
echo "" >> LOG
echo "Input information: complete pathways to reference database" >> LOG

echo "Reftreefasta: $reftreefasta" >> LOG
echo "Reftreemap: $reftreemap" >> LOG
echo "Reftreealigned: $reftreealigned" >> LOG
echo "" >> LOG

# Now, we enter some parameters and variables. There is a prompt first, and then the user may type in their conplete filepath.

	
# Enter the trimming length for fastx_clipper and fastx_trimmer.
echo "Enter trimming length for reads"
read trimlength

#enter the filter_alignment.py paramters. see documentation on QIIME website for more details.
echo "Enter the gap filtering paramter (0.0 to 1.0). recommend 0.90"
read gap

echo "Enter the entropy filtering paramter (0.0 to 1.0). recommend 0.05"
read entropy
	
# Enter Tree building method. Fasttree is fast and the default; raxml seems to be more popular. Although, I've read somewhere they're about the same?
echo "Enter tree building method. Method for tree building. Valid choices are: clustalw, raxml_v730, muscle, fasttree, clearcut"
read treemethod
	
# Enter the minimum count. This is for filtering the OTU table.
echo "Enter minimum number of reads per OTU to include in final OTU table"
read minimumcount
echo "trimlength: $trimlength" >> LOG
echo "minimum entropy: $minimumentropy" >> LOG
echo "Tree method: $treemethod" >> LOG
echo "Minimumcountfraction: $minimumcount" >> LOG
echo "" >> LOG

# COPYING files over and unzipping them
echo "Copying fasta files..."
mkdir ./seq_fasta/

echo "Do not be alarmed if it tells you there is 'no such file or directory'-- this is normal!"
echo "Copying fasta files..."
cp "$fastqlocation"/*fq* ./seq_fasta/
cp "$fastqlocation"/*fastq* ./seq_fasta/
echo "Done copying"

# TODO: wrap each copy command in if statement that only executes if file in question exists
	
echo "Unzipping if necessary"
for f in ./seq_fasta/*.gz; do
    echo "unzipping file: $f ... "
    gunzip "$f"
done

# Printing preamble for LOG. The LOG file will contain intermediate summaries for use of troubleshooting.
echo "" >> LOG
echo "Log and script summaries for $projectname" >> LOG
echo "" >> LOG
	
# Note that all other parameters are default.

#------------------------

# Run multiple_split_libraries.py using default settings
# split_libraries default is:
# 	max_bad_run_length:3
# 	min_per_read_length_fraction: 0.75
#	sequence_max_n:0
# etc
# we set the phred quality threshold to 19, if you want to change it, you can do so below. this shouldn't be necessary, 19 is pretty standard.
# The sample ID indicator means that the stuff before this part will act as the new ID name.

echo "Running split_libraries QIIME script..."

if [ "$projecttype" == 18 ]
then
    touch parameters_mult_split_libraries.txt
    echo "split_libraries_fastq:phred_quality_threshold 19" >> parameters_mult_split_libraries.txt
    multiple_split_libraries_fastq.py -i ./seq_fasta -o multiple_split_libraries_fastq --demultiplexing_method sampleid_by_file -p parameters_mult_split_libraries.txt
    # Output should yield histograms.txt; log; seqs.fna; split_library_log.txt
    echo "multiple_split_libraries script finished..."
    # Still in 'home' directory
fi

###############

if [ "$projecttype" == 16 ]
then
	cd ./seq_fasta/
    ourindex=`ls | grep -v *R1* | grep -v *R2*` #get the name of the index file
    echo "preparing sequences with QIIME..."	
    split_libraries_fastq.py -i *R1* -o ../multiple_split_libraries_fastq --barcode_read_fps "$ourindex" --mapping_fps "$mappingfile" --phred_quality_threshold 19 --barcode_type 12
    echo "split_libraries_fastq script finished..."
    cd ..
fi

# Trimming files to $trimlength bp using fastx.

mkdir Trimmed_Quality_Filtered_MSL
cd multiple_split_libraries_fastq

# Using the $trimlength inputted before, use fastx trimmer to trip the seq file from multiple split libraries.
echo "Trimming and clipping seqeunces using fastx..."
fastx_trimmer -l "$trimlength" -i seqs.fna | fastx_clipper -l "$trimlength" -o ../Trimmed_Quality_Filtered_MSL/seqs_trim_clip.fna
echo "Trimming/Clipping Complete"
echo "" >> ../LOG

# Counting sequences
echo "POST-FASTX SEQUENCE COUNT" >> ../LOG
echo "seqs_trim_clip.fna sequence count" | tee -a ../LOG
grep -c ">" ../Trimmed_Quality_Filtered_MSL/seqs_trim_clip.fna | tee -a ../LOG
echo "" >> ../LOG
#back to project home dir
cd ..
#------------------------

# In home directory.
# Now, we want to alter the seq.fna file so that the format fits with MED requirements.
# MED likes it to be SampleXXX_ReadXXX, so this is what we'll do.
# Note that this ONLY works for halifax data so far-- you'll likely have to change this for other formats. 
# The halifax data gets us fasta headers look like this:
# >JulyZosInnerH_1972165 M02352:23:000000000-ANBEB:1:2119:15710:25291 1:N:0:0 orig_bc=ACCTACTTGTCT new_bc=ACCTACTTGTCT bc_diffs=0
# CCGCGGTAATTCCAGCTCCAATAGCGTATATTAAAGTTGTTGTGGTTAAAAAGCTCGTAGTTGGATCTCAACAGCCTTGAAGCGGTTAACTTTGTAGTTTTACTGCTTTATAGGCTGTGCTTTTGCTGGAGGTTTAGGGGTGCTCTTGATCGAGTGTCTCTGGATGCCAGCAGGTTTACTTTGAAAAAATTAGAGTGCTCAAAGCAGGCTAAAATGCCTGAATATTCGTGCATGGAATAATAGAATAGGA
# where JulyZosInnerH is the sample name and 1972165 is the read number
# In order to prep the sequence headers for MED, we want to get rid of everything after the first space, and add "Read" before the read number

cd ./Trimmed_Quality_Filtered_MSL
echo "Prepping fasta headers for MED..."
cat seqs_trim_clip.fna | awk -F' ' '{ st = index($0," ");print $1}' | awk -F_ '{ print ($2!="" ? $1"_Read"$2 : $1) }' > seqs_MED.fna
echo "done!"
# TODO: the "Read" part of the above string is actually not mandatory, but this wasn't known at the time the script was written. Currently this has no effect on the output, so it's not a high priority. Should be fixed eventually.
# the command above splits each line on " " and then returns the first element only, then splits that first element on _ and returns both but adds "_Read" between them. 
# it doesn't modify the sequence lines since they don't satisfy any of the splitting criteria, and the insertion in the second part of the command is conditional on a second element existing

#checking output fasta for correctly formatted headers and dumping output to LOG file for review
echo "Recording sample and read counts in LOG file..."
echo "SAMPLE AND READ COUNTS FOR MED" >> ../LOG
o-get-sample-info-from-fasta seqs_MED.fna >> ../LOG
echo "" >> ../LOG
echo "done!"

# Cleaning up and moving output to its own directory
echo "Cleaning up workspace..."
rm seqs_trim_clip.fna
mkdir ../MED
mv seqs_MED.fna ../MED/seqs_MED.fna

cd ..
# In the home directory

#------------------------

# SECTION THREE: MED decompose and transposing
# This section decomposes the fasta files using the $minimumentropy input
# Also, transposes the resulting MATRIX-COUNT.txt file into a QIIME-friendly format
# The last part changes the format of NOE-REPRESENTATIVES.fasta to be compatible with QIIME
# (explained below)

#------------------------


# Enter the minimum entropy. This is the only variable in this script that changes for MED (Minimum Entropy Decomposition)

echo "Enter Minimum Substantive Abundance (for MED, see MED documentation for additional details) -- the default is total reads divided by 5000."
read minimumentropy
echo "minimum substantive abundance: $minimumentropy" >> LOG
echo "" >> LOG


# Now, decompose fasta file. This is a single line.

echo "MED is decomposing fasta file..."

decompose ./MED/seqs_MED.fna -M "$minimumentropy" --gen-html -o ./MED/decompose-$minimumentropy

echo "MED run completed!"

# After decomposition, we need to transpose the MATRIX-COUNT.txt file because it is in the wrong format.
# I was too lazy to write my own script for transposing so I copied one from stackoverflow and tested it before including it below
# This part of the script was copy and pasted from stackoverflow

echo "Transposing MATRIX_COUNT.txt..."
awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j];
        }
        print str
    }
}' ./MED/decompose-$minimumentropy/MATRIX-COUNT.txt > MATRIX-COUNT_tmp.txt

# A few more alterations-- change the top left name from samples to #OTU ID, removing leading zeros from OTU IDs

sed 's/samples/#OTU ID/g' MATRIX-COUNT_tmp.txt | sed 's/^0*//' > MATRIX-COUNT_transposed.txt

rm MATRIX-COUNT_tmp.txt # Deletes intermediate file
rm ./MED/seqs_MED.fna

# Also, remove the pesky "|size:XXX" at the end of every OTU ID in the NODE_REPRESENTATIVES.fasta file.
# This is necessary because:
# 1. Adding metadata later on requires that the observation_metadata doesn't have the "|size:XXX"
# 2. The tree cannot have "|size:XXX" otherwise the tree tips will not match up with the OTU table itself.

echo "Modifying MED fasta headers for use with QIIME..."
sed 's/|size:[0-9]*$//g' ./MED/decompose-$minimumentropy/NODE-REPRESENTATIVES.fasta | sed 's/^>0*/>/' > NODE_REP_FORDOWNSTREAM.fasta
echo "done!"
#------------------------

# SECTION FOUR: Making a tree and OTU Table
# This last section goes back to QIIME to make a phylotree and the OTU table
# You can specify which tree you want
# It also records some of the intermediate results in the LOG
# Finally, it adds observation metadata to the file.
# The resulting OTU table and tree have been tested to be compatible for: alpha_rarefaction,beta_diversity_through_plots,summarize_taxa_through_plots 

#------------------------

# Assign Taxonomy using previous filepaths of references database. Need to use the unaligned database.
echo "QIIME STEPS BEGIN HERE"
echo "Assigning taxonomy..."
assign_taxonomy.py -i NODE_REP_FORDOWNSTREAM.fasta -o assign_taxonomy -r "$reftreefasta" -t "$reftreemap"

# Align sequences. This step uses the aligned reference database you supplied.

echo "Aligning sequences..."
align_seqs.py -i NODE_REP_FORDOWNSTREAM.fasta -t "$reftreealigned" -o aligned_seqs

# Now, writing the newly aligned sequences and failures to the log for easy access later.

echo "" >> LOG
echo "Aligned failures: " >> LOG
less ./aligned_seqs/*failures.fasta >> LOG
echo "" >> LOG

# Filtering alignment

echo "Filtering alignment..."
filter_alignment.py -i ./aligned_seqs/*aligned.fasta -s -e $entropy -g $gap -o "filter_alignment_G${gap}_E${entropy}"

# Make phylogenetic tree

echo "Making phylogenetic tree..."
make_phylogeny.py -i ./filter_alignment_G${gap}_E${entropy}/*.fasta -o makephylo_fastree.tre -t "$treemethod"

# Make OTU table

echo "Making OTU Table..."
biom convert -i ./MATRIX-COUNT_transposed.txt -o OTU_Table.biom --to-json --table-type="OTU table"

echo "####################" >> LOG
echo "OTU Table Summary:" >> LOG
echo "####################" >> LOG
biom summarize-table -i OTU_Table.biom >> LOG
echo "" >> LOG

# Filter OTUs by removing singles and OTUs that have very few reads

echo "Filtering OTUs..."
filter_otus_from_otu_table.py -i OTU_Table.biom -o removed_singles_few_reads.biom --min_count "$minimumcount"

echo "" >> LOG
echo "OTU Table Summary After Filtering:" >> LOG
biom summarize-table -i removed_singles_few_reads.biom >> LOG
echo "" >> LOG

echo "Adding obs metadata..."

biom add-metadata -i removed_singles_few_reads.biom -o OTU_Table_wtaxa.biom --observation-metadata-fp ./assign_taxonomy/*.txt --observation-header OTUID,taxonomy,confidence --sc-separated taxonomy

biom summarize-table -i ./OTU_Table_wtaxa.biom | tee -a LOG

echo "DONE!"
echo "DONE" >> LOG
