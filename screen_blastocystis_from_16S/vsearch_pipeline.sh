#!/bin/bash
# Quality filter and process raw fastq files, screen for eukaryotic taxa 

####################################################
#VSEARCH screen for eukaryotes from 16S V4 data
####################################################
# Uses only R1 to simplify

RAW=/Volumes/histolytica/wildAnimalBlast/raw #location where your raw gziped fastq files are
BACREF=/Users/mann/refDB/ezbiocloud/ezbiocloud_qiime_full.fasta #change this to your bacterial reference db taxonomy file
BACTAX=/Users/mann/refDB/ezbiocloud/ezbiocloud_id_taxonomy.txt
EUKTAX=/Users/mann/refDB/silva128/99_SILVA_128_taxa_map_7_level.w_blastocystis_entamoeba.txt #change this to your eukaryotic reference db taxonomy file
EUKREF=/Users/mann/refDB/silva128/99_SILVA_128_euks_rep_set.fasta #change this to your eukaryotic reference db fasta file
NCBIDB=/parfreylab/shared/databases/NCBI_NT/NT_06.08.2018/nt

#first need to trim primers, quality filter
cd $RAW 
mkdir trim 
#change below to match pattern in your R1 fastq
ls *_1.fastq.gz | sed 's/_1.fastq.gz//' | parallel 'cutadapt -a GTGYCAGCMGCCGCGGTAA --trim-n -m 50 -q 30,30 -o trim/{}.trim.fastq {}_1.fastq.gz' 1> trim/cutadapt.stats

#convert to fasta
cd trim
ls *fastq | sed 's/.fastq//' | parallel 'fastq_to_fasta -i {}.fastq -o {}.fasta -Q 33'
rm *fastq 

#add headers and concatenate, clean up intermediate files
ls *fasta | sed 's/.trim.fasta//' | while read line; do awk '/^>/{print ">'$line'."++i; next} {print}' < $line.trim.fasta > $line.headers.fasta; done
rm *trim.fasta 
cat *headers* > ../combinedseq.fa
rm *headers*
cd ..

#dereplicate and cluster to get denovo reference dataset (remove singletons)
vsearch --derep_fulllength combinedseq.fa -sizeout -relabel OTU -output uniques.fa
vsearch --cluster_unoise uniques.fa --sizein --sizeout --uc clustered.uc --centroids clustered.fa --minsize 2
gzip uniques.fa

#assign taxonomy using bacterial database
assign_taxonomy.py -i clustered.fa -o ../assigntax -t $BACTAX -r $BACREF
cd ..

#pull anything that wasn't confidently assiged to bacteria
grep -v "Bacteria;" assigntax/clustered_tax_assignments.txt | grep -v "Archaea;" | awk '{print $1}' > assigntax/query.ids
filter_fasta.py -f raw/clustered.fa -o query.fa -s assigntax/query.ids

#quick screen
#assign_taxonomy.py -i query.fa -o assigntax/ -t $EUKTAX -r $EUKREF

#get frequency table for all otus
vsearch -usearch_global raw/combinedseq.fa -db query.fa -id 0.99 -otutabout otu_counts.txt

#blast screen
blastn -evalue 1e-10 -max_target_seqs 100 -perc_identity 0.90 -db $NCBIDB -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -query query.fa -out blast.out