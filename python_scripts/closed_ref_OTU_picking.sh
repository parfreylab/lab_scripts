#assign taxonomy at most stringent
assign_taxonomy.py -i closed_ref_OTU_picking/NODE-REPRESENTATIVES.DOWNSTREAM.failures.fasta -t ../databases/99_SILVA_128_taxa_map_7_level.w_blastocystis_entamoeba.txt -r ../databases/99_SILVA_128_euks_rep_set.fasta -m uclust -o ./assign_taxonomy_sim100 --similarity 1.0 --uclust_max_accepts 1
grep -v "Unassigned" assign_taxonomy_sim100/NODE-REPRESENTATIVES.DOWNSTREAM.failures_tax_assignments.txt > assign_taxonomy_sim100/sim_100.txt 
grep "Unassigned" assign_taxonomy_sim100/NODE-REPRESENTATIVES.DOWNSTREAM.failures_tax_assignments.txt | awk '{print $1}' > assign_taxonomy_sim100/unassigned_100.txt 
filter_fasta.py -f closed_ref_OTU_picking/NODE-REPRESENTATIVES.DOWNSTREAM.failures.fasta -o closed_ref_OTU_picking/for_sim99_taxa_assignments.fasta -s assign_taxonomy_sim100/unassigned_100.txt

#assign taxonomy at 99% similarity
assign_taxonomy.py -i closed_ref_OTU_picking/for_sim99_taxa_assignments.fasta -t ../databases/99_SILVA_128_taxa_map_7_level.w_blastocystis_entamoeba.txt -r ../databases/99_SILVA_128_euks_rep_set.fasta -m uclust -o ./assign_taxonomy_sim99 --similarity 0.99 --uclust_max_accepts 1
grep -v "Unassigned" assign_taxonomy_sim99/for_sim99_taxa_assignments_tax_assignments.txt > assign_taxonomy_sim99/sim_99.txt
grep "Unassigned" assign_taxonomy_sim99/for_sim99_taxa_assignments_tax_assignments.txt | awk '{print $1}' > assign_taxonomy_sim99/unassigned_99.txt
filter_fasta.py -f closed_ref_OTU_picking/for_sim99_taxa_assignments.fasta -o closed_ref_OTU_picking/for_default_taxa_assignments.fasta -s assign_taxonomy_sim99/unassigned_99.txt

#assign taxonomy at default similarity
assign_taxonomy.py -i closed_ref_OTU_picking/for_default_taxa_assignments.fasta -t ../databases/99_SILVA_128_taxa_map_7_level.w_blastocystis_entamoeba.txt -r ../databases/99_SILVA_128_euks_rep_set.fasta -m uclust -o ./assign_taxonomy_defaults

#convert closed ref biom output to text, extract IDs, filter taxa map for inherited taxonomies
biom convert -i closed_ref_OTU_picking/otu_table.biom -o closed_ref_OTU_picking/otu_table.txt --to-tsv
tail -n +3 closed_ref_OTU_picking/otu_table.txt | awk '{print $1}' > closed_ref_OTU_picking/SILVA_accessions.txt #fetch accessions #the file we get them from has two header lines, skipping them with -n flag
grep -f closed_ref_OTU_picking/SILVA_accessions.txt ../databases/99_SILVA_128_taxa_map_7_level.w_blastocystis_entamoeba.txt > closed_ref_OTU_picking/inherited_taxonomies.txt

#now combine above files into complete taxonomy assignment file
cat closed_ref_OTU_picking/inherited_taxonomies.txt assign_taxonomy_sim100/sim_100.txt assign_taxonomy_sim99/sim_99.txt assign_taxonomy_defaults/*.txt > closed_ref_OTU_picking/complete_taxonomy_assignments.txt
