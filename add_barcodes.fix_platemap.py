import re
import argparse

####a script to create an analysis-ready mapping file and platemap with descriptive names from 96 well plate metadata files (metadata file, plate map, barcode plate spreadsheet -> mapping file w metadata & barcodes, plate map with descriptive names)####
##IMPORTANT: all inputs should be tab separated text files
##inputs:
###1. visual layout of 96 well plate (grid format). user must make this file. (plate map, has swab IDs in place of sample names)
###2. spreadsheet with barcodes in each primer plate. each plate has a unique barcode sequence
###3. metadata file that links swab IDs to sample names, plus other metadata.

#########
##USAGE##
#########
#python scripts/add_barcocdes.fix_platemap.py plate_map_with_swab_ids.txt barcode_platemap_spreadsheet.txt metadata_file_no_barcodes.txt
#test example: python add_barcodes.fix_platemap.py -g Jordan_Hamden_Mouse_Cortisol_16s_plate_map.txt -b 515f_806r_illumina_primers_515barcoded.txt -m Jordan_Hamden_Mouse_Cortisol_16s.txt
#IMPORTANT: see README for rules on file formatting

#PERL TEST--
#perl add_barcodes.fix_platemap.pl Jordan_Hamden_Mouse_Cortisol_16s_plate_map.txt 515f_806r_illumina_primers_515barcoded.txt Jordan_Hamden_Mouse_Cortisol_16s.txt

print "run: python add_barcodes.fix_platemap.py -h for help."
print "\n\n"

# Need to have better descriptions/help !!!
parser = argparse.ArgumentParser(
	description = "Mapping file with barcodes & metadata, and platemap with descriptive names") 
parser.add_argument( #gridfile argument. plate_map.txt
	"-g",
	"--gridfile",
	help = "Path to 96 well plate map containing sample IDs, string",
	required = True)
parser.add_argument( #barcodefile. 515f_806r_illumina_primers_515barcoded.txt
	"-b",
	"--barcodefile",
	help = "Path to barcode file, string",
	required = True)
parser.add_argument( #metadatafile
	"-m",
	"--metadatafile",
	help = "Path to metadata file, string",
	required = True)

args = parser.parse_args()

gridfile = args.gridfile
barcodefile = args.barcodefile
metadatafile = args.metadatafile

header = []
samplemap = {}
wellmap = {}

with open(gridfile) as GRID:
	for i, line in enumerate(GRID): #i is the counter iterating through each line, line is the data stored in each line
		data = line.rstrip() #remove any trailing whitespace
		if i == 0: #line with plate name
			(name, pname) = data.split() #splits by whitespace
			platename = pname #two variables, split by tabs horizontally. take plate name from split. 
		elif i == 1: #line with plate number
			(plate, platenum) = data.split() 
			plateno = platenum.strip() #take plate number from split (second line of file)
		elif i == 2: #column header line
			header = data.split('\t')
		else:			
			row_data = data.split('\t')
			rowname = row_data[0] #grab the first element of the line as the row name
			swabIDs = row_data[1:] #grab the rest of the elements as data
			j = 1
			for ID in swabIDs: #for each sample ID, create a hash entry with the row and column names as the keys
				samplemap[rowname, header[j]] = ID
				well = rowname + header[j]
				wellmap[ID] = well
				j += 1
print "Getting data for plate number: %s" % plateno

#read barcode spreadsheet
##header##	plate	well	name	illumina_5_adapter	golay_barcode	forward_primer	forward_primer_linker	515f_fw_primer	pcr_primer
barcodemap = {}
counter = 0
with open(barcodefile) as BARCODES:
	for i, line in enumerate(BARCODES):
		data = line.rstrip() 
		all_in_line = data.split('\t') 
		plate, well, name, illum5adapter, golaybarcode, fprimer, fprimerlinker, fw515primer, pcrprimer = all_in_line[0:9]
		rest = all_in_line[9:] #we don't need any of this stuff now, but we could modify this script later to do other things with these extra fields if we needed to
		if plate == plateno: #if this line is for the plate we are interested in
			match = re.match(r"([a-z])([0-9]+)", well, re.I) #performs case-insensitive matching of letters and numbers
			if match:
				wellRC = match.groups() #split the matched groups into letters (well row) and numbers (well column) (ex: 'A12' becomes A & 12)
				wellR, wellC = wellRC 
				barcodemap[wellR, wellC] = golaybarcode #dictionary is ready to be accessed by column and row names, same as above loop for grid file
			counter += 1
BARCODES.close()
print "there are %d entries in the barcode spreadsheet for plate %s" % (counter, plateno) 

#third part of the script: read in metadata, add barcodes to it, print out mapping file with barcodes in, new platemap with sample IDs in.
#read metadata file without barcodes. needs the following columns:
##header## sampleID swabID plate project_name person_responsible.
idmapping = {}
header = ""
linkerprimerseq = 'na'
OUTFILE1 = open(platename + ".mapping_file.txt", "w") 
with open(metadatafile) as METADATA:
	for i, line in enumerate(METADATA):
		data = line.strip() 
		if i == 0:
			ID, restofheader = data.split('\t', 1) #split only by first tab; ID is the first element of the header
			OUTFILE1.write("#SampleID \t BarcodeSequence \t LinkerPrimerSequence \t" + restofheader + "\n")
		else:
			sampleID, swabID, plate, project_name, person_responsible, resttoprint = data.split('\t', 5) #split the first five elements as necessary, store the rest in resttoprint
			idmapping[swabID] = sampleID
			well = wellmap[sampleID]
			match = re.match(r"([a-z])([0-9]+)", well, re.I) #performs case-insensitive matching of letters and numbers
			if match:
				wellRC = match.groups() #split the matched groups into letters (well row) and numbers (well column) (ex: 'A12' becomes A & 12)
				wellR, wellC = wellRC 
			barcode = barcodemap[wellR, wellC]
			OUTFILE1.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % 
				(sampleID, barcode, linkerprimerseq, swabID, plate, well, project_name, person_responsible, resttoprint))
		
METADATA.close()

print "finished parsing metadata file"
print "finished printing mapping file"

#now print new barcode plate with descriptive sample names
OUTFILE2 = open(platename + ".platemap.w_sample_names.txt", "w") 
OUTFILE2.write("Plate:\t" + platename + "\nPlate_#:\t" + plateno + "\n") #write in the plate info
colnames = ["1","2","3","4","5","6","7","8","9","10","11","12"]
rownames = ["A","B","C","D","E","F","G","H"]
colnamestoprint = "\t".join(colnames)
OUTFILE2.write("\t" + colnamestoprint + "\n") #print the column names
for row in rownames:
	OUTFILE2.write(row + "\t")
	row_data = []
	for col in colnames:
		try: #filling the list with data
			row_data.append(samplemap[row, col]) 
		except KeyError: #if try to call a well with no information, don't crash just continue
			continue
	toprint = "\t".join(row_data)
	OUTFILE2.write(toprint + "\n")
		
OUTFILE2.close()
print "finished printing new platemap with sampleIDs"

print "done!"
OUTFILE1.close()


