import re
import argparse
import os.path

__author__ = "Kevin Chan"
__maintainer__ = "Evan Morien"
__contact_email__ = "morien@zoology.ubc.ca"

####a script to create an analysis-ready mapping file and platemap with descriptive names from 96 well plate metadata files (metadata file, plate map, barcode plate spreadsheet -> mapping file w metadata & barcodes, plate map with descriptive names)####
##IMPORTANT: all inputs should be tab separated text files
##inputs:
###1. visual layout of 96 well plate (grid format). user must make this file. (plate map, has swab IDs in place of sample names)
###2. spreadsheet with barcodes in each primer plate. each plate has a unique barcode sequence
###3. metadata file that links swab IDs to sample names, plus other metadata.

##outputs:
###1. A mapping file in plain text format, now with barcodes, primer info, project name, and person responsible, plus other existing metadata (named PROJECTNAME.PLATENAME.mapping_file.txt, where PROJECTNAME is the user-defined name of the project, and PLATENAME is the name of the plate as entered in your input platemap). 
###2. A copy of the platemap (grid format) with sampleIDs instead of swabIDs. (named PROJECTNAME.PLATENAME.platemap.w_sample_names.txt)
##Output files will appear in the same directory as the input platemap.

#########
##USAGE##
#########
#python /path/to/generate.plate_map.add_barcodes_to_metadata.py -g plate_map_with_swab_ids.txt -b barcode_platemap_spreadsheet.txt -m metadata_file_no_barcodes.txt -n project_name -p person_responsible
#test example: python agenerate.plate_map.add_barcodes_to_metadata.py -g 16s_plate_map.txt -b 515f_806r_illumina_primers_515barcoded.txt -m 16s_mapping_file.txt -n my_16s_project -p my_name
#IMPORTANT: see README for rules on file formatting

print "run: python add_barcodes.fix_platemap.py -h for help."
print "\n"

#Description of arguments needed
parser = argparse.ArgumentParser(
	description = "Three inputs required: plate_map.txt barcode_file.txt metadata.txt") 
requiredargs = parser.add_argument_group("required arguments")
requiredargs.add_argument( #gridfile argument. user defined.
	"-g",
	"--gridfile",
	help = "Path to 96 well plate map containing sample IDs, string. See README for this script for formatting instructions.",
	required = True)
requiredargs.add_argument( #barcodefile. 515f_806r_illumina_primers_515barcoded.txt
	"-b",
	"--barcodefile",
	help = "Path to barcode file, string. This is the lab's barcode spreadsheet. Can be found on the lab's github repo. Ask on the lab Slack if you don't know where to find it.",
	required = True)
requiredargs.add_argument( #metadatafile. user defined.
	"-m",
	"--metadatafile",
	help = "Path to metadata file, string. Samples in rows, metadata in columns. You can add whatever metadata you like, but the first few columns must fit what's described in the README file for this script.",
	required = True)
requiredargs.add_argument( #project name argument. user defined.
	"-n",
	"--projectname",
	help = "Name of the project.",
	required = True)
requiredargs.add_argument( #person responsible argument. user defined.
	"-p",
	"--personresponsible",
	help = "Name of the person responsible for this project.",
	required = True)
parser.add_argument( #outputpath. user defined. optional.
	"-o",
	"--outputpath",
	help = "Path to where the output files will be saved, string. If unspecified, default is current working directory.")

args = parser.parse_args()

gridfile = args.gridfile
barcodefile = args.barcodefile
metadatafile = args.metadatafile
projectname = args.projectname
personresponsible = args.personresponsible
outputpath = args.outputpath

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
			j = 1 #we count from 1 and not 0 because the 0th element of list 'header' is empty
			for ID in swabIDs: #for each sample ID
				samplemap[rowname, header[j]] = ID #create a double dictionary storing the sample IDs as info with the plate well and column as keys
				well = rowname + header[j] #join plate well and column for single well ID
				wellmap[ID] = well #create single dictionary linking well ID to sample ID
				j += 1 #increment j

print "Getting data for plate number: %s" % plateno
#read barcode spreadsheet
##header##	plate	well	name	golay_barcode	forward_primer	reverse_primer	forward_primer_name	reverse_primer_name	rest...
barcodemap = {}
counter = 0
fprimer = ""
rpirmer = ""
fprimername = ""
rprimername = ""
with open(barcodefile) as BARCODES: #open barcode file
	for i, line in enumerate(BARCODES): #for each line
		data = line.rstrip() #strip whitespace
		all_in_line = data.split('\t') #split on tab
		plate, well, name, golaybarcode, fprimer, rprimer, fprimername, rprimername = all_in_line[0:8] 
		rest = all_in_line[9:] #we don't need any of this stuff now, but we could modify this script later to do other things with these extra fields if we needed to
		if plate == plateno: #if this line is for the plate we are interested in
			match = re.match(r"([a-z])([0-9]+)", well, re.I) #performs case-insensitive matching of letters and numbers to collect the well IDs
			if match: #if we have a well ID
				wellRC = match.groups() #split the matched groups into letters (well row) and numbers (well column) (ex: 'A12' becomes A & 12)
				wellR, wellC = wellRC #list wellRC becomes strings wellR and wellC
				#print "here is our barcode well: %s%s" % (wellR.lower(),wellC)
				barcodemap[wellR.lower(), wellC] = golaybarcode #dictionary is ready to be accessed by column and row names, same structure as above loop for grid file #ensure that we're using lower case letters for the barcode wells
			counter += 1 #increment a counter so we know how many entries we've got in the barcode file that match our plate. should be 96 for each plate.
BARCODES.close()
print "there are %d entries in the barcode spreadsheet for plate %s" % (counter, plateno)
print "if there are fewer than 96 entries for the above plate then there is a problem!\n"

#third part of the script: read in metadata, add barcodes to it, print out mapping file with barcodes in, new platemap with sample IDs in.
#read metadata file without barcodes. needs the following columns:
##header## sampleID swabID plate

idmapping = {}
header = ""
linkerprimerseq = 'na'
try:
	completeName = os.path.join(outputpath, projectname + "." + platename + ".mapping_file.txt")
	OUTFILE1 = open(completeName, "w") #create an outfile in the specified directory.
except AttributeError:
	OUTFILE1 = open(projectname + "." + platename + ".mapping_file.txt", "w") #create an outfile, will be created in the working directory if unspecified.
with open(metadatafile) as METADATA: #read in metadata file
	for i, line in enumerate(METADATA): #for each line
		data = line.strip() #strip whitespace
		if i == 0: #if it's the first line
			ID, restofheader = data.split('\t', 1) #split only by first tab; ID is the first element of the header
			OUTFILE1.write("#SampleID \t BarcodeSequence \t LinkerPrimerSequence \t project_name \t person_responsible \t forward_primer \t reverse_primer \t forward_primer_name \t reverse_primer_name \t" + restofheader + "\n") #write the header for the output file
		else:
			sampleID, swabID, plate, barcode_well, resttoprint = data.split('\t', 4) #split the first five elements as necessary, store the rest in resttoprint
			idmapping[swabID] = sampleID #link swab IDs to sample IDs with a dictionary
			try: #try to retrieve well ID from dictionary and write to the outfile
				well = wellmap[swabID] #retrieve well ID from wellmap dictionary
				match = re.match(r"([a-z])([0-9]+)", well, re.I) #performs case-insensitive matching of letters and numbers
				if match:
					wellRC = match.groups() #split the matched groups into letters (well row) and numbers (well column) (ex: 'A12' becomes A & 12)
					wellR, wellC = wellRC #list wellRC becomes strings wellR and wellC
				barcode = barcodemap[wellR.lower(), wellC] #retrieve barcode from dictionary 'barcodemap' using the well column and row
				#print "trying to match well %s on plate %s with swabID %s and sampleID %s, here is the match: %s and here is the barcode: %s" % (barcode_well, plate, swabID, sampleID, wellRC, barcode)
				OUTFILE1.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % 
					(sampleID, barcode, linkerprimerseq, projectname, personresponsible, fprimer, rprimer, fprimername, rprimername, swabID, plate, well, resttoprint)) #write all collected info to the output file
			except KeyError: #if we encounter a sample ID in the metadata file but not in the platemap, continue running and do not crash.
				print "sample %s not on platemap" % sampleID
				continue
METADATA.close()
print "finished parsing metadata file"
OUTFILE1.close()
print "finished printing mapping file"

#now print new barcode plate with descriptive sample names
try: #this works if the user defines an output path when running the script
	completeName = os.path.join(outputpath, projectname + "." + platename + ".platemap.w_sample_names.txt")
	OUTFILE2 = open(completeName, "w") #create an outfile in the specified directory.
except AttributeError:
	OUTFILE2 = open(projectname + "." + platename + ".platemap.w_sample_names.txt", "w") #create an outfile, will be created in the working directory.
OUTFILE2.write("Plate:\t" + platename + "\nPlate_#:\t" + plateno + "\n") #write in the plate info on two header lines
colnames = ["1","2","3","4","5","6","7","8","9","10","11","12"]
rownames = ["a","b","c","d","e","f","g","h"]
colnamestoprint = "\t".join(colnames)
OUTFILE2.write("\t" + colnamestoprint + "\n") #print the column names on the next line
for row in rownames: #for each row in the file
	OUTFILE2.write(row + "\t") #print the row name
	row_data = []
	for col in colnames: #for each column
		try: #filling the list with data
			row_data.append(idmapping[samplemap[row, col]]) #append values from dictionary 'idmapping' (using keys row and column of dictionary 'samplemap') to the list 'row_data'
		except KeyError: #if we try to call a well with no information, don't crash just continue. it just means there is empty space in the platemap.
			continue
	toprint = "\t".join(row_data) #create tab separated string from collected row data
	OUTFILE2.write(toprint + "\n") #and write it to the output file
		
OUTFILE2.close()
print "finished printing new platemap with sampleIDs"
print "done!"
