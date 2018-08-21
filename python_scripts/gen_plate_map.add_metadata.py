"""
script to generate plate map(s) and metadata file with barcodes, barcode plate number, barcode well number added 
metadata file -> plate map(s), metadata w/ barcodes, barcode plate #, barcode well # added (QIIME readable)
IMPORTANT: input files must be tab-delimited

inputs:
1. metadata file with sample IDs, swab IDs, project name, person responsible, and other metadata (in that order)
2. spreadsheet with barcodes in each primer plate. each plate has a unique barcode sequence
3. (optional) starting plate number (1-10)
4. (optional) starting well (A1-H12)
5. (optional) assigning barcodes by row rather than column
6. (optional) path to output folder

outputs:
1. a mapping file in plain text format, now with plate, barcode well, and barcodes columns added and populated with values (named METADATA.barcodes_added.txt, where METADATA is the inputted metadata filename)
2. a platemap with sample IDs corresponding to the generated metadata file. if more than one plate is required, additional platemaps will be generated (named project_name_sample_plate_num_X.txt, where X is the plate number)
3. a platemap with swab IDs corresponding to the generated metadata file. if more than one plate is required, additional platemaps will be generated (named project_name_swab_plate_num_X.txt, where X is the plate number)

#########
# USAGE #
#########
# example usage: python /path/to/gen_plate_map.add_metadata.py -m /path/to/metadatafile.txt -b /path/to/barcodefile.txt -p 2,7 -w D4 --row -o out/
IMPORTANT: see README for rules on file formatting
"""

import argparse
import os
import re
from sys import exit
from collections import defaultdict

__author__ = "Kevin Chan"
__email__ = "kevchan1@alumni.ubc.ca"

print "run: plate_map.add_metadata.py -h for help."
print "\n"

parser = argparse.ArgumentParser(
	description = "Two inputs required: metadata.txt barcode_file.txt")
requiredargs = parser.add_argument_group("required arguments")
requiredargs.add_argument( # metadatafile. user defined.
	"-m",
	"--metadatafile",
	help = "Path to metadata file for barcodes, plates and wells to be added to, string.",
	required = True)
requiredargs.add_argument( # barcodefile. 515f_806r_illumina_primers_515barcoded.txt
	"-b",
	"--barcodefile",
	help = "Path to barcode file, string. This is the lab's barcode spreadsheet. Can be found on in the botany shared disk space. Ask on the lab Slack if you don't know where to find it.",
	required = True)
parser.add_argument(       # outputpath. user defined. optional.
	"-o",
	"--outputpath",
	default = "./",
	help = "Path to where the output files will be saved, string. If unspecified, default is current working directory.")
parser.add_argument( 	   # platenumber. user defined. optional.
	"-p",
	"--platenumber",
	metavar = "[0-10]",
	help = "Plate number(s) for the generated platemap(s). If multiple plates, numbers must be separated by commas. If unspecified, default is 1. Must be between 0-10.")
parser.add_argument(	   # well number. user defined. optional.
	"-w",
	"--well",
	metavar = "[A1-H12]",
	help = "Well number for the generated platemap, string. If unspecified, default is 'A1'. Must be between A1-H12.")
parser.add_argument(	   # specify to assign barcodes by row. user defined. optional.
	"--row",
	action = 'store_true',
	help = "If specified, barcodes will be assigned by rows (ex. A1, A2, etc.). Default assigning is by columns.")

args = parser.parse_args()

metadatafile = args.metadatafile
barcodefile = args.barcodefile
out = args.outputpath
platenum = args.platenumber
start_well = args.well
row_specification = args.row

# update starting plate number and well (if user specified)
# inputs:	platenum (string)	 the plate number to start on
#			start_well (string)  the well to start on
# outputs:  the updated plate and well numbers
def updater(platenum, start_well):
	if platenum is None:					   # if platenum was not specified, set it to 1
		platenum = str(1)					   # set as string for consistency
	if start_well is None:					   # if start_well was not specified, set it to A1
		start_well = 'A1'
	if ',' in platenum:						   # indicates multiple plates delimited by commas as input
		platenum = platenum.split(',')		   # split by commas into a list of plates
	return platenum, start_well		   		   # return platenum, well in a tuple

if start_well == 'D6':
	print "Warning: well D6 was selected, D7 will be used instead."
elif start_well == 'H12':
	print "Warning: well H12 was selected, A1 on plate %d will be used instead." % (int(updater(platenum, start_well)[0]) + 1)

colnames = map(str, range(1, 13)) 			   # list of column names (exclusive of second item, hence 13 instead of 12). format as string
rownames = ["A","B","C","D","E","F","G","H"]   # list of row names                    

# update the well position in the list of wells from A1-H12
# input:	well (string)	the initial well specified by user. same input from updater() function
#			row  (boolean)  true if assigning barcodes by rows
# outputs: 	an updated well position in the list of wells. position range 0-95.
#			a list of all wells in a 96 well plate (A1-H12)
def update_well_pos(well, row):
	barcode_wells = []		# if assigning barcodes by column
	barcode_wells_row = []  # if assigning barcodes by row
	# list to store wells by row
	if row:					
		for row in rownames:							   # opposite order for rows
			for col in colnames:
				barcode_well_row = row + col           	   # create a barcode well 
				barcode_wells_row.append(barcode_well_row) # append each barcode well to a list. the list consists of wells from A1-H12. order: A1, A2, ..., C4, C5, ..., H11, H12				
		for well_pos in [well_pos for well_pos, well_content in enumerate(barcode_wells_row) if well_content == well]: # fancy list comprehension
			return well_pos, barcode_wells_row	
	# list to store wells by column
	for col in colnames:
		for row in rownames:
			barcode_well = row + col           # create a barcode well 
			barcode_wells.append(barcode_well) # append each barcode well to a list. the list consists of wells from A1-H12. order: A1, B1, C1, ..., G12, H12
	# get the position for the specified well in barcode_wells. barcode_wells has 96 elements (0-95)
	for well_pos in [well_pos for well_pos, well_content in enumerate(barcode_wells) if well_content == well]: # fancy list comprehension
		return well_pos, barcode_wells
print "creating wells"

well_pos = update_well_pos(updater(platenum, start_well)[1], row_specification)[0] 	    # the first value from the function update_well_pos is the specified well (default A1)
barcode_wells = update_well_pos(updater(platenum, start_well)[1], row_specification)[1] # the second value from the function is the barcode_wells or barcode_wells_row list (A1-H12)

# get the filename of the metadatafile
metadata_filename = os.path.basename(metadatafile)
metadata_filename = os.path.splitext(metadata_filename)[0]

# read barcode spreadsheet
##header##	plate	well	name	illumina_5_adapter	golay_barcode	forward_primer	forward_primer_linker	515f_fw_primer	pcr_primer
barcodemap = {}
with open(barcodefile, "U") as BARCODES: # open barcode file
	for i, line in enumerate(BARCODES): # for each line
		data = line.rstrip() # strip whitespace
		all_in_line = data.split('\t') # split on tab
		plate, well, name, illum5adapter, golaybarcode, fprimer, fprimerlinker, fw515primer, pcrprimer = all_in_line[0:9] # first 9 columns are useful or may be useful in the future
		rest = all_in_line[9:] # we don't need any of this stuff now, but we could modify this script later to do other things with these extra fields if we needed to
		match = re.match(r"([a-z])([0-9]+)", well, re.I) # performs case-insensitive matching of letters and numbers to collect the well IDs
		if match: # if we have a well ID
			wellRC = match.groups() # split the matched groups into letters (well row) and numbers (well column) (ex: 'A12' becomes A & 12)
			wellR, wellC = wellRC # list wellRC becomes strings wellR and wellC. required for concatenation in next line
			barcodemap[wellR + wellC + "_" + plate] = golaybarcode # create a dictionary with keys composed of the well and plate number, separated by an underscore ("_")
print "finished parsing barcode file"

linkerprimerseq = 'na'
well_pos = update_well_pos(updater(platenum, start_well)[1], row_specification)[0]
plate_num = updater(platenum, start_well)[0]  # counter to start new plates after each full iteration of barcode_wells list
if type(plate_num) is str: # change str to int to iterate through plates
	plate_num = int(plate_num)
pcount = 0     # count for comma-delimited list of plates
samplemap = {} # to store sample IDs 
swabmap = {}   # to store swab IDs
proj_map = defaultdict(list) # to store unique project names with plate numbers as keys. stores values as a list
proj_title = [] # create a new list to store unique project names in metadata file
oneTime = True
# start by opening metadata file w/ sample IDs. populate the metadata file w/ barcodes, barcode plate #, barcode well #
completeName = os.path.join(out, metadata_filename + ".barcodes_added.txt")
MDout = open(completeName, "w")      						   # create an outfile in the specified output directory
with open(metadatafile, "U") as METADATA:
	for i, line in enumerate(METADATA):
		data = line.strip() # remove all whitespace
		all_in_line = data.split("\t")
		try:
			sampleID, swabID, project_name, person_responsible = all_in_line[0:4] # first four columns for each row
		except ValueError: # exception for blank last line in excel generated tab-delimited files. continue running
			continue
		resttoprint = all_in_line[4:] # all other metadata. will be an empty list if no other metadata present
		if i == 0: # get all headers from metadata file
			restofheader = all_in_line[2:] # rest of header column to print starting with project_name and person_responsible
			MDout.write("#SampleID \t BarcodeSequence \t LinkerPrimerSequence \t swabID \t plate \t barcode_well \t" + "\t".join(restofheader) + "\n")
		else:
			if oneTime: # append the first project name only once. the rest of the project names will be appended when a new plate starts
				if type(plate_num) is int:
					proj_map[plate_num].append(project_name)
				if type(plate_num) is list:
					proj_map[plate_num[pcount]].append(project_name)
				proj_title.append(project_name)
				oneTime = False
			if type(plate_num) is int and not project_name in proj_map[plate_num]: # if the project name isn't already in the associated plate
				proj_map[plate_num].append(project_name)
			if type(plate_num) is list and not project_name in proj_map[plate_num[pcount]]: # same for if plate_num is as list
				proj_map[plate_num[pcount]].append(project_name)
			blank_md = ["unknown" for md in resttoprint] # metadata columns for the blank rows. will be filled in with 'unknown'
			if barcode_wells[well_pos] == 'D6': # if D6 is the current well
				# write the blank1 line for D6
				if isinstance(plate_num, int): # if the plate number is an int
					MDout.write("blank1" +"\t"+ barcodemap[barcode_wells[well_pos] + "_" + str(plate_num)] +"\t"+ linkerprimerseq +"\t"+ 'none' +"\t"+ str(plate_num) +"\t"+ barcode_wells[well_pos] +"\t"+ project_name +"\t"+ person_responsible +"\t"+ "\t".join(blank_md) + "\n")
					# write the line for E6 continuing with the appropriate sample ID
					MDout.write(sampleID +"\t"+ barcodemap[barcode_wells[well_pos + 1] + "_" + str(plate_num)] +"\t"+ linkerprimerseq +"\t"+ swabID +"\t"+ str(plate_num) +"\t"+ barcode_wells[well_pos + 1] +"\t"+ project_name +"\t"+ person_responsible +"\t"+ "\t".join(resttoprint) + "\n")
					samplemap['D6_' + str(plate_num)] = 'blank1' 
					swabmap['D6_' + str(plate_num)] = 'blank1' 					
					if row_specification: # if barcodes are assigned by rows
						samplemap['D7_' + str(plate_num)] = sampleID # append samples/swabs to D7 (same row)
						swabmap['D7_' + str(plate_num)] = swabID
					else:
						samplemap['E6_' + str(plate_num)] = sampleID # append samples/swabs to E6 (same column)
						swabmap['E6_' + str(plate_num)] = swabID
				else: # case where plate number is a list
					MDout.write("blank1" +"\t"+ barcodemap[barcode_wells[well_pos] + "_" + plate_num[pcount]] +"\t"+ linkerprimerseq +"\t"+ 'none' +"\t"+ plate_num[pcount] +"\t"+ barcode_wells[well_pos] +"\t"+ project_name +"\t"+ person_responsible +"\t"+ "\t".join(blank_md) + "\n")
					# write the line for E6 continuing with the appropriate sample ID
					MDout.write(sampleID +"\t"+ barcodemap[barcode_wells[well_pos + 1] + "_" + plate_num[pcount]] +"\t"+ linkerprimerseq +"\t"+ swabID +"\t"+ plate_num[pcount] +"\t"+ barcode_wells[well_pos + 1] +"\t"+ project_name +"\t"+ person_responsible +"\t"+ "\t".join(resttoprint) + "\n")
					samplemap['D6_' + plate_num[pcount]] = 'blank1' 
					swabmap['D6_' + plate_num[pcount]] = 'blank1' 
					if row_specification: # if barcodes are assigned by rows
						samplemap['D7_' + plate_num[pcount]] = sampleID # append samples/swabs to D7 (same row)
						swabmap['D7_' + plate_num[pcount]] = swabID
					else:
						samplemap['E6_' + plate_num[pcount]] = sampleID # append samples/swabs to E6 (same column)
						swabmap['E6_' + plate_num[pcount]] = swabID
				well_pos += 2 # increment by 2 to skip D6 (next line starts at E6 or D8)
			elif well_pos == (len(barcode_wells) - 1): # once the end of barcode_wells is reached (H12)
				if isinstance(plate_num, int):
					# write the blank2 line for well H12
					MDout.write("blank2" +"\t"+ barcodemap[barcode_wells[well_pos] + "_" + str(plate_num)] +"\t"+ linkerprimerseq +"\t"+ 'none' +"\t"+ str(plate_num) +"\t"+ barcode_wells[well_pos] +"\t"+ project_name +"\t"+ person_responsible +"\t"+ "\t".join(blank_md) + "\n")
					samplemap['H12_' + str(plate_num)] = 'blank2'
					swabmap['H12_' + str(plate_num)] = 'blank2'		
					well_pos = 0   	   # reset the counter
					plate_num += 1 	   # starting a new plate
					if plate_num > 10: # if the plate number exceeds the maximum alloted
						exit("The number of plates requested exceeded the maximum (10). Please specify a different starting plate number.")
					MDout.write(sampleID +"\t"+ barcodemap[barcode_wells[well_pos] + "_" + str(plate_num)] +"\t"+ linkerprimerseq +"\t"+ swabID +"\t"+ str(plate_num) +"\t"+ barcode_wells[well_pos] +"\t"+ project_name +"\t"+ person_responsible +"\t"+ "\t".join(resttoprint) + "\n")
					proj_title.append(project_name)
					samplemap[barcode_wells[well_pos] + "_" + str(plate_num)] = sampleID 
					swabmap[barcode_wells[well_pos] + "_" + str(plate_num)] = swabID
				else: # if the plate number is not an int, it will be list
					# write the blank2 line for well H12
					MDout.write("blank2" +"\t"+ barcodemap[barcode_wells[well_pos] + "_" + plate_num[pcount]] +"\t"+ linkerprimerseq +"\t"+ 'none' +"\t"+ plate_num[pcount] +"\t"+ barcode_wells[well_pos] +"\t"+ project_name +"\t"+ person_responsible +"\t"+ "\t".join(blank_md) + "\n")
					proj_title.append(project_name)
					samplemap['H12_' + plate_num[pcount]] = 'blank2'
					swabmap['H12_' + plate_num[pcount]] = 'blank2'		
					well_pos = 0   	   # reset the counter
					pcount += 1 	   # starting a new plate
					try: # try to write the next line by calling plate_num with the index pcount
						MDout.write(sampleID +"\t"+ barcodemap[barcode_wells[well_pos] + "_" + plate_num[pcount]] +"\t"+ linkerprimerseq +"\t"+ swabID +"\t"+ plate_num[pcount] +"\t"+ barcode_wells[well_pos] +"\t"+ project_name +"\t"+ person_responsible +"\t"+ "\t".join(resttoprint) + "\n")
					except IndexError: # occurs if user does not specify enough plates 
						plate_num.append(str(int(plate_num[pcount - 1]) + 1)) # append a new plate (one greater than the last supplied plate)
						print "Warning: not enough plates were entered, additional plate #%s added" % plate_num[-1]
						MDout.write(sampleID +"\t"+ barcodemap[barcode_wells[well_pos] + "_" + plate_num[pcount]] +"\t"+ linkerprimerseq +"\t"+ swabID +"\t"+ plate_num[pcount] +"\t"+ barcode_wells[well_pos] +"\t"+ project_name +"\t"+ person_responsible +"\t"+ "\t".join(resttoprint) + "\n")
						continue
					samplemap[barcode_wells[well_pos] + "_" + plate_num[pcount]] = sampleID 
					swabmap[barcode_wells[well_pos] + "_" + plate_num[pcount]] = swabID
				well_pos += 1
			else:
				try: # conditional depending on the type of plate_num. same idea, different syntax
					if isinstance(plate_num, int):
						barcode = barcodemap[barcode_wells[well_pos] + "_" + str(plate_num)] # get the appropriate barcode specific to each well and plate number
					else:
						barcode = barcodemap[barcode_wells[well_pos] + "_" + plate_num[pcount]] # get the appropriate barcode specific to each well and plate number
				except KeyError:
					continue
				if isinstance(plate_num, int):
					MDout.write(sampleID +"\t"+ barcode +"\t"+ linkerprimerseq +"\t"+ swabID +"\t"+ str(plate_num) +"\t"+ barcode_wells[well_pos] +"\t"+ project_name +"\t"+ person_responsible +"\t"+ "\t".join(resttoprint) + "\n")
					samplemap[barcode_wells[well_pos] + "_" + str(plate_num)] = sampleID # create a dictionary holding sample IDs with keys as wells and plate number separated by an underscore ("_")
					swabmap[barcode_wells[well_pos] + "_" + str(plate_num)] = swabID     # create a dictionary holding swab IDs with keys as wells and plate number separated by an underscore ("_")
				else:
					MDout.write(sampleID +"\t"+ barcode +"\t"+ linkerprimerseq +"\t"+ swabID +"\t"+ plate_num[pcount] +"\t"+ barcode_wells[well_pos] +"\t"+ project_name +"\t"+ person_responsible +"\t"+ "\t".join(resttoprint) + "\n")
					samplemap[barcode_wells[well_pos] + "_" + plate_num[pcount]] = sampleID # create a dictionary holding sample IDs with keys as wells and plate number separated by an underscore ("_")
					swabmap[barcode_wells[well_pos] + "_" + plate_num[pcount]] = swabID     # create a dictionary holding swab IDs with keys as wells and plate number separated by an underscore ("_")					
				well_pos += 1 # increment well position by 1
print "finished parsing metadata file"
print "finished printing new metadata file w/ barcodes, barcode plate #, and barcode well # added"

# print new barcoded plates (sample IDs and swab IDs). each plate will have its own tab-delimited file. each plate will be printed out twice - once for sample IDs, once for swab IDs
# each line is duplicated for sample and swab ID generation
# conditional has same idea, different syntax depending on type of plate_num
proj_iterator = 0
if type(plate_num) is int:
	for num in range(int(updater(platenum, start_well)[0]), plate_num + 1): # for each plate starting on the plate number specified (default 1)
		completeName = os.path.join(out, "%s_sample_plate_num_%d.txt" % (proj_title[proj_iterator], num))
		SAMPLEPLATEout = open(completeName, "w")   # create a platemap (to be filled in w/ sample IDs)
		completeName = os.path.join(out, "%s_swab_plate_num_%d.txt" % (proj_title[proj_iterator], num))
		SWABPLATEout = open(completeName, "w")     # create a platemap (to be filled in w/ swab IDs)
		
		SAMPLEPLATEout.write("Project_name:\t" + ",".join(proj_map[num]) + "\nPlate_#:\t" + str(num) + "\n")  # formatting first line to include plate name + plate number
		SAMPLEPLATEout.write("\t"+ "\t".join(colnames) + "\n")							 # formatting second line to include all column numbers 1-12
		SWABPLATEout.write("Project_name:\t" + ",".join(proj_map[num]) + "\nPlate_#:\t" + str(num) + "\n")
		SWABPLATEout.write("\t"+ "\t".join(colnames) + "\n")
		for row in rownames:
			SAMPLEPLATEout.write(row + "\t") # print each row name
			SWABPLATEout.write(row + "\t")
			row_data_sample = []
			row_data_swab = []
			for col in colnames:
				try: # try to find a sample ID with the current well
					row_data_sample.append(samplemap[row + col + "_" + str(num)]) # append sample IDs corresponding to the barcode well and plate number in a list
					row_data_swab.append(swabmap[row + col + "_" + str(num)])	  # append swab IDs corresponding to barcode well and plate number in a list
				except KeyError: 			   # if no sample ID is found
					row_data_sample.append('') # append an empty string
					row_data_swab.append('')	  				
					continue				   # continue running
			toprint_sample = "\t".join(row_data_sample) # create tab separated rows based on list data
			toprint_swab = "\t".join(row_data_swab)
			SAMPLEPLATEout.write(toprint_sample + "\n") # writing to output files
			SWABPLATEout.write(toprint_swab + "\n")
		proj_iterator += 1
else:
	for num in plate_num: # for each plate starting on the plate number specified (default 1)
		completeName = os.path.join(out, "%s_sample_plate_num_%s.txt" % (proj_title[proj_iterator], num))
		SAMPLEPLATEout = open(completeName, "w")   # create a platemap (to be filled in w/ sample IDs)
		completeName = os.path.join(out, "%s_swab_plate_num_%s.txt" % (proj_title[proj_iterator], num))
		SWABPLATEout = open(completeName, "w")     # create a platemap (to be filled in w/ swab IDs)

		SAMPLEPLATEout.write("Project_name:\t" + ",".join(proj_map[plate_num[proj_iterator]]) + "\nPlate_#:\t" + num + "\n")  # formatting first line to include plate name + plate number
		SAMPLEPLATEout.write("\t"+ "\t".join(colnames) + "\n")							 # formatting second line to include all column numbers 1-12
		SWABPLATEout.write("Project_name:\t" + ",".join(proj_map[plate_num[proj_iterator]]) + "\nPlate_#:\t" + num + "\n")
		SWABPLATEout.write("\t"+ "\t".join(colnames) + "\n")
		for row in rownames:
			SAMPLEPLATEout.write(row + "\t") # print each row name
			SWABPLATEout.write(row + "\t")
			row_data_sample = []
			row_data_swab = []
			for col in colnames:
				try: # try to find a sample ID with the current well
					row_data_sample.append(samplemap[row + col + "_" + num]) # append sample IDs corresponding to the barcode well and plate number in a list
					row_data_swab.append(swabmap[row + col + "_" + num])	  # append swab IDs corresponding to barcode well and plate number in a list
				except KeyError: 			   # if no sample ID is found
					row_data_sample.append('') # append an empty string
					row_data_swab.append('')	  				
					continue				   # continue running
			toprint_sample = "\t".join(row_data_sample) # create tab separated rows based on list data
			toprint_swab = "\t".join(row_data_swab)
			SAMPLEPLATEout.write(toprint_sample + "\n") # writing to output files
			SWABPLATEout.write(toprint_swab + "\n")
		proj_iterator += 1
print "finished printing new platemap with sample IDs"
print "finished printing new platemap with swab IDs too!"
print "done!"
