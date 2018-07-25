"""
A script to concatenate all run stats for each sample from MetaPathways into a single tab-delimited text file named "total_runstats.txt", located in the mp_run_stats/ folder.

input: path to mp_run_stats/ folder under the MetaPathways results folder 

known annoying bug/feature: last line of total_runstats.txt file will have "total_runstats.txt" appended.

#########
# USAGE #
#########
python /path/to/get_metapath_runstats.py -p /path/to/mp_run_stats
test example: python get_metapath_runstats.py -p /Users/parfrey/Desktop/Kevin/weak_cat_lanes_qual_filtered_samples 
"""

import argparse
import os, os.path

print "run: python get_metapath_runstats.py -h for help.\n\nThis script concatenates all the run stats from a MetaPathways results folder into a single file titled 'total_runstats.txt' under mp_run_stats/."
print "\n"

parser = argparse.ArgumentParser(
	description = "One input required: Path to mp_run_stats/ in MetaPathways results folder")
required_arg = parser.add_argument_group("required argument")
required_arg.add_argument(
	"-p",
	"--resultspath",
	help = "Path to MetaPathways results folder, string.",
	required = True)

args = parser.parse_args()

mpres = args.resultspath
one_time = True # used for writing first line in outfile later on

# for each file in the mp_run_stats/ folder. reads through each results file from each sample
for f in os.listdir(mpres):
	if not f.startswith("."): # skip all hidden files in the folder (had an issue with .DS_Store)
		try:
			list_info = []      # list for where categories will be stored
			hits_toprint = {}   # dictionary with keys as sample file names, values as the numbers (results) associated with each category
			hits_from_file = [] # list storing all number (results) associated with each category. values for the above dictionary
			with open(os.path.join(mpres, f)) as fh:
				for i, line in enumerate(fh):
					data = line.strip()
					junk, info, hits = data.split('\t', 3)    # split each input file into 3 variables
					info = info.strip('-')                    # stripping '-' character for better view in excel
					hits_from_file.append(hits)
					list_info.append(info)
			if one_time:
				OUTFILE = open(os.path.join(mpres, "total_runstats.txt"), "w")
				print "printing categories..."
				OUTFILE.write("\t" + "\t".join(list_info) + "\n") # write the first line of outfile. contains all categories/metrics from results files
				one_time = False                                  # one time use. after execution, this conditional will no longer run
			hits_toprint[f] = hits_from_file                          # assign list of values from each file to their respective file names
			print "writing in data for each sample..."
			OUTFILE.write(f + "\t" + "\t".join(hits_toprint[f]) + "\n")	
		except ValueError:
			print "located script in mp_run_stats/ folder, skipping this input..."
			continue
OUTFILE.close()
print "done!"						
