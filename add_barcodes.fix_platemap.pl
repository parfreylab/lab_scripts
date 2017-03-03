#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;

####a script to create an analysis-ready mapping file and platemap with descriptive names from 96 well plate metadata files (metadata file, plate map, barcode plate spreadsheet -> mapping file w metadata & barcodes, plate map with descriptive names)####
##IMPORTANT: all inputs should be tab separated text files
##inputs:
###1. visual layout of 96 well plate (grid format). user must make this file. (plate map, has swab IDs in place of sample names)
###2. spreadsheet with barcodes in each primer plate. each plate has a unique barcode sequence
###3. metadata file that links swab IDs to sample names, plus other metadata.

#########
##USAGE##
#########
#perl scripts/add_barcocdes.fix_platemap.pl plate_map_with_swab_ids.txt barcode_platemap_spreadsheet.txt metadata_file_no_barcodes.txt
#test example: perl scripts/add_barcocdes.fix_platemap.pl test/SEA3platemap.txt 515f_806r_illumina_primers_515barcoded.txt test/SEA3.Seagrass_Dec2016_map.no_barcodes.txt
#IMPORTANT: see README for rules on file formatting 

my $gridfile = $ARGV[0];
chomp $gridfile;
my $barcodefile = $ARGV[1];
chomp $barcodefile;
my $metadatafile = $ARGV[2];
chomp $metadatafile;

#read visual layout file

my $plateno; #declare outside of read file loop, check later against full spreadsheet. we only need the data for the plate in question, the rest is useless to us
my $platename; #declare here to use as part of the output filename
my $i=0; #use for discriminating between plate metadata and actual well data lines
my @header;
my %samplemap;
my %wellmap;
open (GRID, "<$gridfile") or die ("cannot open $gridfile\n");
while(my $line = <GRID>) {
    chomp $line; #remove trailing whitespace if any
    if ($i == 0) { #line with plate name
	my ($foo, $bar) = split(/\t/, $line);
	$platename=$bar;
    }
    elsif ($i==1) { #line with plate number
	my ($foo, $bar) = split(/\t/, $line);
	$plateno = $bar;
    }
    elsif ($i==2){ #column header line
	@header = split(/\t/, $line);
    }
    else {
	my ($rowname, @swabIDs) = split(/\t/, $line);
	my $j = 1; #something to iterate over the header array with. counting starts at 1, don't care about element zero
	foreach my $ID (@swabIDs) { #for each sample ID, create a hash entry with the row and column names as the keys
	    #print STDERR "$rowname, $ID, $header[$j], $j th header item\n";
	    $samplemap{$rowname}{$header[$j]}=$ID;
	    my $well = $rowname . $header[$j];
	    $wellmap{$ID}=$well;
	    $j++;
	}
    }
$i++;    
}
print STDERR "Getting data for plate number: $plateno\n";

#read barcode spreadsheet
##header##	plate	well	name	illumina_5_adapter	golay_barcode	forward_primer	forward_primer_linker	515f_fw_primer	pcr_primer
my $counter = 0;
my %barcodemap;
open (BARCODES, "<$barcodefile") or die ("cannot open $barcodefile\n");
while (my $line = <BARCODES>) {
    chomp $line;
    my ($plate, $well, $name, $illum5adapter, $golaybarcode, $fprimer, $fprimerlinker, $fw515primer, $pcrprimer, @rest) = split(/\t/, $line); #we don't need any of this stuff now, but we could modify this script later to do other things with these extra fields if we needed to
    if ($plate eq $plateno) { #if this line is for the plate we are interested in
	my ($wellR, $wellC) = split("", $well, 2); #split on characters (example: A1 and B6 become A & 1, B & 6, respectively). don't split into more than two elements, so we don't accidentally end up with 3 strings for A12, H10, etc.
	$barcodemap{$wellR}{$wellC} = $golaybarcode; #hash is ready to be accessed by column and row names, same as above loop for grid file
	$counter++;
    }
}
close BARCODES;
print STDERR "there are $counter entries in the barcode spreadsheet for plate $plateno\n";

#third part of the script: read in metadata, add barcodes to it, print out mapping file with barcodes in, new platemap with sample IDs in.
#read metadata file without barcodes. needs the following columns:
##header## sampleID swabID plate project_name person_responsible.
my %idmapping;
my $counter2=0;
my $header;
my $linkerprimerseq="na";
open(OUTFILE1, ">${platename}.mapping_file.txt") or die ("cannot open output file!");
open (METADATA, "<$metadatafile") or die ("cannot open $metadatafile\n");
while (my $line = <METADATA>) {
    chomp $line;
    if ($counter2 == 0) {
	my @header = split(/\t/, $line);
	my $id = shift @header;
	my $restofheader = join "\t", @header;
	#print STDERR $id, "\t", "barcode", "\t", $restofheader, "\n"; #sanity check
	print OUTFILE1 "#SampleID", "\t", "BarcodeSequence", "\t", "LinkerPrimerSequence", "\t", $restofheader, "\n";
    }
    else {
	my ($sampleID, $swabID, $plate, $project_name, $person_responsible, @rest) = split(/\t/, $line); #we don't need any of this stuff now, but we could modify this script later to do other things with these extra fields if we needed to
	my $resttoprint = join "\t", @rest;
	$idmapping{$swabID}=$sampleID;
	my $well = $wellmap{$sampleID};
	my ($wellR, $wellC) = split("", $well, 2); #split on characters (example: A1 and B6 become A & 1, B & 6, respectively). don't split into more than two elements, so we don't accidentally end up with 3 strings for A12, H10, etc.
	#create a new plate map hash with the sample IDs instead of swab IDs
	my $barcode = $barcodemap{$wellR}{$wellC}; #hash is ready to be accessed by column and row names, same as above loop for grid file
	print OUTFILE1 "$sampleID\t$barcode\t$linkerprimerseq\t$swabID\t$plate\t$well\t$project_name\t$person_responsible\t$resttoprint\n";	
    }
    $counter2++;
}
close METADATA;
print STDERR "finished parsing metadata file\n";
print STDERR "finished printing mapping file\n";

#now print new barcode plate with descriptive sample names
open(OUTFILE2, ">${platename}.platemap.w_sample_names.txt") or die ("cannot open output file!");
print OUTFILE2 "Plate:\t$platename\nPlate_#:\t$plateno\n";#print the plate info
my @colnames = ("1","2","3","4","5","6","7","8","9","10","11","12");
my @rownames = ("A","B","C","D","E","F","G","H");
my $colnamestoprint = join "\t", @colnames;
print OUTFILE2 "\t$colnamestoprint\n"; #print the column names
foreach my $row (@rownames) {
    print OUTFILE2 "$row\t";
    my @row_data;
    foreach my $col (@colnames) {
	push(@row_data, $samplemap{$row}{$col});
    }
    my $toprint = join "\t", @row_data;
    print OUTFILE2 "$toprint\n";
}
close OUTFILE2;
print STDERR "finished printing new platemap with sampleIDs\n";

print STDERR "done!\n";
close OUTFILE1;

