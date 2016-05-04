#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;

####a script to create a sample mapping file from 96 well plate metadata files (grid file, barcode spreadsheet -> sample barcode mapping file)####
##IMPORTANT: all inputs should be tab separated text files
##inputs: 
###1. spreadsheet with barcodes in each primer plate. each plate has a unique barcode sequence
###2. visual layout of 96 well plate (grid format). user must make this file.

my $gridfile = $ARGV[0];
chomp $gridfile;
my $barcodefile = $ARGV[1];
chomp $barcodefile;

#read visual layout file

my $plateno; #declare outside of read file loop, check later against full spreadsheet. we only need the data for the plate in question, the rest is useless to us
my $platename; #declare here to use as part of the output filename
my $i=0; #use for discriminating between plate metadata and actual well data lines
my @header;
my %samplemap;
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
	my ($rowname, @sampleIDs) = split(/\t/, $line);
	my $j = 1; #something to iterate over the header array with. counting starts at 1, don't care about element zero
	foreach my $ID (@sampleIDs) { #for each sample ID, create a hash entry with the row and column names as the keys
	    #print STDERR "$rowname, $ID, $header[$j], $j th header item\n";
	    $samplemap{$rowname}{$header[$j]}=$ID;
	    $j++;
	}
    }
$i++;    
}
close GRID;
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

#third part of the script: use the well IDs and the hashes you created to print out all the sample IDs and barcode data in a new spreadsheet
#open and output file that has a useful name based on the plate name. can do this however. this way seems okay for now.
open(OUTFILE, ">$gridfile.${platename}.mapping_file.txt") or die ("cannot open output file!");
print STDERR "printing outfile...\n";
print OUTFILE "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\n"; #print out a header for the output file
foreach my $a (sort keys %barcodemap) { #access both sets of keys and iterate over them
    foreach my $b (sort keys %{$barcodemap{$a}}) {
	#since the two hashes (barcodemap and samplemap) are structured the same way with the same keys, we can just print what we want now:
	print OUTFILE "$samplemap{$a}{$b}\t$barcodemap{$a}{$b}\tNA\t\n"; 
    }
}
#then you're done!
print STDERR "done!\n";
close OUTFILE;

