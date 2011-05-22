####
####
##
## Project PairsDB
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: mask_parser.pl,v 1.1.1.1 2002/07/02 10:46:58 heger Exp $
##
##
####
####


#	parse through fasta-output for masked sequences
#	and write out in format
#
#	hid\tfrom\tto\tmethod 
#
#	input: 
#	>hid
#          sequence

use strict;

my ($method) = (@ARGV);

if (!$method) {
  warn "No method supplied!\n";
  exit(0);
}

my $info = "";

#---------------------------------------------------------------------
my $sequence = "";
my $hid = "";
while (<STDIN>) {
  chop;
  if (/^>/) {
    PrintRanges($hid, $sequence) if ($hid);
    ($hid) = /^>(\S+)/;
    $sequence = "";
  } else {
    $sequence .= $_;
  }
}

sub PrintRanges {
  my ($hid, $sequence) = (@_);

  $sequence =~ tr/[A-Z]/[a-z]/;

  my (@sequence) = split (//,$sequence); 
  
  my $i = 1;
  my $from = 1;

  my $in_masked_region;
  if ($sequence[0] eq 'x') {
    $in_masked_region = 1; 
  } else {
    $in_masked_region = 0;
  } 
  my @ranges;
  while (@sequence) {
    $i++;
    my $thischar = shift(@sequence);

    if ($thischar eq 'x') {

      if (!$in_masked_region) {
	$from = $i;
	$in_masked_region = 1;
      } 

    } else {

      if ($in_masked_region) {
	print join ("\t", $hid, $from, $i-1, $method, $info) . "\n";
	$in_masked_region = 0;
      }
    }
  }
}
