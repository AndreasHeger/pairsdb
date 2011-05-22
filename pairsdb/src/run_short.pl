####
####
##
## Project PairsDB
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: run_short.pl,v 1.1.1.1 2002/07/02 10:46:58 heger Exp $
##
##
####
####


## Run Short for all sequences in fasta file.
# This can be sped up a lot!

$|=1;

use strict;

use lib '/net/cpp-data/backup/andreas/projects/pairsdb/src';

use Short;

my $nid = "";
my $sequence = "";
my $method   = 4;
my $verbose =  0;
my $REPORT_STEP = 1000;
my $MAX_SEQLEN = 10000;

my $iterations = 0;

while(<STDIN>) {
  chop;

  if (/^>/) {
    my ($new_nid)=/^>(\S+)/;

    if ($sequence) {
      $sequence =~ s/\s//g;
      if (length($sequence) > $MAX_SEQLEN) {
	warn "run_short: cowardly refusing to process sequence $nid\n";
      } else {
	  PrintShort($sequence,$nid);
     }
    }

    $nid = $new_nid;
    $sequence = "";

    $iterations ++;
    
    if (($iterations % $REPORT_STEP) == 0) {
	warn "run_short: processing $nid at iteration $iterations\n";
    }

    next;
  }

  $sequence .= $_;
}

PrintShort( $sequence, $nid ) if ($sequence);

############################################################################
############################################################################
# Title: DoShort
# Function: looks for short repeats, i.e. in the range of ktup-size
#---------------------------------------------------------------------------
sub PrintShort {
  my($seq,$hid)=@_;
  my($result)=''; my($irepeat)=1;
  my($flag)=1;

  while($flag) {
      $flag=0;
    my($x, $lrepeat, $score, $diagonal, $mali, $from, $to)=Short::short($seq,$hid); 
    last if($x<1);
    my(@mali)=split(/\n/,$mali); shift(@mali);
    
    my ($min_from, $max_to) = (9999999,0);
    
    foreach (@mali) {
      my($from,$to)=/^\S+\/(\d+)\-(\d+)/;
      next if($to-$from>$lrepeat*2);
      substr($seq,$from-1+$[,$to-$from+1)=~s/./X/g;
      $flag=1;
      $min_from = $from if ($min_from > $from);
      $max_to   = $to   if ($max_to   < $to);
    }
    
    # score is discarded, because info is only five characters long!
    if ( ($max_to >= 1 ) && ($min_from < 9999999) && ($min_from < $max_to) ) {
      print join("\t", $hid, $min_from, $max_to, $method, $lrepeat) . "\n";
    }
  }
}
  









