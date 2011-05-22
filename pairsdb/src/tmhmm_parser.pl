####
####
##
## Project PairsDB
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: tmhmm_parser.pl,v 1.1.1.1 2002/07/02 10:46:58 heger Exp $
##
##
####
####


# parse tmhmm output  
# into format
# hid@from@to@method

use strict;

my ($method) = (@ARGV);
my $type = "";

if (!$method) {
  warn "No method supplied!\n";
  exit(0);
}

my $hid;
while (<STDIN>) {
  if (/^>/) {
    ($hid) = /^>(\S+)/;
    next;
  }
  if (/^%pred/) {
    my ($pred) = /^\%pred NB\(0\): (.*)/;
    my @data = split (/,/, $pred);
    foreach (@data) {
      s/^\s+//;
      my ($mode, $from, $to) = split (/\s+/);
      if ($mode eq 'M') {
	print join ("\t", $hid, $from, $to, $method, $type) . "\n";
      }
    }
    next;
  }
}
