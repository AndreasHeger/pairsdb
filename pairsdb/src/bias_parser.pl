####
####
##
## Project PairsDB
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: bias_parser.pl,v 1.1.1.1 2002/07/02 10:46:58 heger Exp $
##
##
####
####


# parse -rdb option output from biasdb.pl 
# into format
# hid@from@to@method

use strict;

my ($method) = (@ARGV);

if (!$method) {
  warn "No method supplied!\n";
  exit(0);
}

my $i; for ($i = 1; $i <=6; $i++) { $_ = <STDIN>; }

while (<STDIN>) {
  my ($hid, $type, $from, $to, @dummy) = split (/\s+/);
  print join ("\t", $hid, $from, $to, $method, $type) . "\n";
}
