#!/bin/env perl5
####
####
##
## Project PairsDB
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: biasdb.pl,v 1.1.1.1 2002/07/02 10:46:58 heger Exp $
##
##
####
####


#
#     this script detects COMPOSITION BIASed regions in proteins
#     using the algorithm 'biasdb' 
#    
#     Project: biasdb                     Author:  G.Casari Dec 1994
#	
#	4.12.2003	Andreas Heger
#			there is the chance of an endless-loop, when
#			the score for an amino acid is the largest, but
#			the amino acid itself is not present. This can
#			happen because of off-diagonal positive entries
#			in the substitution matrix.
#			Workaround: check, if a residue is actually present
#				in the segment to mask.
#			(checking for last run/excluding last aa will not work,
#			because there might be several regions worth masking)
#			
#======================================================================
$Cutoff = 40;
@AllAA = ( V, L, I, M, F, W, Y, G, A, P, S, T, C, H, R, K, Q, E, N, D);
{
    local($opt);
    for($opt=0;$opt<=$#ARGV;){
	if($ARGV[$opt] eq "-fasta"){
	    $Fasta = 1;
	    splice(@ARGV,$opt,1);
	    next;
	}elsif($ARGV[$opt] eq "-rdb"){
	    $Rdbflag = 1;
	    splice(@ARGV,$opt,1);
	    next;
	}elsif($ARGV[$opt] eq "-m"){ 
	    $Matrixfile = $ARGV[$opt+1];
	    splice(@ARGV,$opt,2);
	    next;
	}elsif($ARGV[$opt] eq "-cut"){
	    die "invalid argument, give numeric value following -cut option" 
	      unless ($ARGV[$opt+1] =~ /^[\d.]+$/);
	    $Cutoff = $ARGV[$opt+1];
	    splice(@ARGV,$opt,2);
	    next;
	}elsif($ARGV[$opt] eq "-gq"){
	    splice(@ARGV,$opt,2);
	    next;
	}elsif($ARGV[$opt] eq "-help" || $ARGV[$opt] eq "-h"){
	    print "usage: $0 sequencefile [options]\n";
	    print "        options: -fasta   ... sequence in fasta format\n";
	    print "                 -help    ... print this text\n";
	    print "                 -rdb     ... print biased regions as rdb file\n";
	    print "                 -cut nn  ... set the cutoff score for reported regions\n";
	    print "                              default is $Cutoff\n";
	    print "                 -m file  ... use matrix in 'file' as comparison matrix\n";
	    print "                              format of matrix: blast comparison matrix\n";
	    print "                              default uses internal blosum62\n";
	    exit;
	}else{
	    $opt++;
	}
    }
}

{				# read comparison matrix
    local($values)=0;
    local($a,@b,@x,$i);

    if($Matrixfile){
	open(DATA,$Matrixfile) || die "could not open matrix file '$Matrixfile'\n";
    }

    while(<DATA>){
	chop;
	if($values){
	    ($a,@x)= split(/\s+/,$_);
	    for($i=0;$i<=$#x;$i++){
		$Matrix{$a.$b[$i]}=$x[$i];
	    }
	}elsif(/^\#/){
	    next;
	}else{
	    @b = grep(!/^\s*$/,split(/\s+/,$_));
	    $values = 1;
	}
    }
}

				# read sequence
#if ($Fasta){
#    die "cannot read sequence" unless (($firstline,@sequence)= <>);
#}else{
#    die "cannot read sequence" unless (($firstline, $secondline,@sequence)= <>);
#}
#$Sequence = join('',@sequence);
#$Sequence =~ s/\s//g;
#$Sequence =~ s/\d//g;
#$Sequence =~ s/^(.*\*).*$/$1/;
#$Sequence =~ tr/a-z/A-Z/;
&prepareRdb('dummy') if($Rdbflag);

#while(1){			# detect all biased regions
# read multisequence fasta-file from STDIN
while(<>) {
  if(/^>/) { 
   $line=$_; 
   if($toprint) { &printit; }
    # new sequence
    $toprint=1;
    $firstline=$line;
    $Sequence='';
  } else {
    $Sequence.=$_;
  } 
}

if($toprint) { &printit; }

sub printit {
    # print old
$Sequence =~ s/\s//g;
$Sequence =~ s/\d//g;
$Sequence =~ s/^(.*\*).*$/$1/;
$Sequence =~ tr/a-z/A-Z/;
$OriginalSequence=$Sequence;
    if($firstline =~ /^>P1;\s*(\S+)\s*/){
	$Seqname = $1;
    }else{
	($Seqname) = ($firstline =~ /^>\s*(\S+)\s*/);
    }
   while(1) {
    $Mscore = 0;
    $maa = '';
    local($Seqlength) = length($Sequence);
    foreach $AA (@AllAA){
	($score,$from,$to,$segment)=&checkaa($AA,$Sequence,$Seqlength);
	if ($score >= $Cutoff){
	    if($Mscore < $score && substr($Sequence, $from, $to+1-$from) =~ /$AA/) {	
		# AH: add check for presence of current AA 
		# note: check full sequence and not just segment (why?)
		# otherwise: infinite loop for 172868
		$maa = $AA;
		$mfrom = $from;
		$mto = $to;
		$Mscore = $score;
	    }	    
	}
    }
    if($maa){
	if($Rdbflag){
	    &reportRdb($maa,$mfrom, $mto,$Mscore,$Seqlength);
	}else{
	    print STDERR "$maa-rich region ", $mfrom+1 ," to ",$mto+1," corrected\n";
	}
	substr($Sequence,$mfrom, $mto+1-$mfrom) =~ s/$maa/X/g;
    }
    last unless $maa;
   }
unless($Rdbflag){
    print $firstline;
    unless ($Fasta){
	print $secondline;
    }
    for ($pos=0;$pos<length($Sequence);$pos+=60){
	print substr($Sequence,$pos,60),"\n";
    }
}
}

sub checkaa{
# --------------------------------------------------------------------------
    local($aa,$sequence,$len)=@_;
    local($s,$i,$im);
    local(@scores) = &run_SW_VA($aa,$sequence,$len);
    $s = 0;
    for($i = 0;$i<$len;$i++){
	if($s <= $scores[$i]){
	    $s = $scores[$i];
	    $TO = $i;
	}
    }
    for($i=$TO-1;$i>=0;$i--){last if($scores[$i]<=0);}
    $FR = $i+1;
    return ($scores[$TO],$FR,$TO,substr($sequence,$FR,$FR-$TO+1));
}





sub run_SW_VA{
#-----------------------------------------------------------------------------
    local($aa,$sequence,$len)=@_;
    local($seqcopy);
    local(@scores,$pos,$s);
    $#scores = $len-1; # initialize?

				# forward
    $scores[0]= $s = $Matrix{$aa.substr($sequence,0,1) };
    if ($s < 0){ $s = 0;}

    for($pos = 1;$pos < $len; $pos++){
	$scores[$pos]= $s += $Matrix{$aa . substr($sequence,$pos,1)};
	$s = 0 if($s < 0);
    }

    return @scores;
}
	
sub prepareRdb{
#--------------------------------------------------------------------------
    ($OriginalSequence) = @_;
    local(@rdbheaders)= ("sequence","bias","from","to","score","fragment","seqleng","nres","fragleng");
    local(%rdbformats)= ("sequence", "20S",
			 "bias",     "4S",
			 "from",     "8N",
			 "to",       "8N",
			 "score",    "8N",
			 "seqleng",  "8N",
			 "fragleng", "8N",
			 "nres",     "4N",
			 "fragment", "10S");
    local($rdbcomment)= "\# Perl-RDB\n".
	"\#composition biased regions detected with program biasdb \n".
	    "\#      pleas cite    .....   G.Casari 1995\n".
		"\#\n";
    local(@format,$key);
    foreach $key (@rdbheaders){ push(@format,$rdbformats{$key});}
    $RdbTop = $rdbcomment . join("\t",@rdbheaders)."\n".
		join("\t",@format)."\n";
}

sub reportRdb{
#------------------------------------------------------------------------
    local($maa,$mfrom, $mto,$mscore,$len)=@_;
    local($nx);
    if($RdbTop){
	print $RdbTop;
	$RdbTop = '';
    }
    $nx= substr($OriginalSequence,$mfrom, $mto+1-$mfrom);
    eval("\$nx =~ s/[^$maa]//g");
    $nx = length($nx);
    print join("\t",$Seqname,$maa,$mfrom+1,$mto +1,$mscore,
#	       substr($OriginalSequence,$mfrom, $mto+1-$mfrom),
		,$len,$nx,$mto+1-$mfrom),
		"\n";
}

1;

__END__
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1 -.3 -4 
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0 -.3 -4 
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
X -.3 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2 -.3 0 -2 -1 -1 -1 -1 -1 -4 
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 




