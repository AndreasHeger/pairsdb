package Short;

use strict;

use pretty;
use suf;
use align;

#used by Bias
my $Cutoff = 40;
my @AllAA = ( 'V','L', 'I', 'M', 'F', 'W', 'Y', 'G', 'A', 'P', 'S', 'T', 'C', 'H', 'R', 'K', 'Q', 'E','N','D');

# used by Short
my $verbose      =  0;
my $gapopen      =-11; 
my $gapelon      = -4;
my $blastube     =  5;
my $mintubescore = 60;
my $ktup         = 3;
my $maxdiagonal  = 20;
my $execpath     = '/net/cpp-data/backup/andreas/projects/pairsdb/src/annabelle';

#-----------------------Code----------------------------------------------
sub Dump {
  print "---- Dump: Short Version 1.0 -----\n";
  print "Cutoff         = $Cutoff\n";
  print "mintubescore   = $mintubescore\n";
  print "gapopen        = $gapopen\n";
  print "gapelon        = $gapelon\n";
  print "blastube       = $blastube\n";
  print "ktup           = $ktup\n";
  print "maxdiagonal    = $maxdiagonal\n";
  print "verbose        = $verbose\n";
  print "---- End Dump -----\n";
}


############################################################################
############################################################################
# Title: DoShort
# Function: looks for short repeats, i.e. in the range of ktup-size
#---------------------------------------------------------------------------
sub DoShort {
	my($seq,$hid)=@_;
	my($result)=''; my($irepeat)=1;
	my($flag)=1;
	while($flag) {
		$flag=0;
		my($x, $lrepeat, $score, $diagonal, $mali, $from, $to)=short($seq,$hid); 
		last if($x<1);
		my(@mali)=split(/\n/,$mali); shift(@mali);
		$result.= sprintf("---------------------------------------------------------------------------\n");
		$result.= sprintf("No. of Repeats|Total Score|Length  |Diagonal| BW-From|   BW-To|   Level\n");
		$result.= sprintf("%14i|%11.2f|%8i|%8i|%8i|%8i|   Short\n", $#mali + 1, $score, $lrepeat, $diagonal, $from, $to);
		$result.= sprintf("---------------------------------------------------------------------------\n");
		foreach (@mali) {
		  $result.= $_ . "\n";


		  my($from,$to)=/^\S+\/(\d+)\-(\d+)/;
		  if ($verbose) { 
		    print "Masking $from-$to\n";
		  }
		  next if($to-$from>$lrepeat*2);
		  substr($seq,$from-1+$[,$to-$from+1)=~s/./X/g;
		  $flag=1;
		}

		if ($verbose) {
		  print $result;
		}

	}
	return($seq,$result);
}

sub short {
  my($seq,$hid)=@_;
  my($icol,@w,%x,$word,@occ,$keep,%hash);

  undef(%hash); 
  my(@seq)=split(//,$seq); unshift(@seq,'?');
  my $lseq=$#seq+1-$[;
  my(@w);

  # 1. create hash of words of length ktup 
  foreach (1..$ktup-1) { push(@w,shift(@seq)); }
  foreach $icol (1..$lseq-$ktup+1) {
    shift(@w); push(@w,shift(@seq));
    $word=join('',@w);
    next if($word=~/[BXZ\W\.]/);	# skip gaps, ambiguous chars
    $hash{$word}.="$icol ";		# create hashs for words of length ktup
    $occ[$icol]++;			# words per column
  }

  # 2. insert information for each word: column:#occurances; 
  # there can be several occuarences of the same word in the same column, for exampel if ktup is 4 and sequence is AAAAAAA
  foreach $word (keys %hash) {		# word column counts 
    foreach $icol (split(/\s+/,$hash{$word})) {
      $x{$icol}++;
    }
    my $x=''; my $keep=0;
    foreach $icol (sort { $a <=> $b } keys %x) {
      my $f="$x{$icol} ";
      $x.="$icol:$f"; 
      $keep=1;
    }
    if($keep) { $hash{$word}=$x; } else { delete $hash{$word}; }
    undef(%x); 
  } 
  # 3. create histogram of words
  my(%hist,$x,$y,$d,@diagonals,$k,$l);
  undef(%hist);
  foreach (keys %hash) {
    my(@x)=sort { $a <=> $b } split(/\s+/,$hash{$_});
    while($x[$[]) {
      ($x,$k)=split(/:/,shift(@x));
      foreach (@x) {
	($y,$l)=split(/:/);
	$d=$y-$x; last if($d>$maxdiagonal);
	$hist{$d}+=$k*$l;
      }
    }
  }
  
  foreach (sort {$hist{$b}<=>$hist{$a}} keys %hist) {
    print "diagonal\t$_\t$hist{$_}\n" if ($verbose); 
    push(@diagonals,$_);
  }
  
  # test best diagonal for tubescore
  my(@profile)=suf::blosumprofile($seq); 
  my(%dot); undef(%dot);
  my($diag)=shift(@diagonals);
  my $d1=$diag/4; if($d1>$blastube) { $d1=$blastube; }
  $d1=0; # ungapped diagonal

  # get all dots in diagonal
  my($n)=&fulldiag(1,$lseq,\%dot,$lseq,$diag-$d1,$diag+$d1,$seq);
  print "$n dots by fulldiag\n" if ($verbose);

  # score dots in diagonal against profile
  my(@col,@row,@evalue);
  ($n)=align::massagedots_profile_seq(\@col,\@row,\@evalue,\%dot,\@profile,$seq,1,1,$lseq,$lseq,-$lseq,$lseq);
  print "$n dots after massage\n" if ($verbose); 
#  my $i;
#  for $i (sort {$row[$a] <=> $row[$b]} (0..$#row) ) {
#    printf ("%5i %5i %5.2f|", $col[$i], $row[$i], $evalue[$i]); 
#  }
#  print "\n";
  
  my ($from, $to) = FindPositiveRegion( \@col, \@row, \@evalue, $lseq, $diag);
  print "Best positive region: $from-$to\n" if ($verbose);
  
  
  # align around dots
  my($tubescore,@ali)=align_abc(\@row,\@col,\@evalue,0,$from,$to,$from,$to,$gapopen,$gapelon,$gapopen,$gapelon); 
  print "tubediag $diag score1 $tubescore ($mintubescore)\n" if ($verbose);
  if($tubescore<$mintubescore) { return(0); }

  #  foreach $i (1..$#ali) {  printf ("%5i %5i|", $i ,$ali[$i]) }; printf ("\n");
  
  # build multiple alignment
  my(@ali1)=foldali_short($diag,$lseq,@ali);			# because no gaps

  # my $i; foreach $i (1..$#ali1){  printf ("%5i %5i|", $i ,$ali1[$i]) }; printf ("\n");
  my($mali)=pretty::pretty($seq,$hid,1,$diag,1,@ali1);
  
  return(1, $diag, $tubescore, $diag, $mali, $from, $to);
}

sub fulldiag {
  my($from2,$to2,$dot,$lseq,$mind,$maxd,$seq)=@_;
  my($ires,$jres,$n,$key);
  #print "* * * fulldiag $from2 $to2 $mind $maxd * * * \n";
  $n=0; 
  my(@x)=split(//,$seq); if($[<1) { unshift(@x,'?'); }
  foreach $ires ($from2..$to2) {
    next if($x[$ires] =~ /X/); # masked
    foreach $jres ($ires+$mind..$ires+$maxd) { 
      next if($jres<1); next if($jres>$lseq);
      next if($x[$jres] =~ /X/); # masked
      $n++;
      $key="$ires $jres";
      $dot->{$key}=1;
    }
  }
  return($n);
}

sub foldali_short {
        # assumes @ali is internally consistent, from tube alignment
        my($window,$nres,@ali)=@_;
        my($n,%a,$last,%left,%rite,%order,$i,$j,$k,$where);
        $n=0; undef(%a); undef(%left); undef(%rite); $last=0;
        # stop at long unaligned region
        $i=1; while($ali[$i]<1&&$i<=$nres) { $i++; } $j=$i; if($i>=$nres) { return(-1); } # no more alis
        while($i-$j<$window && $i<$nres) {
	  $j=$i; while($ali[$j]>0&&$j<=$nres) { $j++; }
	  $i=$j; while($ali[$i]<1&&$i<=$nres) { $i++; }
        } $where=$i-1; 
	# print "foldali: profile for range 1..$where ($i $j)\n";
        foreach $j (1..$where) {
                next if($ali[$j]<1); 
                $i=$ali[$j]; 
                if($a{$j}>0) { # add onto old column
                        $a{$i}=$a{$j};
                } else { # insert new column
                        $n++; 
                        $a{$i}=$n; $a{$j}=$n; 
                        # insert ali[i]->i-->find leftside defined column
                        $k=$j-1; while($a{$k}<1 && $k>0) { $k--; } 
                        if($k>0) { $last=$a{$k}; $k=$rite{$last}; } 
                        else { $last=undef; $k=undef; }
                        $left{$n}=$last; if($last>0) { $rite{$last}=$n; }
                        $rite{$n}=$k; if($k>0) { $left{$k}=$n; }
                }
                # print "$i -> column $a{$i} left $left{$a{$i}} rite $rite{$a{$i}} last $last\n";
	      }
        $i=$n; while($rite{$i}>0) { $i=$rite{$i}; } $k=$n;
	# print "n = $n\n";

        while($i>0 && $k>0) {
	  # print "i=$i k=$k\n";
	  $order{$i}=$k; $k--;
	  $i=$left{$i}; 
        }
        foreach $i (1..$nres) {
	  $j=$order{$a{$i}}; 
	  if($j<1) { $j=0; }
	  # print "$i $a{$i} $order{$a{$i}}\n";
	  $ali[$i]=$j;
        }
        return(@ali);
}

sub align_abc {
  
  my ($row, $col, $score, $roll, $from1, $to1, $from2, $to2, $gop_i, $gep_i, $gop_j, $gep_j ) = @_;
  
  open (OUT, ">$$.dot") or die "Could not open $$.dot\n";
  my $count = $#$row + 1;

  print OUT pack ("i$count", @$row);
  print OUT pack ("i$count", @$col);
  print OUT pack ("d$count", @$score);

  close (OUT);     

  if ($verbose) {
    print "$execpath -A $gop_i -B $gep_i -C $gop_j -D $gep_j -G $from1 -H $to1 -I $from2 -J $to2 -K $roll -M dotalign -P $$.dot\n";
  }

  my (@lines) = `$execpath -A $gop_i -B $gep_i -C $gop_j -D $gep_j -G $from1 -H $to1 -I $from2 -J $to2 -K $roll -M dotalign -P $$.dot`;

  if ($verbose) {
    print @lines;
  }

  $_ = shift (@lines);
  my ($swscore, $lali, $ngaps) = /Alignment: score=\s*([\d\.]+) length=(\d+) gaps=(\d+)/;
  my (@ali, @scores);
  my @ali;
  my @scores;
  foreach (@lines) {
    my ($i, $j, $x) = /^\s+(\d+)\s+(\d+)\s+([\-\d\.]+)/;
    $ali[$i]    = $j;                         # remap profile segment, changed 8.1.99, because change in dotalign;
  }
  system("rm -f $$.dot");
  
  return ($swscore, @ali );
} 


sub FindPositiveRegion {
  
  my ($col, $row, $evalue, $lsequence, $diag) = @_;

  # find region of positive window
  # note: evalue is look-ahead-value!!!, i.e. gives the score i -> i+diag

  my @scores;
  my $i;

  foreach $i (0..$#$col) { $scores[$col->[$i]] = $evalue->[$i]; }  
  
  #  foreach $i (1..$#scores) { printf("%5i %5.2f|", $i, $scores[$i]); } printf("\n");
  my $bestscore = 0;
  my $s = 0;
  my $bestto;
  my $bestfrom = 0;
  my $currentfrom = 1;

  for $i (1..$diag-1) { $s += $scores[$i] }; # calculate score of first window
  my $totalscore = $s;

  my $inregion;
  if ($s > 0) { $inregion = 1; } else { $inregion = 0; }

  $i = 1;
  while ($i <= $lsequence - $diag) {			# i is starting point of window
    $s += $scores[$i + $diag - 1 ] - $scores[$i - 1];   # calculate score of current window starting at i
#    print "$i $s $totalscore $scores[$i]";
    if ($inregion) {
      if ($s < 0) {					
#	print "leaving reagion ";
	$inregion = 0;
	if ($totalscore > $bestscore) {
	  $bestfrom = $currentfrom + 1;  # hmm, looks nicer
	  $bestto   = $i + $diag - 1;	 # include complete region pointed at by previous complete window
	  $bestscore= $totalscore;
#	  print "saving region $bestfrom $bestto $bestscore";
	}
      } else {
	$totalscore += $scores[$i];
      }
    } else {
      if ($s >= 0) { # starting a new window 
	$inregion = 1;
	$currentfrom = $i;
	$totalscore  = 0;
#	print "entering region starting $currentfrom"
      }
    }

    $i++;
#    print "\n";
  }
	
  return( $bestfrom, $bestto);
}



