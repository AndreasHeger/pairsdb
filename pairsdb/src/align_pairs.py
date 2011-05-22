####
####
##
## Project PairsDB
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id$
##
##
####
####

USAGE="""python align_pairs.py [OPTIONS] < stdin

given a list of sequence pairs, compute a pairwise alignment between them. 

stdin: the list of pairs as tab-separated identifiers. Lines starting with
        '#' are ignored.
"""

import sys, re, string, optparse, time, base64, md5, os, tempfile, shutil

import alignlib
import IndexedFasta
import Experiment

def outputResult( result, 
                  id1, seq1, 
                  id2, seq2, 
                  options ):
    """output the result."""

    pid = alignlib.calculatePercentIdentity( result, seq1, seq2) * 100.0

    if options.output_format == "emissions":
        options.stdout.write( "%s\t%s\t%5.2f\t%s\n" % (id1, id2, pid, str(alignlib.AlignmentFormatEmissions( result ) ) ) )
    elif options.output_format == "blocks":
        options.stdout.write( "%s\t%s\t%5.2f\t%s\n" % (id1, id2, pid, str(alignlib.AlignmentFormatBlocks( result ) ) ) )
    elif options.output_format == "explicit":
        options.stdout.write( "%s\t%s\t%5.2f\n%s\n" % (id1, id2, pid, str(alignlib.AlignmentFormatExplicit( result, seq1, seq2 ) ) ) )
    elif options.output_format == "pairsdb":
        f = alignlib.AlignmentFormatEmissions()
        f.fill( result )
        options.stdout.write( "%s\t%i\t%i\t%s\t%s\t%i\t%i\t%s\t%5.2f\t%5.2f\n" %\
                              (id1, f.mRowFrom + 1, f.mRowTo, f.mRowAlignment,
                               id2, f.mColFrom + 1, f.mColTo, f.mColAlignment,
                               result.getScore(), pid ))
    options.stdout.flush()

if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id$", usage=USAGE )

    parser.add_option( "-f", "--filename-fasta", dest="input_filename_fasta", type="string" ,
                       help="INPUT filename of with indexed fasta sequences.")

    parser.add_option( "-a", "--alignator", dest="alignator", type="choice",
                       choices=("groupies", ),
                       help="alignator to use for alignment." )
    parser.add_option( "-o", "--output-format", dest="output_format", type="choice",
                       choices=("blocks", "explicit", "emissions", "pairsdb" ),
                       help="alignment format for output." )
    parser.add_option( "--output-filename-unaligned", dest="output_filename_unaligned", type="string",
                       help="output filename for unalignable pairs." )
    parser.add_option( "--output-filename-notfound", dest="output_filename_notfound", type="string",
                       help="output filename for ids not found in sequence database." )

    parser.set_defaults( input_filename_fasta = None,
                         output_filename_unaligned = "unaligned.pairs",
                         output_filename_notfound = "notfound.ids",
                         output_format = "emissions",
                         alignator = "groupies",
                         tempdir = None,
                         )

    (options, args) = Experiment.Start( parser, add_pipe_options = True )

    fasta = IndexedFasta.IndexedFasta( options.input_filename_fasta )

    if options.alignator == "groupies":
        alignator = alignlib.makeAlignatorGroupies()
    else:
        raise "not implemented."

    result = alignlib.makeAlignmentVector()

    nnotfound1, nnotfound2 = 0, 0
    ninput, noutput, nunaligned = 0, 0, 0

    unaligned_pairs = []
    not_found = set()

    for line in sys.stdin:
        
        if line[0] == "#": continue

        ninput += 1

        id1, id2 = line[:-1].split( "\t" )

        try:
            seq1 = alignlib.makeSequence( fasta.getSequence( id1 ) )
        except KeyError:
            not_found.add( id1 )
            nnotfound1 += 1
            continue
        try:
            seq2 = alignlib.makeSequence( fasta.getSequence( id2 ) )
        except KeyError:
            not_found.add( id2 )
            nnotfound2 += 1
            continue
        
        if options.loglevel >= 3:
            options.stdlog.write("# aligning:\n# %s\t%s\n# %s\t%s\n" % (id1, seq1, id2, seq2) )
            options.stdlog.flush()

        alignator.align( result, seq1, seq2 )

        if result.getLength() == 0:
            nunaligned += 1
            unaligned_pairs.append( (id1, id2) )
            continue

        outputResult( result, id1, seq1, id2, seq2, options )
        noutput += 1

    if options.output_filename_unaligned and unaligned_pairs:
        outfile = open( options.output_filename_unaligned, "w" )
        for a,b in unaligned_pairs:
            outfile.write( "%s\t%s\n" % (a,b))
        outfile.close()

    if options.output_filename_notfound and not_found:
        outfile = open( options.output_filename_notfound, "w" )
        for a in not_found:
            outfile.write("%s\n" % (a))
        outfile.close()
        
    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nunaligned=%i, nnotfound1=%i, nnotfound2=%i\n" % (ninput, noutput, nunaligned, nnotfound1, nnotfound2) )

    Experiment.Stop()
