####
####
##
## Project PairsDB
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: Pairsdb_100x40.py,v 1.1.1.1 2002/07/02 10:46:57 heger Exp $
##
##
####
####

USAGE="""python map_alignments.py [OPTIONS] map_a2b map_b2c 

Map alignments given in map_a2b using alignments in 
map_b2c. The mapping has to be unique. If it is not,
an error is recorded and no mapping is performed.

For example, in the pairsdb project creating pairsdb_100x40 
applies the command:

python map_alignments.py \
        --output-format=emissions \
        --allow-partial-map \
        pairsdb_100x90.table pairsdb_90x40.table \
> pairsdb_100x40.table.

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
                              (id2, f.mColFrom + 1, f.mColTo, f.mColAlignment,
                               id1, f.mRowFrom + 1, f.mRowTo, f.mRowAlignment,
                               result.getScore(), pid ))
    options.stdout.flush()

def readMap( infile, input_format, options, invert = False ):
    """read and index mapping information."""

    if options.loglevel >= 1:
        options.stdlog.write("# reading mapping information and building index.\n" )
        options.stdlog.flush()

    t1 = time.time()

    map_pair2pos = {}
    map_a2b = {}

    while 1:

        pos = infile.tell()
        line = infile.readline()
        if not line: break

        if line[0] == "#": continue

        data = line[:-1].split("\t")

        if input_format == "pairsdb":
            a,b = data[4], data[0]
        else:
            a,b = data[0], data[1]

        pair = (min(a,b), max(a,b))
        if invert: 
            a,b = b,a
            if a not in map_a2b:
                map_a2b[a] = []
            map_a2b[a].append( b )
        else:
            if a in map_a2b: raise "duplicate: %s" % a
            map_a2b[a] = b

        map_pair2pos[ pair ] = pos

    if options.loglevel >= 1:
        options.stdlog.write("# reading mapping information finished in %i seconds - %i entries.\n" % (time.time() - t1, len(map_a2b)))
        options.stdlog.flush()

    return map_a2b, map_pair2pos

def getMap( a, b, infile, map_pair2pos, input_format ):
    """return alignment between a and b."""
    key = (min(a,b), max(a,b))
    if key not in map_pair2pos:
        return None
    else:
        infile.seek( map_pair2pos[key] )

    line = infile.readline()

    map_a2b = alignlib.makeAlignmentVector()
    
    if input_format == "pairsdb":
        rep, rep_start, rep_end, rep_ali, mem, mem_start, mem_end, mem_ali = line[:-1].split("\t")[:8]

#        if rep == b:
#            rep, rep_start, rep_end, rep_ali, mem, mem_start, mem_end, mem_ali =\
#                mem, mem_start, mem_end, mem_ali, rep, rep_start, rep_end, rep_ali

        f = alignlib.AlignmentFormatEmissions()
        try:
            f.mRowFrom = int(mem_start) - 1
            f.mColFrom = int(rep_start) - 1
        except ValueError:
            return None
        f.mRowAlignment = mem_ali
        f.mColAlignment = rep_ali
        f.copy( map_a2b )
    else:
        raise "unimplemented"

    return map_a2b

def getRepresentative( a, map_a2b, map_b2c ):
    """return identifier c that a maps to."""
    if a not in map_a2b:
        return None, None
    b = map_a2b[a]
    if b not in map_b2c:
        return b, None
    return b, map_b2c[b]

if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id$", usage=USAGE )

    parser.add_option( "-f", "--filename-fasta", dest="input_filename_fasta", type="string" ,
                       help="INPUT filename of with indexed fasta sequences.")

    parser.add_option( "-o", "--output-format", dest="output_format", type="choice",
                       choices=("blocks", "explicit", "emissions", "pairsdb" ),
                       help="alignment format for output." )

    parser.add_option( "-i", "--input-format", dest="input_format", type="choice",
                       choices=("blocks", "explicit", "emissions", "pairsdb" ),
                       help="alignment format for output." )

    parser.add_option( "--output-filename-unaligned", dest="output_filename_unaligned", type="string",
                       help="output filename for unalignable pairs." )

    parser.add_option( "--output-filename-notfound", dest="output_filename_notfound", type="string",
                       help="output filename for ids not found in sequence database." )

    parser.add_option( "--allow-partial-map", dest="partial_map", action="store_true",
                       help="if set, allow partial matches, i.e., output a->b, if no b->c found." )

    parser.set_defaults( input_filename_fasta = None,
                         output_filename_unaligned = "unaligned.table",
                         output_filename_notfound = "notfound.ids",
                         output_format = "emissions",
                         input_format = "emissions",
                         tempdir = None,
                         partial_map = False,
                         )

    (options, args) = Experiment.Start( parser, add_pipe_options = True )

    if len(args) != 2:
        print USAGE
        print "please supply two file names"
        sys.exit(1)

    fasta = IndexedFasta.IndexedFasta( options.input_filename_fasta )

    nmapped, nskipped, nerrors, nreps = (0,0,0,0)

    iteration = 0

    nnotfound1, nnotfound2 = 0, 0
    ninput, noutput, nunaligned = 0, 0, 0
    nmapped = 0

    unaligned_pairs = []
    not_found = set()

    infile_a2b = open( args[0] )
    infile_b2c = open( args[1] )

    map_a2b, map_a2pos = readMap( infile_a2b , options.input_format, options )
    map_b2c, map_b2pos = readMap( infile_b2c , options.input_format, options )
    ali_a2c = alignlib.makeAlignmentVector()
    
    infile_a2b = open( args[0] )
    infile_b2c = open( args[1] )

    for a in map_a2b.keys():
        
        ninput += 1
        iteration += 1

        if options.loglevel >= 2:
            options.stdlog.write( "# iteration %i: processing a=%s\n" % (iteration, a) )
            sys.stdout.flush()

        ## retrieve b and c
        b,c = getRepresentative( a, map_a2b, map_b2c )

        ###################################################
        # a->b
        if b == None:
            if options.loglevel >= 1:
                options.stdlog.write( "# skipping, no 'b' found for a=%s\n" % a )
            nskipped += 1
            continue

        ## obtain the alignment a2b
        ali_a2b = getMap( a, b, infile_a2b, map_a2pos, options.input_format )
        if ali_a2b == None:
            options.stdlog.write("# error for a=%s, b=%s: could not retrieve alignment from map_a2b\n" )
            sys.stdout.flush()
            nerrors += 1
            continue

        if options.loglevel >= 3:
            options.stdlog.write( "# a->b: %s\n" % str(alignlib.AlignmentFormatEmissions( ali_a2b )))
                                                                   
        ###################################################
        # b->c
        if c == None:
            if options.partial_map:
                # map a2b straight through as a2c
                ali_a2c = ali_a2b
                c = b
            else:
                if options.loglevel >= 1:
                    options.stdlog.write( "# skipping, no 'c' found for a=%s, b=%s\n" % (a,b))
                    nskipped += 1
                continue
        elif a == c:
            if options.loglevel >= 1:
                options.stdlog.write( "# skipping, because a==c\n" )
            nreps += 1
            continue
        else:
            ## obtain the alignment b2c
            ali_b2c = getMap( b, c, infile_b2c, map_b2pos, options.input_format )

            if options.loglevel >= 3:
                options.stdlog.write( "# b->c: %s\n" % str(alignlib.AlignmentFormatEmissions( ali_b2c )))

            if ali_b2c == None:
                options.stdlog.write("# error for b=%s, c=%s: could not retrieve alignment from map_b2c\n" )
                sys.stdout.flush()
                nerrors += 1
                continue

            ## combine a2b and b2c
            alignlib.combineAlignment( ali_a2c, ali_a2b, ali_b2c, alignlib.CR)
            nmapped += 1

            if options.loglevel >= 3:
                options.stdlog.write( "# a->c: %s\n" % str(alignlib.AlignmentFormatEmissions( ali_a2c )))
        
        if ali_a2c.isEmpty():
            options.stdlog.write("# error for a=%s, b=%s, c=%s: empty alignment map_a2c\n" % ( a, b, c))
            sys.stdout.flush()
            nunaligned += 1
            unaligned_pairs.append( (a,c) )
            continue

        try:
            seq_a = alignlib.makeSequence( fasta.getSequence( a ) )
        except KeyError:
            not_found.add( a )
            nnotfound1 += 1
            continue
        try:
            seq_c = alignlib.makeSequence( fasta.getSequence( c ) )
        except KeyError:
            not_found.add( c )
            nnotfound2 += 1
            continue

        outputResult( ali_a2c, a, seq_a, c, seq_c, options )
        
        noutput += 1
        
    infile_a2b.close()
    infile_b2c.close()

    if options.output_filename_unaligned and unaligned_pairs:
        outfile = open( options.output_filename_unaligned, "w" )
        for a,b in unaligned_pairs:
            outfile.write( "%s\t%s\n" % (a,b))
        outfile.close()

    if options.loglevel >= 1:
        options.stdlog.write( "# ninput=%i, noutput=%i, nskipped=%i, nmapped=%i, nerrors=%i, nunaligned=%i, nnotfound1=%i, nnotfound2=%i, nnomap=%i\n" %\
                              (ninput, noutput, nskipped, nmapped, nerrors, nunaligned, nnotfound1, nnotfound2, nreps) )
        
    Experiment.Stop()


