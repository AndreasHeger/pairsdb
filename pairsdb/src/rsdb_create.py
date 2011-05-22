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

USAGE="""python rsdb_create.py [OPTIONS] blast_graph

Toy implementation for rsdb create.

"""

import sys, re, string, optparse, time, base64, md5, os, tempfile, shutil

import Experiment
import IndexedFasta

def buildIndex( infile ):

    map_nid2index = {}

    last_nid = None
    while 1:
        pos = infile.tell()

        line = infile.readline()
        if not line: break
        
        nid = line[:-1].split("\t")[0]
        if nid != last_nid:
            map_nid2index[nid] = pos
            last_nid = nid
        
    return map_nid2index

if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id$", usage=USAGE )

    parser.add_option( "-g", "--filename-graph", dest="input_filename_graph", type="string" ,
                       help="INPUT filename of pairsdb graph.")

    parser.add_option( "-f", "--filename-fasta", dest="input_filename_fasta", type="string" ,
                       help="INPUT filename of indexed fasta file.")

    
    parser.set_defaults( input_filename_graph = "pairsdb_90x90.table",
                         input_filename_fasta = None,
                         level = 40,
                         output_filename_pairsdb = "pairsdb_90x40.table",
                         output_filename_rsdb = "nrdb40.table",
                         min_coverage = 0.9,
                         max_missing = 30,
                         )

    (options, args) = Experiment.Start( parser, add_pipe_options = True )

    infile_graph = open(options.input_filename_graph, "r" )

    map_nid2index = buildIndex( infile_graph )

    fasta = IndexedFasta.IndexedFasta( options.input_filename_fasta )

    map_nid2length = fasta.getContigSizes()

    reps = map_nid2index.keys()
    reps.sort( lambda x,y: cmp(map_nid2length[x], map_nid2length[y]) )
    reps.reverse()

    reps_set = set(reps)

    outfile_reps = open( options.output_filename_rsdb, "w" )
    outfile_pairsdb = open( options.output_filename_pairsdb, "w" )

    nreps, nmems = 0, 0

    for rep_nid in reps:

        if rep_nid not in reps_set: continue

        nreps += 1
        outfile_reps.write( "%s\n" % str(rep_nid ))
        reps_set.remove(rep_nid)
        
        infile_graph.seek( map_nid2index[rep_nid] )

        while 1:
            line = infile_graph.readline()
            if not line: break
            
            (query_nid, sbjct_nid, evalue, query_from, query_to, query_ali,
             sbjct_from, sbjct_to, sbjct_ali, score, pid) = line[:-1].split("\t")
            
            if query_nid != rep_nid: 
                break
            
            if query_nid == sbjct_nid: continue
            if sbjct_nid not in reps_set: continue
            if float(pid) < options.level: continue
            
            sbjct_length = map_nid2length[sbjct_nid]
            l = int(sbjct_to) - int(sbjct_from)
            if float(l) / float(sbjct_length) < options.min_coverage or \
                    sbjct_length - l > options.max_missing:
                continue
            
            outfile_pairsdb.write( "\t".join( (\
                        query_nid, query_from, query_to, query_ali,
                        sbjct_nid, sbjct_from, sbjct_to, sbjct_ali,
                        score, pid )) + "\n" )

            reps_set.remove( sbjct_nid )
            nmems += 1

    if options.loglevel >= 1:
        options.stdlog.write( "# reps=%i, mems=%i\n" % (nreps, nmems) )

    Experiment.Stop()
