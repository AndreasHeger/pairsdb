####
####
##
## Project PairsDB
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: NrdbClassification.py,v 1.1.1.1 2002/07/02 10:46:57 heger Exp $
##
##
####
####

#----------------------------------------------------------------
import sys, re, string, tempfile, shutil, gzip, time, os, subprocess

from types import *

import SaryFasta
import IndexedFasta
import FastaIterator
import Tools

import alignlib
import map_alignments

class Counts:
    def __init__( self ):
        self.mNUnits = 0
        # shoud get converted to counts during construction.
        self.mNSequences = set()
        self.mNResidues = 0
        
class NrdbAnnotations:
    """class for mapping annotations onto nrdb sequences.
    """ 

    mName = "generic"

    def __init__ (self, options = None):

        self.mTempdirLocation = None
        self.mTempdir = None

        if options:
            self.mOptions = options
            self.mLogLevel = options.loglevel
            self.mStdout = options.stdout
            self.mStdlog = options.stdlog
            self.mStderr = options.stderr        
            
            self.mOutputFilenameAnnotations = options.output_filename_annotations
            self.mOutputFilenameDescriptions = options.output_filename_descriptions
            self.mInputFilenameMap100x40 = options.input_filename_map100x40
            self.mInputFilenameMap100x90 = options.input_filename_map100x90
            
            ## whether or not to test
            self.mTest = options.test

            ## whether or not to check for previous entries
            self.mCheck = 1

            ## whether or not to use blast to find sequences
            self.mPerformBlast = 0

            ## whether or not to ignore pids to find sequences
            self.mIgnorePid = False

            ## whether or not to split discontinous domains
            self.mSplitDiscontinuousDomains = False

            ## if splitting discontinuous domains, do so at gaps of size #
            self.mSplitAtGaps = 20

            ## do not keep empty families
            self.mKeepEmptyFamilies = False

            ## whether or not to use multiple entries if they exist
            self.mUseMultipleEntries = False

            ## whether or not to take first match
            self.mFirstMatch = True

            ## whether or not to realign domain to sequence
            self.mPerformAlignment = False

            self.mMinCoverage = 0.75
            
            ## location of the temporary directory
            self.mTempdirLocation = options.tempdir
            
            ## filenames for files containing optional information
            self.mFilenameReferences = options.input_filename_references
            self.mFilenameNrdb = options.input_filename_nrdb
            
            if options.input_filename_sary:
                self.mSary = SaryFasta.SaryFasta( options.input_filename_sary )
            else:
                self.mSary = None

            if options.input_filename_fasta:
                self.mFasta = IndexedFasta.IndexedFasta( options.input_filename_fasta )
            else:
                self.mFasta = None

            ## set to true, if annotations shall be propagated up to nrdb100
            self.mPropagateUp = False

            ## options for makeNonRedundantClone():
            ## shorten overlapping domains
            self.mNRShortenDomains = False

            ## minimum overlap between adjacent annotations
            self.mNRMinOverlap = 0

            ## mininum size of an annotation to be retained
            self.mNRMinAnnotationLength = 10 

        ## level and suffix to choose for nrdb100 level
        self.mParams100 = {'level': '', 'suffix': '', 'name' : self.mName }

        ## when to output progress report
        self.mReportStep = 1000

        ## maps
        self.mMapHid2Nid = None
        self.mMapPid2Nid = None
        self.mMapAcc2Nid = None
        self.mMapNid2Length = None

    #--------------------------------------------------------------------------------------
    def __del__(self):
        if self.mTempdir:
            shutil.rmtree( self.mTempdir )

    def addOptions( self, parser ):
        """add options to commandline parser."""

        parser.add_option( "-t", "--tempdir", dest="tempdir", type="string" ,
                           help="temporary directory use. Default is system default [/tmp].")

        parser.add_option( "--test", dest="test", type="int" ,
                           help="stop after x entries (for testing purposes).")

        parser.add_option( "--output-filename-annotations", dest="output_filename_annotations", type="string" ,
                           help="output filename for annotations.")

        parser.add_option( "--output-filename-descriptions", dest="output_filename_descriptions", type="string" ,
                           help="output filename for descriptions.")

        parser.add_option( "--input-filename-sary", dest="input_filename_sary", type="string" ,
                           help="input file with sary database.")

        parser.add_option( "--input-filename-fasta", dest="input_filename_fasta", type="string" ,
                           help="input file with indexed fasta file.")

        parser.add_option( "--input-filename-map100x90", dest="input_filename_map100x90", type="string" ,
                           help="input file with mapping information for mapping 100 -> 90.")

        parser.add_option( "--input-filename-map100x40", dest="input_filename_map100x40", type="string" ,
                           help="input file with mapping information for mapping 100 -> 40.")

        parser.add_option( "--input-filename-nrdb", dest="input_filename_nrdb", type="string" ,
                           help="input filename of the nrdb table.")

        parser.add_option( "--input-filename-references", dest="input_filename_references", type="string" ,
                           help="input file with references to look up sequences by accession number or id.")

        parser.add_option( "--skip-non-redundant", dest="make_non_redundant", action="store_false",
                           help="skip redundancy removal step.")

        parser.set_defaults(
            make_non_redundant = True,
            input_filename_nrdb = "nrdb.table.gz",
            input_filename_fasta = None,
            input_filename_sary = None,
            input_filename_references = "references.table.gz",
            input_filename_map100x40 = "pairsdb_100x40.table",
            input_filename_map100x90 = "pairsdb_100x90.table",
            output_filename_annotations  = "nrdb%(level)s_%(name)s_domains%(suffix)s.table.gz",
            output_filename_descriptions = "nrdb%(level)s_%(name)s_families%(suffix)s.table.gz",
            tempdir = "/tmp",
            test = None )

    #--------------------------------------------------------------------------------------
    def prepare( self ):
        pass

    #--------------------------------------------------------------------------------------
    def finish( self ):
        pass

    #--------------------------------------------------------------------------------------
    def resetStats( self ):
        """reset summary statistics for parsing."""
        self.mFound = {}

        self.mFromPid = 0
        self.mFromHid = 0
        self.mFromBlast = 0
        self.mFromSary = 0
        self.mSkipped = 0
        self.mNotFound = 0
        self.mFromLookup = 0
        self.mFailed = 0
        self.mWarnings = 0

        self.mStartTime = time.time()

    #--------------------------------------------------------------------------------------
    def printStats( self ):
        """print summary statistics for parsing."""

        if self.mLogLevel >= 1:
            self.mStdlog.write( "# lookup statistics: time=%i, pid=%i, hid=%i, sary=%i, blast=%i, skipped=%i, notfound=%i, lookup=%i, failed=%i, warnings=%i\n" %\
                  ( time.time() - self.mStartTime ,
                    self.mFromPid, self.mFromHid, self.mFromSary, self.mFromBlast, self.mSkipped,
                    self.mNotFound, self.mFromLookup, self.mFailed, self.mWarnings) )

    #--------------------------------------------------------------------------------------
    def readHids( self ):
        """build map of hid 2 nid from nrdb table.
        """

        if self.mLogLevel >= 1:
            self.mStdlog.write("# reading hids from %s\n" % self.mFilenameNrdb) 
            self.mStdlog.flush()
            
        t1 = time.time()

        infile = self.openInputFile( self.mFilenameNrdb )
        
        self.mMapHid2Nid = {}

        ninput, ninconsistencies, nerrors = 0, 0, 0

        for line in infile:
            if line[0] == "#": continue
            
            data = line[:-1].split("\t")
            
            try:
                nid, filter = map( int, (data[1], data[10]))
                hid = data[2]
            except ValueError:
                self.mStderr.write( "parsing error in file %s: line %s\n" % ( self.mFilenameNrdb, line[:-1] ) )
                nerrors += 1
                continue
            
            if filter == 0:
                continue

            if hid in self.mMapHid2Nid:
                self.mStdlog.write( "# warning: duplicate and inconsistent hid '%s'.\n" % hid)
                ninconsistencies += 1

            self.mMapHid2Nid[hid] =  nid

        infile.close()

        if self.mLogLevel >= 1:
            self.mStdlog.write("# read hids (%i hids) in %i seconds.\n" % ( len(self.mMapHid2Nid),
                                                                            time.time() - t1))
            self.mStdlog.flush()

    #--------------------------------------------------------------------------------------
    def readLengths( self ):
        """build map of nid2length from nrdb table.
        """

        if self.mLogLevel >= 1:
            self.mStdlog.write("# reading lengths from %s\n" % self.mFilenameNrdb) 
            self.mStdlog.flush()
            
        t1 = time.time()

        infile = self.openInputFile( self.mFilenameNrdb )
        
        self.mMapNid2Length = {}

        for line in infile:
            if line[0] == "#": continue
            
            data = line[:-1].split("\t")
            
            try:
                nid, length, filter = map( int, (data[1], data[9], data[10]))
            except ValueError:
                self.mStderr.write( "parsing error in file %s: line %s\n" % ( self.mFilenameNrdb, line[:-1] ) )
                continue
            
            if filter == 0:
                continue

            if nid in self.mMapNid2Length:
                self.mStdlog.write( "# warning: duplicate and inconsistent nid '%i'.\n" % nid)

            self.mMapNid2Length[nid] =  length

        infile.close()

        if self.mLogLevel >= 1:
            self.mStdlog.write("# read lenghts (%i nids) in %i seconds.\n" % ( len(self.mMapNid2Length),
                                                                               time.time() - t1))
            self.mStdlog.flush()

    #--------------------------------------------------------------------------------------
    def readCrossReferences( self, dbs = None):
        """read cross references from file for name mapping."""
        
        if self.mLogLevel >= 1:
            self.mStdlog.write("# reading references from %s\n" % self.mFilenameReferences) 
            self.mStdlog.flush()
            
        t1 = time.time()

        infile = self.openInputFile( self.mFilenameReferences )
        
        self.mMapPid2Nid = {}
        self.mMapAcc2Nid = {}

        ninput, ninconsistencies, nerrors = 0, 0, 0

        for line in infile:

            if line[0] == "#": continue

            ninput += 1
            try:
                nid, db, acc, pid = line[:-1].split("\t")[:4]
            except ValueError:
                self.mStderr.write( "parsing error in file %s: line %s\n" % ( self.mFilenameReferences, line[:-1] ) )
                nerrors += 1
                continue
            
            nid = int(nid)

            if dbs and int(db) not in dbs: continue

            if pid in self.mMapPid2Nid:
                x = self.mMapPid2Nid[pid]
                if nid != x:
                    if self.mLogLevel >= 1:
                        self.mStdlog.write( "# warning: duplicate and inconsistent id '%s' in cross-references.\n" % pid)
                    ninconsistencies += 1
                    if type(x) != ListType:
                        self.mMapPid2Nid[pid] = [x,nid]
                    else:
                        x.append( nid )
            else:
                self.mMapPid2Nid[pid] = nid

            if acc in self.mMapAcc2Nid:
                x = self.mMapAcc2Nid[acc]
                if nid != x:
                    if self.mLogLevel >= 1:
                        self.mStdlog.write( "# warning: duplicate and inconsistent acc '%s' in cross-references.\n" % acc)
                    ninconsistencies += 1
                    if type(x) != ListType:
                        self.mMapAcc2Nid[acc] = [x, nid]
                    else:
                        x.append( nid )
            else:
                self.mMapAcc2Nid[acc] = nid

        infile.close()

        if self.mLogLevel >= 1:
            self.mStdlog.write("# references: input=%i, ninconsistencies=%i, nerrors=%i.\n" % ( ninput,
                                                                                                ninconsistencies,
                                                                                                nerrors ))

            self.mStdlog.write("# read references (%i ids, %i acc) in %i seconds.\n" % ( len(self.mMapPid2Nid),
                                                                                         len(self.mMapAcc2Nid),
                                                                                         time.time() - t1))
            self.mStdlog.flush()

    #--------------------------------------------------------------------------------------        
    def mapPid2Nids( self, pid ):
        """map pid to references."""
        nids = []

        if not self.mMapPid2Nid: return nids

        if type(pid) is TupleType:
            pids = pid
        else:
            pids = [pid,]
                
        for p in pids:
            if p in self.mMapPid2Nid:
                nn = self.mMapPid2Nid[p]
                if type(nn) is ListType:
                    nids += nn
                else:
                    nids.append( nn )

        return list(set(nids))
        
    #--------------------------------------------------------------------------------------        
    def getCountsPerFamily( self, filename_input, is_nr = False ):
        """build family counts from a domains table.
        """
        
        infile = self.openInputFile( filename_input )
        
        counts_per_family = {}
        
        for line in infile:
            
            if is_nr:
                nid, start, end, family = line[:-1].split("\t")[:4]
            else:
                nid, start, end, ali, domain_id, domain_start, domain_end, domain_ali, family = line[:-1].split("\t")[:9]

            if family not in counts_per_family:
                counts_per_family[family] = Counts()
            
            c = counts_per_family[family]
            c.mNUnits += 1
            c.mNSequences.add( nid )
            c.mNResidues += int(end) - int(start) + 1

        for key, c in counts_per_family.items():
            c.mNSequences = len(c.mNSequences)
            
        infile.close()

        return counts_per_family

    #--------------------------------------------------------------------------------------
    def openOutputFile( self, filename ):
        """open a file for output. If it already exists, None is returned.
        """
        if os.path.exists(filename):
            if self.mLogLevel >= 1:
                self.mStdlog.write("# file %s already exists - step will be skipped.\n" % filename )
            return None
        else:
            if filename[-3:] == ".gz":
                return gzip.open( filename, "w" )
            else:
                return open( filename, "w" )

    #--------------------------------------------------------------------------------------
    def openInputFile( self, filename ):
        """open a file for input. Checks if the file exists. If it is compressed,
        open with gzip.
        """
        if not os.path.exists(filename):
            if self.mLogLevel >= 1:
                self.mStdlog.write("# file %s does not exist!" )
            raise "file %s not found" % filename
        else:
            if filename[-3:] == ".gz":
                return gzip.open(filename, "r" )
            else:
                return open(filename, "r" )
            
    #--------------------------------------------------------------------------------------
    def propagate( self ):
        """propagate annotations up or down the hierarchy."""

        params_100 = { 'level' : "", 'suffix' : "", 'name' : self.mName }
        params_100_up = { 'level' : "100", 'suffix' : "", 'name' : self.mName }
        params_90 = { 'level' : "90", 'suffix' : "", 'name' : self.mName }
        params_90_nr = { 'level' : "90", 'suffix' : "_nr", 'name' : self.mName }
        params_40 = { 'level' : "40", 'suffix' : "", 'name' : self.mName }
        params_40_nr = { 'level' : "40", 'suffix' : "_nr", 'name' : self.mName }

        ##########################################################
        ##########################################################
        ##########################################################
        ## propagate to 90 level
        if self.mInputFilenameMap100x90:
            self.propagateDown( self.mOutputFilenameAnnotations % params_100,
                                self.mOutputFilenameDescriptions % params_100,
                                self.mInputFilenameMap100x90,
                                self.mOutputFilenameAnnotations % params_90,
                                self.mOutputFilenameDescriptions % params_90)

            if self.mOptions.make_non_redundant:
                self.makeNonRedundant( 
                    self.mOutputFilenameAnnotations % params_90,
                    self.mOutputFilenameDescriptions % params_90,
                    self.mOutputFilenameAnnotations % params_90_nr,
                    self.mOutputFilenameDescriptions % params_90_nr )


        if self.mInputFilenameMap100x40:
            self.propagateDown( self.mOutputFilenameAnnotations % params_100,
                                self.mOutputFilenameDescriptions % params_100,
                                self.mInputFilenameMap100x40,
                                self.mOutputFilenameAnnotations % params_40,
                                self.mOutputFilenameDescriptions % params_40 )

            ## make non redundant
            if self.mOptions.make_non_redundant:
                self.makeNonRedundant( 
                    self.mOutputFilenameAnnotations % params_40,
                    self.mOutputFilenameDescriptions % params_40,
                    self.mOutputFilenameAnnotations % params_40_nr,
                    self.mOutputFilenameDescriptions % params_40_nr )

                ## propagate to level nrdb100
                if self.mPropagateUp:
                    self.propagateUp( self.mOutputFilenameAnnotations % params_40_nr,
                                      self.mOutputFilenameDescriptions % params_40_nr,
                                      self.mInputFilenameMap100x40,
                                      self.mOutputFilenameAnnotations % params_100_up,
                                      self.mOutputFilenameDescriptions % params_100_up )

    #--------------------------------------------------------------------------------------
    def propagateDown( self, 
                   input_filename_annotations, 
                   input_filename_descriptions,
                   input_filename_map,
                   output_filename_annotations,
                   output_filename_descriptions
                   ):
        """propagate annotations using a map.."""
        
        if self.mLogLevel >= 1:
            self.mStdlog.write( "# propagating %s to %s with %s\n" %\
                                    (input_filename_annotations,
                                     output_filename_annotations,
                                     input_filename_map ))
            sys.stdout.flush()
            

        if input_filename_map[-3:] == ".gz":
            raise "can not used compressed files - need random access for map."

        infile_map = self.openInputFile( input_filename_map )
        
        map_mem2rep, map_pair2pos = map_alignments.readMap( infile_map,
                                                            "pairsdb",
                                                            self.mOptions )
        
        ninput, noutput, nerrors, nmapped, nskipped, nsplits = 0, 0, 0, 0, 0, 0

        f = alignlib.AlignmentFormatEmissions()
        
        infile = self.openInputFile( input_filename_annotations )
        outfile = self.openOutputFile( output_filename_annotations )

        if outfile:

            for line in infile:
                if line[0] == "#": continue

                data = line[:-1].split("\t")
                
                if self.mLogLevel >= 4:
                    self.mStdlog.write( "# processing line: %s\n" % line[:-1] )
                    self.mStdlog.flush()

                mem_nid, start, end, ali, domain_id, domain_start, domain_end, domain_ali = data[:8]

                start, end = int(start), int(end)
                domain_start, domain_end = int(domain_start), int(domain_end)
                ninput += 1

                if self.mLogLevel >= 1 and (ninput % self.mReportStep) == 0:
                    self.mStdlog.write("# processing: iteration=%6i, mem_nid=%s, domain_id=%s\n" % (ninput, mem_nid, domain_id))
                    self.mStdlog.flush()

                if mem_nid not in map_mem2rep:
                    if self.mLogLevel >= 4:
                        self.mStdlog.write( "# no mapping necessary - directly output.\n" )
                        self.mStdlog.flush()
                    outfile.write( line )
                    noutput += 1
                    continue

                rep_nid = map_mem2rep[mem_nid]

                if end - start <= 0:
                    if self.mLogLevel >= 1:
                        self.mStdlog.write("# invalid range in annotation: mapping from %s failed for %s\n" % (rep_nid, line[:-1]))
                    nerrors += 1
                    continue
                
                f.mRowFrom, f.mRowTo = start - 1, end
                f.mColFrom, f.mColTo = domain_start - 1, domain_end
                f.mRowAlignment = ali
                f.mColAlignment = domain_ali
                ali_mem2domain = alignlib.makeAlignmentVector()
                f.copy( ali_mem2domain )

                ali_mem2rep = map_alignments.getMap( mem_nid, 
                                                     rep_nid, 
                                                     infile_map, 
                                                     map_pair2pos,
                                                     "pairsdb" )
                
                if ali_mem2rep == None:
                    if self.mLogLevel >= 1:
                        self.mStdlog.write( "# could not find map from %s to %s\n" % (mem_nid, rep_nid ))
                    nerrors += 1
                    continue

                if self.mLogLevel >= 4:
                    self.mStdlog.write( "# using  map: %s\n" % str(f) )
                    self.mStdlog.write( "# applying map: %s\n" % str(alignlib.AlignmentFormatEmissions( ali_mem2rep) ))
                    self.mStdlog.flush()

                ali_rep2domain = alignlib.makeAlignmentVector()
                alignlib.combineAlignment( ali_rep2domain, 
                                           ali_mem2rep,
                                           ali_mem2domain,
                                           alignlib.RR )

                if ali_rep2domain.isEmpty():
                    nerrors += 1
                    if self.mLogLevel >= 1:
                        self.mStdlog.write("# mapping to rep %s failed for %s\n" % (rep_nid, line[:-1]))
                    continue

                noutput += 1
                nmapped += 1

                if self.mSplitDiscontinuousDomains:
                    alis = alignlib.splitAlignment( ali_rep2domain, self.mSplitAtGaps )
                    if len(alis) > 1:
                        nsplits += 1
                else:
                    alis = (ali_rep2domain,)
                        
                for ali in alis:
                    f.fill( ali )
                    outfile.write( "\t".join( map(str, (
                                                  rep_nid, 
                                                  f.mRowFrom + 1,
                                                  f.mRowTo,
                                                  f.mRowAlignment,
                                                  domain_id,
                                                  f.mColFrom + 1,
                                                  f.mColTo,
                                                  f.mColAlignment)) + data[8:] ) + "\n" )
            outfile.close()

            if self.mLogLevel >= 1:
                self.mStdlog.write("# propagating domains: ninput=%i, noutput=%i, nmapped=%i, nskipped=%i, nerrors=%i, nsplits=%i\n" % (ninput, noutput, nmapped, nskipped, nerrors, nsplits) )
                

        infile.close()
        
        self.createDescriptions( output_filename_annotations,
                                 output_filename_descriptions,
                                 input_filename_descriptions )

    #--------------------------------------------------------------------------------------
    def propagateUp( self, 
                     input_filename_annotations, 
                     input_filename_descriptions,
                     input_filename_map,
                     output_filename_annotations,
                     output_filename_descriptions
                     ):
        """propagate annotations using a map."""
        
        if self.mLogLevel >= 1:
            self.mStdlog.write( "# propagating %s to %s with %s\n" %\
                                    (input_filename_annotations,
                                     output_filename_annotations,
                                     input_filename_map ))
            sys.stdout.flush()

        if input_filename_map[-3:] == ".gz":
            raise "can not used compressed files - need random access for map."

        infile_map = self.openInputFile( input_filename_map )
        
        map_rep2mem, map_pair2pos = map_alignments.readMap( infile_map,
                                                            "pairsdb",
                                                            self.mOptions, 
                                                            invert = True )
        
        ninput, noutput, nerrors, nmapped, nskipped, nsplits = 0, 0, 0, 0, 0, 0

        f = alignlib.AlignmentFormatEmissions()
        
        infile = self.openInputFile( input_filename_annotations )
        outfile = self.openOutputFile( output_filename_annotations )

        if outfile:

            for line in infile:
                if line[0] == "#": continue

                data = line[:-1].split("\t")
                
                if self.mLogLevel >= 4:
                    self.mStdlog.write( "# processing line: %s\n" % line[:-1] )
                    self.mStdlog.flush()

                rep_nid, start, end, family = data[:4]
                info = data[4:]

                start, end = int(start), int(end)
                ninput += 1

                if end - start <= 0:
                    if self.mLogLevel >= 1:
                        self.mStdlog.write("# invalid range in annotation: mapping from %s failed for %s\n" % (rep_nid, line[:-1]))
                    nerrors += 1
                    continue
                
                if self.mLogLevel >= 1 and (ninput % self.mReportStep) == 0:
                    self.mStdlog.write("# processing: iteration=%6i, mem_nid=%s, domain_id=%s\n" % (ninput, rep_nid, family))
                    self.mStdlog.flush()
                
                ## output the representative itself
                outfile.write( line )
                noutput += 1
                
                if rep_nid not in map_rep2mem: continue

                for mem_nid in map_rep2mem[rep_nid]:

                    ali_rep2domain = alignlib.makeAlignmentVector()
                    alignlib.addDiagonal2Alignment( ali_rep2domain, start, end, -start )

                    ali_mem2rep = map_alignments.getMap( rep_nid, 
                                                         mem_nid, 
                                                         infile_map, 
                                                         map_pair2pos,
                                                         "pairsdb" )

                    if ali_mem2rep == None:
                        if self.mLogLevel >= 1:
                            self.mStdlog.write( "# could not find map from %s to %s\n" % (rep_nid, mem_nid ))
                        nerrors += 1
                        continue

                    if self.mLogLevel >= 4:
                        self.mStdlog.write( "# applying map: %s\n" % str(alignlib.AlignmentFormatEmissions( ali_mem2rep) ))
                        self.mStdlog.flush()

                    ali_mem2domain = alignlib.makeAlignmentVector()
                    alignlib.combineAlignment( ali_mem2domain, 
                                               ali_mem2rep,
                                               ali_rep2domain,
                                               alignlib.CR )

                    if ali_rep2domain.isEmpty():
                        nerrors += 1
                        if self.mLogLevel >= 1:
                            self.mStdlog.write("# mapping to rep %s failed for %s\n" % (rep_nid, line[:-1]))
                        continue

                    noutput += 1
                    nmapped += 1

                    if self.mSplitDiscontinuousDomains:
                        alis = alignlib.splitAlignment( ali_mem2domain, self.mSplitAtGaps )
                        if len(alis) > 1:
                            nsplits += 1
                    else:
                        alis = (ali_mem2domain,)

                    for ali in alis:
                        f.fill( ali )
                        outfile.write( "\t".join( map(str, (
                                                      mem_nid, 
                                                      f.mRowFrom + 1,
                                                      f.mRowTo,
                                                      family ) ) ) +\
                                           "\t".join( info ) +\
                                           "\n" )
            outfile.close()

            if self.mLogLevel >= 1:
                self.mStdlog.write("# propagating domains: ninput=%i, noutput=%i, nmapped=%i, nskipped=%i, nerrors=%i, nsplits=%i\n" % (ninput, noutput, nmapped, nskipped, nerrors, nsplits) )
                

        infile.close()
        
        self.createDescriptions( output_filename_annotations,
                                 output_filename_descriptions,
                                 input_filename_descriptions,
                                 is_nr = True)

    #--------------------------------------------------------------------------------------
    def getTemporaryFilename( self, fn ):
        """return a temporary filename."""
        if not self.mTempdir:
            self.mTempdir = tempfile.mkdtemp( dir = self.mTempdirLocation )
            
        return self.mTempdir + "/" + fn 

    #--------------------------------------------------------------------------------------
    def makeNonRedundant( self, 
                          input_filename_annotations, 
                          input_filename_descriptions,
                          output_filename_annotations,
                          output_filename_descriptions
                          ):
        """make a non-redundant clone of a database. 

        Overlapping annotations of the same type are merged. 
        The non-redundant clone loses information of the individual 
        domains.
        """
        
        if self.mLogLevel >= 1:
            self.mStdlog.write( "# making non-redundant clone of %s to %s\n" %\
                                    (input_filename_annotations,
                                     output_filename_annotations ))
            sys.stdout.flush()

        ninput, noutput, nerrors, nmapped, nskipped = 0, 0, 0, 0, 0

        infile = self.openInputFile( input_filename_annotations )
        outfile = self.openOutputFile( output_filename_annotations )

        if outfile:
            tempfile_name = self.getTemporaryFilename( "make_non_redundant_%s_%s" % (input_filename_annotations, output_filename_annotations) )
            
            tmpfile = open(tempfile_name, "w")

            if self.mLogLevel >= 1:
                self.mStdlog.write("# creating temporary data in %s\n" % tempfile_name )

            for line in infile:
                if line[0] == "#": continue

                ninput += 1

                data = line[:-1].split("\t")

                nid, start, end, ali, domain_id, domain_start, domain_end, domain_ali, family = data[:9]
                tmpfile.write( "\t".join( (nid, family, domain_start, domain_end ) ) + "\n" )
                
            tmpfile.close()

            infile.close()
            
            if self.mLogLevel >= 1:
                self.mStdlog.write("# starting sorting of %i entries\n" % ninput )

            cmd = "sort -T%s -k1,1n -k2,2 -k3,3n %s" % (self.mTempdir, tempfile_name)

            (last_nid, last_family, first_from, last_info) = (None, None, None, None )
            last_info = ()

            last_to  = 0
            
            bufsize = 1024

            for line in subprocess.Popen( cmd, shell=True, bufsize=bufsize, stdout=subprocess.PIPE).stdout:
                
                nid, family, start, end = line[:-1].split("\t")

                nid, start, end = map( int, (nid, start, end))

                ## skip small domains
                if end - start < self.mNRMinAnnotationLength:
                    continue

                ## new entry if different nid
                if last_nid != nid:
                    if last_nid:
                        outfile.write( string.join( map(str,
                                                        (last_nid, first_from, last_to, last_family) + last_info),
                                                    "\t") + "\n" )
                        noutput += 1
                    last_to = 0
                    first_from = start

                ## new entry if different family
                elif last_family != family:

                    # check if domain is contained in current domain. If so, skip this one.
                    # if start >= first_from and end <= last_to:
                    #    continue

                    # write current domain as there is a new assignment
                    outfile.write( string.join( map(str,
                                                    (last_nid, first_from, last_to, last_family) + last_info),
                                                "\t") + "\n" )
                    noutput += 1
                    if self.mNRShortenDomains:
                        if start <= last_to:
                            print "--> shortened domain in %i due to %s (%i-%i)" % (last_nid,
                                                                                    str(last_family),
                                                                                    first_from,
                                                                                    last_to), r                
                            first_from = last_to + 1
                        else:
                            first_from = start
                    else:
                        first_from = start

                    last_to = 0

                ## new entry if no overlap
                elif last_to - self.mNRMinOverlap < start:
                    # write current domain, as there is a gap in between two domains
                    outfile.write( string.join( map(str,
                                                    (last_nid, first_from, last_to, last_family) + last_info),
                                                "\t") + "\n" )
                    noutput += 1
                    first_from = start
                    last_to = 0

                last_nid = nid
                last_family = family
                # last_info = info
                last_to = max( last_to, end )
                
            if last_nid:
                outfile.write( string.join( map(str,
                                                (last_nid, first_from, last_to, last_family) + last_info),
                                            "\t") + "\n" )
                noutput += 1

            if self.mLogLevel >= 1:
                self.mStdlog.write("# removing redundant domains: ninput=%i, noutput=%i, nmapped=%i, nskipped=%i, nerrors=%i\n" % (ninput, noutput, nmapped, nskipped, nerrors) )

            outfile.close()                

        self.createDescriptions( 
            output_filename_annotations,
            output_filename_descriptions,
            input_filename_descriptions,
            is_nr = True )

    #--------------------------------------------------------------------------------------
    def createDescriptions( self, 
                            input_filename_annotations,
                            output_filename_descriptions,
                            input_filename_descriptions = None,
                            is_nr = False ):

        """create descriptions from an annotations file and an
        optional descriptions file.

        The optional descriptions file is used to map extra columns.
        """

        if self.mLogLevel >= 1:
            self.mStdlog.write("# creating descriptions for %s.\n" % input_filename_annotations )
            self.mStdlog.write("# descriptions from %s to %s.\n" % (input_filename_descriptions, output_filename_descriptions ))

        outfile = self.openOutputFile( output_filename_descriptions )
        if outfile:

            tstart = time.time()
            counts = self.getCountsPerFamily( input_filename_annotations, is_nr = is_nr )

            if input_filename_descriptions:
                map_family2info = {}
                infile = self.openInputFile( input_filename_descriptions )
                for line in infile:
                    if line[0] == "#": continue
                    data = line[:-1].split("\t")
                    if len(data) <= 5:
                        ## no extra info, skip
                        map_family2info = None
                        break

                    map_family2info[data[0]] = data[5:]
            else:
                map_family2info = None
                
            if map_family2info:
                families = map_family2info.keys()
            else:
                families = counts.keys()

            ninput, noutput, nskipped, nempty = 0, 0, 0, 0
            info = ""

            for id in families:
                
                ninput += 1

                if map_family2info:
                    if id not in map_family2info:
                        if self.mLogLevel >= 1:
                            self.mStdlog.write("# no info for family %s - skipped\n" % id)
                        nskipped += 1
                        continue
                    info = "\t" + "\t".join( map_family2info[id] )
                    
                if id not in counts:
                    nempty += 1
                        
                    if self.mKeepEmptyFamilies:
                        outfile.write( "\t".join( map(str, (id, 
                                                            0,
                                                            0,
                                                            0,
                                                            0))) +\
                                           "%s\n" % info )
                        noutput += 1
                    else:
                        if self.mLogLevel >= 1:
                            self.mStdlog.write("# empty family %s - skipped\n" % id)
                        nskipped += 1
                        continue
                else:
                    c = counts[id]
                    outfile.write( "\t".join( map(str, (id, 
                                                        c.mNUnits,
                                                        c.mNSequences,
                                                        c.mNResidues,
                                                        c.mNResidues / c.mNUnits ))) +\
                                       "%s\n" % info )
                    noutput += 1

            outfile.close()

            if self.mLogLevel >= 1:
                self.mStdlog.write("# created descriptions: ninput=%i, noutput=%i, nskipped=%i, nempty=%i, time=%i\n" % (ninput, noutput, nskipped, nempty, time.time() - tstart ) )

    #--------------------------------------------------------------------------------------
    #--------------------------------------------------------------------------------------
    def __findNids( self, 
                    pid = None, 
                    sequence = None, 
                    perform_blast = False ):
        """find an nid for a pid or a sequence. It returns
        the first pid found.
        """

        if not pid and not sequence:
            raise "either pid or sequence have to be specified."
        
        nids = []
        # 1. lookup identifier in cross_references
        if pid and self.mMapPid2Nid:
            nids = self.mapPid2Nids( pid )

        if nids: return nids, "PID"
        if not sequence: return None, "MISSED"

        sequence = string.replace( sequence, "\n", "")
        sequence = string.upper( sequence ) 
        first_res = 0

        # 2. if unsuccessful, lookup hid in nrdb
        if self.mMapHid2Nid:
            hid = Tools.calculateHID( sequence )
            if hid in self.mMapHid2Nid:
                return (self.mMapHid2Nid[hid],), "HID"

        # 3. if unsuccessfull, try to find it using suffix array
        if self.mSary and len(sequence) > 10:
            matches = self.mSary.search( sequence )
            if matches:
                return map( int, matches), "SARY"

        # 4. if unsuccessful, do BLAST-search against nrdb100
        if perform_blast:
            b = WrapperBlast()
            b.Set_Database("rsdb100")
            ## perform search
            # 1. no masking: -F F
            # 2. no multiple hits: -P 1
            # 3. speed up search: -W 3 -f 1000
            result = b.Search_Sequence( sequence, " -F F -P 1 -W 3 -f 1000 -v 10 -b 10" )     

            if result:
                x = result[0]                                           # just look at first entry
                if x[:13] == 'No-hits-found':
                    return 0, "MISSED"

                (values) = string.split(x, "\t")
                pide = values[2]            
                nid = string.atoi(string.split(values[0], "|")[0])

                if string.atoi(pide) >= min_pide:                              # allow up to 5% mismatches
                    return nid, "BLAST"

        return 0, "MISSED"

    #--------------------------------------------------------------------------------------
    def findNids( self, pid = None, sequence = None ):
        """find an nid given a protein identifier and/or a sequence. 
        """

        if self.mLogLevel >= 1:
            if pid:
                self.mStdlog.write( "# --> processing pid: %s\n" % pid )
                self.mStdlog.flush()

        if self.mIgnorePid: pid = None
            
        if pid:
            if not self.mFound.has_key(pid):
                nids, method = self.__findNids( pid = pid,
                                                sequence = sequence,
                                                perform_blast = self.mPerformBlast )
                self.mFound[pid] = nids            
            else:
                nids = self.mFound[pid]
                if nids == None:
                    method = "MISSED LOOKUP"
                else:
                    method = "LOOKUP"
        else:
            nids, method = self.__findNids( pid, 
                                            sequence,
                                            perform_blast = self.mPerformBlast )
            
        if self.mLogLevel >= 2:
            if nids:
                self.mStdlog.write("# found: %s -> %s, method=%s\n" % (pid, str(nids), method ))
            else:
                self.mStdlog.write("# not found: %s\n" % (pid))
                                 
        if method == "PID":
            self.mFromPid += 1
        elif method == "HID":
            self.mFromHid += 1
        elif method == "BLAST":
            self.mFromBlast += 1
        elif method == "LOOKUP":
            self.mFromLookup += 1
        elif method == "SARY":
            self.mFromSary += 1
        else:
            self.mNotFound += 1

        return nids

    #--------------------------------------------------------------------------------------        
    def outputAnnotationWithAlignment( self, outfile, nid, sequence, domain_sequence ):
        """align domain_sequence to sequence and output.

        return True, if successfully written.
        """

        s1 = alignlib.makeSequence( domain_sequence )
        s2 = alignlib.makeSequence( sequence )

        alignator    = alignlib.makeAlignatorGroupies()
        map_domain2nrdb = alignlib.makeAlignmentVector()
        alignator.align( map_domain2nrdb, s1, s2 )

        if map_domain2nrdb.isEmpty():
            if self.mLogLevel >= 1:
                self.mStdlog.write("# empty alignment for nid %s\n" % nid )
            self.mFailed += 1
            return False

        domain_from  = map_domain2nrdb.getRowFrom() + 1
        domain_to = map_domain2nrdb.getRowTo()
        nrdb_from = map_domain2nrdb.getColFrom() + 1
        nrdb_to = map_domain2nrdb.getColTo()

        coverage = float(domain_to - domain_from + 1) / float(len(domain_sequence))

        if self.mLogLevel >= 3:
            self.mStdlog.write("# processing nid %s\n" % str(nid) )
            self.mStdlog.write("# --> coverage of domain= %5.2f, domain length=%i, sequence length=%i\n" %\
                               (coverage, len(domain_sequence), len(sequence) ))

        if coverage < self.mMinCoverage:
            if self.mLogLevel >= 1:
                self.mStdlog.write( "# -----> coverage low for nid=%s: from=%i, to=%i, coverage of domain=%5.2f, domain length=%i, sequence length=%i\n" %\
                                    ( str(nid), domain_from, domain_to, coverage, len(domain_sequence), len(sequence) ) )
            self.mFailed += 1
            return False

        if self.mLogLevel >= 4:
            self.mStdlog.write( "# --> map_domain2nrdb: %s\n"  % str(alignlib.AlignmentFormatEmissions( map_domain2nrdb ) ))
            sys.stdout.flush()

        f = alignlib.AlignmentFormatEmissions( map_domain2nrdb )

        return self.outputAnnotation( outfile, nid,
                                      nrdb_from, nrdb_to, f.mColAlignment,
                                      domain_from, domain_to, f.mRowAlignment )

    
    #--------------------------------------------------------------------------------------        
    def processFastaEntry( self, outfile,
                           description_line, domain_sequence = None, nids = None):
        """process a fasta entry.

        If nids are given, no new lookup is performed. This allows this routine to be
        called with a sequence of fragments.
        """

        # try different methods to locate a nid for a given domain_sequence with domain-identifier
        if domain_sequence:
            domain_sequence = re.sub("\s", "", domain_sequence)
            domain_sequence = string.upper(domain_sequence)
            domain_sequence = re.sub("[^A-Z]","", domain_sequence)

        
        ## find pids using sequence
        if not nids:
            # parse description line
            if not self.parseDescriptionLine( description_line ):
                return False

            nids = self.findPid( pid = self.mPid, sequence = domain_sequence )

            if not nids: return False
            
        nsuccess = 0
        for nid in nids:
            ## if domain_sequence is given, align sequence to it
            if domain_sequence and self.mPerformAlignment and self.mFasta:

                sequence = self.mFasta.getSequence( str(nid) )

                if self.mLogLevel >= 5:
                    self.mStdlog.write("# domain=%s\n# nrdb=%s\n" % (domain_sequence, sequence))

                success = self.outputAnnotationWithAlignment( outfile, nid, sequence, domain_sequence )
            else:
                success = self.outputAnnotation( outfile, nid, 0, 0, "", 0, 0, "" )
                
            if success: nsuccess += 1

        return nsuccess > 0
    #--------------------------------------------------------------------------------------
    def outputAnnotation( self, outfile, nid, nrdb_from, nrdb_to, nrdb_ali, domain_from, domain_to, domain_ali ):
        """write a line to outputfile for loading. 

        This default implemenation just writes
        mDomainId and mDomainClass without coordinates.

        return True, if successfully written.
        """
    
        outfile.write( string.join ( (
                    str(nid),
                    str(nrdb_from),
                    str(nrdb_to),
                    nrdb_ali,
                    self.mDomainId,
                    str(domain_from),
                    str(domain_to),
                    domain_ali,
                    self.mDomainClass), "\t") + "\n" )

        outfile.flush()

        return True
    #--------------------------------------------------------------------------------------    
    def createAnnotationsFromFasta( self, input_filename, output_filename_annotations ):
        """parse a fasta file and create annotations from it. 

        This method calls parseDescriptionLine for every description line 
        and outputAnnotation.
        """

        if self.mLogLevel >= 1:
            self.mStdlog.write("# starting parsing fasta file %s\n" % input_filename )

        outfile = self.openOutputFile( output_filename_annotations )            

        if outfile:

            self.resetStats()

            infile = self.openInputFile( input_filename )
            
            iterator = FastaIterator.FastaIterator( infile )

            ninput, noutput, nskipped = 0, 0, 0            
            tstart = time.time()

            while 1:
                cur_record = iterator.next()
                if cur_record is None: break
                ninput += 1
                if self.processFastaEntry( outfile, cur_record.title, cur_record.sequence ):
                    noutput += 1
                else:
                    nskipped += 1

            if self.mLogLevel >= 1:
                self.printStats()
                self.mStdlog.write("# parsing from fasta finished: ninput=%i, noutput=%i, nskipped=%i, time=%i\n" % (ninput, noutput, nskipped, time.time() - tstart ) )

            infile.close()
            outfile.close()

    #--------------------------------------------------------------------------------------
    def checkIfOld( self ):
        """check, if entry exists. This default implementation requires mDomainId and mDomainClass
        members to be set
        """
        
        if self.mTableAnnotations.HasEntry( self.mDomainId, self.mDomainClass ):
            return 1
        else:
            # if just scop_class has changed, all entries for this domain have to be deleted.
            self.mTableAnnotations.RemoveEntry( self.mDomainId )
            return 0
        
