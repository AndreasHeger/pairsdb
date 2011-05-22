####
####
##
## Project PairsDB
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: Nrdb_Interpro.py,v 1.1.1.1 2002/07/02 10:46:57 heger Exp $
##
##
####
####

#----------------------------------------------------------------
import sys, re, string, gzip, optparse, time

import Experiment

USAGE="""

SCOP annotations are taken from the cross-references table. 
"""

from NrdbAnnotations import NrdbAnnotations

class Scop (NrdbAnnotations):

    mName = "scop"

    def __init__ (self, options = None):

	NrdbAnnotations.__init__( self, options )

        if options:
            self.mFilenameAnnotations = options.input_filename_annotations
            self.mFilenameDescriptions = options.input_filename_descriptions

        # parser for ddd classification id
        self.mClassParser = re.compile( "^\S+\s+(\S\.[0-9.]+)")          
        self.mClassIdentifierLength = 3

        ## set flags for lookup
        self.mUseMultipleEntries = True
        self.mPerformAlignment = True
        self.mSplitDiscontinuousDomains = True
        self.mKeepEmptyFamilies = True
        self.mWithFragments = True
        self.mIgnorePid = False
        self.mPropagateUp = True

    #--------------------------------------------------------------------------------------
    def run( self ):
        
        ## read references for scop only
        self.readCrossReferences( dbs= (3,) )
        
        ## create annotations
        self.createAnnotationsFromFasta( self.mFilenameAnnotations,
                                         self.mOutputFilenameAnnotations % self.mParams100 )
        
        ## create descriptions
        self.addDescriptions( 
            self.mFilenameDescriptions,
            self.mOutputFilenameAnnotations % self.mParams100,
            self.mOutputFilenameDescriptions % self.mParams100 )
        
        ## propagate to lower levels
        self.propagate()
        
    #--------------------------------------------------------------------------------------
    def processFastaEntry( self, outfile, description_line, domain_sequence):
        """process a fasta entry.

        Discontinuous domains contain an X in the middle. This is
        the "genetic format". For more information, see

        http://astral.berkeley.edu/scopseq-os-1.73.html
        """

        if self.mWithFragments:
            fragments = string.split( domain_sequence, "X")
        else:
            fragments = [domain_sequence,]

        ## find nid based on longest fragment
        max_fragment = ""
        l = 0
        for fragment in fragments:
            if len(fragment) > l:
                max_fragment = fragment
                l = len(fragment)

        ## parse description line
        if not self.parseDescriptionLine( description_line ):
            if self.mLogLevel >= 2: 
                self.mStdlog.write("# error in parsing description: %s\n" % description_line )
            return False

        ## find pid using sequence
        nids = self.findNids( pid = self.mPid, sequence = max_fragment )

        if not nids: 
            if self.mLogLevel >= 2: 
                self.mStdlog.write("# no matches found for pid %s\n" % str(self.mPid) )
            return False

        self.mFragmentNo = 0

        n = 0
        for fragment in fragments:
            self.mFragmentNo += 1
            if NrdbAnnotations.processFastaEntry( self, outfile, description_line, fragment, nids):
                n += 1
        return n > 0

    #--------------------------------------------------------------------------------------        
    def parseDescriptionLine( self, line):
        """extract information from description line.
        Convert a scop-class from 1.1.1.1 to 001001001001.
        """

        id = string.upper(line[1:5])
        chain = string.upper(line[5:6])

        self.mDomainId  = line[1:7]

        if chain != "_":
            self.mPdbId = "%s-%s" % (id, chain)
        else:
            self.mPdbId = id
        
        try:
            class_id = self.mClassParser.search( line ).groups()[0]
        except AttributeError:
            raise "parsing error in line %s" % line
            
        ids = class_id.split(".")
        self.mFamily = string.join( map( lambda x: "0" * (self.mClassIdentifierLength - len(x)) + x, ids),"")

        self.mPid = self.mPdbId
        
        return True
        
    #--------------------------------------------------------------------------------------
    def outputAnnotation( self, outfile, nid, nrdb_from, nrdb_to, nrdb_ali, pdb_from, pdb_to, pdb_ali ):
        """write a line to outputfile for loading.
        the pdb-identifier is combined to pdb-chain, which is conform then to table
        cross_references. Note, that id and chain are upcasted, which introduces a
        a problem in structures with more than 26 chains (so far only 1fnt, whcih is not in scop)
        """

        outfile.write( string.join ( (
            str(nid),
            str(nrdb_from),
            str(nrdb_to),
            nrdb_ali,
            self.mDomainId,
            str(pdb_from),
            str(pdb_to),
            pdb_ali,
            self.mFamily,
            self.mPdbId,
            str(self.mFragmentNo)), "\t") + "\n" )

        return True

    def parseIdFromLine( self, line ):
        """parse family identifier and description from a scop hierarchy file."""

        (id, level, domain_class, scopid, description) = string.split(line[:-1], "\t")

        if level not in ("cl", "cf", "sf", "fa"):
            return None, None

        ids = string.split(domain_class, ".")
        id = string.join( map( lambda x: "0" * (self.mClassIdentifierLength - len(x)) + x, ids),"")

        return id, description

    #--------------------------------------------------------------------------------------
    def addDescriptions( self, 
                         input_filename_descriptions,
                         input_filename_annotations,
                         output_filename_descriptions,
                         ):

        """parser description from interpro-file.

        call this method after parseDescriptions().
        """

        if self.mLogLevel >= 1:
            self.mStdlog.write("# parsing descriptions.\n" )
            
        infile = self.openInputFile( input_filename_descriptions )
        outfile = self.openOutputFile( output_filename_descriptions )

        if outfile:

            counts = self.getCountsPerFamily( input_filename_annotations )

            ninput, noutput, nskipped, nempty = 0, 0, 0, 0

            for line in infile:

                if line[0] == "#": continue

                ninput += 1

                id, description = self.parseIdFromLine( line )
                
                if id == None: continue

                if id not in counts:
                    nempty += 1
                    if self.mKeepEmptyFamilies:
                        outfile.write( "\t".join( map(str, (id, 
                                                            0, 
                                                            0,
                                                            0,
                                                            0,
                                                            description,
                                                            ))) + "\n")
                        noutput += 1
                    else:
                        if self.mLogLevel >= 1:
                            self.mStdlog.write("# skipping empty family: %s \n" % id)
                        nskipped += 1
                        continue
                else:
                    outfile.write( "\t".join( map(str, (id, 
                                                        counts[id].mNUnits,
                                                        counts[id].mNSequences,
                                                        counts[id].mNResidues,
                                                        counts[id].mNResidues / counts[id].mNUnits,
                                                        description,
                                                        ))) + "\n")
                    noutput += 1

            outfile.close()
            infile.close

            if self.mLogLevel >= 1:
                self.mStdlog.write("# parsing descriptions: ninput=%i, noutput=%i, nskipped=%i, nempty=%i\n" % (ninput, noutput, nskipped, nempty) )

#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id$", usage=USAGE )

    parser.add_option( "--input-filename-annotations", dest="input_filename_annotations", type="string" ,
                       help="INPUT filename of for scop: ASTRAL domain sequences.")

    parser.add_option( "--input-filename-descriptions", dest="input_filename_descriptions", type="string" ,
                       help="INPUT filename of with SCOP annotations: hierarchy.")

    parser.set_defaults( 
        input_filename_annotations = "scop.fasta",
        input_filename_descriptions = "scop.hierarchy",
        )

    Scop().addOptions( parser )

    (options, args) = Experiment.Start( parser )

    scop = Scop( options )

    scop.run()
    
    Experiment.Stop()
        



        




