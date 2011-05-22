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

PDB annotations are taken from the cross-references table. 
"""

from NrdbAnnotations import NrdbAnnotations

class Pdb (NrdbAnnotations):

    mName = "pdb"

    def __init__ (self, options = None):

	NrdbAnnotations.__init__( self, options )
        
    #--------------------------------------------------------------------------------------
    def run( self ):
        
        ## read references for pdb only
        self.readCrossReferences( dbs= (3,) )
        
        ## read length map
        self.readLengths()

        ## create annotations
        self.createAnnotations()
        
        ## create descriptions
        self.createDescriptions( self.mOutputFilenameAnnotations % self.mParams100,
                                 self.mOutputFilenameDescriptions % self.mParams100 )
        
        self.propagate()

    #--------------------------------------------------------------------------------------        
    def createAnnotations( self ):
        """create annotations.
        """
        
        if self.mLogLevel >= 1:
            self.mStdlog.write("# parsing annotations - output goes to %s.\n" % (self.mOutputFilenameAnnotations % self.mParams100))
        
        outfile = self.openOutputFile( self.mOutputFilenameAnnotations % self.mParams100 )

        self.resetStats()
        
        if outfile:

            tstart = time.time()
            ninput, noutput, nskipped = 0, 0, 0

            for acc, nid in self.mMapAcc2Nid.items():

                ninput += 1
                
                if nid not in self.mMapNid2Length:
                    if self.mLogLevel >= 1:
                        self.mStdlog.write("# can not find length for nid %i (acc=%s)\n" % (nid, acc))
                    nskipped += 1
                    continue
                
                length = self.mMapNid2Length[nid]
                self.mFromPid += 1
                outfile.write( "\t".join( map(str, \
                                                  (nid, 
                                                   1, length, "+%i" % length,
                                                   acc,
                                                   1, length, "+%i" % length,
                                                   acc) ) ) + "\n" )
                noutput += 1
            
            outfile.close()

            if self.mLogLevel >= 1:
                self.printStats()
                self.mStdlog.write("# creating pdb: ninput=%i, noutput=%i, nskipped=%i, time=%i\n" % (ninput, noutput, nskipped, time.time() - tstart ) )
            
        
#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id$", usage=USAGE )

    parser.set_defaults( )

    Pdb().addOptions( parser )

    (options, args) = Experiment.Start( parser )

    pdb = Pdb( options )

    pdb.run()
    
    Experiment.Stop()
        



        




