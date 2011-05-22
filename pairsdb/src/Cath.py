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
"""

from NrdbAnnotations import NrdbAnnotations
from Scop import Scop

class Cath (Scop):

    mName = "cath"

    def __init__ (self, options = None):

	Scop.__init__( self, options )

        # parser for classification id
        self.mClassParser = re.compile( "^\S+\s+(\S+)")          
        self.mClassIdentifierLength = 4

    def parseIdFromLine( self, line ):
        """parse family identifier and description from a scop hierarchy file."""

        data = re.split("\s+", line[:-1] )

        level, rep = data[0], data[1]
        description = string.join(data[2:], " ")[1:]
        ids = string.split(level, ".")

        id = string.join( map( lambda x: "0" * (self.mClassIdentifierLength - len(x)) + x, ids),"")

        return id, description

#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id$", usage=USAGE )

    parser.add_option( "--input-filename-annotations", dest="input_filename_annotations", type="string" ,
                       help="INPUT filename of cath domain sequences in fasta format.")

    parser.add_option( "--input-filename-descriptions", dest="input_filename_descriptions", type="string" ,
                       help="INPUT filename of with CATH annotations: hierarchy.")

    parser.set_defaults( 
        input_filename_annotations = "cath.fasta",
        input_filename_descriptions = "cath.hierarchy",
        )

    Cath().addOptions( parser )

    (options, args) = Experiment.Start( parser )

    cath = Cath( options )

    cath.run()
    
    Experiment.Stop()
        



        




