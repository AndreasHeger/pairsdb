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

USAGE="""python check_requirements.py [OPTIONS]

check if required binaries exist for python.
"""

import os, sys, re, string, optparse
import Experiment

def checkBinary( binary, paths ):
    
    print "checking %s..." % binary,

    for path in paths:
        if os.path.exists( path + "/" + binary ):
            print "success"
            return True

    return False

if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id$", usage=USAGE )

    parser.set_defaults()

    (options, args) = Experiment.Start( parser, add_pipe_options = True )

    binaries = ("blastall", "coils", "cd-hit", "decodeanhmm", "sary", "mksary", "formatdb" )

    paths = os.environ['PATH'].split(":")

    nerrors = 0

    for binary in binaries:
        nerrors += checkBinary( binary, paths )

    print "importing alignlib...",

    try:
        import alignlib
        x = alignlib.AlignmentFormatExplicit()
        print "success"
    except ImportError:
        nerrors += 1
        print "failed"

    Experiment.Stop()
    
    if nerrors == 0:
        sys.exit(0)
    else:
        sys.exit(1)
