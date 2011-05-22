################################################################################
#
#   PairsDB
#
#   $Id$
#
#   Copyright (C) 2007 Andreas Heger 
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
import os, sys, string, re, optparse

USAGE="""python %s [OPTIONS] < in.fasta > out.fasta
""" % sys.argv[0]

import sys, re, optparse

import Experiment

class Masker:

    mMaskChar = "X"

    mMapMethod2Id = {
        'bias' : 1,
        'tmhmm' : 2,
        'coils' : 3,
        'short' : 4 }
        
    def __init__(self, infile, masks, options ):

        self.mErrors = 0
        self.mMasks = {}

        self.mStdlog = options.stdlog
        self.mLogLevel = options.loglevel
        
        take = set()
        for x in masks:
            if x == "all":
                take = set( self.mMapMethod2Id.values() )
                break
            else:
                take.add( self.mMapMethod2Id[x] )
                
        with_residue = set( [1,] )
                  
        for line in infile:
            if line[0] == "#": continue

            data = line[:-1].split("\t" )

            if len(data) == 5:
                name, first_res, last_res, method, residue = data
            elif len(data) == 4:
                name, first_res, last_res, method = data
                residue = None
            else:
                self.mStdlog.write("# parsing error in line: %s\n" % line[:-1] )
                self.mErrors += 1
                continue
            
            try:
                method = int(method)
            except ValueError:
                self.mStdlog.write("# parsing error in line: %s\n" % line[:-1] )
                self.mErrors += 1
                continue
            
            if method not in take: continue
            if method not in with_residue:
                residue = None
            
            if name not in self.mMasks:
                self.mMasks[ name ] = []
                
            self.mMasks[name].append( (int(first_res)-1, int(last_res), residue) )

        if self.mLogLevel >= 1 and self.mErrors > 0:
            self.mStdlog.write("# number of parsing errors: %i\n" % (self.mErrors) )
            
    def __call__( self, peptide_sequence, id ):
        """mask peptide sequence
        """

        if id not in self.mMasks:
            return peptide_sequence, 0
        
        s = list(peptide_sequence.upper())

        nmasked = 0
        for first_res, last_res, residue in self.mMasks[id]:
            if residue:
                for x in range(max(0,first_res), min((len(peptide_sequence),last_res))):
                    if s[x] == residue:
                        s[x] = self.mMaskChar
                        nmasked += 1
            else:
                for x in range(max(0,first_res), min((len(peptide_sequence),last_res))):                
                    s[x] = self.mMaskChar
                nmasked += last_res - first_res

        masked_sequence = "".join(s)
        return masked_sequence, nmasked
    
##------------------------------------------------------------
if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id$", usage = USAGE)

    parser.add_option("-a", "--filename-masks", dest="filename_masks", type="string",
                      help="""filename with mask information."""  )

    parser.add_option("-m", "--masks", dest="methods", type="choice", action="append",
                      choices=("bias", "short", "tmhmm", "coils", "all"),
                      help="sections to mask. Supply more than one or 'all' for all." )

    parser.set_defaults(
        filename_masks = None,
        methods=[],
        )

    options, args = Experiment.Start( parser )

    if len(options.methods) == 0:
        raise "please supply at least one masking method."
    
    if not options.filename_masks:
        raise "please supply filename with masking information."
    
    masker = Masker( open(options.filename_masks, "r" ),
                     options.methods,
                     options )

    nmasked_sequences, nmasked_residues = 0, 0
    ninput, noutput, nresidues, nmasked = 0, 0, 0, 0
    
    def processEntry( header, sequence ):

        global ninput, noutput, nmasked_sequences, nmasked_residues, nmasked, nresidues

        ninput += 1
        id = re.match( "(\S+)", header ).groups()[0]
        
        sequence, nmasked = masker( sequence, id )

        if nmasked:
            nmasked_sequences += 1
            if options.loglevel >= 2:
                options.stdlog.write("# %s: %i residues masked\n" % (id, nmasked))
            nmasked_residues += nmasked
        
        options.stdout.write( ">%s\n%s\n" % (header, sequence ))
        nresidues += len(sequence)
        noutput += 1
        
    header = None
    s = []
    for line in sys.stdin:
        if line[0] == "#": continue
        if line[0] == ">":
            if header:
                processEntry( header, "".join(s) )
            header = line[1:-1]
            s = []
            continue
        s.append( line[:-1] )

    if header:
        processEntry( header, "".join(s) )            

    if options.loglevel >= 1:
        options.stdlog.write("# ninput=%i, noutput=%i, nresidues=%i, nmasked_sequences=%i, nmasked_residues=%i\n" %\
                             (ninput, noutput, nresidues, nmasked_sequences, nmasked_residues))
            
    Experiment.Stop()
    
