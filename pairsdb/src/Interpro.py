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
import sys, re, string, gzip, optparse

import Experiment

USAGE="""
"""

from NrdbAnnotations import NrdbAnnotations

databases = {
    'PRINTS' : 1,
    'PROSITE' : 2,
    'PFAM' : 3,
    'PRODOM' : 4,
    'PROFILE' : 5,
    'PREFILE' : 6,
    'SMART' : 7,
    'TIGRFAMs' : 8,
    'SSF' : 9,
    'PIRSF' : 10,
    'PANTHER' : 11,
    'GENE3D' : 12,
    }
		 
class Interpro (NrdbAnnotations):

    mName = "interpro"

    def __init__ (self, options = None):

	NrdbAnnotations.__init__( self, options )
        
        if options:
            self.mFilenameAnnotations = options.input_filename_annotations
            self.mFilenameDescriptions = options.input_filename_descriptions
            
    #--------------------------------------------------------------------------------------
    def run( self ):
        
        ## read references for uniprot only
        self.readCrossReferences( dbs= (1,) )
        
        self.createAnnotations()

        self.addDescriptions( self.mOutputFilenameAnnotations % self.mParams100,
                              self.mOutputFilenameDescriptions % self.mParams100 )
        self.propagate()

    #--------------------------------------------------------------------------------------        
    def createAnnotations( self ):
        """parse annotations from interpro xml file.
        """

        if self.mLogLevel >= 1:
            self.mStdlog.write("# parsing annotations.\n" )
        
        infile = self.openInputFile( self.mFilenameAnnotations )
        outfile = self.openOutputFile( self.mOutputFilenameAnnotations % self.mParams100 )
        
        if outfile:

            ninput, noutput, nskipped = 0, 0, 0

            self.resetStats()

            for line in infile:

                x = re.search("<protein id=\"(\S+)\" name=\"(\S+)\" length=\"(\d+)\"", line)
                if x:

                    if self.mTest and ninput > self.mTest: break

                    ninput += 1

                    (acc, id, length) = x.groups()

                    if "-" in acc: acc = acc[:acc.find("-")]
                    if "-" in id: id = id[:id.find("-")]
                    
                    length = string.atoi(length)
                    
                    if id not in self.mMapPid2Nid and acc not in self.mMapAcc2Nid:
                        self.mNotFound += 1
                        if self.mLogLevel >= 1:
                            self.mStdlog.write( "# no nids found for acc=%s id=%s in line: %s\n" % (acc, id, line[:-1]))
                        nid = None
                    else:
                        if id in self.mMapPid2Nid:
                            nid = self.mMapPid2Nid[id]
                        else:
                            nid = self.mMapAcc2Nid[acc]

                    continue

                x = re.search("<ipr id=\"(IPR\S+)\" ", line)
                if x: 
                    interpro_id = x.groups()[0]
                    continue

                x = re.search('<match id="(\S+)" name="(.+)" dbname="(\S+)" status="(\S+)" evd="(\S+)"',line)
                if x: 
                    (db_key, db_annotation, db_name, status, evidence) = x.groups()
                    if self.mLogLevel >= 4:
                        self.mStdlog.write("# found match for %s (acc=%s,nid=%i): name=%s, db=%s, status=%s, evd=%s\n" %\
                                               (id, acc, nid, db_key, db_name, status, evidence))
                    continue

                if re.search("<lcn", line):

                    if not nid: continue

                    x = re.search('<lcn start="(\d+)" end="(\d+)" score="(.+)"',line)
                    if x:
                        (start, end, score) = x.groups()
                    else:
                        if self.mLogLevel >= 1:
                            self.mStdlog.write( "# no location found in line: %s\n" % (line[:-1]))
                        self.mFailed += 1
                        continue
                    
                    if start == "0" and end == "0":
                        start = "1"
                        end = str(length)
                        # some entries do not have a location, for example:
                        # <location status="N" evidence="Man_class" />
                        # these are not included here.

                    if not databases.has_key(db_name):
                        if self.mLogLevel >= 1:
                            self.mStdlog.write( "# db_name %s not found for: %s\n" (db_name, ",".join( map(str, nid,id,name,db_key,db_annotation,start,end,status,evidence))))
                        self.mFailed += 1
                        continue
                    else:
                        self.mFromPid += 1

                        db_code = databases[db_name]

                        ali = "+%i" % (string.atoi(end)-string.atoi(start)+1)
                        if nid != None:
                            outfile.write( string.join((
                                str(nid),
                                start,
                                end,
                                ali,
                                "%s.%i.%s-%s" % (interpro_id, db_code, start, end),
                                start,
                                end,
                                ali,
                                interpro_id,
                                str(db_code),
                                db_annotation,
                                db_key,
                                status,
                                evidence), "\t") + "\n")
                            noutput += 1
                            outfile.flush()
                        else:
                            nskipped += 1
                            

            outfile.close()

            if self.mLogLevel >= 1:
                self.printStats()
                self.mStdlog.write("# parsing domains: ninput=%i, noutput=%i, nskipped=%i\n" % (ninput, noutput, nskipped) )
            
    #--------------------------------------------------------------------------------------
    def addDescriptions( self, 
                         input_filename_annotations,
                         output_filename_descriptions,
                         ):
        """parser description from interpro-file.

        call this method after parseDescriptions().
        """

        if self.mLogLevel >= 1:
            self.mStdlog.write("# parsing descriptions.\n" )
            
        infile = self.openInputFile( self.mFilenameDescriptions )
        outfile = self.openOutputFile( output_filename_descriptions )

        if outfile:
            counts = self.getCountsPerFamily( input_filename_annotations )

            ninput, noutput, nskipped, ntypes = 0, 0, 0, 0

            for line in infile:

                if line[:3] == "IPR":
                    ninput += 1
                    id = line[:9]
                    description = line[10:-1]

                    if id not in counts:
                        if self.mLogLevel >= 1:
                            self.mStdlog.write("# missing family in annotations: %s\n" % id)
                        nskipped += 1
                        continue

                    outfile.write( "\t".join( map(str, (id, 
                                                        counts[id].mNUnits,
                                                        counts[id].mNSequences,
                                                        counts[id].mNResidues,
                                                        counts[id].mNResidues / counts[id].mNUnits,
                                                        description,
                                                        type))) + "\n")
                    noutput += 1
                else:
                    ntypes += 1
                    type = line[:-1]

            outfile.close()
            infile.close

            if self.mLogLevel >= 1:
                self.printStats()
                self.mStdlog.write("# parsing descriptions: ninput=%i, noutput=%i, nskipped=%i, ntypes=%i\n" % (ninput, noutput, nskipped, ntypes) )

#---------------------------------------------------------------------------------------------------

if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id$", usage=USAGE )

    parser.add_option( "--input-filename-annotations", dest="input_filename_annotations", type="string" ,
                       help="INPUT filename of with interpro data (matches).")

    parser.add_option( "--input-filename-descriptions", dest="input_filename_descriptions", type="string" ,
                       help="INPUT filename of with interpro data (descriptions).")

    parser.set_defaults( 
        input_filename_annotations = "match_complete.xml.gz",
        input_filename_descriptions = "entry.list",
        )

    Interpro().addOptions( parser )

    (options, args) = Experiment.Start( parser )

    interpro = Interpro( options )

    interpro.run()
    
    Experiment.Stop()
        



        




