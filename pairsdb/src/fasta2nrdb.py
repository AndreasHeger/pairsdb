####
####
##
## Project PairsDB
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: Nrdb.py,v 1.1.1.1 2002/07/02 10:46:57 heger Exp $
##
##
####
####

USAGE="""python fasta2nrdb.py [OPTIONS] filename1 filename2 > nrdb

add one or several fasta files to an nrdb table.

Possible styles:
   swiss:  format is acc|id description - species
   trembl: format is acc|id description - species   
   generic: format is db|acc|id description [species]

Output:

nrdb table: list of non-redundant sequences. 

Sequences are only added to this table, never removed. New sequences 
are appended, while obsolete sequences are flagged (filter = 0). The
fields in this table are:

updated         data of last update of this entry
nid             numerical identifier of the sequence
hid             hash code of the sequence
sequence        sequence
sequencedbs_id  sequence database for preferred reference
accession       accession number of preferred reference
identifier      identifier of preferred reference
description     description of preferred reference
created         time of first creation of this entry
length          sequence length
filter          0=obsolete sequence

references table: list of alternative identifiers for sequences

nid             numerical identifier of sequence
sequencedbs_id  sequence database of reference
accession       accession number of reference
identifier      identifier of reference
description     description associated with reference

"""

import sys, re, string, optparse, time, base64, md5, os, tempfile, shutil

import Experiment

class NrdbEntry:

    def __init__(self):
        self.mNid = 0
        self.mHid = None
        self.mCreated = None
        self.mUpdated = None
        self.mDescription = None
        self.mDatabase = None
        self.mIdentifer = None
        self.mAccession = None
        self.mSequence = None
        self.mFilter = None
        
    def readFromLine(self, line):

        try:
            (self.mUpdated,
             self.mNid,
             self.mHid,
             self.mSequence,
             self.mDatabase,
             self.mAccession,                                    
             self.mIdentifier,
             self.mDescription,
             self.mCreated,
             l,
             self.mFilter ) = line[:-1].split("\t")
        except ValueError:
            raise "parsing error in line: %s" % line[:-1]
        
        self.mNid = int(self.mNid)
        self.mFilter = int(self.mFilter)
        self.mDatabase = int(self.mDatabase)

    def __str__(self):

        return "\t".join( map(str, (self.mUpdated,
                                    self.mNid,
                                    self.mHid,
                                    self.mSequence,
                                    self.mDatabase,
                                    self.mAccession,                                    
                                    self.mIdentifier,
                                    self.mDescription,
                                    self.mCreated,
                                    len(self.mSequence),
                                    self.mFilter )))

class Reference:
    def __init__(self, nid=0, db=0, acc=None, id=None, desc=None):
        self.mNid = nid
        self.mDatabase = db
        self.mAccession = acc
        self.mIdentifier = id
        self.mDescription = desc
        
    def __str__(self):
        return "\t".join( map(str,(self.mNid, self.mDatabase, self.mAccession, self.mIdentifier, self.mDescription)))

    def readFromLine(self, line):
        self.mNid, self.mDatabase, self.mAccession, self.mIdentifier, self.mDescription = line[:-1].split("\t")
        self.mNid = int(self.mNid)
        self.mDatabase = int(self.mDatabase)
        
class Species:
    def __init__(self):
        self.mNid = 0
        self.mSpecies = 0

    def _str__(self):
        return "\t".join( map(str,(self.mNid, self.mSpecies)))

class Nrdb:

    name		= "Nrdb"	 # name of sender 
    lastline		= ''		 # LAST line while stepping through file
    source_file		= 0		 # file-handle of nrdb

    requirements        = ()             # modules that should have finished before this module is executed

    def __init__ (self, options):

	self.mLogLevel = options.loglevel
        self.mStdout = options.stdout
        self.mStdlog = options.stdlog
        self.mStderr = options.stderr        
        
        self.mInputFilenameNrdb = options.input_filename_nrdb
        self.mInputFilenameSynonyms = options.input_filename_synonyms
        self.mInputFilenameTaxonomy = options.input_filename_taxonomy

        self.mOutputFilenameReferences = options.output_filename_references
        self.mOutputFilenameTaxonomy = options.output_filename_taxonomy
        self.mOutputFilenameQuality = options.output_filename_quality

        if options.tempdir:
            self.mTempdir = tempfile.mkdtemp(dir=options.tempdir)            
        else:
            self.mTempdir = tempfile.mkdtemp()
        if self.mLogLevel >= 2:
            self.mStdlog.write("# temporary directory used: %s\n" % self.mTempdir)
        
        self.mOutputFilenameNrdbNew = self.mTempdir + "/new"
        self.mOutputFilenameNrdbUpdated = self.mTempdir + "/updated"

        self.mOutputFilenameNrdb = options.output_filename_nrdb
        self.mOutputFilenameNrdbFasta = options.output_filename_nrdb_fasta

        ## output files: these are opened while adding new fasta files to the
        ## existing nrdb table.

        self.mOutfileReferences = None
        self.mOutfileTaxonomy = None
        self.mOutfileNew = None
        
        ## set to true, if old entries shall not be marked obsolete
        self.mRescan = options.rescan

        self.readSynonyms()
        self.readTaxonomy()
        
        # legitimate alphabet for sequences (currently restricted by TMHMM)
        self.mAlphabet = "ACDEFGHIKLMNPQRSTVWYXBZUOJ\-"        # UOJ added. Are replaced in Nrdb90_masks

        self.mAliveNids = set()

    #-------------------------------------------------------------------------------------------------------
    def readSynonyms(self):
        """read synonyms from file, alternatively, use hard-coded table."""
        
        if not self.mInputFilenameSynonyms:

            if self.mLogLevel >= 1:
                self.mStdlog.write("# WARNING: using hard-coded map of database names.\n" )
                self.mStdlog.flush()
                
            self.mMapId2Synonyms = {
                1 : ("swall","swiss","swissnew","sptrembl","tremblnew","uniprot", "trembl"),
                2 : ("remtrembl",),
                3 : ("pdb",),
                4 : ("gp,gpnew,gi",),
                5 : ("pir,pironly",),
                6 : ("worm",),
                7 : ("ensembl",),
                8 : ("yeast",),
                9 : ("ncbi",),    
                10 : ("ncbi",),
                11 : ("peptide",),
                12 : ("refseq",),
                }

            self.mMapSynonym2Id = {}
            for key, values in self.mMapId2Synonyms.items():
                for value in values:
                    self.mMapSynonym2Id[value] = key

        else:
            if self.mLogLevel >= 1:
                self.mStdlog.write("# WARNING: reading database names from %s.\n" % self.mInputFilenameSynonyms)
                self.mStdlog.flush()

            infile = open(self.mInputFilenameSynonyms, "r")
            
            for line in infile:
                if line[0] == "#": continue
                key, synonym = line[:-1].split("\t")
                self.mMapSynonym2Id[synonym] = int(key)
                
            infile.close()

        self.mSequencedbsOrder = [3, 1, 12, 5, 9, 4, 6, 7, 10, 8, 2, 11]
        
        self.mMapDatabase2Preference = {}
        for x,y in enumerate(self.mSequencedbsOrder):
            self.mMapDatabase2Preference[y] = x

        ## curated databases: pdb, refseq
        self.mCuratedDatabases = set( [3, 12,] )

        if self.mLogLevel >= 1:
            for key, value in self.mMapSynonym2Id.items():
                self.mStdlog.write("# %s: %i\n" % (key, value))

    #-------------------------------------------------------------------------------------------------------
    def readTaxonomy( self ):
        """read taxonomy information.
        
        Reads in also synonyms and common names.
        """
        
        infile = open(self.mInputFilenameTaxonomy, "r")
        self.mMapSpecies2TaxId = {}
        for line in infile:

            (tax_id, name, full_name, name_class) = re.split("\t\|\t", line[:-3])[:4]
            self.mMapSpecies2TaxId[name.lower()] = int(tax_id)
            
        infile.close()

        if self.mLogLevel >= 1:
            self.mStdlog.write("# read %i species names\n" % (len(self.mMapSpecies2TaxId)))

    #-------------------------------------------------------------------------------------------------------
    def calculateHID ( self, sequence ):
        """calculate HID for a sequence."""
        # do the encryption
        h = md5.new(sequence).digest()

        # map to printable letters: hid has length 22, so the padded '=' are
        # truncated. You have to add them, if you ever want to decode,
        # but who would do such a thing :=)

        r = base64.encodestring(h)[0:22]

        # finally substitute some characters:
        # '/' for '_', so we have legal file names
        # '[' for '+' and ']' for '=' for internet-applications

        hid = string.replace(r  , '/', '_') 
        hid = string.replace(hid, '+', '[') 
        hid = string.replace(hid, '=', ']') 

        return hid
    
    #-------------------------------------------------------------------------------------------------------
    def openFile( self, filename ):
        """open a fasta file for reading."""
        try:
            self.mInfile = open(filename, "r" )
        except IOError:
            self.mStderr.write("# FATAL: could not open file %s\n" % filename)
            sys.exit(1)
            
        self.mLastLine = None

    #-------------------------------------------------------------------------------------------------------
    def closeFile( self ):
        """close currently opened fasta file."""
        self.mInfile.close()
        
    #-------------------------------------------------------------------------------------------------------
    def getNextEntry( self ):
        """read next entry in source file.

	note: assuming no-hobo fasta-format, i.e. no special characters and no spaces, just
	plain sequences. Since this routine is called often, I do not want to do any pattern
	matching.
        """

	# if we are not on the beginning of a fasta entry, go there
	if (not self.lastline) or (self.lastline[0] != ">"):
            while 1:
                line = self.mInfile.readline()
                if not line: break
                if line[0] == ">": break
	else:
	    line = self.lastline
	    
	if not line: return ('','')
	    
	description = line[1:-1]

        sequence = []

	while 1:
	    line = self.mInfile.readline();
	    if not line: break
	    if line[:1] == ">": break
            sequence.append( re.sub("\s", "", line[:-1] ) )
	    
	self.lastline = line

	return (description, "".join(sequence))
    
#------------------------------------< Methods for updating/inserting/deleting>-----------------------------

    def parseDescriptionLine( self, header, style ):

        parsed = []
        unparsed = []

        if style == "generic":
            data = string.split( header, '//:')	# split at /: and skip the three characters ('>/:')

            for t in data:
                m  = re.compile( '^(\S+)\|(\S+)\|(\S+) ([\S\s]*)' ).search(t)
                if m:
                    (db, ac, id, de) = m.groups()
                else:
                    m  = re.compile( '^(\S+)\|(\S+)\|(\S+)' ).search(t)
                    if m:
                        (db, ac, id) = m.groups()
                        de = ""
                    else:
                        unparsed.append(t)
                        continue

                de = string.strip(de)				# remove white spaces and newlines
                if len(de) > 0 and de[len(de) - 1] == '.':      # remove trailing . 
                    de = de[:-1]
                de = re.sub("\t", "", de)
                if db in self.mMapSynonym2Id:
                    parsed.append( Reference(0, db=self.mMapSynonym2Id[db], acc=ac, id=id, desc=de) )

                else:
                    if self.mLogLevel >= 1:
                        self.mStdlog.write("# WARNING: unknown database: %s\n" % db)
                    unparsed.append(t)

        elif style in ("swiss", "trembl", "uniprot" ):

            db = style
            m  = re.compile( '^(\S+)\|(\S+) ([\S\s]*)' ).search(header)
            if m:
                (ac, id, de) = m.groups()
            else:
                m  = re.compile( '^(\S+)\|(\S+)' ).search(t)
                if m:
                    (ac, id) = m.groups()
                    de = ""
                else:
                    unparsed.append(header)
                    return parsed, unparsed
                
            de = string.strip(de)			    # remove white spaces and newlines
            if len(de) > 0 and de[len(de) - 1] == '.':      # remove trailing . 
                de = de[:-1]
            de = re.sub("\t", "", de)
            
            if db in self.mMapSynonym2Id:
                parsed.append( Reference(0, db=self.mMapSynonym2Id[db], acc=ac, id=id, desc=de) )

            else:
                if self.mLogLevel >= 1:
                    self.mStdlog.write("# WARNING: unknown database: %s\n" % db)
                unparsed.append(header)

        else:
            raise "unknown style %s" % style
        
        return parsed, unparsed

    #-------------------------------------------------------------------------------------------------------
    def writeReferences( self, nid, references ):
	"""write references to outputfile.
        """
        if not self.mOutfileReferences: return
        for reference in references:
            reference.mNid = nid
            self.mOutfileReferences.write( "%s\n" % str(reference))

    #-------------------------------------------------------------------------------------------------------
    def writeTaxonomy( self, nid, references ):
	"""write references to outputfile.
        """
        if not self.mOutfileTaxonomy: return False

        species_list = []
        
        for reference in references:

            if reference.mDatabase == 1:
                ## uniprot: only one species, but with synonyms
                data = reference.mDescription.split("-")
                if len(data) != 2:
                    continue
                else:
                    t = data[1].strip()
                    if "(" not in t:
                        species_list.append( t )
                    else:
                        x = re.match( "([^()]+) (\([^()]+\))*", t)
                        if not x:
                            if self.mLogLevel >= 2:
                                self.mStdlog.write("# ERROR: could not extract taxa from %s\n" % reference.mDescription )
                        else:
                            taxa = x.groups()
                            species_list.append( taxa[0] )
                            for x in taxa[1:]:
                                if x:
                                    # remove brackets and space
                                    species_list.append( x[1:-1] )
            else:
                ## others: species in [], possibly multiple species
                x = re.search("\[([^\[\]]+)\]", reference.mDescription)
                if x:
                    species_list += x.groups()
                
        tax_ids = set()
        missed = set()

        for species in species_list:
            if species.lower() not in self.mMapSpecies2TaxId:
                missed.add( species )
                continue
            tax_id = self.mMapSpecies2TaxId[species.lower()]
            tax_ids.add( tax_id )
            
        if len(tax_ids) == 0 and len(missed) > 0:
            if self.mLogLevel >= 2:
                self.mStdlog.write("# ERROR: could not find any taxonomic assignment for nid=%i, species=%s\n" % (nid, str(missed)))

        for tax_id in tax_ids:
            self.mOutfileTaxonomy.write( "%i\t%i\n" % (nid, tax_id) )

        return len(tax_ids) > 0
    #-------------------------------------------------------------------------------------------------------
    def getPreferredReference( self, references ):
        """get preferred reference from a list of references."""
        l = []
        for x in references:
            l.append( (self.mMapDatabase2Preference[x.mDatabase], x ) )

        l.sort()
        return l[0][1]
    
    #-------------------------------------------------------------------------------------------------------
    def addFasta( self, filename, style = "generic"  ):
        """read through fasta formatted file and add sequences.
        """

        if self.mLogLevel >= 1:
            options.stdlog.write("# adding %s to nrdb table\n" % filename )
            options.stdlog.flush()

        time_0 = time.time()
        
        self.openFile( filename )
 	date = time.strftime("%Y-%m-%d", time.localtime(time.time()))

        rx = re.compile("[^" + self.mAlphabet + "]")

        nnoid, nunknown, ndatabase = 0,0,0
        nmasked, ninput, nadded, nredundant = 0,0,0,0
        noutput_with_taxonomy = 0

	while 1:

	    # iterate through file
	    (identifier, sequence) = self.getNextEntry()

	    if not identifier: break
            
            if self.mLogLevel >= 3:
                self.mStdlog.write( "# processing %s\n" % identifier )
            
            ninput += 1

            # check alphabet
            if rx.search(sequence):
                nmasked += 1
                if self.mLogLevel >= 1:
                    options.stdlog.write("# WARNING: wrong alphabet in sequence %s - substituting illegal characters with X\n" % (identifier))
                    options.stdlog.write("# original sequence: %s\n" % (sequence))                    
                sequence = rx.sub( "X", sequence )
                
            (references, unparsed) = self.parseDescriptionLine( identifier, style )

            ## check if unparsed entries
            if unparsed:
                self.mStderr.write( "# parsing error in line: %s\n" % identifier)
                self.mStderr.write( "# unparsed: %s\n" % unparsed )
                sys.exit(1)

            ## check databases
            if len(references) == 0:
                if self.mLogLevel >= 1:
                    self.mStdlog.write("# WARNING: refused to add sequence %s, because no identifier can be retrieved.\n" % identifier)
                nnoid += 1
                continue

	    ## calculate hid of entry and get nid
	    hid = self.calculateHID( sequence )

            if hid not in self.mMapHid2Nid:

                nadded += 1
                ## new entry, write to file
                self.mMaxNid += 1
                nid = self.mMaxNid
                self.mMapHid2Nid[hid] = nid

                entry = NrdbEntry()

                entry.mNid = nid
                entry.mHid = hid
                entry.mSequence = sequence

                # get preferred reference
                pr = self.getPreferredReference( references )
                
                (entry.mDatabase, entry.mAccession, entry.mIdentifier, entry.mDescription) =\
                                  (pr.mDatabase, pr.mAccession, pr.mIdentifier, pr.mDescription)

                entry.mUpdated = time.strftime( "%Y%m%d%H%M%S", time.localtime() )
                entry.mCreated = date
                entry.mFilter = 100
                self.mOutfileNew.write( "%s\n" % str( entry ) )
            else:
                
                nredundant += 1
                ## old entry
                nid = self.mMapHid2Nid[hid]

            self.writeReferences( nid, references )
            if self.writeTaxonomy( nid, references ):
                noutput_with_taxonomy += 1

            ## mark nid as alive
            self.mAliveNids.add( nid )

	# end of while-----------------------------------------------------------
        self.closeFile()
        
        if self.mLogLevel >= 1:
            self.mStdlog.write("# %s: added %i entries to nrdb table in %i seconds\n" % (filename, nadded, time.time() - time_0))
            self.mStdlog.write("# %s: ninput=%i, nadded=%i, nredundant=%i, noutput_with_taxonomy=%i\n" % (filename, ninput, nadded, nredundant, noutput_with_taxonomy))
            self.mStdlog.write("# %s: WARNINGS: nmasked=%i, ndatabase=%i\n" % (filename, nmasked, ndatabase))
            self.mStdlog.write("# %s: ERRORS: noids=%i\n" % (filename, nnoid))

    #-------------------------------------------------------------------------------------------------------
    def readHids( self ):
        """read hids from nrdb table."""
        
        infile = open(self.mInputFilenameNrdb, "r" )

        self.mMapHid2Nid = {}
        self.mAliveNids = set()
        self.mMaxNid = 0
        if self.mLogLevel >= 1:
            self.mStdlog.write( "# reading hids from file %s\n" % self.mInputFilenameNrdb )
            
        for line in infile:
            data = line[:-1].split("\t") 
            nid, hid = data[1], data[2]
            nid = int(nid)
            self.mMapHid2Nid[hid] = nid
            self.mMaxNid = max(self.mMaxNid, nid )
            
        infile.close()

        if self.mLogLevel >= 1:
            self.mStdlog.write( "# read %i hids\n" % len( self.mMapHid2Nid) )

    #-------------------------------------------------------------------------------------------------------
    def prepare(self):
        """prepare for building a new nrdb database."""

        if self.mOutputFilenameReferences:
            self.mOutfileReferences = open(self.mOutputFilenameReferences, "w")
            if self.mLogLevel >= 1:
                self.mStdlog.write("# output with references goes to: %s\n" % self.mOutputFilenameReferences)
        else:
            self.mOutfileReferences = None
            
        if self.mOutputFilenameTaxonomy:
            self.mOutfileTaxonomy = open(self.mOutputFilenameTaxonomy, "w")
            if self.mLogLevel >= 1:
                self.mStdlog.write("# output with taxonomic assignments goes to: %s\n" % self.mOutputFilenameTaxonomy)
        else:
            self.mOutfileTaxonomy = None

        self.mOutfileNew = open( self.mOutputFilenameNrdbNew, "w" )
        
        self.readHids()

    #-------------------------------------------------------------------------------------------------------
    def indexReferences( self ):
        """index references on filesystem."""

        if self.mLogLevel >= 1:
            self.mStdlog.write("# indexing references.\n" )
            self.mStdlog.flush()
            
        os.system( "sort -T %s -k1,1n %s > %s_tmp" % (self.mTempdir, self.mOutputFilenameReferences, self.mOutputFilenameReferences ))
        os.system( "mv %s_tmp %s" % (self.mOutputFilenameReferences, self.mOutputFilenameReferences ))

        infile = open(self.mOutputFilenameReferences, "r" )
        self.mMapNid2Pos = {}
        reference = Reference()
        while 1:
            pos = infile.tell()
            line = infile.readline()
            if not line: break

            reference.readFromLine( line )
            self.mMapNid2Pos[reference.mNid] = pos
        infile.close()

        if self.mLogLevel >= 1:
            self.mStdlog.write("# indexed references for %i nids.\n" % len(self.mMapNid2Pos))
            self.mStdlog.flush()

    #-------------------------------------------------------------------------------------------------------
    def updateNrdb( self ):
        """update an old nrdb file."""
        
        self.mOutfileUpdated = open( self.mOutputFilenameNrdbUpdated, "w" )
        
        self.indexReferences()
        infile_references = open(self.mOutputFilenameReferences, "r")
        
        infile = open(self.mInputFilenameNrdb, "r" )
        nrdb = NrdbEntry()

        if self.mRescan:
            if self.mOutputFilenameTaxonomy:
                self.mOutfileTaxonomy = open(self.mOutputFilenameTaxonomy, "w")
                if self.mLogLevel >= 1:
                    self.mStdlog.write("# output with taxonomic assignments goes to: %s\n" % self.mOutputFilenameTaxonomy)
            else:
                self.mOutfileTaxonomy = None

        ninput, nobsolete, noutput, noutput_with_taxonomy = 0, 0, 0, 0
        for line in infile:
            ninput += 1
            nrdb.readFromLine( line )

            ## test whether this entry shall be processed
            is_obsolete = False
            if self.mRescan:
                if nrdb.mFilter == 0:
                    is_obsolete = True
            elif nrdb.mNid not in self.mAliveNids:
                nrdb.mFilter = 0
                is_obsolete = True
            else:
                nrdb.mFilter = 100
                is_obsolete = False
                
            ## obsolete entries are simply output
            if is_obsolete:
                self.mOutfileUpdated.write( "%s\n" % (str(nrdb) ) )
                nobsolete += 1
                continue

            ## get/set preferred reference
            infile_references.seek( self.mMapNid2Pos[nrdb.mNid] )
            references = []
            while 1:
                line = infile_references.readline()
                if not line: break
                r = Reference()
                r.readFromLine( line )
                if r.mNid != nrdb.mNid: break
                references.append( r )

            rp = self.getPreferredReference( references )
            nrdb.mDatabase, nrdb.mIdentifier, nrdb.mAccession, nrdb.mDescription = \
                            rp.mDatabase, rp.mIdentifier, rp.mAccession, rp.mDescription

            ## write taxonomic assignments
            if self.mRescan:
                if self.writeTaxonomy( nrdb.mNid, references ):
                    noutput_with_taxonomy += 1
                
            ## write updated entry
            self.mOutfileUpdated.write( "%s\n" % (str(nrdb) ) )
            noutput +=1

        self.mOutfileUpdated.close()
        
        if self.mLogLevel >= 1:
            self.mStdlog.write("# updated file %s\n" % self.mInputFilenameNrdb )
            self.mStdlog.write("# ninput=%i, nupdated=%i, nobsolete=%i, noutput_with_taxonomy=%i\n" % (ninput, noutput, nobsolete, noutput_with_taxonomy) )

        if self.mOutfileTaxonomy:
            self.mOutfileTaxonomy.close()

    #-------------------------------------------------------------------------------------------------------
    def outputQuality( self, outfile_quality, nrdb ):
        """output quality status for an nrdb entry.
        """
        txt = nrdb.mDescription.lower()

        is_trusted= False
        is_complete = not re.search( "fragment", txt )
        is_nonhypothetical = not re.search( "hypothetical", txt )
        is_nonartificial = not re.search( "artificial", txt )
        is_nonsynthetic = not (re.search( "synthetic", txt ) and not \
            (re.search("biosynthetic",txt) or  re.search("photosynthetic", txt) ))
        is_genuine = is_complete and is_nonhypothetical and is_nonartificial

        ## curated: check for swissprot identifier and/or databas
        is_curated = nrdb.mDatabase in self.mCuratedDatabases or re.search("_", nrdb.mIdentifier)

        def __yn( x ):
            if x: return "Y"
            return "N"
        
        outfile_quality.write( "\t".join( [str(nrdb.mNid)] +\
                                              map( __yn, \
                                                       ( is_complete,
                                                         is_trusted,
                                                         is_genuine,
                                                         is_nonhypothetical,
                                                         is_nonsynthetic,
                                                         is_nonartificial,
                                                         is_curated ) ) ) + "\n" )

    #-------------------------------------------------------------------------------------------------------
    def createNrdb( self ):
        """create new nrdb files.

        The following files are created:

        an nrdb table containing all entries (alive and obsolete).
        an nrdb fasta file containing all alive entries.
        an nrdb quality file containing all alive entries.

        This also servers as a quality control step.
        """

        if self.mLogLevel >= 1:
            self.mStdlog.write("# creating table and fasta file of nrdb.\n")
            self.mStdlog.write("# output goes to: %s and %s\n" % (self.mOutputFilenameNrdb, self.mOutputFilenameNrdbFasta ))
            
        outfile_nrdb = open( self.mOutputFilenameNrdb, "w" )
        outfile_fasta = open( self.mOutputFilenameNrdbFasta, "w" )
        outfile_quality = open( self.mOutputFilenameQuality, "w" )

        nrdb = NrdbEntry()
        
        def __write(infile):
            nobsolete, nalive = 0, 0
            for line in infile:
                nrdb.readFromLine( line )
                outfile_nrdb.write( "%s\n" % str(nrdb) )                
                if nrdb.mFilter == 0:
                    nobsolete += 1
                    continue
                if not self.mRescan:
                    if nrdb.mFilter != 100:
                        print str(nrdb)
                    assert( nrdb.mFilter == 100 )

                self.outputQuality( outfile_quality, nrdb )

                nalive += 1
                outfile_fasta.write(">%s\n%s\n" % (nrdb.mNid, nrdb.mSequence ) )
            return nobsolete, nalive
        
        infile = open(self.mOutputFilenameNrdbUpdated, "r" )
        nobsolete1, nupdated = __write(infile)
        infile.close()
        
        if os.path.exists( self.mOutputFilenameNrdbNew ):
            infile = open(self.mOutputFilenameNrdbNew, "r" )
            nobsolete2, nadded = __write(infile)
            infile.close()
        else:
            nobsolete2, nadded = 0, 0

        outfile_fasta.close()
        outfile_nrdb.close()

        if self.mLogLevel >= 1:
            self.mStdlog.write("# finished creating nrdb table and fasta file.\n")
            self.mStdlog.write("# ntotal=%i, noriginal=%i, nobsolete=%i, nupdated=%i, nadded=%i\n" % \
                               (nupdated + nadded,
                                nupdated + nobsolete1,
                                nobsolete1, 
                                nupdated, nadded))
                
    #-------------------------------------------------------------------------------------------------------
    def finish(self):
        """perform cleanup entries.
        """
        if self.mOutfileReferences:
            self.mOutfileReferences.close()
        if self.mOutfileTaxonomy:
            self.mOutfileTaxonomy.close()
        if self.mOutfileNew:
            self.mOutfileNew.close()
        
        ## update previous entries
        self.updateNrdb()

        ## create new nrdb table
        self.createNrdb()

        ## clean up temporary files
        shutil.rmtree( self.mTempdir )

#--------------------------------------< end of class definition >-----------------------

if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id$", usage=USAGE )

    parser.add_option( "-n", "--input-filename-nrdb", dest="input_filename_nrdb", type="string" ,
                       help="INPUT filename of an old nrdb table.")

    parser.add_option( "--input-filename-taxonomy", dest="input_filename_taxonomy", type="string" ,
                       help="INPUT filename mapping species names to taxonomy ids.")

    parser.add_option( "-t", "--tempdir", dest="tempdir", type="string" ,
                       help="temporary directory use. Default is system default [/tmp].")

    parser.add_option( "--output-filename-references", dest="output_filename_references", type="string" ,
                       help="OUTPUT filename of cross-references.")

    parser.add_option( "--output-filename-preferred-references", dest="output_filename_preferred_references", type="string" ,
                       help="OUTPUT filename of preferred references.")

    parser.add_option( "--output-filename-taxonomy", dest="output_filename_taxonomy", type="string" ,
                       help="OUTPUT filename of taxonomy assignments.")

    parser.add_option( "--output-filename-nrdb", dest="output_filename_nrdb", type="string" ,
                       help="OUTPUT filename of the new nrdb table.")

    parser.add_option( "--output-filename-nrdb-fasta", dest="output_filename_nrdb_fasta", type="string" ,
                       help="OUTPUT filename with updated nrdb entries.")

    parser.add_option( "--rescan", dest="rescan", action="store_true",
                       help="rescan an existing nrdb table.""")

    parser.set_defaults( input_filename_nrdb = None,
                         input_filename_synonyms = None,
                         input_filename_taxonomy = "names.dmp",
                         output_filename_references = "references.table",
                         output_filename_taxonomy = "taxonomy_assignments.table",
                         output_filename_quality = "nrdb_quality.table",
                         output_filename_nrdb = "nrdb.table",
                         output_filename_nrdb_fasta = "nrdb.fasta",
                         tempdir = None,
                         rescan = False,
                         )

    (options, args) = Experiment.Start( parser, add_pipe_options = True )

    nrdb = Nrdb( options )

    if options.rescan:
        nrdb.finish()
    else:
        nrdb.prepare()

        if len(args) % 2 != 0:
            raise "please supply style and fasta information."

        for x in range(0, len(args), 2):
            filename = args[x]
            style = args[x+1]
            nrdb.addFasta( filename, style )

        nrdb.finish()

    Experiment.Stop()
        
    
    

    
