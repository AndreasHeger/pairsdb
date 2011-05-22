import unittest, sys, os, re, gzip
from pprint import pprint


class rsdbTest(unittest.TestCase):
    """
    A test class for the feedparser module.
    """

    def setUp(self):
        """
        set up data used in the tests.
        setUp is called before each test function execution.
        """
        
        self.mLevels = ( 30, 40, 50, 60, 70, 80, 90 )

        self.mPatternRsdb = "nrdb%i.table"
        self.mPatternPairsdb = "pairsdb_100x%i.table"
            
    def readRsdb( self, filename ):
        """read nids"""
        
        infile = self.openFile( filename )
        
        nids = []
        for line in infile:
            if line[0] == "#": continue
            nid = int(line.split( "\t" )[0])
            nids.append( nid )
            
        infile.close()

        return nids
            
    def readPairsdb( self, filename ):
        
        infile = self.openFile( filename )

        map_mem2rep, map_rep2mem = {}, {}

        for line in infile:
            if line[0] == "#": continue
            data = line.split( "\t" )
            rep_nid, mem_nid = map( int, [ data[x] for x in (0,4) ] )
            if rep_nid not in map_rep2mem: map_rep2mem[ rep_nid ] = []
            map_rep2mem[ rep_nid ].append( mem_nid )
            if mem_nid in map_mem2rep:
                raise "member nid %i occurs more than once" % mem_nid
            map_mem2rep[ mem_nid ] = rep_nid
        infile.close()

        return map_rep2mem, map_mem2rep
        
    def openFile( self, filename ):
        """open a file for reading.
        
        checks if a zipped file is present.
        """

        if os.path.exists( filename ):
            infile = open( filename, "r" )
        elif os.path.exists( filename + ".gz" ):
            infile = gzip.open( filename + ".gz", "r" )
        else:
            raise "not found"
        return infile

    def runTests( self, filename_rsdb, filename_pairsdb ):
        """check rsdb and pairsdb for consistency.

        The conditions are:

        - reps in pairsdb are a subset of reps in rsdb
        - there is no overlap between reps in rsdb and mems in 
        pairsdb
        """
        
        reps_rsdb = set( self.readRsdb( filename_rsdb ) )

        map_rep2mem, map_mem2rep = self.readPairsdb( filename_pairsdb )
        
        reps_pairsdb = set( map_rep2mem.keys() )
        mems_pairsdb = set( map_mem2rep.keys() )

        ## reps in pairsdb not in rsdb
        extra_reps = reps_pairsdb.difference( reps_rsdb )
        
        mems_in_reps = reps_rsdb.intersection( mems_pairsdb )

        self.assertEqual( len( extra_reps ), 0,
                          msg = "%i reps in %s and not in %s" %\
                              ( len( extra_reps ),
                                filename_pairsdb, 
                                filename_rsdb ))
        
        self.assertEqual( len( mems_in_reps), 0,
                          msg = "%i members in %s also reps in %s" %\
                              (len( mems_in_reps ),
                               filename_pairsdb, 
                               filename_rsdb ))

    def tearDown(self):
        """
        tear down any data used in tests
        tearDown is called after each test function execution.
        """
        pass

    def testRsdb90(self):
        """test on 90 level.
        """
        self.runTests( \
            self.mPatternRsdb % 90,
            self.mPatternPairsdb % 90 )

    def testRsdb40(self):
        """test on 40 level.
        """
        self.runTests( \
            self.mPatternRsdb % 40,
            self.mPatternPairsdb % 40 )

 
if __name__ == '__main__':
    unittest.main()

