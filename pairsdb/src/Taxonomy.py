################################################################################
#
#   Taxonomy table mirror
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
import sys, string, re, os, optparse, collections

USAGE="""python Taxonomy.py --action=action [OPTIONS]

build and query ncbi taxonomy. The NCBI taxonomy can be downloaded from

wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xvzf taxdump.tar.gz
"""

import Experiment

def trueFalse( v ): 
    if v: return "t" 
    else: return "f"

class Node:
    def __init__(self,
                 tax_id, 
                 parent, 
                 rank,
                 division,
                 name,
                 complete_genome_flag,
                 node_id,
                 min_node_id,
                 num_children,
                 num_leaves ):
        self.mTaxId = tax_id
        self.mParent = parent 
        self.mRank = rank
        self.mDivision = division
        self.mName = name
        self.mIsComplete = complete_genome_flag
        self.mNodeId = node_id
        self.mMinNodeId = min_node_id
        self.mNumChildren = num_children
        self.mNumLeaves = num_leaves

    def __str__(self):
        return "\t".join( map(str, (self.mTaxId, 
                                    self.mParent, 
                                    self.mRank,
                                    self.mDivision,
                                    self.mName,
                                    trueFalse( self.mIsComplete), 
                                    self.mNodeId,
                                    self.mMinNodeId,
                                    self.mNumChildren,
                                    self.mNumLeaves ) ) )


class Taxonomy:
    """wrapper around a taxonomy structure kept in a postgres database.
    """
    def __init__ ( self ): 
        pass

    ##--------------------------------------------------------------------
    def createTable( self ):
        self.mTree = {}

    ##--------------------------------------------------------------------
    def getChildrenTaxId( self, tax_id ):
        return [ self.mNodes[x].mTaxId for x in self.mMapParent2Children[ tax_id ] ]

    ##--------------------------------------------------------------------
    def setNode( self, tax_id, node_id ):
        """update node statistics for a node.
        """
        children = [ self.mNodes[x] for x in self.mMapParent2Children[ tax_id ] ]
        t = self.mNodes[tax_id]
        t.mNodeId = node_id
        t.mNumChildren = len(children)
        t.mMinNodeId = min( [x.mNodeId for x in children ] )
        t.mNLeaves = sum( [x.mNumLeaves for x in children ] )

    ##--------------------------------------------------------------------
    def setLeaf( self, tax_id, node_id ):
        t = self.mNodes[tax_id]
        t.mNodeId = node_id
        t.mMinNodeId = node_id
        t.mNumChildren = 1
        t.mNumLeaves = 1

    ##--------------------------------------------------------------------
    def getRoot( self ):
        """returns root of the tree.
        """
        return 1

    ##--------------------------------------------------------------------        
    def calculateTree( self ):
        """compute node ids in tree.
        
        The node-ids impose a topological sort order on the nodes allowing
        to query all children of a node in the tree. This is useful for storing
        the tree in a database.

        Simply does a depth-first traversal of taxonomy tree using a stack structure. 
        """

        stack = []
        stack.append( self.getRoot() )
        stack.append( self.getRoot() )         
        
        node = 1
        visited = {}
        
        while stack:
            
            tax_id = stack.pop()
            
            if visited.has_key(tax_id):
                self.setNode( tax_id, node )
                node+=1
                continue
            else:
                visited[tax_id] = 1
                
            children = self.getChildrenTaxId( tax_id )

            if not children:
                stack.pop()
                self.setLeaf( tax_id, node )
                node+=1                
            else:
                for child in children:
                    stack.append( child )
                    stack.append( child )                    
                
    ##--------------------------------------------------------------------------
    def loadFromFiles( self, infile_nodes, infile_divisions, infile_names ):
        """Load taxonomy from NCBI files using python dictionaries.
        """
        
        ## read divisions
        divisions = {}
        for line in infile_divisions:
            (id, abbreviation, explicit) = re.split("\t\|\t", line[:-1])[:3]
            divisions[id] = abbreviation

        ## read names
        names = {}
        for line in infile_names:

            (tax_id, name, full_name, name_class) = re.split("\t\|\t", line[:-3])[:4]
            
            if name_class == "scientific name":
                names[int(tax_id)] = name

        self.mNodes = {}
        self.mMapParent2Children = collections.defaultdict( list )
        for line in infile_nodes:
            
            fields = re.split("\t\|\t", line[:-3])
            (tax_id, parent, rank,
             embl_code,
             division, div_flag,
             gc_id, gc_flag,
             mgc_id, mgc_flag,
             hidden_gb_flag,
             hidden_subtree_flag, comment) = fields
            
            tax_id, parent = map(int, (tax_id, parent) )
            self.mNodes[ tax_id ] = Node( tax_id, 
                                          parent, 
                                          rank,
                                          divisions[division],
                                          names[tax_id],
                                          False,
                                          0, 0, 0, 0)

            self.mMapParent2Children[parent].append( tax_id )

        self.calculateTree()

    def writeToFile( self, outfile ):
        
        for tax_id in sorted(self.mNodes.keys()):
            outfile.write( str(self.mNodes[tax_id]) + "\n" )

class TaxonomyPostgres( Taxonomy ):
    """ """
    def __init__( self, dbhandle, table_name, *args, **kwargs ):
        Table.__init__( *args, **kwargs )

        self.mDbHandle = dbhandle
        
        self.name   =  table_name
        self.fields = (
            ('tax_id', 'INT NOT NULL'),
            ('parent_id', 'INT NOT NULL'),
            ('rank', 'VARCHAR(20) NOT NULL'),
            ('division', 'CHAR(3) NOT NULL'),
            ("scientific_name", "VARCHAR(100) NOT NULL DEFAULT ''"),
            ('complete_genome_flag', 'CHAR(1) NOT NULL'),
            ('node_id', 'INT NOT NULL'),
            ('min_node_id', 'INT NOT NULL'),
            ('num_children', 'INT NOT NULL'),
            ('num_leaves', 'INT NOT NULL'),
            )

        self.indices = ( 'tax_id',
                        'parent_id',
                        'node_id',
                        'rank',
                        )
           
        import pgdb

    ##--------------------------------------------------------------------
    def createTable( self ):
        """create a new table."""
        
        statement = "CREATE TABLE " + self.name + ' ( '
                
        f = map (string.join, self.fields )
        statement = statement + string.join( f, ',')			       	# add fields

        statement = statement + ' ) '
        
        self.execute(statement)

        n = 0
        for index in self.indices:
            n += 1
            statement = "CREATE INDEX %s_index%i ON %s (%s)" % ( re.sub( "[.]", "_", self.name),
                                                                 n,
                                                                 self.name,
                                                                 index )
            
            self.execute(statement)

        self.mDbHandle.commit()
        
        return True

    ##--------------------------------------------------------------------
    def dropTable( self ):
        """drop table."""
        self.execute( "DROP TABLE %s" % self.name )
                     
    ##--------------------------------------------------------------------
    def execute( self, statement ):
        """execute statement and return a cursor with result."""
        cc = self.mDbHandle.cursor()
        cc.execute(statement)
        cc.close()

    ##--------------------------------------------------------------------
    def getCursor( self, statement ):
        """execute statement and return a cursor with result."""
        cc = self.mDbHandle.cursor()
        cc.execute(statement)
        return cc
        
    ##--------------------------------------------------------------------
    def getChildrenTaxId( self, tax_id ):
        statement = "SELECT tax_id FROM %s WHERE parent_id = %i" % (self.name, tax_id)
        cc = self.getCursor( statement )
        result = map(lambda x:x[0], cc.fetchall())
        cc.close()
        return result

    ##--------------------------------------------------------------------
    def setLeaf( self, tax_id, node ):
        statement = """UPDATE %s
        SET node_id = %i, min_node_id = %i, num_children = 1, num_leaves = 1
        WHERE tax_id = %i""" % (self.name, node, node, tax_id)
        return self.execute(statement)

    ##--------------------------------------------------------------------
    def setNode( self, tax_id, node_id ):
        """update node statistics for a node.
        """
        statement = """SELECT MIN(min_node_id), SUM(num_children), SUM(num_leaves)
        FROM %s WHERE parent_id = %i""" % (self.name, tax_id)
        (min_node_id, num_children, num_leaves) = self.getCursor(statement).fetchone()
        
        statement = """UPDATE %s
        SET node_id = %i, min_node_id = %i, num_children = %i, num_leaves = %i
        WHERE tax_id = %i
        """ % (self.name, node_id, min_node_id, num_children+1, num_leaves, tax_id)
        
        return self.execute(statement)

    ##--------------------------------------------------------------------
    def getChildren( self, tax_id ):
        """return children of node."""
        return self.getCursor( "SELECT min_node_id, node_id FROM %s WHERE tax_id = %i" % (self.name, tax_id)).fetchone()

    ##--------------------------------------------------------------------
    def getNodeId( self, tax_id ):
        return self.getCursor( "SELECT node_id FROM %s WHERE tax_id = %i" % (self.name, tax_id)).fetchone()[0]

    ##--------------------------------------------------------------------
    def getDomainForNodeId( self, node_id ):
        result = self.getCursor( """SELECT
        scientific_name FROM %s WHERE rank = 'superkingdom' AND
        %i BETWEEN min_node_id AND node_id""" % (self.name, node_id)).fetchone()
        if result: return result[0]
        else: return None
    

    ##--------------------------------------------------------------------------
    def loadFromFiles( self, infile_nodes, infile_divisions, infile_names ):
        """Load taxonomy from NCBI files."""
        
        ## read divisions
        divisions = {}
        for line in infile_divisions:
            (id, abbreviation, explicit) = re.split("\t\|\t", line[:-1])[:3]
            divisions[id] = abbreviation

        ## read names
        names = {}
        for line in infile_names:

            (tax_id, name, full_name, name_class) = re.split("\t\|\t", line[:-3])[:4]
            
            if name_class == "scientific name":
                names[tax_id] = name

        self.dropTable()
        self.createTable()
        
        template = """INSERT INTO %s (tax_id, parent_id, rank, division, scientific_name,
        complete_genome_flag, node_id, min_node_id, num_children, num_leaves) VALUES (E'%%s')""" % self.name
        
        for line in infile_nodes:
            
            fields = re.split("\t\|\t", line[:-3])
            (tax_id, parent, rank,
             embl_code,
             division, div_flag,
             gc_id, gc_flag,
             mgc_id, mgc_flag,
             hidden_gb_flag,
             hidden_subtree_flag, comment) = fields

            statement =  template % "',E'".join( (tax_id, parent, rank,
                                                  divisions[division].encode("string_escape"),
                                                  names[tax_id].encode("string_escape"),
                                                  "f",
                                                  "0", "0", "0", "0") )

            self.execute( statement )

        self.mDbHandle.commit()
        
        self.calculateTree()
    
    ##--------------------------------------------------------------------
    def getTaxIdForName( self, name, rank = None):
        statement = "SELECT tax_id FROM %s WHERE scientific_name LIKE '%%%s%%'" % (self.name, name)
        if rank: statement += " AND rank = '%s'" % rank
        result = self.getCursor(statement).fetchone()
        if not result:
            return None
        return result[0]

    ##--------------------------------------------------------------------
    def getNameMap( self ):
        """return a map of identifiers to tax_ids."""

        statement = "SELECT scientific_name, tax_id FROM %s" % self.name

        result= {}
        for scientific_name, tax_id in self.getCursor(statement).fetchall():
            result[scientific_name] = tax_id

        return result

    ##--------------------------------------------------------------------    
    def getDomainForNid( self, nid, table_name_assignments = "taxonomy_assignments" ):

        tax_id = self.getCursor( "SELECT tax_id FROM %s WHERE nid = %i" % (table_name_assignments, nid)).fetchone()
        if not tax_id: return None
        tax_id = tax_id[0]

        node_id = self.getNodeId( tax_id )
        if not node_id: return None
        
        return self.GetDomainForNodeId( node_id )

    ##--------------------------------------------------------------------
    def getNameForTaxId( self, tax_id ):
        statement = "SELECT scientific_name FROM %s WHERE tax_id = %i" % (self.name, tax_id)
        return self.getCursor(statement).fetchone()[0]
    
    ##--------------------------------------------------------------------
    def guessTaxIdForName( self, name):
        """do an approximate match for finding a certain taxonomic
        id for a given string.
        """

        words = string.split(name, " ")

        ## retrieve species
        statement = "SELECT tax_id, scientific_name FROM %s WHERE scientific_name LIKE E'%s %%'" % (self.name, words[0].encode("string_escape"))
        
        result = self.getCursor(statement).fetchall()
        
        if not result:
            return None

        if len(result) == 1:
            return result[0][0]
        
        ## match other words in strain specifier
        matches = [0] * len(result)
        for word in words[1:]:
            w = re.escape(word)
            for x in range(len(result)):
                if re.search(w, result[x][1]):
                    matches[x] += 1

        best  = 0
        bestx = 0
        for x in range(len(matches)):
            if best < matches[x]:
                best = matches[x]
                bestx = x

        ## if no further words match, return tax_id with complete match
        if bestx == 0:
            statement = "SELECT tax_id FROM %s WHERE scientific_name = E'%s'" % (self.name, words[0].encode("string_escape"))
            
            result = self.getCursor(statement).fetchone()

            if not result:
                return None
            else:
                return result[0]
        else:
            return result[bestx][0]


##------------------------------------------------------------------------
if __name__ == "__main__":

    parser = optparse.OptionParser( version = "%prog version: $Id$", usage = USAGE)

    parser.add_option("-a", "--action", dest="action", type="choice",
                      choices=("create", "map-species", "build-pairsdb" ),
                      help="action to perform. create: create taxonomy database, build-pairsdb: create table dump for pairsdb, map-species: map a list of species identifiers to taxonomic ids from stdin."  )

    parser.add_option( "--input-filename-names", dest="input_filename_names", type="string" ,
                       help="INPUT filename of taxonomy assignments: path to names.dmp file.")

    parser.add_option( "--input-filename-division", dest="input_filename_divisions", type="string" ,
                       help="INPUT filename of taxonomy assignments: path to division.dmp file.")

    parser.add_option( "--input-filename-nodes", dest="input_filename_nodes", type="string" ,
                       help="INPUT filename of taxonomy assignments: path to nodes.dmp file.")

    parser.set_defaults(
        input_filename_names = "names.dmp",
        input_filename_divisions = "division.dmp" ,
        input_filename_nodes = "nodes.dmp",
        action = None,
        table_name = "taxonomy",
        filename_species = None,
        )
        
    (options, args) = Experiment.Start( parser, add_psql_options = True )
    
    # cmd = "wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
    # dbhandle = pgdb.connect( options.connection )
    # taxonomy = TaxonomyPostgres( dbhandle, options.table_name )
    taxonomy = Taxonomy()
    
    if options.action == "create":
        if options.loglevel >= 1:
            options.stdlog.write( "# creating taxonomy table %s\n" % options.table_name )
            options.stdlog.flush()
            
        taxonomy.loadFromFiles( open(options.input_filename_nodes),
                                open(options.input_filename_divisions),
                                open(options.input_filename_names) )

    elif options.action == "build-pairsdb":
        if options.loglevel >= 1:
            options.stdlog.write( "# creating taxonomy table for pairsdb %s\n" % options.table_name )
            options.stdlog.flush()
            
        taxonomy.loadFromFiles( open(options.input_filename_nodes),
                                open(options.input_filename_divisions),
                                open(options.input_filename_names) )

        taxonomy.writeToFile( options.stdout )

    elif options.action == "map-species":

        options.stdout.write( "species\ttax_id\n" )

        missed = []
        ninput, noutput = 0, 0

        name_map = taxonomy.getNameMap()
        
        for line in sys.stdin:

            if line[0] == "#": continue
            ninput += 1
            species_name = line[:-1].split("\t")[0]

            if species_name in name_map:
                tax_id = name_map[species_name]
            else:
                tax_id = taxonomy.guessTaxIdForName( species_name )

            if tax_id:
                options.stdout.write( "\t".join(map(str, (species_name, tax_id))) + "\n" )
                options.stdout.flush()
                noutput += 1
            else:
                missed.append( species_name )
                if options.loglevel >= 2:
                    options.stdlog.write( "# missed: %s\n" % species_name )
                    options.stdlog.flush()
                    
        if options.loglevel >= 1:
            options.stdlog.write( "# ninput=%i, noutput=%i, nskipped=%i\n" % (ninput, noutput, len(missed) ) )
                    
    Experiment.Stop()
