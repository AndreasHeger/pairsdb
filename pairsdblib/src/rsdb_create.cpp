//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: rsdb_create.cpp,v 1.2 2003/01/06 14:32:06 aheger Exp $
//--------------------------------------------------------------------------------    

/**
 creating rsdb 

 input file format:
 query_nid
 sbjct_nid
 evalue (usually ignored)
 query_from
 query_to
 query_ali
 sbjct_from
 sbjct_to
 sbjct_ali
 score
 pide
 
 output file format:
 rep_nid
 rep_from
 rep_to
 rep_ali
 mem_nid
 mem_from
 mem_ali
 score
 new pide

 21.3.2001	Order alignments according to pide for each nid. This ensures, that
 if a sequence aligns multiply to another due to repeats, the best one
 is chosen first.

 10.12.2002   use index file instead of database

 27.10.2003	check for internal gaps. If there is an insert in the sbjct of more
 than the miminum domain size, then keep it.

 @author Andreas Heger
 @version $Id: rsdb_create.cpp,v 1.2 2003/01/06 14:32:06 aheger Exp $

 */

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include <map>
#include <set>
#include <vector>
#include <functional>
#include <iterator>
#include <algorithm>
#include <alignlib.h>
#include <alignlib/AlignlibIndex.h>

#include "pairsdblib.h"
#include "Connection.h"
#include "Row.h"
#include "Query.h"
#include "HelpersAlignandum.h"
#include "SQLQueries.h"

#include <fcntl.h>
#include <cstdio>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
typedef po::variables_map Options;

using namespace std;
using namespace pairsdblib;
using namespace alignlib;

#define SEPARATOR '\t'
#define MAX_LINE_LENGTH 10000

typedef std::set<Nid> NidSet;
typedef std::vector<int> LevelVector;

//----------------------------------------
// abstract retrieval of neighbors
struct NeighborGetter
{
  NeighborGetter() {};
  virtual ~NeighborGetter() {};
  virtual HNeighborLinkVector operator()( Nid & nid ) {};
};

struct NeighborGetterFromFile : public NeighborGetter
{
  typedef alignlib::Index< Nid, alignlib::RecorderTable< Nid > > Index;
  
  NeighborGetterFromFile( const Options & options ) : NeighborGetter()
  {
    debug_func_cerr( 5 );
    
    int loglevel = options["verbose"].as<int>();

    if (loglevel >= 1)
      std::cout << "## using neighbour lists from "
		<< options["filename-neighbors"].as<std::string>() << endl;
    
    mFile = alignlib::openFileForRead( options["filename-neighbors"].as<std::string>() );
    mIndex.setData( mFile );

    if (options.count("filename-index"))
      {
	// open links file: read indices
	if (loglevel >=1)
	  std::cout << "## reading index file from " << options["filename-index"].as<std::string>() << std::endl;
	    
	FILE * index_file = alignlib::openFileForRead( options["filename-neighbors"].as<std::string>() );
	mIndex.load( index_file );
	fclose( index_file );
      }
    else
      {
	if (loglevel >=1)
	  std::cout << "## building index..." << std::endl;
	mIndex.create();
      }
  }

  virtual ~NeighborGetterFromFile()
  {
    fclose( mFile );
  }
  
  virtual HNeighborLinkVector operator()(Nid & nid )
  {
    debug_func_cerr( 5 );
    
    HNeighborLinkVector result(new NeighborLinkVector() );
    
    // iterate over all neighbors in graph
    if (!mIndex.hasToken( nid ))
	return result;
	    
    mIndex.goTo( nid );

    char buffer[MAX_LINE_LENGTH+1];
    while (!feof(mFile))
      {
		
	long query_nid, sbjct_nid;
	Position query_from, query_to, sbjct_from, sbjct_to;
	char query_ali[MAX_LINE_LENGTH+1];
	char sbjct_ali[MAX_LINE_LENGTH+1];
	int score;
	float pid;
	float evalue;
		
	fscanf(mFile, "%ld\t%ld\t%f\t%i\t%i\t%s\t%i\t%i\t%s\t%i\t%f",
	       &query_nid, &sbjct_nid, &evalue, &query_from,
	       &query_to, query_ali, &sbjct_from, &sbjct_to,
	       sbjct_ali, &score, &pid);

	fgets(buffer, MAX_LINE_LENGTH, mFile);
		
	// if the sequence contains no links, the position is
	// undefined and query_nid will not be *it_reps
	if (query_nid != nid)
	  break;

	result->push_back( NeighborLink( query_nid,
					 sbjct_nid,
					 evalue,
					 query_from - OFFSET_CORRECTION,
					 query_ali,
					 sbjct_from - OFFSET_CORRECTION,
					 sbjct_ali,
					 score,
					 pid ) );

      }

    return result;
  }
  
  Index mIndex;
  FILE * mFile;
};

struct NeighborGetterFromDatabase : public NeighborGetter
{
  NeighborGetterFromDatabase( const Options & options,
			      HConnection & connection ) :
    NeighborGetter(),
    mConnection(connection),
    mTableName( options["tablename-neighbors"].as<std::string>() )
  {
    debug_func_cerr( 5 );
  }

  virtual ~NeighborGetterFromDatabase()
  {
  }
  
  virtual HNeighborLinkVector operator()(Nid & nid )
  {
    debug_func_cerr( 5 );
    return SQL_GetNeighbourLinks( (*mConnection), nid, mTableName );
  }
  std::string mTableName;
  HConnection mConnection;
};





//-------------------------------------> parameter parsing <----------------------------------
void parseArguments( po::variables_map & vm,
		     int argc, char *argv[])
{

    try 
    {
      // note: default_length option needs to be specified due
      // to a bug in boost 1.35 and 1.36 - otherwise undefined
      // value at the linking stage
      #define DEFAULT_LINE_LENGTH 100
      po::options_description generic("Generic options", DEFAULT_LINE_LENGTH);
        
        generic.add_options()                              
	  ("help,h", "produce help message.")            
	  ("version", "print Version.")
	  ("verbose,v", po::value<int>()->default_value(1),        
	   "enable verbosity (optionally specify level)");

	po::options_description specific("Specific options", DEFAULT_LINE_LENGTH);
        
        specific.add_options()                              
            ("min-domain-size,d", po::value<int>()->default_value(40),
                  "the minimum domain size.")
            ("filename-neighbors,n", po::value< std::string >(),
                  "filename with neighbors - if not given, neigbors will be retrieved from the database.")
            ("filename-index,i", po::value< std::string >(), 
            	  "filename with index of neighbor file. If not given, the index will be created on-the-fly.")
            ("filename-nids,q", po::value< std::string >(),
            		"filename with nids to analyse (+ sequence lengths). If not given, nids will be retrieved from the database.")
	  ("prefix-rsdb,p", po::value< std::string >()->default_value("rsdb"), 
            		"prefix for rsdb output files. The level will be appended." )
	  ("prefix-groupies,o", po::value< std::string >()->default_value("groupies"), 
	     "prefix for groupies output files. The level will be appended." )
	    ("report-step,r", po::value<int>()->default_value(10000),        
	     "report step")
	    ("check-insertions,r", po::value<bool>()->default_value(false),        
	     "check for insertions in sbjct")
	    ("use-evalue,e", po::value<bool>()->default_value(false),        
	     "use the evalue by filtering")
	    ("level,l", po::value< vector<int> >(),
	     "rsdb level to create. Multiple level can be given.")
	  ("tablename-neighbors,n", po::value< std::string >()->default_value( "pairsdb_90x90" ),
	   "tablename with neighbors.")
	  ("tablename-nids,n", po::value< std::string >()->default_value( "nrdb90" ),
	   "tablename with nids to use as a starting set.")
	  ;
	
	      po::options_description database("Database options", DEFAULT_LINE_LENGTH);

	database.add_options()
	  ("database,D", po::value< std::string >()->default_value("pairsdb"), 
	   "database")
	  ("user,U", po::value< std::string >()->default_value("test"), 
	   "user")
	  ("password,A", po::value< std::string >()->default_value(""), 
	   "password")
	  ("port,P", po::value< int >()->default_value(3306), 
	   "port")
	  ("host,H", po::value< std::string >()->default_value("localhost"), 
	   "host")
	  ("socket,S", po::value< std::string >()->default_value("/tmp/mysql.sock"), 
	   "socket")
	  ;

	// combine all options
	po::options_description cmdline_options(DEFAULT_LINE_LENGTH);
	cmdline_options.add(generic).add(specific).add(database);
	
	po::positional_options_description p;
	p.add("level", -1);
	
        po::store(po::command_line_parser(argc, argv).
                  options(cmdline_options).positional(p).run(), vm);
        po::notify(vm);

        if (vm.count("help")) 
        {
	  std::cout << "Usage: options_description [options]\n";
	  std::cout << cmdline_options;
	  exit(EXIT_SUCCESS);
        }

	if (vm.count("version")) 
        {
	  std::cout << "Version $Id$\n";
	  exit(EXIT_SUCCESS);
        }

    }
    catch(exception& e)
    {
        cout << e.what() << "\n";
        exit(EXIT_FAILURE);
	
    }
    
    return;
}

//---------------------------> end of parameter parsing <--------------------------------------------

template<class T>
ostream& operator<<(ostream& os, const vector<T>& v)
{

  copy(v.begin(), v.end(), ostream_iterator<T>(cout, " "));

  return os;

}

HConnection getDatabaseConnection( const po::variables_map & options )
{
  debug_func_cerr( 5 );
  
  HConnection connection(new Connection( 
					options["host"].as<std::string>(),
					options["user"].as<std::string>(),
					options["password"].as<std::string>(),
					options["port"].as<int>(),			
					options["socket"].as<std::string>() ) );
	
  connection->Connect( options["database"].as<std::string>() );
  return connection;		
}

// sort both by length an nid to produce reproducible results
class ComparatorLength : public std::binary_function<Nid,Nid,bool>
{
public:
  ComparatorLength( MapNid2Length & map ) : mMap(map) {}
  bool operator()(const Nid & a, const Nid & b) const {
	  if (mMap[a] == mMap[b])
		  return a < b;
	  else
		  return mMap[a] > mMap[b]; 
		 }
private:
  MapNid2Length & mMap;
};

/*--------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{

  Options options;

  parseArguments(options, argc, argv);
  
  // extract options
  int loglevel = options["verbose"].as<int>();
  
  bool param_use_evalue = options["use-evalue"].as<bool>();
  bool param_check_insertions = options["check-insertions"].as<bool>();
  int param_report_step = options["report-step"].as<int>();
  int param_min_domain_length = options["min-domain-size"].as<int>();
  std::vector<int> param_levels = options["level"].as< std::vector<int> >();
  
  //---------------------------------------------------------------
  // do we need database connection?
  HConnection connection;
  if (!options.count("filename-neighbors") || !options.count("filename-nids") )
    {
      connection = getDatabaseConnection( options );
    }
  
  //------------------------------------------------------------------
  // how to get neighbors: from file or from database
  NeighborGetter * get_neighbors = NULL;
  if (options.count("filename-neighbors"))
    {
      get_neighbors = new NeighborGetterFromFile( options );
      }
    else
      {
	get_neighbors = new NeighborGetterFromDatabase( options, connection );
      }

  //----------------------------------------------------------
  // read list of nids to work with. Also build map of nid to length
  HMapNid2Length map_nid2length;
  HNidVector nids;

  if (options.count("filename-nids"))
    {
      map_nid2length = HMapNid2Length(new MapNid2Length);
      nids = HNidVector(new NidVector);
      
      //------------------------------------------------------------------
      // open links file: read indices
      if (loglevel >=1)
	std::cout << "## reading nids and sequence lengths from "
		  << options["filename-nids"].as<std::string>() << endl;
  
      {
	std::ifstream fin(options["filename-nids"].as<std::string>().c_str());
	  
	while (!fin.eof())
	  {
	    Nid nid;
	    Length i;

	    fin >> nid >> i;
	    (*map_nid2length)[nid] = i;
	    (*nids).push_back(nid);
	  }
	fin.close();
      }

      if (loglevel >=1)
	std::cout << "## read nids: " << map_nid2length->size() << std::endl;
    }
  else
    {
      map_nid2length = SQL_GetLengthMap( *connection );
      nids = SQL_GetNids( *connection, options["tablename-nids"].as<std::string>() );
    }

  //----------------------------------------------------------
  // sort nids by length
  std::sort( nids->begin(), nids->end(), ComparatorLength( *map_nid2length ) );
    //----------------------------------------------------------
    // select nids to cluster
    // open new file each round

    for (unsigned int x = 0; x < param_levels.size(); x++)
      {

	unsigned int cutoff = param_levels[x];
	std::ostringstream groupies_file;
	groupies_file << options["prefix-groupies"].as<std::string>() << cutoff << std::ends;
	ofstream outfile_groupies(groupies_file.str().c_str());

	if (loglevel >= 1)
	  std::cout << "## saving alignments on level " << cutoff << " in file "
		    << groupies_file.str() << std::endl;

	// build sets of nids selected as representatives (empty at the beginning)
	// and nids available as representatives (full at the beginning)
	NidSet selected_reps;
	NidSet putative_reps;
	{
	  NidVector::iterator it(nids->begin()), end(nids->end());
	  for (; it!=end; ++it)
	    putative_reps.insert( *it);
	}

	if (loglevel >= 1)
	  std::cout << "## processing " << cutoff << " with " << nids->size()
		    << " entries " << std::endl;

	NidVector::iterator it_reps(nids->begin()), end_reps(nids->end());
	
	unsigned int iteration = 0;
	unsigned int eliminated = 0;

	// ------------------------------------------
	// ------------------------------------------
	// ------------------------------------------
	// Loop over representatives starts
	// ------------------------------------------
	// ------------------------------------------
	// ------------------------------------------
	Length last_length = std::numeric_limits<Length>::max();
	
	for (; it_reps!=end_reps; ++it_reps)
	  {
	    
	    iteration++;
	    Nid rep_nid = *it_reps;
	    
	    Length query_length = (*map_nid2length)[rep_nid];
	    assert( last_length >= query_length );
	    last_length = query_length;
	    
	    if (loglevel >= 1)
	      if (!(iteration % param_report_step))
		std::cout << "level=" << cutoff << " iteration="
			  << iteration << " query_nid=" << rep_nid
			  << " length=" << query_length << " eliminated="
			  << eliminated << std::endl;

	    if (loglevel >= 2)
	      std::cout << "## working on " << rep_nid << " length=" << query_length << std::endl;

	    // skip if not in putative reps
	    if (putative_reps.find(rep_nid) == putative_reps.end())
	      {
		if (loglevel >= 2)
		  std::cout << " eliminated" << endl;
		continue;
	      }

	    // add new as a new rep and remove from putative reps
	    selected_reps.insert(rep_nid);
	    putative_reps.erase(rep_nid);
	    
	    unsigned int nneighbors = 0;
	    unsigned int neliminated = 0;
	    
	    // ------------------------------------------
	    // ------------------------------------------
	    // ------------------------------------------
	    // Loop over neighbors starts
	    // ------------------------------------------
	    // ------------------------------------------
	    // ------------------------------------------
	    
	    HNeighborLinkVector neighbors = (*get_neighbors)( rep_nid );

	    if (loglevel >= 2)
	      std::cout << " found " << neighbors->size() << " neighbors "
			<< std::endl;

	    NeighborLinkVector::iterator it(neighbors->begin()), end(neighbors->end());
	    
	    for (; it!=end; ++it)
	      {
		// next if already selected rep
		if (selected_reps.find(it->mSbjctNid) != selected_reps.end())
		  continue;
		
		// next if not a putative rep
		if (putative_reps.find(it->mSbjctNid) == putative_reps.end())
		  continue;

		assert( map_nid2length->find( it->mSbjctNid ) != map_nid2length->end() );
		Length sbjct_length = (*map_nid2length)[it->mSbjctNid];

		Pid pid = it->mPid;
		
		// sbjct has been found in set putative reps, check, if we should keep it.
		if (loglevel >= 3)
		  std::cout << "## processing query_nid=" << it->mQueryNid << "("
			    << it->mMapQuery2Sbjct.mRowFrom << "-" << it->mMapQuery2Sbjct.mRowTo << ")"
			    << " sbjct_nid=" << it->mSbjctNid << "(" << it->mMapQuery2Sbjct.mColFrom
			    << "-" << it->mMapQuery2Sbjct.mColTo << ")" << " query_length="
			    << query_length << " sbjct_length=" << sbjct_length
			    << " pid= " << it->mPid << std::endl;
		
		if (loglevel >= 4)
		  {
		    std::cout << "## query_ali=" << it->mMapQuery2Sbjct.mRowAlignment << std::endl;
		    std::cout << "## sbjct_ali=" << it->mMapQuery2Sbjct.mColAlignment << std::endl;
		  }

		// calculate the length cutoff. If the shorter sequence is only aligned by residues 
		// less than max(length -40; length * 0.9),
		// do not look at the pid, because this sequence will not be a groupie.
		Length reference_length;
		Length aligned_region;

		// choose the smaller length sequence as reference length. The valid length for calculating 
		// the minimum length is the number
		// of aligned residues in the shorter sequence, i.e. not the overlap length, that is
		// used for calculating the percent identity.

		if (sbjct_length < query_length)
		  {
		    aligned_region = it->mMapQuery2Sbjct.mColTo - it->mMapQuery2Sbjct.mColFrom ;
		    reference_length = sbjct_length;
		  }
		else
		  {
		    aligned_region = it->mMapQuery2Sbjct.mRowTo - it->mMapQuery2Sbjct.mRowFrom;
		    reference_length = query_length;
		  }

		if ( (aligned_region < (reference_length - param_min_domain_length)) ||
		     (aligned_region < ((double)reference_length * 0.9)))
		  continue;
		
		// calculate new alignment length
		Length ali_length = 0;
		Length increment = 0;
		Length internal_gaps = 0;
		
		// use query_ali, so that we can check for large insertions in sbjct.
		// calculate length of alignment, parse a string like '+10-10+10-10+10' = 50
		std::istringstream is_row(it->mMapQuery2Sbjct.mRowAlignment);
		while (is_row >> increment)
		  if (increment > 0)
		    ali_length += increment;
		  else
		    {
		      internal_gaps -= increment;
		      ali_length -= increment;
		    }
		
		if (loglevel >= 4)
		  {
		    std::cout << "## igaps=" << internal_gaps << " lali="
			      << ali_length << std::endl;
		  }
		
		// this is a change of algorithm
		// checkpoint: do not remove alignment, if number of internal gaps is to large
		if (param_check_insertions && internal_gaps
		    > param_min_domain_length)
		  continue;
		
		if (param_use_evalue)
		  {
		    pid = -it->mEValue;
		  }
		else
		  {
		    
		    // calculate new percent identity. The percent identity is based on the total residues in the 
		    // shorter sequence. By this, residues not part of the alignment of the shorter sequence are counted 
		    // as dissimilar.
		    
		    pid *= (float)ali_length / (float)reference_length;
		    
		    if (loglevel >= 3)
		      std::cout << "## new pid= " << pid << " ref_l= "
				<< reference_length << " ali_l= " << ali_length
				<< " igaps=" << internal_gaps << " ali="
				<< it->mMapQuery2Sbjct.mRowAlignment << endl;
		  }
		
		// if percent identity of neighbour is larger than cutoff, eliminate
		// neighbour from list of putative representatives and output association
		if (pid > cutoff)
		  {
		    putative_reps.erase(it->mSbjctNid);
		    outfile_groupies << it->mQueryNid << SEPARATOR
				     << it->mMapQuery2Sbjct.mRowFrom << SEPARATOR
				     << it->mMapQuery2Sbjct.mRowTo << SEPARATOR
				     << it->mMapQuery2Sbjct.mRowAlignment << SEPARATOR
				     << it->mSbjctNid << SEPARATOR
				     << it->mMapQuery2Sbjct.mColFrom << SEPARATOR
				     << it->mMapQuery2Sbjct.mColTo << SEPARATOR
				     << it->mMapQuery2Sbjct.mColAlignment << SEPARATOR
				     << it->mScore << SEPARATOR
				     << it->mPid << SEPARATOR << endl;
		    ++neliminated;
		  }
		
	      } // end of loop over all neighbors
	    
	    eliminated += neliminated;

	    if (loglevel >= 1)
	      std::cout << "nid=" << rep_nid
			<< " len=" << query_length
			<< " nneighbors=" << neighbors->size()
			<< " neliminated=" << neliminated
			<< std::endl;
	    
	    // if (iteration > 10) break;

	  } // end of loop over all nids

	// write set of representatives
	{
	  std::ostringstream rsdb_file;
	  rsdb_file << options["prefix-rsdb"].as<std::string>() << cutoff << ends;
	  ofstream outfile_rsdb(rsdb_file.str().c_str());
	  
	  if (loglevel >= 1)
	    std::cout << "## saving << " << selected_reps.size()
		      << " representatives on level " << cutoff << " in "
		      << rsdb_file.str() << endl;
	  
	  std::copy(selected_reps.begin(), selected_reps.end(),
		    std::ostream_iterator< Nid>(outfile_rsdb, "\n"));
	  outfile_rsdb.close();
	}
	
	putative_reps.clear();
	selected_reps.clear();
				
	outfile_groupies.close();

      } // end of loop over cutoffs
    
    delete get_neighbors;

  if (!options.count("filename-neighbors") || !options.count("filename-nids") )
      connection->Disconnect();
    
}

