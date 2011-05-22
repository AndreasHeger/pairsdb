//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: get_conservation.cpp,v 1.1.1.1 2002/07/03 11:17:13 heger Exp $
//--------------------------------------------------------------------------------    

/** Retrieve a multiple alignment from Picasso and write to stdout the 
    consensus string. Then put out for each member of the alignment on
    stderr:

    nid query_column sbjct_column 
    
    for all positions in the alignment, where the aligned residue is 
    equal to the consensus string.

*/

#include <iostream>
#include <iomanip>
#include <sstream>

#include <alignlib.h>

#include "pairsdblib.h"
#include "Connection.h"
#include "Row.h"
#include "Query.h"
#include "HelpersMultipleAlignment.h"

using namespace std;
using namespace pairsdblib;
using namespace alignlib; 

/* Global parameters that are can be set by command line arguments */
static Nid param_nid = 0;
static Filter  param_filter_lower = DEFAULT_LOWER_FILTER;
static Filter param_filter_upper = DEFAULT_UPPER_FILTER;
static unsigned int param_max_lines = DEFAULT_MAX_LINES;
static EValue param_max_evalue = DEFAULT_MAX_EVALUE;
static double param_cutoff = 0.51;

static std::string param_host("localhost");
static std::string param_user("test");
static std::string param_password("");
static int param_port = 3306;
static std::string param_database("pairsdb");

static int param_loglevel = 1;

static std::string param_table_name("pairsdb_90x90");

#include <unistd.h>
const char * my_progname = "get_conservation";
const char * SYSTEM_TYPE = "..";
const char * MACHINE_TYPE = "..";


static void print_version() {
  cout << my_progname << " Version ... for ... at ..." << endl;
}    

static void usage()
{
  print_version();
  cout << "Usage: " << my_progname << " [OPTIONS] nid" << endl;
  cout << "calculates conservation for pairsdb_90x90 multiple alignment" << endl;
  cout << "-D		database to use [default=pairsdb]." << endl;
  cout << "-H		host[localhost]." << endl;
  cout << "-U		username[test]." << endl;
  cout << "-P		port[3306]." << endl;
  cout << "-A		password[]." << endl;
  cout << "-c		conservation cutoff (conserved, if max frequency > cutoff)." << endl;
  cout << "-e		maximum evalue." << endl;
  cout << "-f		minimum filter (in rsdb)." << endl;
  cout << "-g		maximum filter (in rsdb)." << endl;
  cout << "-l		max num lines (of multiple alignment)." << endl;
  cout << "-t		tablename [pairsdb_90x90]." << std::endl;
  cout << "-v		loglevel " << endl;
  cout << "-V		print version and exit." << endl;
}

void ParseArguments (int argc, char *argv[]) {

  int c;  

  extern char * optarg;

  while ((c=getopt(argc, argv, "V?D:c:e:f:g:l:v:t:H:U:P:A:")) != EOF) {
    switch(c) {
    case 'D':
      param_database = optarg; break;
    case 'H':
      param_host = optarg; break;
    case 'U':
      param_user = optarg; break;
    case 'A':
      param_password = optarg; break;
    case 'P':
      param_port = atoi(optarg); break;
    case 'c':
	param_cutoff = atof(optarg); break;
    case 'e':
      param_max_evalue = atof(optarg); break;
    case 'f':
      param_filter_lower = atoi(optarg); break;
    case 'g':
      param_filter_upper = atoi(optarg); break;
    case 'l':
      param_max_lines = atoi(optarg); break;
    case 'v':
      param_loglevel = atoi(optarg); break;
    case 't':
      param_table_name = optarg; break;
    case 'V':
      print_version(); exit(0);
    case '?':
      usage(); exit(0);
    }
  }

  // set pointers to end of options
  (argc)-=optind;
  (argv)+=optind;

  if (argc == 0) {
    usage();
    exit(-1);
  }     

  param_nid = atol(argv[0]); 
  
}


/*--------------------------------------------------------------------------*/
int main (int argc, char *argv[]) {

    ParseArguments( argc, argv );

    pairsdblib::Connection connection( param_host.c_str(), 
				       param_user.c_str(), 
				       param_password.c_str(), 
				       param_port );
  
    connection.Connect( param_database.c_str());
    
    if (!connection.isConnected()) {
      cerr << "Could not connect to database " << param_database << std::endl;
      exit(EXIT_FAILURE);
    }

    HMultipleAlignment m( makeMultipleAlignmentNeighbours( connection, 
							   param_nid, 
							   param_table_name.c_str(),
							   param_max_lines,
							   0,
							   param_filter_upper,
							   param_filter_lower,
							   param_max_evalue ) );
  int mali_width = m->get();
  if (! m->isEmpty())
    {

      std::string consensus = calculateConservation( m, param_cutoff);
    
      cout << consensus << endl;

    Position length = m->getLength();
    char gap_char = alignlib::getDefaultTranslator()->getGapChar();
    
    // retrieve the neighbours. They have to be in the same order as in makeMultipleAlignmentNeighbours90
    Query query_sbjcts( connection );
    {
      ostringstream query_str;
      
      query_str << "SELECT p.query_from, p.query_ali, p.sbjct_nid, p.sbjct_from, p.sbjct_ali "
		<< " FROM " << param_table_name << " AS p, nrdb n "
		<< " WHERE p.query_nid = " << param_nid 
		<< " AND n.nid = p.sbjct_nid ";
      
      if (param_filter_lower != DEFAULT_LOWER_FILTER)
	query_str << " AND nrdb.filter >= " << param_filter_lower;
      if (param_filter_upper != DEFAULT_UPPER_FILTER)
	query_str << " AND nrdb.filter <= " << param_filter_upper;
      if (param_max_evalue != DEFAULT_MAX_EVALUE)
	query_str << " AND p.evalue <= " << param_max_evalue;

      query_str << " ORDER BY p.evalue";

      if (param_max_lines != DEFAULT_MAX_LINES)
	query_str << " LIMIT " << param_max_lines;

      query_str << ends;
      
      query_sbjcts.startQuery( query_str.str());
    }
    
    typedef enum { FIELD_query_from, FIELD_query_ali, 
		   FIELD_sbjct_nid, FIELD_sbjct_from, FIELD_sbjct_ali } MY_FIELDS ;

    Row row;
    

    // there has to be an offset of 1, since the first row in the multiple alignment is the representative
    int current_row = 1;		

    assert( (mali_width-1) == query_sbjcts.getResultNumRows());
    
    // iterate over each row in the multiple alignment
    while (query_sbjcts.fetchRow( row )) {
      
      assert(current_row < mali_width);

      Nid sbjct_nid = atol(row[FIELD_sbjct_nid]);
      const std::string row_ali = row[FIELD_query_ali];
      const std::string col_ali = row[FIELD_sbjct_ali];
      
      // create alignment
      // row = query
      // col = sbjct

      alignlib::Alignment * ali = alignlib::fillAlignmentCompressed( alignlib::makeAlignmentVector(), 
								   atoi(row[FIELD_query_from]),
								   row_ali,
								   atoi(row[FIELD_sbjct_from]),
								   col_ali 
								   );

      const std::string & alignment_string = (*m)[current_row];
#ifdef DEBUG
      cerr << "-----------------------------" << endl;
      cerr << sbjct_nid << endl << alignment_string << endl << consensus << endl;
#endif
      
      for (int pos = 0; pos < length; pos++) {
#ifdef DEBUG
	cerr << pos << "\t" << alignment_string[pos] << "\t" 
	     << consensus[pos] << endl; 

#endif
	  if (gap_char != alignment_string[pos] && alignment_string[pos] == consensus[pos]) {
	    cerr << sbjct_nid << "\t" 
		 << pos + 1   << "\t" 
		 << ali->mapRowToCol( pos + 1) << "\t" 
		 << endl;
	  }
      }


      current_row ++;
      delete ali;
    }
    
    query_sbjcts.endQuery();
  } 

  delete m;

  connection.Disconnect();

}     


