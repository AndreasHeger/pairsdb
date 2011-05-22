//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: align_groupies.cpp,v 1.4 2004/10/14 23:33:05 aheger Exp $
//--------------------------------------------------------------------------------    

/**
   align all sequences in pairsdb_100x90 without an alignment

   @author Andreas Heger
   @version $Id: align_groupies.cpp,v 1.4 2004/10/14 23:33:05 aheger Exp $

*/

#include <iostream>
#include <iomanip>
#include <sstream>

#include <alignlib.h>

#include "pairsdblib.h"
#include "Connection.h"
#include "Row.h"
#include "Query.h"

using namespace std;
using namespace pairsdblib;
using namespace alignlib;

#define SEPARATOR '\t'
#define MIN_DOMAIN_SIZE 30

//-------------------------------------> parameter parsing <----------------------------------
#include <unistd.h>

const char * my_progname = "align_groupies";
const char * SYSTEM_TYPE = "..";
const char * MACHINE_TYPE = "..";

static void print_version() {
  cout << my_progname << " Version ... for ... at ..." << endl;
}    

static void usage()
{
  print_version();
  cout << "Usage: " << my_progname << "[OPTIONS] tablename" << endl;
  cout << "-D   database to use[pairsdb]." << endl;
  cout << "-H   host[localhost]." << endl;
  cout << "-U   username[test]." << endl;
  cout << "-P   port[3306]." << endl;
  cout << "-A   password[]." << endl;
  cout << "-r	report steps[10000]." << endl;
  cout << "-o	gap opening penalty [default=-10]" << endl;
  cout << "-e	gap extension penalty [default=-2]" << endl;
  cout << "-v	loglevel " << endl;
  cout << "-V   print version and exit." << endl;
}

static std::string param_host("localhost");
static std::string param_user("test");
static std::string param_password("");
static int param_port = 3306;

static int param_report_step = 10000;

static std::string param_database("pairsdb");
static std::string param_table_name("");
static alignlib::Score param_gop = -10;
static alignlib::Score param_gep = -2;
static int param_loglevel = 1;

void ParseArguments (int argc, char *argv[]) {

  int c;  

  extern char * optarg;

  while ((c=getopt(argc, argv, "V?D:t:o:e:v:r:H:U:P:A:")) != EOF) {
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
    case 'r':
      param_report_step = atoi(optarg); break;
    case 't':
      param_table_name = optarg; break;
    case 'o':
      param_gop = atof(optarg); break;
    case 'e':
      param_gep = atof(optarg); break;
    case 'v':
      param_loglevel = atoi(optarg); break;
    case 'V':
      print_version(); exit(EXIT_SUCCESS);
    case '?':
      usage(); exit(EXIT_SUCCESS);
    }
  }

  // set pointers to end of options
  (argc)-=optind;
  (argv)+=optind;

  if (argc != 1) {
    usage();
    exit(EXIT_FAILURE);
  }     

  param_table_name=argv[0];

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

    std::ostringstream query_str; 
  
    query_str << "SELECT a.sequence, a.nid, b.sequence, b.nid "
	      << " FROM nrdb a, nrdb b, " << param_table_name << " AS c" 
	      << " WHERE  a.nid = c.rep_nid AND " 
	      << "        b.nid = c.mem_nid AND "
	      << "        c.mem_ali = ''" 
	      << ends; 

    Query query( connection );
    query.startQuery( query_str.str() );
    
    Row row; 

    if (param_loglevel >= 1) 
      std::cerr << "aligning " << query.getResultNumRows() << " pairs from " << param_table_name << endl;

    // set up the objects
    HAlignment result(makeAlignmentVector());
    alignlib::setDefaultSubstitutionMatrix( alignlib::makeSubstitutionMatrixPam30() );

   int iterations = 0;

   HAlignator alignator_gaps(makeAlignatorDPFull( ALIGNMENT_GLOBAL, param_gop, param_gep ));
   
    // The Result class has a read-only Random Access Iterator 
    while (query.fetchRow( row )) 
      {
	iterations++;
	
	if (param_loglevel >= 2)
	  if (!(iterations % param_report_step)) 
	    std::cerr << "----> at iteration " << iterations << std::endl;

	alignlib::HAlignandum seq1(alignlib::makeSequence(row[0]));
	Nid nid1 = atol(row[1]);
	alignlib::HAlignandum seq2(alignlib::makeSequence(row[2]));
	Nid nid2 = atol(row[3]);

	if (param_loglevel >= 3) {
	  cout << *seq1 << "....." << seq1->getLength() << endl;
	  cout << *seq2 << "....." << seq2->getLength() << endl;
	}

	Position lseq1 = seq1->getLength();
	// Position lseq2 = seq2->getLength();
      
	alignlib::HAlignator dottor( alignlib::makeAlignatorTuples( 3 ) );
	alignlib::HAlignment dots(alignlib::makeAlignmentMatrixUnsorted());

	dottor->align( dots, seq1, seq2 );

	if (param_loglevel >= 3) 
	  cout << "-> there are " << dots->getLength() << " dots " << endl;
	  
	// calculate maximum diagonal
	alignlib::AlignmentIterator it(dots->begin());
	alignlib::AlignmentIterator it_end(dots->end());

	std::vector<int> diagonals(lseq1*2,0);
	for (; it != it_end; ++it) {
	  Position diagonal = it->mCol - it->mRow + lseq1;
	  diagonals[diagonal] += 1;
	}

	int max_diag = 0;
	for (int x = 0; x < lseq1*2; x++)
	  {
	    if (diagonals[max_diag] < diagonals[x])
	      max_diag = x;
	  }
      
	max_diag -= lseq1;
	
	alignlib::HAlignment new_dots(alignlib::makeAlignmentMatrixRow());
	
	alignlib::copyAlignment( new_dots, dots, 0,0,0,0, 
				max_diag - MIN_DOMAIN_SIZE, 
				max_diag + MIN_DOMAIN_SIZE );
	
	if (param_loglevel >= 3) {
	  cout << "-> there are " << new_dots->getLength() << " dots after filtering" << endl;
	  cout << "-> tube used: " << max_diag - MIN_DOMAIN_SIZE << "-" <<  max_diag + MIN_DOMAIN_SIZE << endl;
	}

	if (param_loglevel >= 4 )
	  std::cout << "Dots=" << std::endl << *new_dots;
	
	alignlib::HAlignator p(alignlib::makeAlignatorPrebuilt( new_dots ));
	alignlib::HAlignator alignator(alignlib::makeAlignatorDots( p, param_gop, param_gep ));
	
	alignator->align( result, seq1, seq2 );

	if (result->getLength() > 0)
	  {

	    // fill gaps
	    alignlib::fillAlignmentGaps( result, alignator_gaps, seq1, seq2 );
	    

	    alignlib::AlignmentFormatEmissions f( result );
	    
	    cout << nid1   << SEPARATOR ;
	    cout << result->getRowFrom() << SEPARATOR << result->getRowTo() << SEPARATOR;
	    cout << f.mRowAlignment << SEPARATOR;
	    
	    cout << nid2   << SEPARATOR ;
	    cout << result->getColFrom() << SEPARATOR << result->getColTo() << SEPARATOR;
	    cout << f.mColAlignment << SEPARATOR;
	
	    cout << result->getScore() << SEPARATOR;
	    cout << calculatePercentIdentity( result, seq1, seq2) * 100;
	    cout << endl;
	    
	    if (param_loglevel >= 3)
	      {
		std::cout << *result;
		std::cout << AlignmentFormatExplicit( result, seq1, seq2 );
		std::cout << std::endl;
	    }
	  }
      } 

    query.endQuery();
    
    connection.Disconnect();
    
    return EXIT_SUCCESS;
}


