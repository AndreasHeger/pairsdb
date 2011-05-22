//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: get_mali.cpp,v 1.2 2003/01/06 14:32:06 aheger Exp $
//--------------------------------------------------------------------------------    

#include <stdlib.h>
#include <stdio.h>

#include <iostream>
#include <iomanip>
#include <string>
#include <alignlib.h>

#include "pairsdblib.h"
#include "Connection.h"
#include "Row.h"
#include "Query.h"
#include "HelpersMultipleAlignment.h"

using namespace std;
using namespace pairsdblib;
using namespace alignlib; 

//-------------------------------------> parameter parsing <----------------------------------

#include <unistd.h>

const char * my_progname = "get_mali";
const char * SYSTEM_TYPE = "..";
const char * MACHINE_TYPE = "..";

static void print_version() {
  cout << my_progname << " Version ... for ... at ..." << endl;
}    

static void usage()
{
  print_version();
  cout << "Usage: " << my_progname << "[OPTIONS] nid|did [repeat_id]\n" << endl;
  cout << "-D   database to use." << endl;
  cout << "-H   host[localhost]." << endl;
  cout << "-U   username[test]." << endl;
  cout << "-P   port[3306]." << endl;
  cout << "-S   socket[3306]." << endl;
  cout << "-A   password[]." << endl;
  cout << "-r   renderer to use." << endl;
  cout << "-c   write consensus string." << endl;
  cout << "-e   maximum evalue." << endl;
  cout << "-l	maximum number of lines." << endl;
  cout << "-o	offset of alignment [lines]." << endl;
  cout << "-g	write full expanded alignment (domain alignments)." << endl;
  cout << "-i	maximum number of insertions per sequence (domain alignments)." << endl;  
  cout << "-f   apply filter: upper range." << endl;
  cout << "-m   apply filter: lower range." << endl;
  cout << "-t	table that contains alignments." << endl;
  cout << "-i   table that contains alignment information." << endl;
  cout << "-y   table type (0=pairsdb,1=groupies,2=radar,3=domains)." << endl;
  cout << "-s from,to  only print out region from from to to (counting starts at 1)." << std::endl;
  cout << "-v	loglevel " << endl;
  cout << "-V   print version and exit." << endl;
}

/* Global parameters that are can be set by command line arguments */
static int  param_renderer = 0;
static Did param_nid = NULL;
static int param_table_type = 0;
static std::string param_table_name1("");
static RepeatId param_repeat_id = 0;
static bool param_consensus = false;
static bool param_compress_alignment = true;
static int param_max_insertion_length = -1;

static int param_max_lines = DEFAULT_MAX_LINES;
static int param_offset_lines = 0;
static double param_max_evalue = DEFAULT_MAX_EVALUE;
static Filter param_filter_upper = DEFAULT_UPPER_FILTER;
static Filter param_filter_lower = DEFAULT_LOWER_FILTER;


static std::string param_host("localhost");
static std::string param_user("test");
static std::string param_password("");
static int param_port = 3306;
static std::string param_database("pairsdb");
static std::string param_socket("/var/lib/mysql/mysql.sock");

static int param_loglevel = 1;

void ParseArguments (int argc, char *argv[]) {

  int c;  

  extern char * optarg;

  while ((c=getopt(argc, argv, "gcV?D:e:f:t:i:l:m:n:o:r:y:v:H:U:P:A:S:s:")) != EOF) {
    switch(c) {
    case 'D':
      param_database = optarg; break;
    case 'H':
      param_host = optarg; break;
    case 'U':
      param_user = optarg; break;
    case 'A':
      param_password = optarg; break;
    case 'S':
      param_socket = optarg; break;
    case 'P':
      param_port = atoi(optarg); break;
    case 'c':
	param_consensus = true; break;
    case 'r':
      param_renderer = atoi(optarg); break;
    case 't':
      param_table_name1 = optarg; break;
    case 'e':
      param_max_evalue = atof(optarg); break;
    case 'y':
      param_table_type = atoi(optarg); break;
    case 'l':
      param_max_lines = atoi(optarg); break;
    case 'o':
      param_offset_lines = atoi(optarg); break;
    case 'f':
      param_filter_upper = atoi(optarg); break;
    case 'm':
      param_filter_lower = atoi(optarg); break;
    case 'g':
      param_compress_alignment = false; break;
    case 'i':
      param_max_insertion_length = atoi(optarg); break;
    case 'v':
      param_loglevel = atoi(optarg); break;
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

  param_nid = argv[0]; 
  
  if (argc > 1)
    param_repeat_id = atoi(argv[1]);

}

//---------------------------> end of parameter parsing <--------------------------------------------

/*--------------------------------------------------------------------------*/
int main (int argc, char *argv[]) {

    ParseArguments( argc, argv );

    pairsdblib::Connection connection( param_host.c_str(), 
				       param_user.c_str(), 
				       param_password.c_str(), 
				       param_port,
				       param_socket.c_str() );
  
    connection.Connect( param_database.c_str());
    
    if (!connection.isConnected()) {
      cerr << "Could not connect to database " << param_database << std::endl;
      exit(EXIT_FAILURE);
    }

  //----------------------------------------------------------  
  // create the multiple alignment

    HMultipleAlignment m;

  switch (param_table_type) {
  case 0: m = makeMultipleAlignmentNeighbours( connection, 
					       atol(param_nid), 
					       param_table_name1.c_str(), 
					       param_max_lines,
					       param_offset_lines,
					       param_filter_upper, 
					       param_filter_lower,
					       param_max_evalue); 
  break;
  case 1: m = makeMultipleAlignmentGroupies( connection, 
					     atol(param_nid), 
					     param_table_name1.c_str(),
					     param_max_lines,
					     param_offset_lines); 
  break;
  case 2: m = makeMultipleAlignmentRadar( connection, 
					  atol(param_nid), 
					  param_repeat_id, 
					  param_table_name1.c_str(),
					  param_max_lines,
					  param_offset_lines); 
  break;
  case 3: m = makeMultipleAlignmentDomains( connection, 
					    param_nid, 
					    param_table_name1.c_str(), 
					    param_max_lines,
					    param_offset_lines,
					    param_compress_alignment,
					    param_max_insertion_length); 
  break;
  default: cerr << "Unknown table type" << endl; exit(-1);
  }    

  std::string consensus( calculateConservation( m, 0.51) );

  // write consensus string at top
  if (param_consensus) 
    cout << consensus << endl;

  // write multiple alignment
  switch (param_renderer) 
  {
  	case 0:
  		std::cout << alignlib::MultipleAlignmentFormatPlain( m );
  		break;
  	case 1: 
  		std::cout << alignlib::MultipleAlignmentFormatHTML( m, alignlib::makePaletteMView());
  		break;
  }    
  connection.Disconnect();
}     





