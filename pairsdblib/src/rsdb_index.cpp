//--------------------------------------------------------------------------------
// Project adda
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id$
//--------------------------------------------------------------------------------    

#include <math.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include <fcntl.h>
#include <cstdio>
#include <cstdlib>

#include <alignlib.h>
#include <alignlib/AlignlibIndex.h>


#include "pairsdblib.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;
typedef po::variables_map Options;

using namespace std;
using namespace pairsdblib;
using namespace alignlib;

#define MAX_LINE_LENGTH 65536

/*--------------------------------------------------------------*/
#include <unistd.h>

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
	("command,c", po::value< std::string >(),
	 "commands to execute")
	("filename-neighbors,n", po::value< std::string >(),
	 "filename with neighbors")
	("filename-index,i", po::value< std::string >(), 
	 "filename with index of neighbor file. If not given, the index will be created on-the-fly.")
	;
      
      // combine all options
      po::options_description cmdline_options(DEFAULT_LINE_LENGTH);
      cmdline_options.add(generic).add(specific);

      po::positional_options_description p;
      p.add("command", 1);

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
/*--------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{

  Options options;

  parseArguments(options, argc, argv);

  std::string command = options["command"].as<std::string> ();

  // extract options
  int loglevel = options["verbose"].as<int>();
  if (options.count("filename-neighbors") == 0)
    {
      throw("no filename neighbours");
    }

  typedef alignlib::Index< Nid, alignlib::RecorderTable< Nid > > Index;
  
  FILE * infile = alignlib::openFileForRead( options["filename-neighbors"].as<std::string>() );
  Index index;

  index.setData( infile );

  if (command == std::string("build") )
    {
      if (loglevel >= 1)
	std::cout << "## building index..." << std::endl;
      index.create();
      if (loglevel >= 1)
	std::cout << "## saving index..." << std::endl;
      
      FILE * outfile = alignlib::openFileForWrite( options["filename-index"].as<std::string>() );
      index.save( outfile );
    }
  else if (command == std::string("check") )
    {
      if (loglevel >= 1)
	std::cout << "## checking index..." << std::endl;
      
      FILE * index_file = alignlib::openFileForRead( options["filename-index"].as<std::string>() );
      index.load( index_file );
      fclose( index_file );
      
      Index::MapToken2IndexIterator it(index.begin());

      for (; it != index.end(); ++it)
	{
	  Nid nid = it->first;
	  index.goTo( nid );
	  
	}
    }
  else
    {
      throw "unknown command";
    }
	
  if (loglevel >= 1)
    std::cout << "## done" << std::endl;

  exit( EXIT_SUCCESS );
  
}






