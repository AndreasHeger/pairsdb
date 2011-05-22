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

#include "adda.h"

using namespace std;

/*--------------------------------------------------------------*/
#include <unistd.h>

static std::string param_file_name_neighbours = "neighbours.in";
static std::string param_file_name_index = "index.in";

static unsigned int param_report_step = 10000;
static unsigned int param_loglevel = 0;
static bool param_create = true;

const char * my_progname = "index";
const char * SYSTEM_TYPE = "..";
const char * MACHINE_TYPE = "..";

static void print_version() {
  cout << my_progname << " Version $Id$ for  at ..." << endl;
}    

static void usage()
{
  print_version();
  cout << "Usage: " << my_progname << "[OPTIONS] links \n" << endl;
  cout << "-v #		loglevel [" << param_loglevel << "]" << endl;
  cout << "-r #		loglevel [" << param_report_step << "]" << endl;
  cout << "-n #		file of sorted neighbours [" << param_file_name_neighbours << "]" << endl;
  cout << "-f #		file of indices for neighbours [" << param_file_name_index << "]" << endl;
  cout << "-c 		check index file versus neighbours" << endl;

}

void ParseArguments (int argc, char *argv[]) {

  int c;  
  
  extern char * optarg;

  while ((c=getopt(argc, argv, "r:v:n:f:c")) != EOF) {
    switch(c) {
    case 'v':
      param_loglevel = atoi(optarg); break;
    case 'r':
      param_report_step = atoi(optarg); break;
    case 'n':
      param_file_name_neighbours = optarg; break;
    case 'f':
      param_file_name_index = optarg; break;
    case 'c':
      param_create = false; break;
    }
  }

  // set pointers to end of options
  (argc)-=optind;
  (argv)+=optind;

  if (argc != 0) {
    usage();
    exit(EXIT_FAILURE);
  }     

}


//--------------------------------------------------------------------------------
int main (int argc, char *argv[]) {

  ParseArguments( argc, argv );

  FILE * infile = fopen(param_file_name_neighbours.c_str(), "r"); 
  
  if (infile == NULL) {
    std::cerr << "could not open filename with neighbours: " << param_file_name_neighbours << std::endl;
    exit(EXIT_FAILURE);
  }

  FILE * outfile;
  if (param_create) 
    outfile = fopen(param_file_name_index.c_str(), "w"); 
  else
    outfile = fopen(param_file_name_index.c_str(), "r"); 
  
  if (outfile == NULL) {
    std::cerr << "could not open filename for indices: " << param_file_name_index << std::endl;
    exit(EXIT_FAILURE);
  }

  char buffer[MAX_LINE_LENGTH+1];
  unsigned int iteration = 0;
  Nid nid;
  FileIndex index;

  if (param_create) {

    int last_nid = 0;
  
    while(!feof(infile)) { 
      
      iteration ++;
      
      if (param_loglevel >= 2) {
	if (!(iteration % param_report_step)) {
	  std::cout << "# line=" << iteration << " last_nid=" << last_nid << std::endl;
	}
      }
      
      fgetpos( infile, &index );    
      if (fscanf(infile, "%ld", &nid) != 1) break;
      
      fgets( buffer, MAX_LINE_LENGTH, infile );

      if (last_nid != nid) {
	fwrite(&nid,sizeof(Nid),1,outfile);
	fwrite(&index,sizeof(FileIndex),1,outfile); 
	last_nid = nid;
      }
      
    }

  } else {

    if (param_loglevel >= 1) {
      std::cout << "# checking index " << param_file_name_index << std::endl;
      std::cout << "# against neighbours in " << param_file_name_neighbours << std::endl;
    }

    while(!feof(outfile)) { 

      iteration ++;
    
      fread(&nid,sizeof(Nid), 1, outfile);
      if (feof(outfile)) break;
      fread(&index, sizeof(FileIndex), 1, outfile);

      if (param_loglevel >= 2) {
	if (!(iteration % param_report_step)) {
	  std::cout << "# line=" << iteration << " nid=" << nid << std::endl;
	}
      }

      Nid check_nid;
      fsetpos( infile, &index );    
      if (fscanf(infile, "%ld", &check_nid) != 1) {
	std::cerr << "error for nid " << nid << ": incorrect file position gives" 
		  << check_nid << std::endl;
	exit(EXIT_FAILURE);
      }
    }
  }
  
  fclose(infile);  
  fclose(outfile);  

}








