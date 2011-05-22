//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: rsdb_filter.cpp,v 1.1.1.1 2002/07/03 11:17:13 heger Exp $
//--------------------------------------------------------------------------------    

/**
 filtering a list of sequences for redundancy:

 -> Criterea for filtering:
 1. full coverage (more than param_min_coverage or less than param_min_domain_length residues missing)
 2. E-Value large enough

 @author Andreas Heger
 @version $Id: rsdb_filter.cpp,v 1.1.1.1 2002/07/03 11:17:13 heger Exp $

 */

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include <alignlib.h>

#include <set>
#include <vector>

#include "pairsdblib.h"
#include "Connection.h"
#include "Row.h"
#include "Query.h"
#include "HelpersAlignandum.h"

using namespace std;
using namespace pairsdblib;
using namespace alignlib;

#define SEPARATOR '@'

//-------------------------------------> parameter parsing <----------------------------------

#include <unistd.h>

const char * my_progname = "rsdb_create";
const char * SYSTEM_TYPE = "..";
const char * MACHINE_TYPE = "..";

static void print_version()
{
	cout << my_progname << " Version ... for ... at ..." << endl;
}

static void usage()
{
	print_version();
	cout << "Usage: " << my_progname << "[OPTIONS] \n" << endl;
	cout << "-D   database to use[pairsdb]." << endl;
	cout << "-H   host[localhost]." << endl;
	cout << "-U   username[test]." << endl;
	cout << "-P   port[3306]." << endl;
	cout << "-A   password[]." << endl;
	cout << "-r	  report steps[10000]." << endl;
	cout << "-d	  minimum domain length." << endl;
	cout << "-c	  minimum coverage." << endl;
	cout << "-f   filename for rsdb output." << endl;
	cout << "-g   filename for nids output." << endl;
	cout << "-t   tablename of nids to cluster." << endl;
	cout << "-v	  loglevel " << endl;
	cout << "-V   print version and exit." << endl;
}

/* Global parameters that are can be set by command line arguments */
static std::string param_table_name = "nrdb40";

// minimum coverage
static double param_min_coverage = 0.9;

// minimum domain length
static int param_min_domain_length = 40;

// database to use
static std::string param_database = "pairsdb";

// evalue cutoff
static double param_max_evalue = -460.0;

// outputs
static std::string param_output_rsdb = "rsdb.out";
static std::string param_output_nids = "nids.out";

static std::string param_host("localhost");
static std::string param_user("test");
static std::string param_password("");
static int param_port = 3306;

static int param_report_step = 10000;
static int param_loglevel = 1;

void ParseArguments(int argc, char *argv[])
{

	int c;

	extern char * optarg;

	while ((c=getopt(argc, argv, "V?D:d:c:e:f:g:t:r:H:U:P:A:v:")) != EOF)
	{
		switch (c)
		{
		case 'D':
			param_database = optarg;
			break;
		case 'H':
			param_host = optarg;
			break;
		case 'U':
			param_user = optarg;
			break;
		case 'A':
			param_password = optarg;
			break;
		case 'P':
			param_port = atoi(optarg);
			break;
		case 'r':
			param_report_step = atoi(optarg);
			break;
		case 'f':
			param_output_rsdb = optarg;
			break;
		case 'g':
			param_output_nids = optarg;
			break;
		case 't':
			param_table_name = optarg;
			break;
		case 'd':
			param_min_domain_length = atoi(optarg);
			break;
		case 'c':
			param_min_coverage = atof(optarg);
			break;
		case 'e':
			param_max_evalue = atof(optarg);
			break;
		case 'v':
			param_loglevel = atoi(optarg);
			break;
		case 'V':
			print_version();
			exit(0);
		case '?':
			usage();
			exit(0);
		}
	}

	// set pointers to end of options
	(argc)-=optind;
	(argv)+=optind;

	if (argc > 0)
	{
		usage();
		exit(-1);
	}

}

//---------------------------> end of parameter parsing <--------------------------------------------

/*--------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{

	ParseArguments(argc, argv);

	pairsdblib::Connection connection(param_host.c_str(), param_user.c_str(),
			param_password.c_str(), param_port);

	connection.Connect(param_database.c_str());

	if (!connection.isConnected())
	{
		cerr << "Could not connect to database " << param_database << std::endl;
		exit(EXIT_FAILURE);
	}

	ofstream out_rsdb(param_output_rsdb.c_str() );
	ofstream out_nids(param_output_nids.c_str() );

	//----------------------------------------------------------
	// select nids to cluster
	Query query_nids(connection);
	{

		std::ostringstream query_str;

		query_str << "SELECT n.nid, nrdb.length FROM "
				<< param_table_name.c_str() << " AS n "
				<< ", nrdb WHERE nrdb.nid = n.nid "
				<< " ORDER BY length DESC, nid ASC " << ends;

		query_nids.startQuery(query_str.str() );
	}

	cout << "processing " << param_max_evalue << " with "
			<< query_nids.getResultNumRows() << " entries from "
			<< param_table_name.c_str() << endl;

	// use hash implementation, if it becomes part of the standard
	typedef std::set<Nid> NID_SET;

	NID_SET groupies_nids; // set of assigned groupies, gets filled
	NID_SET reps_nids; // set of current representatives, gets filled

	Row row_nids;
	long iterations = 0;

	while (query_nids.fetchRow(row_nids))
	{

		iterations ++;

		Nid query_nid = atol(row_nids[0]);
		Nid query_length = atol(row_nids[1]);

		if (param_loglevel >= 3)
			cout << "Working on " << query_nid << " " << query_length << endl;

		if (param_loglevel >= 1)
			if (!(iterations % param_report_step))
				cout << "-> at iteration " << iterations
						<< " working on: query_nid=" << query_nid << " length="
						<< query_length << endl;

		if (groupies_nids.find(query_nid) == groupies_nids.end())
		{

			// add query to hash of representatives
			reps_nids.insert(query_nid);

			out_nids << query_nid << endl;

			/* check if neighbour to previous sequence
			 and go through neighbours and eliminate all those above cutoff. Only select those neighbours, which
			 had been a representative on the previous level. Sort the neighbours by length and then query_nid, 
			 so that if a new rep will be found, it is on top of the list. This means, there will be no need 
			 to remap.
			 */

			Query query_neighbours(connection);
			typedef enum
			{ FIELD_sbjct_nid, FIELD_evalue, FIELD_length, FIELD_query_from, FIELD_query_to, FIELD_sbjct_from, FIELD_sbjct_to };

			{
				std::ostringstream query_str;

				query_str << "SELECT p.sbjct_nid, p.evalue, nrdb.length,"
						<< " p.query_from, p.query_to, p.sbjct_from, p.sbjct_to "
						<< " FROM pairsdb_90x90 p, nrdb,  "
						<< param_table_name.c_str() << " AS n"
						<< " WHERE p.sbjct_nid = n.nid "
						<< " AND nrdb.length <= " << query_length
						<< " AND p.sbjct_nid = nrdb.nid "
						<< " AND p.query_nid = " << query_nid
						<< " AND p.sbjct_nid != p.query_nid "
						<< " AND evalue <= " << param_max_evalue
						<< " ORDER BY nrdb.length DESC, p.sbjct_nid ASC"
						<< ends;

				query_neighbours.startQuery(query_str.str() );
			}

			Row row_neighbour;

			Nid last_nid = 0;

			std::vector< bool> aligned_residues;

			while (query_neighbours.fetchRow(row_neighbour))
			{

				Nid sbjct_nid = atol(row_neighbour[FIELD_sbjct_nid]);
				Position sbjct_length = atoi(row_neighbour[FIELD_length]);
				Position sbjct_from = atoi(row_neighbour[FIELD_sbjct_from]);
				Position sbjct_to = atoi(row_neighbour[FIELD_sbjct_to]);
				EValue evalue = atol(row_neighbour[FIELD_evalue]);

				// skip, if current representative has already been groupie
				if (groupies_nids.find(sbjct_nid) != groupies_nids.end())
					continue;

				if (sbjct_length == 0) // consistency check
					continue;

				if (sbjct_length > query_length)
					cout << "--> error in length " << sbjct_length << " > "
							<< query_length << endl;

				if (last_nid != sbjct_nid)
				{
					if (last_nid)
					{
						// calculate aligned positions
						Position total_aligned = 0;
						for (Position x = 1; x <= sbjct_length; x++)
							if (aligned_residues[x])
								total_aligned ++;
						if (param_loglevel >= 3)
							cout << "--> lali= " << total_aligned << " ref_l= "
									<< sbjct_length << " evalue=" << evalue
									<< endl;

						if (evalue <= param_max_evalue && (total_aligned
								>= (sbjct_length - param_min_domain_length))
								&& total_aligned
										>= ((double)sbjct_length * 0.9))
						{

							// mark sbjct as a groupie		
							groupies_nids.insert(sbjct_nid);
							out_rsdb << query_nid << SEPARATOR << sbjct_nid
									<< SEPARATOR << evalue << endl;
						}
					}

					if (param_loglevel >= 3)
						cout << "--> analysing query_nid=" << query_nid
								<< " sbjct_nid=" << sbjct_nid
								<< " query_length=" << query_length
								<< " sbjct_length=" << sbjct_length << endl;

					aligned_residues.clear();
					aligned_residues.resize(sbjct_length + 1, false);
				}

				if (param_loglevel >= 3)
					cout << "--> masking region from " << sbjct_from << " to "
							<< sbjct_to << endl;

				for (Position x = sbjct_from; x <= sbjct_to; x++)
					aligned_residues[x] = true;

				last_nid = sbjct_nid;

			} // end of loop over all neighbours

			query_neighbours.endQuery();
		}

	} // end of loop over all nids
	query_nids.endQuery();

	connection.Disconnect();

	out_rsdb.close();
	out_nids.close();
}

