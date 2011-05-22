//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: HelpersAlignandum.cpp,v 1.3 2004/06/02 09:03:58 aheger Exp $
//--------------------------------------------------------------------------------    


#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <alignlib.h>

#include "HelpersAlignandum.h"
#include "pairsdblib.h"
#include "Connection.h"
#include "Query.h"
#include "Row.h"
#include "PairsdblibDebug.h"

using namespace std;
using namespace alignlib;


namespace pairsdblib
{
  
/** add a limit clause to a query */
static void AddLimit( ostringstream & query_str, unsigned int max_lines = DEFAULT_MAX_LINES, unsigned int offset_lines = 0) {

  if (offset_lines > 0) {
    query_str << " LIMIT " << offset_lines << ", " << max_lines;
  } else {
    if (max_lines != DEFAULT_MAX_LINES)
      query_str << " LIMIT " << max_lines;
  }

}

HAlignandum makeSequenceFromPairsdb( 
		Connection & connection, 
		Nid nid )
{
    
	debug_func_cerr(5);
	
    Query query( connection );
    
    { 
        ostringstream query_str;
	
        query_str << "SELECT sequence FROM nrdb WHERE nid = " << nid << ends;
        query.startQuery( query_str.str());
    }
  
    Row row;

    query.fetchRow( row );

    HAlignandum result(alignlib::makeSequence( row[0] ));

    query.endQuery();
    
    return result;
}

//-------------------------------------------------------------------------------------------------
void fillProfileNeighbours( 
		HAlignandum & source, 
		Connection & connection, 
		Nid nid, 
		const char * table_name, 
		unsigned int max_lines,
		unsigned int offset_lines,
		Filter filter_upper_range,
		Filter filter_lower_range,
		EValue max_evalue)
{

	debug_func_cerr(5);

	// retrieve representative

	HProfile profile(toProfile(source));
  
	Query query_representative( connection );
  
	{ 
		ostringstream query_str;
    
		query_str << "SELECT sequence, length FROM nrdb WHERE nid = " << nid << ends;
		query_representative.startQuery( query_str.str());
	}
  
	Row row;

	query_representative.fetchRow( row );

	Length rep_length = atoi(row[1]);

	profile->resize( rep_length );

	HAlignment map_sequence2profile( makeAlignmentVector() );

	// add counts from representative
	map_sequence2profile->addDiagonal( 0, rep_length );

	HAlignandum sequence(makeSequence( row[0] ));
  
	profile->add( sequence, map_sequence2profile );
	map_sequence2profile->clear();

	query_representative.endQuery();

	// now collect neighbours and add them

	Query query_sbjcts( connection );
	{
		ostringstream query_str;
		
		query_str << " SELECT p.query_from, p.query_ali, "
	      	<< " p.sbjct_from, p.sbjct_ali, "
	      	<< " nrdb.sequence "
	      	<< " FROM " << table_name << " AS p, nrdb "
	      	<< " WHERE p.query_nid = " << nid << " AND nrdb.nid = p.sbjct_nid ";

		if (filter_lower_range != DEFAULT_LOWER_FILTER)
			query_str << " AND nrdb.filter >= " << filter_lower_range;
		if (filter_upper_range != DEFAULT_UPPER_FILTER)
			query_str << " AND nrdb.filter <= " << filter_upper_range;
		if (max_evalue != DEFAULT_MAX_EVALUE)
			query_str << " AND p.evalue <= " << max_evalue;

		AddLimit( query_str, max_lines, offset_lines );

		query_str << ends;
      
		query_sbjcts.startQuery( query_str.str());
	}
	
	typedef enum { FIELD_query_from, FIELD_query_ali, 
		 	FIELD_sbjct_from, FIELD_sbjct_ali, 
		 	FIELD_sequence } MY_FIELDS ;
  
	while (query_sbjcts.fetchRow( row ))
    {

		const std::string sbjct_ali = row[FIELD_sbjct_ali];
		const std::string query_ali = row[FIELD_query_ali];

		AlignmentFormatEmissions( 
				atoi(row[FIELD_sbjct_from]) - OFFSET_CORRECTION,
				sbjct_ali,
				atoi(row[FIELD_query_from]) - OFFSET_CORRECTION,
				query_ali).copy(map_sequence2profile);
    
		// sometimes there are errors probably due to blast parsing script. As
		// a result the alignment is longer than the profile. Skip these entries.
		if (map_sequence2profile->getColTo() <= rep_length)
		{
			HAlignandum sequence(makeSequence( row[FIELD_sequence] ));
			profile->add( sequence, map_sequence2profile );
		}

		map_sequence2profile->clear();
      
    }

	query_sbjcts.endQuery();
}

//-------------------------------------------------------------------------------------------------
void fillProfileDomains( HAlignandum & source, 
			 Connection & connection, 
			 Did did, 
			 const char * table_name_alignments,
			 const char * table_name_domains,
			 unsigned int max_lines,
			 unsigned int offset_lines )
{

	HProfile profile(toProfile(source));
    
  // retrieve length of domain
  Query query_length( connection );
  { 
    ostringstream query_str;
    
    query_str << "SELECT length FROM " << table_name_domains << " WHERE did = " << did << ends;
    query_length.startQuery( query_str.str());
  }
  
  Row row;
  query_length.fetchRow( row );
  
  Length  rep_length = atoi( row[0] );

  query_length.endQuery();

  //--------------------------------------------------------------------------------------------
  // build query
  Query query_domains( connection );
  {
    ostringstream query_str;
    
    query_str << "SELECT master_from, master_ali, domain_from, domain_ali, "
	      << " nrdb.sequence "
	      << " FROM " << table_name_alignments << " AS p, nrdb " 
	      << " WHERE master_did = " << did
	      << " AND nrdb.nid = p.domain_nid ";
      
    AddLimit( query_str, max_lines, offset_lines );
    
    query_str << ends;
    
    query_domains.startQuery( query_str.str());
  }

  typedef enum { FIELD_master_from, FIELD_master_ali, 
		 FIELD_domain_from, FIELD_domain_ali, 
		 FIELD_sequence } MY_FIELDS;

  // add counts
  profile->resize( rep_length );

  HAlignment map_sequence2profile(makeAlignmentVector());

  while (query_domains.fetchRow( row )) {
    
    const std::string master_ali = row[FIELD_master_ali];
    const std::string domain_ali = row[FIELD_domain_ali];

    AlignmentFormatEmissions( atoi(row[FIELD_domain_from]),
			      domain_ali,
			      atoi(row[FIELD_master_from]),
			      master_ali).copy(map_sequence2profile);
    
    HAlignandum sequence(makeSequence( row[FIELD_sequence] ));
    profile->add(sequence, map_sequence2profile );
    
    map_sequence2profile->clear();
    
  }
  
  query_domains.endQuery();
}



} // namespace pairsdblib
















