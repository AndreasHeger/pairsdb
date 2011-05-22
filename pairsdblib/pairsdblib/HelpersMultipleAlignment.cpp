//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: HelpersMultipleAlignment.cpp,v 1.6 2004/06/02 09:03:58 aheger Exp $
//--------------------------------------------------------------------------------    


#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <alignlib.h>

#include "pairsdblib.h"
#include "PairsdblibDebug.h"
#include "Connection.h"
#include "Query.h"
#include "Row.h"
#include "HelpersMultipleAlignment.h"
#include "HelpersAlignatum.h"

using namespace std;
using namespace alignlib;


namespace pairsdblib {
  
  /** build an alignment of nrdb90 groupies */
  /** forward declarations */

  //-------------------------------------------------------------------------------------------------
  // get sequence of representative and load it 
  // do it directly via AlignatumSequence, because otherwise you will translate the sequence back and forth

/** add a limit clause to a query */
static void AddLimit( ostringstream & query_str, unsigned int max_lines = DEFAULT_MAX_LINES, unsigned int offset_lines = 0) {

  if (offset_lines > 0) {
    query_str << " LIMIT " << offset_lines << ", " << max_lines;
  } else {
    if (max_lines != DEFAULT_MAX_LINES)
      query_str << " LIMIT " << max_lines;
  }

}

HMultipleAlignment makeMultipleAlignmentDomains( Connection & connection, 
						  Did did,
						  const char * table_name,
						  unsigned int max_lines,
						  unsigned int offset_lines,
						  bool compress_alignment,
						  int max_insertion_length)
{

  debug_func_cerr(5);

  HMultipleAlignment result(makeMultipleAlignmentDots( compress_alignment, max_insertion_length ));

  Query query_representative( connection );
  
  typedef enum { FIELD_rep_nid, FIELD_sequence, 
		 FIELD_identifier, FIELD_description,
		 FIELD_rep_from, FIELD_rep_ali, 
		 FIELD_mali_from, FIELD_mali_ali,
		 FIELD_score, FIELD_zscore } MY_FIELDS;


  Query query( connection );
  
  {
    ostringstream query_str;
    
    query_str << "SELECT rep_nid, sequence, identifier, description, " 
	      << " rep_from, rep_ali, domain_from, domain_ali"
	      << " FROM " << table_name << ", nrdb AS n "
	      << " WHERE family = '" << did << "'"
	      << " AND n.nid = rep_nid " ;

    AddLimit( query_str, max_lines, offset_lines );
    
    query_str << ends;

    query.startQuery( query_str.str());
  }
  
  if ( query.isEmpty() ) {
    query.endQuery();
    return result;
  }

  Row row;

  while (query.fetchRow( row ))
    {

    const std::string rep_ali = row[FIELD_rep_ali];
    const std::string mali_ali = row[FIELD_mali_ali];
    
    HAlignment map_mali2rep(makeAlignmentVector());

    AlignmentFormatEmissions( atoi(row[FIELD_mali_from]),
			      mali_ali,
			      atoi(row[FIELD_rep_from]),
			      rep_ali ).copy( map_mali2rep );

    HAlignatum s ( makeAlignatumDomain( row[FIELD_sequence], 
					row[FIELD_identifier], 
					row[FIELD_description],
					atol(row[FIELD_rep_nid]))); 

    result->add( s, map_mali2rep, true, false, true, true, false);      
    
  }
  
  query.endQuery();
  return result;
}

//-------------------------------------------------------------------------------------------------
// get sequence of representative and load it 
// do it directly via AlignatumSequence, because otherwise you will translate the sequence back and forth

HMultipleAlignment makeMultipleAlignmentGroupies( 
		Connection & connection, 
		Nid nid, 
		const char * table_name,
		unsigned int max_lines,
		unsigned int offset_lines)
{
  debug_func_cerr(5);

  HMultipleAlignment result(makeMultipleAlignment());
    
  Query query_representative( connection );
  { 
    ostringstream query_str;
    
    query_str << "SELECT sequence, identifier, description FROM nrdb WHERE nid = " << nid << ends;
    query_representative.startQuery( query_str.str());
	
  }
  Row row;
  while (query_representative.fetchRow( row )) {				       // should be just one row
    HAlignatum s(makeAlignatumGroupie( row[0], row[1], row[2], nid ));
    result->add( s );
#ifdef DEBUG
    cout << "Representative is " << *s << endl;
#endif
    }
    
  query_representative.endQuery();

    //--------------------------------------------------------------------------------------------
    // get groupies and alignments of representatives and load them from the corresponding table
	
    Query query_groupies( connection );
    {
	ostringstream query_str;
	
	query_str << "SELECT rep_from, rep_to, rep_ali, mem_nid, mem_from, mem_to, mem_ali, score, percent_identity, "
		  << " nrdb.sequence, nrdb.identifier, nrdb.description " 
		  << " FROM " << table_name << " AS p, nrdb WHERE rep_nid = " << nid
		  << " AND nrdb.nid = mem_nid ORDER BY percent_identity DESC";

	AddLimit( query_str, max_lines, offset_lines );
	
	query_str << ends;

	query_groupies.startQuery( query_str.str());
    }
    typedef enum { FIELD_rep_from, FIELD_rep_to, FIELD_rep_ali, 
		   FIELD_mem_nid, FIELD_mem_from, FIELD_mem_to, 
		   FIELD_mem_ali, FIELD_score, FIELD_pide, FIELD_sequence, 
		   FIELD_identifier, FIELD_description } MY_FIELDS;

    while (query_groupies.fetchRow( row )) {
      
      /* skip unaligned entries */
      if (atoi(row[FIELD_rep_from]) == 0) 
	continue;
	
      const std::string rep_ali = row[FIELD_rep_ali];
      const std::string mem_ali = row[FIELD_mem_ali];

      HAlignment map_rep2mem(makeAlignmentVector());

      AlignmentFormatEmissions( 
    		  	atoi(row[FIELD_rep_from]) - OFFSET_CORRECTION,
				rep_ali,
				atoi(row[FIELD_mem_from]) - OFFSET_CORRECTION,
				mem_ali ).copy( map_rep2mem );
      
      HAlignatum s(makeAlignatumGroupie( row[FIELD_sequence], 
    		  row[FIELD_identifier], 
    		  row[FIELD_description],
    		  atol(row[FIELD_mem_nid]), 
    		  atof(row[FIELD_pide]), 
    		  atof(row[FIELD_score])));

      result->add( s, map_rep2mem, true, false, true, true, false );

    }
    
    query_groupies.endQuery();

    return result;
}

//-------------------------------------------------------------------------------------------------
HMultipleAlignment makeMultipleAlignmentNeighbours( Connection & connection, 
						    Nid nid, 
						    const char * table_name, 
						    unsigned int max_lines,
						    unsigned int offset_lines,
						    Filter filter_upper_range,
						    Filter filter_lower_range,
						    EValue max_evalue)
{
  debug_func_cerr(5);

  HMultipleAlignment result( makeMultipleAlignment() );

  Query query_representative( connection );
  
  { 
    ostringstream query_str;
    
    query_str << "SELECT sequence, identifier, description, filter FROM nrdb WHERE nid = " << nid << ends;
    query_representative.startQuery( query_str.str());

  }
  
  Row row;

  Filter rep_filter = 100;
  
  while (query_representative.fetchRow( row )) {				       // should be just one row
    rep_filter = atoi(row[3]);
    HAlignatum s(makeAlignatumNeighbour( row[0], row[1], row[2], 0, 1.0, 0.0, 0.0, rep_filter  ));
    result->add( s );
    debug_cerr( 5, "Query is " << *s );
  }
 
  query_representative.endQuery();

  //--------------------------------------------------------------------------------------------
  // get groupies and alignments of representatives and load them from the corresponding table
	
  Query query_sbjcts( connection );
  {
    ostringstream query_str;
    
    query_str << " SELECT p.query_from, p.query_to, p.query_ali, "
	      << " p.sbjct_nid, p.sbjct_from, p.sbjct_to, p.sbjct_ali, "
	      << " p.score, p.percent_identity, p.evalue, " 
	      << " nrdb.sequence, nrdb.filter, nrdb.identifier, nrdb.description "
	      << " FROM " << table_name << " p, nrdb "
	      << " WHERE p.query_nid = " << nid << " AND nrdb.nid = p.sbjct_nid ";

    if (filter_lower_range != DEFAULT_LOWER_FILTER)
      query_str << " AND nrdb.filter >= " << filter_lower_range;
    if (filter_upper_range != DEFAULT_UPPER_FILTER)
      query_str << " AND nrdb.filter <= " << filter_upper_range;
    if (max_evalue != DEFAULT_MAX_EVALUE)
      query_str << " AND p.evalue <= " << max_evalue;
    
    query_str << " ORDER BY p.evalue";
    
    AddLimit( query_str, max_lines, offset_lines );
    
    query_str << ends;
    
    query_sbjcts.startQuery( query_str.str());
  }
  
  typedef enum { FIELD_query_from, FIELD_query_to, FIELD_query_ali, 
		 FIELD_sbjct_nid, FIELD_sbjct_from, FIELD_sbjct_to, FIELD_sbjct_ali, 
		 FIELD_score, FIELD_pide, FIELD_evalue, FIELD_sequence, FIELD_filter, FIELD_identifier, 
		 FIELD_description } MY_FIELDS ;
  
  int counter = 0;
  
  while (query_sbjcts.fetchRow( row ))
    {
      counter++;
      const std::string row_ali = row[FIELD_query_ali];
      const std::string col_ali = row[FIELD_sbjct_ali];
    
      HAlignment map_query2sbjct(makeAlignmentVector());
      
      AlignmentFormatEmissions( 
    		  atoi(row[FIELD_query_from]) - OFFSET_CORRECTION,
    		  row_ali,
    		  atoi(row[FIELD_sbjct_from]) - OFFSET_CORRECTION,
    		  col_ali ).copy( map_query2sbjct );

      
      HAlignatum s (makeAlignatumNeighbour( row[FIELD_sequence], 
					    row[FIELD_identifier], 
					    row[FIELD_description],
					    atol(row[FIELD_sbjct_nid]), 
					    atof(row[FIELD_pide]), 
					    atof(row[FIELD_score]),
					    atof(row[FIELD_evalue]),
					    atoi(row[FIELD_filter])
					    ));

      debug_cerr( 5, "adding neighbour " << counter );
    
      result->add( s, map_query2sbjct, true, false, true, true, false );
    
      //        if (counter > 400)
      //  	break;
	
    }
  
  query_sbjcts.endQuery();

  return result;
}


//----------------------------------------------------------------------------------------------------
Position CalculateAlignmentLength( const std::string & alignment )
{
  debug_func_cerr(5);

  std::istringstream is( alignment.c_str() );   
  
  Position result = 0, d = 0;
  while (is >> d) {
    if (d < 0) result -= d;
    if (d > 0) result += d;
  }
  return result;
}

//-----------------------------------------------------------------------------------------------------
HMultipleAlignment makeMultipleAlignmentRadar( Connection & connection, 
					       Nid nid,
					       RepeatId repeat_id,
					       const char * table_name,
					       unsigned int max_lines,
					       unsigned int offset_lines)
{
  debug_func_cerr(5);

  HMultipleAlignment result(makeMultipleAlignment());

  Query query_representative( connection );
  
  typedef enum { FIELD_sequence, FIELD_identifier, FIELD_description, 
		 FIELD_bw_from, FIELD_bw_to, FIELD_nrepreats, FIELD_diagonal,
		 FIELD_length, FIELD_score, FIELD_level } MY_FIELDS1 ;
  { 
    ostringstream query_str;
    
    query_str << "SELECT n.sequence " 
	      << " FROM nrdb AS n WHERE n.nid = " << nid 
	      << ends;
    query_representative.startQuery( query_str.str());
    
  }
  
  Row row;
  query_representative.fetchRow( row );
  
  std::string sequence(row[FIELD_sequence]);
  
  query_representative.endQuery();

  Query query_repeatunits( connection );
  
  {
    ostringstream query_str;
    
    query_str << "SELECT rep_from, rep_to, rep_ali, score, zscore "
	      << " FROM " << table_name 
	      << " WHERE family = \"" << nid << '-' << repeat_id << "\""
	      << " ORDER BY rep_from ";

    AddLimit( query_str, max_lines, offset_lines );
    
    query_str << ends;

    query_repeatunits.startQuery( query_str.str());
  }
  
  if ( query_repeatunits.isEmpty() ) {
    query_repeatunits.endQuery();
    return result;
  }

  typedef enum { FIELD2_repeat_from, FIELD2_repeat_to, FIELD2_repeat_ali, 
		 FIELD2_score, FIELD2_zscore } MY_FIELDS2 ;
  

  std::string row_ali; 
  
  Position repeat_length = 0;
  while (query_repeatunits.fetchRow( row ))
    {
    
      const std::string col_ali = row[FIELD2_repeat_ali];

      // the first time, calculate the repeat length and add first repeat explicitely
      if (repeat_length == 0)
	{

	  repeat_length = CalculateAlignmentLength( col_ali );

	  char * buffer = new char[10];

	  // this is a bit ugly. I rely on the fact, that the adding of pairs to the alignment stops, when
	  // no more characters can be emitted from the shorter aligned string.
	  
	  sprintf( buffer, "+%i", repeat_length);
	  row_ali = buffer;
	  delete [] buffer;
	
	  //-----------------------------------------------------------
	  // now add first repeat unit
	  HAlignment map_profile2repeat(makeAlignmentVector());

	  AlignmentFormatEmissions( 
			  1 - OFFSET_CORRECTION,
			  row_ali,
			  atoi(row[FIELD2_repeat_from]) - OFFSET_CORRECTION,
			  col_ali ).copy(map_profile2repeat);

      // create mapping of master to new mali
	HAlignment map_profile2mali(makeAlignmentVector());
	HAlignment map_repeat2mali(makeAlignmentVector());
	
	expandAlignment( map_profile2mali, 
			 map_repeat2mali, 
			 map_profile2repeat, 
			 false,
			 true,
			 true,
			 false,
			 repeat_length );
      
	HAlignatum s( makeAlignatumRadar( sequence.c_str(), 
					  atof(row[FIELD2_score]), 
					  atof(row[FIELD2_zscore])));

	s->mapOnAlignment( map_repeat2mali, repeat_length ); 

	// add to multiple alignment
	result->add( s );

	continue;
	
    }

      //-----------------------------------------------------------
      // create alignment
      // row = query
      // col = sbjct
      HAlignment map_profile2repeat(makeAlignmentVector());
      AlignmentFormatEmissions( 
    		  1 - OFFSET_CORRECTION,
    		  row_ali,
    		  atoi(row[FIELD2_repeat_from]) - OFFSET_CORRECTION,
    		  col_ali ).copy(map_profile2repeat);


      HAlignatum s( makeAlignatumRadar( sequence.c_str(), 
					atof(row[FIELD2_score]), 
					atof(row[FIELD2_zscore])));

    
      result->add( s, map_profile2repeat, true, false, true, true, false);      
    
  }
  
  query_repeatunits.endQuery();
  return result;
}


  
} // namespace pairsdblib
















