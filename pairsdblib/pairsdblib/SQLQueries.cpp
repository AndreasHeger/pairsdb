//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: SQLQueries.cpp,v 1.3 2004/06/02 09:03:58 aheger Exp $
//--------------------------------------------------------------------------------    


#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdio.h>

#include <alignlib/alignlib.h>

#include "SQLQueries.h"

#include "pairsdblib.h"
#include "Connection.h"
#include "Query.h"
#include "Row.h"


using namespace std;
using namespace alignlib;


namespace pairsdblib 
{
  
  /** build an alignment of nrdb90 groupies */
  /** forward declarations */

  //-------------------------------------------------------------------------------------------------
  // get sequence of representative and load it 
  // do it directly via AlignatumSequence, because otherwise you will translate the sequence back and forth

  void SQL_GetSequenceDomainBoundaries( Connection & connection,
					std::vector< SQLType_Boundaries > & result,
					Nid nid,
					const char * table_name) {
    
    Query query( connection ); 
    std::ostringstream query_str; 
    
    query_str << "SELECT domain_from, domain_to  "
	      << " FROM " << table_name << " WHERE domain_nid = " << nid
	      << " ORDER BY domain_from ASC " 
	      << ends;
    
    query.startQuery( query_str.str() );

    pairsdblib::Row row; 
    
    while (query.fetchRow( row))
      result.push_back( SQLType_Boundaries(atoi(row[0]), atoi(row[1])) );
    
    query.endQuery();
  }

  //----------------------------------------------------------------------------------------
  void SQL_GetSequenceMaskBoundaries( Connection & connection,
				      std::vector< SQLType_Boundaries > & result,
				      Nid nid,
				      Method method,
				      const char * table_name) {
    
    Query query( connection ); 
    std::ostringstream query_str; 
    
    query_str << "SELECT first_res, last_res "
	      << " FROM " << table_name << " WHERE nid = " << nid
              << " AND method = " << method
	      << " ORDER BY first_res ASC " 
	      << ends;
    
    query.startQuery( query_str.str() );

    pairsdblib::Row row; 
    
    while (query.fetchRow( row))
      result.push_back( SQLType_Boundaries(atoi(row[0]), atoi(row[1])) );
    
    query.endQuery();
  }

//----------------------------------------------------------------------------------------
Position SQL_GetSequenceLength( Connection & connection,
					   Nid nid,
					   const char * table_name) {
    
    Query query( connection ); 
    std::ostringstream query_str; 
    
    query_str << "SELECT length "
    	      << " FROM " << table_name
	      << " WHERE nid = " << nid 
	      << ends;
    
    query.startQuery( query_str.str() );

    pairsdblib::Row row; 

    Position length = 0;

    if (query.fetchRow( row))
      length = atol(row[0]);

    query.endQuery();

    return length;
  }


//----------------------------------------------------------------------------------------
void SQL_GetNeighbours( Connection & connection,
			std::vector< Nid > & result,
			Nid nid,
			const char * table_name) {
  
  Query query( connection ); 
  std::ostringstream query_str; 
  
  query_str << "SELECT sbjct_nid "
	    << " FROM " << table_name 
	    << " WHERE query_nid = " << nid 
	    << ends;
    
  query.startQuery( query_str.str() );

  pairsdblib::Row row; 
  
  while (query.fetchRow( row)) 
    result.push_back( atol( row[0]) );
  
  query.endQuery();

}

//--------------------------------------------------------------------------------------------------
Nid SQL_GetRepresentative( Connection & connection, Nid nid, 
				Filter from_level, Filter to_level ) {
  
  char buffer[100];
  sprintf( buffer, "pairsdb_%ix%i", from_level,to_level );

  Nid return_nid = 0;

  Query query_groupies( connection );
  
  ostringstream query_str;
  
  query_str << "SELECT rep_nid FROM " << buffer << " WHERE mem_nid = " << nid << ends;
  query_groupies.startQuery( query_str.str());
    
  if (query_groupies.getResultNumRows() > 0) {
    Row row;
    query_groupies.fetchRow( row );
    
    return_nid = atol(row[0]);
  } else {
    return_nid = 0;
  }
  
  query_groupies.endQuery();
  return return_nid;
}

  //----------------------------------------------------------------------------------------
  HMapNid2Length SQL_GetLengthMap( 
			Connection & connection )

  {
    debug_func_cerr(5);

    HMapNid2Length result(new MapNid2Length() );
    
    Query query( connection ); 
    std::ostringstream query_str; 
  
    query_str << "SELECT nid, length FROM nrdb" << ends;
  
    query.startQuery( query_str.str() );
  
    pairsdblib::Row row; 
  
    while (query.fetchRow( row)) 	
      {
	(*result)[ atol( row[0] ) ] = atoi(row[1]);
      }
  
    query.endQuery();
  return result;
}

    
//----------------------------------------------------------------------------------------
HNeighborLinkVector SQL_GetNeighbourLinks( 
			   Connection & connection,
			   const Nid & nid,
			   const std::string & table_name) 
{
  debug_func_cerr(5);
  
  HNeighborLinkVector result(new NeighborLinkVector() );
	
  Query query( connection ); 
  std::ostringstream query_str; 
  
  query_str << "SELECT sbjct_nid, evalue, query_from, query_ali, sbjct_from, sbjct_ali, score, percent_identity "
	    << " FROM " << table_name 
	    << " WHERE query_nid = " << nid 
	    << ends;
  
  query.startQuery( query_str.str() );
  
  pairsdblib::Row row; 
  
  while (query.fetchRow( row)) 	
    {
      result->push_back( NeighborLink( nid,
				      atol( row[0] ),
				      atol( row[1] ),
				      atoi( row[2] ) - OFFSET_CORRECTION,
				      row[3],
				      atoi( row[4] ) - OFFSET_CORRECTION,
				      row[5],
				      atof(row[6]),
				      atof(row[7]) ) );
    }
  
  query.endQuery();

  return result;
}

  HNidVector SQL_GetNids( Connection & connection,
			  const std::string & table_name )
  {
    debug_func_cerr(5);

    HNidVector result(new NidVector() );
    
    Query query( connection ); 
    std::ostringstream query_str; 
  
    query_str << "SELECT DISTINCT nid FROM " << table_name << std::ends;

    query.startQuery( query_str.str() );
  
    pairsdblib::Row row; 
  
    while (query.fetchRow( row)) 	
      {
	(*result).push_back(atol(row[0]));
      }
  
    query.endQuery();
    return result;
  }
  

				
} // namespace pairsdblib







