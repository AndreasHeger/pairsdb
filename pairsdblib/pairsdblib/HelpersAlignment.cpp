//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: HelpersAlignment.cpp,v 1.4 2004/06/02 09:03:58 aheger Exp $
//--------------------------------------------------------------------------------    


#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include <alignlib/alignlib.h>
#include <alignlib/Alignment.h>
#include <alignlib/HelpersAlignment.h>

#include "HelpersAlignment.h"
#include "SQLQueries.h"
#include "pairsdblib.h"
#include "Connection.h"
#include "Query.h"
#include "Row.h"

using namespace std;
using namespace alignlib;

namespace pairsdblib 
{
  
//----------------------------------------------------------------------------------------------------
void fillAlignmentGroupie( HAlignment & dest, 
			   Connection & connection, 
			   Nid rep_nid, 
			   Nid mem_nid, 
			   const char * table_name )
{
    
  Query query_groupies( connection );
  {
    ostringstream query_str;
    
    query_str << "SELECT rep_from, rep_ali, mem_from, mem_ali "
	      << " FROM " << table_name << " WHERE rep_nid = " << rep_nid
	      << " AND mem_nid = " << mem_nid << ends;
    
    query_groupies.startQuery( query_str.str());
  }
  
  typedef enum { FIELD_rep_from, FIELD_rep_ali, 
		 FIELD_mem_from, FIELD_mem_ali } MY_FIELDS;
    
    
    Row row;
    query_groupies.fetchRow( row );

    std::string rep_ali (row[FIELD_rep_ali]);
    std::string mem_ali (row[FIELD_mem_ali]);

    AlignmentFormatEmissions( 
    		atoi(row[FIELD_mem_from]) - OFFSET_CORRECTION,
    		mem_ali,
    		atoi(row[FIELD_rep_from]) - OFFSET_CORRECTION,
    		rep_ali).copy(dest);
    

    query_groupies.endQuery();
}



//--------------------------------------------------------------------------------------------------
void fillAlignment100x( HAlignment & dest, 
			Connection & connection, 
			Nid nrdb_nid,
			Filter level)
{
  
  assert( level < 90);

  dest->clear();

  /** fill Alignment with alignment 100x90 */
  Nid nrdb90_nid = SQL_GetRepresentative(connection, nrdb_nid, 100, 90 );
  if (nrdb90_nid == 0)
    nrdb90_nid = nrdb_nid;

  Nid nrdb40_nid = SQL_GetRepresentative(connection, nrdb90_nid, 90, level);
  if (nrdb40_nid == 0)
    nrdb40_nid = nrdb90_nid;
  
  // The procedure below could be optimized by exiting, but the version below is easier to understand
  HAlignment map_90x40(makeAlignmentVector());
  HAlignment map_100x90(makeAlignmentVector());

  // cout << " nrdb40_nid=" << nrdb40_nid << " nrdb90_nid=" << nrdb90_nid << " nrdb_nid=" << nrdb_nid << endl;
  /** now build the alignments and combine them: */
  if (nrdb90_nid == nrdb40_nid)
    {
      Position length = SQL_GetSequenceLength( connection, nrdb90_nid );
      map_90x40->addDiagonal(0, length );
    }
  else
    {
      char tbl[100];
      sprintf( tbl, "pairsdb_90x%i", level);
      fillAlignmentGroupie( map_90x40, connection, nrdb40_nid, nrdb90_nid, tbl);
  }
  
  if (nrdb_nid == nrdb90_nid)
    {
      Position length = SQL_GetSequenceLength( connection, nrdb_nid );
      map_100x90->addDiagonal(0,length);
    }
  else
    {
      fillAlignmentGroupie( map_100x90, connection, nrdb90_nid, nrdb_nid, "pairsdb_100x90");
  }

  combineAlignment( dest, map_100x90, map_90x40, CR);

}



} // namespace pairsdblib
















