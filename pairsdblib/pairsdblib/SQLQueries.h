//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: SQLQueries.h,v 1.2 2003/01/06 14:34:22 aheger Exp $
//--------------------------------------------------------------------------------    


#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef SQL_QUERIES_H
#define SQL_QUERIES_H 1

#include <iosfwd>
#include "pairsdblib.h"
#include <alignlib/alignlib.h>
#include <vector>
#include <utility>

namespace pairsdblib
{

    /** factory functions for creating alignatm-objects
	
	@author Andreas Heger
	@version $Id: SQLQueries.h,v 1.2 2003/01/06 14:34:22 aheger Exp $
	@short protocol class for aligned objects

     */

  typedef std::pair< alignlib::Position, alignlib::Position> SQLType_Boundaries;

  //--------------------------------------------------------------------------------------
  /** Data structure mapping query to sbjct
   */
  struct NeighborLink
  {
  NeighborLink() : 
    mQueryNid(0), mSbjctNid(0), mEValue(0),
      mMapQuery2Sbjct(),
      mPid(0), mScore(0) {};
  NeighborLink( const Nid & query_nid,
		const Nid & sbjct_nid,
		const EValue & evalue,
		const alignlib::Position & query_from,
		const std::string & query_ali,
		const alignlib::Position & sbjct_from,
		const std::string & sbjct_ali,
		const Score & score,
		const Pid & pid ) :
    mQueryNid( query_nid ),
      mSbjctNid( sbjct_nid ),
      mEValue ( evalue ),
      mMapQuery2Sbjct( query_from, query_ali, sbjct_from, sbjct_ali ),
      mPid ( pid ),
      mScore( score )
    {
    }
						
		
    Nid mQueryNid;
    Nid mSbjctNid;
    EValue mEValue;
    alignlib::AlignmentFormatEmissions mMapQuery2Sbjct;
    Pid mPid;
    Score mScore;
  };
  

  /** Data structure mapping query to sbjct
   */
  

  
  class Connection;
  
  void SQL_GetSequenceDomainBoundaries( Connection & connection, 
					std::vector< SQLType_Boundaries> & dest,
					const Nid nid,
					const char * table_name = "domain_alignments" );

  void SQL_GetSequenceMaskBoundaries( Connection & connection, 
				      std::vector< SQLType_Boundaries> & dest,
				      const Nid nid,
				      const Method method,
				      const char * table_name = "nrdb90_masks" );

  alignlib::Position SQL_GetSequenceLength( Connection & connection,
						       const Nid nid,
						       const char * table_name = "nrdb" );

   /** retrieve list of neighbours for nid */
  void SQL_GetNeighbours( 
		  	Connection & connection,
		  	std::vector< Nid > & result,
		  	Nid nid,
		  	const char * table_name = "pairsdb_90x90" );
  
  /** retrieve list of neighbours for nid */
  HNeighborLinkVector SQL_GetNeighbourLinks( 
					     Connection & connection,
					     const Nid & nid,
					     const std::string & table_name = "pairsdb_90x90" );
  

  /** @short return map of nid to length
      @param connection the connection
      @result the result
  */
  HMapNid2Length SQL_GetLengthMap( Connection & connection );

  /** @short return nids in a tabse
      @param connection the connection
      @param tablen_ame the table name
  */
  HNidVector SQL_GetNids( Connection & connection,
			  const std::string & table_name = "nrdb90" );

  /** return nid of representative for member given by nid. 
      return 0 if not found. */
  Nid SQL_GetRepresentative( Connection & connection, Nid nid, Filter from_level, Filter to_level );
  
}; /* namespace */
#endif 

