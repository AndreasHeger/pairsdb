//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: test_Query.cpp,v 1.1.1.1 2002/07/02 14:16:54 heger Exp $
//--------------------------------------------------------------------------------

/** Test the MultipleAlignment - object
*/

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <strstream>

#ifdef WITH_DMALLOC
#include <dmalloc.h>
#endif

#include "Connection.h"
#include "Query.h"
#include "Row.h"

#include "tests.h"
#include "test_tools.h"

using namespace std;
using namespace pairsdblib;

BOOST_AUTO_TEST_CASE( test_Query )
{
	Connection connection = connect();

	Query query( connection );
    
	{ 
		ostrstream query_str;
	
		query_str << "SELECT nid, hid, sequence, length FROM nrdb WHERE nid = 1" << ends;
		query.startQuery( query_str.str());
    
		query_str.freeze( false );
	}

	Row row;

	query.fetchRow( row );
    
  	BOOST_CHECK_EQUAL( atoi(row[0]), 1);
  	BOOST_CHECK_EQUAL( row[1], "OVo48fOuViJtTRAyfLEjKA" );
  	BOOST_CHECK_EQUAL( strlen(row[2]), atoi(row[3]) );
 
  	query.endQuery();

  	connection.Disconnect();

}



