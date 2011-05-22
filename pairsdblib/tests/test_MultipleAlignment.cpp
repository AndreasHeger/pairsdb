//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: test_MultipleAlignment.cpp,v 1.3 2004/01/20 12:03:33 aheger Exp $
//--------------------------------------------------------------------------------

/** Test the MultipleAlignment - object
*/

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>

#include <alignlib.h>

#include "Connection.h"
#include "HelpersMultipleAlignment.h"

#include "test_tools.h"

using namespace std;
using namespace pairsdblib;

BOOST_AUTO_TEST_CASE( test_MultipleAlignmentNeighbours )
{
	Connection connection = connect();
	
	alignlib::HMultipleAlignment m = makeMultipleAlignmentNeighbours( connection, 1, "pairsdb_90x90");
	BOOST_CHECK_EQUAL( m->getNumSequences(), 6 );
}

BOOST_AUTO_TEST_CASE( test_MultipleAlignmentGroupies )
{
	Connection connection = connect();
	
	alignlib::HMultipleAlignment m = makeMultipleAlignmentGroupies( connection, 28185, "pairsdb_100x90");
	BOOST_CHECK_EQUAL( m->getNumSequences(), 5 );

}

BOOST_AUTO_TEST_CASE( test_MultipleAlignmentRadar )
{
	Connection connection = connect();

	{
		alignlib::HMultipleAlignment m = makeMultipleAlignmentRadar( connection, 64, 1, "nrdb40_radar_domains" );
		BOOST_CHECK_EQUAL( m->getNumSequences(), 2 );
	}	
	
	{
		alignlib::HMultipleAlignment m = makeMultipleAlignmentRadar( connection, 64, 2, "nrdb40_radar_domains" );
		BOOST_CHECK_EQUAL( m->getNumSequences(), 2 );
	}	

}


/* makeMultipleAlignmentDomains not tested, used to be
    m = makeMultipleAlignmentDomains( connection, 1, "pairsdb.modules_alignments", "pairsdb.modules" );
*/















