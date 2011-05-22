//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: test_FillAlignata.cpp,v 1.2 2004/01/20 12:03:33 aheger Exp $
//--------------------------------------------------------------------------------

/** Test the MultipleAlignment - object
*/

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>

#include <alignlib.h>

#include "test_tools.h"
#include "Connection.h"
#include "HelpersAlignment.h"

using namespace std;
using namespace pairsdblib;

BOOST_AUTO_TEST_CASE( test_FillAlignataGroupie )
{
	Connection connection = connect();
	
	alignlib::HAlignment a = alignlib::makeAlignmentVector();

    fillAlignmentGroupie( a, connection, 3, 4, "pairsdb_100x90" );
    
    BOOST_CHECK_EQUAL( a->getRowFrom(), 0);
    BOOST_CHECK_EQUAL( a->getColFrom(), 0);

}

BOOST_AUTO_TEST_CASE( test_FillAlignata100x )
{
	Connection connection = connect();
	
	alignlib::HAlignment a = alignlib::makeAlignmentVector();

    fillAlignment100x( a, connection, 20588, 40 );
    
    BOOST_CHECK_EQUAL( a->getRowFrom(), 0);
    BOOST_CHECK_EQUAL( a->getColFrom(), 0);
}



