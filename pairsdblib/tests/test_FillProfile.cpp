//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: test_FillProfile.cpp,v 1.2 2004/01/20 12:03:33 aheger Exp $
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
#include "HelpersAlignandum.h"

#include "test_tools.h"
        
using namespace std;
using namespace pairsdblib;

BOOST_AUTO_TEST_CASE( test_FillProfileNeighbours )
{
	Connection connection = connect();

	alignlib::HAlignandum p(alignlib::makeProfile());

	fillProfileNeighbours( p, connection, 3, "pairsdb_90x90" );
    p->prepare();
}


