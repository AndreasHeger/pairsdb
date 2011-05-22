//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: test_Connection.cpp,v 1.1.1.1 2002/07/02 14:16:54 heger Exp $
//--------------------------------------------------------------------------------

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>

#include <cstdlib>
#include "Connection.h"

#include "test_tools.h"

using namespace std;
using namespace pairsdblib;

BOOST_AUTO_TEST_CASE( test_ConnectionPort )
{
	if (PORT)
	{
		Connection connection( HOST, USER, PASSWORD, PORT );
		connection.Connect( DATABASE );
		BOOST_CHECK_EQUAL(connection.isConnected(), true);
	}
}

BOOST_AUTO_TEST_CASE( test_ConnectionSocket )
{
	if (SOCKET)
	{
		Connection connection( HOST, USER, PASSWORD, SOCKET );
		connection.Connect( DATABASE );
		BOOST_CHECK_EQUAL(connection.isConnected(), true);
	}
}


