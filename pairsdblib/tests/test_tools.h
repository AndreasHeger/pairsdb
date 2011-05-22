#ifndef TEST_TOOLS_H_
#define TEST_TOOLS_H_

#include "tests.h"
#include "pairsdblib.h"
#include <cstdlib>
#include <iostream>

using namespace pairsdblib;

#define BOOST_TEST_MODULE
#include <boost/test/included/unit_test.hpp>
using boost::unit_test::test_suite;

Connection connect()
{
  
	Connection connection( HOST, USER, PASSWORD, PORT, SOCKET );
 
	connection.Connect( DATABASE );

	if (!connection.isConnected()) 
	{
		std::cerr << "could not connect to test database: host=" << HOST
			<< " database=" << DATABASE
			<< " user=" << USER
			<< " password=" << PASSWORD
			<< " port=" << PORT
			<< " socket=" << SOCKET
			<< std::endl;
		exit(EXIT_FAILURE);
	}

	return connection;
}
  
#endif /*TEST_TOOLS_H_*/
