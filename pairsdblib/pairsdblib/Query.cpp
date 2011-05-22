//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: Query.cpp,v 1.2 2004/06/02 09:03:58 aheger Exp $
//--------------------------------------------------------------------------------


#include <iostream>
#include <iomanip>
#include <fstream>
#include "PairsdblibException.h"
#include "PairsdblibDebug.h"
#include "Connection.h"
#include "Query.h"
#include "pairsdblib.h"
#include <stdio.h>
#include <time.h>

using namespace std;

namespace pairsdblib
{

//---------------------------------------------------------< constructors and destructors >----------
Query::Query ( Connection & connection) : mConnection(connection), mResult(NULL)
{
}

Query::~Query ()
{
  endQuery();
}

Query::Query (const Query & src ) :
  mConnection( src.mConnection), mResult(NULL)
{
}


//------------------------------------------------------------------------------------------------------
int Query::getResultNumRows() const
{
  if (mResult)
    return mysql_num_rows( mResult);
  else
    return 0;
}

//------------------------------------------------------------------------------------------------------
int Query::getResultNumColumns() const
{
    if (mResult)
	return mNumFields;
    else
	return 0;
}

//------------------------------------------------------------------------------------------------------
bool Query::isEmpty() const
{
  return (getResultNumRows() == 0);
}


//------------------------------------------------------------------------------------------------------
bool Query::startQuery(const std::string & query, bool client_side )
{

  debug_func_cerr(5);

  debug_cerr(5, "query=" << query );

  endQuery();

  if (mysql_query( mConnection.getConnection(), query.c_str() ) != 0)
    {
      std::string message = std::string( "Query failed in Query.cpp: ") + query;
      message += "\nError message:";
      message +=  mysql_error(mConnection.getConnection());
      throw PairsdblibException( message.c_str() );
    }

  if (client_side)
    {
      mResult = mysql_store_result( mConnection.getConnection() );
    }
  else
    {
      mResult = mysql_use_result( mConnection.getConnection() );
    }

  if (mResult)
    mNumFields = mysql_num_fields(mResult);
  else
    mNumFields = 0;

  return true;
}

//------------------------------------------------------------------------------------------------------
void Query::endQuery()
{
  if (mResult)
    {
      mysql_free_result( mResult );
      mResult = NULL;
    }
}

//------------------------------------------------------------------------------------------------------
bool Query::fetchRow( Row& row)
{

  if (mResult == NULL)
    return false;

  MYSQL_ROW new_row = mysql_fetch_row( mResult );

  if (!new_row)
    return false;

  row.Set( new_row, mNumFields );
  return true;
}

//---------------------------------------------------------< Input/Output routines >---------------------------------------------

std::istream & operator>> (std::istream & input, Query & src)
{
  return input;
}

}
