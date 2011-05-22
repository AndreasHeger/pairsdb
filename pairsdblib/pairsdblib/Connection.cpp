//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: Connection.cpp,v 1.2 2004/06/02 09:03:58 aheger Exp $
//--------------------------------------------------------------------------------    


#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "PairsdblibException.h"
#include "Connection.h"
#include "pairsdblib.h"

#include <mysql.h>
#include <stdio.h>
#include <time.h>

#include "Query.h"
#include "Row.h"

using namespace std;

namespace pairsdblib
{

//---------------------------------------------------------< constructors and destructors >----------------------
Connection::Connection () : 
  mConnection(NULL), 
  mHost(""), mUser(""), mPassword(""), 
  mSocket(""), mPort(0), mClientFlag(0) 
{
}

Connection::Connection( const std::string & host, const std::string & user, 
			const std::string & password, const std::string & socket ) : 
  mConnection(NULL), 
  mHost(host), mUser(user), mPassword(password), 
  mSocket(socket), mPort(0), mClientFlag(0) 
{
  debug_func_cerr(5);
}

Connection::Connection( 
		const std::string & host, 
		const std::string & user, 
		const std::string & password, 
		int port) : 
  mConnection(NULL), 
  mHost(host), mUser(user), mPassword(password), 
  mSocket(""), mPort(port), mClientFlag(0) 
{
  debug_func_cerr(5);
}

Connection::Connection( const std::string & host, const std::string & user, 
			const std::string & password, int port,
			const std::string & socket ) : 
  mConnection(NULL), 
  mHost(host), mUser(user), mPassword(password), 
  mSocket(socket), mPort(port), mClientFlag(0) 
{
  debug_func_cerr(5);
}

Connection::~Connection () 
{
  debug_func_cerr( 5 );
  
  Disconnect();
}

Connection::Connection (const Connection & src ) : 
  mConnection( src.mConnection) 
{
  debug_func_cerr( 5 );
    
  setHost(src.mHost);
  setUser(src.mUser);
  setPassword(src.mPassword);
  setSocket(src.mSocket);
  setPort(src.mPort);
  setClientFlag(src.mClientFlag);
}

//---------------------------------------------------------------------------------------------------------------
bool Connection::Connect(const std::string & database) 
{
  debug_func_cerr( 5 );
  
  mConnection = mysql_init(NULL);
  
  if (mConnection == NULL) 
    throw PairsdblibException( "mysql_init() failed in Connection::Connect");

  debug_cerr( 5, "connecting with host=" << mHost
	      << " user=" << mUser
	      << " pass=" << mPassword
	      << " database=" << database
	      << " port=" << mPort
	      << " socket=" << mSocket
	      << " flag=" << mClientFlag );
  
  if (!mysql_real_connect( mConnection, 
		  ( mHost != "" ? mHost.c_str() : 0 ), 
		  ( mUser != "" ? mUser.c_str() : 0 ), 
		  ( mPassword != "" ? mPassword.c_str() : 0) , 
		  ( database != "" ? database.c_str() : 0), 
		  mPort,
		  ( mSocket != "" ? mSocket.c_str() : 0 ),					  					
		  mClientFlag ) )
    {
      std::cerr << "Failed to connect to database: Error: ---\n" << mysql_error(mConnection) << "\n---" << std::endl;

      throw PairsdblibException( "Connection failed in Connection::Connect");
    }
  
  return true;
}

//------------------------------------------------------------------------------------------------------
bool Connection::Disconnect()
{
  debug_func_cerr( 5 );
  if (isConnected()) 
    mysql_close( mConnection );
  mConnection = NULL;
  return true;
}

//------------------------------------------------------------------------------------------------------
bool Connection::isConnected() const
{
    if (mConnection) 
	return true;
    else 
	return false;
}

//------------------------------------------------------------------------------------------------------
void Connection::setDatabase( const std::string & new_database)
{
  debug_func_cerr( 5 );

  std::ostringstream query_str; 
    
  query_str << "USE " << new_database << ends;
  
  Query query( *this );
  query.startQuery( query_str.str() );
  query.endQuery();
}

//------------------------------------------------------------------------------------------------------
std::string Connection::getDatabase() 
{

  std::ostringstream query_str; 
    
  query_str << "SELECT DATABASE() " << endl;
  
  Query query( *this );
  query.startQuery( query_str.str() );

  Row  row;

  query.fetchRow( row );

  std::string result = row[0];

  query.endQuery();

  return result;
}


  /** member access functions */

  const std::string & Connection::getHost() const { return mHost; }
  void Connection::setHost( const std::string & host) 
  {
	  mHost = host;
  }

  const std::string & Connection::getUser() const { return mUser; }
  void Connection::setUser( const std::string & user) 
  {
	  mUser = user;
  }

  void Connection::setPassword( const std::string & password) 
  {
	  mPassword = password;
  }

  const std::string & Connection::getSocket() const { return mSocket; }
  void Connection::setSocket( const std::string & socket) 
  {
	  mSocket = socket;
  }

  int Connection::getPort() const { return mPort; }
  void Connection::setPort( int port ) {
    mPort = port;
  }

  int Connection::getClientFlag() const { return mClientFlag; }
  void Connection::setClientFlag( int clientflag ) {
    mClientFlag = clientflag;
  }


}



