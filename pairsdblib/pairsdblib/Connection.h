//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: Connection.h,v 1.1.1.1 2002/07/02 14:16:54 heger Exp $
//--------------------------------------------------------------------------------    

/**
   Interface class for connecting to pairsdb
   and mediating sql-calls, etc.

   Input strings for setting database, hostname, etc.
   are copied, the results are pointers to the data in the object
   and not copies, so they are const.

 */

#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef CONNECTION_PICASSO_H
#define CONNECTION_PICASSO_H 1

#include <iosfwd>
#include <string>
#include <cstring>
#include <mysql.h>

namespace pairsdblib 
{
    
class Connection 
{

  // class member functions
 public:

  /** constructor, connect via socket */
  Connection( const std::string & host, const std::string & user, const std::string & password, const std::string & socket );

    /** constructor, connect via port and socket */
  Connection( const std::string & host, const std::string & user, const std::string & password, int port, const std::string & socket );

  /** constructor, connect via port */
  Connection( const std::string & host, const std::string & user, const std::string & password, int port );
    
  Connection  (const Connection &);

  virtual ~Connection ();

  /** is connected */
  bool isConnected() const;

  /** Connect to Pairsdb
   */
  bool Connect(const std::string & database);
    
  /** Disconnect from Pairsdb
   */
  bool Disconnect();

  /** set database
   */
  void setDatabase( const std::string & database );

  /** get database
   */
  std::string getDatabase();
  
  /** export MYSQL connection structure */
  MYSQL * getConnection() { return mConnection;}

  /** member access functions */
  const std::string & getHost() const;
  void setHost( const std::string & host);

  const std::string & getUser() const;
  void setUser( const std::string & user);

  void setPassword( const std::string & password);

  const std::string & getSocket() const;
  void setSocket( const std::string & socket);
  
  int getPort() const;
  void setPort( int port );
  int getClientFlag() const;
  void setClientFlag( int clientflag );

 protected:

 private:
  // member data

  /** empty constructor */
  Connection();

 public:

 protected:
  /** pointer to mysql connection data structure */
  MYSQL * mConnection;

 private:
  std::string mHost;
  std::string mUser;
  std::string mPassword;
  std::string mSocket;
  unsigned int mPort;
  unsigned int mClientFlag;

};

}

#endif 




