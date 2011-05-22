//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: Query.h,v 1.2 2004/06/02 09:03:58 aheger Exp $
//--------------------------------------------------------------------------------    

/**
   Interface class for connecting to pairsdb
   and mediating sql-calls, etc.

 */

#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef QUERY_H
#define QUERY_H 1

#include <iostream>
#include <alignlib/MultipleAlignment.h>
#include <mysql.h>
#include "Row.h"
#include "Connection.h"

namespace pairsdblib 
{
    
class Query {
    // friends
    friend std::istream & operator>>( std::istream &, Query &);
    // class member functions
 public:
    // constructors and desctructors
    Query  (Connection & connection);
    
    Query  (const Query & src);

    virtual ~Query ();

    /** Perform a query and store/use the result */
    bool startQuery( const std::string & query, bool client_side = true ); 

    /** Fetch next row, @ref Row is a smart pointer over MYSQL_ROW. Since I was not
	quite sure how to create null pointers, this function fills in an existing object. 
    */
    bool fetchRow(Row & row);
    
    /** Retrieve size of result set: number of rows */
    int getResultNumRows() const ;
    
    /** Retrieve size of result set: number of columns (fields) */
    int getResultNumColumns() const ;  

    /** tell, if any result was obtained */
    bool isEmpty() const;

    /** ends a query by releasing the memory */
    void endQuery();

 protected:
 private:
    // member data
 public:
 protected:
    Connection & mConnection;
    
    /** pointer to result */
    MYSQL_RES * mResult;

    int mNumFields;

 private:

};

}

#endif /* MULTALIPICASSOH */

