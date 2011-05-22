//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: Row.h,v 1.1.1.1 2002/07/02 14:16:54 heger Exp $
//--------------------------------------------------------------------------------    

/**
   Interface class for rows from query results.
   This class wraps around MYSQL_ROW, so that I
   do not have to use any MYSQL datatypes.
   The implementation is pretty basic, this pointer
   just maps the operator[] to const char *. The user
   still has to do any conversion he needs to do.

   Note that the memory location is managed by the
   @ref Connection object and the mysql-library. This
   class is just for data read access.

   This class does no consistency checking, since the information
   about a row-length is stored in the mysql-result and not in the
   row and out-of-bounds-checking would cause additional overhead.
   
   This class could use some work. Should write some fancy operators
   so that it really behaves like a pointer.

*/

#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef _ROW_H
#define _ROW_H 1

#include <iostream>
#include <mysql.h>
#include "pairsdblib_fwd.h"

namespace pairsdblib {
    
class Row {
    // friends
    friend std::ostream & operator<<( std::ostream &, const Row &);
    friend std::istream & operator>>( std::istream &, Row &);
    // class member functions
 public:
    // constructors and desctructors
    /** empty constructor */
    Row();

    /** Construct a pointer pointing to a MYSQL_ROW location */
    Row( const MYSQL_ROW src);
    
    /** copy constructor */
    Row(const Row &);

    /** destructor */
    virtual ~Row ();

    /** Set a new pointer */
    void Set(const MYSQL_ROW src, int num_fields );

    /** access a field, fields are const char *, you have to do 
	type conversion yourself */
    virtual const char * operator[]( int index ) const;

 protected:
 private:

 public:
 protected:
    MYSQL_ROW mRow;
    int mNumFields;
 private:

};

}

#endif /* MULTALIPICASSOH */

