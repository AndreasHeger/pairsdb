//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: Row.cpp,v 1.2 2004/01/20 12:03:33 aheger Exp $
//--------------------------------------------------------------------------------    


#include <iostream>
#include <iomanip>
#include "pairsdblib.h"
#include "PairsdblibDebug.h"
#include "Row.h"

using namespace std;

namespace pairsdblib 
{

//---------------------------------------------------------< constructors and destructors >--------------------------------------
Row::Row () : mRow(NULL) {
}

Row::Row ( MYSQL_ROW src) : mRow(src) 
{
}
		       
Row::~Row () 
{
	debug_func_cerr( 5 );
}

Row::Row (const Row & src ) : mRow( src.mRow) 
{
}

//---------------------------------------------------------------------------------------------------------
void Row::Set( MYSQL_ROW src, int num_fields) {
  mRow = src;
  mNumFields = num_fields;
}

//---------------------------------------------------------------------------------------------------------
const char * Row::operator[](int index) const {
  return mRow[index];
}

//---------------------------------------------------------< Input/Output routines >-----------------------

std::ostream & operator<< (std::ostream & output, const Row & src) {
    for (int i = 0; i < src.mNumFields; i++) 
	output << "\t" << src[i];

    return output;
}

std::istream & operator>> (std::istream & input, Row & src) {
    return input;
}



}
