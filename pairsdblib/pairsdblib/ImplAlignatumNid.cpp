//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: ImplAlignatumNid.cpp,v 1.1 2004/06/02 10:43:40 aheger Exp $
//--------------------------------------------------------------------------------    


#include <iostream>
#include <iomanip>
#include <string>
#include "ImplAlignatumNid.h"
#include <alignlib.h>
#include "PairsdblibDebug.h"
using namespace std;
using namespace alignlib;

namespace pairsdblib {

  HAlignatum makeAlignatumNid(const Sequence sequence, const Nid nid)
  {
    std::string s(sequence);
    return HAlignatum( new ImplAlignatumNid( s, nid )); 
  }
  
//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplAlignatumNid::ImplAlignatumNid (const std::string & sequence, 
				    const Nid nid ):
  ImplAlignatum( sequence, 1, 0 ),  
  mNid (nid) {
}

ImplAlignatumNid::~ImplAlignatumNid () {
}

ImplAlignatumNid::ImplAlignatumNid (const ImplAlignatumNid & src ) :
  ImplAlignatum( src ),
  mNid ( src.mNid ) {
}


//--------------------------------------------------------------------------------------------------------------------------------
HAlignatum ImplAlignatumNid::getClone() const
{
	debug_func_cerr(5);
	return HAlignatum(new ImplAlignatumNid(*this)); 
}

//--------------------------------------------------------------------------------------------------------------------------------
  void ImplAlignatumNid::write( std::ostream & output ) const
  {

    ImplAlignatum::write( output);		

  output << getFieldSeparator() << mNid;

}

//--------------------------------------------------------------------------------------------------------------------------------

} // namespace pairsdblib
