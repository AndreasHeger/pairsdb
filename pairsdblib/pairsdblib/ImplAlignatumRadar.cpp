//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: ImplAlignatumRadar.cpp,v 1.2 2003/01/06 14:34:22 aheger Exp $
//--------------------------------------------------------------------------------    


#include <iostream>
#include <iomanip>
#include <string>
#include <alignlib.h>

#include "ImplAlignatumRadar.h"
#include "PairsdblibDebug.h"

using namespace std;
using namespace alignlib;

namespace pairsdblib {

    /** desctructor */

  HAlignatum makeAlignatumRadar(const Sequence sequence,
				   const Score score,
				   const ZScore zscore) {
	std::string s(sequence);
	return HAlignatum(new ImplAlignatumRadar( s, score, zscore ));  
    }
    

//---------------------------------------------------------< constructors and destructors >-----------------
ImplAlignatumRadar::ImplAlignatumRadar(const std::string & sequence,  
				       const Score score,
				       const ZScore zscore) :
    ImplAlignatum( sequence, 1, 0 ),  
    mScore( score ),
    mZScore( zscore ) {
}

ImplAlignatumRadar::~ImplAlignatumRadar () {
}

ImplAlignatumRadar::ImplAlignatumRadar (const ImplAlignatumRadar & src ) :
    ImplAlignatum( src),
    mScore( src.mScore ),
    mZScore( src.mZScore ) {
}

//----------------------------------------------------------------------------------------------------
HAlignatum ImplAlignatumRadar::getClone() const
{
  return HAlignatum ( new ImplAlignatumRadar(*this) ); 
}

//----------------------------------------------------------------------------------------------------
  void ImplAlignatumRadar::write( std::ostream & output ) const

{

  ImplAlignatum::write( output );
	
  output << getFieldSeparator() << mScore 
	 << getFieldSeparator() << mZScore;
}    

//---------------------------------------------------------------------------------------------------

} // namespace pairsdblib
