//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: ImplAlignatumDomain.cpp,v 1.2 2003/01/06 14:34:21 aheger Exp $
//--------------------------------------------------------------------------------    


#include <iostream>
#include <iomanip>
#include <string>
#include <alignlib.h>

#include "ImplAlignatumDomain.h"
#include "PairsdblibDebug.h"
using namespace std;
using namespace alignlib;

namespace pairsdblib {

  HAlignatum makeAlignatumDomain(const Sequence sequence,
				  const Identifier identifier, 
				  const Description description,
				  const Nid nid, 
				  const Score score,
				  const ZScore zscore) {
    std::string s(sequence);
    return HAlignatum(new ImplAlignatumDomain( s, identifier, description, nid, score, zscore )); 
  }
    

//---------------------------------------------------------< constructors and destructors >-----------
ImplAlignatumDomain::ImplAlignatumDomain (const std::string & sequence, 
					  const Identifier identifier, 
					  const Description description,
					  const Nid nid, 
					  const Score score,
					  const ZScore zscore ) :
  ImplAlignatum( sequence, 1, 0 ),  
  mNid (nid), 
  mIdentifier(identifier), 
  mDescription(description), 
  mScore(score),
  mZScore( zscore ) {
}

ImplAlignatumDomain::~ImplAlignatumDomain () {
}

ImplAlignatumDomain::ImplAlignatumDomain (const ImplAlignatumDomain & src ) :
  ImplAlignatum( src ),
  mNid ( src.mNid ), 
  mIdentifier (src.mIdentifier), 
  mDescription( src.mDescription),
  mScore( src.mScore ), 
  mZScore( src.mZScore ) {
}

//--------------------------------------------------------------------------------------------------------------------------------
HAlignatum ImplAlignatumDomain::getClone() const
{
	debug_func_cerr( 5 );
  return HAlignatum(new ImplAlignatumDomain(*this)); 
}

//--------------------------------------------------------------------------------------------------------------------------------
  void ImplAlignatumDomain::write( std::ostream & output ) const
  {
	
  ImplAlignatum::write( output );
	
  output << getFieldSeparator() << mNid 
	 << getFieldSeparator() << mScore 
	 << getFieldSeparator() << mZScore
	 << getFieldSeparator() << mIdentifier
	 << getFieldSeparator() << mDescription;
}

//-------------------------------------------------------------------------------------------------

} // namespace pairsdblib
