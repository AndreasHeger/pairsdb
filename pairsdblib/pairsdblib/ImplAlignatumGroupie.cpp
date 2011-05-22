//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: ImplAlignatumGroupie.cpp,v 1.3 2004/06/02 09:03:58 aheger Exp $
//--------------------------------------------------------------------------------    


#include <iostream>
#include <iomanip>
#include <string>
#include "ImplAlignatumGroupie.h"
#include <alignlib.h>

using namespace std;
using namespace alignlib;

namespace pairsdblib {

  HAlignatum makeAlignatumGroupie(const Sequence sequence,
				  const Identifier identifier, 
				  const Description description,
				  const Nid nid, 
				  const Pid pide, 
				  const Score score)
  {
    std::string s(sequence);
    return HAlignatum( new ImplAlignatumGroupie( s, identifier, description, nid, pide, score ));
  }
  
//---------------------------------------------------------< constructors and destructors >--------------------------------------
ImplAlignatumGroupie::ImplAlignatumGroupie (const std::string & sequence, 
					    const Identifier identifier, 
					    const Description description,
					    const Nid nid, 
					    const Pid pide, 
					    const Score score ) :
  ImplAlignatum( sequence, 1, 0 ),  
  mNid (nid), 
  mIdentifier(identifier), 
  mDescription(description), 
  mScore(score), 
  mPercentIdentity (pide ) {
}

ImplAlignatumGroupie::~ImplAlignatumGroupie () {
}

ImplAlignatumGroupie::ImplAlignatumGroupie (const ImplAlignatumGroupie & src ) :
  ImplAlignatum( src ),
  mNid ( src.mNid ), 
  mIdentifier (src.mIdentifier), 
  mDescription( src.mDescription),
  mScore( src.mScore ),
  mPercentIdentity( src.mPercentIdentity)
{
}


//--------------------------------------------------------------------------------------------------------------------------------
HAlignatum ImplAlignatumGroupie::getClone() const {
  debug_func_cerr(5);
  
  return HAlignatum ( new ImplAlignatumGroupie( *this) ); 
}

//-------------------------------------------------------------------------------------------------------------------------------
  void ImplAlignatumGroupie::write( std::ostream & output ) const
  {
    ImplAlignatum::write( output );		

    output << getFieldSeparator() << mNid 
	   << getFieldSeparator() << mScore 
	   << getFieldSeparator() << mPercentIdentity
	   << getFieldSeparator() << mIdentifier
	   << getFieldSeparator() << mDescription;

  }

//--------------------------------------------------------------------------------------------------------------------------------

} // namespace pairsdblib
