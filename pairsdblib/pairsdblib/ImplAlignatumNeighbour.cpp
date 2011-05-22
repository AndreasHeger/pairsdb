//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: ImplAlignatumNeighbour.cpp,v 1.2 2003/01/06 14:34:21 aheger Exp $
//--------------------------------------------------------------------------------    


#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>   /* for exp */
#include <alignlib.h>

#include "ImplAlignatumNeighbour.h"
#include "PairsdblibDebug.h"

using namespace std;
using namespace alignlib;

namespace pairsdblib
{

  HAlignatum  makeAlignatumNeighbour(const Sequence sequence,
				     const Identifier identifier, 
				     const Description description,
				     const Nid nid, 
				     const Pid pide, 
				     const Score score,
				     EValue evalue,
				     Filter filter) {
    std::string s(sequence);
    return HAlignatum (new ImplAlignatumNeighbour( s, identifier, description, nid, pide, score, evalue, filter )); 
  }

//----------------------------------< constructors and destructors >--------------------------------------
ImplAlignatumNeighbour::ImplAlignatumNeighbour (const std::string & sequence, 
						const Identifier identifier, 
						const Description description,
						const Nid nid, 
						const Pid pide, 
						const Score score,
						EValue evalue,
						Filter filter) :
      ImplAlignatum( sequence, 1, 0 ),  
      mNid (nid), 
      mIdentifier(identifier), 
      mDescription(description), 
      mScore(score), 
      mPercentIdentity ( pide ),
      mEvalue(evalue),
      mFilter(filter) {
}

ImplAlignatumNeighbour::~ImplAlignatumNeighbour () {
}

  ImplAlignatumNeighbour::ImplAlignatumNeighbour (const ImplAlignatumNeighbour & src ):
    ImplAlignatum( src ),
    mNid ( src.mNid ), 
    mIdentifier (src.mIdentifier), 
    mDescription( src.mDescription),
    mScore( src.mScore ),
    mPercentIdentity( src.mPercentIdentity),
    mEvalue(src.mEvalue), 
    mFilter(src.mFilter) {
}

//-------------------------------------------------------------------------------------------------
HAlignatum ImplAlignatumNeighbour::getClone() const
{
	debug_func_cerr( 5 );
  return HAlignatum(new ImplAlignatumNeighbour(*this)); 
}

//-------------------------------------------------------------------------------------------------
  void ImplAlignatumNeighbour::write( std::ostream & output ) const
  {
    ImplAlignatum::write( output );		

    output << getFieldSeparator() << mNid 
	   << getFieldSeparator() << mFilter
	 << getFieldSeparator() << mScore
	 << getFieldSeparator() << mPercentIdentity * 100
	 << getFieldSeparator() << mEvalue
	 << getFieldSeparator() << mIdentifier
	 << getFieldSeparator() << mDescription;
    
}

//--------------------------------------------------------------------------------------------

} // namespace pairsdblib
