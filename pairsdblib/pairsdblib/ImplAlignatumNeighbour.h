//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: ImplAlignatumNeighbour.h,v 1.2 2003/01/06 14:34:22 aheger Exp $
//--------------------------------------------------------------------------------    


#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef IMPL_ALIGNATUM_NEIGHBOUR_H
#define IMPL_ALIGNATUM_NEIGHBOUR_H 1

#include <iosfwd>
#include "pairsdblib.h"
#include <alignlib_fwd.h>
#include <alignlib/ImplAlignatum.h>

namespace pairsdblib {

    /** implementation of an aligned sequence, in this case part of a neighbour alignment
	
	@author Andreas Heger
	@version $Id: ImplAlignatumNeighbour.h,v 1.2 2003/01/06 14:34:22 aheger Exp $
	@short protocol class for aligned objects

     */


class ImplAlignatumNeighbour : public alignlib::ImplAlignatum {

 public:
    //--------- constructors and desctructors

  ImplAlignatumNeighbour(const std::string & sequence,  
			 const Identifier identifier, 
			 const Description description,
			 const Nid nid = 0, 
			 const Pid pide = 0.0, 
			 const Score score = 0.0,
			 EValue evalue = 0.0,
			 Filter filter = 90);

  /** desctructor */
  virtual ~ImplAlignatumNeighbour ();
  
  /** copy constructor */
  ImplAlignatumNeighbour (const ImplAlignatumNeighbour & src );

  /* member access functions--------------------------------------------------------------- */
  virtual alignlib::HAlignatum getClone()  const;

  
 protected:
    /** Write a formatted line for use as output in a multiple alignment into a stream. This function gets called by
	GetRow and GetFragment. These two functions provide WriteLine with a sequence fragment and the residue coordinates,
	WriteLine then produces formatted output for this.
	@param output	stream to write to
	@param from	first residue of sequence segment
	@param out		last residue of sequence segment
	@param sequence	sequence segment
    */
  virtual void write( std::ostream & output ) const;
  
    
 private:


    /** nid */
    long   mNid;					

    /** the identifier */
    std::string mIdentifier;

    /** the description */
    std::string mDescription;

    /** alignment score */
    Score mScore;
    
    /** pide */
    Pid  mPercentIdentity;			

    /** evalue */
    EValue mEvalue;
    
    /** Filter */
    Filter mFilter;
    

};



}
#endif /* _ALIGNATUMSEQUENCENEIGHBOURS_H */

