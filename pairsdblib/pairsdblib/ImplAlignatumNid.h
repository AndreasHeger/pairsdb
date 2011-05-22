//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: ImplAlignatumNid.h,v 1.1 2004/06/02 10:43:40 aheger Exp $
//--------------------------------------------------------------------------------    


#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef IMPL_ALIGNATUM_NID_H
#define IMPL_ALIGNATUM_NID_H 1

#include <iosfwd>
#include "pairsdblib.h"
#include <alignlib_fwd.h>
#include <alignlib/ImplAlignatum.h>

namespace pairsdblib {

    /** implementation of an aligned sequence, in this case part of a nid alignment
	
	@author Andreas Heger
	@version $Id: ImplAlignatumNid.h,v 1.1 2004/06/02 10:43:40 aheger Exp $
	@short protocol class for aligned objects

     */


class ImplAlignatumNid : public alignlib::ImplAlignatum {

 public:
    //--------- constructors and desctructors

    ImplAlignatumNid(const std::string & sequence,  
		     const Nid nid);

    /** desctructor */
    virtual ~ImplAlignatumNid ();

    /** copy constructor */
    ImplAlignatumNid (const ImplAlignatumNid & src );

    /* member access functions--------------------------------------------------------------- */
    virtual alignlib::HAlignatum getClone() const;

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
    Nid   mNid;					

};



}
#endif /* _ALIGNATUMSEQUENCENIDS_H */

