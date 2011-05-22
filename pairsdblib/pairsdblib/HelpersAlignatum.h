//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: HelpersAlignatum.h,v 1.3 2004/06/02 09:03:58 aheger Exp $
//--------------------------------------------------------------------------------    


#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef PICASSO_HELPERS_ALIGNATUM_H
#define PICASSO_HELPERS_ALIGNATUM_H 1

#include <iosfwd>
#include <alignlib/alignlib.h>
#include "pairsdblib.h"

namespace alignlib {
  class Alignatum;
};

namespace pairsdblib {

    /** factory functions for creating alignatum-objects
	
	@author Andreas Heger
	@version $Id: HelpersAlignatum.h,v 1.3 2004/06/02 09:03:58 aheger Exp $
	@short protocol class for aligned objects

     */

    class Connection;

    /*--------------------------------------------------------------------*/
    /* Alignatum functions */
    alignlib::HAlignatum makeAlignatumNid( Sequence sequence,  
					    Nid nid = 0 );
    
    alignlib::HAlignatum makeAlignatumGroupie( Sequence sequence,  
						Identifier identifier, 
						Description description,
						Nid nid = 0, 
						Pid pide = 1.0, 
						Score score = 0.0 );
    
    alignlib::HAlignatum makeAlignatumDomain( Sequence sequence,  
					       Identifier identifier, 
					       Description description,
					       Nid nid = 0, 
					       Score score = 0.0,
					       ZScore zscore = 0.0);
    
    alignlib::HAlignatum makeAlignatumNeighbour( Sequence sequence,  
						  Identifier identifier, 
						  Description description,
						  Nid nid = 0, 
						  Pid pide = 1.0, 
						  Score score = 0.0,
						  EValue evalue = 0.0,
						  Filter filter = 90);
    
    alignlib::HAlignatum makeAlignatumRadar( Sequence sequence,
					      Score score = 0.0,
					      ZScore zscore = 0.0);
    
}; /* namespace */
#endif 






