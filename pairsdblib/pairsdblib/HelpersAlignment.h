//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: HelpersAlignment.h,v 1.2 2003/01/06 14:34:21 aheger Exp $
//--------------------------------------------------------------------------------    


#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef PICASSO_HELPERS_ALIGNATA_H
#define PICASSO_HELPERS_ALIGNATA_H 1

#include <iosfwd>
#include <alignlib_fwd.h>
#include "pairsdblib.h"
#include "Connection.h"

namespace pairsdblib
{

    /** factory functions for creating alignatum-objects
	
	@author Andreas Heger
	@version $Id: HelpersAlignment.h,v 1.2 2003/01/06 14:34:21 aheger Exp $
	@short protocol class for aligned objects

     */

    /*--------------------------------------------------------------------*/
    /* Fill functions */
    /* in the fill functions, rep is in col */
  void fillAlignmentGroupie( alignlib::HAlignment & dest,
			Connection & connection, 
			Nid rep_nid, 
			Nid mem_nid, 
			const char * table_name );
  
    /** fill an Alignment between nid and groupie representative, where it does not matter, who the nrdb40_representative is
     */
  void fillAlignment100x( alignlib::HAlignment & dest,
		       Connection & connection,
		       Nid nid,
		       Filter level);

    
}; /* namespace */
#endif 






