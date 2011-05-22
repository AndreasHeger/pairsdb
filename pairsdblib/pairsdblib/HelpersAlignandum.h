//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: HelpersAlignandum.h,v 1.1.1.1 2002/07/02 14:16:54 heger Exp $
//--------------------------------------------------------------------------------    


#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef PICASSO_ALIGNATUM_H
#define PICASSO_ALIGNATUM_H 1

#include <iosfwd>
#include <alignlib/alignlib.h>
#include "pairsdblib.h"

namespace alignlib {
    class Alignandum;
};

namespace pairsdblib {

    /** factory functions for creating alignatum-objects
	
	@author Andreas Heger
	@version $Id: HelpersAlignandum.h,v 1.1.1.1 2002/07/02 14:16:54 heger Exp $
	@short protocol class for aligned objects

     */

    class Connection;


    /** retrieve a sequence from Pairsdb and put it into Alignandum-object 
     */
    alignlib::HAlignandum makeSequenceFromPairsdb( Connection & connection,
						    Nid nid );

    /* methods for filling profiles */

    void fillProfileNeighbours( alignlib::HAlignandum & dest, 
				Connection & connection, 
				Nid nid, 
				const char * table_name,
				unsigned int max_lines = DEFAULT_MAX_LINES,
				unsigned int offset_lines = 0,
				Filter filter_upper_range = DEFAULT_UPPER_FILTER,
				Filter filter_lower_range = DEFAULT_LOWER_FILTER,
				EValue max_evalue = DEFAULT_MAX_EVALUE); 
						  
    void fillProfileDomains( alignlib::HAlignandum & dest, 
			     Connection & connection, 
			     Did did,
			     const char * table_name_alis,
			     const char * table_name_info,
			     unsigned int max_lines = DEFAULT_MAX_LINES,
			     unsigned int offset_lines = 0 );


    
    
}; /* namespace */
#endif 






