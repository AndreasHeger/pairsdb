//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: HelpersMultipleAlignment.h,v 1.3 2004/01/27 12:14:28 aheger Exp $
//--------------------------------------------------------------------------------    


#if HAVE_CONFIG_H
#include <config.h>
#endif

#ifndef PICASSO_HELPERS_MULTIPLEALIGNMENT_H
#define PICASSO_HELPERS_MULTIPLEALIGNMENT_H 1

#include <iosfwd>
#include <alignlib_fwd.h>
#include "pairsdblib.h"

namespace pairsdblib {

    /** factory functions for creating alignatum-objects
	
	@author Andreas Heger
	@version $Id: HelpersMultipleAlignment.h,v 1.3 2004/01/27 12:14:28 aheger Exp $
	@short protocol class for aligned objects

     */

    class Connection;

    /*--------------------------------------------------------------------*/
    /** Generic multiple alignment factory functions */
    alignlib::HMultipleAlignment makeMultipleAlignmentGroupies( Connection & connection, 
								Nid nid, 
								const char * table_name,
								unsigned int max_lines = DEFAULT_MAX_LINES,
								unsigned int offset_lines = 0);

    alignlib::HMultipleAlignment makeMultipleAlignmentNeighbours( 
								 Connection & connection, 
								 Nid nid, 
								 const char * table_name,
								 unsigned int max_lines = DEFAULT_MAX_LINES,
								 unsigned int offset_lines = 0,
								 Filter filter_upper_range = DEFAULT_UPPER_FILTER,
								 Filter filter_lower_range = DEFAULT_LOWER_FILTER,
								 EValue max_evalue = DEFAULT_MAX_EVALUE
								  ); 
    
    alignlib::HMultipleAlignment makeMultipleAlignmentDomains( Connection & connection, 
							       Did did,
							       const char * table_name,
							       unsigned int max_lines = DEFAULT_MAX_LINES,
							       unsigned int offset_lines = 0,
							       bool insert_gaps_representative = false,
							       int max_insertion_length = 10);

    alignlib::HMultipleAlignment makeMultipleAlignmentRadar( Connection & connection, 
							     Nid nid, 
							     RepeatId repeat_id,
							     const char * table_name,
							     unsigned int max_lines = DEFAULT_MAX_LINES,
							     unsigned int offset_lines = 0);


}; /* namespace */
#endif 






