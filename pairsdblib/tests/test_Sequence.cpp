//--------------------------------------------------------------------------------
// Project LibAlign
//
// Copyright (C) 2000 Andreas Heger All rights reserved
//
// Author: Andreas Heger <heger@ebi.ac.uk>
//
// $Id: test_Sequence.cpp,v 1.1.1.1 2002/07/02 14:16:54 heger Exp $
//--------------------------------------------------------------------------------

/** Test the MultipleAlignment - object
*/

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <fstream>

#include <alignlib.h>
#include "test_tools.h"

#include "Connection.h"
#include "HelpersAlignandum.h"

using namespace std;
using namespace pairsdblib;

BOOST_AUTO_TEST_CASE( test_makeSequenceFromPairsdb )
{
	Connection connection = connect();

	alignlib::HAlignandum s( makeSequenceFromPairsdb( connection, 1 ) );  
	
  	BOOST_CHECK_EQUAL( s->asString(), "MELRHTPARDLDKFIEDHLLPNTCFRTQVKEAIDIVCRFLKERCFQGTADPVRVSKVVKGGSSGKGTTLRGRSDADLVVFLTKLTSFEDQLRRRGEFIQEIRRQLEACQREQKFKVTFEVQSPRRENPRALSFVLSSPQLQQEVEFDVLPAFDALGQWTPGYKPNPEIYVQLIKECKSRGKEGEFSTCFTELQRDFLRNRPTKLKSLIRLVKHWYQTCKKTHGNKLPPQYALELLTVYAWEQGSRKTDFSTAQGFQTVLELVLKHQKLCIFWEAYYDFTNPVVGRCMLQQLKKPRPVILDPADPTGNVGGGDTHSWQRLAQEARVWLGYPCCKNLDGSLVGAWTMLQKI" );
}






