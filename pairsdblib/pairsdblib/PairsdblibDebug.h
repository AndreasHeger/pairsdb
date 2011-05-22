/*
  pairsdblib - a library for aligning protein sequences

  $Id: pairsdblib.h,v 1.5 2005/02/24 11:07:25 aheger Exp $

  Copyright (C) 2004 Andreas Heger
  
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#ifndef _PAIRSDBLIB_DEBUG_H
#define _PAIRSDBLIB_DEBUG_H
#include <string>
#include <ostream>

namespace pairsdblib
{
  
  /* this debugging method was inspired by Leo Goodstadt's code
   */
  
#ifdef DEBUG
  // The inclusion of nl_types.h was necessary for this to work on cygwin
#include <nl_types.h>
#if defined __cplusplus
#if __GNUC_PREREQ(2, 6)
#define __ASSERT_FUNCTION __PRETTY_FUNCTION__
#endif
#elif __GNUC_PREREQ(2, 4)
#define __ASSERT_FUNCTION __PRETTY_FUNCTION__
#else
#if defined __STDC_VERSION__ && __STDC_VERSION__ >= 199901L
#define __ASSERT_FUNCTION __func__
#else
#define __ASSERT_FUNCTION ((__const char *) 0)
#endif
#endif

  struct PairsdblibDebug
  {

    PairsdblibDebug(std::ostream& os, unsigned int i, const std::string& s);

    ~PairsdblibDebug();

    static std::string getPrefix();

    static bool showDebug(unsigned int v);

    static void setVerbosity( unsigned int v );
    
    enum {indent = 4 } ;
    std::string mName;
    unsigned int mPriority;
    std::ostream& mStream;
    static unsigned int mRecurse;
    static unsigned int mVerbosity;
  };

  
#define debug_func_cerr(ii) PairsdblibDebug debug(std::cerr, ii, __ASSERT_FUNCTION)
#define debug_func_cout(ii) PairsdblibDebug debug(std::cout, ii, __ASSERT_FUNCTION)
#define debug_prefix PairsdblibDebug::getPrefix()
#define debug_show(ii) PairsdblibDebug::showDebug(ii)
#define debug_cout(ii, text) {if (debug_show(ii)){std::cout << debug_prefix << text << std::endl;}}
#define debug_cout_start(ii, text) {if (debug_show(ii)){std::cout << debug_prefix << text;}}
#define debug_cout_add(ii, text) {if (debug_show(ii)){std::cout << text;}}  
#define debug_cerr(ii, text) {if (debug_show(ii)){std::cerr << debug_prefix << text << std::endl;}}
#define debug_cerr_start(ii, text) {if (debug_show(ii)){std::cerr << debug_prefix << text;}}
#define debug_cerr_add(ii, text) {if (debug_show(ii)){std::cerr << text;}}
#else
#define debug_func_cerr(ii)
#define debug_func_cout(ii)
#define debug_prefix
#define debug_show(ii)
#define debug_cout(ii, text)
#define debug_cerr_start(ii, text)
#define debug_cout_add(ii, text)
#define debug_cerr(ii, text)
#define debug_cerr_start(ii, text)
#define debug_cerr_add(ii, text)

#endif
  
}

#endif  //_ALIGNLIB_DEBUG_H

