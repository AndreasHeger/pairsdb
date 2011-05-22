/*
  alignlib - a library for aligning protein sequences

  $Id: alignlib.h,v 1.5 2005/02/24 11:07:25 aheger Exp $

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

/** Declaration of types */

#ifndef PAIRSDBLIB_DECLS_H_
#define PAIRSDBLIB_DECLS_H_

#include <alignlib_fwd.h>

namespace pairsdblib
{

	/** type of nids */
	typedef long Nid;
	typedef const char * Did;
	typedef double Pid;
	typedef double Score;
	typedef double ZScore;
	typedef long RepeatId;
	typedef const char * Identifier;
	typedef const char * Description;
	typedef const char * Sequence;
	typedef double EValue;
	typedef int Filter;
	typedef int Length;
	typedef int Method;
	typedef const char * InformationString;

	#define DEFAULT_UPPER_FILTER 100
	#define DEFAULT_LOWER_FILTER 0
	#define DEFAULT_MAX_EVALUE 10.0
	#define DEFAULT_MAX_LINES 0

	/** maximum length of a string buffer for users, passwords, etc. */
	#define MAX_BUFFER_LENGTH 100
	
	/* offset correction - subtract one from first position to translate
	 * pairsdb coordinates (1-based, inclusive) to alignlib coordinates 
	 * (0-based, inclusive/exclusive)
	 */
#define OFFSET_CORRECTION 1


} 

#endif /*PAIRSDBLIB_DECLS_H_*/
