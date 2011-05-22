/*
  pairsdblib - a library for aligning protein sequences

  $Id: PairsdblibException.cpp,v 1.2 2004/01/07 14:35:31 aheger Exp $

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



#include <iostream>
#include "PairsdblibException.h"
#include "PairsdblibDebug.h"

using namespace std;

namespace pairsdblib 
{

PairsdblibException::PairsdblibException( const char *msg ) :
	mMessage(msg)
{
	debug_func_cerr(5);
	std::cerr << mMessage;
}

PairsdblibException::~PairsdblibException() throw()
{
};

const char * PairsdblibException::what() const throw()
{
	debug_func_cerr(5);
	return mMessage.c_str();
}

} // namespace pairsdblib
