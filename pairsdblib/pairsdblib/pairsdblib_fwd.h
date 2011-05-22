/*
  alignlib - a library for aligning protein sequences

  $Id$

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


#ifndef PAIRSDBLIB_FWD_H_
#define PAIRSDBLIB_FWD_H_

#include "pairsdblib_types.h"
#include <boost/shared_ptr.hpp>
#include <vector>
namespace pairsdblib
{
	/** actor objects and their handles */
	class Query;
	typedef boost::shared_ptr<Query>HQuery;
	
	class Connection;
	typedef boost::shared_ptr<Connection>HConnection;
	
	class SQLQuery;
	typedef boost::shared_ptr<SQLQuery>HSQLQuery;
	
	class Row;
	typedef boost::shared_ptr<Row>HRow;
	    
	struct NeighborLink;
	typedef std::vector<NeighborLink>NeighborLinkVector;
	typedef boost::shared_ptr<NeighborLinkVector>HNeighborLinkVector;

	typedef  std::map<Nid, int>MapNid2Length;
	typedef boost::shared_ptr<MapNid2Length>HMapNid2Length;

	typedef  std::vector<Nid>NidVector;
	typedef boost::shared_ptr<NidVector>HNidVector;
}



#endif /*PAIRSDBLIB_FWD_H_*/
