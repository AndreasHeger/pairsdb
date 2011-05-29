# Copyright (C) 2000 Andreas Heger <heger@ebi.ac.uk>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

INCLUDES = -I$(top_srcdir)/pairsdblib

AM_CPPFLAGS=$(MYSQL_CFLAGS) $(BOOST_CPPFLAGS)
AM_LDFLAGS=$(MYSQL_LDFLAGS) -L../pairsdblib/.libs $(BOOST_PROGRAM_OPTIONS_LIB)
AM_CXXFLAGS=$(BOOST_CPPFLAGS)

LDADD=../pairsdblib/.libs/libpairsdblib.a

bin_PROGRAMS =  align_groupies \
		get_mali \
		rsdb_filter \
		rsdb_create \
		rsdb_index

align_groupies_SOURCES = align_groupies.cpp 
rsdb_create_SOURCES = rsdb_create.cpp
rsdb_filter_SOURCES = rsdb_filter.cpp
get_mali_SOURCES = get_mali.cpp
rsdb_index_SOURCES = rsdb_index.cpp



