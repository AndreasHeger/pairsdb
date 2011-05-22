# Copyright (C) 2000 Andreas Heger <heger@ebi.ac.uk>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

EXTRA_DIST=data/README data/pairsdblib_data.tgz data/test_rsdb_create.sh

check_PROGRAMS = test_Connection \
		 test_Query \
		 test_Sequence \
		 test_FillAlignata \
		 test_FillProfile \
		 test_MultipleAlignment 

AM_CPPFLAGS=$(MYSQL_CFLAGS)
AM_LDFLAGS=$(MYSQL_LDFLAGS)

TESTS = $(check_PROGRAMS)

LDADD = $(top_srcdir)/pairsdblib/.libs/libpairsdblib.a -lm 

test_FillAlignata_SOURCES = test_FillAlignata.cpp 
test_FillProfile_SOURCES = test_FillProfile.cpp 
test_MultipleAlignment_SOURCES = test_MultipleAlignment.cpp 
test_Sequence_SOURCES = test_Sequence.cpp 
test_Connection_SOURCES = test_Connection.cpp 
test_Query_SOURCES = test_Query.cpp 

INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/pairsdblib 
