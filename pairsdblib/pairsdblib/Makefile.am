# Copyright (C) 2000 Andreas Heger <heger@ebi.ac.uk>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

INCLUDES = -I$(top_srcdir)/alignlib 

AM_CPPFLAGS=$(MYSQL_CFLAGS)
AM_LDFLAGS=$(MYSQL_LDFLAGS)

include_HEADERS = pairsdblib.h

pkginclude_HEADERS = pairsdblib.h \
			pairsdblib_fwd.h \
			pairsdblib_types.h \
			pairsdblib_functions.h \
			pairsdblib_interfaces.h \
			PairsdblibDebug.h \
			Connection.h \
			Row.h \
			Query.h \
			HelpersAlignatum.h \
			HelpersAlignment.h \
			HelpersAlignandum.h \
			HelpersMultipleAlignment.h \
			PairsdblibException.h \
			ImplAlignatumGroupie.h \
			ImplAlignatumNeighbour.h \
			ImplAlignatumRadar.h \
			ImplAlignatumDomain.h \
			ImplAlignatumNid.h \
			SQLQueries.h

PARTS = Connection.cpp Row.cpp Query.cpp \
		PairsdblibDebug.cpp \
		PairsdblibException.cpp \
		ImplAlignatumGroupie.cpp \
		ImplAlignatumNeighbour.cpp \
		ImplAlignatumRadar.cpp \
		ImplAlignatumDomain.cpp \
		ImplAlignatumNid.cpp \
		HelpersMultipleAlignment.cpp \
		HelpersAlignment.cpp \
		HelpersAlignatum.cpp \
		HelpersAlignandum.cpp \
		SQLQueries.cpp

lib_LTLIBRARIES = libpairsdblib.la

libpairsdblib_la_SOURCES = $(PARTS)

LIBOBJS= $(PARTS)
