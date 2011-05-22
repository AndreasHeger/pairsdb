# Copyright (C) 2000 Andreas Heger <heger@ebi.ac.uk>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# using htmldir directly does not work
# all: documentation

if MAKE_HTML
hdir = @htmldir@
else
hdir=
endif

SOURCE=$(wildcard $(top_srcdir)/src/*.h) $(wildcard $(top_srcdir)/src/*.cpp) 
PYTHONSOURCES=$(wildcard $(top_srcdir)/interfaces/python/*.i)

## kept here for old times sake
kdoc: 
	if test $(KDOC); then \
		rm -rf ./html; \
		mkdir ./html; \
		$(KDOC) -d ./html $(top_srcdir)/src/*.h $(top_srcdir)/src/*.cpp ; \
	fi

## use doxygen for building your documentation
doxygen: $(SOURCE)
	if test $(DOXYGEN); then \
		$(DOXYGEN) doxygen.cfg >& $@; \
	fi

happydoc: $(PYTHONSOURCES)
	if test $(HAPPYDOC); then \
		$(HAPPYDOC) -d python $(top_srcdir)/interfaces/python/*.py >& $@; \
	fi

html: doxygen happydoc

install-data-hook:
	if test ! -z "$(hdir)"; then \
		if test -d ./html; then \
			if ! test -e $(hdir); then mkdir $(hdir); fi; \
			cp -dpR html/* $(hdir); \
		fi; \
		if test -d ./python; then \
			if ! test -e $(hdir)/python; then mkdir $(hdir)/python; fi; \
			cp -dpR python/* $(hdir)/python; \
		fi; \
	fi

clean: 
	rm -rf ./html
	rm -f happydoc doxygen











