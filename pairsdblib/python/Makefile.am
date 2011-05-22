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

EXTRA_DIST=setup.py build.py includes.h pairsdblib.cpp modules

python:
	python setup.py build

install-exec-hook:
	python setup.py install --prefix=${prefix}
