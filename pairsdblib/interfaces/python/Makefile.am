# uses to big extent the Makefile.template supplied by SWIG

# Since this module depends also on definitions in the library AlignLib, SWIG
# has to use the runtime-libraries to resolve name conflicts. See 
# http://www.swig.org/Doc1.1/HTML/Advanced.html#n1 for an explanation.

# There are four libraries you need:
# 1. alignlib 
# 2. SWIG runtime library (contained in swig_libdir)

# 20.11.2002 added CXXFLAGS instead of CFLAGS

TARGET     = $(PACKAGE)
INTERFACE  = $(PACKAGE).i
SWIGOPT    = -shadow -noruntime -I$(SWIG_LIB) 		# build shadow classes, use runtime library 
CFLAGS     = -fexceptions				# This is needed for g++
LIBS       = -L$(top_srcdir)/src/.libs -l$(PACKAGE) \
	        $(alignlib_libdir) \
	        $(mysql_libdir) \
	        -L$(swig_libdir) -lswigpy 
INCLUDES   = -I$(top_srcdir) -I$(top_srcdir)/src $(PYTHON_INCLUDES) $(mysql_includedir) $(alignlib_includedir)

IWRAP      = $(INTERFACE:.i=_wrap.i)
ISRCS      = $(IWRAP:.i=.c)
ICXXSRCS   = $(IWRAP:.i=.cxx)
IOBJS      = $(IWRAP:.i=.o)
OBJS       = $(SRCS:.c=.o) $(CXXSRCS:.cxx=.o)
SO	   = .so

PYTHON_INTERFACE = $(INTERFACE:.i=.py)
CPYTHON_INTERFACE = $(INTERFACE:.i=.pyc)

install:
	$(INSTALL_PROGRAM) $(PYTHON_INTERFACE) $(PYTHON_SITE)
	$(INSTALL_PROGRAM) $(CPYTHON_INTERFACE) $(PYTHON_SITE)
	$(INSTALL_PROGRAM) _$(TARGET)$(SO) $(PYTHON_SITE_EXEC)

all:
	$(SWIG) -c++ -python $(SWIGOPT) $(INTERFACE)
	$(CXX) -c $(CCSHARED) $(CXXFLAGS) $(ICXXSRCS) $(SRCS) $(CXXSRCS) $(INCLUDES)
	$(CXX) -shared $(OBJS) $(IOBJS) $(LIBS) -o _$(TARGET)$(SO) 
	python -c 'import py_compile;py_compile.compile("$(PYTHON_INTERFACE)")'

static:
	$(MAKE) -f $(TOP)/Makefile.template SRCS='$(SRCS)' SWIG='$(SWIG)' SWIGOPT='$(SWIGOPT)' \
	CFLAGS='$(CFLAGS)' TARGET='mypython' INTERFACE='$(INTERFACE)' python_static_cpp

clean:
	rm -f *_wrap* *.o *~ *.so *.py mypython *.pyc .~*




