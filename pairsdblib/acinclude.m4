dnl ------------------------------------------------------------------------------------
dnl adl_SEARCH_LIB_PATH(FUNCTION, LIB [, PATH [, VARIABLE-TO-SET
dnl            [, ACTION-IF-NOT-FOUND [, OTHER-LIBRARIES]]]])
dnl Find the path for the library defining FUNC, if it's not already available.
AC_DEFUN([adl_SEARCH_LIB_PATH],
[
[AC_PREREQ(2.13)]
AC_CACHE_CHECK(
		[for library path for $2], 
		[ac_cv_ldflags_$2],
		[	ac_ldflags_search_save_LIBS="$LIBS"
			ac_cv_ldflags_$2="no"
			AC_TRY_LINK_FUNC([$1], [ac_cv_ldflags_$2="none required"])
			test "$ac_cv_ldflags_$2" = "no" && for i in $3; do
				echo Trying "$i"
				LIBS="-L$i -l$2 $6 $ac_ldflags_search_save_LIBS"
				AC_TRY_LINK_FUNC([$1],[ac_cv_ldflags_$2="-L$i -l$2" break])
			done
			LIBS="$ac_ldflags_search_save_LIBS"
		]
)

if test "$ac_cv_ldflags_$2" != "no"; then
	if test "$ac_cv_ldflags_$2" = "none required"; then
		$4=""
	else :
		LIBS="$ac_cv_ldflags_$2 $LIBS"
		$4="$ac_cv_ldflags_$2"
	fi
	AC_SUBST($4)
else :
	$5
fi
]
)  

dnl ------------------------------------------------------------------------------------
dnl adl_SEARCH_LIB_PATH_NOFUNC(HEADERS, PROGRAM BODY, LIB [, PATH [, VARIABLE-TO-SET
dnl            [, ACTION-IF-NOT-FOUND [, OTHER-LIBRARIES [, variable-name]]]]])
dnl variable-name: hack for mysql++, the ++ cause problems in filename
dnl Find the path for the library defining FUNC, if it's not already available.
dnl use a whole program to check if linking succeeds
AC_DEFUN([adl_SEARCH_LIB_PATH_NOFUNC],
[
AC_PREREQ([2.13])
AC_CACHE_CHECK(
	[for library path for $3], 
	[ac_cv_ldflags_$3],
	[	ac_ldflags_search_save_LIBS="$LIBS"
		ac_cv_ldflags_$3="no"

		dnl Try if library is needed
		AC_TRY_LINK([$1], [$2], [ac_cv_ldflags_$3="none required"])

		dnl Try if library is already on path
		LIBS="$7 $ac_ldflags_search_save_LIBS"	
		AC_TRY_LINK([$1], [$2], [ac_cv_ldflags_$3="none required"])

		test "$ac_cv_ldflags_$3" = "no" && for i in $4; do
			dnl echo Trying "$i"
			LIBS="-L$i -l$3 $7 $ac_ldflags_search_save_LIBS"
			AC_TRY_LINK([$1], [$2], [ac_cv_ldflags_$3="-L$i -l$3" break])
		done
		LIBS="$ac_ldflags_search_save_LIBS"
	]
)

if test "$ac_cv_ldflags_$3" != "no"; then
	if test "$ac_cv_ldflags_$3" = "none required"; then
		$5="" 
	else :
		LIBS="$ac_cv_ldflags_$3 $LIBS"
		$5="$ac_cv_ldflags_$3"
	fi
	AC_SUBST($5)
else :
  	$6
fi
]
)  


dnl ------------------------------------------------------------------------------------
dnl adl_SEARCH_LIB_PATH_NOFUNC_NOLIBS(HEADERS, PROGRAM BODY, LIB [, PATH [, VARIABLE-TO-SET
dnl            [, ACTION-IF-NOT-FOUND [, OTHER-LIBRARIES ]]]])
dnl Search for library in different paths. Put path into variable, but do not add to
dnl LIBS. You have to run AC_SUBST again for some reason.
AC_DEFUN([adl_SEARCH_LIB_PATH_NOFUNC_NOLIBS],
[
AC_PREREQ([2.13])
AC_CACHE_CHECK(
	[for library path for $3], 
	[ac_cv_ldflags_$3],
	[	ac_ldflags_search_save_LIBS="$LIBS"
		ac_cv_ldflags_$3="no"
		AC_TRY_LINK([$1], [$2], [ac_cv_ldflags_$3="none required"])
		test "$ac_cv_ldflags_$3" = "no" && for i in $4; do
			LIBS="-L$i -l$3 $7 $ac_ldflags_search_save_LIBS"
			AC_TRY_LINK([$1], [$2], [ac_cv_ldflags_$3="$i"])
		done
		LIBS="$ac_ldflags_search_save_LIBS"
	]
)

if test "$ac_cv_ldflags_$3" != "no"; then
	if test "$ac_cv_ldflags_$3" = "none required"; then
		$5="" 
	else :
		$5="$ac_cv_ldflags_$3"
	fi
	AC_SUBST($5)
else :
  	$6
fi
]
)  

dnl #########################################################################
dnl adl_SEARCH_INCLUDE_PATH(HEADER-FILE [, PATH [, VARIABLE-TO-SET
dnl            [, ACTION-IF-NOT-FOUND ]]])
dnl Find the path for the HEADER-FILE, if it's not already available.
AC_DEFUN([adl_SEARCH_INCLUDE_PATH],
[
AC_PREREQ([2.13])
ac_include_search_save_CPPFLAGS="$CPPFLAGS"
ac_safe=`echo "$1" | sed 'y%./+-%__p_%'`

eval "ac_cv_header_$ac_safe=\"no\""

dnl check whether header is in system headers
AC_CHECK_HEADER($1, eval "ac_cv_header_$ac_safe=\"none required\"")

dnl check whether header already found by current CPPFLAGS
AC_TRY_CPP(
	[#include <$1>],
	[eval "ac_cv_header_$ac_safe=\"none required\""],
	[]
)	

dnl check list of possible locations
eval "test \"`echo '$ac_cv_header_'$ac_safe`\" = \"no\"" && for i in $2; do
	AC_MSG_CHECKING([for $1 in $i])
	CPPFLAGS="-I$i $ac_include_search_save_CPPFLAGS"

	AC_TRY_CPP(
		[#include <$1>],
		[eval "ac_cv_header_$ac_safe=\"-I$i\""
		AC_MSG_RESULT("-I$i")
		break],
		AC_MSG_RESULT(no)
		)	
	done

CPPFLAGS="$ac_include_search_save_CPPFLAGS"
eval "ac_include_search_result=`echo \\\$ac_cv_header_$ac_safe`"

dnl append to CPPFLAGS
if test "$ac_include_search_result" != "no"; then
	if test "$ac_include_search_result" = "none required"; then
	 	$3=""
	else
		CPPFLAGS="$ac_include_search_result $CPPFLAGS"
	 	$3="$ac_include_search_result"
	fi
	AC_SUBST($3)
else :
  	$4
fi
])      

dnl #########################################################################
dnl Check command line option --enable-DEBUG
AC_DEFUN([adl_ENABLE_DEBUG],
[
AC_PROVIDE([$0])
 
ac_cv_enable_debug=no
 
dnl command line options and --help descriptions
AC_ARG_ENABLE( debug,
[  --enable-debug		compile with option DEBUG \[default=no\] ],
ac_cv_enable_debug=$enableval,	
ac_cv_enable_debug=no	
)

if test x"$ac_cv_enable_debug" = xyes; then
	CXXFLAGS="$CXXFLAGS -DDEBUG"
	AC_SUBST(CXXFLAGS)

	AC_MSG_CHECKING("Enable option DEBUG")
	AC_MSG_RESULT(yes)
fi
]
)
   
dnl
dnl
dnl AM_CHECK_PYMOD(MODNAME [,SYMBOL [,ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]]])
dnl Check if a module containing a given symbol is visible to python.
AC_DEFUN([AM_CHECK_PYMOD],
[AC_REQUIRE([AM_PATH_PYTHON])
py_mod_var=`echo $1['_']$2 | sed 'y%./+-%__p_%'`
AC_MSG_CHECKING(for ifelse([$2],[],,[$2 in ])python module $1)
AC_CACHE_VAL(py_cv_mod_$py_mod_var, [
ifelse([$2],[], [prog="
import sys
try:
        import $1
except ImportError:
        sys.exit(1)
except:
        sys.exit(0)
sys.exit(0)"], [prog="
import $1
$1.$2"])
if $PYTHON -c "$prog" 1>&AC_FD_CC 2>&AC_FD_CC
  then
    eval "py_cv_mod_$py_mod_var=yes"
  else
    eval "py_cv_mod_$py_mod_var=no"
  fi
])
py_val=`eval "echo \`echo '$py_cv_mod_'$py_mod_var\`"`
if test "x$py_val" != xno; then
  AC_MSG_RESULT(yes)
  ifelse([$3], [],, [$3
])dnl
else
  AC_MSG_RESULT(no)
  ifelse([$4], [],, [$4
])dnl
fi
])         

dnl a macro to check for ability to create python extensions
dnl  AM_CHECK_PYTHON_HEADERS([ACTION-IF-POSSIBLE], [ACTION-IF-NOT-POSSIBLE])
dnl function also defines PYTHON_INCLUDES
AC_DEFUN([AM_CHECK_PYTHON_HEADERS],
[AC_REQUIRE([AM_PATH_PYTHON])
AC_MSG_CHECKING(for headers required to compile python extensions)
dnl deduce PYTHON_INCLUDES
py_prefix=`$PYTHON -c "import sys; print sys.prefix"`
py_exec_prefix=`$PYTHON -c "import sys; print sys.exec_prefix"`
PYTHON_INCLUDES="-I${py_prefix}/include/python${PYTHON_VERSION}"
if test "$py_prefix" != "$py_exec_prefix"; then
  PYTHON_INCLUDES="$PYTHON_INCLUDES -I${py_exec_prefix}/include/python${PYTHON_VERSION}"
fi
AC_SUBST(PYTHON_INCLUDES)
dnl check if the headers exist:
save_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $PYTHON_INCLUDES"
AC_TRY_CPP([#include <Python.h>],dnl
[AC_MSG_RESULT(found)
$1],dnl
[AC_MSG_RESULT(not found)
$2])
CPPFLAGS="$save_CPPFLAGS"
])
 
# AM_PATH_PYTHON([package, module])
#
 
# Adds support for distributing Python modules or the special form
# of a module called a `package.'  Modules of the first type are
# files ending in `.py' with no '__init__.py' file.  This must be
# placed on the PYTHONPATH, and the default location is PYTHON_SITE,
# or $(prefix)/lib/python$(PYTHON_VERSION)/site-package
#
# A package module is contained in its own directory.  This directory
# is named PACKAGE, which was the name given to automake.  The full
# directory path is PYTHON_SITE_PACKAGE or
#   $(prefix)/lib/python$(PYTHON_VERSION)/site-package/$(PACKAGE)
# where site-package is on the PYTHONPATH.  The `__init__.py' file is
# located in the directory, along with any other submodules which may
# be necessary.
 
 
AC_DEFUN([AM_PATH_PYTHON],
 [
  dnl Find a version of Python.  I could check for python versions 1.4
  dnl or earlier, but the default installation locations changed from
  dnl $prefix/lib/site-python in 1.4 to $prefix/lib/python1.5/site-packages
  dnl in 1.5, and I don't want to maintain that logic.
 
  if test -z "$PYBIN"; then
	AC_CHECK_PROGS(PYTHON, $prefix/bin/python python python2.4 python2.3 python2.2 python2.1 python2.0 python1.6 python1.5 python1.4 python)
  else
	PYTHON="$PYBIN"	
  fi

  AC_MSG_CHECKING([local Python configuration])
 
  dnl Query Python for its version number.  Getting [:3] seems to be
  dnl the best way to do this; it's what "site.py" does in the standard
  dnl library.  Need to change quote character because of [:3]
 
  AC_SUBST(PYTHON_VERSION)
  changequote(<<, >>)dnl
  PYTHON_VERSION=`$PYTHON -c "import sys; print sys.version[:3]"`
  changequote([, ])dnl
 
  dnl Use the values of $prefix and $exec_prefix for the corresponding
  dnl values of PYTHON_PREFIX and PYTHON_EXEC_PREFIX.  These are made
  dnl distinct variables so they can be overridden if need be.  However,
  dnl general consensus is that you shouldn't need this ability.
  AC_MSG_CHECKING(for Python prefix)
  PYTHON_PREFIX=`($PYTHON -c "import sys; print sys.prefix") 2>/dev/null`
  AC_MSG_RESULT($PYTHON_PREFIX)
    
  AC_MSG_CHECKING(for Python exec-prefix)
  PYTHON_EXEC_PREFIX=`($PYTHON -c "import sys; print sys.exec_prefix") 2>/dev/null`
  AC_MSG_RESULT($PYTHON_EXEC_PREFIX)

  dnl At times (like when building shared libraries) you may want
  dnl to know which OS platform Python thinks this is.
 
  AC_SUBST(PYTHON_PLATFORM)
  PYTHON_PLATFORM=`$PYTHON -c "import sys; print sys.platform"`   
  dnl Set up 4 directories:
 
  dnl   pythondir -- location of the standard python libraries
  dnl     Also lets automake think PYTHON means something.
 
  AC_SUBST(pythondir)
  pythondir=$PYTHON_PREFIX"/lib/python"$PYTHON_VERSION
 
  dnl   PYTHON_SITE -- the platform independent site-packages directory
 
  AC_SUBST(PYTHON_SITE)
  PYTHON_SITE=$pythondir/site-packages
 
  dnl   PYTHON_SITE_PACKAGE -- the $PACKAGE directory under PYTHON_SITE
 
  AC_SUBST(PYTHON_SITE_PACKAGE)
  PYTHON_SITE_PACKAGE=$pythondir/site-packages/$PACKAGE
 
  dnl   PYTHON_SITE_EXEC -- platform dependent site-packages dir (eg, for
  dnl        shared libraries)
 
  AC_SUBST(PYTHON_SITE_EXEC)
  PYTHON_SITE_EXEC=$PYTHON_EXEC_PREFIX"/lib/python"$PYTHON_VERSION/site-packages
 
 
  dnl Configure PYTHON_SITE_INSTALL so it points to either the module
  dnl directory or the package subdirectory, depending on the $1
  dnl parameter ("module" or "package").
 
  AC_SUBST(PYTHON_SITE_INSTALL)
  ifelse($1, module, [PYTHON_SITE_INSTALL=$PYTHON_SITE],
         $1, package, [PYTHON_SITE_INSTALL=$PYTHON_SITE_PACKAGE],
   [errprint([Unknown option `$1' used in call to AM_PATH_PYTHON.
Valid options are `module' or `package'
])m4exit(4)])
 
  dnl All done.
 
  AC_MSG_RESULT([looks good])
])   

dnl Handle the --enable-html-doc option.
dnl
dnl The following snippet goes into doc/Makefile.am
dnl ------------------------
dnl CLEANFILES = yourfile.html
dnl
dnl if MAKE_HTML
dnl hdir = @htmldir@
dnl h_DATA = yourfile.html
dnl endif
dnl
dnl SUFFIXES = .html
dnl
dnl .texi.html:
dnl         @cd $(srcdir) && rm -f $@ $@-[0-9] $@-[0-9][0-9]
dnl         cd $(srcdir) && $(MAKEINFO) --html `echo $< | sed 's,.*/,,'`
dnl ------------------------
dnl
dnl Written by Alexandre Duret-Lutz <duret_g@epita.fr>
 
AC_DEFUN([adl_ENABLE_HTML_DOC],
 [AC_ARG_ENABLE([html-doc],
                [--enable-html-doc=DIR       build and install html documentation (install via make install-data)])
  if test "${enable_html_doc-no}" != no; then
    if test "x$enable_html_doc" = xyes; then
      htmldir="$prefix/doc/$PACKAGE"
    else
      htmldir="$enable_html_doc"
    fi
    AC_SUBST([htmldir])
  else
    htmldir='\<none\>'
  fi
  AM_CONDITIONAL(MAKE_HTML, [test "${enable_html_doc-no}" != no])])     


dnl #########################################################################
dnl Check for alignlib headers and libraries
AC_DEFUN([AC_PACKAGE_ALIGNLIB],
[
AC_PROVIDE([$0])
AC_REQUIRE_CPP

AC_LANG_SAVE
dnl Use C++ for compiling
AC_LANG_CPLUSPLUS
 
ac_cv_with_alignlib_dir=""
ac_cv_with_alignlib_lib=""
ac_cv_with_alignlib_inc=""

ALIGNLIB_CFLAGS=""
ALIGNLIB_LIBS=""

dnl command line options and --help descriptions
 AC_ARG_WITH(alignlib-dir,
    [  --with-alignlib-dir           where the root of alignlib is installed (/usr/local)],
    [  ac_cv_with_alignlib_dir="$withval" ]
)
AC_ARG_WITH(alignlib-include,
    [  --with-alignlib-include       where the alignlib headers are. (/usr/local/include) ],
    [  ac_cv_with_alignlib_inc="$withval" ]
)
AC_ARG_WITH(alignlib-lib,
    [  --with-alignlib-lib           where the alignlib library is installed. (/usr/local/lib)],
    [  ac_cv_with_alignlib_lib="$withval" ]
)
 
dnl 1. Check for include path:
adl_SEARCH_INCLUDE_PATH([alignlib.h],
[   $ac_cv_with_alignlib_inc \
    $ac_cv_with_alignlib_dir/include \
    /usr/include                  \
    /usr/local/include                     \
    /usr/local/alignlib/include            \
],
[ALIGNLIB_CFLAGS],
AC_MSG_ERROR([Can\'t find the alignlib headers])
)

adl_SEARCH_LIB_PATH_NOFUNC([ 
	#include <alignlib.h>
],
[ alignlib::makeSequence("AAA") ], 
alignlib,
[  $ac_cv_with_alignlib_dir/lib   \
 	   $ac_cv_with_alignlib_lib       \
	   /usr/local/lib \
    	   /usr/lib/alignlib \
	   /usr/local/lib/alignlib \
    	   /usr/local/alignlib/lib \
    	   /usr/local/alignlib/lib/alignlib \
    	   /usr/share/alignlib/lib \
    	   /usr/share/lib/alignlib \
],
[ALIGNLIB_LIBS],
AC_MSG_ERROR([Can\'t find libalignlib]),
-lm
)

dnl Substitute some variables for use with SWIG
alignlib_includedir="$ALIGNLIB_CFLAGS"
AC_SUBST(alignlib_includedir)
alignlib_libdir="$ALIGNLIB_LIBS"
AC_SUBST(alignlib_libdir)

LIBS="-lalignlib $LIBS"

dnl Restore the compile language
AC_LANG_RESTORE 

]
)

