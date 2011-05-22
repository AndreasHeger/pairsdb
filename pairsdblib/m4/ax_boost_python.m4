##### http://autoconf-archive.cryp.to/ax_boost_python.html
#
# SYNOPSIS
#
#   AX_BOOST_PYTHON
#
# DESCRIPTION
#
#   This macro checks to see if the Boost.Python library is installed.
#   It also attempts to guess the currect library name using several
#   attempts. It tries to build the library name using a user supplied
#   name or suffix and then just the raw library.
#
#   If the library is found, HAVE_BOOST_PYTHON is defined and
#   BOOST_PYTHON_LIB is set to the name of the library.
#
#   This macro calls AC_SUBST(BOOST_PYTHON_LIB).
#
#   In order to ensure that the Python headers are specified on the
#   include path, this macro requires AX_PYTHON to be called.
#
# LAST MODIFICATION
#
#   2007-07-24
#
# COPYLEFT
#
#   Copyright (c) 2007 Michael Tindal
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License as
#   published by the Free Software Foundation; either version 2 of the
#   License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#   General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
#   02111-1307, USA.
#
#   As a special exception, the respective Autoconf Macro's copyright
#   owner gives unlimited permission to copy, distribute and modify the
#   configure scripts that are the output of Autoconf when processing
#   the Macro. You need not follow the terms of the GNU General Public
#   License when using or distributing such scripts, even though
#   portions of the text of the Macro appear in them. The GNU General
#   Public License (GPL) does govern all other use of the material that
#   constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the
#   Autoconf Macro released by the Autoconf Macro Archive. When you
#   make and distribute a modified version of the Autoconf Macro, you
#   may extend this special exception to the GPL to apply to your
#   modified version as well.

AC_DEFUN([AX_BOOST_PYTHON],
[AC_REQUIRE([AC_PYTHON_DEVEL])dnl
AC_CACHE_CHECK(whether the Boost::Python library is available,
ac_cv_boost_python,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 CPPFLAGS_SAVE=$CPPFLAGS
 LDFLAGS_SAVE=$LDFLAGS
 if test ! -z "$PYTHON_LDFLAGS"; then
   LDFLAGS=$PYTHON_LDFLAGS $LDFLAGS
 fi
 if test x$PYTHON_CPPFLAGS != x; then
   CPPFLAGS=$PYTHON_CPPFLAGS $CPPFLAGS
 fi
 AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[
 #include <boost/python/module.hpp>
 using namespace boost::python;
 BOOST_PYTHON_MODULE(test) { throw "Boost::Python test."; }]],
 			   [[return 0;]]),
  			   ac_cv_boost_python=yes, ac_cv_boost_python=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_boost_python" = "yes"; then
  AC_DEFINE(HAVE_BOOST_PYTHON,,[define if the Boost::Python library is available])
  ax_python_lib=boost_python
  AC_ARG_WITH([boost-python],AS_HELP_STRING([--with-boost-python],[specify the boost python library or suffix to use]),
  [if test "x$with_boost_python" != "xno"; then
     ax_python_lib=$with_boost_python
     ax_boost_python_lib=boost_python-$with_boost_python
   fi])
  BOOSTLIBDIR=`echo $BOOST_LDFLAGS | sed -e 's/@<:@^\/@:>@*//'`
  for libextension in `ls $BOOSTLIBDIR/libboost_python*.so* 2>/dev/null  | sed 's,.*/,,' | sed -e 's;^lib\(boost_python.*\)\.so.*$;\1;'` ; do
      ax_lib=${libextension}
      AC_CHECK_LIB($ax_lib, exit,[BOOST_PYTHON_LIB="$ax_lib"; break],[link_python="no"],[$PYTHON_LDFLAGS $PYTHON_EXTRA_LIBS])
  done
  AC_SUBST(BOOST_PYTHON_LIB)
fi
CPPFLAGS=$CPPFLAGS_SAVE
LDFLAGS=$LDFLAGS_SAVE
])dnl