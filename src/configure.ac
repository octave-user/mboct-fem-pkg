dnl  Copyright (C) 2018(-2020) Reinhard <octave-user@a1.net>

dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 3 of the License, or
dnl (at your option) any later version.

dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.

dnl You should have received a copy of the GNU General Public License
dnl along with this program; If not, see <http://www.gnu.org/licenses/>.

AC_PREREQ([2.69])
AC_INIT([mboct-fem-pkg], [0.1.0])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIRS([m4])

AM_INIT_AUTOMAKE

AC_PROG_CXX

AC_ARG_VAR([PKG_CONFIG],[Program to return meta information about installed libraries])
AC_ARG_VAR([NLOPT_INC],[nlopt include dirs])
AC_ARG_VAR([NLOPT_LIBS],[nlopt libraries])
AC_ARG_VAR([DEBUG],[Enable debug macros])
AC_ARG_VAR([OCTINCLUDEDIR],[Octave include directory])
AC_ARG_VAR([OCTLDFLAGS],[Octave linker flags])
AC_ARG_VAR([OCTLIBS],[Octave libraries])
AC_CHECK_PROG([PKG_CONFIG], [pkg-config], [pkg-config])
AC_CHECK_PROG([MKOCTFILE], [mkoctfile], [mkoctfile])
AC_CHECK_PROG([OCTAVE_CONFIG], [octave-config], [octave-config])

if test -z "${MKOCTFILE}"; then
  AC_MSG_ERROR([mkoctfile not found], 1);
fi

if test -z "${OCTAVE_CONFIG}"; then
  AC_MSG_ERROR([octave-config not found], 1);
fi

if test -z "${NLOPT_INC}"; then
  if ! test -z "${PKG_CONFIG}"; then
    NLOPT_INC="`${PKG_CONFIG} --cflags nlopt`"
  fi
fi

if test -z "${NLOPT_LIBS}"; then
  if ! test -z "${PKG_CONFIG}"; then
    NLOPT_LIBS="`${PKG_CONFIG} --libs nlopt`"
  else
    NLOPT_LIBS="-lnlopt"
  fi
fi

AC_LANG(C++)

have_win32_io=no

AC_CHECK_HEADERS([fcntl.h io.h],[have_fcntl=yes])

if test "${have_fcntl}" = "yes"; then
   AC_CHECK_FUNCS([setmode read write fileno],[have_setmode=yes])

   if test "${have_setmode}" = "yes"; then
      AC_CHECK_DECL([_O_BINARY],[have_binary=yes],[],[[#include <fcntl.h>]])

      if test "${have_binary}" = "yes"; then
         AC_DEFINE([HAVE_WIN32_IO],[1],[Define if Win32 IO functions can be used])
         have_win32_io=yes
      fi         
   fi
fi

have_unix_io=no

AC_CHECK_HEADERS([unistd.h],[have_unistd=yes])

if test "${have_unistd}" = "yes"; then
   AC_CHECK_FUNCS([read write fileno],[have_read=yes])

   if test "${have_read}" = "yes"; then
      AC_DEFINE([HAVE_UNIX_IO],[1],[Define if Unix IO functions can be used])
      have_unix_io=yes
   fi
fi

AM_CONDITIONAL([USE_MESH_SIZE],[test "${have_win32_io}" = "yes" || test "${have_unix_io}" = "yes"])

AC_CHECK_HEADERS([algorithm array cassert iostream memory stdexcept string vector], [], AC_MSG_ERROR([standard C++ headers not found], 1))

TCPPFLAGS="${CPPFLAGS}"
TLIBS="${LIBS}"

CPPFLAGS="${NLOPT_INC} ${CPPFLAGS}"
LIBS="${NLOPT_LIBS} ${LIBS}"

AC_LANG(C)
AC_CHECK_HEADER([nlopt.h])

have_nlopt=no

if test "${ac_cv_header_nlopt_h}" = "yes"; then
   AC_CHECK_TYPE([nlopt_result],[],[],[[#include <nlopt.h>]])

   if test "${ac_cv_type_nlopt_result}" = "yes"; then
      AC_CHECK_DECLS([NLOPT_FAILURE],[],[],[[#include <nlopt.h>]])

      if test "${ac_cv_have_decl_NLOPT_FAILURE}" = "yes"; then
            AC_CHECK_DECLS([NLOPT_LD_SLSQP],[],[],[[#include <nlopt.h>]])

            if test "${ac_cv_have_decl_NLOPT_LD_SLSQP}" = "yes"; then
               AC_CHECK_FUNCS([nlopt_create \
                               nlopt_destroy \
                               nlopt_set_lower_bounds \
                               nlopt_set_min_objective \
                               nlopt_add_equality_constraint \
                               nlopt_set_upper_bounds \
                               nlopt_set_maxeval \
                               nlopt_set_xtol_abs1 \
                               nlopt_set_ftol_abs \
                               nlopt_optimize],[have_nlopt=yes])
            fi
      fi
   fi
fi

if test "${have_nlopt}" = "yes"; then
   AC_DEFINE([HAVE_NLOPT],[1],[Define if nlopt can be used])
fi

AM_CONDITIONAL([USE_NLOPT],[test "${have_nlopt}" = "yes"])

CPPFLAGS="${TCPPFLAGS}"
LIBS="${TLIBS}"

OCTINCLUDEDIR="${OCTINCLUDEDIR:-`${MKOCTFILE} -p INCFLAGS`}"
OCTLIBDIR="${OCTLIBDIR:-`${MKOCTFILE} -p OCTLIBDIR`}"
OCTLIBS="${OCTLIBS:-`${MKOCTFILE} -p LIBOCTAVE` `${MKOCTFILE} -p LIBOCTINTERP`}"
CXX="${CXX:-`${MKOCTFILE} -p CXX`}"
CXXFLAGS="${CXXFLAGS} `${MKOCTFILE} -p CXXFLAGS`"
OCTLDFLAGS="${OCTLDFLAGS:--L${OCTLIBDIR} -Wl,-rpath=${OCTLIBDIR} -Wl,--no-undefined}"

TCPPFLAGS="${CPPFLAGS}"
TLDFLAGS="${LDFLAGS}"
TLIBS="${LIBS}"

CPPFLAGS="${OCTINCLUDEDIR} ${CPPFLAGS}"

AC_LANG(C++)
AC_CHECK_HEADERS(octave/oct.h octave/oct-map.h, [], AC_MSG_ERROR([octave header files not found], 1))

LDFLAGS="${OCTLDFLAGS} ${LDFLAGS}"
LIBS="${OCTLIBS} ${LIBS}"

dnl Force mkoctfile to link to Octave's libraries in order to avoid unresolved symbols at runtime
AC_MSG_CHECKING([for octave libraries])
AC_TRY_LINK([
  #include <octave/oct.h>
],[
  octave_value_list ovl;
  ovl.append(octave_value(1.5));
],[use_octlibs=yes],[use_octlibs=no])

AC_MSG_RESULT([${use_octlibs}])

AM_CONDITIONAL(USE_OCTLIBS, [test "${use_octlibs}" = "yes"])

## Simple symbol alternatives of different Octave versions.
OF_OCTAVE_LIST_ALT_SYMS([

[dnl
  [is_integer_type],
  [isinteger],
  [[octave_value ().isinteger ();]],
  [OV_ISINTEGER],
  [],
  []
],

[dnl
  [is_real_type],
  [isreal],
  [[octave_value ().isreal ();]],
  [OV_ISREAL],
  [],
  []
]

],
   [oct-alt-includes.h])

CPPFLAGS="${TCPPFLAGS}"
LIBS="${TLIBS}"
LDFLAGS="${TLDFLAGS}"

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
