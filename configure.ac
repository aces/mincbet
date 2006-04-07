dnl Process this file with autoconf to produce a configure script.

# Require autoconf 2.57 or newer.
AC_PREREQ([2.57])

# The arguments are package name, and package version.
AC_INIT([mincbet],[1.0], [Andrew Janke <a.janke@gmail.com>])
AC_CONFIG_SRCDIR([bet.c])

AM_INIT_AUTOMAKE

# The argument is the name of the generated header file.
# It is recommended to leave it as "config.h".
AC_CONFIG_HEADERS([config.h])

CFLAGS="$CFLAGS -DFDT=\"float\""

# This macro eases building the package; see m4/README.
smr_WITH_BUILD_PATH

AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

mni_REQUIRE_BICPL

dnl Build only static libs by default
AC_DISABLE_SHARED
AC_PROG_LIBTOOL

AC_CONFIG_FILES([Makefile])
AC_OUTPUT(
epm-header)

