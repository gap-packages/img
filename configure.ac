#############################################################################
##
#W  configure.ac                                            Laurent Bartholdi
##
#Y Copyright (C) 2009-2012, Laurent Bartholdi
##
#############################################################################

# autoconf
AC_PREREQ([2.68])
AC_INIT([img],[package],[https://github.com/gap-packages/img/issues],[img],[https://gap-packages.github.io/img])
AC_CONFIG_SRCDIR([src/img.c])
AC_CONFIG_AUX_DIR([cnf])
AC_CONFIG_HEADERS([src/pkgconfig.h:cnf/pkgconfig.h.in])
AC_CONFIG_MACRO_DIR([m4])

# automake
AM_INIT_AUTOMAKE([1.11 -Wall foreign subdir-objects])
## AM_SILENT_RULES([yes])
AM_MAINTAINER_MODE
AM_PROG_AR

# libtool
LT_PREREQ([2.4.2])
LT_INIT([disable-static dlopen win32-dll])

# compilers
AC_PROG_CC
AC_PROG_CXX
AC_CANONICAL_HOST
case $host_os in
  *cygwin* ) CYGWIN=yes;;
	 * ) CYGWIN=no;;
esac
AM_CONDITIONAL([SYS_IS_CYGWIN], [test "$CYGWIN" = "yes"])

# gap
FIND_GAP

AC_CHECK_LIB(m,sincos,AC_DEFINE(HAVE_SINCOS,1,do we have sincos?))

# AC_CHECK_LEVMAR
AC_CHECK_CHOLMOD
# AX_BLAS # too heavy, requires F77
# AC_CHECK_LIBDOGLEG

################################################################
# generate files

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
