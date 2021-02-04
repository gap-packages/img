#############################################################################
##
#W  configure.ac                                            Laurent Bartholdi
##
#Y Copyright (C) 2009-2012, Laurent Bartholdi
##
#############################################################################

# autoconf
AC_PREREQ(2.68)
AC_INIT([img],m4_esyscmd([tr -d '\n' < VERSION]),[https://github.com/gap-packages/img/issues],[https://gap-packages.github.io/img])
AC_CONFIG_SRCDIR([src/img.c])
AC_CONFIG_AUX_DIR([cnf])
AC_CONFIG_HEADER([src/pkgconfig.h:cnf/pkgconfig.h.in])
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
AC_PROG_CC_STDC
AC_PROG_CXX
AX_CC_MAXOPT

# gap
AC_FIND_GAP

# C
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h])
AC_C_CONST
AC_C_INLINE
AC_CHECK_HEADERS(malloc/malloc.h)
AC_CHECK_HEADERS(malloc.h)
AC_CHECK_LIB(m,sincos,AC_DEFINE(HAVE_SINCOS,1,do we have sincos?))

AC_CYGWIN
AM_CONDITIONAL([SYS_IS_CYGWIN], [test "$CYGWIN" = "yes"])
if test "$CYGWIN" = "yes"; then
  AC_DEFINE(SYS_IS_CYGWIN32, 1, are we on CYGWIN?)
else
  AC_DEFINE(SYS_IS_CYGWIN32, 0, are we on CYGWIN?)
fi

# AC_CHECK_LEVMAR
AC_CHECK_CHOLMOD
# AX_BLAS # too heavy, requires F77
# AC_CHECK_LIBDOGLEG

################################################################
# generate files

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
