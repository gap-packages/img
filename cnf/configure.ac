#############################################################################
##
#W  configure.ac                                            Laurent Bartholdi
##
#Y Copyright (C) 2009-2012, Laurent Bartholdi
##
#############################################################################

AC_PREREQ(2.68)
LT_PREREQ([2.4.2])
AC_INIT(img,,laurent.bartholdi@gmail.com)
AC_CONFIG_SRCDIR([src/img_dll.c])
AC_CONFIG_AUX_DIR([cnf])
AC_CONFIG_HEADERS([src/config.h:cnf/config.h.in])
LT_INIT([disable-static dlopen win32-dll])

AC_PREFIX_DEFAULT([${PWD}])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_JAVAC
AC_CANONICAL_TARGET
AX_CC_MAXOPT

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

# Checks for standard C++ stuff
AC_CHECK_HEADERS(malloc/malloc.h)
AC_CHECK_HEADERS(malloc.h)

AC_CHECK_LIB(m,sincos,AC_DEFINE(HAVE_SINCOS,1,do we have sincos?))

# Locates GAP
AC_FIND_GAP

AC_SYS_IS_CYGWIN
AC_SYS_IS_DARWIN

# External libraries configuration
EXTERN="\$(CURDIR)/bin/$TARGET/extern"
MAKE_LIBTARGETS=""

AC_CHECK_LEVMAR
AC_CHECK_CHOLMOD
AC_CHECK_LIBDOGLEG

# External programs configuration

AC_PATH_PROG(APPLETVIEWER,[appletviewer])

if test -z "$APPLETVIEWER"; then
    AC_WARN([Could not find 'appletviewer' (debian package java-6-sdk)... you won't be able to draw spiders])
fi

################################################################
# generate files

mkdir -p bin/$GAPARCH
CONFIG_STATUS=bin/$GAPARCH/config.status

AC_CONFIG_FILES([$GAP_MAKEFILE:cnf/Makefile.in bin/$GAPARCH/files.g:cnf/files.g.in])

if test "$GAP_MAKEFILE" != Makefile; then
    ln -sf "$GAP_MAKEFILE" Makefile
fi

AC_OUTPUT
