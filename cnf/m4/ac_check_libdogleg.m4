# check for libdogleg library
# sets the LIBDOGLEG_CFLAGS, LIBDOGLEG_LDFLAGS and LIBDOGLEG_LIBS appropriately
# set also targets MAKE_LIBTARGETS and LIBDOGLEG_MAKELIB

AC_DEFUN([AC_CHECK_LIBDOGLEG],[
temp_LIBS="$LIBS"
temp_CPPFLAGS="$CPPFLAGS"
temp_LDFLAGS="$LDFLAGS"
LIBDOGLEG=unknown
LIBDOGLEG_CFLAGS=""
LIBDOGLEG_LDFLAGS=""
LIBDOGLEG_LIBS=""

AC_CHECK_CHOLMOD

AC_ARG_WITH(libdogleg,
 [  --with-libdogleg=<location>
    Location at which the libdogleg library, needed for layout, was installed.],
 [if test "$withval" = extern -o "$withval" = yes -o "$withval" = no; then
    LIBDOGLEG="$withval"
  else
    LIBDOGLEG=yes
    LIBDOGLEG_CFLAGS="-I$withval/include"; LIBDOGLEG_LDFLAGS="-L$withval/lib"
  fi]
)

AC_ARG_WITH(libdogleg-include,
 [  --with-libdogleg-include=<location>
    Location at which the libdogleg include files were installed.],
 [LIBDOGLEG=yes; LIBDOGLEG_CFLAGS="-I$withval"]
)

AC_ARG_WITH(libdogleg-lib,
 [  --with-libdogleg-lib=<location>
    Location at which the libdogleg library files were installed.
 ],
 [LIBDOGLEG=yes; LIBDOGLEG_LDFLAGS="-L$withval"]
)

if test "$LIBDOGLEG" != no; then

LIBDOGLEG_LIBS="-lcholmod -ldogleg"

if test "$LIBDOGLEG" != extern; then

AC_LANG_PUSH([C])
temp_status=true
CPPFLAGS="$CPPFLAGS $LIBDOGLEG_CFLAGS"
AC_CHECK_HEADER(dogleg.h,,[temp_status=false])
LDFLAGS="$LDFLAGS $LIBDOGLEG_LDFLAGS"
AC_CHECK_LIB(dogleg,dogleg_optimize,,[temp_status=false],[],[-lcholmod])
AC_LANG_POP([C])

if test "$temp_status" = false; then
    if test "$LIBDOGLEG" = yes; then
        AC_MSG_ERROR([dogleg.h not found. Using --with-libdogleg, specify its location, "extern" to compile it locally, or "no" to disable it.
The package may be downloaded from https://github.com/Oblong/libdogleg.git])
    else
	LIBDOGLEG=extern
    fi
else
    LIBDOGLEG=yes
fi

fi

if test "$LIBDOGLEG" = extern; then

AC_MSG_NOTICE([I will compile libdogleg for you from the extern/ directory])

LIBDOGLEG_MAKELIB=`printf 'libdogleg:
	mkdir -p "$(EXTERN)/include" "$(EXTERN)/lib"
	if [[ ! -f "$(EXTERN)/include/dogleg.h" ]]; then \\
		$(MAKE) CFLAGS="$(CHOLMOD_CFLAGS)" -B -C $(LIBDOGLEG) && \\
		cp $(LIBDOGLEG)/libdogleg.a $(EXTERN)/lib/ && \\
		cp $(LIBDOGLEG)/dogleg.h $(EXTERN)/include/ && \\
		if [[ "$(shell uname)" = Darwin ]]; then install_name_tool -id $(EXTERN)/lib/libdogleg.dylib $(EXTERN)/lib/libdogleg.dylib; fi; \\
	fi\n'`

MAKE_LIBTARGETS="$MAKE_LIBTARGETS libdogleg"

LIBDOGLEG_CFLAGS='-I$(EXTERN)/include'
LIBDOGLEG_LDFLAGS='-L$(EXTERN)/lib'

fi

fi

CPPFLAGS="$temp_CPPFLAGS"
LDFLAGS="$temp_LDFLAGS"
LIBS="$temp_LIBS"

if test "$LIBDOGLEG" != no; then
    AC_DEFINE(USE_LIBDOGLEG)
fi

AC_SUBST(LIBDOGLEG_CFLAGS)
AC_SUBST(LIBDOGLEG_LDFLAGS)
AC_SUBST(LIBDOGLEG_LIBS)
AC_SUBST(LIBDOGLEG_MAKELIB)
AC_SUBST(MAKE_LIBTARGETS)
])
