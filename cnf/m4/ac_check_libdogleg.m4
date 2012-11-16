# check for libdogleg library
# sets the LIBDOGLEG_CFLAGS, LIBDOGLEG_LDFLAGS and LIBDOGLEG_LIBS appropriately

AC_DEFUN([AC_CHECK_LIBDOGLEG],[

extern_libs=false
AC_ARG_WITH(libdogleg,
 [  --with-libdogleg=<location>
    Location at which the libdogleg library, needed for layout, was installed.],
 [if test "$withval" = extern; then
    extern_libs=true
  else
    LIBDOGLEG_CFLAGS="-I$withval/include"; LIBDOGLEG_LDFLAGS="-L$withval/lib"
  fi]
)

AC_ARG_WITH(libdogleg-include,
 [  --with-libdogleg-include=<location>
    Location at which the libdogleg include files were installed.],
 [LIBDOGLEG_CFLAGS="-I$withval"]
)

AC_ARG_WITH(libdogleg-lib,
 [  --with-libdogleg-lib=<location>
    Location at which the libdogleg library files were installed.
 ],
 [LIBDOGLEG_LDFLAGS="-L$withval"]
)

LIBDOGLEG_LIBS="-lcholmod -ldogleg"

if test "$extern_libs" = true; then
AC_MSG_CHECKING([for libdogleg])
AC_MSG_RESULT([extern])
LIBDOGLEG_MAKELIB=`printf 'libdogleg:
	mkdir -p $(EXTERN)/include $(EXTERN)/lib
	if [[ ! -f $(EXTERN)/include/dogleg.h ]]; then \\
		$(MAKE) -B -C $(LIBDOGLEG); \\
		cp $(LIBDOGLEG)/libdogleg.* $(EXTERN)/lib/; \\
		cp $(LIBDOGLEG)/dogleg.h $(EXTERN)/include/; \\
	fi\n'`

MAKE_LIBTARGETS="$MAKE_LIBTARGETS libdogleg"
LIBDOGLEG_CFLAGS='-I$(EXTERN)/include'
LIBDOGLEG_LDFLAGS='-L$(EXTERN)/lib'

else

AC_LANG_PUSH([C])

ldl_CPPFLAGS=$CPPFLAGS
CPPFLAGS="$CPPFLAGS $LIBDOGLEG_CFLAGS"
AC_CHECK_HEADER(dogleg.h,,AC_MSG_ERROR([dogleg.h not found. Specify its location using --with-libdogleg.
The package may be downloaded from https://github.com/Oblong/libdogleg.git]))
CPPFLAGS=$ldl_CPPFLAGS

ldl_LDFLAGS=$LDFLAGS
LDFLAGS="$LDFLAGS $LIBDOGLEG_LDFLAGS"
AC_CHECK_LIB(dogleg,dogleg_optimize,,AC_MSG_ERROR([libdogleg not found. Specify its location using --with-libdogleg.]))
LDFLAGS=$ldl_LDFLAGS

AC_LANG_POP([C])

fi

AC_SUBST(LIBDOGLEG_CFLAGS)
AC_SUBST(LIBDOGLEG_LDFLAGS)
AC_SUBST(LIBDOGLEG_LIBS)
AC_SUBST(LIBDOGLEG_MAKELIB)
AC_SUBST(MAKE_LIBTARGETS)
])
