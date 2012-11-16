# check for cholmod library
# sets the CHOLMOD_CFLAGS, CHOLMOD_LDFLAGS and CHOLMOD_LIBS appropriately

AC_DEFUN([AC_CHECK_CHOLMOD],[

cm_LIBS="$LIBS"

AC_ARG_WITH(cholmod,
 [  --with-cholmod=<location>
    Location at which the cholmod library, needed for layout, was installed.],
 [CHOLMOD_CFLAGS="-I$withval/include"; CHOLMODLIB="$withval/lib"]
)

AC_ARG_WITH(cholmod-include,
 [  --with-cholmod-include=<location>
    Location at which the cholmod include files were installed.],
 [CHOLMOD_CFLAGS="-I$withval"]
)

AC_ARG_WITH(cholmod-lib,
 [  --with-cholmod-lib=<location>
    Location at which the cholmod library files were installed.
 ],
 [CHOLMODLIB="$withval"]
)

AC_CHECK_HEADER(suitesparse/cholmod.h,[CHOLMOD_CFLAGS="-I/usr/include/suitesparse"])

if test "$CHOLMOD_CFLAGS" != ""; then
    CPPFLAGS="$CPPFLAGS $CHOLMOD_CFLAGS"
fi

AC_CHECK_HEADER(cholmod.h,,AC_MSG_ERROR([cholmod.h not found. Specify its location using --with-cholmod.]))

if test "$CHOLMODLIB" != ""; then
    LDFLAGS="$LDFLAGS -L$CHOLMODLIB"
    CHOLMOD_LDFLAGS="-L$CHOLMODLIB"
fi

CHOLMOD_LIBS=""

for cm_extra in -llapack -lcolamd -lsuitesparseconfig -lamd -lrt; do
    CHOLMOD_LIBS="$CHOLMOD_LIBS $cm_extra"
    AC_CHECK_LIB(cholmod,cholmod_allocate_triplet,[cm_found=true],[cm_found=false],[$CHOLMOD_LIBS])
    if test $cm_found = true; then break; fi
done

if test $cm_found = false; then
    AC_MSG_ERROR([libcholmod not found. Specify its location using --with-cholmod.])
fi

CHOLMOD_LIBS="-lcholmod $CHOLMOD_LIBS"

LIBS="$cm_LIBS"
AC_SUBST(CHOLMOD_CFLAGS)
AC_SUBST(CHOLMOD_LDFLAGS)
AC_SUBST(CHOLMOD_LIBS)
])
