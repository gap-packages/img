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
AC_CHECK_LIB(cholmod,cholmod_allocate_triplet,,AC_MSG_ERROR([libcholmod not found. Specify its location using --with-cholmod.]),[-llapack -lcolamd -lsuitesparseconfig -lamd -lrt])

CHOLMOD_LIBS="-lcholmod -llapack -lcolamd -lsuitesparseconfig -lamd -lrt"

LIBS="$cm_LIBS"
AC_SUBST(CHOLMOD_CFLAGS)
AC_SUBST(CHOLMOD_LDFLAGS)
AC_SUBST(CHOLMOD_LIBS)
])
