# check for levmar library
# sets the LEVMAR_CFLAGS, LEVMAR_LDFLAGS and LEVMAR_LIBS appropriately

AC_DEFUN([AC_CHECK_LEVMAR],[
lm_LIBS="$LIBS"

LEVMAR=yes
LEVMAR_CFLAGS=""
LEVMAR_LDFLAGS=""

AC_ARG_WITH(levmar,
 [  --with-levmar=<location>
    Location at which the levmar library, needed for layout, was installed.],
 [if test "$withval" = extern -o "$withval" = yes -o "$withval" = no; then
    LEVMAR="$withval"
  else
    LEVMAR_CFLAGS="-I$withval/include"; LEVMAR_LDFLAGS="-L$withval/lib"
  fi]
)

AC_ARG_WITH(levmar-include,
 [  --with-levmar-include=<location>
    Location at which the levmar include files were installed.],
 [LEVMAR_CFLAGS="-I$withval"]
)

AC_ARG_WITH(levmar-lib,
 [  --with-levmar-lib=<location>
    Location at which the levmar library files were installed.
 ],
 [LEVMAR_LDFLAGS="-L$withval"]
)

LEVMAR_LIBS="-llevmar -llapack -lblas"

if test "$LEVMAR" = extern; then

AC_MSG_CHECKING([for levmar])
AC_MSG_RESULT([extern])
LEVMAR_MAKELIB=`printf 'liblevmar:
	mkdir -p $(EXTERN)/include $(EXTERN)/lib
	if [[ ! -f $(EXTERN)/include/levmar.h ]]; then \\
		cmake $(LEVMAR); \\
		$(MAKE) -C $(LEVMAR) levmar C_FLAGS=-fPIC; \\
		cp $(LEVMAR)/liblevmar.a $(EXTERN)/lib/; \\
		cp $(LEVMAR)/levmar.h $(EXTERN)/include/; \\
	fi\n'`

MAKE_LIBTARGETS="$MAKE_LIBTARGETS liblevmar"

LEVMAR_CFLAGS='-I$(EXTERN)/include'
LEVMAR_LDFLAGS='-L$(EXTERN)/lib'

elif test "$LEVMAR" != no; then

AC_LANG_PUSH([C])

lm_CPPFLAGS=$CPPFLAGS
CPPFLAGS="$CPPFLAGS $LEVMAR_CFLAGS"
AC_CHECK_HEADER(levmar.h,,AC_MSG_ERROR([levmar.h not found. Specify its location using --with-levmar.
The package may be downloaded from http://www.ics.forth.gr/~lourakis/levmar/]))
CPPFLAGS=$lm_CPPFLAGS

lm_LDFLAGS=$LDFLAGS
LDFLAGS="$LDFLAGS $LEVMAR_LDFLAGS"
AC_CHECK_LIB(levmar,dlevmar_dif,,AC_MSG_ERROR([liblevmar not found. Specify its location using --with-levmar.]),[-llapack -lblas -lm])
LDFLAGS=$lm_LDFLAGS

AC_LANG_POP([C])

fi

LIBS="$lm_LIBS"
AC_SUBST(LEVMAR_CFLAGS)
AC_SUBST(LEVMAR_LDFLAGS)
AC_SUBST(LEVMAR_LIBS)
AC_SUBST(LEVMAR_MAKELIB)
AC_SUBST(MAKE_LIBTARGETS)
])
