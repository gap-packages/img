# check for levmar library
# sets the LEVMAR_CFLAGS, LEVMAR_LDFLAGS and LEVMAR_LIBS appropriately

AC_DEFUN([AC_CHECK_LEVMAR],[
temp_LIBS="$LIBS"
temp_CPPFLAGS="$CPPFLAGS"
temp_LDFLAGS="$LDFLAGS"
LEVMAR=unknown
LEVMAR_CFLAGS=""
LEVMAR_LDFLAGS=""
LEVMAR_LIBS=""

AC_ARG_WITH(levmar,
 [  --with-levmar=<location>
    Location at which the levmar library, needed for layout, was installed.],
 [if test "$withval" = extern -o "$withval" = yes -o "$withval" = no; then
    LEVMAR="$withval"
  else
    LEVMAR=yes
    LEVMAR_CFLAGS="-I$withval/include"; LEVMAR_LDFLAGS="-L$withval/lib"
  fi]
)

AC_ARG_WITH(levmar-include,
 [  --with-levmar-include=<location>
    Location at which the levmar include files were installed.],
 [LEVMAR=yes; LEVMAR_CFLAGS="-I$withval"]
)

AC_ARG_WITH(levmar-lib,
 [  --with-levmar-lib=<location>
    Location at which the levmar library files were installed.
 ],
 [LEVMAR=yes; LEVMAR_LDFLAGS="-L$withval"]
)

if test "$LEVMAR" != no; then

LEVMAR_LIBS="-llevmar -llapack -lblas"

if test "$LEVMAR" != extern; then

AC_LANG_PUSH([C])
temp_status=true
CPPFLAGS="$CPPFLAGS $LEVMAR_CFLAGS"
AC_CHECK_HEADER(levmar.h,,[temp_status=false])
LDFLAGS="$LDFLAGS $LEVMAR_LDFLAGS"
AC_CHECK_LIB(levmar,dlevmar_dif,,[temp_status=false],[-llapack -lblas -lm])
AC_LANG_POP([C])

if test "$temp_status" = false; then
    if test "$LEVMAR" = yes; then
        AC_MSG_ERROR([levmar.h not found. Using --with-levmar, specify its location, "extern" to compile it locally, or "no" to disable it.
The package may be downloaded from http://www.ics.forth.gr/~lourakis/levmar/])
    else
	LEVMAR=extern
    fi
else
    LEVMAR=yes
fi

fi

if test "$LEVMAR" = extern; then

AC_MSG_NOTICE([I will compile levmar for you from the extern/ directory])

LEVMAR_MAKELIB=`printf 'liblevmar:
	mkdir -p $(EXTERN)/include $(EXTERN)/lib
	if [[ ! -f $(EXTERN)/include/levmar.h ]]; then \\
		(cd $(LEVMAR) && cmake .) && \\
		$(MAKE) -C $(LEVMAR) levmar C_FLAGS=-fPIC && \\
		cp $(LEVMAR)/liblevmar.a $(EXTERN)/lib/ && \\
		cp $(LEVMAR)/levmar.h $(EXTERN)/include/; \\
	fi\n'`

MAKE_LIBTARGETS="$MAKE_LIBTARGETS liblevmar"

LEVMAR_CFLAGS='-I$(EXTERN)/include'
LEVMAR_LDFLAGS='-L$(EXTERN)/lib'

fi

fi

CPPFLAGS="$temp_CPPFLAGS"
LDFLAGS="$temp_LDFLAGS"
LIBS="$temp_LIBS"

if test "$LEVMAR" != no; then
    AC_DEFINE(USE_LEVMAR)
fi

AC_SUBST(LEVMAR_CFLAGS)
AC_SUBST(LEVMAR_LDFLAGS)
AC_SUBST(LEVMAR_LIBS)
AC_SUBST(LEVMAR_MAKELIB)
AC_SUBST(MAKE_LIBTARGETS)
])
