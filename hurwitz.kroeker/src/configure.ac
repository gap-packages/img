#############################################################################
##
#W  configure.ac                                            Laurent Bartholdi
##                                                              Jakob Kroeker
##
#H   @(#)$Id$
##
#Y Copyright (C) 2009, Laurent Bartholdi
##
#############################################################################


# jk: for AC_OPENMP check (parallelization)   2.62 is required! (2.62 available since  2009)
# AC_PREREQ(2.00)
AC_PREREQ([2.62])

AC_INIT(fr,,laurent.bartholdi@gmail.com)
AC_CONFIG_SRCDIR([src/fr_dll.c])
AC_CONFIG_AUX_DIR(cnf)

# Checks for programs.
AC_PROG_CC


################################################################
### lines added by jk

    AC_PROG_CXX

    # checking compile flag '-std=c++0x' does somehow not work.
    # AX_CHECK_COMPILE_FLAG([-std=c++0x], [CXXFLAGS="$CXXFLAGS -std=c++0x"], [ AX_MSG_ERROR( [need  -std=c++0x flag but it did not work]) ] )  
    #    CXXFLAGS="$CFLAGS  -std=c++0x"

   
    AC_OPENMP 

    # since some systems have broken OMP libraries
    # we also check that the actual package will work
    
    
    # @TODO: following check does not work; what is wrong?
    
#    ac_pkg_openmp=no
#    if test -n "${OPENMP_CFLAGS}"; then
#      AC_MSG_CHECKING([whether OpenMP will work in a package])
#      AC_LANG_CONFTEST(
#      [AC_LANG_PROGRAM([[#include <omp.h>]], [[ return omp_get_num_threads (); ]])])
#      PKG_CFLAGS="${OPENMP_CFLAGS}" PKG_LIBS="${OPENMP_CFLAGS}" "${CURDIR}" CMD SHLIB conftest.c 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && "${CURDIR}" --vanilla -q -e "dyn.load(paste('conftest',.Platform\$dynlib.ext,sep=''))" 1>&AS_MESSAGE_LOG_FD 2>&AS_MESSAGE_LOG_FD && ac_pkg_openmp=yes
#      AC_MSG_RESULT([${ac_pkg_openmp}])
#    fi
#
#    # if ${ac_pkg_openmp} = "yes" then we have OMP, otherwise it will be "no"
#    if test "${ac_pkg_openmp}" = no; then
#      OPENMP_CFLAGS=''
#      # you could put AC_MSG_ERROR here is OpenMP is required
#      AC_MSG_WARN([NO OpenMP for parallelization detected. Use more recent compiler!]) 
#    fi 
    
#    CFLAGS="$CFLAGS $OPENMP_CFLAGS "
    CXXFLAGS="$CXXFLAGS $OPENMP_CFLAGS"
    

# end changes by jk 
################################################################

                        

# Check for -fno-stack-protector, because we link within GAP
AC_CACHE_CHECK([whether $CC accepts -fno-stack-protector],
    [ns_cv_cc__nostackprotector],
    [save_CFLAGS=$CFLAGS
     CFLAGS="$CFLAGS -fno-stack-protector"
     AC_LINK_IFELSE([AC_LANG_PROGRAM([], [])],
                    [ns_cv_cc__nostackprotector=yes],
                    [ns_cv_cc__nostackprotector=no])
     CFLAGS=$save_CFLAGS])
if test $ns_cv_cc__nostackprotector = yes; then
     CFLAGS="$CFLAGS -fno-stack-protector"
fi

################################################################
# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

AC_ARG_VAR([GAPDIR], [Location of the GAP root directory, e.g. ../..])

if test -z "$GAPDIR"; then
    for dir in ../.. /Applications/gap4r5; do
	if test -f $dir/sysinfo.gap; then
	    GAPDIR=$dir
	    break
	fi
    done
fi

AC_ARG_VAR([CONFIGNAME],[Name of GAP build configuration])
if test -n "$CONFIGNAME"; then
    SYSINFO="$GAPDIR/sysinfo.gap-$CONFIGNAME"
    MAKEFILE="Makefile-$CONFIGNAME"
    GAPPROG="$GAPDIR/bin/gap-$CONFIGNAME.sh"
else
    SYSINFO="$GAPDIR/sysinfo.gap"
    MAKEFILE="Makefile"
    GAPPROG="$GAPDIR/bin/gap.sh"
fi

if ! test -f "$SYSINFO"; then
    AC_ERROR([Could not locate the GAP root directory;
        specify its location with './configure GAPDIR=DIR [CONFIGNAME=NAME)]'])
fi

GAPDIR=`cd $GAPDIR && pwd` # make path absolute
AC_SUBST(GAPDIR)

. "$SYSINFO"
TARGET="$GAParch"

echo checking target... "$TARGET"

XTARGET="`cnf/config.guess`-$CC-`echo $TARGET | sed 's/.*-//'`"
if test "$XTARGET" != "$GAParch_system"; then
   AC_WARN([The guessed target $XTARGET is not the gap target $GAParch_system. Cross your fingers])
fi
AC_SUBST(TARGET)

echo checking gap executable ... "$GAPPROG"

if ! test -e "$GAPPROG"; then
    AC_WARN([Could not find GAP executable $GAPPROG; I won't compile the documentation])
fi
AC_SUBST(GAPPROG)

GAC="$GAPDIR/bin/$TARGET/gac"

echo checking gac compiler... $GAC

if ! test -e "$GAC"; then
    AC_ERROR([Could not find GAP compiler $GAC])
fi
AC_SUBST(GAC)

eval `grep '^c_[[a-z_]]*=' $GAC`

AC_SUBST(c_compiler)
AC_SUBST(c_options)
AC_SUBST(c_linker)
AC_SUBST(c_link_options)
AC_SUBST(c_libs)
AC_SUBST(c_dyn_options)
AC_SUBST(c_dyn_linker)
AC_SUBST(c_dyn_linking)
AC_SUBST(c_dynlibs)
AC_SUBST(c_addlibs)


################################################################################################################################
### jk changes

    #    todo: dependencies between flintlib and mpir, mpft
    # prevent parallel make on mpfr, mpc, mpfi before mpfr is compiled
    # if test "$MPFRDIR" == "$EXTERN"; then MPFRDEPEND=mpfrlib; fi       
    
    # probably not required: 
    # AC_SUBST(cxx_compiler)
    # AC_SUBST(cxx_options)
    # AC_SUBST(cxx_linker)
    # AC_SUBST(cxx_link_options)
    # AC_SUBST(cxx_libs)
    # AC_SUBST(cxx_dyn_options)
    # AC_SUBST(cxx_dyn_linker)
    # AC_SUBST(cxx_dyn_linking)
    # AC_SUBST(cxx_dynlibs)
    # AC_SUBST(cxx_addlibs)


    ################################################################
    # mpir configuration

    EXTERN="\$(CURDIR)/bin/$TARGET/extern"

    MPIRDIR="$EXTERN"
    MPIRINCLUDE=""
    MPIRLIB=""

    AC_ARG_WITH(mpir,
     [  --with-mpir=<path>|yes|no|extern
        Location at which the MPIR library was installed.
        If the argument is omitted, the library is assumed to be reachable
        under the standard search path (/usr, /usr/local,...).  Otherwise
        you must give the <path> to the directory which contains the
        library. The special value "extern", which is the default, asks Float
        to compile a version of mpir in the subdirectory extern/.
     ],
     [if test "$withval" != extern; then MPIRDIR="$withval"; fi]
    )

    AC_ARG_WITH(mpir-include,
     [  --with-mpir-include=<location>
        Location at which the MPIR include files were installed.],
     [MPIRINCLUDE="$withval"]
    )

    AC_ARG_WITH(mpir-lib,
     [  --with-mpir-lib=<location>
        Location at which the MPIR library files were installed.],
     [MPIRLIB="$withval"]
    )

    if test "$MPIRDIR" != yes; then
    if test "$MPIRINCLUDE" == ""; then MPIRINCLUDE="$MPIRDIR/include"; fi
    if test "$MPIRLIB" == ""; then MPIRLIB="$MPIRDIR/lib"; fi
    fi

    if test "$MPIRINCLUDE" != ""; then
        CPPFLAGS="$CPPFLAGS -I$MPIRINCLUDE"
        CFLAGS="$CFLAGS -p -I$MPIRINCLUDE"
        CXXFLAGS="$CXXFLAGS -p -I$MPIRINCLUDE"
    fi

    ################################################################

    ################################################################
    # mpfr configuration

    EXTERN="\$(CURDIR)/bin/$TARGET/extern"

    MPFRDIR="$EXTERN"
    MPFRINCLUDE=""
    MPFRLIB=""

    AC_ARG_WITH(mpfr,
     [  --with-mpfr=<path>|yes|no|extern
        Location at which the MPFR library was installed.
        If the argument is omitted, the library is assumed to be reachable
        under the standard search path (/usr, /usr/local,...).  Otherwise
        you must give the <path> to the directory which contains the
        library. The special value "extern", which is the default, asks Float
        to compile a version of mpfr in the subdirectory extern/.
     ],
     [if test "$withval" != extern; then MPFRDIR="$withval"; fi]
    )

    AC_ARG_WITH(mpfr-include,
     [  --with-mpfr-include=<location>
        Location at which the MPFR include files were installed.],
     [MPFRINCLUDE="$withval"]
    )

    AC_ARG_WITH(mpfr-lib,
     [  --with-mpfr-lib=<location>
        Location at which the MPFR library files were installed.],
     [MPFRLIB="$withval"]
    )

    if test "$MPFRDIR" != yes; then
    if test "$MPFRINCLUDE" == ""; then MPFRINCLUDE="$MPFRDIR/include"; fi
    if test "$MPFRLIB" == ""; then MPFRLIB="$MPFRDIR/lib"; fi
    fi

    if test "$MPFRINCLUDE" != ""; then
        CPPFLAGS="$CPPFLAGS -I$MPFRINCLUDE"
        CFLAGS="$CFLAGS -p -I$MPFRINCLUDE"
        CXXFLAGS="$CXXFLAGS -p -I$MPFRINCLUDE"
    fi

    ################################################################

### end changes by jk
################################################################################################################################


# gmp configuration

GMPDIR="$GAPDIR/bin/$TARGET/extern/gmp"
GMPINCLUDE=""
GMPLIB=""

AC_ARG_WITH(gmp,
 [  --with-gmp=<location>|yes|no|gap
    Location at which the GMP library, needed for FR, was installed.
    If the argument is omitted, the library is assumed to be reachable
    under the standard search path (/usr, /usr/local,...).  Otherwise
    you must give the <path> to the directory which contains the
    library. The special value "gap", which is the default, asks FR
    to use the version of gmp included in the GAP distribution.
 ],
 [if test "$withval" != gap; then GMPDIR="$withval"; fi]
)

AC_ARG_WITH(gmp-include,
 [  --with-gmp-include=<location>
    Location at which the GMP include files were installed.],
 [GMPINCLUDE="$withval"]
)

AC_ARG_WITH(gmp-lib,
 [  --with-gmp-lib=<location>
    Location at which the GMP library files were installed.],
 [GMPLIB="$withval"]
)

if test "$GMPDIR" != yes; then
if test "$GMPINCLUDE" == ""; then GMPINCLUDE="$GMPDIR/include"; fi
if test "$GMPLIB" == ""; then GMPLIB="$GMPDIR/lib"; fi
fi

echo using GMP directory... $GMPDIR

if test "$GMPINCLUDE" != ""; then
    CPPFLAGS="$CPPFLAGS -I$GMPINCLUDE"
    CFLAGS="$CFLAGS -p -I$GMPINCLUDE"
fi
AC_CHECK_HEADER(gmp.h,[],[AC_MSG_ERROR([library gmp not found. Specify its location using --with-gmp])],[])

################################################################
# Check for gsl library

EXTERN="\$(CURDIR)/bin/$TARGET/extern"
LIB_TARGET=""

GSLDIR="$EXTERN"
GSLINCLUDE=""
GSLLIB=""
AC_ARG_WITH(gsl,
 [  --with-gsl=<path>|yes|no|extern
    Location at which the GSL library was installed.
    If the argument is omitted, the library is assumed to be reachable
    under the standard search path (/usr, /usr/local,...).  Otherwise
    you must give the <path> to the directory which contains the
    library. The special value "extern", which is the default, asks Float
    to compile a version of gsl in the subdirectory extern/.
 ],
 [if test "$withval" != extern; then GSLDIR="$withval"; fi]
)

AC_ARG_WITH(gsl-include,
 [  --with-gsl-include=<location>
    Location at which the GSL include files were installed.],
 [GSLINCLUDE="$withval"]
)

AC_ARG_WITH(gsl-lib,
 [  --with-gsl-lib=<location>
    Location at which the GSL library files were installed.],
 [GSLLIB="$withval"]
)

if test "$GSLDIR" != yes; then
if test "$GSLINCLUDE" == ""; then GSLINCLUDE="$GSLDIR/include"; fi
if test "$GSLLIB" == ""; then GSLLIB="$GSLDIR/lib"; fi
fi

echo using GSL directory... $GSLDIR

if test "$GSLINCLUDE" != ""; then
    CPPFLAGS="$CPPFLAGS -I$GSLINCLUDE"
    CFLAGS="$CFLAGS -I$GSLINCLUDE"
fi

if test "$GSLDIR" == "$EXTERN"; then
    LIB_TARGET="$LIB_TARGET gsllib"
else
    AC_CHECK_HEADER(gsl/gsl_vector.h,[],[AC_MSG_ERROR([library gsl not found. Specify its location using --with-gsl])],[])
fi

if test "$GSLLIB" != ""; then
    LIBS="$LIBS -L$GSLLIB"
    GACFLAGS="$GACFLAGS -L -L$GSLLIB -L -Wl,-rpath,$GSLLIB -L -lgsl -L -lgslcblas"
    if test "$GSLLIB" != "$EXTERN/lib"; then
        AC_CHECK_LIB([gsl],[gsl_multiroot_fsolver_set],,
            [AC_ERROR([The GSL library could not be found. It is needed for IMG calculations.])],[-lgslcblas])

        AC_CHECK_LIB([gslcblas],[cblas_ctrmv],,
            [AC_ERROR([The GSL CBlas library could not be found. It is needed for IMG calculations.])])
    fi
fi

################################################################
# Check for givaro library

GIVARODIR="$EXTERN"
GIVAROINCLUDE=""
GIVAROLIB=""
AC_ARG_WITH(givaro,
 [  --with-givaro=<path>|yes|no|extern
    Location at which the Givaro library was installed.
    If the argument is omitted, the library is assumed to be reachable
    under the standard search path (/usr, /usr/local,...).  Otherwise
    you must give the <path> to the directory which contains the
    library. The special value "extern", which is the default, asks Float
    to compile a version of givaro in the subdirectory extern/.
 ],
 [if test "$withval" != extern; then GIVARODIR="$withval"; fi]
)

AC_ARG_WITH(givaro-include,
 [  --with-givaro-include=<location>
    Location at which the Givaro include files were installed.],
 [GIVAROINCLUDE="$withval"]
)

AC_ARG_WITH(givaro-lib,
 [  --with-givaro-lib=<location>
    Location at which the Givaro library files were installed.],
 [GIVAROLIB="$withval"]
)

if test "$GIVARODIR" != yes; then
if test "$GIVAROINCLUDE" == ""; then GIVAROINCLUDE="$GIVARODIR/include"; fi
if test "$GIVAROLIB" == ""; then GIVAROLIB="$GIVARODIR/lib"; fi
fi

echo using Givaro directory... $GIVARODIR

if test "$GIVAROINCLUDE" != ""; then
    CPPFLAGS="$CPPFLAGS -I$GIVAROINCLUDE"
    CFLAGS="$CFLAGS -I$GIVAROINCLUDE"
fi

if test "$GIVARODIR" == "$EXTERN"; then
    LIB_TARGET="$LIB_TARGET givarolib"
else
    AC_CHECK_HEADER(givaro-config.h,[],[AC_MSG_ERROR([library givaro not found. Specify its location using --with-givaro])],[])
fi

if test "$GIVAROLIB" != ""; then
    LIBS="$LIBS -L$GIVAROLIB"
    #GACFLAGS="$GACFLAGS -L -L$GIVAROLIB -L -Wl,-rpath,$GIVAROLIB -L -lgivaro"
fi

################################################################
# external programs configuration

AC_PATH_PROGS(DOT,[dot],[],[$PATH$PATH_SEPARATOR/usr/local/graphviz/bin])

if test -z "$DOT"; then
    AC_WARN([Could not find 'dot' (debian package graphviz)... you won't be able to draw automata])
fi

AC_PATH_PROG(DISP,[display])

if test -z "$DISP"; then
    AC_WARN([Could not find 'display' (debian package imagemagick)... you won't be able to draw automata])
fi

AC_PATH_PROG(APPLETVIEWER,[appletviewer])

if test -z "$APPLETVIEWER"; then
    AC_WARN([Could not find 'appletviewer' (debian package java-6-sdk)... you won't be able to draw spiders])
fi

# Check for java compiler
AC_CHECK_PROGS(JAVAC,[javac])

if test -z "$JAVAC"; then
    AC_WARN([Could not find java compiler... you won't be able to draw spiders])
else
    JAVABUILD="java/javaplot.class java/javaview.jar"
fi
AC_SUBST(JAVABUILD)
AC_SUBST(GACFLAGS)
AC_SUBST(LIBS)
AC_SUBST(LIB_TARGET)

################################################################
# generate files

WITHGMP=""

INCLGMP=""
LINKGMP=""

if test "$GMPDIR" != ""; then
   WITHGMP="$WITHGMP --with-gmp=$GMPDIR"
fi


if test "$GMPINCLUDE" != "";  then
    INCLGMP="-I$GMPINCLUDE"
    if test "$GMPLIB" != ""; then
        LINKGMP="-L$GMPLIB"
        # changes by jk: WITHGMP should define either '--with-gmp' or '--with-gmp-include' AND '--with-gmp-lib' !
        WITHGMP="$ --with-gmp-include=$GMPINCLUDE  --with-gmp-lib=$GMPLIB"
    fi
fi

if test "$GMPLIB" != ""; then
    LINKGMP="-L$GMPLIB"
    if test "$GMPINCLUDE" != "";  then
        INCLGMP="-I$GMPINCLUDE"
        # changes by jk: WITHGMP should define either '--with-gmp' or '--with-gmp-include' AND '--with-gmp-lib' !
        WITHGMP="$ --with-gmp-include=$GMPINCLUDE  --with-gmp-lib=$GMPLIB"
    fi
fi

if test "$WITHGMP" = "";  then
    AC_ERROR([ Please set GMPDIR or GMPLIB and GMPINCLUDE! ])
fi;



AC_SUBST(WITHGMP)
AC_SUBST(INCLGMP)
AC_SUBST(LINKGMP)


################################################################
## changes by jk
    WITHMPIR=""
    WITHMPFR=""
    
  
    if test "$MPIRDIR" != ""; then
       WITHMPIR="$WITHMPIR --with-mpir=$MPIRDIR"
    fi
    
    if test "$MPIRLIB" != ""; then
       WITHMPIRLIB="$  --with-mpir-lib=$MPIRLIB"
       LINKMPIR=" -L$MPIRLIB "
    fi

    if test "$MPIRINCLUDE" != ""; then
       WITHMPIRINCLUDE="$  --with-mpir-include=$MPIRINCLUDE"
       INCLMPIR=" -I$MPIRINCLUDE "
    fi

    if test "$WITHMPIR" = "";  then
        AC_ERROR([ MPIR not found! ])
    fi;

  
    if test "$MPFRDIR" != ""; then
       WITHMPFR="$WITHMPFR --with-mpfr=$MPFRDIR"
    fi

    if test "$MPFRLIB" != ""; then
       WITHMPFRLIB="$  --with-mpfr-lib=$MPFRLIB"
       LINKMPFR=" -L$MPFRLIB "
    fi

    if test "$MPFRINCLUDE" != ""; then
       WITHMPFRINCLUDE="$  --with-mpfr-include=$MPFRINCLUDE"
       INCLMPFR=" -I$MPFRINCLUDE "
    fi
    
    if test "$WITHMPFR" = "";  then
        AC_ERROR([ MPFR not found! ])
    fi;

    AC_SUBST(WITHMPFR)
    AC_SUBST(WITHMPIR)
# end changes by jk
################################################################


mkdir -p bin/$TARGET
CONFIG_STATUS=bin/$TARGET/config.status

AC_CONFIG_FILES([$MAKEFILE:cnf/Makefile.in])

if test "$MAKEFILE" != Makefile; then
    ln -sf "$MAKEFILE" Makefile
fi

if test "$GSLDIR" == "$EXTERN" -o "$GIVARODIR" == "$EXTERN"; then
    echo -n "****** Remember to download:"
    if test "$GSLDIR" == "$EXTERN"; then echo -n " gsl-1.15.tar.gz"; fi
    if test "$GIVARODIR" == "$EXTERN"; then echo -n " givaro-3.6.0.tar.gz"; fi
    echo; echo "****** (instructions in extern/GET_LIBRARIES)"
fi

AC_OUTPUT
