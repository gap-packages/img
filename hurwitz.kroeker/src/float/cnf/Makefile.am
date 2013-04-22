#############################################################################
##
#W Makefile.am                                              Laurent Bartholdi
##
#H   @(#)$Id$
##
#Y Copyright (C) 2007, Laurent Bartholdi
##
#############################################################################
##
##  This compiles the module polroots, and creates archives
##
#############################################################################
ACLOCAL_AMFLAGS = -I m4


GAC=@GAC@
GAP=@GAPPROG@
LOCALBIN=bin/@TARGET@
EXTERN=$(CURDIR)/bin/@TARGET@/extern
GACFLAGS=@GACFLAGS@ -p -g @GAC_stack_protector_flag@
GACPP=bin/@TARGET@/gac++
MPFRLIB=mpfr-3.1.0
MPFILIB=mpfi-1.5
MPCLIB=mpc-0.9
FPLLLLIB=libfplll-3.1.3
CXSCLIB=cxsc-2-5-1

.PHONY: all distribute doc clean mrproper wwwdir checkblocks \
	lib mpfrlib mpfilib mpclib fpllllib cxsclib  

all: @LIB_TARGET@ @DLL_TARGET@

lib: mpfrlib mpfilib mpclib fpllllib cxsclib


LINK = $(LIBTOOL) --tag=CC $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) \
	--mode=link $(CCLD) $(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) \
	$(LDFLAGS) -o $@
	

extern/$(MPFRLIB).tar.bz2: 
	echo "I can't find $(MPFRLIB), so I'm going to download it"
	(cd extern; @GETBIN@ @GETBINPARAM@  http://www.mpfr.org/$(MPFRLIB)/$(MPFRLIB).tar.bz2 )

extern/$(MPFILIB).tar.bz2: 
	echo "I can't find $(MPFILIB), so I'm going to download it"
	(cd extern; @GETBIN@ @GETBINPARAM@  https://gforge.inria.fr/frs/download.php/27345/mpfi-1.5.tar.bz2   )

extern/$(MPCLIB).tar.gz: 
	echo "I can't find $(MPCLIB), so I'm going to download it"
	(cd extern; @GETBIN@ @GETBINPARAM@  http://www.multiprecision.org/mpc/download/$(MPCLIB).tar.gz  )

extern/$(CXSCLIB).tar.gz: 
	echo "I can't find $(CXSCLIB), so I'm going to download it"
	(cd extern; @GETBIN@ @GETBINPARAM@  http://www2.math.uni-wuppertal.de/~xsc/xsc/cxsc/$(CXSCLIB).tar.gz  )

extern/$(FPLLLLIB).tar.gz: 
	echo "I can't find $(FPLLLLIB), so I'm going to download it"
	(cd extern; @GETBIN@ @GETBINPARAM@  http://perso.ens-lyon.fr/xavier.pujol/fplll/$(FPLLLLIB).tar.gz  )


mpfrlib: extern/$(MPFRLIB).tar.bz2
	if ! test -r $(EXTERN)/include/mpfr.h; then \
	    cd extern && \
	    tar -x -f $(MPFRLIB).tar.bz2 -j && \
	    cd $(MPFRLIB) && \
	    ./configure @CONFIGUREPARAMS@ @WITHGMP@ --prefix=$(EXTERN) && \
	    $(MAKE) install; \
	fi

mpfilib: @MPFRDEPEND@ extern/$(MPFILIB).tar.bz2
	if ! test -r $(EXTERN)/include/mpfi.h; then \
	    cd extern && \
	    tar -x -f $(MPFILIB).tar.bz2 -j && \
	    cd $(MPFILIB) && \
	    ./configure  @CONFIGUREPARAMS@ @WITHGMP@ @WITHMPFR@ --prefix=$(EXTERN) && \
	    $(MAKE) install; \
	fi

mpclib: @MPFRDEPEND@ extern/$(MPCLIB).tar.gz
	if ! test -r $(EXTERN)/include/mpc.h; then \
	    cd extern && \
	    tar -x -f $(MPCLIB).tar.gz -z && \
	    cd $(MPCLIB) && \
	    ./configure  @CONFIGUREPARAMS@  @WITHGMP@ @WITHMPFR@ --prefix=$(EXTERN) && \
	    $(MAKE) install; \
	fi

fpllllib: @MPFRDEPEND@ extern/$(FPLLLLIB).tar.gz
	if ! test -r $(EXTERN)/include/fplll.h; then \
	    cd extern && \
	    tar -x -f $(FPLLLLIB).tar.gz -z && \
	    cd $(FPLLLLIB) && \
	    ./configure  @CONFIGUREPARAMS@ CPPFLAGS="@INCLGMP@ @INCLMPFR@ $(CPPFLAGS)" LDFLAGS="@LINKGMP@ @LINKMPFR@ $(LDFLAGS)" --prefix=$(EXTERN) && \
	    $(MAKE) install; \
	fi

cxsclib: extern/$(CXSCLIB).tar.gz
	mkdir -p $(EXTERN)
	rm -f $(EXTERN)/cxsc
	ln -s . $(EXTERN)/cxsc
# a dirty hack: we'll pretend our home directory
# is $(EXTERN), and the installation will go into $(EXTERN)/cxsc
	if ! test -r $(EXTERN)/include/real.hpp; then \
	    cd extern && \
	    tar -x -f $(CXSCLIB).tar.gz -z && \
	    cd $(CXSCLIB) && \
	    (echo "yes"; for i in 1 2 3 4 5 6 7 8 9 10; do echo ""; done) | \
		HOME=$(EXTERN) ./install_cxsc; \
	fi

distribute: doc wwwdir doc tarballs
	rsync -arvp --delete www/ laurent@elanpc00.uni-math.gwdg.de:public_html/Float/

$(LOCALBIN):
	mkdir -p $(LOCALBIN)

$(GACPP):
	sed -e 's/"gcc"/"g++"/g' -e 's/\([* ]\.\)c/\1C/g' < $(GAC) > $(GACPP)
	chmod a+x $(GACPP)

$(LOCALBIN)/%.o: src/%.C $(GACPP)
	$(GACPP) $(GACFLAGS) -d -o $@ -c $<

$(LOCALBIN)/%.o: src/%.c
	$(GAC) $(GACFLAGS) -d -o $@ -c $<

$(LOCALBIN)/mp_float.so: @MP_FLOAT_O@ $(GACPP)
	$(GACPP) $(GACFLAGS) -d -o $@ @MP_FLOAT_O@ @MP_FLOAT_LIB@

$(LOCALBIN)/cxsc_float.so: @CXSC_FLOAT_O@ $(GACPP)
	$(GACPP) $(GACFLAGS) -d -o $@ @CXSC_FLOAT_O@ @CXSC_FLOAT_LIB@

$(LOCALBIN)/mp_float.o: src/mp_float.c src/mp_float.h

$(LOCALBIN)/mpfr.o: src/mpfr.c src/mp_float.h

$(LOCALBIN)/mpfi.o: src/mpfi.c src/mp_float.h

$(LOCALBIN)/mpc.o: src/mpc.c src/mp_float.h

$(LOCALBIN)/fplll.o: src/fplll.C

$(LOCALBIN)/mp_poly.o: src/cpoly.C src/mp_poly.C

$(LOCALBIN)/cxsc_poly.o: src/cpoly.C src/cxsc_poly.C src/cxsc_float.h

clean:
	rm -rf .version $(LOCALBIN) `find doc -type l`

configure: cnf/Makefile.in cnf/configure.ac
	(cd cnf; autoconf; mv -f configure ..)

mrproper: clean
	rm Makefile

.version: PackageInfo.g
	grep '^Version :=' $< | awk -F'"' '{print $$2}' > $@

wwwdir: .version tarballs
	mkdir -p www
	rm -f `find www -type l`
	cp README www/README.float
	cp PackageInfo.g www/PackageInfo.g
	cp doc/chap0.html www/index.html
	ln -s float-`cat .version`.tar.gz www/float.tar.gz
	cp doc/manual.pdf www/manual.pdf
	(cd doc; for i in *.html; do cp $$i ../www/$$i; done)
	cp doc/manual.css www/manual.css

doc: doc/chap0.html

doc/chap0.html: doc/float.xml lib/float.gd
	echo 'LoadPackage("float"); MAKEDOC@Float();' | $(GAP)

checkblocks:
	grep '<#GAPDoc' lib/*.g[di] | awk -F'"' '{print $$2}' | sort > @@-blocks
	grep '<#Include' doc/float.xml | awk -F'"' '{print $$2}' | sort > @@-in
	comm -3 @@-blocks @@-in
	rm @@-blocks @@-in

tarballs: .version doc
	mkdir -p www
	tar cfz www/float-`cat .version`.tar.gz --exclude '*~' --exclude sandbox --exclude bin --exclude www --exclude CVS --exclude .version --exclude cnf/autom4te.cache --exclude 'Makefile*' --exclude config.log --exclude 'extern/*-[0-9]*' -C .. float

#E Makefile . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
