ACLOCAL_AMFLAGS = -I m4

BINARCHDIR = bin/$(GAPARCH)
GAPINSTALLLIB = $(BINARCHDIR)/img.so

LEVMAR=extern/levmar-2.6
LIBDOGLEG=extern/libdogleg-0.15.4

lib_LTLIBRARIES = img.la
img_la_DEPENDENCIES = liblevmar
img_la_SOURCES = src/img.c src/p1.c src/rpoly.c
img_la_CPPFLAGS = $(GAP_CPPFLAGS) -I$(LEVMAR)
img_la_CFLAGS = $(GAP_CFLAGS)
img_la_LDFLAGS = $(GAP_LDFLAGS) -module -avoid-version -L$(BINARCHDIR)/lib -llevmar -llapack -lblas
if SYS_IS_CYGWIN
img_la_LDFLAGS += -no-undefined -version-info 0:0:0 -Wl,$(GAPROOT)/bin/$(GAPARCH)/gap.dll
endif

bin_PROGRAMS = layout hsolve
layout_DEPENDENCIES = libdogleg
layout_SOURCES = src/layout.C
layout_CPPFLAGS = -I$(LIBDOGLEG) $(CHOLMOD_CFLAGS)
layout_LDFLAGS = -L$(BINARCHDIR)/lib $(CHOLMOD_LDFLAGS)
layout_LDADD = -ldogleg -lcholmod -lsuitesparseconfig -llapack -lblas

hsolve_DEPENDENCIES = liblevmar
hsolve_SOURCES = src/hsolve.C
hsolve_CPPFLAGS = -I$(LEVMAR)
hsolve_LDFLAGS = -L$(BINARCHDIR)/lib
hsolve_LDADD = -llevmar -llapack -lblas

all-local: liblevmar libdogleg $(GAPINSTALLLIB) $(BINARCHDIR)/layout $(BINARCHDIR)/hsolve

liblevmar:
	mkdir -p $(BINARCHDIR)/include $(BINARCHDIR)/lib
	if [[ ! -f $(BINARCHDIR)/include/levmar.h ]]; then \
		(cd $(LEVMAR) && cmake .) && \
		$(MAKE) -B -C $(LEVMAR) levmar C_DEFINES="-O3 -fPIC" && \
		cp $(LEVMAR)/liblevmar.a $(BINARCHDIR)/lib/ && \
		cp $(LEVMAR)/levmar.h $(BINARCHDIR)/include/; \
	fi

libdogleg:
	mkdir -p "$(BINARCHDIR)/include" "$(BINARCHDIR)/lib"
	if [[ ! -f "$(BINARCHDIR)/include/dogleg.h" ]]; then \
	    $(MAKE) -B -C $(LIBDOGLEG) CFLAGS="-std=c99 -O3 $(CHOLMOD_CFLAGS)" LDFLAGS="$(CHOLMOD_LDFLAGS)" && \
	    cp $(LIBDOGLEG)/libdogleg.a $(BINARCHDIR)/lib/ && \
	    cp $(LIBDOGLEG)/dogleg.h $(BINARCHDIR)/include/; \
	fi

$(GAPINSTALLLIB): img.la
	$(mkdir_p) $(top_srcdir)/$(BINARCHDIR)
if SYS_IS_CYGWIN
	cp .libs/img.dll $(GAPINSTALLLIB)
else
	cp .libs/img.so $(GAPINSTALLLIB)
endif

$(BINARCHDIR)/layout: layout
	mv $< $@

$(BINARCHDIR)/hsolve: hsolve
	mv $< $@

clean-local:
	rm -f $(GAPINSTALLLIB)

distclean-local:
	rm -rf bin/*
	(cd doc ; ./clean)

checkblocks:
	grep '<#GAPDoc' PackageInfo.g gap/*d | awk -F'"' '{print $$2}' | sort > @@-blocks
	grep '<#Include' doc/img.xml | awk -F'"' '{print $$2}' | sort > @@-in
	comm -3 @@-blocks @@-in
	@rm @@-blocks @@-in

doc:
	(echo 'Read("makedoc.g");' | $(GAPROOT)/bin/gap.sh -A -q)

.PHONY: doc
