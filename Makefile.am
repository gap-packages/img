ACLOCAL_AMFLAGS = -I m4

BINARCHDIR = bin/$(GAPARCH)
GAPINSTALLLIB = $(BINARCHDIR)/img.so

LEVMAR=extern/levmar-2.6
LIBDOGLEG=extern/libdogleg-0.09

lib_LTLIBRARIES = img.la
img_la_SOURCES = src/img.c src/p1.c src/rpoly.c
img_la_CPPFLAGS = $(GAP_CPPFLAGS) -DCONFIG_H
# Note that the latter is only for GAP 4.4.12
img_la_LDFLAGS = -module -avoid-version -L$(BINARCHDIR)/lib -llevmar -llapack -lblas
if SYS_IS_CYGWIN
img_la_LDFLAGS += -no-undefined -version-info 0:0:0 -Wl,$(GAPROOT)/bin/$(GAPARCH)/gap.dll
endif

bin_PROGRAMS = layout hsolve
layout_DEPENDENCIES = libdogleg
layout_SOURCES = src/layout.C
layout_CPPFLAGS = -I$(LIBDOGLEG)
layout_LDFLAGS = -L$(BINARCHDIR)/lib
layout_LDADD = -ldogleg -lcholmod

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
		$(MAKE) -B -C $(LEVMAR) levmar C_FLAGS=-fPIC && \
		cp $(LEVMAR)/liblevmar.a $(BINARCHDIR)/lib/ && \
		cp $(LEVMAR)/levmar.h $(BINARCHDIR)/include/; \
	fi

libdogleg:
	mkdir -p "$(BINARCHDIR)/include" "$(BINARCHDIR)/lib"
	if [[ ! -f "$(BINARCHDIR)/include/dogleg.h" ]]; then \
	    $(MAKE) -B -C $(LIBDOGLEG) CFLAGS=$(CHOLMOD_CFLAGS) && \
	    cp $(LIBDOGLEG)/libdogleg.a $(BINARCHDIR)/lib/ && \
	    cp $(LIBDOGLEG)/dogleg.h $(BINARCHDIR)/include/ && \
	    if [[ "$(shell uname)" = Darwin ]]; then \
		cp $(LIBDOGLEG)/libdogleg.dylib $(BINARCHDIR)/lib/libdogleg.dylib; \
		install_name_tool -id $(BINARCHDIR)/lib/libdogleg.dylib $(BINARCHDIR)/lib/libdogleg.dylib; \
	    fi; \
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
