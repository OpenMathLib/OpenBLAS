TOPDIR	= ..

include ../Makefile.system

ifndef EXPRECISION
EXPRECISION	= 0
endif

ifndef NO_CBLAS
NO_CBLAS	= 0
endif

ifndef NO_LAPACK
NO_LAPACK	= 0
endif

ifndef NO_LAPACKE
NO_LAPACKE	= 0
endif

ifeq ($(OSNAME), WINNT)
ifeq ($(F_COMPILER), GFORTRAN)
EXTRALIB += -lgfortran
endif
endif

ifeq ($(OSNAME), CYGWIN_NT)
ifeq ($(F_COMPILER), GFORTRAN)
EXTRALIB += -lgfortran
endif
endif

all::

libs::

prof::

hpl:: libgoto_hpl.$(LIBSUFFIX)

hpl_p:: libgoto_hpl_p.$(LIBSUFFIX)

libgoto_hpl.$(LIBSUFFIX) : ../$(LIBNAME)
	rm -f $(@F)
	$(LD) -r $(LDFLAGS) -o goto.$(SUFFIX) --whole-archive $< --no-whole-archive
	$(AR) cq $(@F) goto.$(SUFFIX)
	$(RANLIB) libgoto_hpl.$(LIBSUFFIX)

libgoto_hpl_p.$(LIBSUFFIX) :  ../$(LIBNAME_P)
	rm -f $(@F)
	$(LD) -r  $(LDFLAGS) -o goto.$(PSUFFIX)  --whole-archive $< --no-whole-archive
	$(AR) cq $(@F) goto.$(PSUFFIX)
	$(RANLIB) libgoto_hpl_p.$(LIBSUFFIX)

libgoto_hpl.dll : libgoto_hpl.$(LIBSUFFIX) dllinit.$(SUFFIX) libgoto_hpl.def
	$(DLLWRAP) -o $(@F) --def libgoto_hpl.def --entry _dllinit -s dllinit.$(SUFFIX) --dllname libgoto_hpl.dll libgoto_hpl.$(LIBSUFFIX)
	lib /machine:X64 /def:libgoto_hpl.def

dyn : $(LIBDYNNAME)

zip : dll
	zip $(LIBZIPNAME) $(LIBDLLNAME) $(LIBNAME)

dll  : ../$(LIBDLLNAME)
#libgoto2.dll

dll2 : libgoto2_shared.dll

../$(LIBDLLNAME) : ../$(LIBNAME) libopenblas.def dllinit.$(SUFFIX)
	$(RANLIB) ../$(LIBNAME)
ifeq ($(BINARY32), 1)
	$(DLLWRAP) -o ../$(LIBDLLNAME) --def libopenblas.def \
	--entry _dllinit@12 -s dllinit.$(SUFFIX) --dllname $(@F) ../$(LIBNAME) $(FEXTRALIB)
	-lib /machine:i386 /def:libopenblas.def
else
	$(DLLWRAP) -o ../$(LIBDLLNAME) --def libopenblas.def \
	--entry $(FU)dllinit    -s dllinit.$(SUFFIX) --dllname $(@F) ../$(LIBNAME) $(FEXTRALIB)
	-lib /machine:X64 /def:libopenblas.def 
endif

libgoto2_shared.dll : ../$(LIBNAME) libgoto2_shared.def
	$(CC) $(CFLAGS) libgoto2_shared.def -shared -o $(@F) \
	-Wl,--whole-archive ../$(LIBNAME) -Wl,--no-whole-archive \
	-Wl,--out-implib,libgoto2_shared.lib $(FEXTRALIB)

libopenblas.def : gensymbol
	perl ./gensymbol win2k    $(ARCH) dummy $(EXPRECISION) $(NO_CBLAS) $(NO_LAPACK) $(NO_LAPACKE) > $(@F)

libgoto2_shared.def : gensymbol
	perl ./gensymbol win2k    $(ARCH) dummy $(EXPRECISION) $(NO_CBLAS) $(NO_LAPACK) $(NO_LAPACKE) > $(@F)

libgoto_hpl.def : gensymbol
	perl ./gensymbol win2khpl $(ARCH) dummy $(EXPRECISION) $(NO_CBLAS) $(NO_LAPACK) $(NO_LAPACKE) > $(@F)

$(LIBDYNNAME) : ../$(LIBNAME) osx.def
	$(CC) $(CFLAGS) -all_load -headerpad_max_install_names -install_name $(CURDIR)/../$(LIBDYNNAME) -dynamiclib -o ../$(LIBDYNNAME) $< -Wl,-exported_symbols_list,osx.def  $(FEXTRALIB)

symbol.$(SUFFIX) : symbol.S
	$(CC) $(CFLAGS) -c -o $(@F) $^

dllinit.$(SUFFIX) : dllinit.c
	$(CC) $(CFLAGS) -c -o $(@F) -s $<

ifeq ($(OSNAME), Linux)

so : ../$(LIBSONAME)

../$(LIBSONAME) : ../$(LIBNAME) linux.def linktest.c
	$(CC) $(CFLAGS) -shared -o ../$(LIBSONAME) \
	-Wl,--whole-archive ../$(LIBNAME) -Wl,--no-whole-archive \
	-Wl,--retain-symbols-file=linux.def -Wl,-soname,$(LIBPREFIX).so.$(MAJOR_VERSION) $(EXTRALIB)
	$(CC) $(CFLAGS) -w -o linktest linktest.c ../$(LIBSONAME) $(FEXTRALIB) && echo OK.
	rm -f linktest

endif

ifeq ($(OSNAME), FreeBSD)

so : ../$(LIBSONAME)

../$(LIBSONAME) : ../$(LIBNAME) linux.def linktest.c
	$(CC) $(CFLAGS)  -shared -o ../$(LIBSONAME) \
	-Wl,--whole-archive ../$(LIBNAME) -Wl,--no-whole-archive \
	-Wl,--retain-symbols-file=linux.def $(EXTRALIB)
	$(CC) $(CFLAGS) -w -o linktest linktest.c ../$(LIBSONAME) $(FEXTRALIB) && echo OK.
	rm -f linktest

endif

ifeq ($(OSNAME), OSF1)

so : ../$(LIBSONAME)

../$(LIBSONAME) :
	$(CC) -shared -o ../$(LIBSONAME) ../$(LIBNAME)
endif

ifeq ($(OSNAME), SunOS)

so : ../$(LIBSONAME)
	$(CC) $(CFLAGS)  -shared -o ../$(LIBSONAME) \
	-Wl,--whole-archive ../$(LIBNAME) -Wl,--no-whole-archive $(EXTRALIB)
	$(CC) $(CFLAGS) -w -o linktest linktest.c ../$(LIBSONAME) $(FEXTRALIB) && echo OK.
	rm -f linktest

endif

ifeq ($(OSNAME), AIX)

ifeq ($(COMPILER_F77), xlf)

goto32.$(SUFFIX) : ../$(LIBNAME) aix.def
	ld -o $(@F) ../$(LIBNAME) -bE:aix.def -bM:SRE -bnoexpall -bnoentry -L$(HOME)/misc/lib -lxlf90 -lc -lm  -lpthread

goto64.$(SUFFIX) : ../$(LIBNAME) aix.def
	ld -b64 -o $(@F) ../$(LIBNAME) -bE:aix.def -bM:SRE -bnoexpall -bnoentry -L$(HOME)/misc/lib/ppc64 -lxlf90 -lc -lm -lpthread
else
goto32.$(SUFFIX) : ../$(LIBNAME) aix.def
	ld -o $(@F) ../$(LIBNAME) -bE:aix.def -bM:SRE -bnoexpall -bnoentry -L$(HOME)/misc/lib -lg2c -lc -lm 

goto64.$(SUFFIX) : ../$(LIBNAME) aix.def
	ld -b64 -o $(@F) ../$(LIBNAME) -bE:aix.def -bM:SRE -bnoexpall -bnoentry -L$(HOME)/misc/lib/ppc64 -lg2c -lc -lm 
endif
endif

static : ../$(LIBNAME)
	$(LD) $(LDFLAGS) -r -o goto.$(SUFFIX) \
	--whole-archive ../$(LIBNAME) --no-whole-archive
	rm -f ../$(LIBNAME)
	$(AR) -cq ../$(LIBNAME) goto.$(SUFFIX)
	rm -f goto.$(SUFFIX)

linux.def : gensymbol ../Makefile.system ../getarch.c
	perl ./gensymbol linux $(ARCH) $(BU) $(EXPRECISION) $(NO_CBLAS)  $(NO_LAPACK) $(NO_LAPACKE) > $(@F)

osx.def : gensymbol ../Makefile.system ../getarch.c
	perl ./gensymbol osx $(ARCH) $(BU) $(EXPRECISION) $(NO_CBLAS)  $(NO_LAPACK) $(NO_LAPACKE) > $(@F)

aix.def : gensymbol ../Makefile.system ../getarch.c
	perl ./gensymbol aix $(ARCH) $(BU) $(EXPRECISION) $(NO_CBLAS)  $(NO_LAPACK) $(NO_LAPACKE) > $(@F)

symbol.S : gensymbol
	perl ./gensymbol win2kasm noarch dummy $(EXPRECISION) $(NO_CBLAS)  $(NO_LAPACK) $(NO_LAPACKE) > symbol.S

test : linktest.c
	$(CC) $(CFLAGS) -w -o linktest linktest.c ../$(LIBSONAME) -lm && echo OK.
	rm -f linktest

linktest.c : gensymbol ../Makefile.system ../getarch.c
	perl ./gensymbol linktest  $(ARCH) $(BU) $(EXPRECISION) $(NO_CBLAS) $(NO_LAPACK) $(NO_LAPACKE) > linktest.c

clean ::
	@rm -f *.def *.dylib __.SYMDEF*

include ../Makefile.tail


