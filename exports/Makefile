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

ifndef NEED2UNDERSCORES
NEED2UNDERSCORES=0
endif

ifndef ONLY_CBLAS
ONLY_CBLAS	= 0
endif

ifeq ($(OSNAME), WINNT)
ifeq ($(F_COMPILER), GFORTRAN)
EXTRALIB += -lgfortran
endif
ifeq ($(USE_OPENMP), 1)
ifeq ($(C_COMPILER), GCC)
EXTRALIB += -lgomp
endif
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

# On Windows, we only generate a DLL without a version suffix.  This is because
# applications which link against the dynamic library reference a fixed DLL name
# in their import table.  By instead using a stable name it is possible to
# upgrade between library versions, without needing to re-link an application.
# For more details see: https://github.com/xianyi/OpenBLAS/issues/127.
../$(LIBDLLNAME) : ../$(LIBNAME) libopenblas.def dllinit.$(SUFFIX)
	$(RANLIB) ../$(LIBNAME)
	$(CC) $(CFLAGS) $(LDFLAGS) libopenblas.def dllinit.$(SUFFIX) \
	-shared -o ../$(LIBDLLNAME) -Wl,--out-implib,../$(LIBDLLNAME).a \
	-Wl,--whole-archive ../$(LIBNAME) -Wl,--no-whole-archive $(FEXTRALIB) $(EXTRALIB)

libopenblas.def : gensymbol
	perl ./gensymbol win2k    $(ARCH) dummy $(EXPRECISION) $(NO_CBLAS) $(NO_LAPACK) $(NO_LAPACKE) $(NEED2UNDERSCORES) $(ONLY_CBLAS) > $(@F)

libgoto_hpl.def : gensymbol
	perl ./gensymbol win2khpl $(ARCH) dummy $(EXPRECISION) $(NO_CBLAS) $(NO_LAPACK) $(NO_LAPACKE) $(NEED2UNDERSCORES) $(ONLY_CBLAS) > $(@F)

$(LIBDYNNAME) : ../$(LIBNAME) osx.def
	$(FC) $(FFLAGS) -all_load -headerpad_max_install_names -install_name $(CURDIR)/../$(LIBDYNNAME) -dynamiclib -o ../$(LIBDYNNAME) $< -Wl,-exported_symbols_list,osx.def  $(FEXTRALIB)

dllinit.$(SUFFIX) : dllinit.c
	$(CC) $(CFLAGS) -c -o $(@F) -s $<

ifeq ($(OSNAME), Linux)

so : ../$(LIBSONAME)

../$(LIBSONAME) : ../$(LIBNAME) linktest.c
ifneq ($(C_COMPILER), LSB)
	$(CC) $(CFLAGS) $(LDFLAGS) -shared -o ../$(LIBSONAME) \
	-Wl,--whole-archive ../$(LIBNAME) -Wl,--no-whole-archive \
	-Wl,-soname,$(LIBPREFIX).so.$(MAJOR_VERSION) $(EXTRALIB)
	$(CC) $(CFLAGS) $(LDFLAGS) -w -o linktest linktest.c ../$(LIBSONAME) $(FEXTRALIB) && echo OK.
else
#for LSB
	env LSBCC_SHAREDLIBS=gfortran $(CC) $(CFLAGS) $(LDFLAGS) -shared -o ../$(LIBSONAME) \
	-Wl,--whole-archive ../$(LIBNAME) -Wl,--no-whole-archive \
	-Wl,-soname,$(LIBPREFIX).so.$(MAJOR_VERSION) $(EXTRALIB)
	$(FC) $(CFLAGS) $(LDFLAGS) -w -o linktest linktest.c ../$(LIBSONAME) $(FEXTRALIB) && echo OK.
endif
	rm -f linktest

endif

#http://stackoverflow.com/questions/7656425/makefile-ifeq-logical-or
ifeq ($(OSNAME), $(filter $(OSNAME),FreeBSD NetBSD))

so : ../$(LIBSONAME)

../$(LIBSONAME) : ../$(LIBNAME) linktest.c
	$(CC) $(CFLAGS) $(LDFLAGS)  -shared -o ../$(LIBSONAME) \
	-Wl,--whole-archive ../$(LIBNAME) -Wl,--no-whole-archive \
	$(FEXTRALIB) $(EXTRALIB)
	$(CC) $(CFLAGS) $(LDFLAGS) -w -o linktest linktest.c ../$(LIBSONAME) $(FEXTRALIB) && echo OK.
	rm -f linktest

endif

ifeq ($(OSNAME), OSF1)

so : ../$(LIBSONAME)

../$(LIBSONAME) :
	$(CC) $(CFLAGS) $(LDFLAGS) -shared -o ../$(LIBSONAME) ../$(LIBNAME)
endif

ifeq ($(OSNAME), SunOS)

so : ../$(LIBSONAME)
	$(CC) $(CFLAGS) $(LDFLAGS)  -shared -o ../$(LIBSONAME) \
	-Wl,--whole-archive ../$(LIBNAME) -Wl,--no-whole-archive $(EXTRALIB)
	$(CC) $(CFLAGS) $(LDFLAGS) -w -o linktest linktest.c ../$(LIBSONAME) $(FEXTRALIB) && echo OK.
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

osx.def : gensymbol ../Makefile.system ../getarch.c
	perl ./gensymbol osx $(ARCH) $(BU) $(EXPRECISION) $(NO_CBLAS)  $(NO_LAPACK) $(NO_LAPACKE) $(NEED2UNDERSCORES) $(ONLY_CBLAS) > $(@F)

aix.def : gensymbol ../Makefile.system ../getarch.c
	perl ./gensymbol aix $(ARCH) $(BU) $(EXPRECISION) $(NO_CBLAS)  $(NO_LAPACK) $(NO_LAPACKE) $(NEED2UNDERSCORES) $(ONLY_CBLAS) > $(@F)

test : linktest.c
	$(CC) $(CFLAGS) $(LDFLAGS) -w -o linktest linktest.c ../$(LIBSONAME) -lm && echo OK.
	rm -f linktest

linktest.c : gensymbol ../Makefile.system ../getarch.c
	perl ./gensymbol linktest  $(ARCH) $(BU) $(EXPRECISION) $(NO_CBLAS) $(NO_LAPACK) $(NO_LAPACKE) $(NEED2UNDERSCORES) $(ONLY_CBLAS) > linktest.c

clean ::
	@rm -f *.def *.dylib __.SYMDEF*

include ../Makefile.tail


