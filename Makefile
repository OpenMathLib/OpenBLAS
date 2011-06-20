TOPDIR	= .
include ./Makefile.system

BLASDIRS = interface driver/level2 driver/level3 driver/others

ifndef DYNAMIC_ARCH
BLASDIRS += kernel 
endif

ifdef UTEST_CHECK
SANITY_CHECK = 1
endif

ifdef SANITY_CHECK
BLASDIRS += reference
endif

ifndef PREFIX
PREFIX = /opt/OpenBLAS
endif

SUBDIRS	= $(BLASDIRS)
ifneq ($(NO_LAPACK), 1)
SUBDIRS	+= lapack
endif

SUBDIRS_ALL = $(SUBDIRS) test ctest utest exports benchmark ../laswp ../bench

.PHONY : all libs netlib test ctest shared
.NOTPARALLEL : all libs prof lapack-test

all :: libs netlib tests shared
	@echo
	@echo " GotoBLAS build complete."
	@echo
	@echo "  OS               ... $(OSNAME)             "
	@echo "  Architecture     ... $(ARCH)               "
ifndef BINARY64
	@echo "  BINARY           ... 32bit                 "
else
	@echo "  BINARY           ... 64bit                 "
endif
	@echo "  C compiler       ... $(C_COMPILER)  (command line : $(CC))"
	@echo "  Fortran compiler ... $(F_COMPILER)  (command line : $(FC))"
ifneq ($(OSNAME), AIX)
	@echo -n "  Library Name     ... $(LIBNAME)"
else
	@echo "  Library Name     ... $(LIBNAME)"
endif

ifndef SMP
	@echo " (Single threaded)  "
else
	@echo " (Multi threaded; Max num-threads is $(NUM_THREADS))"
endif
	@echo

shared :
ifeq ($(OSNAME), Linux)
	$(MAKE) -C exports so
	-ln -fs $(LIBSONAME) libopenblas.so
endif
ifeq ($(OSNAME), FreeBSD)
	$(MAKE) -C exports so
	-ln -fs $(LIBSONAME) libopenblas.so
endif
ifeq ($(OSNAME), NetBSD)
	$(MAKE) -C exports so
	-ln -fs $(LIBSONAME) libopenblas.so
endif
ifeq ($(OSNAME), Darwin)
	$(MAKE) -C exports dyn
	-ln -fs $(LIBDYNNAME) libopenblas.dylib
endif
ifeq ($(OSNAME), WINNT)
	$(MAKE) -C exports dll
#	-ln -fs $(LIBDLLNAME) libopenblas.dll
endif
ifeq ($(OSNAME), CYGWIN_NT)
	$(MAKE) -C exports dll
	-ln -fs $(LIBDLLNAME) libopenblas.dll
endif

tests :
ifndef NOFORTRAN
ifndef TARGET
ifndef CROSS
	touch $(LIBNAME)
ifndef NO_FBLAS
	$(MAKE) -C test all
ifdef UTEST_CHECK
	$(MAKE) -C utest all
endif
endif
ifndef NO_CBLAS
	$(MAKE) -C ctest all
endif
endif
endif
endif

libs :
ifeq ($(CORE), UNKOWN)
	$(error OpenBLAS: Detecting CPU failed. Please set TARGET explicitly, e.g. make TARGET=your_cpu_target. Please read README for the detail.)
endif
	-ln -fs $(LIBNAME) libopenblas.$(LIBSUFFIX)
	for d in $(SUBDIRS) ; \
	do if test -d $$d; then \
	  $(MAKE) -C $$d $(@F) || exit 1 ; \
	fi; \
	done
ifdef DYNAMIC_ARCH
	  $(MAKE) -C kernel commonlibs || exit 1
	for d in $(DYNAMIC_CORE) ; \
	do  $(MAKE) GOTOBLAS_MAKEFILE= -C kernel TARGET_CORE=$$d kernel || exit 1 ;\
	done
endif
	touch lib.grd

prof : prof_blas prof_lapack

prof_blas :
	ln -fs $(LIBNAME_P) libopenblas_p.$(LIBSUFFIX)
	for d in $(SUBDIRS) ; \
	do if test -d $$d; then \
	  $(MAKE) -C $$d prof || exit 1 ; \
	fi; \
	done
ifdef DYNAMIC_ARCH
	  $(MAKE) -C kernel commonprof || exit 1
endif

blas :
	ln -fs $(LIBNAME) libopenblas.$(LIBSUFFIX)
	for d in $(BLASDIRS) ; \
	do if test -d $$d; then \
	  $(MAKE) -C $$d libs || exit 1 ; \
	fi; \
	done

hpl : 
	ln -fs $(LIBNAME) libopenblas.$(LIBSUFFIX)
	for d in $(BLASDIRS) ../laswp exports ; \
	do if test -d $$d; then \
	  $(MAKE) -C $$d $(@F) || exit 1 ; \
	fi; \
	done
ifdef DYNAMIC_ARCH
	  $(MAKE) -C kernel commonlibs || exit 1
	for d in $(DYNAMIC_CORE) ; \
	do  $(MAKE) GOTOBLAS_MAKEFILE= -C kernel TARGET_CORE=$$d kernel || exit 1 ;\
	done
endif

hpl_p :
	ln -fs $(LIBNAME_P) libopenblas_p.$(LIBSUFFIX)
	for d in $(SUBDIRS) ../laswp exports ; \
	do if test -d $$d; then \
	  $(MAKE) -C $$d $(@F) || exit 1 ; \
	fi; \
	done

ifeq ($(NO_LAPACK), 1)
netlib : 

else
netlib : lapack-3.1.1 patch.for_lapack-3.1.1 lapack-3.1.1/make.inc
ifndef NOFORTRAN
	-@$(MAKE) -C lapack-3.1.1 lapacklib
endif
endif

prof_lapack : lapack-3.1.1 lapack-3.1.1/make.inc
	-@$(MAKE) -C lapack-3.1.1 lapack_prof

lapack-3.1.1/make.inc :
ifndef NOFORTRAN
	-@echo "FORTRAN   = $(FC)" > lapack-3.1.1/make.inc
	-@echo "OPTS      = $(FFLAGS)" >> lapack-3.1.1/make.inc
	-@echo "POPTS     = $(FPFLAGS)" >> lapack-3.1.1/make.inc
	-@echo "NOOPT     = $(FFLAGS) -O0" >> lapack-3.1.1/make.inc
	-@echo "PNOOPT     = $(FPFLAGS) -O0" >> lapack-3.1.1/make.inc
	-@echo "LOADOPTS  = $(FFLAGS) $(EXTRALIB)" >> lapack-3.1.1/make.inc
	-@echo "ARCH      = $(AR)" >> lapack-3.1.1/make.inc
	-@echo "RANLIB    = $(RANLIB)" >> lapack-3.1.1/make.inc
	-@echo "LAPACKLIB = ../$(LIBNAME)" >> lapack-3.1.1/make.inc
	-@echo "LAPACKLIB_P = ../$(LIBNAME_P)" >> lapack-3.1.1/make.inc
	-@echo "SUFFIX     = $(SUFFIX)" >> lapack-3.1.1/make.inc
	-@echo "PSUFFIX    = $(PSUFFIX)" >> lapack-3.1.1/make.inc
#	-@echo "CEXTRALIB  = $(CEXTRALIB)" >> lapack-3.1.1/make.inc
	-@cat  make.inc >> lapack-3.1.1/make.inc
endif

lapack-3.1.1 : lapack-3.1.1.tgz
ifndef NOFORTRAN
	@if test `$(MD5SUM) lapack-3.1.1.tgz | $(AWK) '{print $$1}'` = 00b21551a899bcfbaa7b8443e1faeef9; then \
		echo $(TAR) zxf $< ;\
		$(TAR) zxf $< && (cd lapack-3.1.1; $(PATCH) -p1 < ../patch.for_lapack-3.1.1) ;\
	else \
		echo "	lapack-3.1.1.tgz check sum is wrong (Please use orignal)." ;\
		rm -rf lapack-3.1.1 ;\
	fi
endif

lapack-3.1.1.tgz :
ifndef NOFORTRAN
	-wget http://www.netlib.org/lapack/lapack-3.1.1.tgz
endif

large.tgz : 
ifndef NOFORTRAN
	-wget http://www.netlib.org/lapack/timing/large.tgz
endif

timing.tgz :
ifndef NOFORTRAN
	-wget http://www.netlib.org/lapack/timing/timing.tgz
endif

lapack-timing : lapack-3.1.1 large.tgz timing.tgz
ifndef NOFORTRAN
	(cd lapack-3.1.1; $(TAR) zxf ../timing.tgz TIMING)
	(cd lapack-3.1.1/TIMING; $(TAR) zxf ../../large.tgz )
	make -C lapack-3.1.1 tmglib
	make -C lapack-3.1.1/TIMING
endif


lapack-test :
	$(MAKE) -C lapack-3.1.1 tmglib
	$(MAKE) -C lapack-3.1.1/TESTING xeigtstc xeigtstd xeigtsts xeigtstz xlintstc xlintstd xlintstds xlintsts xlintstz xlintstzc
	@rm	-f lapack-3.1.1/TESTING/*.out
	$(MAKE) -j 1 -C lapack-3.1.1/TESTING
	$(GREP) failed lapack-3.1.1/TESTING/*.out

dummy :

lib.grd :
	$(error OpenBLAS: Please run "make" firstly)

install : 	lib.grd
	@-mkdir -p $(PREFIX)
	@echo Generating openblas_config.h in $(PREFIX)
#for inc 
	@echo \#ifndef OPENBLAS_CONFIG_H > $(PREFIX)/openblas_config.h
	@echo \#define OPENBLAS_CONFIG_H >> $(PREFIX)/openblas_config.h
	@cat config.h >> $(PREFIX)/openblas_config.h
	@echo \#define VERSION \" OpenBLAS $(VERSION) \" >> $(PREFIX)/openblas_config.h
	@cat openblas_config_template.h >> $(PREFIX)/openblas_config.h
	@echo \#endif >> $(PREFIX)/openblas_config.h

	@echo Generating f77blas.h in $(PREFIX)
	@echo \#ifndef OPENBLAS_F77BLAS_H > $(PREFIX)/f77blas.h
	@echo \#define OPENBLAS_F77BLAS_H >> $(PREFIX)/f77blas.h
	@echo \#include \"openblas_config.h\" >> $(PREFIX)/f77blas.h
	@cat common_interface.h >> $(PREFIX)/f77blas.h
	@echo \#endif >> $(PREFIX)/f77blas.h

	@echo Generating cblas.h in $(PREFIX)
	@sed 's/common/openblas_config/g' cblas.h > $(PREFIX)/cblas.h

#for install static library 
	@echo Copy the static library to $(PREFIX)
	@cp $(LIBNAME) $(PREFIX)
	@-ln -fs $(PREFIX)/$(LIBNAME) $(PREFIX)/libopenblas.$(LIBSUFFIX)
#for install shared library 
	@echo Copy the shared library to $(PREFIX)
ifeq ($(OSNAME), Linux)
	-cp $(LIBSONAME) $(PREFIX)
	-ln -fs $(PREFIX)/$(LIBSONAME) $(PREFIX)/libopenblas.so
endif
ifeq ($(OSNAME), FreeBSD)
	-cp $(LIBSONAME) $(PREFIX)
	-ln -fs $(PREFIX)/$(LIBSONAME) $(PREFIX)/libopenblas.so
endif
ifeq ($(OSNAME), NetBSD)
	-cp $(LIBSONAME) $(PREFIX)
	-ln -fs $(PREFIX)/$(LIBSONAME) $(PREFIX)/libopenblas.so
endif
ifeq ($(OSNAME), Darwin)
	-cp $(LIBDYNNAME) $(PREFIX)
	-ln -fs $(PREFIX)/$(LIBDYNNAME) $(PREFIX)/libopenblas.dylib
endif
ifeq ($(OSNAME), WINNT)
	-cp $(LIBDLLNAME) $(PREFIX)
	-ln -fs $(PREFIX)/$(LIBDLLNAME) $(PREFIX)/libopenblas.dll
endif
ifeq ($(OSNAME), CYGWIN_NT)
	-cp $(LIBDLLNAME) $(PREFIX)
	-ln -fs $(PREFIX)/$(LIBDLLNAME) $(PREFIX)/libopenblas.dll
endif

	@echo Install OK!

clean ::
	@for d in $(SUBDIRS_ALL) ; \
	do if test -d $$d; then \
	  $(MAKE) -C $$d $(@F) || exit 1 ; \
	fi; \
	done
ifdef DYNAMIC_ARCH
	@$(MAKE) -C kernel clean
endif
	@rm -f *.$(LIBSUFFIX) *.so *~ *.exe getarch getarch_2nd *.dll *.lib *.$(SUFFIX) *.dwf libopenblas.$(LIBSUFFIX) libopenblas_p.$(LIBSUFFIX) *.lnk myconfig.h
	@rm -f Makefile.conf config.h Makefile_kernel.conf config_kernel.h st* *.dylib
	@if test -d lapack-3.1.1; then \
	echo deleting lapack-3.1.1; \
	rm -rf lapack-3.1.1 ;\
	fi
	@rm -f *.grd
	@echo Done.