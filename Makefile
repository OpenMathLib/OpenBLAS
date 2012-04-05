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

.PHONY : all libs netlib test ctest shared install 
.NOTPARALLEL : all libs prof lapack-test install

all :: libs netlib tests shared
	@echo
	@echo " OpenBLAS build complete."
	@echo
	@echo "  OS               ... $(OSNAME)             "
	@echo "  Architecture     ... $(ARCH)               "
ifndef BINARY64
	@echo "  BINARY           ... 32bit                 "
else
	@echo "  BINARY           ... 64bit                 "
endif
ifdef INTERFACE64
	@echo "  Use 64 bits int    (equivalent to \"-i8\" in Fortran)      "
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

ifeq ($(USE_OPENMP), 1)
	@echo
	@echo " Use OpenMP in the multithreading. Becasue of ignoring OPENBLAS_NUM_THREADS and GOTO_NUM_THREADS flags, "
	@echo " you should use OMP_NUM_THREADS environment variable to control the number of threads."
	@echo
endif

ifeq ($(OSNAME), Darwin)
	@echo "WARNING: If you plan to use the dynamic library $(LIBDYNNAME), you must run:"
	@echo
	@echo "\"make PREFIX=/your_installation_path/ install\"."
	@echo
	@echo "(or set PREFIX in Makefile.rule and run make install."
	@echo "If you want to move the .dylib to a new location later, make sure you change"
	@echo "the internal name of the dylib with:"
	@echo
	@echo "install_name_tool -id /new/absolute/path/to/$(LIBDYNNAME) $(LIBDYNNAME)"
endif
	@echo
	@echo "To install the library, you can run \"make PREFIX=/path/to/your/installation install\"."
	@echo

shared :
ifeq ($(OSNAME), Linux)
	$(MAKE) -C exports so
	-ln -fs $(LIBSONAME) $(LIBPREFIX).so
	-ln -fs $(LIBSONAME) $(LIBPREFIX).so.$(MAJOR_VERSION)
endif
ifeq ($(OSNAME), FreeBSD)
	$(MAKE) -C exports so
	-ln -fs $(LIBSONAME) $(LIBPREFIX).so
endif
ifeq ($(OSNAME), NetBSD)
	$(MAKE) -C exports so
	-ln -fs $(LIBSONAME) $(LIBPREFIX).so
endif
ifeq ($(OSNAME), Darwin)
	$(MAKE) -C exports dyn
	-ln -fs $(LIBDYNNAME) $(LIBPREFIX).dylib
endif
ifeq ($(OSNAME), WINNT)
	$(MAKE) -C exports dll
	-ln -fs $(LIBDLLNAME) $(LIBPREFIX).dll
endif
ifeq ($(OSNAME), CYGWIN_NT)
	$(MAKE) -C exports dll
	-ln -fs $(LIBDLLNAME) $(LIBPREFIX).dll
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
ifeq ($(NOFORTRAN), 1)
	$(error OpenBLAS: Detecting fortran compiler failed. Please install fortran compiler, e.g. gfortran, ifort, openf90.)
endif
	-ln -fs $(LIBNAME) $(LIBPREFIX).$(LIBSUFFIX)
	for d in $(SUBDIRS) ; \
	do if test -d $$d; then \
	  $(MAKE) -C $$d $(@F) || exit 1 ; \
	fi; \
	done
#Save the config files for installation
	cp Makefile.conf Makefile.conf_last
	cp config.h config_last.h
ifdef QUAD_PRECISION
	echo "#define QUAD_PRECISION">> config_last.h
endif
ifeq ($(EXPRECISION), 1)
	echo "#define EXPRECISION">> config_last.h
endif
## 
ifdef DYNAMIC_ARCH
	  $(MAKE) -C kernel commonlibs || exit 1
	for d in $(DYNAMIC_CORE) ; \
	do  $(MAKE) GOTOBLAS_MAKEFILE= -C kernel TARGET_CORE=$$d kernel || exit 1 ;\
	done
	echo DYNAMIC_ARCH=1 >> Makefile.conf_last
endif
	touch lib.grd

prof : prof_blas prof_lapack

prof_blas :
	ln -fs $(LIBNAME_P) $(LIBPREFIX)_p.$(LIBSUFFIX)
	for d in $(SUBDIRS) ; \
	do if test -d $$d; then \
	  $(MAKE) -C $$d prof || exit 1 ; \
	fi; \
	done
ifdef DYNAMIC_ARCH
	  $(MAKE) -C kernel commonprof || exit 1
endif

blas :
	ln -fs $(LIBNAME) $(LIBPREFIX).$(LIBSUFFIX)
	for d in $(BLASDIRS) ; \
	do if test -d $$d; then \
	  $(MAKE) -C $$d libs || exit 1 ; \
	fi; \
	done

hpl : 
	ln -fs $(LIBNAME) $(LIBPREFIX).$(LIBSUFFIX)
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
	ln -fs $(LIBNAME_P) $(LIBPREFIX)_p.$(LIBSUFFIX)
	for d in $(SUBDIRS) ../laswp exports ; \
	do if test -d $$d; then \
	  $(MAKE) -C $$d $(@F) || exit 1 ; \
	fi; \
	done

ifeq ($(NO_LAPACK), 1)
netlib : 

else
netlib : lapack-3.4.0 patch.for_lapack-3.4.0 lapack-3.4.0/make.inc
ifndef NOFORTRAN
	-@$(MAKE) -C lapack-3.4.0 lapacklib
endif
endif

prof_lapack : lapack-3.4.0 lapack-3.4.0/make.inc
	-@$(MAKE) -C lapack-3.4.0 lapack_prof

lapack-3.4.0/make.inc :
ifndef NOFORTRAN
	-@echo "FORTRAN   = $(FC)" > lapack-3.4.0/make.inc
	-@echo "OPTS      = $(FFLAGS)" >> lapack-3.4.0/make.inc
	-@echo "POPTS     = $(FPFLAGS)" >> lapack-3.4.0/make.inc
	-@echo "NOOPT     = $(FFLAGS) -O0" >> lapack-3.4.0/make.inc
	-@echo "PNOOPT     = $(FPFLAGS) -O0" >> lapack-3.4.0/make.inc
	-@echo "LOADOPTS  = $(FFLAGS) $(EXTRALIB)" >> lapack-3.4.0/make.inc
	-@echo "ARCH      = $(AR)" >> lapack-3.4.0/make.inc
	-@echo "RANLIB    = $(RANLIB)" >> lapack-3.4.0/make.inc
	-@echo "LAPACKLIB = ../$(LIBNAME)" >> lapack-3.4.0/make.inc
	-@echo "LAPACKLIB_P = ../$(LIBNAME_P)" >> lapack-3.4.0/make.inc
	-@echo "SUFFIX     = $(SUFFIX)" >> lapack-3.4.0/make.inc
	-@echo "PSUFFIX    = $(PSUFFIX)" >> lapack-3.4.0/make.inc
#	-@echo "CEXTRALIB  = $(CEXTRALIB)" >> lapack-3.4.0/make.inc
	-@cat  make.inc >> lapack-3.4.0/make.inc
endif

lapack-3.4.0 : lapack-3.4.0.tgz
ifndef NOFORTRAN
ifndef NO_LAPACK
	@if test `$(MD5SUM) lapack-3.4.0.tgz | $(AWK) '{print $$1}'` = 02d5706ec03ba885fc246e5fa10d8c70; then \
		echo $(TAR) zxf $< ;\
		$(TAR) zxf $< && (cd lapack-3.4.0; $(PATCH) -p1 < ../patch.for_lapack-3.4.0) ;\
	else \
		rm -rf lapack-3.4.0 ;\
		echo "	Cannot download lapack-3.4.0.tgz or the MD5 check sum is wrong (Please use orignal)."; \
		exit 1; \
	fi
endif
endif

LAPACK_URL=http://www.netlib.org/lapack/lapack-3.4.0.tgz

lapack-3.4.0.tgz :
ifndef NOFORTRAN
ifeq ($(OSNAME), Darwin)
	curl -O $(LAPACK_URL)
else
	wget $(LAPACK_URL)
endif
endif

large.tgz : 
ifndef NOFORTRAN
	-wget http://www.netlib.org/lapack/timing/large.tgz
endif

timing.tgz :
ifndef NOFORTRAN
	-wget http://www.netlib.org/lapack/timing/timing.tgz
endif

lapack-timing : lapack-3.4.0 large.tgz timing.tgz
ifndef NOFORTRAN
	(cd lapack-3.4.0; $(TAR) zxf ../timing.tgz TIMING)
	(cd lapack-3.4.0/TIMING; $(TAR) zxf ../../large.tgz )
	make -C lapack-3.4.0 tmglib
	make -C lapack-3.4.0/TIMING
endif


lapack-test :
	$(MAKE) -C lapack-3.4.0 tmglib
	$(MAKE) -C lapack-3.4.0/TESTING xeigtstc xeigtstd xeigtsts xeigtstz xlintstc xlintstd xlintstds xlintsts xlintstz xlintstzc
	@rm	-f lapack-3.4.0/TESTING/*.out
	$(MAKE) -j 1 -C lapack-3.4.0/TESTING
	$(GREP) failed lapack-3.4.0/TESTING/*.out

dummy :

install :
	$(MAKE) -f Makefile.install install

clean ::
	@for d in $(SUBDIRS_ALL) ; \
	do if test -d $$d; then \
	  $(MAKE) -C $$d $(@F) || exit 1 ; \
	fi; \
	done
#ifdef DYNAMIC_ARCH
	@$(MAKE) -C kernel clean
#endif
	@$(MAKE) -C reference clean
	@rm -f *.$(LIBSUFFIX) *.so *~ *.exe getarch getarch_2nd *.dll *.lib *.$(SUFFIX) *.dwf $(LIBPREFIX).$(LIBSUFFIX) $(LIBPREFIX)_p.$(LIBSUFFIX) $(LIBPREFIX).so.$(MAJOR_VERSION) *.lnk myconfig.h
	@rm -f Makefile.conf config.h Makefile_kernel.conf config_kernel.h st* *.dylib
	@if test -d lapack-3.4.0; then \
	echo deleting lapack-3.4.0; \
	rm -rf lapack-3.4.0 ;\
	fi
	@rm -f *.grd Makefile.conf_last config_last.h
	@echo Done.