include ../make.inc

all:
	cd include && cp cblas_mangling_with_flags.h cblas_mangling.h
	cd src && $(MAKE) all


clean: cleanlib

cleanlib:
	cd src && $(MAKE) clean

cleanexe:
	cd testing && $(MAKE) cleanexe

cleanall: clean cleanexe
	rm -f $(CBLASLIB)
	cd examples && rm -f *.o cblas_ex1 cblas_ex2

cblas_testing:
	cd testing && $(MAKE) all

runtst:
	cd testing && $(MAKE) run
	
example: all
	cd examples && make all
