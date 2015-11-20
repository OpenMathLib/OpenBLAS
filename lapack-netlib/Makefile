#
#  Top Level Makefile for LAPACK
#  Version 3.4.1
#  April 2012
#

include make.inc

all: lapack_install lib blas_testing lapack_testing 

lib: lapacklib tmglib
#lib: blaslib variants lapacklib tmglib

clean: cleanlib cleantesting cleanblas_testing cleancblas_testing

lapack_install:
	( cd INSTALL; $(MAKE); )

blaslib:
	( cd BLAS/SRC; $(MAKE) )

cblaslib:
	( cd CBLAS/src; $(MAKE) )

lapacklib:	lapack_install
	( cd SRC; $(MAKE) )

lapackelib: lapacklib
	( cd LAPACKE; $(MAKE) )

cblas_example: cblaslib blaslib
	( cd CBLAS/examples; $(MAKE) )

lapacke_example: lapackelib
	( cd LAPACKE/example; $(MAKE) )

variants:
	( cd SRC/VARIANTS ; $(MAKE))

tmglib:
	( cd TESTING/MATGEN; $(MAKE) )

lapack_testing:	lib
	( cd TESTING ; $(MAKE) )
	./lapack_testing.py

variants_testing: lib variants
	( cd TESTING ; rm -f xlintst* ; $(MAKE)  VARLIB='SRC/VARIANTS/LIB/cholrl.a' ; \
	mv stest.out stest_cholrl.out ; mv dtest.out dtest_cholrl.out ; mv ctest.out ctest_cholrl.out ; mv ztest.out ztest_cholrl.out )
	( cd TESTING ; rm -f xlintst* ; $(MAKE)  VARLIB='SRC/VARIANTS/LIB/choltop.a' ; \
	mv stest.out stest_choltop.out ; mv dtest.out dtest_choltop.out ; mv ctest.out ctest_choltop.out ; mv ztest.out ztest_choltop.out )
	( cd TESTING ; rm -f xlintst* ; $(MAKE)  VARLIB='SRC/VARIANTS/LIB/lucr.a' ; \
	mv stest.out stest_lucr.out ; mv dtest.out dtest_lucr.out ; mv ctest.out ctest_lucr.out ; mv ztest.out ztest_lucr.out )
	( cd TESTING ;  rm -f xlintst* ; $(MAKE)  VARLIB='SRC/VARIANTS/LIB/lull.a' ; \
	mv stest.out stest_lull.out ; mv dtest.out dtest_lull.out ; mv ctest.out ctest_lull.out ; mv ztest.out ztest_lull.out )
	( cd TESTING ;  rm -f xlintst* ; $(MAKE)  VARLIB='SRC/VARIANTS/LIB/lurec.a' ; \
	mv stest.out stest_lurec.out ; mv dtest.out dtest_lurec.out ; mv ctest.out ctest_lurec.out ; mv ztest.out ztest_lurec.out )
	( cd TESTING ;  rm -f xlintst* ; $(MAKE)  VARLIB='SRC/VARIANTS/LIB/qrll.a' ; \
	mv stest.out stest_qrll.out ; mv dtest.out dtest_qrll.out ; mv ctest.out ctest_qrll.out ; mv ztest.out ztest_qrll.out )
	
blas_testing:
	( cd BLAS/TESTING; $(MAKE) -f Makeblat1 )
	( cd BLAS; ./xblat1s > sblat1.out    ; \
	           ./xblat1d > dblat1.out    ; \
	           ./xblat1c > cblat1.out    ; \
	           ./xblat1z > zblat1.out    ) 
	( cd BLAS/TESTING; $(MAKE) -f Makeblat2 )
	( cd BLAS; ./xblat2s < sblat2.in     ; \
	           ./xblat2d < dblat2.in     ; \
	           ./xblat2c < cblat2.in     ; \
	           ./xblat2z < zblat2.in     )
	( cd BLAS/TESTING; $(MAKE) -f Makeblat3 )
	( cd BLAS; ./xblat3s < sblat3.in     ; \
	           ./xblat3d < dblat3.in     ; \
	           ./xblat3c < cblat3.in     ; \
	           ./xblat3z < zblat3.in     ) 

cblas_testing: blaslib
	( cd CBLAS ; $(MAKE) cblas_testing)
	( cd CBLAS ; $(MAKE) runtst)



html:
	@echo "LAPACK HTML PAGES GENRATION with Doxygen"
	doxygen DOCS/Doxyfile
	@echo "=================="	
	@echo "LAPACK HTML PAGES GENRATED in DOCS/explore-html"
	@echo "Usage: open DOCS/explore-html/index.html"
	@echo "Online version available at http://www.netlib.org/lapack/explore-html/"
	@echo "=================="

man:
	@echo "LAPACK MAN PAGES GENRATION with Doxygen"
	doxygen DOCS/Doxyfile_man
	@echo "=================="
	@echo "LAPACK MAN PAGES GENRATED in DOCS/MAN"
	@echo "Set your MANPATH env variable accordingly"
	@echo "Usage: man dgetrf.f"
	@echo "=================="

cleanlib:
	( cd INSTALL; $(MAKE) clean )
	( cd BLAS/SRC; $(MAKE) clean )
	( cd CBLAS; $(MAKE) clean )
	( cd SRC; $(MAKE) clean )
	( cd SRC/VARIANTS; $(MAKE) clean )
	( cd TESTING/MATGEN; $(MAKE) clean )
	( cd LAPACKE; $(MAKE) clean )


cleanblas_testing:	
	( cd BLAS/TESTING; $(MAKE) -f Makeblat1 clean )
	( cd BLAS/TESTING; $(MAKE) -f Makeblat2 clean )
	( cd BLAS/TESTING; $(MAKE) -f Makeblat3 clean )
	( cd BLAS; rm -f xblat* )

cleancblas_testing:
	( cd CBLAS; $(MAKE) cleanexe )

cleantesting:
	( cd TESTING/LIN; $(MAKE) clean )
	( cd TESTING/EIG; $(MAKE) clean )
	( cd TESTING; rm -f xlin* xeig* )

cleanall: cleanlib cleanblas_testing cleancblas_testing cleantesting 
	rm -f *.a TESTING/*.out INSTALL/test*  BLAS/*.out

