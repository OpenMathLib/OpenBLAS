#######################################################################
#  This makefile runs the test programs for the linear equation routines
#  and the eigenvalue routines in LAPACK.  The test output files
#  are grouped as follows:
#
#       SLINTST,SEIGTST  -- Single precision real test routines
#       CLINTST,CEIGTST  -- Single precision complex test routines
#       DLINTST,DEIGTST  -- Double precision real test routines
#       ZLINTST,ZEIGTST  -- Double precision complex test routines
#
#  Test programs can be executed for all or some of the four different
#  precisions.  Enter 'make' followed by one or more of the data types
#  desired.
#  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments runs all eight test programs.
#  The executable files are called:
#       xlintsts, xlintstd, xlintstc, and xlintstz for LIN
#       xeigtsts, xeigtstd, xeigtstc, and xeigtstz for EIG
#  and exist in the current directory level.
#
#  To remove the output files after the tests have been run, enter
#       make clean
#
#  To re-run specific tests after a make, enter (for example):
#       'rm ssvd.out; make'  or:
#       'make ssvd.out' or:
#       'touch svd.in; make' (to re-run the single precision SVD tests.)
#
#       'rm *svd.out; make'  (to re-run all the SVD tests.)
#
#######################################################################

include ../make.inc

ifneq ($(strip $(VARLIB)),)
    LAPACKLIB :=  $(VARLIB) ../$(LAPACKLIB)
endif


all:     single complex double complex16 singleproto doubleproto complexproto complex16proto

SEIGTST= snep.out \
         ssep.out \
         ssvd.out \
         sec.out \
         sed.out \
         sgg.out \
         sgd.out \
         ssb.out \
	 ssg.out \
         sbal.out \
         sbak.out \
         sgbal.out \
         sgbak.out \
         sbb.out \
         sglm.out \
         sgqr.out \
         sgsv.out \
         scsd.out \
         slse.out

CEIGTST= cnep.out \
         csep.out \
         csvd.out \
         cec.out \
         ced.out \
         cgg.out \
         cgd.out \
         csb.out \
	 csg.out \
         cbal.out \
         cbak.out \
         cgbal.out \
         cgbak.out \
         cbb.out \
         cglm.out \
         cgqr.out \
         cgsv.out \
         ccsd.out \
         clse.out

DEIGTST= dnep.out \
         dsep.out \
         dsvd.out \
         dec.out \
         ded.out \
         dgg.out \
         dgd.out \
         dsb.out \
	 dsg.out \
         dbal.out \
         dbak.out \
         dgbal.out \
         dgbak.out \
         dbb.out \
         dglm.out \
         dgqr.out \
         dgsv.out \
         dcsd.out \
         dlse.out

ZEIGTST= znep.out \
         zsep.out \
         zsvd.out \
         zec.out \
         zed.out \
         zgg.out \
         zgd.out \
         zsb.out \
	 zsg.out \
         zbal.out \
         zbak.out \
         zgbal.out \
         zgbak.out \
         zbb.out \
         zglm.out \
         zgqr.out \
         zgsv.out \
         zcsd.out \
         zlse.out


SLINTST= stest.out

SLINTSTPROTO= stest_rfp.out

CLINTST= ctest.out

CLINTSTPROTO= ctest_rfp.out

DLINTST= dtest.out

DLINTSTPROTO= dstest.out dtest_rfp.out

ZLINTST= ztest.out

ZLINTSTPROTO= zctest.out ztest_rfp.out

single:         $(SLINTST) $(SEIGTST)
complex:        $(CLINTST) $(CEIGTST)
double:         $(DLINTST) $(DEIGTST)
complex16:      $(ZLINTST) $(ZEIGTST)
singleproto:    $(SLINTSTPROTO)
complexproto:   $(CLINTSTPROTO)
doubleproto:    $(DLINTSTPROTO)
complex16proto: $(ZLINTSTPROTO)

#
# ======== SINGLE LIN TESTS ===========================

stest.out: stest.in xlintsts
	@echo Testing REAL LAPACK linear equation routines
	./xlintsts < stest.in > $@ 2>&1
#
# ======== COMPLEX LIN TESTS ==========================

ctest.out: ctest.in xlintstc
	@echo Testing COMPLEX LAPACK linear equation routines
	./xlintstc < ctest.in > $@ 2>&1
#
# ======== DOUBLE LIN TESTS ===========================

dtest.out: dtest.in xlintstd
	@echo Testing DOUBLE PRECISION LAPACK linear equation routines
	./xlintstd < dtest.in > $@ 2>&1
#
# ======== COMPLEX16 LIN TESTS ========================

ztest.out: ztest.in xlintstz
	@echo Testing COMPLEX16 LAPACK linear equation routines
	./xlintstz < ztest.in > $@ 2>&1
#
# ======== SINGLE-DOUBLE PROTO LIN TESTS ==============

dstest.out: dstest.in xlintstds
	@echo Testing SINGLE-DOUBLE PRECISION LAPACK prototype linear equation routines
	./xlintstds < dstest.in > $@ 2>&1
#
# ======== COMPLEX-COMPLEX16 LIN TESTS ========================

zctest.out: zctest.in xlintstzc
	@echo Testing COMPLEX-COMPLEX16 LAPACK protoype linear equation routines
	./xlintstzc < zctest.in > $@ 2>&1
#
# ======== SINGLE RFP LIN TESTS ========================

stest_rfp.out: stest_rfp.in xlintstrfs
	@echo Testing REAL LAPACK RFP protoype linear equation routines
	./xlintstrfs < stest_rfp.in > $@ 2>&1
#
# ======== COMPLEX16 RFP LIN TESTS ========================

dtest_rfp.out: dtest_rfp.in xlintstrfd
	@echo Testing DOUBLE PRECISION LAPACK RFP protoype linear equation routines
	./xlintstrfd < dtest_rfp.in > $@ 2>&1
#
# ======== COMPLEX16 RFP LIN TESTS ========================

ctest_rfp.out: ctest_rfp.in xlintstrfc
	@echo Testing COMPLEX LAPACK RFP protoype linear equation routines
	./xlintstrfc < ctest_rfp.in > $@ 2>&1
#
# ======== COMPLEX16 RFP LIN TESTS ========================

ztest_rfp.out: ztest_rfp.in xlintstrfz
	@echo Testing COMPLEX16 LAPACK RFP protoype linear equation routines
	./xlintstrfz < ztest_rfp.in > $@ 2>&1
#
#
# ======== SINGLE EIG TESTS ===========================
#

snep.out: nep.in xeigtsts
	@echo NEP: Testing Nonsymmetric Eigenvalue Problem routines
	./xeigtsts < nep.in > $@ 2>&1

ssep.out: sep.in xeigtsts
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtsts < sep.in > $@ 2>&1

ssvd.out: svd.in xeigtsts
	@echo SVD: Testing Singular Value Decomposition routines
	./xeigtsts < svd.in > $@ 2>&1

sec.out: sec.in xeigtsts
	@echo SEC: Testing REAL Eigen Condition Routines
	./xeigtsts < sec.in > $@ 2>&1

sed.out: sed.in xeigtsts
	@echo SEV: Testing REAL Nonsymmetric Eigenvalue Driver
	./xeigtsts < sed.in > $@ 2>&1

sgg.out: sgg.in xeigtsts
	@echo SGG: Testing REAL Nonsymmetric Generalized Eigenvalue Problem routines
	./xeigtsts < sgg.in > $@ 2>&1

sgd.out: sgd.in xeigtsts
	@echo SGD: Testing REAL Nonsymmetric Generalized Eigenvalue Problem driver routines
	./xeigtsts < sgd.in > $@ 2>&1

ssb.out: ssb.in xeigtsts
	@echo SSB: Testing REAL Symmetric Eigenvalue Problem routines
	./xeigtsts < ssb.in > $@ 2>&1

ssg.out: ssg.in xeigtsts
	@echo SSG: Testing REAL Symmetric Generalized Eigenvalue Problem routines
	./xeigtsts < ssg.in > $@ 2>&1

sbal.out: sbal.in xeigtsts
	@echo SGEBAL: Testing the balancing of a REAL general matrix
	./xeigtsts < sbal.in > $@ 2>&1

sbak.out: sbak.in xeigtsts
	@echo SGEBAK: Testing the back transformation of a REAL balanced matrix
	./xeigtsts < sbak.in > $@ 2>&1

sgbal.out: sgbal.in xeigtsts
	@echo SGGBAL: Testing the balancing of a pair of REAL general matrices
	./xeigtsts < sgbal.in > $@ 2>&1

sgbak.out: sgbak.in xeigtsts
	@echo SGGBAK: Testing the back transformation of a pair of REAL balanced matrices
	./xeigtsts < sgbak.in > $@ 2>&1

sbb.out: sbb.in xeigtsts
	@echo SBB:  Testing banded Singular Value Decomposition routines
	./xeigtsts < sbb.in > $@ 2>&1

sglm.out: glm.in xeigtsts
	@echo GLM: Testing Generalized Linear Regression Model routines
	./xeigtsts < glm.in > $@ 2>&1

sgqr.out: gqr.in xeigtsts
	@echo GQR: Testing Generalized QR and RQ factorization routines
	./xeigtsts < gqr.in > $@ 2>&1

sgsv.out: gsv.in xeigtsts
	@echo GSV: Testing Generalized Singular Value Decomposition routines
	./xeigtsts < gsv.in > $@ 2>&1

scsd.out: csd.in xeigtsts
	@echo CSD: Testing CS Decomposition routines
	./xeigtsts < csd.in > $@ 2>&1

slse.out: lse.in xeigtsts
	@echo LSE: Testing Constrained Linear Least Squares routines
	./xeigtsts < lse.in > $@ 2>&1
#
# ======== COMPLEX EIG TESTS ===========================

cnep.out: nep.in xeigtstc
	@echo NEP: Testing Nonsymmetric Eigenvalue Problem routines
	./xeigtstc < nep.in > $@ 2>&1

csep.out: sep.in xeigtstc
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstc < sep.in > $@ 2>&1

csvd.out: svd.in xeigtstc
	@echo SVD: Testing Singular Value Decomposition routines
	./xeigtstc < svd.in > $@ 2>&1

cec.out: cec.in xeigtstc
	@echo CEC: Testing COMPLEX Eigen Condition Routines
	./xeigtstc < cec.in > $@ 2>&1

ced.out: ced.in xeigtstc
	@echo CES: Testing COMPLEX Nonsymmetric Schur Form Driver
	./xeigtstc < ced.in > $@ 2>&1

cgg.out: cgg.in xeigtstc
	@echo CGG: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem routines
	./xeigtstc < cgg.in > $@ 2>&1

cgd.out: cgd.in xeigtstc
	@echo CGD: Testing COMPLEX Nonsymmetric Generalized Eigenvalue Problem driver routines
	./xeigtstc < cgd.in > $@ 2>&1

csb.out: csb.in xeigtstc
	@echo CHB: Testing Hermitian Eigenvalue Problem routines
	./xeigtstc < csb.in > $@ 2>&1

csg.out: csg.in xeigtstc
	@echo CSG: Testing Symmetric Generalized Eigenvalue Problem routines
	./xeigtstc < csg.in > $@ 2>&1

cbal.out: cbal.in xeigtstc
	@echo CGEBAL: Testing the balancing of a COMPLEX general matrix
	./xeigtstc < cbal.in > $@ 2>&1

cbak.out: cbak.in xeigtstc
	@echo CGEBAK: Testing the back transformation of a COMPLEX balanced matrix
	./xeigtstc < cbak.in > $@ 2>&1

cgbal.out: cgbal.in xeigtstc
	@echo CGGBAL: Testing the balancing of a pair of COMPLEX general matrices
	./xeigtstc < cgbal.in > $@ 2>&1

cgbak.out: cgbak.in xeigtstc
	@echo CGGBAK: Testing the back transformation of a pair of COMPLEX balanced matrices
	./xeigtstc < cgbak.in > $@ 2>&1

cbb.out: cbb.in xeigtstc
	@echo CBB:  Testing banded Singular Value Decomposition routines
	./xeigtstc < cbb.in > $@ 2>&1

cglm.out: glm.in xeigtstc
	@echo GLM: Testing Generalized Linear Regression Model routines
	./xeigtstc < glm.in > $@ 2>&1

cgqr.out: gqr.in xeigtstc
	@echo GQR: Testing Generalized QR and RQ factorization routines
	./xeigtstc < gqr.in > $@ 2>&1

cgsv.out: gsv.in xeigtstc
	@echo GSV: Testing Generalized Singular Value Decomposition routines
	./xeigtstc < gsv.in > $@ 2>&1

ccsd.out: csd.in xeigtstc
	@echo CSD: Testing CS Decomposition routines
	./xeigtstc < csd.in > $@ 2>&1

clse.out: lse.in xeigtstc
	@echo LSE: Testing Constrained Linear Least Squares routines
	./xeigtstc < lse.in > $@ 2>&1
#
# ======== DOUBLE EIG TESTS ===========================

dnep.out: nep.in xeigtstd
	@echo NEP: Testing Nonsymmetric Eigenvalue Problem routines
	./xeigtstd < nep.in > $@ 2>&1

dsep.out: sep.in xeigtstd
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstd < sep.in > $@ 2>&1

dsvd.out: svd.in xeigtstd
	@echo SVD: Testing Singular Value Decomposition routines
	./xeigtstd < svd.in > $@ 2>&1

dec.out: dec.in xeigtstd
	@echo DEC: Testing DOUBLE PRECISION Eigen Condition Routines
	./xeigtstd < dec.in > $@ 2>&1

ded.out: ded.in xeigtstd
	@echo DEV: Testing DOUBLE PRECISION Nonsymmetric Eigenvalue Driver
	./xeigtstd < ded.in > $@ 2>&1

dgg.out: dgg.in xeigtstd
	@echo DGG: Testing DOUBLE PRECISION Nonsymmetric Generalized Eigenvalue Problem routines
	./xeigtstd < dgg.in > $@ 2>&1

dgd.out: dgd.in xeigtstd
	@echo DGD: Testing DOUBLE PRECISION Nonsymmetric Generalized Eigenvalue Problem driver routines
	./xeigtstd < dgd.in > $@ 2>&1

dsb.out: dsb.in xeigtstd
	@echo DSB: Testing DOUBLE PRECISION Symmetric Eigenvalue Problem routines
	./xeigtstd < dsb.in > $@ 2>&1

dsg.out: dsg.in xeigtstd
	@echo DSG: Testing DOUBLE PRECISION Symmetric Generalized Eigenvalue Problem routines
	./xeigtstd < dsg.in > $@ 2>&1

dbal.out: dbal.in xeigtstd
	@echo DGEBAL: Testing the balancing of a DOUBLE PRECISION general matrix
	./xeigtstd < dbal.in > $@ 2>&1

dbak.out: dbak.in xeigtstd
	@echo DGEBAK:  Testing the back transformation of a DOUBLE PRECISION balanced matrix
	./xeigtstd < dbak.in > $@ 2>&1

dgbal.out: dgbal.in xeigtstd
	@echo DGGBAL: Testing the balancing of a pair of DOUBLE PRECISION general matrices
	./xeigtstd < dgbal.in > $@ 2>&1

dgbak.out: dgbak.in xeigtstd
	@echo DGGBAK: Testing the back transformation of a pair of DOUBLE PRECISION balanced matrices
	./xeigtstd < dgbak.in > $@ 2>&1

dbb.out: dbb.in xeigtstd
	@echo DBB:  Testing banded Singular Value Decomposition routines
	./xeigtstd < dbb.in > $@ 2>&1

dglm.out: glm.in xeigtstd
	@echo GLM: Testing Generalized Linear Regression Model routines
	./xeigtstd < glm.in > $@ 2>&1

dgqr.out: gqr.in xeigtstd
	@echo GQR: Testing Generalized QR and RQ factorization routines
	./xeigtstd < gqr.in > $@ 2>&1

dgsv.out: gsv.in xeigtstd
	@echo GSV: Testing Generalized Singular Value Decomposition routines
	./xeigtstd < gsv.in > $@ 2>&1

dcsd.out: csd.in xeigtstd
	@echo CSD: Testing CS Decomposition routines
	./xeigtstd < csd.in > $@ 2>&1

dlse.out: lse.in xeigtstd
	@echo LSE: Testing Constrained Linear Least Squares routines
	./xeigtstd < lse.in > $@ 2>&1
#
# ======== COMPLEX16 EIG TESTS ===========================

znep.out: nep.in xeigtstz
	@echo NEP: Testing Nonsymmetric Eigenvalue Problem routines
	./xeigtstz < nep.in > $@ 2>&1

zsep.out: sep.in xeigtstz
	@echo SEP: Testing Symmetric Eigenvalue Problem routines
	./xeigtstz < sep.in > $@ 2>&1

zsvd.out: svd.in xeigtstz
	@echo SVD: Testing Singular Value Decomposition routines
	./xeigtstz < svd.in > $@ 2>&1

zec.out: zec.in xeigtstz
	@echo ZEC: Testing COMPLEX16 Eigen Condition Routines
	./xeigtstz < zec.in > $@ 2>&1

zed.out: zed.in xeigtstz
	@echo ZES: Testing COMPLEX16 Nonsymmetric Schur Form Driver
	./xeigtstz < zed.in > $@ 2>&1

zgg.out: zgg.in xeigtstz
	@echo ZGG: Testing COMPLEX16 Nonsymmetric Generalized Eigenvalue Problem routines
	./xeigtstz < zgg.in > $@ 2>&1

zgd.out: zgd.in xeigtstz
	@echo ZGD: Testing COMPLEX16 Nonsymmetric Generalized Eigenvalue Problem driver routines
	./xeigtstz < zgd.in > $@ 2>&1

zsb.out: zsb.in xeigtstz
	@echo ZHB: Testing Hermitian Eigenvalue Problem routines
	./xeigtstz < zsb.in > $@ 2>&1

zsg.out: zsg.in xeigtstz
	@echo ZSG: Testing Symmetric Generalized Eigenvalue Problem routines
	./xeigtstz < zsg.in > $@ 2>&1

zbal.out: zbal.in xeigtstz
	@echo ZGEBAL: Testing the balancing of a COMPLEX16 general matrix
	./xeigtstz < zbal.in > $@ 2>&1

zbak.out: zbak.in xeigtstz
	@echo ZGEBAK: Testing the back transformation of a COMPLEX16 balanced matrix
	./xeigtstz < zbak.in > $@ 2>&1

zgbal.out: zgbal.in xeigtstz
	@echo ZGGBAL: Testing the balancing of a pair of COMPLEX general matrices
	./xeigtstz < zgbal.in > $@ 2>&1

zgbak.out: zgbak.in xeigtstz
	@echo ZGGBAK: Testing the back transformation of a pair of COMPLEX16 balanced matrices
	./xeigtstz < zgbak.in > $@ 2>&1

zbb.out: zbb.in xeigtstz
	@echo ZBB: Testing banded Singular Value Decomposition routines
	./xeigtstz < zbb.in > $@ 2>&1

zglm.out: glm.in xeigtstz
	@echo GLM: Testing Generalized Linear Regression Model routines
	./xeigtstz < glm.in > $@ 2>&1

zgqr.out: gqr.in xeigtstz
	@echo GQR: Testing Generalized QR and RQ factorization routines
	./xeigtstz < gqr.in > $@ 2>&1

zgsv.out: gsv.in xeigtstz
	@echo GSV: Testing Generalized Singular Value Decomposition routines
	./xeigtstz < gsv.in > $@ 2>&1

zcsd.out: csd.in xeigtstz
	@echo CSD: Testing CS Decomposition routines
	./xeigtstz < csd.in > $@ 2>&1

zlse.out: lse.in xeigtstz
	@echo LSE: Testing Constrained Linear Least Squares routines
	./xeigtstz < lse.in > $@ 2>&1
# ==============================================================================

xlintsts:	../$(LAPACKLIB) ../$(TMGLIB) $(FRCLIN) $(FRC)
	cd LIN ; $(MAKE) single

xlintstc:	../$(LAPACKLIB) ../$(TMGLIB) $(FRCLIN) $(FRC)
	cd LIN ; $(MAKE) complex

xlintstd:	../$(LAPACKLIB) ../$(TMGLIB) $(FRCLIN) $(FRC)
	cd LIN ; $(MAKE) double

xlintstz:	../$(LAPACKLIB) ../$(TMGLIB) $(FRCLIN) $(FRC)
	cd LIN ; $(MAKE) complex16

xlintstrfs:	../$(LAPACKLIB) ../$(TMGLIB) $(FRCLIN) $(FRC)
	cd LIN ; $(MAKE) proto-single

xlintstrfc:	../$(LAPACKLIB) ../$(TMGLIB) $(FRCLIN) $(FRC)
	cd LIN ; $(MAKE) proto-complex

xlintstrfd:	../$(LAPACKLIB) ../$(TMGLIB) $(FRCLIN) $(FRC)
	cd LIN ; $(MAKE) proto-double

xlintstrfz:	../$(LAPACKLIB) ../$(TMGLIB) $(FRCLIN) $(FRC)
	cd LIN ; $(MAKE) proto-complex16

xlintstds:	../$(LAPACKLIB) ../$(TMGLIB) $(FRCLIN) $(FRC)
	cd LIN ; $(MAKE) proto-double

xlintstzc:	../$(LAPACKLIB) ../$(TMGLIB) $(FRCLIN) $(FRC)
	cd LIN ; $(MAKE) proto-complex16

xeigtsts:	../$(LAPACKLIB) ../$(TMGLIB) $(FRCEIG) $(FRC)
	cd EIG ; $(MAKE) single

xeigtstc:	../$(LAPACKLIB) ../$(TMGLIB) $(FRCEIG) $(FRC)
	cd EIG ; $(MAKE) complex

xeigtstd:	../$(LAPACKLIB) ../$(TMGLIB) $(FRCEIG) $(FRC)
	cd EIG ; $(MAKE) double

xeigtstz:	../$(LAPACKLIB) ../$(TMGLIB) $(FRCEIG) $(FRC)
	cd EIG ; $(MAKE) complex16

clean:
	rm -f *.out core

cleanup:
	rm -f x* *.out core

FRCLIN:
	@FRCLIN=$(FRCLIN)

FRCEIG:
	@FRCEIG=$(FRCEIG)

FRC:
	@FRC=$(FRC)
