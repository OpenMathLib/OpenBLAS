###############################################################################
# Copyright (c) 2025, The OpenBLAS Project
# All rights reserved.
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in
#    the documentation and/or other materials provided with the
#    distribution.
# 3. Neither the name of the OpenBLAS project nor the names of
#    its contributors may be used to endorse or promote products
#    derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
###############################################################################
TOPDIR	= ..
include ../Makefile.system
ifeq ($(F_COMPILER),GFORTRAN)
ifneq (, $(filter $(CORE),LOONGSON3R3 LOONGSON3R4))
	override FFLAGS := $(filter_out(-O2 -O3,$(FFLAGS))) -O0
endif
        override FFLAGS += -fno-tree-vectorize
endif
# override CFLAGS += -std=c11 -Wall -Werror

SUPPORT_GEMM3M = 0
XERBLA_TEST = xerbla_test.$(SUFFIX)

define CHECK_TEST_SUMMARY
@status=0; $(GREP) -q -e FATAL -e FAILED "$(1)" || status=$$?; \
case $$status in \
  0) cat "$(1)"; exit 1 ;; \
  1) ;; \
  *) echo "Unable to read test summary: $(1)"; exit $$status ;; \
esac
endef

ifeq ($(ARCH), x86)
SUPPORT_GEMM3M = 1
endif

ifeq ($(ARCH), x86_64)
SUPPORT_GEMM3M = 1
endif

ifeq ($(ARCH), ia64)
SUPPORT_GEMM3M = 1
endif

ifeq ($(ARCH), MIPS)
SUPPORT_GEMM3M = 1
endif

ifeq ($(NOFORTRAN),1)
all ::
else
all :: level1 level2 level3
endif

ifeq ($(BUILD_SINGLE),1)
S1=sblat1$(EXE)
endif
ifeq ($(BUILD_DOUBLE),1)
D1=dblat1$(EXE)
endif
ifeq ($(BUILD_COMPLEX),1)
C1=cblat1$(EXE)
endif
ifeq ($(BUILD_COMPLEX16),1)
Z1=zblat1$(EXE)
endif

level1: $(S1) $(D1) $(C1) $(Z1)

ifneq ($(CROSS), 1)
ifeq ($(BUILD_SINGLE),1)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./sblat1$(EXE)
endif
ifeq ($(BUILD_DOUBLE),1)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./dblat1$(EXE)
endif
ifeq ($(BUILD_COMPLEX),1)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./cblat1$(EXE)
endif
ifeq ($(BUILD_COMPLEX16),1)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./zblat1$(EXE)
endif
ifdef SMP
ifeq ($(USE_OPENMP), 1)
ifeq ($(BUILD_SINGLE),1)
	OMP_NUM_THREADS=2 ./sblat1$(EXE)
endif
ifeq ($(BUILD_DOUBLE),1)
	OMP_NUM_THREADS=2 ./dblat1$(EXE)
endif
ifeq ($(BUILD_COMPLEX),1)
	OMP_NUM_THREADS=2 ./cblat1$(EXE)
endif
ifeq ($(BUILD_COMPLEX16),1)
	OMP_NUM_THREADS=2 ./zblat1$(EXE)
endif
else
ifeq ($(BUILD_SINGLE),1)
	OPENBLAS_NUM_THREADS=2 ./sblat1$(EXE)
endif
ifeq ($(BUILD_DOUBLE),1)
	OPENBLAS_NUM_THREADS=2 ./dblat1$(EXE)
endif
ifeq ($(BUILD_COMPLEX),1)
	OPENBLAS_NUM_THREADS=2 ./cblat1$(EXE)
endif
ifeq ($(BUILD_COMPLEX16),1)
	OPENBLAS_NUM_THREADS=2 ./zblat1$(EXE)
endif
endif
endif
endif

ifeq ($(BUILD_HFLOAT16), 1)
SH2 = test_shgemv$(EXE)
endif
ifeq ($(BUILD_BFLOAT16), 1)
BB2 = test_bgemv$(EXE)
B2 = test_sbgemv$(EXE)
endif
ifeq ($(BUILD_SINGLE),1)
S2=sblat2$(EXE)
endif
ifeq ($(BUILD_DOUBLE),1)
D2=dblat2$(EXE)
endif
ifeq ($(BUILD_COMPLEX),1)
C2=cblat2$(EXE)
endif
ifeq ($(BUILD_COMPLEX16),1)
Z2=zblat2$(EXE)
endif

level2: $(SH2) $(BB2) $(B2) $(S2) $(D2) $(C2) $(Z2)


ifneq ($(CROSS), 1)
	rm -f ?BLAT2.SUMM
ifeq ($(BUILD_BFLOAT16),1)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./test_bgemv$(EXE) > BBLAT2.SUMM
	$(call CHECK_TEST_SUMMARY,BBLAT2.SUMM)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./test_sbgemv$(EXE) > SBBLAT2.SUMM
	$(call CHECK_TEST_SUMMARY,SBBLAT2.SUMM)
endif
ifeq ($(BUILD_HFLOAT16),1)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./test_shgemv$(EXE) > SHBLAT2.SUMM
	$(call CHECK_TEST_SUMMARY,SHBLAT2.SUMM)
endif
ifeq ($(BUILD_SINGLE),1)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./sblat2$(EXE) < ./sblat2.dat
	$(call CHECK_TEST_SUMMARY,SBLAT2.SUMM)
endif
ifeq ($(BUILD_DOUBLE),1)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./dblat2$(EXE) < ./dblat2.dat
	$(call CHECK_TEST_SUMMARY,DBLAT2.SUMM)
endif
ifeq ($(BUILD_COMPLEX),1)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./cblat2$(EXE) < ./cblat2.dat
	$(call CHECK_TEST_SUMMARY,CBLAT2.SUMM)
endif
ifeq ($(BUILD_COMPLEX16),1)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./zblat2$(EXE) < ./zblat2.dat
	$(call CHECK_TEST_SUMMARY,ZBLAT2.SUMM)
endif
ifdef SMP
	rm -f ?BLAT2.SUMM
ifeq ($(USE_OPENMP), 1)
ifeq ($(BUILD_BFLOAT16),1)
	OMP_NUM_THREADS=2 ./test_bgemv$(EXE) > BBLAT2.SUMM
	$(call CHECK_TEST_SUMMARY,BBLAT2.SUMM)
	OMP_NUM_THREADS=2 ./test_sbgemv$(EXE) > SBBLAT2.SUMM
	$(call CHECK_TEST_SUMMARY,SBBLAT2.SUMM)
endif
ifeq ($(BUILD_HFLOAT16),1)
	OMP_NUM_THREADS=2 ./test_shgemv$(EXE) > SHBLAT2.SUMM
	$(call CHECK_TEST_SUMMARY,SHBLAT2.SUMM)
endif
ifeq ($(BUILD_SINGLE),1)
	OMP_NUM_THREADS=2 ./sblat2$(EXE) < ./sblat2.dat
	$(call CHECK_TEST_SUMMARY,SBLAT2.SUMM)
endif
ifeq ($(BUILD_DOUBLE),1)
	OMP_NUM_THREADS=2 ./dblat2$(EXE) < ./dblat2.dat
	$(call CHECK_TEST_SUMMARY,DBLAT2.SUMM)
endif
ifeq ($(BUILD_COMPLEX),1)
	OMP_NUM_THREADS=2 ./cblat2$(EXE) < ./cblat2.dat
	$(call CHECK_TEST_SUMMARY,CBLAT2.SUMM)
endif
ifeq ($(BUILD_COMPLEX16),1)
	OMP_NUM_THREADS=2 ./zblat2$(EXE) < ./zblat2.dat
	$(call CHECK_TEST_SUMMARY,ZBLAT2.SUMM)
endif
else
ifeq ($(BUILD_BFLOAT16),1)
	OMP_NUM_THREADS=2 ./test_bgemv$(EXE) > BBLAT2.SUMM
	$(call CHECK_TEST_SUMMARY,BBLAT2.SUMM)
	OMP_NUM_THREADS=2 ./test_sbgemv$(EXE) > SBBLAT2.SUMM
	$(call CHECK_TEST_SUMMARY,SBBLAT2.SUMM)
endif
ifeq ($(BUILD_HFLOAT16),1)
	OMP_NUM_THREADS=2 ./test_shgemv$(EXE) > SHBLAT2.SUMM
	$(call CHECK_TEST_SUMMARY,SHBLAT2.SUMM)
endif
ifeq ($(BUILD_SINGLE),1)
	OPENBLAS_NUM_THREADS=2 ./sblat2$(EXE) < ./sblat2.dat
	$(call CHECK_TEST_SUMMARY,SBLAT2.SUMM)
endif
ifeq ($(BUILD_DOUBLE),1)
	OPENBLAS_NUM_THREADS=2 ./dblat2$(EXE) < ./dblat2.dat
	$(call CHECK_TEST_SUMMARY,DBLAT2.SUMM)
endif
ifeq ($(BUILD_COMPLEX),1)
	OPENBLAS_NUM_THREADS=2 ./cblat2$(EXE) < ./cblat2.dat
	$(call CHECK_TEST_SUMMARY,CBLAT2.SUMM)
endif
ifeq ($(BUILD_COMPLEX16),1)
	OPENBLAS_NUM_THREADS=2 ./zblat2$(EXE) < ./zblat2.dat
	$(call CHECK_TEST_SUMMARY,ZBLAT2.SUMM)
endif
endif
endif
endif

ifeq ($(BUILD_BFLOAT16),1)
BF3= test_bgemm$(EXE)
B3 = test_sbgemm$(EXE)
endif
ifeq ($(BUILD_HFLOAT16),1)
H3 = test_shgemm$(EXE)
endif
ifeq ($(BUILD_SINGLE),1)
S3=sblat3$(EXE)
endif
ifeq ($(BUILD_DOUBLE),1)
D3=dblat3$(EXE)
endif
ifeq ($(BUILD_COMPLEX),1)
ifeq ($(SUPPORT_GEMM3M),1)
C3=cblat3$(EXE)
C3M=cblat3_3m$(EXE)
else
C3=cblat3$(EXE)
C3M=
endif
endif
ifeq ($(BUILD_COMPLEX16),1)
ifeq ($(SUPPORT_GEMM3M),1)
Z3=zblat3$(EXE)
Z3M=zblat3_3m$(EXE)
else
Z3=zblat3$(EXE)
Z3M=
endif
endif


ifeq ($(SUPPORT_GEMM3M),1)
level3: $(BF3) $(B3) $(H3) $(S3) $(D3) $(C3) $(Z3) level3_3m
else
level3: $(BF3) $(B3) $(H3) $(S3) $(D3) $(C3) $(Z3)
endif

ifneq ($(CROSS), 1)
	rm -f ?BLAT3.SUMM
ifeq ($(BUILD_BFLOAT16),1)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./test_sbgemm$(EXE) > SBBLAT3.SUMM
	$(call CHECK_TEST_SUMMARY,SBBLAT3.SUMM)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./test_bgemm$(EXE) > BBLAT3.SUMM
	$(call CHECK_TEST_SUMMARY,BBLAT3.SUMM)
endif
ifeq ($(BUILD_HFLOAT16),1)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./test_shgemm > SHBLAT3.SUMM
	$(call CHECK_TEST_SUMMARY,SHBLAT3.SUMM)
endif
ifeq ($(BUILD_SINGLE),1)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./sblat3$(EXE) < ./sblat3.dat
	$(call CHECK_TEST_SUMMARY,SBLAT3.SUMM)
endif
ifeq ($(BUILD_DOUBLE),1)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./dblat3$(EXE) < ./dblat3.dat
	$(call CHECK_TEST_SUMMARY,DBLAT3.SUMM)
endif
ifeq ($(BUILD_COMPLEX),1)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./cblat3$(EXE) < ./cblat3.dat
	$(call CHECK_TEST_SUMMARY,CBLAT3.SUMM)
endif
ifeq ($(BUILD_COMPLEX16),1)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./zblat3$(EXE) < ./zblat3.dat
	$(call CHECK_TEST_SUMMARY,ZBLAT3.SUMM)
endif
ifdef SMP
	rm -f ?BLAT3.SUMM
ifeq ($(USE_OPENMP), 1)
ifeq ($(BUILD_BFLOAT16),1)
	OMP_NUM_THREADS=2 ./test_sbgemm$(EXE) > SBBLAT3.SUMM
	$(call CHECK_TEST_SUMMARY,SBBLAT3.SUMM)
	OMP_NUM_THREADS=2 ./test_bgemm$(EXE) > BBLAT3.SUMM
	$(call CHECK_TEST_SUMMARY,BBLAT3.SUMM)
endif
ifeq ($(BUILD_HFLOAT16),1)
	OMP_NUM_THREADS=2 ./test_shgemm > SHBLAT3.SUMM
	$(call CHECK_TEST_SUMMARY,SHBLAT3.SUMM)
endif

ifeq ($(BUILD_SINGLE),1)
	OMP_NUM_THREADS=2 ./sblat3$(EXE) < ./sblat3.dat
	$(call CHECK_TEST_SUMMARY,SBLAT3.SUMM)
endif
ifeq ($(BUILD_DOUBLE),1)
	OMP_NUM_THREADS=2 ./dblat3$(EXE) < ./dblat3.dat
	$(call CHECK_TEST_SUMMARY,DBLAT3.SUMM)
endif
ifeq ($(BUILD_COMPLEX),1)
	OMP_NUM_THREADS=2 ./cblat3$(EXE) < ./cblat3.dat
	$(call CHECK_TEST_SUMMARY,CBLAT3.SUMM)
endif
ifeq ($(BUILD_COMPLEX16),1)
	OMP_NUM_THREADS=2 ./zblat3$(EXE) < ./zblat3.dat
	$(call CHECK_TEST_SUMMARY,ZBLAT3.SUMM)
endif
else
ifeq ($(BUILD_BFLOAT16),1)
	OPENBLAS_NUM_THREADS=2 ./test_sbgemm$(EXE) > SBBLAT3.SUMM
	$(call CHECK_TEST_SUMMARY,SBBLAT3.SUMM)
	OPENBLAS_NUM_THREADS=2 ./test_bgemm$(EXE) > BBLAT3.SUMM
	$(call CHECK_TEST_SUMMARY,BBLAT3.SUMM)
endif
ifeq ($(BUILD_HFLOAT16),1)
	OPENBLAS_NUM_THREADS=2 ./test_shgemm > SHBLAT3.SUMM
	$(call CHECK_TEST_SUMMARY,SHBLAT3.SUMM)
endif
ifeq ($(BUILD_SINGLE),1)
	OPENBLAS_NUM_THREADS=2 ./sblat3$(EXE) < ./sblat3.dat
	$(call CHECK_TEST_SUMMARY,SBLAT3.SUMM)
endif
ifeq ($(BUILD_DOUBLE),1)
	OPENBLAS_NUM_THREADS=2 ./dblat3$(EXE) < ./dblat3.dat
	$(call CHECK_TEST_SUMMARY,DBLAT3.SUMM)
endif
ifeq ($(BUILD_COMPLEX),1)
	OPENBLAS_NUM_THREADS=2 ./cblat3$(EXE) < ./cblat3.dat
	$(call CHECK_TEST_SUMMARY,CBLAT3.SUMM)
endif
ifeq ($(BUILD_COMPLEX16),1)
	OPENBLAS_NUM_THREADS=2 ./zblat3$(EXE) < ./zblat3.dat
	$(call CHECK_TEST_SUMMARY,ZBLAT3.SUMM)
endif
endif
endif
endif


level3_3m: $(C3M) $(Z3M)
ifneq ($(CROSS), 1)
	rm -f ?BLAT3_3M.SUMM
ifeq ($(BUILD_COMPLEX),1)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./cblat3_3m$(EXE) < ./cblat3_3m.dat
	$(call CHECK_TEST_SUMMARY,CBLAT3_3M.SUMM)
endif
ifeq ($(BUILD_COMPLEX16),1)
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./zblat3_3m$(EXE) < ./zblat3_3m.dat
	$(call CHECK_TEST_SUMMARY,ZBLAT3_3M.SUMM)
endif
ifdef SMP
	rm -f ?BLAT3_3M.SUMM
ifeq ($(USE_OPENMP), 1)
ifeq ($(BUILD_COMPLEX),1)
	OMP_NUM_THREADS=2 ./cblat3_3m$(EXE) < ./cblat3_3m.dat
	$(call CHECK_TEST_SUMMARY,CBLAT3_3M.SUMM)
endif
ifeq ($(BUILD_COMPLEX16),1)
	OMP_NUM_THREADS=2 ./zblat3_3m$(EXE) < ./zblat3_3m.dat
	$(call CHECK_TEST_SUMMARY,ZBLAT3_3M.SUMM)
endif
else
ifeq ($(BUILD_COMPLEX),1)
	OPENBLAS_NUM_THREADS=2 ./cblat3_3m$(EXE) < ./cblat3_3m.dat
	$(call CHECK_TEST_SUMMARY,CBLAT3_3M.SUMM)
endif
ifeq ($(BUILD_COMPLEX16),1)
	OPENBLAS_NUM_THREADS=2 ./zblat3_3m$(EXE) < ./zblat3_3m.dat
	$(call CHECK_TEST_SUMMARY,ZBLAT3_3M.SUMM)
endif
endif
endif
endif




FLDFLAGS = $(FFLAGS:-fPIC=) $(LDFLAGS)
CLDFLAGS = $(CFLAGS) $(LDFLAGS)


ifeq ($(USE_OPENMP), 1)
ifeq ($(F_COMPILER), GFORTRAN)
ifeq ($(C_COMPILER), CLANG)
EXTRALIB += -lomp
endif
endif
ifeq ($(F_COMPILER), NAG)
EXTRALIB = -lgomp
endif
ifeq ($(F_COMPILER), IBM)
ifeq ($(C_COMPILER), GCC)
EXTRALIB += -lgomp
endif
ifeq ($(C_COMPILER), CLANG)
EXTRALIB += -lomp
endif
endif
endif

$(XERBLA_TEST): xerbla_test.c
	$(CC) $(CFLAGS) -c $< -o $@

ifeq ($(BUILD_SINGLE),1)
sblat1$(EXE) : sblat1.$(SUFFIX) ../$(LIBNAME)
	$(FC) $(FLDFLAGS) -o $@ sblat1.$(SUFFIX) ../$(LIBNAME) $(EXTRALIB) $(CEXTRALIB)

sblat2$(EXE) : sblat2.$(SUFFIX) $(XERBLA_TEST) ../$(LIBNAME)
	$(FC) $(FLDFLAGS) -o $@ sblat2.$(SUFFIX) $(XERBLA_TEST) ../$(LIBNAME) $(EXTRALIB) $(CEXTRALIB)

sblat3$(EXE) : sblat3.$(SUFFIX) $(XERBLA_TEST) ../$(LIBNAME)
	$(FC) $(FLDFLAGS) -o $@ sblat3.$(SUFFIX) $(XERBLA_TEST) ../$(LIBNAME) $(EXTRALIB) $(CEXTRALIB)
endif

ifeq ($(BUILD_DOUBLE),1)
dblat1$(EXE) : dblat1.$(SUFFIX) ../$(LIBNAME)
	$(FC) $(FLDFLAGS) -o $@ dblat1.$(SUFFIX) ../$(LIBNAME) $(EXTRALIB) $(CEXTRALIB)

dblat2$(EXE) : dblat2.$(SUFFIX) $(XERBLA_TEST) ../$(LIBNAME)
	$(FC) $(FLDFLAGS) -o $@ dblat2.$(SUFFIX) $(XERBLA_TEST) ../$(LIBNAME) $(EXTRALIB) $(CEXTRALIB)

dblat3$(EXE) : dblat3.$(SUFFIX) $(XERBLA_TEST) ../$(LIBNAME)
	$(FC) $(FLDFLAGS) -o $@ dblat3.$(SUFFIX) $(XERBLA_TEST) ../$(LIBNAME) $(EXTRALIB) $(CEXTRALIB)
else
dblat2$(EXE):
dblat3$(EXE):
endif


qblat1$(EXE) : qblat1.$(SUFFIX) ../$(LIBNAME)
	$(FC) $(FLDFLAGS) -o $@ qblat1.$(SUFFIX) ../$(LIBNAME) $(EXTRALIB) $(CEXTRALIB)

ifeq ($(BUILD_COMPLEX),1)
cblat1$(EXE) : cblat1.$(SUFFIX) ../$(LIBNAME)
	$(FC) $(FLDFLAGS) -o $@ cblat1.$(SUFFIX) ../$(LIBNAME) $(EXTRALIB) $(CEXTRALIB)

cblat2$(EXE) : cblat2.$(SUFFIX) $(XERBLA_TEST) ../$(LIBNAME)
	$(FC) $(FLDFLAGS) -o $@ cblat2.$(SUFFIX) $(XERBLA_TEST) ../$(LIBNAME) $(EXTRALIB) $(CEXTRALIB)

cblat3$(EXE) : cblat3.$(SUFFIX) $(XERBLA_TEST) ../$(LIBNAME)
	$(FC) $(FLDFLAGS) -o $@ cblat3.$(SUFFIX) $(XERBLA_TEST) ../$(LIBNAME) $(EXTRALIB) $(CEXTRALIB)
endif

ifeq ($(BUILD_COMPLEX16),1)
zblat1$(EXE) : zblat1.$(SUFFIX) ../$(LIBNAME)
	$(FC) $(FLDFLAGS) -o $@ zblat1.$(SUFFIX) ../$(LIBNAME) $(EXTRALIB) $(CEXTRALIB)

zblat2$(EXE) : zblat2.$(SUFFIX) $(XERBLA_TEST) ../$(LIBNAME)
	$(FC) $(FLDFLAGS) -o $@ zblat2.$(SUFFIX) $(XERBLA_TEST) ../$(LIBNAME) $(EXTRALIB) $(CEXTRALIB)

zblat3$(EXE) : zblat3.$(SUFFIX) $(XERBLA_TEST) ../$(LIBNAME)
	$(FC) $(FLDFLAGS) -o $@ zblat3.$(SUFFIX) $(XERBLA_TEST) ../$(LIBNAME) $(EXTRALIB) $(CEXTRALIB)
endif

ifeq ($(BUILD_BFLOAT16),1)
test_bgemm$(EXE) : compare_sgemm_bgemm.c test_helpers.h ../$(LIBNAME)
	$(CC) $(CLDFLAGS) -DIBFLOAT16 -DOBFLOAT16 -o test_bgemm compare_sgemm_bgemm.c ../$(LIBNAME) $(EXTRALIB) $(FEXTRALIB) $(CEXTRALIB)

test_bgemv$(EXE) : compare_sgemv_bgemv.c ../$(LIBNAME)
	$(CC) $(CLDFLAGS) -DIBFLOAT16 -DOBFLOAT16 -o test_bgemv compare_sgemv_bgemv.c ../$(LIBNAME) $(EXTRALIB) $(FEXTRALIB) $(CEXTRALIB)

test_sbgemm$(EXE) : compare_sgemm_sbgemm.c test_helpers.h ../$(LIBNAME)
	$(CC) $(CLDFLAGS) -DIBFLOAT16 -o test_sbgemm compare_sgemm_sbgemm.c ../$(LIBNAME) $(EXTRALIB) $(FEXTRALIB) $(CEXTRALIB)

test_sbgemv$(EXE) : compare_sgemv_sbgemv.c ../$(LIBNAME)
	$(CC) $(CLDFLAGS) -DIBFLOAT16 -o test_sbgemv compare_sgemv_sbgemv.c ../$(LIBNAME) $(EXTRALIB) $(FEXTRALIB) $(CEXTRALIB)
endif

ifeq ($(BUILD_HFLOAT16),1)
test_shgemm$(EXE) : compare_sgemm_shgemm.c test_helpers.h ../$(LIBNAME)
	$(CC) $(CLDFLAGS) -DIHFLOAT16 -o test_shgemm compare_sgemm_shgemm.c ../$(LIBNAME) $(EXTRALIB) $(FEXTRALIB) $(CEXTRALIB)

test_shgemv$(EXE) : compare_sgemv_shgemv.c ../$(LIBNAME)
	$(CC) $(CLDFLAGS) -o test_shgemv compare_sgemv_shgemv.c ../$(LIBNAME) $(EXTRALIB) $(FEXTRALIB) $(CEXTRALIB)
endif


ifeq ($(BUILD_COMPLEX),1)
cblat3_3m$(EXE) : cblat3_3m.$(SUFFIX) $(XERBLA_TEST) ../$(LIBNAME)
	$(FC) $(FLDFLAGS) -o $@ cblat3_3m.$(SUFFIX) $(XERBLA_TEST) ../$(LIBNAME) $(EXTRALIB) $(CEXTRALIB)
endif

ifeq ($(BUILD_COMPLEX16),1)
zblat3_3m$(EXE) : zblat3_3m.$(SUFFIX) $(XERBLA_TEST) ../$(LIBNAME)
	$(FC) $(FLDFLAGS) -o $@ zblat3_3m.$(SUFFIX) $(XERBLA_TEST) ../$(LIBNAME) $(EXTRALIB) $(CEXTRALIB)
endif



clean:
	@rm -f *.$(SUFFIX) *.$(PSUFFIX) gmon.$(SUFFIX)ut *.SUMM *.cxml *.exe *.pdb *.dwf \
	sblat1$(EXE) dblat1$(EXE) cblat1$(EXE) zblat1$(EXE) \
	sblat2$(EXE) dblat2$(EXE) cblat2$(EXE) zblat2$(EXE) \
	test_bgemm$(EXE) test_bgemv$(EXE) test_sbgemm$(EXE) test_sbgemv$(EXE) test_shgemm$(EXE) test_shgemv$(EXE) sblat3$(EXE) dblat3$(EXE) cblat3$(EXE) zblat3$(EXE) \
	sblat1p$(EXE) dblat1p$(EXE) cblat1p$(EXE) zblat1p$(EXE) \
	sblat2p$(EXE) dblat2p$(EXE) cblat2p$(EXE) zblat2p$(EXE) \
	sblat3p$(EXE) dblat3p$(EXE) cblat3p$(EXE) zblat3p$(EXE) \
	zblat3_3m$(EXE) zblat3_3mp$(EXE) \
	cblat3_3m$(EXE) cblat3_3mp$(EXE) \
	*.stackdump *.dll

libs:

prof:

quick :
	$(MAKE) -C $(TOPDIR) libs

# include ../Makefile.tail
