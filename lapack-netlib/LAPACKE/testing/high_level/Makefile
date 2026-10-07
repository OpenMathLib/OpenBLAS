#
# The Makefile compiles c wrappers and testers for CBLAS.
#

TOPSRCDIR = ../../..
include $(TOPSRCDIR)/make.inc

ifeq ($(BUILD_SINGLE),1)
SOURCESS = sgels.c  sgeqrf.c  sgetrf.c  sgetri.c  sgetrs.c  slange.c  smatgen.c  spotrf.c  spotrs.c
STARGET = slapacke_high_level$(EXE)
endif

ifeq ($(BUILD_DOUBLE),1)
SOURCESD = dgels.c  dgeqrf.c  dgetrf.c  dgetri.c  dgetrs.c  dlange.c  dmatgen.c  dpotrf.c  dpotrs.c
DTARGET = dlapacke_high_level$(EXE)
endif

ifeq ($(BUILD_COMPLEX),1)
SOURCESC = cgels.c  cgeqrf.c  cgetrf.c  cgetri.c  cgetrs.c  clange.c  cmatgen.c  cpotrf.c  cpotrs.c
CTARGET = clapacke_high_level$(EXE)
endif

ifeq ($(BUILD_COMPLEX16),1)
SOURCESZ = zgels.c  zgeqrf.c  zgetrf.c  zgetri.c  zgetrs.c  zlange.c  zmatgen.c  zpotrf.c  zpotrs.c
ZTARGET = zlapacke_high_level$(EXE)
endif

TESTSRCS = lapacke_test.c lapacke_test_alloc.c lapacke_test_nan_sweep.c lapacke_test_matgen_s.c slapacke_main.c
TESTSRCD = lapacke_test.c lapacke_test_alloc.c lapacke_test_nan_sweep.c lapacke_test_matgen_d.c dlapacke_main.c
TESTSRCC = lapacke_test.c lapacke_test_alloc.c lapacke_test_nan_sweep.c lapacke_test_matgen_c.c clapacke_main.c
TESTSRCZ = lapacke_test.c lapacke_test_alloc.c lapacke_test_nan_sweep.c lapacke_test_matgen_z.c zlapacke_main.c

all :: $(STARGET) $(DTARGET) $(CTARGET) $(ZTARGET)

slapacke_high_level$(EXE):
	$(CC) $(CFLAGS) -I../../include $(SOURCESS) $(TESTSRCS) -o $@ ../$(BLASLIB) $(EXTRALIB) $(CEXTRALIB)

dlapacke_high_level$(EXE):
	$(CC) $(CFLAGS) -I../../include $(SOURCESD) $(TESTSRCD) -o $@ ../$(BLASLIB) $(EXTRALIB) $(CEXTRALIB)

clapacke_high_level$(EXE):
	$(CC) $(CFLAGS) -I../../include $(SOURCESC) $(TESTSRCC) -o $@ ../$(BLASLIB) $(EXTRALIB) $(CEXTRALIB)

zlapacke_high_level$(EXE):
	$(CC) $(CFLAGS) -I../../include $(SOURCESZ) $(TESTSRCZ) -o $@ ../$(BLASLIB) $(EXTRALIB) $(CEXTRALIB)

clean ::
	rm -f *.o $(STARGET) $(DTARGET) $(CTARGET) $(ZTARGET)

