# This Makefile compiles the CBLAS routines
#
include ../../make.inc

#
# Erase all object and archive files
#
all: cblaslib

clean:
	rm -f *.o a.out core 

# Error handling routines for level 2 & 3

errhand = cblas_globals.o cblas_xerbla.o xerbla.o

# Object files of all routines
alev = $(alev1) $(alev2) $(alev3) $(errhand) 
#
#
# CBLAS routines
#
# Level 1
#
# 

#
# All object files for single real precision
#
slev1 = cblas_srotg.o cblas_srotmg.o cblas_srot.o cblas_srotm.o  \
         cblas_sswap.o cblas_sscal.o  cblas_scopy.o cblas_saxpy.o \
         cblas_sdot.o  cblas_sdsdot.o cblas_snrm2.o cblas_sasum.o \
         cblas_isamax.o   sdotsub.o   sdsdotsub.o snrm2sub.o sasumsub.o \
         isamaxsub.o 
#
# All object files for double real precision
#
dlev1 = cblas_drotg.o cblas_drotmg.o cblas_drot.o cblas_drotm.o  \
         cblas_dswap.o cblas_dscal.o  cblas_dcopy.o cblas_daxpy.o \
         cblas_ddot.o  cblas_dsdot.o cblas_dnrm2.o cblas_dasum.o  \
         cblas_idamax.o   ddotsub.o   dsdotsub.o       dnrm2sub.o \
         dasumsub.o idamaxsub.o

#
# All object files for single complex precision
#
clev1 = cblas_cswap.o cblas_cscal.o cblas_csscal.o cblas_ccopy.o \
         cblas_caxpy.o cblas_cdotu_sub.o cblas_cdotc_sub.o \
         cblas_icamax.o cdotcsub.o cdotusub.o icamaxsub.o 

#
# All object files for double complex precision
#
zlev1 = cblas_zswap.o cblas_zscal.o cblas_zdscal.o cblas_zcopy.o \
         cblas_zaxpy.o cblas_zdotu_sub.o cblas_zdotc_sub.o cblas_dznrm2.o \
         cblas_dzasum.o cblas_izamax.o zdotcsub.o zdotusub.o      \
         dzasumsub.o  dznrm2sub.o  izamaxsub.o 

#
# Common files for single / complex precision
#
sclev1 = cblas_scasum.o scasumsub.o cblas_scnrm2.o scnrm2sub.o

#
# All object files
#
alev1 = $(slev1) $(dlev1) $(clev1) $(zlev1) $(sclev1)


#
# Make an archive file
#

# Single real precision
slib1: $(slev1) $(sclev1)
	$(ARCH) $(ARCHFLAGS) $(CBLASLIB) $(slev1) $(sclev1)
	$(RANLIB) $(CBLASLIB)

# Double real precision
dlib1: $(dlev1) 
	$(ARCH) $(ARCHFLAGS) $(CBLASLIB) $(dlev1)
	$(RANLIB) $(CBLASLIB)

# Single complex precision
clib1: $(clev1) $(sclev1)
	$(ARCH) $(ARCHFLAGS) $(CBLASLIB) $(clev1) $(sclev1)
	$(RANLIB) $(CBLASLIB)

# Double complex precision
zlib1: $(zlev1)
	$(ARCH) $(ARCHFLAGS) $(CBLASLIB) $(zlev1)
	$(RANLIB) $(CBLASLIB)

# All precisions
all1:  $(alev1)
	$(ARCH) $(ARCHFLAGS)  $(CBLASLIB) $(alev1)
	$(RANLIB) $(CBLASLIB)

#
#
# CBLAS routines
#
# Level 2
#
# 

#
# All object files for single real precision
#
slev2 = cblas_sgemv.o cblas_sgbmv.o cblas_sger.o cblas_ssbmv.o  cblas_sspmv.o \
        cblas_sspr.o  cblas_sspr2.o cblas_ssymv.o cblas_ssyr.o  cblas_ssyr2.o \
        cblas_stbmv.o cblas_stbsv.o cblas_stpmv.o cblas_stpsv.o  cblas_strmv.o \
        cblas_strsv.o  
 
#
# All object files for double real precision
#
dlev2 = cblas_dgemv.o cblas_dgbmv.o cblas_dger.o cblas_dsbmv.o  cblas_dspmv.o \
        cblas_dspr.o  cblas_dspr2.o cblas_dsymv.o cblas_dsyr.o  cblas_dsyr2.o \
        cblas_dtbmv.o cblas_dtbsv.o cblas_dtpmv.o cblas_dtpsv.o  cblas_dtrmv.o \
        cblas_dtrsv.o  

#
# All object files for single complex precision
#
clev2 = cblas_cgemv.o cblas_cgbmv.o cblas_chemv.o cblas_chbmv.o cblas_chpmv.o \
        cblas_ctrmv.o cblas_ctbmv.o cblas_ctpmv.o cblas_ctrsv.o cblas_ctbsv.o \
        cblas_ctpsv.o cblas_cgeru.o cblas_cgerc.o cblas_cher.o  cblas_cher2.o \
        cblas_chpr.o  cblas_chpr2.o 

#
# All object files for double complex precision
#
zlev2 = cblas_zgemv.o cblas_zgbmv.o cblas_zhemv.o cblas_zhbmv.o cblas_zhpmv.o \
        cblas_ztrmv.o cblas_ztbmv.o cblas_ztpmv.o cblas_ztrsv.o cblas_ztbsv.o \
        cblas_ztpsv.o cblas_zgeru.o cblas_zgerc.o cblas_zher.o  cblas_zher2.o \
        cblas_zhpr.o  cblas_zhpr2.o 
#
# All object files
#
alev2 = $(slev2) $(dlev2) $(clev2) $(zlev2)

#
# Make an archive file
#

# Single real precision
slib2: $(slev2) $(errhand)
	$(ARCH) $(ARCHFLAGS) $(CBLASLIB) $(slev2) $(errhand)
	$(RANLIB) $(CBLASLIB)

# Double real  precision
dlib2: $(dlev2) $(errhand)
	$(ARCH) $(ARCHFLAGS) $(CBLASLIB) $(dlev2) $(errhand)
	$(RANLIB) $(CBLASLIB)

# Single complex precision
clib2: $(clev2) $(errhand)
	$(ARCH) $(ARCHFLAGS) $(CBLASLIB) $(clev2) $(errhand)
	$(RANLIB) $(CBLASLIB)

# Double complex precision
zlib2: $(zlev2) $(errhand)
	$(ARCH) $(ARCHFLAGS) $(CBLASLIB) $(zlev2) $(errhand)
	$(RANLIB) $(CBLASLIB)

# All precisions
all2:  $(alev2) $(errhand)
	$(ARCH) $(ARCHFLAGS)  $(CBLASLIB) $(alev2) $(errhand)
	$(RANLIB) $(CBLASLIB)
#
#
# CBLAS routines
#
# Level 3
#
# 

#
# All object files for single real precision
#
slev3 = cblas_sgemm.o cblas_ssymm.o cblas_ssyrk.o cblas_ssyr2k.o cblas_strmm.o\
        cblas_strsm.o
 
#
# All object files for double real precision
#
dlev3 = cblas_dgemm.o cblas_dsymm.o cblas_dsyrk.o cblas_dsyr2k.o cblas_dtrmm.o\
        cblas_dtrsm.o  

#
# All object files for single complex precision
#
clev3 = cblas_cgemm.o cblas_csymm.o cblas_chemm.o cblas_cherk.o\
        cblas_cher2k.o cblas_ctrmm.o cblas_ctrsm.o  cblas_csyrk.o\
        cblas_csyr2k.o
#
# All object files for double complex precision
#
zlev3 = cblas_zgemm.o cblas_zsymm.o cblas_zhemm.o cblas_zherk.o\
        cblas_zher2k.o cblas_ztrmm.o cblas_ztrsm.o  cblas_zsyrk.o\
        cblas_zsyr2k.o
#
# All object files
#
alev3 = $(slev3) $(dlev3) $(clev3) $(zlev3)

#
# Make an archive file
#

# Single real precision
slib3: $(slev3) $(errhand)
	$(ARCH) $(ARCHFLAGS) $(CBLASLIB) $(slev3) $(errhand)
	$(RANLIB) $(CBLASLIB)

# Double real precision
dlib3: $(dlev3) $(errhand)
	$(ARCH) $(ARCHFLAGS) $(CBLASLIB) $(dlev3) $(errhand)
	$(RANLIB) $(CBLASLIB)

# Single complex precision
clib3: $(clev3) $(errhand)
	$(ARCH) $(ARCHFLAGS) $(CBLASLIB) $(clev3) $(errhand)
	$(RANLIB) $(CBLASLIB)

# Single complex precision
zlib3: $(zlev3) $(errhand)
	$(ARCH) $(ARCHFLAGS) $(CBLASLIB) $(zlev3) $(errhand)
	$(RANLIB) $(CBLASLIB)

# All precisions
all3:  $(alev3) $(errhand)
	$(ARCH) $(ARCHFLAGS) $(CBLASLIB) $(alev3) 
	$(RANLIB) $(CBLASLIB)

# All levels and precisions
cblaslib: $(alev)  
	$(ARCH) $(ARCHFLAGS)  $(CBLASLIB) $(alev) 
	$(RANLIB) $(CBLASLIB)

FRC:
	@FRC=$(FRC)

.c.o:
	$(CC) -c $(CFLAGS) -I ../include -o $@ $<
	
.f.o: 
	$(FORTRAN) $(OPTS) -c $< -o $@
