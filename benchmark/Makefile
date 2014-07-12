TOPDIR	= ..
include $(TOPDIR)/Makefile.system

# ACML standard
ACML=/opt/acml5.3.1/gfortran64_mp/lib
LIBACML = -fopenmp $(ACML)/libacml_mp.a -lgfortran -lm

# ACML custom
#ACML=/opt/pb/acml-5-3-1-gfortran-64bit/gfortran64_fma4_mp/lib
#LIBACML = -fopenmp $(ACML)/libacml_mp.a -lgfortran -lm

# Atlas Ubuntu 
#ATLAS=/usr/lib/atlas-base
#LIBATLAS = -fopenmp $(ATLAS)/liblapack_atlas.a  $(ATLAS)/libptcblas.a  $(ATLAS)/libptf77blas.a  $(ATLAS)/libatlas.a -lgfortran -lm

# Atlas RHEL and Fedora
ATLAS=/usr/lib64/atlas
LIBATLAS = -fopenmp $(ATLAS)/liblapack.a  $(ATLAS)/libptcblas.a  $(ATLAS)/libptf77blas.a  $(ATLAS)/libatlas.a -lgfortran -lm

# Intel standard
MKL=/opt/intel/mkl/lib/intel64
LIBMKL = -L$(MKL) -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm

# Intel custom
#MKL=/home/saar/intel_mkl
#LIBMKL = -L$(MKL) -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm



goto :: slinpack.goto dlinpack.goto clinpack.goto zlinpack.goto \
       scholesky.goto dcholesky.goto ccholesky.goto zcholesky.goto \
       sgemm.goto dgemm.goto cgemm.goto zgemm.goto \
       strmm.goto dtrmm.goto ctrmm.goto ztrmm.goto \
       strsm.goto dtrsm.goto ctrsm.goto ztrsm.goto \
       ssyrk.goto dsyrk.goto csyrk.goto zsyrk.goto \
       ssyr2k.goto dsyr2k.goto csyr2k.goto zsyr2k.goto \
       chemm.goto zhemm.goto \
       cherk.goto zherk.goto \
       cher2k.goto zher2k.goto \
       ssymm.goto dsymm.goto csymm.goto zsymm.goto

acml :: slinpack.acml dlinpack.acml clinpack.acml zlinpack.acml \
       scholesky.acml dcholesky.acml ccholesky.acml zcholesky.acml \
       sgemm.acml dgemm.acml cgemm.acml zgemm.acml \
       strmm.acml dtrmm.acml ctrmm.acml ztrmm.acml \
       strsm.acml dtrsm.acml ctrsm.acml ztrsm.acml \
       ssyrk.acml dsyrk.acml csyrk.acml zsyrk.acml \
       ssyr2k.acml dsyr2k.acml csyr2k.acml zsyr2k.acml \
       chemm.acml zhemm.acml \
       cherk.acml zherk.acml \
       cher2k.acml zher2k.acml \
       ssymm.acml dsymm.acml csymm.acml zsymm.acml

atlas :: slinpack.atlas dlinpack.atlas clinpack.atlas zlinpack.atlas \
       scholesky.atlas dcholesky.atlas ccholesky.atlas zcholesky.atlas \
       sgemm.atlas dgemm.atlas cgemm.atlas zgemm.atlas \
       strmm.atlas dtrmm.atlas ctrmm.atlas ztrmm.atlas \
       strsm.atlas dtrsm.atlas ctrsm.atlas ztrsm.atlas \
       ssyrk.atlas dsyrk.atlas csyrk.atlas zsyrk.atlas \
       ssyr2k.atlas dsyr2k.atlas csyr2k.atlas zsyr2k.atlas \
       chemm.atlas zhemm.atlas \
       cherk.atlas zherk.atlas \
       cher2k.atlas zher2k.atlas \
       ssymm.atlas dsymm.atlas csymm.atlas zsymm.atlas

mkl :: slinpack.mkl dlinpack.mkl clinpack.mkl zlinpack.mkl \
       scholesky.mkl dcholesky.mkl ccholesky.mkl zcholesky.mkl \
       sgemm.mkl dgemm.mkl cgemm.mkl zgemm.mkl \
       strmm.mkl dtrmm.mkl ctrmm.mkl ztrmm.mkl \
       strsm.mkl dtrsm.mkl ctrsm.mkl ztrsm.mkl \
       ssyrk.mkl dsyrk.mkl csyrk.mkl zsyrk.mkl \
       ssyr2k.mkl dsyr2k.mkl csyr2k.mkl zsyr2k.mkl \
       chemm.mkl zhemm.mkl \
       cherk.mkl zherk.mkl \
       cher2k.mkl zher2k.mkl \
       ssymm.mkl dsymm.mkl csymm.mkl zsymm.mkl

all :: goto atlas acml mkl

##################################### Slinpack ####################################################
slinpack.goto : slinpack.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

slinpack.acml : slinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

slinpack.atlas : slinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

slinpack.mkl : slinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Dlinpack ####################################################
dlinpack.goto : dlinpack.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

dlinpack.acml : dlinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dlinpack.atlas : dlinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dlinpack.mkl : dlinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Clinpack ####################################################

clinpack.goto : clinpack.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

clinpack.acml : clinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

clinpack.atlas : clinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

clinpack.mkl : clinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Zlinpack ####################################################

zlinpack.goto : zlinpack.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

zlinpack.acml : zlinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zlinpack.atlas : zlinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zlinpack.mkl : zlinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Scholesky ###################################################

scholesky.goto : scholesky.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

scholesky.acml : scholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

scholesky.atlas : scholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

scholesky.mkl : scholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Dcholesky ###################################################

dcholesky.goto : dcholesky.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

dcholesky.acml : dcholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dcholesky.atlas : dcholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dcholesky.mkl : dcholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Ccholesky ###################################################

ccholesky.goto : ccholesky.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

ccholesky.acml : ccholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

ccholesky.atlas : ccholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

ccholesky.mkl : ccholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zcholesky.goto : zcholesky.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

##################################### Zcholesky ###################################################

xcholesky.goto : xcholesky.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

zcholesky.acml : zcholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zcholesky.atlas : zcholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zcholesky.mkl : zcholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)


##################################### Sgemm ####################################################
sgemm.goto : sgemm.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

sgemm.acml : sgemm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

sgemm.atlas : sgemm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

sgemm.mkl : sgemm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Dgemm ####################################################
dgemm.goto : dgemm.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

dgemm.acml : dgemm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dgemm.atlas : dgemm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dgemm.mkl : dgemm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Cgemm ####################################################

cgemm.goto : cgemm.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

cgemm.acml : cgemm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

cgemm.atlas : cgemm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

cgemm.mkl : cgemm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Zgemm ####################################################

zgemm.goto : zgemm.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

zgemm.acml : zgemm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zgemm.atlas : zgemm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zgemm.mkl : zgemm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Ssymm ####################################################
ssymm.goto : ssymm.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

ssymm.acml : ssymm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

ssymm.atlas : ssymm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

ssymm.mkl : ssymm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Dsymm ####################################################
dsymm.goto : dsymm.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

dsymm.acml : dsymm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dsymm.atlas : dsymm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dsymm.mkl : dsymm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Csymm ####################################################

csymm.goto : csymm.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

csymm.acml : csymm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

csymm.atlas : csymm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

csymm.mkl : csymm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Zsymm ####################################################

zsymm.goto : zsymm.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

zsymm.acml : zsymm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zsymm.atlas : zsymm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zsymm.mkl : zsymm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Strmm ####################################################
strmm.goto : strmm.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

strmm.acml : strmm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

strmm.atlas : strmm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

strmm.mkl : strmm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Dtrmm ####################################################
dtrmm.goto : dtrmm.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

dtrmm.acml : dtrmm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dtrmm.atlas : dtrmm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dtrmm.mkl : dtrmm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Ctrmm ####################################################

ctrmm.goto : ctrmm.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

ctrmm.acml : ctrmm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

ctrmm.atlas : ctrmm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

ctrmm.mkl : ctrmm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Ztrmm ####################################################

ztrmm.goto : ztrmm.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

ztrmm.acml : ztrmm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

ztrmm.atlas : ztrmm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

ztrmm.mkl : ztrmm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)


##################################### Strsm ####################################################
strsm.goto : strsm.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

strsm.acml : strsm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

strsm.atlas : strsm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

strsm.mkl : strsm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Dtrsm ####################################################
dtrsm.goto : dtrsm.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

dtrsm.acml : dtrsm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dtrsm.atlas : dtrsm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dtrsm.mkl : dtrsm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Ctrsm ####################################################

ctrsm.goto : ctrsm.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

ctrsm.acml : ctrsm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

ctrsm.atlas : ctrsm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

ctrsm.mkl : ctrsm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Ztrsm ####################################################

ztrsm.goto : ztrsm.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

ztrsm.acml : ztrsm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

ztrsm.atlas : ztrsm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

ztrsm.mkl : ztrsm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Ssyrk ####################################################
ssyrk.goto : ssyrk.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

ssyrk.acml : ssyrk.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

ssyrk.atlas : ssyrk.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

ssyrk.mkl : ssyrk.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Dsyrk ####################################################
dsyrk.goto : dsyrk.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

dsyrk.acml : dsyrk.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dsyrk.atlas : dsyrk.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dsyrk.mkl : dsyrk.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Csyrk ####################################################

csyrk.goto : csyrk.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

csyrk.acml : csyrk.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

csyrk.atlas : csyrk.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

csyrk.mkl : csyrk.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Zsyrk ####################################################

zsyrk.goto : zsyrk.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

zsyrk.acml : zsyrk.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zsyrk.atlas : zsyrk.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zsyrk.mkl : zsyrk.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)


##################################### Ssyr2k ####################################################
ssyr2k.goto : ssyr2k.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

ssyr2k.acml : ssyr2k.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

ssyr2k.atlas : ssyr2k.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

ssyr2k.mkl : ssyr2k.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Dsyr2k ####################################################
dsyr2k.goto : dsyr2k.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

dsyr2k.acml : dsyr2k.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dsyr2k.atlas : dsyr2k.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dsyr2k.mkl : dsyr2k.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Csyr2k ####################################################

csyr2k.goto : csyr2k.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

csyr2k.acml : csyr2k.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

csyr2k.atlas : csyr2k.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

csyr2k.mkl : csyr2k.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Zsyr2k ####################################################

zsyr2k.goto : zsyr2k.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

zsyr2k.acml : zsyr2k.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zsyr2k.atlas : zsyr2k.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zsyr2k.mkl : zsyr2k.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Chemm ####################################################

chemm.goto : chemm.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

chemm.acml : chemm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

chemm.atlas : chemm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

chemm.mkl : chemm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Zhemm ####################################################

zhemm.goto : zhemm.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

zhemm.acml : zhemm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zhemm.atlas : zhemm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zhemm.mkl : zhemm.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Cherk ####################################################

cherk.goto : cherk.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

cherk.acml : cherk.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

cherk.atlas : cherk.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

cherk.mkl : cherk.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Zherk ####################################################

zherk.goto : zherk.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

zherk.acml : zherk.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zherk.atlas : zherk.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zherk.mkl : zherk.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Cher2k ####################################################

cher2k.goto : cher2k.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

cher2k.acml : cher2k.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

cher2k.atlas : cher2k.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

cher2k.mkl : cher2k.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

##################################### Zher2k ####################################################

zher2k.goto : zher2k.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

zher2k.acml : zher2k.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zher2k.atlas : zher2k.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBATLAS) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zher2k.mkl : zher2k.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

###################################################################################################

slinpack.$(SUFFIX) : linpack.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -UDOUBLE -o $(@F) $^

dlinpack.$(SUFFIX) : linpack.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -DDOUBLE -o $(@F) $^

clinpack.$(SUFFIX) : linpack.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -UDOUBLE -o $(@F) $^

zlinpack.$(SUFFIX) : linpack.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -DDOUBLE -o $(@F) $^

scholesky.$(SUFFIX) : cholesky.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -UDOUBLE -o $(@F) $^

dcholesky.$(SUFFIX) : cholesky.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -DDOUBLE -o $(@F) $^

ccholesky.$(SUFFIX) : cholesky.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -UDOUBLE -o $(@F) $^

zcholesky.$(SUFFIX) : cholesky.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -DDOUBLE -o $(@F) $^

sgemm.$(SUFFIX) : gemm.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -UDOUBLE -o $(@F) $^

dgemm.$(SUFFIX) : gemm.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -DDOUBLE -o $(@F) $^

cgemm.$(SUFFIX) : gemm.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -UDOUBLE -o $(@F) $^

zgemm.$(SUFFIX) : gemm.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -DDOUBLE -o $(@F) $^

ssymm.$(SUFFIX) : symm.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -UDOUBLE -o $(@F) $^

dsymm.$(SUFFIX) : symm.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -DDOUBLE -o $(@F) $^

csymm.$(SUFFIX) : symm.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -UDOUBLE -o $(@F) $^

zsymm.$(SUFFIX) : symm.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -DDOUBLE -o $(@F) $^

strmm.$(SUFFIX) : trmm.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -UDOUBLE -o $(@F) $^

dtrmm.$(SUFFIX) : trmm.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -DDOUBLE -o $(@F) $^

ctrmm.$(SUFFIX) : trmm.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -UDOUBLE -o $(@F) $^

ztrmm.$(SUFFIX) : trmm.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -DDOUBLE -o $(@F) $^

strsm.$(SUFFIX) : trsm.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -UDOUBLE -o $(@F) $^

dtrsm.$(SUFFIX) : trsm.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -DDOUBLE -o $(@F) $^

ctrsm.$(SUFFIX) : trsm.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -UDOUBLE -o $(@F) $^

ztrsm.$(SUFFIX) : trsm.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -DDOUBLE -o $(@F) $^

ssyrk.$(SUFFIX) : syrk.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -UDOUBLE -o $(@F) $^

dsyrk.$(SUFFIX) : syrk.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -DDOUBLE -o $(@F) $^

csyrk.$(SUFFIX) : syrk.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -UDOUBLE -o $(@F) $^

zsyrk.$(SUFFIX) : syrk.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -DDOUBLE -o $(@F) $^

ssyr2k.$(SUFFIX) : syr2k.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -UDOUBLE -o $(@F) $^

dsyr2k.$(SUFFIX) : syr2k.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -DDOUBLE -o $(@F) $^

csyr2k.$(SUFFIX) : syr2k.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -UDOUBLE -o $(@F) $^

zsyr2k.$(SUFFIX) : syr2k.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -DDOUBLE -o $(@F) $^

chemm.$(SUFFIX) : hemm.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -UDOUBLE -o $(@F) $^

zhemm.$(SUFFIX) : hemm.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -DDOUBLE -o $(@F) $^

cherk.$(SUFFIX) : herk.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -UDOUBLE -o $(@F) $^

zherk.$(SUFFIX) : herk.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -DDOUBLE -o $(@F) $^

cher2k.$(SUFFIX) : her2k.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -UDOUBLE -o $(@F) $^

zher2k.$(SUFFIX) : her2k.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -DDOUBLE -o $(@F) $^




clean ::
	@rm -f *.goto *.mkl *.acml *.atlas

include $(TOPDIR)/Makefile.tail

