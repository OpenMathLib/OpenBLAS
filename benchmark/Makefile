TOPDIR	= ..
include $(TOPDIR)/Makefile.system

CULA_INC = -I/usr/local/cula/include
CULA_LIB = -L/usr/local/cula/lib64 -Wl,-rpath,/usr/local/cula/lib64 -lcula_fortran -lcula -lcublas

all :: dlinpack.goto dlinpack.mkl dlinpack.acml dcholesky.goto dcholesky.mkl dcholesky.acml
	./dlinpack.goto  4000 4000 1
	-./dlinpack.mkl  4000 4000 1
	-./dlinpack.acml 4000 4000 1
	./dcholesky.goto  4000 4000 1
	-./dcholesky.mkl  4000 4000 1
	-./dcholesky.acml 4000 4000 1

slinpack.goto : slinpack.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

dlinpack.goto : dlinpack.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

qlinpack.goto : qlinpack.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

clinpack.goto : clinpack.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

zlinpack.goto : zlinpack.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

xlinpack.goto : xlinpack.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

scholesky.goto : scholesky.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

dcholesky.goto : dcholesky.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

qcholesky.goto : qcholesky.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

ccholesky.goto : ccholesky.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

zcholesky.goto : zcholesky.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

xcholesky.goto : xcholesky.$(SUFFIX) ../$(LIBNAME)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CEXTRALIB) $(EXTRALIB) -lm

slinpack.mkl : slinpack.$(SUFFIX)
	-$(CC) -static $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dlinpack.mkl : dlinpack.$(SUFFIX)
	-$(CC) -static $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

clinpack.mkl : clinpack.$(SUFFIX)
	-$(CC) -static $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zlinpack.mkl : zlinpack.$(SUFFIX)
	-$(CC) -static $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

scholesky.mkl : scholesky.$(SUFFIX)
	-$(CC) -static $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dcholesky.mkl : dcholesky.$(SUFFIX)
	-$(CC) -static $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

ccholesky.mkl : ccholesky.$(SUFFIX)
	-$(CC) -static $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zcholesky.mkl : zcholesky.$(SUFFIX)
	-$(CC) -static $(CFLAGS) -o $(@F) $^ $(LIBMKL) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

slinpack.acml : slinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dlinpack.acml : dlinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

clinpack.acml : clinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zlinpack.acml : zlinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

scholesky.acml : scholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dcholesky.acml : dcholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

ccholesky.acml : ccholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zcholesky.acml : zcholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBACML) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

slinpack.flame : slinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBFLAME) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dlinpack.flame : dlinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBFLAME) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

clinpack.flame : clinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBFLAME) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zlinpack.flame : zlinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBFLAME) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

scholesky.flame : scholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBFLAME) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dcholesky.flame : dcholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBFLAME) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

ccholesky.flame : ccholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBFLAME) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zcholesky.flame : zcholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBFLAME) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

slinpack.sun : slinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBSUNPERF) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dlinpack.sun : dlinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBSUNPERF) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

clinpack.sun : clinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBSUNPERF) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zlinpack.sun : zlinpack.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBSUNPERF) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

scholesky.sun : scholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBSUNPERF) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

dcholesky.sun : dcholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBSUNPERF) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

ccholesky.sun : ccholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBSUNPERF) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

zcholesky.sun : zcholesky.$(SUFFIX)
	-$(CC) $(CFLAGS) -o $(@F) $^ $(LIBSUNPERF) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

slinpack.cula : slinpack.$(SUFFIX) cula_wrapper.$(SUFFIX)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CULA_LIB) ../$(LIBNAME) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

clinpack.cula : clinpack.$(SUFFIX) cula_wrapper.$(SUFFIX)
	$(CC) $(CFLAGS) -o $(@F) $^ $(CULA_LIB) ../$(LIBNAME) $(CEXTRALIB) $(EXTRALIB) $(FEXTRALIB)

cula_wrapper.$(SUFFIX) : cula_wrapper.c
	$(CC) $(CFLAGS) -c $(CULA_INC) -o $(@F) $^

slinpack.$(SUFFIX) : linpack.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -UDOUBLE -o $(@F) $^

dlinpack.$(SUFFIX) : linpack.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -DDOUBLE -o $(@F) $^

qlinpack.$(SUFFIX) : linpack.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -DXDOUBLE -o $(@F) $^

clinpack.$(SUFFIX) : linpack.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -UDOUBLE -o $(@F) $^

zlinpack.$(SUFFIX) : linpack.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -DDOUBLE -o $(@F) $^

xlinpack.$(SUFFIX) : linpack.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -DXDOUBLE -o $(@F) $^

scholesky.$(SUFFIX) : cholesky.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -UDOUBLE -o $(@F) $^

dcholesky.$(SUFFIX) : cholesky.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -DDOUBLE -o $(@F) $^

qcholesky.$(SUFFIX) : cholesky.c
	$(CC) $(CFLAGS) -c -UCOMPLEX -DXDOUBLE -o $(@F) $^

ccholesky.$(SUFFIX) : cholesky.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -UDOUBLE -o $(@F) $^

zcholesky.$(SUFFIX) : cholesky.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -DDOUBLE -o $(@F) $^

xcholesky.$(SUFFIX) : cholesky.c
	$(CC) $(CFLAGS) -c -DCOMPLEX -DXDOUBLE -o $(@F) $^

clean ::
	@rm -f *.goto *.mkl *.acml *.sun *.cula

include $(TOPDIR)/Makefile.tail
