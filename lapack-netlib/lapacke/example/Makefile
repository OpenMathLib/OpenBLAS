include ../../make.inc

all: xexample_DGESV_rowmajor xexample_ZGESV_rowmajor

xexample_DGESV_rowmajor: example_DGESV_rowmajor.o ../../$(LAPACKLIB) ../../$(LAPACKELIB)
	$(LOADER) $(LOADOPTS)  example_DGESV_rowmajor.o \
         ../../$(LAPACKELIB) $(CEXTRALIB) -o $@
	./$@

xexample_ZGESV_rowmajor: example_ZGESV_rowmajor.o ../../$(LAPACKLIB) ../../$(LAPACKELIB)
	$(LOADER) $(LOADOPTS)  example_ZGESV_rowmajor.o \
         ../../$(LAPACKELIB) $(CEXTRALIB) -o $@
	./$@

.c.o:
	$(CC) -c $(CFLAGS) -I ../include -o $@ $<

clean:
	rm -f *.o xexample_DGESV_rowmajor xexample_ZGESV_rowmajor