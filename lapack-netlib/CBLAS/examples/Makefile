include ../../make.inc

all: example1 example2

example1:
	$(CC) $(CFLAGS) -I../include -c cblas_example1.c
	$(LOADER) $(LOADOPTS) -o cblas_ex1 cblas_example1.o $(CBLASLIB) $(BLASLIB)

example2:
	$(CC) $(CFLAGS) -I../include -c cblas_example2.c
	$(LOADER) $(LOADOPTS) -o cblas_ex2 cblas_example2.o $(CBLASLIB) $(BLASLIB)

cleanall:
	rm -f *.o cblas_ex1 cblas_ex2
