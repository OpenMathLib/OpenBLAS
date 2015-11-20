include ../../make.inc

all: example1 example2

example1: 
	$(CC) -c $(CFLAGS) -I../include cblas_example1.c
	$(LOADER) -o cblas_ex1 cblas_example1.o $(CBLASLIB) $(BLASLIB)

example2: 
	$(CC) -c $(CFLAGS) -I../include cblas_example2.c
	$(LOADER) -o cblas_ex2 cblas_example2.o $(CBLASLIB) $(BLASLIB)
   
cleanall:
	rm -f *.o cblas_ex1 cblas_ex2
