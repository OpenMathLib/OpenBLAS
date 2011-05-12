UTEST_CHECK = 1
TOPDIR	= ..
include $(TOPDIR)/Makefile.system

TARGET=openblas_utest
CUNIT_LIB=/usr/local/lib/libcunit.a

OBJS=main.o test_rot.o test_swap.o test_axpy.o test_dotu.o test_rotmg.o test_dsdot.o

all : run_test

$(TARGET): $(OBJS)
	$(FC) -o $@ $^ ../$(LIBNAME) $(CUNIT_LIB) $(EXTRALIB)

run_test: $(TARGET)
	./$(TARGET)

clean:
	rm -f *.o $(TARGET)

libs:

