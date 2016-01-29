UTEST_CHECK = 1
TOPDIR	= ..

TARGET=openblas_utest

.PHONY : all
.NOTPARALLEL : all run_test $(TARGET)

include $(TOPDIR)/Makefile.system

OBJS=utest_main.o test_amax.o
#test_rot.o test_swap.o test_axpy.o test_dotu.o test_rotmg.o test_dsdot.o test_fork.o

all : run_test

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ ../$(LIBNAME) $(EXTRALIB)

run_test: $(TARGET)
	./$(TARGET)

clean:
	-rm -f *.o $(TARGET)

libs:

