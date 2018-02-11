UTEST_CHECK = 1
TOPDIR	= ..

UTESTBIN=openblas_utest

.PHONY : all
.NOTPARALLEL : all run_test $(UTESTBIN)

include $(TOPDIR)/Makefile.system

OBJS=utest_main.o test_amax.o test_rotmg.o
#test_rot.o test_swap.o test_axpy.o test_dotu.o test_dsdot.o test_fork.o

ifneq ($(NO_LAPACK), 1)
#OBJS += test_potrs.o
endif

ifndef USE_OPENMP
ifndef OS_WINDOWS
OBJS += test_fork.o
else
ifdef OS_CYGWIN_NT
OBJS += test_fork.o
endif
endif
endif

all : run_test

$(UTESTBIN): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ ../$(LIBNAME) $(EXTRALIB) $(FEXTRALIB)

run_test: $(UTESTBIN)
ifndef CROSS
	./$(UTESTBIN)
endif

clean:
	-rm -f *.o $(UTESTBIN)

libs:

