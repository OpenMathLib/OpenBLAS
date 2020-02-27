UTEST_CHECK = 1
TOPDIR	= ..

override TARGET_ARCH=
override TARGET_MACH=

UTESTBIN=openblas_utest

.PHONY : all
.NOTPARALLEL : all run_test $(UTESTBIN)

include $(TOPDIR)/Makefile.system

OBJS=utest_main.o test_amax.o test_ismin.o test_rotmg.o test_axpy.o test_dotu.o test_dsdot.o test_swap.o test_rot.o
#test_rot.o test_swap.o test_axpy.o test_dotu.o test_dsdot.o test_fork.o

ifneq ($(NO_LAPACK), 1)
OBJS += test_potrs.o
ifneq ($(NO_CBLAS), 1)
ifneq ($(NO_LAPACKE), 1)
OBJS += test_kernel_regress.o
endif
endif
endif

#this does not work with OpenMP nor with native Windows or Android threads
# FIXME TBD if this works on OSX, SunOS, POWER and zarch
ifndef USE_OPENMP
ifeq ($(OSNAME), $(filter $(OSNAME),Linux CYGWIN_NT))
OBJS += test_fork.o
endif
endif

ifeq ($(C_COMPILER), PGI)
OBJS = utest_main2.o
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
