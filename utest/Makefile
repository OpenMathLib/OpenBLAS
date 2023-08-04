UTEST_CHECK = 1
TOPDIR	= ..

override TARGET_ARCH=
override TARGET_MACH=

UTESTBIN=openblas_utest

.PHONY : all
.NOTPARALLEL : all run_test $(UTESTBIN)

include $(TOPDIR)/Makefile.system

OBJS=utest_main.o test_min.o test_amax.o test_ismin.o test_rotmg.o test_axpy.o test_dotu.o test_dsdot.o test_swap.o test_rot.o test_dnrm2.o
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
ifeq ($(OSNAME), $(filter $(OSNAME),Linux CYGWIN_NT))
ifneq ($(USE_OPENMP), 1)
OBJS += test_fork.o
endif
OBJS += test_post_fork.o
endif

ifeq ($(C_COMPILER), PGI)
OBJS = utest_main2.o
endif
ifeq ($(C_COMPILER), SUN)
OBJS = utest_main2.o
endif
ifeq ($(OSNAME), AIX)
OBJS = utest_main2.o
endif

all : run_test

$(UTESTBIN): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ ../$(LIBNAME) $(EXTRALIB) $(FEXTRALIB)

run_test: $(UTESTBIN)
ifneq ($(CROSS), 1)
	./$(UTESTBIN)
endif

clean:
	-rm -f *.o $(UTESTBIN)

libs:
