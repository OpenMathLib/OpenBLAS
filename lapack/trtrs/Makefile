TOPDIR	= ../..
include ../../Makefile.system

SBLASOBJS = strtrs_UNU_single.$(SUFFIX) strtrs_UNN_single.$(SUFFIX) strtrs_UTU_single.$(SUFFIX) strtrs_UTN_single.$(SUFFIX) strtrs_LNU_single.$(SUFFIX) strtrs_LNN_single.$(SUFFIX) strtrs_LTU_single.$(SUFFIX) strtrs_LTN_single.$(SUFFIX)
DBLASOBJS = dtrtrs_UNU_single.$(SUFFIX) dtrtrs_UNN_single.$(SUFFIX) dtrtrs_UTU_single.$(SUFFIX) dtrtrs_UTN_single.$(SUFFIX) dtrtrs_LNU_single.$(SUFFIX) dtrtrs_LNN_single.$(SUFFIX) dtrtrs_LTU_single.$(SUFFIX) dtrtrs_LTN_single.$(SUFFIX)
QBLASOBJS = qtrtrs_UNU_single.$(SUFFIX) qtrtrs_UNN_single.$(SUFFIX) qtrtrs_UTU_single.$(SUFFIX) qtrtrs_UTN_single.$(SUFFIX) qtrtrs_LNU_single.$(SUFFIX) qtrtrs_LNN_single.$(SUFFIX) qtrtrs_LTU_single.$(SUFFIX) qtrtrs_LTN_single.$(SUFFIX)
CBLASOBJS = cgetrs_N_single.$(SUFFIX) cgetrs_T_single.$(SUFFIX) cgetrs_R_single.$(SUFFIX) cgetrs_C_single.$(SUFFIX)
ZBLASOBJS = zgetrs_N_single.$(SUFFIX) zgetrs_T_single.$(SUFFIX) zgetrs_R_single.$(SUFFIX) zgetrs_C_single.$(SUFFIX)
XBLASOBJS = xgetrs_N_single.$(SUFFIX) xgetrs_T_single.$(SUFFIX) xgetrs_R_single.$(SUFFIX) xgetrs_C_single.$(SUFFIX)

ifdef SMP
SBLASOBJS += sgetrs_N_parallel.$(SUFFIX) sgetrs_T_parallel.$(SUFFIX)
DBLASOBJS += dgetrs_N_parallel.$(SUFFIX) dgetrs_T_parallel.$(SUFFIX)
QBLASOBJS += qgetrs_N_parallel.$(SUFFIX) qgetrs_T_parallel.$(SUFFIX)
CBLASOBJS += cgetrs_N_parallel.$(SUFFIX) cgetrs_T_parallel.$(SUFFIX) cgetrs_R_parallel.$(SUFFIX) cgetrs_C_parallel.$(SUFFIX)
ZBLASOBJS += zgetrs_N_parallel.$(SUFFIX) zgetrs_T_parallel.$(SUFFIX) zgetrs_R_parallel.$(SUFFIX) zgetrs_C_parallel.$(SUFFIX)
XBLASOBJS += xgetrs_N_parallel.$(SUFFIX) xgetrs_T_parallel.$(SUFFIX) xgetrs_R_parallel.$(SUFFIX) xgetrs_C_parallel.$(SUFFIX)
endif

strtrs_UNU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -UDOUBLE -UUPLO -UTRANS -UDIAG $< -o $(@F)

strtrs_UNN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -UDOUBLE -UUPLO -UTRANS -DDIAG $< -o $(@F)

strtrs_UTU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -UDOUBLE -UUPLO -DTRANS -UDIAG $< -o $(@F)

strtrs_UTN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -UDOUBLE -UUPLO -DTRANS -DDIAG $< -o $(@F)

strtrs_LNU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -UDOUBLE -DUPLO -UTRANS -UDIAG $< -o $(@F)

strtrs_LNN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -UDOUBLE -DUPLO -UTRANS -DDIAG $< -o $(@F)

strtrs_LTU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -UDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

strtrs_LTN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -UDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

strtrs_UNU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -UDOUBLE -UUPLO -UTRANS -UDIAG $< -o $(@F)

strtrs_UNN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -UDOUBLE -UUPLO -UTRANS -DDIAG $< -o $(@F)

strtrs_UTU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -UDOUBLE -UUPLO -DTRANS -UDIAG $< -o $(@F)

strtrs_UTN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -UDOUBLE -UUPLO -DTRANS -DDIAG $< -o $(@F)

strtrs_LNU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -UDOUBLE -DUPLO -UTRANS -UDIAG $< -o $(@F)

strtrs_LNN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -UDOUBLE -DUPLO -UTRANS -DDIAG $< -o $(@F)

strtrs_LTU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -UDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

strtrs_LTN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -UDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

dtrtrs_UNU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DDOUBLE -UUPLO -UTRANS -UDIAG $< -o $(@F)

dtrtrs_UNN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DDOUBLE -UUPLO -UTRANS -DDIAG $< -o $(@F)

dtrtrs_UTU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DDOUBLE -UUPLO -DTRANS -UDIAG $< -o $(@F)

dtrtrs_UTN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DDOUBLE -UUPLO -DTRANS -DDIAG $< -o $(@F)

dtrtrs_LNU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DDOUBLE -DUPLO -UTRANS -UDIAG $< -o $(@F)

dtrtrs_LNN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DDOUBLE -DUPLO -UTRANS -DDIAG $< -o $(@F)

dtrtrs_LTU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

dtrtrs_LTN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

dtrtrs_UNU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DDOUBLE -UUPLO -UTRANS -UDIAG $< -o $(@F)

dtrtrs_UNN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DDOUBLE -UUPLO -UTRANS -DDIAG $< -o $(@F)

dtrtrs_UTU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DDOUBLE -UUPLO -DTRANS -UDIAG $< -o $(@F)

dtrtrs_UTN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DDOUBLE -UUPLO -DTRANS -DDIAG $< -o $(@F)

dtrtrs_LNU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DDOUBLE -DUPLO -UTRANS -UDIAG $< -o $(@F)

dtrtrs_LNN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DDOUBLE -DUPLO -UTRANS -DDIAG $< -o $(@F)

dtrtrs_LTU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

dtrtrs_LTN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

qtrtrs_UNU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DXDOUBLE -UUPLO -UTRANS -UDIAG $< -o $(@F)

qtrtrs_UNN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DXDOUBLE -UUPLO -UTRANS -DDIAG $< -o $(@F)

qtrtrs_UTU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DXDOUBLE -UUPLO -DTRANS -UDIAG $< -o $(@F)

qtrtrs_UTN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DXDOUBLE -UUPLO -DTRANS -DDIAG $< -o $(@F)

qtrtrs_LNU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DXDOUBLE -DUPLO -UTRANS -UDIAG $< -o $(@F)

qtrtrs_LNN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DXDOUBLE -DUPLO -UTRANS -DDIAG $< -o $(@F)

qtrtrs_LTU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DXDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

qtrtrs_LTN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DXDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

qtrtrs_UNU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DXDOUBLE -UUPLO -UTRANS -UDIAG $< -o $(@F)

qtrtrs_UNN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DXDOUBLE -UUPLO -UTRANS -DDIAG $< -o $(@F)

qtrtrs_UTU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DXDOUBLE -UUPLO -DTRANS -UDIAG $< -o $(@F)

qtrtrs_UTN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DXDOUBLE -UUPLO -DTRANS -DDIAG $< -o $(@F)

qtrtrs_LNU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DXDOUBLE -DUPLO -UTRANS -UDIAG $< -o $(@F)

qtrtrs_LNN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DXDOUBLE -DUPLO -UTRANS -DDIAG $< -o $(@F)

qtrtrs_LTU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DXDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

qtrtrs_LTN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DXDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

ctrtrs_UNU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -UDOUBLE -UUPLO -UTRANS -UDIAG $< -o $(@F)

ctrtrs_UNN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -UDOUBLE -UUPLO -UTRANS -DDIAG $< -o $(@F)

ctrtrs_UTU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -UDOUBLE -UUPLO -DTRANS -UDIAG $< -o $(@F)

ctrtrs_UTN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -UDOUBLE -UUPLO -DTRANS -DDIAG $< -o $(@F)

ctrtrs_LNU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -UDOUBLE -DUPLO -UTRANS -UDIAG $< -o $(@F)

ctrtrs_LNN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -UDOUBLE -DUPLO -UTRANS -DDIAG $< -o $(@F)

ctrtrs_LTU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -UDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

ctrtrs_LTN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -UDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

ctrtrs_UNU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -UDOUBLE -UUPLO -UTRANS -UDIAG $< -o $(@F)

ctrtrs_UNN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -UDOUBLE -UUPLO -UTRANS -DDIAG $< -o $(@F)

ctrtrs_UTU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -UDOUBLE -UUPLO -DTRANS -UDIAG $< -o $(@F)

ctrtrs_UTN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -UDOUBLE -UUPLO -DTRANS -DDIAG $< -o $(@F)

ctrtrs_LNU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -UDOUBLE -DUPLO -UTRANS -UDIAG $< -o $(@F)

ctrtrs_LNN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -UDOUBLE -DUPLO -UTRANS -DDIAG $< -o $(@F)

ctrtrs_LTU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -UDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

ctrtrs_LTN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -UDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

ztrtrs_UNU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DDOUBLE -UUPLO -UTRANS -UDIAG $< -o $(@F)

ztrtrs_UNN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DDOUBLE -UUPLO -UTRANS -DDIAG $< -o $(@F)

ztrtrs_UTU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DDOUBLE -UUPLO -DTRANS -UDIAG $< -o $(@F)

ztrtrs_UTN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DDOUBLE -UUPLO -DTRANS -DDIAG $< -o $(@F)

ztrtrs_LNU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DDOUBLE -DUPLO -UTRANS -UDIAG $< -o $(@F)

ztrtrs_LNN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DDOUBLE -DUPLO -UTRANS -DDIAG $< -o $(@F)

ztrtrs_LTU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

ztrtrs_LTN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

ztrtrs_UNU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DDOUBLE -UUPLO -UTRANS -UDIAG $< -o $(@F)

ztrtrs_UNN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DDOUBLE -UUPLO -UTRANS -DDIAG $< -o $(@F)

ztrtrs_UTU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DDOUBLE -UUPLO -DTRANS -UDIAG $< -o $(@F)

ztrtrs_UTN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DDOUBLE -UUPLO -DTRANS -DDIAG $< -o $(@F)

ztrtrs_LNU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DDOUBLE -DUPLO -UTRANS -UDIAG $< -o $(@F)

ztrtrs_LNN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DDOUBLE -DUPLO -UTRANS -DDIAG $< -o $(@F)

ztrtrs_LTU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

ztrtrs_LTN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

qtrtrs_UNU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DXDOUBLE -UUPLO -UTRANS -UDIAG $< -o $(@F)

qtrtrs_UNN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DXDOUBLE -UUPLO -UTRANS -DDIAG $< -o $(@F)

qtrtrs_UTU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DXDOUBLE -UUPLO -DTRANS -UDIAG $< -o $(@F)

qtrtrs_UTN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DXDOUBLE -UUPLO -DTRANS -DDIAG $< -o $(@F)

qtrtrs_LNU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DXDOUBLE -DUPLO -UTRANS -UDIAG $< -o $(@F)

qtrtrs_LNN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DXDOUBLE -DUPLO -UTRANS -DDIAG $< -o $(@F)

qtrtrs_LTU_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DXDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

qtrtrs_LTN_single.$(SUFFIX) : trtrs_single.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DXDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

qtrtrs_UNU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DXDOUBLE -UUPLO -UTRANS -UDIAG $< -o $(@F)

qtrtrs_UNN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DXDOUBLE -UUPLO -UTRANS -DDIAG $< -o $(@F)

qtrtrs_UTU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DXDOUBLE -UUPLO -DTRANS -UDIAG $< -o $(@F)

qtrtrs_UTN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DXDOUBLE -UUPLO -DTRANS -DDIAG $< -o $(@F)

qtrtrs_LNU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DXDOUBLE -DUPLO -UTRANS -UDIAG $< -o $(@F)

qtrtrs_LNN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DXDOUBLE -DUPLO -UTRANS -DDIAG $< -o $(@F)

qtrtrs_LTU_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DXDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

qtrtrs_LTN_parallel.$(SUFFIX) : trtrs_parallel.c
	$(CC) -c $(CFLAGS) -DCOMPLEX -DXDOUBLE -DUPLO -DTRANS -DDIAG $< -o $(@F)

include ../../Makefile.tail
