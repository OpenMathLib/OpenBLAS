TOPDIR	= ../..
include ../../Makefile.system

SBLASOBJS = slaed3_single.$(SUFFIX)
DBLASOBJS = dlaed3_single.$(SUFFIX)

ifdef SMP
SBLASOBJS += slaed3_parallel.$(SUFFIX)
DBLASOBJS += dlaed3_parallel.$(SUFFIX)
endif

ifeq "$(or $(BUILD_SINGLE),$(BUILD_DOUBLE))" ""
SBLASOBJS=
endif
ifneq ($(BUILD_DOUBLE),1)
DBLASOBJS=
endif

slaed3_single.$(SUFFIX) : laed3_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -UDOUBLE $< -o $(@F)
dlaed3_single.$(SUFFIX) : laed3_single.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DDOUBLE $< -o $(@F)
slaed3_parallel.$(SUFFIX) : laed3_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -UDOUBLE $< -o $(@F)
dlaed3_parallel.$(SUFFIX) : laed3_parallel.c
	$(CC) -c $(CFLAGS) -UCOMPLEX -DDOUBLE $< -o $(@F)
slaed3_single.$(PSUFFIX) : laed3_single.c
	$(CC) -c $(PFLAGS) -UCOMPLEX -UDOUBLE $< -o $(@F)
dlaed3_single.$(PSUFFIX) : laed3_single.c
	$(CC) -c $(PFLAGS) -UCOMPLEX -DDOUBLE $< -o $(@F)
slaed3_parallel.$(PSUFFIX) : laed3_parallel.c
	$(CC) -c $(PFLAGS) -UCOMPLEX -UDOUBLE $< -o $(@F)
dlaed3_parallel.$(PSUFFIX) : laed3_parallel.c
	$(CC) -c $(PFLAGS) -UCOMPLEX -DDOUBLE $< -o $(@F)

include ../../Makefile.tail
