ifdef TARGET_CORE
TARGET = $(TARGET_CORE)
endif

TOPDIR	= ..
include $(TOPDIR)/Makefile.system

AVX2OPT = 
ifeq ($(C_COMPILER), GCC)
# AVX2 support was added in 4.7.0
GCCVERSIONCHECK := $(GCCVERSIONGT4)$(GCCVERSIONGTEQ4)$(GCCMINORVERSIONGTEQ7)
ifeq ($(GCCVERSIONCHECK), $(filter $(GCCVERSIONCHECK), 011 110 111))
   AVX2OPT = -mavx2
  endif
endif
ifeq ($(C_COMPILER), CLANG)
# Any clang posing as gcc 4.2 should be new enough (3.4 or later)
  GCCVERSIONCHECK := $(GCCVERSIONGT4)$(GCCVERSIONGTEQ4)$(GCCMINORVERSIONGTEQ2)
  ifeq ($(GCCVERSIONCHECK), $(filter $(GCCVERSIONCHECK), 011 110 111))
   AVX2OPT = -mavx2 -mfma
  endif
endif
ifdef NO_AVX2
 AVX2OPT=
endif


ifdef TARGET_CORE
ifeq ($(TARGET_CORE), ARMV9SME)
 override CFLAGS += -DBUILD_KERNEL -DTABLE_NAME=gotoblas_$(TARGET_CORE) -march=armv9-a+sve2+sme+sme-f64f64
 ifdef OS_WINDOWS
 ifeq ($(C_COMPILER), CLANG)
   override CFLAGS += --aarch64-stack-hazard-size=0
 endif
 endif
endif
ifeq ($(TARGET_CORE), VORTEXM4)
  ifeq ($(C_COMPILER), GCC)
    override CFLAGS += -DBUILD_KERNEL -DTABLE_NAME=gotoblas_$(TARGET_CORE) -UHAVE_SME -march=armv8.4-a
  else
    override CFLAGS += -DBUILD_KERNEL -DTABLE_NAME=gotoblas_$(TARGET_CORE) -march=armv8.4-a+sme+sme-f64f64
#   ifneq ($(APPLECLANG),1)
#     override LDFLAGS += -lclang_rt_builtins-aarch64
#   endif
    ifdef OS_WINDOWS
      ifeq ($(C_COMPILER), CLANG)
        override CFLAGS += --aarch64-stack-hazard-size=0
      endif
    endif
  endif
endif
ifeq ($(TARGET_CORE), SAPPHIRERAPIDS)
 override CFLAGS += -DBUILD_KERNEL -DTABLE_NAME=gotoblas_$(TARGET_CORE)
 ifeq (1, $(filter 1,$(GCCVERSIONGTEQ11) $(CLANGVERSIONGTEQ12)))
  override CFLAGS += -march=sapphirerapids
 else 
  override CFLAGS += -march=skylake-avx512 -mavx512f
 endif 
 ifeq ($(OSNAME), CYGWIN_NT)
  override CFLAGS += -fno-asynchronous-unwind-tables
 endif
 ifeq ($(OSNAME), WINNT)
  ifeq ($(C_COMPILER), GCC)
   override CFLAGS += -fno-asynchronous-unwind-tables
  endif
 endif
else ifeq ($(TARGET_CORE), COOPERLAKE)
 override CFLAGS += -DBUILD_KERNEL -DTABLE_NAME=gotoblas_$(TARGET_CORE)
 ifeq (1, $(filter 1,$(GCCVERSIONGTEQ10) $(CLANGVERSIONGTEQ9)))
  override CFLAGS += -march=cooperlake
 else 
  override CFLAGS += -march=skylake-avx512 -mavx512f
 endif 
 ifeq ($(OSNAME), CYGWIN_NT)
  override CFLAGS += -fno-asynchronous-unwind-tables
 endif
 ifeq ($(OSNAME), WINNT)
  ifeq ($(C_COMPILER), GCC)
   override CFLAGS += -fno-asynchronous-unwind-tables
  endif
 endif
else ifeq ($(TARGET_CORE), SKYLAKEX)
 override CFLAGS += -DBUILD_KERNEL -DTABLE_NAME=gotoblas_$(TARGET_CORE) -march=skylake-avx512 -mavx512f
 ifeq ($(OSNAME), CYGWIN_NT)
  override CFLAGS += -fno-asynchronous-unwind-tables
 endif
 ifeq ($(OSNAME), WINNT)
  ifeq ($(C_COMPILER), GCC)
   override CFLAGS += -fno-asynchronous-unwind-tables
  endif
 endif
else ifeq ($(TARGET_CORE), HASWELL)
 override CFLAGS += -DBUILD_KERNEL -DTABLE_NAME=gotoblas_$(TARGET_CORE) $(AVX2OPT)
else ifeq ($(TARGET_CORE), ZEN)
 override CFLAGS += -DBUILD_KERNEL -DTABLE_NAME=gotoblas_$(TARGET_CORE) $(AVX2OPT)
else ifeq ($(TARGET_CORE), LOONGSON3R4)
 override CFLAGS += -DBUILD_KERNEL -DTABLE_NAME=gotoblas_$(TARGET_CORE) $(MSA_FLAGS)
else ifneq ($(filter NEOVERSEN2 NEOVERSEV1 ARMV8SVE A64FX, $(TARGET_CORE)),)
 ifeq ($(C_COMPILER), PGI)
  override CFLAGS += -DBUILD_KERNEL -DTABLE_NAME=gotoblas_$(TARGET_CORE) -Msve_intrinsics
 else
  override CFLAGS += -DBUILD_KERNEL -DTABLE_NAME=gotoblas_$(TARGET_CORE) -march=armv8.4-a+sve+bf16
 endif
else
 override CFLAGS += -DBUILD_KERNEL -DTABLE_NAME=gotoblas_$(TARGET_CORE)
endif
BUILD_KERNEL = 1
KDIR =
TSUFFIX = _$(TARGET_CORE)
else
TARGET_CORE = $(CORE)
KDIR =
TSUFFIX =
ifeq ($(TARGET_CORE), LOONGSON3R4)
  override CFLAGS += $(MSA_FLAGS)
endif
endif

-include $(KERNELDIR)/KERNEL.$(TARGET_CORE)

# RDNA2/RDNA4 HIP kernels. Built whenever USE_HIP is set, whatever the CPU core
# is: the CPU side keeps that core's own kernels and these add the GPU side, so
# an ordinary build (TARGET=zen, TARGET=haswell, ...) gets both. A GPU-core
# build (TARGET=GFX1201) is just the special case where there is no separate CPU
# core. Defined before Makefile.L3 so SKERNELOBJS picks the objects up (SBLASOBJS
# is expanded from it at parse time).
ifdef USE_HIP
USE_RDNA2 = 1
endif

ifdef USE_RDNA2
HIPCC ?= hipcc
# $(SUFFIX) keeps these out of the profiled object list, which would otherwise
# compile the same file a second time with -pg.
RDNA2_OBJS = rdna2/rdna2_kernels.$(SUFFIX) rdna2/rdna2_oblas.$(SUFFIX)
# GPU_ARCH selects the device architecture(s) independently of the CPU core;
# a GPU-core target supplies it implicitly. One --offload-arch per entry, so
# GPU_ARCH="gfx1031 gfx1201" produces a fat object that runs on both families.
ifneq ($(strip $(GPU_ARCH)),)
RDNA2_GPU_ARCH = $(GPU_ARCH)
else
RDNA2_GPU_ARCH = $(if $(filter GFX%,$(TARGET_CORE)),gfx$(patsubst GFX%,%,$(TARGET_CORE)),gfx1201)
endif
RDNA2_OFFLOAD = $(addprefix --offload-arch=,$(RDNA2_GPU_ARCH))
# The kernels have a single source of truth, gpu/rdna2/, shared with the
# standalone test and benchmark harnesses.
RDNA2_HIP_SRC = $(TOPDIR)/gpu/rdna2/rdna2_kernels.hip
# The architecture list is baked into the device object but make cannot see
# flags, so record it in a stamp the object depends on: changing GPU_ARCH then
# relinks the library instead of silently keeping the previous cards' code.
RDNA2_ARCH_STAMP = rdna2/rdna2_gpu_arch.stamp
$(RDNA2_ARCH_STAMP): FORCE
	@printf '%s\n' '$(RDNA2_GPU_ARCH)' | cmp -s - $@ || printf '%s\n' '$(RDNA2_GPU_ARCH)' > $@
rdna2/rdna2_kernels.$(SUFFIX): $(RDNA2_HIP_SRC) $(RDNA2_ARCH_STAMP)
	$(HIPCC) $(HIPFLAGS) $(RDNA2_OFFLOAD) -x hip -c $< -o $@
rdna2/rdna2_oblas.$(SUFFIX): rdna2/rdna2_oblas.c
	$(CC) $(CFLAGS) -c $< -o $@
SKERNELOBJS += $(RDNA2_OBJS)
endif

# The device objects carry an -offload-arch list inside them, so they must not
# outlive a clean: obj/o are spelled out because $(SUFFIX) is empty here.
clean ::
	@rm -f rdna2/rdna2_gpu_arch.stamp \
		rdna2/rdna2_kernels.obj rdna2/rdna2_kernels.o \
		rdna2/rdna2_oblas.obj rdna2/rdna2_oblas.o

.PHONY: FORCE
FORCE:

include $(KERNELDIR)/KERNEL

include Makefile.L1

include Makefile.L2

include Makefile.L3

include Makefile.LA

HPLOBJS = \
	dgemm_kernel.$(SUFFIX) \
	$(DGEMMINCOPYOBJ) $(DGEMMITCOPYOBJ) \
	$(DGEMMONCOPYOBJ) $(DGEMMOTCOPYOBJ) \
	dtrsm_kernel_LN.$(SUFFIX) dtrsm_kernel_LT.$(SUFFIX) \
	dtrsm_kernel_RN.$(SUFFIX) dtrsm_kernel_RT.$(SUFFIX) \
	daxpy_k.$(SUFFIX) dcopy_k.$(SUFFIX) ddot_k.$(SUFFIX) \
	dger_k.$(SUFFIX)  dscal_k.$(SUFFIX) idamax_k.$(SUFFIX) \
	dgemv_n.$(SUFFIX) dgemv_t.$(SUFFIX) dgemm_beta.$(SUFFIX) \
	dtrsm_iunucopy.$(SUFFIX) dtrsm_iunncopy.$(SUFFIX) \
	dtrsm_ilnucopy.$(SUFFIX) dtrsm_ilnncopy.$(SUFFIX) \
	dtrsm_iutucopy.$(SUFFIX) dtrsm_iutncopy.$(SUFFIX) \
	dtrsm_iltucopy.$(SUFFIX) dtrsm_iltncopy.$(SUFFIX) \
	dtrsm_ounucopy.$(SUFFIX) dtrsm_ounncopy.$(SUFFIX) \
	dtrsm_olnucopy.$(SUFFIX) dtrsm_olnncopy.$(SUFFIX) \
	dtrsm_outucopy.$(SUFFIX) dtrsm_outncopy.$(SUFFIX) \
	dtrsm_oltucopy.$(SUFFIX) dtrsm_oltncopy.$(SUFFIX)

COMMONOBJS	+= lsame.$(SUFFIX) scabs1.$(SUFFIX) dcabs1.$(SUFFIX)

ifeq ($(DYNAMIC_ARCH), 1)
SBLASOBJS	+= setparam$(TSUFFIX).$(SUFFIX)
CCOMMON_OPT	+= -DTS=$(TSUFFIX)
endif

KERNEL_INTERFACE = ../common_level1.h ../common_level2.h ../common_level3.h
ifneq ($(NO_LAPACK), 1)
KERNEL_INTERFACE += ../common_lapack.h
endif

ifeq ($(ARCH), x86)
COMMONOBJS	+= cpuid.$(SUFFIX)
endif

ifdef EXPRECISION
COMMONOBJS	+= qconjg.$(SUFFIX) qcabs1.$(SUFFIX)
endif

ifdef QUAD_PRECISION
COMMONOBJS	+= qconjg.$(SUFFIX) qcabs1.$(SUFFIX)
endif

all : libs

scabs1.$(SUFFIX): $(KERNELDIR)/$(SCABS_KERNEL)
	$(CC) -c $(CFLAGS) -DCOMPLEX -UDOUBLE -DF_INTERFACE $< -o $(@F)

dcabs1.$(SUFFIX): $(KERNELDIR)/$(DCABS_KERNEL)
	$(CC) -c $(CFLAGS) -DCOMPLEX -DDOUBLE -DF_INTERFACE $< -o $(@F)

qcabs1.$(SUFFIX): $(KERNELDIR)/$(QCABS_KERNEL)
	$(CC) -c $(CFLAGS) -DCOMPLEX -DXDOUBLE -DF_INTERFACE $< -o $(@F)

qconjg.$(SUFFIX): $(KERNELDIR)/qconjg.S
	$(CC) -c $(CFLAGS) -DCOMPLEX -DXDOUBLE -DF_INTERFACE $< -o $(@F)

lsame.$(SUFFIX): $(KERNELDIR)/$(LSAME_KERNEL)
	$(CC) -c $(CFLAGS) -DF_INTERFACE $< -o $(@F)

setparam$(TSUFFIX).$(SUFFIX): setparam$(TSUFFIX).c kernel$(TSUFFIX).h
ifeq ($(USE_GEMM3M), 1)
	$(CC) -c $(CFLAGS) -DUSE_GEMM3M $< -o $@
else
	$(CC) -c $(CFLAGS) $< -o $@
endif

setparam$(TSUFFIX).c : setparam-ref.c
	sed 's/TS/$(TSUFFIX)/g' $< > $(@F)

kernel$(TSUFFIX).h : $(KERNEL_INTERFACE)
	sed 's/\ *(/$(TSUFFIX)(/g' $^ > $(@F)


cpuid.$(SUFFIX): $(KERNELDIR)/cpuid.S
	$(CC) -c $(CFLAGS) $< -o $(@F)

scabs1.$(PSUFFIX): $(KERNELDIR)/$(SCABS_KERNEL)
	$(CC) -c $(PFLAGS) -DCOMPLEX -UDOUBLE -DF_INTERFACE $< -o $(@F)

dcabs1.$(PSUFFIX): $(KERNELDIR)/$(DCABS_KERNEL)
	$(CC) -c $(PFLAGS) -DCOMPLEX -DDOUBLE -DF_INTERFACE $< -o $(@F)

qcabs1.$(PSUFFIX): $(KERNELDIR)/$(QCABS_KERNEL)
	$(CC) -c $(PFLAGS) -DCOMPLEX -DXDOUBLE -DF_INTERFACE $< -o $(@F)

qconjg.$(PSUFFIX): $(KERNELDIR)/qconjg.S
	$(CC) -c $(PFLAGS) -DCOMPLEX -DXDOUBLE -DF_INTERFACE $< -o $(@F)

lsame.$(PSUFFIX): $(KERNELDIR)/$(LSAME_KERNEL)
	$(CC) -c $(PFLAGS) -DF_INTERFACE $< -o $(@F)

cpuid.$(PSUFFIX): $(KERNELDIR)/cpuid.S
	$(CC) -c $(PFLAGS) $< -o $(@F)

#ifdef DYNAMIC_ARCH
clean ::
	@rm -f setparam_*.c kernel_*.h setparam.h kernel.h

#endif

include $(TOPDIR)/Makefile.tail
