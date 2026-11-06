# Build the RDNA2 / RDNA4 int8/f32 kernels and their correctness test.
#
#   make                          build test_rdna2 for gfx1031
#   make test                     build and run it (GPU if present, else CPU)
#   make HIP_ARCH=gfx1030 test    build for gfx1030 instead
#   make HIP_ARCH=gfx1201 test    build for gfx1201 (RX 9070 / 9070 XT)
#   make ARCHS="gfx1031 gfx1201" test
#                                 one fat binary covering several GPUs
#   make multi                    fat binary for gfx1030 gfx1031 gfx1200 gfx1201
#   make oblas-test               run the level-3 bridge test (see below)
#   make clean
#
# Targets: gfx1031 = RX 6700 XT / 6750 XT, gfx1030 = RX 6800 / 6900,
#          gfx1201 = RX 9070 / 9070 XT (Navi 48), gfx1200 = RX 9060 XT (Navi 44).
# gfx1200 and gfx1201 are not binary compatible; a fat binary compiled for both
# (or gfx12-generic) is what runs on either.
#
# Override the compiler if your ROCm install names it differently:
#   make HIPCC=/opt/rocm/llvm/bin/clang++ INC=-I/opt/rocm/include
#
# Changing HIP_ARCH/ARCHS does not rebuild existing objects: run make clean
# first, or use "make multi", which forces a fresh build for its arch set.
# On gfx12, if the toolchain rejects v_dot4_i32_i8 ("needs target feature
# dot1-insts"), rebuild with HIPFLAGS="-O2 -DRDNA2_NO_DOT1=1" to use the scalar
# sign-extend/multiply fallback instead.
HIPCC   ?= amdclang++
HIP_ARCH ?= gfx1031
ARCHS    ?= $(HIP_ARCH)
HIPFLAGS ?= -O2
INC      ?=
LDFLAGS  ?=

CXXFLAGS_COMMON = $(HIPFLAGS) -std=c++17 -Wall -Wextra $(INC)
OFFLOAD = $(addprefix --offload-arch=,$(ARCHS))

all: test_rdna2

rdna2_kernels.o: rdna2_kernels.hip rdna2_gemm_ref.h rdna2_kernels_math.h
	$(HIPCC) -x hip $(OFFLOAD) $(CXXFLAGS_COMMON) -c $< -o $@

test_rdna2: test_rdna2_host.cpp rdna2_kernels.o rdna2_gemm_ref.h rdna2_kernels_math.h
	$(HIPCC) -x hip $(OFFLOAD) $(CXXFLAGS_COMMON) \
		test_rdna2_host.cpp -x none rdna2_kernels.o $(LDFLAGS) -o $@

test: test_rdna2
	./test_rdna2

# CPU vs GPU comparison. CPU side is this file's blocked -O3 gemm plus the
# OpenBLAS DLL shipped with numpy when OB_DLL points at it (or is autodetected).
OB_DLL ?= $(shell python -c "import glob,os,numpy;print((glob.glob(os.path.join(os.path.dirname(numpy.__file__),'..','numpy.libs','*openblas*.dll'))+[''])[0])" 2>/dev/null)

bench: rdna2_kernels.o
	$(HIPCC) -x hip $(OFFLOAD) $(CXXFLAGS_COMMON) \
		bench_rdna2.cpp -x none rdna2_kernels.o $(LDFLAGS) -o bench_rdna2

bench-run: bench
	./bench_rdna2 "$(OB_DLL)"

# Level-3 bridge test against the library built in this tree: checks the
# interface-side kernel plumbing (SGEMM, SYMM, TRMM, TRSM, SYRK and SYR2K where
# the dispatch is compiled in) against an independent double-precision reference,
# including the contract that the unstored triangle and the ldb/ldc padding are
# untouched, over every layout and over zero, negative and large scalars.
# The same binary is useful either way -- the device path is exercised only when
# OBLAS_LIB is a USE_HIP build, while the reference checks still run on a
# CPU-only build. Point OBLAS_LIB/LDFLAGS at another tree or a ROCm runtime
# outside the default search path if needed:
#   make oblas-test OBLAS_LIB=/path/libopenblas.a LDFLAGS=-L/opt/rocm/lib
CC          ?= cc
AR          ?= ar
OBLAS_INC   ?= ../..
OBLAS_LIB   ?= ../../libopenblas.a
OBLAS_CFLAGS?= -O2 -Wall -I$(OBLAS_INC)
# A USE_HIP library carries the device objects and needs the HIP runtime at link
# time; a plain CPU library does not, so the runtime stays out of the link line
# on machines (and CI) without ROCm instead of failing there.
OBLAS_HIP   := $(shell $(AR) t $(OBLAS_LIB) 2>/dev/null | grep -q rdna2_kernels && echo 1)
ifeq ($(OBLAS_HIP),1)
# HIP_LIB_PATH comes from the top-level build (Makefile.system exports it) when
# this is run through "make level3-test" there.
OBLAS_LIBS  ?= $(if $(HIP_LIB_PATH),-L$(subst \,/,$(HIP_LIB_PATH))) -lamdhip64 -lstdc++ -lpthread -lm
else
OBLAS_LIBS  ?= -lpthread -lm
endif

test_oblas_gpu: test_oblas_gpu.c
	$(CC) $(OBLAS_CFLAGS) $< $(OBLAS_LIB) $(LDFLAGS) $(OBLAS_LIBS) -o $@

oblas-test: test_oblas_gpu
	./test_oblas_gpu

# Every target this tree is shaped for, in one binary. -B because the object
# files are otherwise reused from whatever arch the last build used.
multi:
	$(MAKE) -B ARCHS="gfx1030 gfx1031 gfx1200 gfx1201" test

clean:
	rm -f rdna2_kernels.o test_rdna2 bench_rdna2 test_oblas_gpu

.PHONY: all test bench bench-run oblas-test multi clean
