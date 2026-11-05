# Build the RDNA2 / RDNA4 int8/f32 kernels and their correctness test.
#
#   make                          build test_rdna2 for gfx1031
#   make test                     build and run it (GPU if present, else CPU)
#   make HIP_ARCH=gfx1030 test    build for gfx1030 instead
#   make HIP_ARCH=gfx1201 test    build for gfx1201 (RX 9070 / 9070 XT)
#   make ARCHS="gfx1031 gfx1201" test
#                                 one fat binary covering several GPUs
#   make multi                    fat binary for gfx1030 gfx1031 gfx1200 gfx1201
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

# Every target this tree is shaped for, in one binary. -B because the object
# files are otherwise reused from whatever arch the last build used.
multi:
	$(MAKE) -B ARCHS="gfx1030 gfx1031 gfx1200 gfx1201" test

clean:
	rm -f rdna2_kernels.o test_rdna2

.PHONY: all test multi clean
