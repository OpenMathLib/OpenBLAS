include_directories(${CMAKE_SOURCE_DIR})

if (${CORE} STREQUAL "PPC440")
  set(MEMORY memory_qalloc.c)
else ()
  set(MEMORY memory.c)
endif ()

add_library(COMMON_OBJS OBJECT
  ${MEMORY}
  xerbla.c
  abs.c # TODO: this is split into c_abs (DOUBLE unset) and z_abs (DOUBLE set) in the Makefile
  openblas_set_num_threads.c
  openblas_get_config.c
  openblas_get_parallel.c
  openblas_error_handle.c
)

#ifdef SMP
#COMMONOBJS	+= blas_server.$(SUFFIX) divtable.$(SUFFIX) blasL1thread.$(SUFFIX)
#ifndef NO_AFFINITY
#COMMONOBJS	+= init.$(SUFFIX)
#endif
#endif
#
#ifeq ($(DYNAMIC_ARCH), 1)
#COMMONOBJS	+=  dynamic.$(SUFFIX)
#else
#COMMONOBJS	+=  parameter.$(SUFFIX)
#endif
#
#ifdef EXPRECISION
#COMMONOBJS	+= x_abs.$(SUFFIX) qlamch.$(SUFFIX) qlamc3.$(SUFFIX)
#endif
#
#ifdef QUAD_PRECISION
#COMMONOBJS	+= addx.$(SUFFIX) mulx.$(SUFFIX)
#endif
#
#ifdef USE_CUDA
#COMMONOBJS	+= cuda_init.$(SUFFIX)
#endif
#
#ifdef FUNCTION_PROFILE
#COMMONOBJS	+= profile.$(SUFFIX)
#endif

list(APPEND DBLAS_OBJS "COMMON_OBJS")

#LIBOTHERS = libothers.$(LIBSUFFIX)
#
#ifeq ($(USE_OPENMP), 1)
#BLAS_SERVER = blas_server_omp.c
#else
#ifeq ($(OSNAME), WINNT)
#BLAS_SERVER = blas_server_win32.c
#endif
#ifeq ($(OSNAME), CYGWIN_NT)
#BLAS_SERVER = blas_server_win32.c
#endif
#ifeq ($(OSNAME), Interix)
#BLAS_SERVER = blas_server_win32.c
#endif
#endif
#
#ifndef BLAS_SERVER
#BLAS_SERVER = blas_server.c
#endif
#
#ifeq ($(DYNAMIC_ARCH), 1)
#HPLOBJS = memory.$(SUFFIX) xerbla.$(SUFFIX) dynamic.$(SUFFIX)
#else
#HPLOBJS = memory.$(SUFFIX) xerbla.$(SUFFIX) parameter.$(SUFFIX)
#endif

set(DBLAS_OBJS ${DBLAS_OBJS} PARENT_SCOPE) # list append removes the scope from DBLAS_OBJS

