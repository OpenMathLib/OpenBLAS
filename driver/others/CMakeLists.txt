include_directories(${CMAKE_SOURCE_DIR})

if (${CORE} STREQUAL "PPC440")
  set(MEMORY memory_qalloc.c)
else ()
  set(MEMORY memory.c)
endif ()

if (SMP)

  if (USE_OPENMP)
    set(BLAS_SERVER blas_server_omp.c)
  elseif (${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    set(BLAS_SERVER blas_server_win32.c)
  endif ()

  if (NOT DEFINED BLAS_SERVER)
    set(BLAS_SERVER blas_server.c)
  endif ()

  set(SMP_SOURCES
    ${BLAS_SERVER}
    divtable.c # TODO: Makefile has -UDOUBLE
    blas_l1_thread.c
  )

  if (NOT NO_AFFINITY)
    list(APPEND SMP_SOURCES init.c)
  endif ()
endif ()

set(COMMON_SOURCES
  xerbla.c
  openblas_set_num_threads.c
  openblas_error_handle.c
)

# these need to have NAME/CNAME set, so use GenerateNamedObjects
GenerateNamedObjects("abs.c" "" "" "c_abs")
GenerateNamedObjects("abs.c" "" "DOUBLE" "z_abs")
GenerateNamedObjects("openblas_get_config.c;openblas_get_parallel.c")

if (DYNAMIC_ARCH)
  list(APPEND COMMON_SOURCES dynamic.c)
else ()
  list(APPEND COMMON_SOURCES parameter.c)
endif ()

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

add_library(COMMON_OBJS OBJECT
  ${MEMORY}
  ${SMP_SOURCES}
  ${COMMON_SOURCES}
)

list(APPEND DBLAS_OBJS "COMMON_OBJS")

#LIBOTHERS = libothers.$(LIBSUFFIX)

#ifeq ($(DYNAMIC_ARCH), 1)
#HPLOBJS = memory.$(SUFFIX) xerbla.$(SUFFIX) dynamic.$(SUFFIX)
#else
#HPLOBJS = memory.$(SUFFIX) xerbla.$(SUFFIX) parameter.$(SUFFIX)
#endif

set(DBLAS_OBJS ${DBLAS_OBJS} PARENT_SCOPE) # list append removes the scope from DBLAS_OBJS

