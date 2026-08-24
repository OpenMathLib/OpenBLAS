include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR})

enable_language(CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -DADD${BU} -DCBLAS")

set(CPP_THREAD_SAFETY_LIBS ${OpenBLAS_LIBNAME})
find_package(Threads REQUIRED)
list(APPEND CPP_THREAD_SAFETY_LIBS Threads::Threads)
add_definitions(-DOPENBLAS_USE_GENERATED_CBLAS_H)

if ((CPP_THREAD_SAFETY_TEST OR CPP_THREAD_SAFETY_GEMV) AND CPP_THREAD_SAFETY_USE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
  list(APPEND CPP_THREAD_SAFETY_LIBS OpenMP::OpenMP_CXX)
  add_definitions(-DCPP_THREAD_SAFETY_USE_OPENMP)
endif()

set(CPP_THREAD_SAFETY_DGEMM_ARGS "" CACHE STRING "Arguments passed to the DGEMM thread safety test")
set(CPP_THREAD_SAFETY_DGEMM_MIXED_ARGS "" CACHE STRING "Arguments passed to the mixed DGEMM thread safety test")
set(CPP_THREAD_SAFETY_DGEMV_ARGS "" CACHE STRING "Arguments passed to the DGEMV thread safety test")
set(CPP_THREAD_SAFETY_SHUTDOWN_ARGS "" CACHE STRING "Arguments passed to the DGEMM shutdown safety test")

if (CPP_THREAD_SAFETY_TEST)
  message(STATUS "building thread safety test")
  add_executable(dgemm_thread_safety dgemm_thread_safety.cpp)
  target_link_libraries(dgemm_thread_safety ${CPP_THREAD_SAFETY_LIBS})
  add_test(NAME dgemm_thread_safety COMMAND ${CMAKE_CURRENT_BINARY_DIR}/dgemm_thread_safety ${CPP_THREAD_SAFETY_DGEMM_ARGS})

  add_executable(dgemm_thread_safety_mixed dgemm_thread_safety_mixed.cpp)
  target_link_libraries(dgemm_thread_safety_mixed ${CPP_THREAD_SAFETY_LIBS})
  add_test(NAME dgemm_thread_safety_mixed COMMAND ${CMAKE_CURRENT_BINARY_DIR}/dgemm_thread_safety_mixed ${CPP_THREAD_SAFETY_DGEMM_MIXED_ARGS})

  if (WIN32)
    add_executable(dgemm_thread_safety_shutdown dgemm_thread_safety_shutdown.cpp)
    target_link_libraries(dgemm_thread_safety_shutdown ${CPP_THREAD_SAFETY_LIBS})
    add_test(NAME dgemm_thread_safety_shutdown COMMAND ${CMAKE_CURRENT_BINARY_DIR}/dgemm_thread_safety_shutdown ${CPP_THREAD_SAFETY_SHUTDOWN_ARGS})
    set_tests_properties(dgemm_thread_safety_shutdown PROPERTIES TIMEOUT 900)
  endif()
  if (USE_THREAD AND (USE_OPENMP OR (NOT WIN32 AND NOT CYGWIN)))
    add_test(NAME dgemm_thread_safety_mixed_callback COMMAND ${CMAKE_CURRENT_BINARY_DIR}/dgemm_thread_safety_mixed ${CPP_THREAD_SAFETY_DGEMM_MIXED_ARGS} --callback)
  endif()
endif()


if (CPP_THREAD_SAFETY_TEST OR CPP_THREAD_SAFETY_GEMV)
  add_executable(dgemv_thread_safety dgemv_thread_safety.cpp)
  target_link_libraries(dgemv_thread_safety ${CPP_THREAD_SAFETY_LIBS})
  add_test(NAME dgemv_thread_safety COMMAND ${CMAKE_CURRENT_BINARY_DIR}/dgemv_thread_safety ${CPP_THREAD_SAFETY_DGEMV_ARGS})
endif()
