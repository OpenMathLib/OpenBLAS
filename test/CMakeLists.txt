include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_BINARY_DIR})

enable_language(Fortran)
if (CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-tree-vectorize")
endif()
if (CMAKE_Fortran_COMPILER_ID STREQUAL Flang)
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2")
endif()

if (BUILD_SINGLE)
	list( APPEND OpenBLAS_Tests sblat1 sblat2 sblat3)
endif()
if (BUILD_DOUBLE)
	list (APPEND OpenBLAS_Tests dblat1 dblat2 dblat3)
endif()
if (BUILD_COMPLEX)
	list (APPEND OpenBLAS_Tests cblat1 cblat2 cblat3)
endif()
if (BUILD_COMPLEX16)
	list (APPEND OpenBLAS_Tests zblat1 zblat2 zblat3)
endif()

if (USE_GEMM3M)
	if (BUILD_COMPLEX)
		list (APPEND OpenBLAS_Tests cblat3_3m)
	endif ()
	if (BUILD_COMPLEX16)
		list (APPEND OpenBLAS_Tests zblat3_3m)
	endif ()
endif ()

add_library(openblas_test_xerbla OBJECT xerbla_test.c)
# Static test executables provide their own XERBLA and must continue to test
# link-time interposition. Shared-only tests use the explicit handler because
# calls originating in a DLL or dylib cannot portably bind to that definition.
if (BUILD_SHARED_LIBS AND NOT BUILD_STATIC_LIBS)
	target_compile_definitions(openblas_test_xerbla PRIVATE
		OPENBLAS_TEST_USE_XERBLA_HANDLER)
endif()

foreach(test_bin ${OpenBLAS_Tests})
	add_executable(${test_bin} ${test_bin}.f)
	if (test_bin MATCHES "blat(2|3)(_3m)?$")
		target_sources(${test_bin} PRIVATE
			$<TARGET_OBJECTS:openblas_test_xerbla>)
	endif()
	target_link_libraries(${test_bin} ${OpenBLAS_LIBNAME})
endforeach()

if (BUILD_BFLOAT16)
	add_executable(test_bgemm compare_sgemm_bgemm.c)
	target_compile_definitions(test_bgemm PUBLIC -DIBFLOAT16 -DOBFLOAT16)
	target_link_libraries(test_bgemm ${OpenBLAS_LIBNAME})
	add_executable(test_bgemv compare_sgemv_bgemv.c)
	target_compile_definitions(test_bgemv PUBLIC -DIBFLOAT16 -DOBFLOAT16)
	target_link_libraries(test_bgemv ${OpenBLAS_LIBNAME})
	add_executable(test_sbgemm compare_sgemm_sbgemm.c)
	target_compile_definitions(test_sbgemm PUBLIC -DIBFLOAT16)
	target_link_libraries(test_sbgemm ${OpenBLAS_LIBNAME})
	add_executable(test_sbgemv compare_sgemv_sbgemv.c)
	target_compile_definitions(test_sbgemv PUBLIC -DIBFLOAT16)
	target_link_libraries(test_sbgemv ${OpenBLAS_LIBNAME})
endif()

if (BUILD_HFLOAT16)
	add_executable(test_shgemm compare_sgemm_shgemm.c)
	target_link_libraries(test_shgemm ${OpenBLAS_LIBNAME})
	add_executable(test_shgemv compare_sgemv_shgemv.c)
	target_link_libraries(test_shgemv ${OpenBLAS_LIBNAME})
endif()

# $1 exec, $2 input, $3 output_result
if(WIN32)
FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/test_helper.ps1
"[Console]::InputEncoding = New-Object Text.UTF8Encoding $false\n"
"$ErrorActionPreference = \"Stop\"\n"
"if (Test-Path -LiteralPath $args[2]) { Remove-Item -LiteralPath $args[2] -Force }\n"
"Get-Content -LiteralPath $args[1] | & $args[0]\n"
"$testExitCode = $LASTEXITCODE\n"
"if ($testExitCode -ne 0) { exit $testExitCode }\n"
"if (-not (Test-Path -LiteralPath $args[2])) {\n"
"  Write-Host \"Missing test summary: $($args[2])\"\n"
"  exit 1\n"
"}\n"
"if (Select-String -LiteralPath $args[2] -Pattern \"FATAL|FAILED\" -Quiet) {\n"
"  Write-Host \"Error in $($args[1]):\"\n"
"  Get-Content -LiteralPath $args[2]\n"
"  exit 1\n"
"}\n"
"exit 0\n"
)
set(helper_prefix powershell -ExecutionPolicy Bypass "${CMAKE_CURRENT_BINARY_DIR}/test_helper.ps1")
else()
FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/test_helper.sh
"rm -f \"$3\"\n"
"remove_status=$?\n"
"if [ $remove_status -ne 0 ]; then\n"
"  echo \"Unable to remove stale test summary: $3\"\n"
"  exit $remove_status\n"
"fi\n"
"\"$1\" < \"$2\"\n"
"test_exit_code=$?\n"
"if [ $test_exit_code -ne 0 ]; then exit $test_exit_code; fi\n"
"if [ ! -f \"$3\" ]; then\n"
"  echo \"Missing test summary: $3\"\n"
"  exit 1\n"
"fi\n"
"grep -q -e FATAL -e FAILED \"$3\"\n"
"summary_status=$?\n"
"case $summary_status in\n"
"  0) echo \"Error in $2:\"; cat \"$3\"; exit 1 ;;\n"
"  1) exit 0 ;;\n"
"  *) echo \"Unable to read test summary: $3\"; exit $summary_status ;;\n"
"esac\n"
)
set(helper_prefix sh "${CMAKE_CURRENT_BINARY_DIR}/test_helper.sh")
endif()

#set(float_types s d c z)
if (BUILD_SINGLE)
	list (APPEND float_types s)
endif()
if (BUILD_DOUBLE)
	list (APPEND float_types d)
endif()
if (BUILD_COMPLEX)
	list (APPEND float_types c)
endif()
if (BUILD_COMPLEX16)
	list (APPEND float_types z)
endif()
foreach(float_type ${float_types})
string(TOUPPER ${float_type} float_type_upper)
add_test(NAME "${float_type}blas1" 
  COMMAND $<TARGET_FILE:${float_type}blat1>)
add_test(NAME "${float_type}blas2" 
  COMMAND ${helper_prefix} $<TARGET_FILE:${float_type}blat2> "${PROJECT_SOURCE_DIR}/test/${float_type}blat2.dat" ${float_type_upper}BLAT2.SUMM)
add_test(NAME "${float_type}blas3" 
  COMMAND ${helper_prefix} $<TARGET_FILE:${float_type}blat3> "${PROJECT_SOURCE_DIR}/test/${float_type}blat3.dat" ${float_type_upper}BLAT3.SUMM)
if (USE_GEMM3M)
if ((${float_type} STREQUAL "c") OR (${float_type} STREQUAL "z"))
add_test(NAME "${float_type}blas3_3m" 
  COMMAND ${helper_prefix} $<TARGET_FILE:${float_type}blat3_3m> "${PROJECT_SOURCE_DIR}/test/${float_type}blat3_3m.dat" ${float_type_upper}BLAT3_3M.SUMM)
endif()
endif()
endforeach()

if (BUILD_BFLOAT16)
	add_test(NAME "bgemm" 
	COMMAND $<TARGET_FILE:test_bgemm>)
	add_test(NAME "bgemv" 
	COMMAND $<TARGET_FILE:test_bgemv>)
	add_test(NAME "sbgemm" 
	COMMAND $<TARGET_FILE:test_sbgemm>)
	add_test(NAME "sbgemv" 
	COMMAND $<TARGET_FILE:test_sbgemv>)
endif()

if (BUILD_HFLOAT16)
	add_test(NAME "shgemm" 
	COMMAND $<TARGET_FILE:test_shgemm>)
	add_test(NAME "shgemv" 
	COMMAND $<TARGET_FILE:test_shgemv>)
endif()
