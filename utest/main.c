/*****************************************************************************
Copyright (c) 2011, Lab of Parallel Software and Computational Science,ICSAS
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.
   3. Neither the name of the ISCAS nor the names of its contributors may 
      be used to endorse or promote products derived from this software 
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

**********************************************************************************/

#include <stdio.h>
#include <string.h>

#include "common_utest.h"
#include <CUnit/Basic.h>

CU_TestInfo test_level1[]={
	{"Testing srot when incx || incy == 0",test_srot_inc_0},
	{"Testing drot when incx || incy == 0",test_drot_inc_0},
	{"Testing csrot when incx || incy == 0",test_csrot_inc_0},
	{"Testing zdrot when incx || incy == 0",test_zdrot_inc_0},

	{"Testing sswap with incx || incy == 0",test_sswap_inc_0},
	{"Testing dswap with incx || incy == 0",test_dswap_inc_0},
	{"Testing cswap with incx || incy == 0",test_cswap_inc_0},
	{"Testing zswap with incx || incy == 0",test_zswap_inc_0},

	{"Testing saxpy with incx || incy == 0",test_saxpy_inc_0},
	{"Testing daxpy with incx || incy == 0",test_daxpy_inc_0},
	{"Testing caxpy with incx || incy == 0",test_caxpy_inc_0},
	{"Testing zaxpy with incx || incy == 0",test_zaxpy_inc_0},

	{"Testing zdotu with n == 1",test_zdotu_n_1},
	{"Testing zdotu with input x & y offset == 1",test_zdotu_offset_1},

	{"Testing drotmg",test_drotmg},

	{"Testing dsdot with n == 1",test_dsdot_n_1},

	{"Testing samax", test_samax},
	CU_TEST_INFO_NULL,
};

CU_SuiteInfo suites[]={
	{"Level1 Test Suite", NULL,NULL,test_level1},
	CU_SUITE_INFO_NULL,
};

int main()
{
	CU_ErrorCode error;
	if (CUE_SUCCESS != CU_initialize_registry())
		return CU_get_error();
	
	error=CU_register_suites(suites);
	
	if (error != CUE_SUCCESS) {
		perror(CU_get_error_msg());
		CU_cleanup_registry();
		return CU_get_error();
		
	}
	

	
	printf("Seting OK\n");
	fflush(stdout);
	
	/* Run all tests using the CUnit Basic interface */
	CU_basic_set_mode(CU_BRM_VERBOSE);
	
	CU_basic_run_tests();
	
	CU_cleanup_registry();
	
	return CU_get_error();
	
}

