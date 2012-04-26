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

#ifndef COMMON_UTEST_H_
#define COMMON_UTEST_H_
#include <CUnit/CUnit.h>

#include <common.h>

#define CHECK_EPS 0.00002

//Testcase list
void test_drot_inc_0(void);
void test_srot_inc_0(void);
void test_zdrot_inc_0(void);
void test_csrot_inc_0(void);

void test_dswap_inc_0(void);
void test_zswap_inc_0(void);
void test_sswap_inc_0(void);
void test_cswap_inc_0(void);

void test_daxpy_inc_0(void);
void test_zaxpy_inc_0(void);
void test_saxpy_inc_0(void);
void test_caxpy_inc_0(void);

void test_zdotu_n_1(void);
void test_zdotu_offset_1(void);

void test_drotmg(void);

void test_dsdot_n_1(void);

void test_samax(void);

#endif
