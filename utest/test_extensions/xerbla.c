/*****************************************************************************
Copyright (c) 2023, The OpenBLAS Project
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
   3. Neither the name of the OpenBLAS project nor the names of
      its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

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

#include <limits.h>
#include "utest/openblas_utest.h"
#include "common.h"

static int handler_installed = FALSE;
static int lerr, _info, ok;
static char *rout;

static void test_xerbla(const char *srname, const blasint *vinfo, size_t length)
{
   blasint info = *vinfo;
   int name_length = length > (size_t)INT_MAX ? INT_MAX : (int)length;

   if (rout != NULL &&
       (length != strlen(rout) || memcmp(rout, srname, length) != 0)) {
      printf("***** XERBLA WAS CALLED WITH AN UNEXPECTED SRNAME INSTEAD OF <%s> *******\n",
             rout);
      ok = FALSE;
   }

   if (info != _info){
      printf("***** XERBLA WAS CALLED WITH INFO = %lld INSTEAD OF %d in %.*s *******\n",
             (long long)info, _info, name_length, srname == NULL ? "" : srname);
      lerr = TRUE;
      ok = FALSE;
   } else lerr = FALSE;
}

static void alternate_xerbla(const char *srname, const blasint *vinfo,
                             size_t length)
{
   (void) srname;
   (void) vinfo;
   (void) length;
}

CTEST(openblas_extensions, xerbla_handler_registration)
{
   openblas_xerbla_handler original;
   openblas_xerbla_handler previous;
   openblas_xerbla_handler restored;
   openblas_xerbla_handler default_handler;

   original = openblas_set_xerbla(test_xerbla);
   previous = openblas_set_xerbla(alternate_xerbla);
   restored = openblas_set_xerbla(NULL);
   default_handler = openblas_set_xerbla(original);

   ASSERT_TRUE(previous == test_xerbla);
   ASSERT_TRUE(restored == alternate_xerbla);
   ASSERT_TRUE(default_handler != NULL);
}

int check_error(void) {
   if (lerr == TRUE ) {
       printf("***** ILLEGAL VALUE OF PARAMETER NUMBER %d NOT DETECTED BY %s *****\n", _info, rout);
      ok = FALSE;
   }
   lerr = TRUE;
   return ok;
}

void set_xerbla(char* current_rout, int expected_info){
   if (!handler_installed) {
      openblas_set_xerbla(test_xerbla);
      handler_installed = TRUE;
   }

   ok = TRUE;
   lerr = TRUE;
   _info = expected_info;
   rout = current_rout;
}
