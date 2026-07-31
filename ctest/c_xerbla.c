#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>
#include <string.h>
#include "common.h"
#include "cblas_test.h"

void cblas_xerbla(blasint info, char *rout, char *form, ...)
{
   extern int cblas_lerr, cblas_info, cblas_ok;
   extern int RowMajorStrg;
   extern char *cblas_rout;

   if (cblas_rout != NULL && strcmp(cblas_rout, rout) != 0){
      printf("***** XERBLA WAS CALLED WITH SRNAME = <%s> INSTEAD OF <%s> *******\n", rout, cblas_rout);
      cblas_ok = FALSE;
   }

   if (RowMajorStrg)
   {
      /* To properly check leading dimension problems in cblas__gemm, we
       * need to do the following trick. When cblas__gemm is called with
       * CblasRowMajor, the arguments A and B switch places in the call to
       * f77__gemm. Thus when we test for bad leading dimension problems
       * for A and B, lda is in position 11 instead of 9, and ldb is in
       * position 9 instead of 11.
       */
      if (strstr(rout,"gemm") != 0)
      {
         if      (info == 5 ) info =  4;
         else if (info == 4 ) info =  5;
         else if (info == 11) info =  9;
         else if (info == 9 ) info = 11;
      }
      else if (strstr(rout,"symm") != 0 || strstr(rout,"hemm") != 0)
      {
         if      (info == 5 ) info =  4;
         else if (info == 4 ) info =  5;
      }
      else if (strstr(rout,"trmm") != 0 || strstr(rout,"trsm") != 0)
      {
         if      (info == 7 ) info =  6;
         else if (info == 6 ) info =  7;
      }
      else if (strstr(rout,"gemv") != 0)
      {
         if      (info == 4)  info = 3;
         else if (info == 3)  info = 4;
      }
      else if (strstr(rout,"gbmv") != 0)
      {
         if      (info == 4)  info = 3;
         else if (info == 3)  info = 4;
         else if (info == 6)  info = 5;
         else if (info == 5)  info = 6;
      }
      else if (strstr(rout,"ger") != 0)
      {
         if      (info == 3) info = 2;
         else if (info == 2) info = 3;
         else if (info == 8) info = 6;
         else if (info == 6) info = 8;
      }
      else if ( ( strstr(rout,"her2") != 0 || strstr(rout,"hpr2") != 0 )
               && strstr(rout,"her2k") == 0 )
      {
         if      (info == 8) info = 6;
         else if (info == 6) info = 8;
      }
   }

   if (info != cblas_info){
      printf("***** XERBLA WAS CALLED WITH INFO = %lld INSTEAD OF %lld in %s *******\n",
             (long long)info, (long long)cblas_info, rout);
      cblas_lerr = PASSED;
      cblas_ok = FALSE;
   } else cblas_lerr = FAILED;
}

static void cblas_test_xerbla(const char *srname, const blasint *info,
                              size_t length)
{
   extern int cblas_ok;
   char rout[] = {'c','b','l','a','s','_','\0','\0','\0','\0','\0','\0','\0'};
   blasint i;

   if (length < 6) {
      printf("***** XERBLA WAS CALLED WITH AN INVALID ROUTINE NAME LENGTH *******\n");
      cblas_ok = FALSE;
      return;
   }

   for(i=0;  i  < 6; i++) rout[i+6] = tolower((unsigned char)srname[i]);
   for(i=11; i >= 9; i--) if (rout[i] == ' ') rout[i] = '\0';

   /* We increment *info by 1 since the CBLAS interface adds one more
    * argument to all level 2 and 3 routines.
    */
   cblas_xerbla(*info+1,rout,"");
}

void cblas_test_set_xerbla(void) {
   openblas_set_xerbla(cblas_test_xerbla);
}

void cblas_test_fail(void) {
   exit(EXIT_FAILURE);
}
