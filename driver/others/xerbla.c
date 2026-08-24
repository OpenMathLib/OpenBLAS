/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* All rights reserved.                                              */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/*   1. Redistributions of source code must retain the above         */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer.                                                  */
/*                                                                   */
/*   2. Redistributions in binary form must reproduce the above      */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer in the documentation and/or other materials       */
/*      provided with the distribution.                              */
/*                                                                   */
/*    THIS  SOFTWARE IS PROVIDED  BY THE  UNIVERSITY OF  TEXAS AT    */
/*    AUSTIN  ``AS IS''  AND ANY  EXPRESS OR  IMPLIED WARRANTIES,    */
/*    INCLUDING, BUT  NOT LIMITED  TO, THE IMPLIED  WARRANTIES OF    */
/*    MERCHANTABILITY  AND FITNESS FOR  A PARTICULAR  PURPOSE ARE    */
/*    DISCLAIMED.  IN  NO EVENT SHALL THE UNIVERSITY  OF TEXAS AT    */
/*    AUSTIN OR CONTRIBUTORS BE  LIABLE FOR ANY DIRECT, INDIRECT,    */
/*    INCIDENTAL,  SPECIAL, EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES    */
/*    (INCLUDING, BUT  NOT LIMITED TO,  PROCUREMENT OF SUBSTITUTE    */
/*    GOODS  OR  SERVICES; LOSS  OF  USE,  DATA,  OR PROFITS;  OR    */
/*    BUSINESS INTERRUPTION) HOWEVER CAUSED  AND ON ANY THEORY OF    */
/*    LIABILITY, WHETHER  IN CONTRACT, STRICT  LIABILITY, OR TORT    */
/*    (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY WAY OUT    */
/*    OF  THE  USE OF  THIS  SOFTWARE,  EVEN  IF ADVISED  OF  THE    */
/*    POSSIBILITY OF SUCH DAMAGE.                                    */
/*                                                                   */
/* The views and conclusions contained in the software and           */
/* documentation are those of the authors and should not be          */
/* interpreted as representing official policies, either expressed   */
/* or implied, of The University of Texas at Austin.                 */
/*********************************************************************/

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"

#if defined(OS_WINDOWS) && (defined(__MINGW32__) || defined(__MINGW64__))
#include <conio.h>
#undef  printf
#define printf	_cprintf
#endif

#define MSGFMT " ** On entry to %6.*s parameter number %2lld had an illegal value\n"

static size_t openblas_xerbla_name_length(const char *message,
                                          blasint length) {
  const char *terminator;

  if (message == NULL || length <= 0) return 0;

  terminator = memchr(message, '\0', (size_t)length);
  if (terminator != NULL) return (size_t)(terminator - message);

  return (size_t)length;
}

static void openblas_xerbla_default(const char *message, const blasint *info,
                                    size_t length) {
  int precision = length > INT_MAX ? INT_MAX : (int)length;

  printf(MSGFMT, precision, message == NULL ? "" : message,
         (long long)*info);
}

static openblas_xerbla_handler openblas_xerbla = openblas_xerbla_default;
static volatile BLASULONG openblas_xerbla_lock = 0;

openblas_xerbla_handler
openblas_set_xerbla(openblas_xerbla_handler handler) {
  openblas_xerbla_handler previous;

  if (handler == NULL) handler = openblas_xerbla_default;

  blas_lock(&openblas_xerbla_lock);
  previous = openblas_xerbla;
  openblas_xerbla = handler;
  blas_unlock(&openblas_xerbla_lock);

  return previous;
}

static int openblas_xerbla_dispatch(char *message, blasint *info,
                                    blasint length) {
  openblas_xerbla_handler handler;
  size_t name_length = openblas_xerbla_name_length(message, length);

  blas_lock(&openblas_xerbla_lock);
  handler = openblas_xerbla;
  blas_unlock(&openblas_xerbla_lock);

  handler(message, info, name_length);
  return 0;
}

#ifdef __ELF__
int __xerbla(char *message, blasint *info, blasint length) {
  return openblas_xerbla_dispatch(message, info, length);
}

int BLASFUNC(xerbla)(char *, blasint *, blasint)
  __attribute__ ((weak, alias ("__xerbla")));
#else
int BLASFUNC(xerbla)(char *message, blasint *info, blasint length) {
  return openblas_xerbla_dispatch(message, info, length);
}
#endif
