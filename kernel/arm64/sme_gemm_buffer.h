/***************************************************************************
Copyright (c) 2026, The OpenBLAS Project
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
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

/* Packing buffers of the SME GEMM kernels sme_[sdcz]gemm_kernel.c, private to the calling thread. */

#ifndef SME_GEMM_BUFFER_H
#define SME_GEMM_BUFFER_H

#include <stdlib.h>

#define SME_BUF_CACHE (4UL << 20) /* largest buffer kept per thread; larger ones come from blas_memory_alloc */

#if !defined(_WIN32)
#include <pthread.h>
static pthread_key_t sme_buf_key;
static pthread_once_t sme_buf_once = PTHREAD_ONCE_INIT;
static void sme_buf_free(void *p) { free(p); }
static void sme_buf_key_init(void) { pthread_key_create(&sme_buf_key, sme_buf_free); }

static void *sme_buf_cached(size_t bytes) {
  pthread_once(&sme_buf_once, sme_buf_key_init);
  size_t *p = (size_t *)pthread_getspecific(sme_buf_key);
  if (p == NULL || p[0] < bytes) {
    free(p);
    void *q = NULL;
    if (posix_memalign(&q, 256, SME_BUF_CACHE + 256)) q = NULL;
    p = (size_t *)q;
    if (p) p[0] = SME_BUF_CACHE;
    pthread_setspecific(sme_buf_key, p);
    if (p == NULL) return NULL;
  }
  return (char *)p + 256;
}
#else
static void *sme_buf_cached(size_t bytes) { (void)bytes; return NULL; }
#endif

/* At least bytes (at most BUFFER_SIZE) for the calling thread; *pool is the buffer to return with sme_buf_put. */
static void *sme_buf_get(size_t bytes, void **pool) {
  void *p = bytes <= SME_BUF_CACHE ? sme_buf_cached(bytes) : NULL;
  *pool = NULL;
  if (p == NULL) p = *pool = blas_memory_alloc(0);
  return p;
}

static void sme_buf_put(void *pool) {
  if (pool) blas_memory_free(pool);
}

#endif
