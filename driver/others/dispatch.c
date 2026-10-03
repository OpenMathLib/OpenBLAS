/*****************************************************************************
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
   3. Neither the name of the OpenBLAS project nor the names of its contributors may
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
*****************************************************************************/

/* The arrays of per-core kernel tables that OPENBLAS_DISPATCH(group)
   indexes (see common_param.h), one per kernel group. */

#include "common.h"

#define DISPATCH_EXTERN(core, group) \
  extern const openblas_##group##_dispatch_t openblas_##group##_dispatch_##core;
#define DISPATCH_ENTRY(core, group) \
  [OPENBLAS_CORE_##core] = &openblas_##group##_dispatch_##core,

/* openblas_<group>_dispatch[OPENBLAS_CORE_<CORE>] points to
   openblas_<group>_dispatch_<CORE>, which kernel/setparam-ref.c defines for
   that core. */
#define DEFINE_DISPATCH_TABLE(group) \
  OPENBLAS_CORE_LIST(DISPATCH_EXTERN, group) \
  const openblas_##group##_dispatch_t *const openblas_##group##_dispatch[OPENBLAS_NUM_CORES] = { \
    OPENBLAS_CORE_LIST(DISPATCH_ENTRY, group) \
  };
