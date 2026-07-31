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

#include <stddef.h>
#include "common.h"

/* Cooperative cancellation of in-flight OpenBLAS operations.
 *
 * Every thread owns a pointer-sized generation slot in thread-local
 * storage:
 *
 *   - openblas_cancel_token() returns the address of the calling thread's
 *     slot.  The address is stable for the lifetime of the thread.
 *   - Instrumented compute drivers advance the slot to a fresh even value
 *     at operation entry (on the issuing thread) and poll it at block
 *     granularity: any change - the cancel bit or a later generation -
 *     makes the operation abandon its remaining work and return early.
 *   - openblas_cancel(token, loaded_token) requests cancellation of the
 *     operation that was in flight when the caller loaded *token: it
 *     atomically sets the cancel bit iff the slot still holds
 *     loaded_token.  A racing or stale request (the operation completed;
 *     a new one may have started) either fails the compare-exchange or
 *     dirties an already-dead generation - both harmless.  May be called
 *     from any thread.
 *
 * A cancelled operation returns quickly but leaves its output buffers in
 * an unspecified, partially-updated state; the caller is responsible for
 * discarding the result.  The library itself remains in a consistent
 * state and can service further calls.
 *
 * These symbols are deliberately exported without any SYMBOLPREFIX /
 * SYMBOLSUFFIX decoration, like the other openblas_* utility functions
 * that are not per-precision entry points.
 */

/* Storage for the per-thread generation slot.  Compiler-emitted TLS is
 * used where the target supports it (clang reports per-target support
 * via __has_feature(tls); e.g. 32-bit iOS with an old deployment minimum
 * lacks it).  Where it does not, the slot falls back to a
 * pthread_getspecific() lookup of a lazily heap-allocated slot: slightly
 * more expensive at operation entry, but the compute drivers cache the
 * slot's address, so the polling fast path is identical. */
#if defined(_MSC_VER) && !defined(__clang__)
#define OPENBLAS_CANCEL_TLS __declspec(thread)
#elif defined(__clang__)
#if __has_feature(tls)
#define OPENBLAS_CANCEL_TLS __thread
#endif
#elif defined(__GNUC__) || defined(__SUNPRO_C) || defined(__xlC__)
#define OPENBLAS_CANCEL_TLS __thread
#endif

#ifdef OPENBLAS_CANCEL_TLS

/* This thread's generation slot.  Bit 0 is the cancel flag; the remaining
 * bits count operations issued by this thread. */
static OPENBLAS_CANCEL_TLS size_t openblas_cancel_slot = 0;

/* Address of the calling thread's generation slot. */
static size_t *openblas_cancel_slot_addr(void) {
  return (size_t *)&openblas_cancel_slot;
}

#else

#include <pthread.h>
#include <stdlib.h>

static pthread_key_t openblas_cancel_key;
static pthread_once_t openblas_cancel_key_once = PTHREAD_ONCE_INIT;
static int openblas_cancel_key_ok = 0;

static void openblas_cancel_key_init(void) {
  openblas_cancel_key_ok = (pthread_key_create(&openblas_cancel_key, free) == 0);
}

/* Address of the calling thread's generation slot, allocated on first
 * use and freed by the key destructor at thread exit.  Returns NULL if
 * per-thread state cannot be set up; cancellation is then unavailable
 * on this thread (all other entry points accept a NULL slot). */
static size_t *openblas_cancel_slot_addr(void) {
  size_t *slot;
  pthread_once(&openblas_cancel_key_once, openblas_cancel_key_init);
  if (!openblas_cancel_key_ok) return NULL;
  slot = (size_t *)pthread_getspecific(openblas_cancel_key);
  if (slot == NULL) {
    slot = (size_t *)calloc(1, sizeof(size_t));
    if (slot == NULL) return NULL;
    if (pthread_setspecific(openblas_cancel_key, slot) != 0) {
      free(slot);
      return NULL;
    }
  }
  return slot;
}

#endif

#if defined(__GNUC__) || defined(__clang__)
#define OPENBLAS_CANCEL_LOAD(PTR) __atomic_load_n((PTR), __ATOMIC_RELAXED)
#define OPENBLAS_CANCEL_STORE(PTR, VAL) __atomic_store_n((PTR), (VAL), __ATOMIC_RELEASE)
#define OPENBLAS_CANCEL_CAS(PTR, EXPP, VAL) \
  __atomic_compare_exchange_n((PTR), (EXPP), (VAL), 0, __ATOMIC_RELEASE, __ATOMIC_RELAXED)
#else
/* Fallback for compilers without the GNU atomic builtins.  The slot is
 * only advanced by its owning thread; the cancel request degrades to a
 * check-then-store, which can at worst dirty a dead generation. */
#define OPENBLAS_CANCEL_LOAD(PTR) (*(volatile size_t *)(PTR))
#define OPENBLAS_CANCEL_STORE(PTR, VAL) do { WMB; *(volatile size_t *)(PTR) = (VAL); } while (0)
static int openblas_cancel_cas_fallback(size_t *ptr, size_t *expp, size_t val) {
  if (*(volatile size_t *)ptr != *expp) return 0;
  WMB;
  *(volatile size_t *)ptr = val;
  return 1;
}
#define OPENBLAS_CANCEL_CAS(PTR, EXPP, VAL) openblas_cancel_cas_fallback((PTR), (EXPP), (VAL))
#endif

/* Address of the calling thread's generation slot (NULL if unavailable). */
size_t *openblas_cancel_token(void) {
  return openblas_cancel_slot_addr();
}

/* Request cancellation of the operation that was in flight on the token's
 * thread when the caller loaded LOADED_TOKEN from it.  May be called from
 * any thread. */
void openblas_cancel(size_t *token, size_t loaded_token) {
  if (token == NULL) return;
  OPENBLAS_CANCEL_CAS(token, &loaded_token, loaded_token | (size_t)1);
}

/* --- internal helpers (not part of the public API) --------------------- */

/* Begin an instrumented operation on the calling thread: advance the slot
 * past any stale cancel bit to a fresh even generation and return it.
 * Returns 0 if the thread has no slot (paired with a NULL slot from
 * openblas_cancel_self, so polling never reports cancellation). */
size_t openblas_cancel_begin(void) {
  size_t *slot = openblas_cancel_slot_addr();
  size_t gen;
  if (slot == NULL) return 0;
  gen = (*slot | (size_t)1) + 1;
  OPENBLAS_CANCEL_STORE(slot, gen);
  return gen;
}

/* The calling thread's slot address (same as openblas_cancel_token, for
 * internal use without going through the exported symbol). */
size_t *openblas_cancel_self(void) {
  return openblas_cancel_slot_addr();
}

/* Cheap poll usable from compute drivers: nonzero once GEN is no longer
 * the live value of SLOT (cancelled, or superseded).  SLOT may be NULL. */
int openblas_cancel_poll(size_t *slot, size_t gen) {
  if (slot == NULL) return 0;
  return OPENBLAS_CANCEL_LOAD(slot) != gen;
}
