/*****************************************************************************
 * Replaceable allocator for LAPACKE-internal allocations.
 *
 * The LAPACKE_malloc / LAPACKE_free macros default to the functions below,
 * which forward to a run-time replaceable allocator (malloc / free until
 * changed via LAPACKE_set_alloc). Like the nancheck flag, the allocator is
 * shared between all APIs, so this file must be compiled exactly once,
 * without an API suffix.
 *****************************************************************************/

#include <stdlib.h>

#include "lapacke.h"

static void *(*lapacke_malloc_fn)(size_t) = malloc;
static void (*lapacke_free_fn)(void *) = free;

/**
 * \brief Allocate through the installed allocator (malloc by default).
 *
 * Default of the LAPACKE_malloc macro.
 *
 * \param[in] size Allocation size in bytes.
 * \return The allocation, or NULL on failure.
 */
void *LAPACKE_malloc_proxy(size_t size)
{
    return lapacke_malloc_fn(size);
}

/**
 * \brief Release through the installed deallocator (free by default).
 *
 * Default of the LAPACKE_free macro.
 *
 * \param[in] ptr Allocation obtained through LAPACKE_malloc_proxy.
 */
void LAPACKE_free_proxy(void *ptr)
{
    lapacke_free_fn(ptr);
}

/**
 * \brief Install a custom allocator for LAPACKE-internal allocations.
 *
 * Subsequent LAPACKE-internal allocations (workspaces, row-major
 * transposition buffers) are made with malloc_fn and released with
 * free_fn. The two must form a matching pair: memory allocated with
 * malloc_fn is released through free_fn, so pairing an allocator with an
 * incompatible deallocator is undefined behavior. The setter mutates
 * global state shared by all threads and is not thread-safe; call it
 * during setup, like LAPACKE_set_nancheck.
 *
 * \param[in] malloc_fn Replacement for malloc; NULL restores malloc.
 * \param[in] free_fn   Replacement for free; NULL restores free.
 */
void LAPACKE_set_alloc(void *(*malloc_fn)(size_t), void (*free_fn)(void *))
{
    lapacke_malloc_fn = (malloc_fn != NULL) ? malloc_fn : malloc;
    lapacke_free_fn = (free_fn != NULL) ? free_fn : free;
}
