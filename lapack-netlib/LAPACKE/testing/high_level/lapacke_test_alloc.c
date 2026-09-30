/*****************************************************************************
 * Allocation failure injection and leak tracking for the dedicated
 * high-level LAPACKE tests. See lapacke_test.h.
 *
 * The runner main installs the allocator below for LAPACKE-internal
 * allocations via LAPACKE_set_alloc, so all tests run against the stock
 * LAPACKE library. When a failure is scheduled, the (countdown+1)-th
 * allocation fails once and the schedule clears itself. Every live
 * allocation is tracked so that lapacke_test_check_leaks can verify that
 * LAPACKE released all its buffers, in particular on error paths (nonzero
 * info, injected allocation failures).
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "lapacke_test.h"

#define LAPACKE_TEST_MAX_LIVE_ALLOCS 64

static long mallocs_until_failure = -1;
static void *live_allocs[LAPACKE_TEST_MAX_LIVE_ALLOCS];
static size_t live_alloc_count = 0;

/**
 * \brief The counting/failing test allocator.
 *
 * Returns NULL (once) when the scheduled countdown reaches zero; otherwise
 * forwards to malloc and records the allocation in the live table.
 *
 * \param[in] size Allocation size in bytes.
 * \return The allocation, or NULL on (injected) failure.
 */
void *lapacke_test_malloc(size_t size)
{
    if (mallocs_until_failure >= 0) {
        if (mallocs_until_failure == 0) {
            mallocs_until_failure = -1;
            return NULL;
        }
        mallocs_until_failure--;
    }

    void *ptr = malloc(size);
    if (ptr != NULL) {
        if (live_alloc_count < LAPACKE_TEST_MAX_LIVE_ALLOCS) {
            live_allocs[live_alloc_count] = ptr;
            live_alloc_count++;
        } else {
            lapacke_test_checks++;
            lapacke_test_failures++;
            printf("FAIL allocation tracking table overflow\n");
        }
    }
    return ptr;
}

/**
 * \brief The deallocator matching lapacke_test_malloc.
 *
 * Removes the pointer from the live table and frees it; freeing a pointer
 * that is not a live allocation (double free, foreign pointer) fails the
 * test.
 *
 * \param[in] ptr The allocation to release; NULL is ignored.
 */
void lapacke_test_free(void *ptr)
{
    if (ptr == NULL) return;

    for (size_t i = 0; i < live_alloc_count; i++) {
        if (live_allocs[i] == ptr) {
            live_alloc_count--;
            live_allocs[i] = live_allocs[live_alloc_count];
            free(ptr);
            return;
        }
    }
    lapacke_test_checks++;
    lapacke_test_failures++;
    printf("FAIL free of a pointer that is not a live allocation"
           " (double or foreign free?)\n");
}

/**
 * \brief Schedule the (countdown+1)-th LAPACKE allocation to fail once.
 *
 * The schedule clears itself when the failure fires, so at most one
 * allocation fails per schedule. Used to test the LAPACK_WORK_MEMORY_ERROR
 * / LAPACK_TRANSPOSE_MEMORY_ERROR paths.
 *
 * \param[in] countdown Number of allocations to let through before failing.
 */
void lapacke_test_schedule_malloc_failure(long countdown)
{
    mallocs_until_failure = countdown;
}

/**
 * \brief Record one check that a call consumed the scheduled countdown exactly.
 *
 * Call directly after a call with a failure scheduled at its exact
 * expected number of allocations: more allocations would have failed the
 * call, fewer leave the countdown positive. Clears the schedule
 * afterwards.
 *
 * \param[in] name Name of the check.
 */
void lapacke_test_check_alloc_count(const char *name)
{
    lapacke_test_checks++;
    if (mallocs_until_failure == 0) {
        printf("ok   %s: allocations match the scheduled countdown\n", name);
    } else if (mallocs_until_failure < 0) {
        lapacke_test_failures++;
        printf("FAIL %s: more allocations than expected (the scheduled"
               " failure fired)\n",
               name);
    } else {
        lapacke_test_failures++;
        printf("FAIL %s: %ld fewer allocations than expected\n", name,
               mallocs_until_failure);
    }
    mallocs_until_failure = -1;
}

/**
 * \brief Record one check that every tracked allocation has been freed.
 *
 * LAPACKE must release its buffers on error paths as well; the test runner
 * calls this after each test. Resets the allocation tracking, releases any
 * leaked allocations and fails an additional check if a scheduled
 * allocation failure never triggered.
 *
 * \param[in] name Name of the test the runner just executed.
 */
void lapacke_test_check_leaks(const char *name)
{
    lapacke_test_checks++;
    if (live_alloc_count != 0) {
        lapacke_test_failures++;
        printf("FAIL %s: %d allocation(s) not freed\n", name,
               (int)live_alloc_count);

        for (size_t i = 0; i < live_alloc_count; i++) {
            free(live_allocs[i]);
        }
        live_alloc_count = 0;
    } else {
        printf("ok   %s: no leaked allocations\n", name);
    }

    if (mallocs_until_failure >= 0) {
        lapacke_test_checks++;
        lapacke_test_failures++;
        printf("FAIL %s: allocation failure scheduled but never triggered\n",
               name);
        mallocs_until_failure = -1;
    }
}
