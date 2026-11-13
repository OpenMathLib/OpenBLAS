/*
 * Rare finite-input range recovery for the AArch64 NRM2 kernels.
 * Distributed under the OpenBLAS BSD license (see LICENSE).
 */
#include "common.h"
#include <fenv.h>
#include <stdint.h>
#include <string.h>

/* A binary64 square needs at most 4198 bits in units of 2^-2150.
 * Fewer than 2^64 components (including complex components) need at most
 * 4262 bits. 67 limbs suffice, including carries. Two low guard bits make
 * squares of rounding midpoints integers too. No floating-point arithmetic
 * or precision-dependent long double is used to decide the result.
 */
#define LIMBS 67

static void add_limb(uint64_t *s, int i, uint64_t v)
{
    while (v) {
        uint64_t old = s[i];
        s[i++] += v;
        v = s[i - 1] < old;
    }
}

static void add_square(uint64_t *s, uint64_t m, int shift)
{
    /* m is at most 54 bits, including the midpoint's extra bit. */
    uint64_t a = (uint32_t)m, b = m >> 32;
    uint64_t low_product = a * a, cross = a * b;
    uint64_t lo = low_product + (cross << 33);
    uint64_t hi = b * b + (cross >> 31) + (lo < low_product);
    int i = shift / 64, r = shift % 64;
    add_limb(s, i, lo << r);
    add_limb(s, i + 1, hi << r);
    if (r) {
        add_limb(s, i + 1, lo >> (64 - r));
        add_limb(s, i + 2, hi >> (64 - r));
    }
}

static int compare(const uint64_t *a, const uint64_t *b)
{
    int i;
    for (i = LIMBS - 1; i >= 0; --i) {
        if (a[i] != b[i]) return a[i] > b[i] ? 1 : -1;
    }
    return 0;
}

static uint64_t range_significand(uint64_t bits, int fraction)
{
    uint64_t implicit = UINT64_C(1) << fraction;
    return (bits & (implicit - 1)) | (bits >> fraction ? implicit : 0);
}

static int square_shift(uint64_t bits, int fraction)
{
    int e = (int)(bits >> fraction);
    return 2 * (e ? e - 1 : 0) + 2;
}

static int compare_square(const uint64_t *sum, uint64_t bits, int fraction,
                          int midpoint)
{
    uint64_t square[LIMBS] = {0};
    uint64_t m = range_significand(bits, fraction);
    int shift = square_shift(bits, fraction);
    if (midpoint) {
        m = 2 * m + 1;
        shift -= 2;
    }
    add_square(square, m, shift);
    return compare(sum, square);
}

static uint64_t range_bits(BLASLONG n, const void *vx, BLASLONG incx,
                           int components, int single)
{
    uint64_t sum[LIMBS] = {0};
    const unsigned char *x = vx;
    int fraction = single ? 23 : 52;
    int bytes = single ? 4 : 8;
    uint64_t infinity = single ? UINT64_C(0x7f800000) : UINT64_C(0x7ff0000000000000);
    uint64_t sign = infinity | (infinity - 1);
    uint64_t lo = 0, hi = infinity - 1, result;
    BLASLONG i;
    int c, have_inf = 0, mode, exact, exceptions = 0;

    for (i = 0; i < n; ++i) {
        for (c = 0; c < components; ++c) {
            uint64_t bits;
            if (single) {
                uint32_t u;
                memcpy(&u, x + c * bytes, sizeof(u));
                bits = u;
            } else {
                memcpy(&bits, x + c * bytes, sizeof(bits));
            }
            bits &= sign;
            /* Normally excluded by the assembly's finite-scale guard.
             * A later NaN can nevertheless leave a finite final scale.
             */
            if (bits > infinity) return bits | (UINT64_C(1) << (fraction - 1));
            if (bits == infinity) have_inf = 1;
            else add_square(sum, range_significand(bits, fraction), square_shift(bits, fraction));
        }
        /* Avoid forming a pointer outside the object after the last item. */
        if (i + 1 < n) x += incx * components * bytes;
    }
    if (have_inf) return infinity;

    /* Find the greatest finite representable value with square <= sum.
     * Positive IEEE encodings have the same ordering as their values.
     */
    while (lo < hi) {
        uint64_t mid = lo + (hi - lo + 1) / 2;
        if (compare_square(sum, mid, fraction, 0) >= 0) lo = mid;
        else hi = mid - 1;
    }
    exact = compare_square(sum, lo, fraction, 0) == 0;
    result = lo;
    mode = fegetround();
    if (!exact) {
        exceptions = FE_INEXACT;
        if (mode == FE_UPWARD) ++result;
        else if (mode == FE_TONEAREST) {
            int cmp = compare_square(sum, lo, fraction, 1);
            if (cmp > 0 || (cmp == 0 && (lo & 1))) ++result;
        }
        /* Directed rounding can overflow while returning max finite.
         * infinity's encoding is interpreted here as the next *finite*
         * binade, i.e. 2^1024 or 2^128, for this squared comparison only.
         */
        if (result == infinity ||
            (lo == infinity - 1 && compare_square(sum, infinity, fraction, 0) >= 0))
            exceptions |= FE_OVERFLOW;
        if (result < (UINT64_C(1) << fraction)) exceptions |= FE_UNDERFLOW;
        feraiseexcept(exceptions);
    }
    return result;
}

#if defined(__GNUC__) && !defined(_WIN32)
__attribute__((visibility("hidden")))
#endif
double openblas_dnrm2_range(BLASLONG n, const double *x, BLASLONG incx, int components)
{
    uint64_t bits = range_bits(n, x, incx, components, 0);
    double result;
    memcpy(&result, &bits, sizeof(result));
    return result;
}

#if defined(__GNUC__) && !defined(_WIN32)
__attribute__((visibility("hidden")))
#endif
float openblas_snrm2_range(BLASLONG n, const float *x, BLASLONG incx, int components)
{
    uint32_t bits = (uint32_t)range_bits(n, x, incx, components, 1);
    float result;
    memcpy(&result, &bits, sizeof(result));
    return result;
}
