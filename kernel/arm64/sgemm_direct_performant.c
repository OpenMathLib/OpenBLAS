#include "common.h"
#include "sme2_gemm_detect.h"
/* helper for the direct sgemm code adapted from Arjan van der Ven's x86_64 version */

int CNAME(BLASLONG M, BLASLONG N, BLASLONG K)
{
	if (M < 3 || N <= 0 || K <= 0)
		return 0;

	unsigned long long mnk = (unsigned long long)M * (unsigned long long)N * (unsigned long long)K;
#ifdef HAVE_SME2_GEMM
	/* the SME2 kernel behind SME_SGEMM_KERNEL is faster from about 64^3 (Apple M4) */
	if (mnk >= 64ULL * 64ULL * 64ULL && s2_usable())
		return 0;
#endif
	/* benchmark performance on M4 peaks around 512 and crosses the graph of the NEON SGEMM at about 3100  */
	if (mnk >= 3100ULL * 3100ULL * 3100ULL)
		return 0;

	return 1;
}
