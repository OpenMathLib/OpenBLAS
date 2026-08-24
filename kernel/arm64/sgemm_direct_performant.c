#include "common.h"
/* helper for the direct sgemm code adapted from Arjan van der Ven's x86_64 version */

int CNAME(BLASLONG M, BLASLONG N, BLASLONG K)
{
	if (M < 3 || N <= 0 || K <= 0)
		return 0;

	unsigned long long mnk = (unsigned long long)M * (unsigned long long)N * (unsigned long long)K;
	/* benchmark performance on M4 peaks around 512 and crosses the graph of the NEON SGEMM at about 3100  */
	if (mnk >= 3100ULL * 3100ULL * 3100ULL)
		return 0;

	return 1;
}
