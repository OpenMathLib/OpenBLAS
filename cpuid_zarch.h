#include <stdlib.h>

#define CPU_GENERIC     0
#define CPU_Z13         1
#define CPU_Z14         2
#define CPU_Z15         3

static char *cpuname[] = {
  "ZARCH_GENERIC",
  "Z13",
  "Z14",
  "Z15"
};

static char *cpuname_lower[] = {
  "zarch_generic",
  "z13",
  "z14",
  "z15"
};

// Guard the use of getauxval() on glibc version >= 2.16
#ifdef __GLIBC__
#include <features.h>
#if __GLIBC_PREREQ(2, 16)
#include <sys/auxv.h>
#define HAVE_GETAUXVAL 1

static unsigned long get_hwcap(void)
{
	unsigned long hwcap = getauxval(AT_HWCAP);
	char *maskenv;

	// honor requests for not using specific CPU features in LD_HWCAP_MASK
	maskenv = getenv("LD_HWCAP_MASK");
	if (maskenv)
		hwcap &= strtoul(maskenv, NULL, 0);

	return hwcap;
	// note that a missing auxval is interpreted as no capabilities
	// available, which is safe.
}

#else // __GLIBC_PREREQ(2, 16)
#warn "Cannot detect SIMD support in Z13 or newer architectures since glibc is older than 2.16"

static unsigned long get_hwcap(void) {
	// treat missing support for getauxval() as no capabilities available,
	// which is safe.
	return 0;
}
#endif // __GLIBC_PREREQ(2, 16)
#endif // __GLIBC

static int detect(void)
{
	unsigned long hwcap = get_hwcap();

	if ((hwcap & HWCAP_S390_VX) && (hwcap & HWCAP_S390_VXE))
		return CPU_Z14;

	if (hwcap & HWCAP_S390_VX)
		return CPU_Z13;

	return CPU_GENERIC;
}

