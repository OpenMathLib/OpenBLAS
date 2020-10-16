#include "common.h"
#include <stdbool.h>

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

extern gotoblas_t gotoblas_ZARCH_GENERIC;
#ifdef DYN_Z13
extern gotoblas_t gotoblas_Z13;
#endif
#ifdef DYN_Z14
extern gotoblas_t gotoblas_Z14;
#endif

#define NUM_CORETYPES 4

extern void openblas_warning(int verbose, const char* msg);

static char* corename[] = {
	"unknown",
	"Z13",
	"Z14",
	"ZARCH_GENERIC",
};

char* gotoblas_corename(void) {
#ifdef DYN_Z13
	if (gotoblas == &gotoblas_Z13)	return corename[1];
#endif
#ifdef DYN_Z14
	if (gotoblas == &gotoblas_Z14)	return corename[2];
#endif
	if (gotoblas == &gotoblas_ZARCH_GENERIC) return corename[3];

	return corename[0];
}

#ifndef HWCAP_S390_VXE
#define HWCAP_S390_VXE 8192
#endif

/**
 * Detect the fitting set of kernels by retrieving the CPU features supported by
 * OS from the auxiliary value AT_HWCAP and choosing the set of kernels
 * ("coretype") that exploits most of the features and can be compiled with the
 * available gcc version.
 * Note that we cannot use vector registers on a z13 or newer unless supported
 * by the OS kernel (which needs to handle them properly during context switch).
 */
static gotoblas_t* get_coretype(void) {

	unsigned long hwcap __attribute__((unused)) = get_hwcap();

#ifdef DYN_Z14
	// z14 and z15 systems: exploit Vector Facility (SIMD) and
	// Vector-Enhancements Facility 1 (float SIMD instructions), if present.
	if ((hwcap & HWCAP_S390_VX) && (hwcap & HWCAP_S390_VXE))
		return &gotoblas_Z14;
#endif

#ifdef DYN_Z13
	// z13: Vector Facility (SIMD for double)
	if (hwcap & HWCAP_S390_VX)
		return &gotoblas_Z13;
#endif

	// fallback in case of missing compiler support, systems before z13, or
	// when the OS does not advertise support for the Vector Facility (e.g.,
	// missing support in the OS kernel)
	return &gotoblas_ZARCH_GENERIC;
}

static gotoblas_t* force_coretype(char* coretype) {

	int i;
	int found = -1;
	char message[128];

	for (i = 0; i < NUM_CORETYPES; i++)
	{
		if (!strncasecmp(coretype, corename[i], 20))
		{
			found = i;
			break;
		}
	}

	if (found == 1) {
#ifdef DYN_Z13
		return &gotoblas_Z13;
#else
		openblas_warning(1, "Z13 support not compiled in");
		return NULL;
#endif
	} else if (found == 2) {
#ifdef DYN_Z14
		return &gotoblas_Z14;
#else
		openblas_warning(1, "Z14 support not compiled in");
		return NULL;
#endif
	} else if (found == 3) {
		return &gotoblas_ZARCH_GENERIC;
	}

	snprintf(message, 128, "Core not found: %s\n", coretype);
	openblas_warning(1, message);
	return NULL;
}

void gotoblas_dynamic_init(void) {

	char coremsg[128];
	char coren[22];
	char* p;


	if (gotoblas) return;

	p = getenv("OPENBLAS_CORETYPE");
	if (p)
	{
		gotoblas = force_coretype(p);
	}
	else
	{
		gotoblas = get_coretype();
	}

	if (gotoblas == NULL)
	{
		snprintf(coremsg, 128, "Failed to detect system, falling back to generic z support.\n");
		openblas_warning(1, coremsg);
		gotoblas = &gotoblas_ZARCH_GENERIC;
	}

	if (gotoblas && gotoblas->init) {
		strncpy(coren, gotoblas_corename(), 20);
		sprintf(coremsg, "Core: %s\n", coren);
		openblas_warning(2, coremsg);
		gotoblas->init();
	}
	else {
		openblas_warning(0, "OpenBLAS : Architecture Initialization failed. No initialization function found.\n");
		exit(1);
	}
}

void gotoblas_dynamic_quit(void) {
	gotoblas = NULL;
}
