#include "common.h"
#include <stdbool.h>

// Gate kernels for z13 and z14 on gcc version
#if (__GNUC__ == 5 && __GNUC_MINOR__ >= 2) || __GNUC__ >= 6 ||           \
    /* RHEL 7 since 7.3: */                                              \
    (__GNUC__ == 4 && __GNUC_MINOR__ == 8 && __GNUC_PATCHLEVEL__ == 5 && \
     __GNUC_RH_RELEASE__ >= 11)
#define HAVE_Z13_SUPPORT
#endif

#if __GNUC__ >= 7
#define HAVE_Z14_SUPPORT
#endif

extern gotoblas_t gotoblas_ZARCH_GENERIC;
#ifdef HAVE_Z13_SUPPORT
extern gotoblas_t gotoblas_Z13;
#endif
#ifdef HAVE_Z14_SUPPORT
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
#ifdef HAVE_Z13_SUPPORT
	if (gotoblas == &gotoblas_Z13)	return corename[1];
#endif
#ifdef HAVE_Z14_SUPPORT
	if (gotoblas == &gotoblas_Z14)	return corename[2];
#endif
	if (gotoblas == &gotoblas_ZARCH_GENERIC) return corename[3];

	return corename[0];
}

// __builtin_cpu_is is not supported by zarch
static gotoblas_t* get_coretype(void) {
	FILE* infile;
	char buffer[512], * p;

	p = (char*)NULL;
	infile = fopen("/proc/sysinfo", "r");
	while (fgets(buffer, sizeof(buffer), infile)) {
		if (!strncmp("Type", buffer, 4)) {
			p = strchr(buffer, ':') + 2;
#if 0
			fprintf(stderr, "%s\n", p);
#endif
			break;
		}
	}

	fclose(infile);

#ifdef HAVE_Z13_SUPPORT
	if (strstr(p, "2964") || strstr(p, "2965")) return &gotoblas_Z13;
#endif

	// Z14 and Z15 systems
	if (strstr(p, "3906") || strstr(p, "3907") || strstr(p, "8561") ||
	    strstr(p, "8562"))
#ifdef HAVE_Z14_SUPPORT
		return &gotoblas_Z14;
#else
		return &gotoblas_Z13;
#endif

	// unknown system or compiler too old? use generic code for z architecture
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

	switch (found)
	{
#ifdef HAVE_Z13_SUPPORT
	case  1: return (&gotoblas_Z13);
#endif
#ifdef HAVE_Z14_SUPPORT
	case  2: return (&gotoblas_Z14);
#endif
	case  3: return (&gotoblas_ZARCH_GENERIC);
	default: return NULL;
	}
	snprintf(message, 128, "Core not found: %s\n", coretype);
	openblas_warning(1, message);
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
