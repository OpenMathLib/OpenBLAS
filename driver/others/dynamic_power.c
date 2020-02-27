
#include "common.h"

extern gotoblas_t gotoblas_POWER6;
extern gotoblas_t gotoblas_POWER8;
#if (!defined C_GCC) || (GCC_VERSION >= 60000)
extern gotoblas_t gotoblas_POWER9;
#endif

extern void openblas_warning(int verbose, const char *msg);

static char *corename[] = {
	"unknown",
	"POWER6",
	"POWER8",
	"POWER9"
};

#define NUM_CORETYPES 4

char *gotoblas_corename(void) {
	if (gotoblas == &gotoblas_POWER6)	return corename[1];
	if (gotoblas == &gotoblas_POWER8)	return corename[2];
#if (!defined C_GCC) || (GCC_VERSION >= 60000)
	if (gotoblas == &gotoblas_POWER9)	return corename[3];
#endif
	return corename[0];
}

static gotoblas_t *get_coretype(void) {

	if (__builtin_cpu_is("power6") || __builtin_cpu_is("power6x"))
		return &gotoblas_POWER6;
	if (__builtin_cpu_is("power8"))
		return &gotoblas_POWER8;
#if (!defined C_GCC) || (GCC_VERSION >= 60000)
	if (__builtin_cpu_is("power9"))
		return &gotoblas_POWER9;
#endif
	return NULL;
}

static gotoblas_t *force_coretype(char * coretype) {

	int i ;
	int found = -1;
	char message[128];

	for ( i = 0 ; i < NUM_CORETYPES; i++)
	{
		if (!strncasecmp(coretype, corename[i], 20))
		{
			found = i;
			break;
		}
	}

	switch (found)
	{
	case  1: return (&gotoblas_POWER6);
	case  2: return (&gotoblas_POWER8);
#if (!defined C_GCC) || (GCC_VERSION >= 60000)
	case  3: return (&gotoblas_POWER9);
#endif
	default: return NULL;
	}
	snprintf(message, 128, "Core not found: %s\n", coretype);
	openblas_warning(1, message);
}

void gotoblas_dynamic_init(void) {

	char coremsg[128];
	char coren[22];
	char *p;


	if (gotoblas) return;

	p = getenv("OPENBLAS_CORETYPE");
	if ( p )
	{
		gotoblas = force_coretype(p);
	}
	else
	{
		gotoblas = get_coretype();
	}

	if (gotoblas == NULL)
	{
		snprintf(coremsg, 128, "Falling back to POWER8 core\n");
		openblas_warning(1, coremsg);
		gotoblas = &gotoblas_POWER8;
	}

	if (gotoblas && gotoblas -> init) {
		strncpy(coren,gotoblas_corename(),20);
		sprintf(coremsg, "Core: %s\n",coren);
		openblas_warning(2, coremsg);
		gotoblas -> init();
	} else {
		openblas_warning(0, "OpenBLAS : Architecture Initialization failed. No initialization function found.\n");
		exit(1);
	}
}

void gotoblas_dynamic_quit(void) {
	gotoblas = NULL;
}
