
#include "common.h"

extern gotoblas_t gotoblas_Z13;
extern gotoblas_t gotoblas_Z14;
//extern gotoblas_t gotoblas_Z15;
//#if (!defined C_GCC) || (GCC_VERSION >= 60000)
//extern gotoblas_t gotoblas_Z14;
//#endif

#define NUM_CORETYPES 4

extern void openblas_warning(int verbose, const char* msg);

static char* corename[] = {
	"unknown",
	"Z13",
	"Z14",
//	"Z15",
	"ZARCH_GENERIC",
};

char* gotoblas_corename(void) {
	if (gotoblas == &gotoblas_Z13)	return corename[1];
	if (gotoblas == &gotoblas_Z14)	return corename[2];
//	if (gotoblas == &gotoblas_Z15)	return corename[3];
//#if (!defined C_GCC) || (GCC_VERSION >= 60000)
//	if (gotoblas == &gotoblas_POWER9)	return corename[3];
//#endif
	return corename[0]; // try generic?
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

	if (strstr(p, "2964")) return &gotoblas_Z13;
	if (strstr(p, "2965")) return &gotoblas_Z13;
	if (strstr(p, "3906")) return &gotoblas_Z14;
	if (strstr(p, "3907")) return &gotoblas_Z14;
	if (strstr(p, "8561")) return &gotoblas_Z14;        // fallback z15 to z14
	if (strstr(p, "8562")) return &gotoblas_Z14;        // fallback z15 to z14

	return NULL; // should be ZARCH_GENERIC
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
	case  1: return (&gotoblas_Z13);
	case  2: return (&gotoblas_Z14);
//	case  3: return (&gotoblas_Z15);
//#if (!defined C_GCC) || (GCC_VERSION >= 60000)
//	case  3: return (&gotoblas_POWER9);
//#endif
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
		snprintf(coremsg, 128, "Falling back to Z14 core\n");
		openblas_warning(1, coremsg);
		gotoblas = &gotoblas_Z14;
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
