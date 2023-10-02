
#include "common.h"

extern gotoblas_t gotoblas_POWER6;
extern gotoblas_t gotoblas_POWER8;
#if (!defined __GNUC__) || ( __GNUC__ >= 6)
extern gotoblas_t gotoblas_POWER9;
#endif
#ifdef HAVE_P10_SUPPORT
extern gotoblas_t gotoblas_POWER10;
#endif

extern void openblas_warning(int verbose, const char *msg);

static char *corename[] = {
    "unknown",
    "POWER6",
    "POWER8",
    "POWER9",
    "POWER10"
};

#define NUM_CORETYPES 5

char *gotoblas_corename(void) {
    if (gotoblas == &gotoblas_POWER6)    return corename[1];
    if (gotoblas == &gotoblas_POWER8)    return corename[2];
#if (!defined __GNUC__) || ( __GNUC__ >= 6)
    if (gotoblas == &gotoblas_POWER9)    return corename[3];
#endif
#ifdef HAVE_P10_SUPPORT
    if (gotoblas == &gotoblas_POWER10)   return corename[4];
#endif
    return corename[0];
}

#ifdef _AIX
#include <sys/systemcfg.h>

#define CPU_UNKNOWN  0
#define CPU_POWER6   6
#define CPU_POWER7   7
#define CPU_POWER8   8
#define CPU_POWER9   9
#define CPU_POWER10 10

static int cpuid(void)
{
    int arch = _system_configuration.implementation;
#ifdef POWER_6
    if (arch == POWER_6) return CPU_POWER6;
#endif
#ifdef POWER_7
    else if (arch == POWER_7) return CPU_POWER7;
#endif
#ifdef POWER_8
    else if (arch == POWER_8) return CPU_POWER8;
#endif
#ifdef POWER_9
    else if (arch == POWER_9) return CPU_POWER9;
#endif
#ifdef POWER_10
    else if (arch == POWER_10) return CPU_POWER10;
#endif
    return CPU_UNKNOWN;
}

#ifndef __BUILTIN_CPU_SUPPORTS__
static int __builtin_cpu_supports(char* arg)
{
    static int ipinfo = -1;
    if (ipinfo < 0) {
        ipinfo = cpuid();
    }
    if (ipinfo >= CPU_POWER10) {
        if (!strcmp(arg, "power10")) return 1;
    }
    if (ipinfo >= CPU_POWER9) {
        if (!strcmp(arg, "power9")) return 1;
    }
    if (ipinfo >= CPU_POWER8) {
        if (!strcmp(arg, "power8")) return 1;
    }
    if (ipinfo >= CPU_POWER6) {
        if (!strcmp(arg, "power6")) return 1;
    }
    return 0;
}
#endif

static gotoblas_t *get_coretype(void) {

    if (__builtin_cpu_supports("power6"))
        return &gotoblas_POWER6;
    if (__builtin_cpu_supports("power8"))
        return &gotoblas_POWER8;
#if (!defined __GNUC__) || ( __GNUC__ >= 6)
    if (__builtin_cpu_supports("power9"))
        return &gotoblas_POWER9;
#endif
#ifdef HAVE_P10_SUPPORT
#ifdef _AIX
    if (__builtin_cpu_supports("power10"))
#else
    if (__builtin_cpu_supports("arch_3_1") && __builtin_cpu_supports("mma"))
#endif
        return &gotoblas_POWER10;
#endif
    /* Fall back to the POWER9 implementation if the toolchain is too old or the MMA feature is not set */
#if (!defined __GNUC__) || ( __GNUC__ < 11) || (__GNUC__ == 10 && __GNUC_MINOR__ < 2)
    if (__builtin_cpu_supports("power10"))
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
#if (!defined __GNUC__) || ( __GNUC__ >= 6)
    case  3: return (&gotoblas_POWER9);
#endif
#ifdef HAVE_P10_SUPPORT
    case  4: return (&gotoblas_POWER10);
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
