/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
/* Copyright 2023 The OpenBLAS Project                               */
/* All rights reserved.                                              */
/*                                                                   */
/* Redistribution and use in source and binary forms, with or        */
/* without modification, are permitted provided that the following   */
/* conditions are met:                                               */
/*                                                                   */
/*   1. Redistributions of source code must retain the above         */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer.                                                  */
/*                                                                   */
/*   2. Redistributions in binary form must reproduce the above      */
/*      copyright notice, this list of conditions and the following  */
/*      disclaimer in the documentation and/or other materials       */
/*      provided with the distribution.                              */
/*                                                                   */
/*    THIS  SOFTWARE IS PROVIDED  BY THE  UNIVERSITY OF  TEXAS AT    */
/*    AUSTIN  ``AS IS''  AND ANY  EXPRESS OR  IMPLIED WARRANTIES,    */
/*    INCLUDING, BUT  NOT LIMITED  TO, THE IMPLIED  WARRANTIES OF    */
/*    MERCHANTABILITY  AND FITNESS FOR  A PARTICULAR  PURPOSE ARE    */
/*    DISCLAIMED.  IN  NO EVENT SHALL THE UNIVERSITY  OF TEXAS AT    */
/*    AUSTIN OR CONTRIBUTORS BE  LIABLE FOR ANY DIRECT, INDIRECT,    */
/*    INCIDENTAL,  SPECIAL, EXEMPLARY,  OR  CONSEQUENTIAL DAMAGES    */
/*    (INCLUDING, BUT  NOT LIMITED TO,  PROCUREMENT OF SUBSTITUTE    */
/*    GOODS  OR  SERVICES; LOSS  OF  USE,  DATA,  OR PROFITS;  OR    */
/*    BUSINESS INTERRUPTION) HOWEVER CAUSED  AND ON ANY THEORY OF    */
/*    LIABILITY, WHETHER  IN CONTRACT, STRICT  LIABILITY, OR TORT    */
/*    (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY WAY OUT    */
/*    OF  THE  USE OF  THIS  SOFTWARE,  EVEN  IF ADVISED  OF  THE    */
/*    POSSIBILITY OF SUCH DAMAGE.                                    */
/*                                                                   */
/* The views and conclusions contained in the software and           */
/* documentation are those of the authors and should not be          */
/* interpreted as representing official policies, either expressed   */
/* or implied, of The University of Texas at Austin.                 */
/*********************************************************************/

#include "common.h"
#if (defined OS_LINUX || defined OS_ANDROID)
#include <asm/hwcap.h>
#include <sys/auxv.h>
#endif

extern gotoblas_t  gotoblas_ARMV8;
#ifdef DYNAMIC_LIST
#ifdef DYN_CORTEXA53
extern gotoblas_t  gotoblas_CORTEXA53;
#else
#define gotoblas_CORTEXA53 gotoblas_ARMV8
#endif
#ifdef DYN_CORTEXA57
extern gotoblas_t  gotoblas_CORTEXA57;
#else
#define gotoblas_CORTEXA57 gotoblas_ARMV8
#endif
#ifdef DYN_CORTEXA72
extern gotoblas_t  gotoblas_CORTEXA72;
#else
#define gotoblas_CORTEXA72 gotoblas_ARMV8
#endif
#ifdef DYN_CORTEXA73
extern gotoblas_t  gotoblas_CORTEXA73;
#else
#define gotoblas_CORTEXA73 gotoblas_ARMV8
#endif
#ifdef DYN_FALKOR
extern gotoblas_t  gotoblas_FALKOR;
#else
#define gotoblas_FALKOR gotoblas_ARMV8
#endif
#ifdef DYN_TSV110
extern gotoblas_t  gotoblas_TSV110;
#else
#define gotoblas_TSV110 gotoblas_ARMV8
#endif
#ifdef DYN_THUNDERX
extern gotoblas_t  gotoblas_THUNDERX;
#else
#define gotoblas_THUNDERX gotoblas_ARMV8
#endif
#ifdef DYN_THUNDERX2T99
extern gotoblas_t  gotoblas_THUNDERX2T99;
#else
#define gotoblas_THUNDERX2T99 gotoblas_ARMV8
#endif
#ifdef DYN_THUNDERX3T110
extern gotoblas_t  gotoblas_THUNDERX3T110;
#else
#define gotoblas_THUNDERX3T110 gotoblas_ARMV8
#endif
#ifdef DYN_EMAG8180
extern gotoblas_t  gotoblas_EMAG8180;
#else
#define gotoblas_EMAG8180 gotoblas_ARMV8
#endif
#ifdef DYN_NEOVERSEN1
extern gotoblas_t  gotoblas_NEOVERSEN1;
#else
#define gotoblas_NEOVERSEN1 gotoblas_ARMV8
#endif
#ifdef DYN_NEOVERSEV1
extern gotoblas_t  gotoblas_NEOVERSEV1;
#else
#define gotoblas_NEOVERSEV1 gotoblas_ARMV8
#endif
#ifdef DYN_NEOVERSEN2
extern gotoblas_t gotoblas_NEOVERSEN2;
#else
#define gotoblas_NEOVERSEN2 gotoblas_ARMV8
#endif
#ifdef DYN_ARMV8SVE
extern gotoblas_t gotoblas_ARMV8SVE;
#else
#define gotoblas_ARMV8SVE gotoblas_ARMV8
#endif
#ifdef DYN_CORTEX_A55
extern gotoblas_t  gotoblas_CORTEXA55;
#else
#define gotoblas_CORTEXA55 gotoblas_ARMV8
#endif
#else
extern gotoblas_t  gotoblas_CORTEXA53;
extern gotoblas_t  gotoblas_CORTEXA57;
extern gotoblas_t  gotoblas_CORTEXA72;
extern gotoblas_t  gotoblas_CORTEXA73;
extern gotoblas_t  gotoblas_FALKOR;
extern gotoblas_t  gotoblas_THUNDERX;
extern gotoblas_t  gotoblas_THUNDERX2T99;
extern gotoblas_t  gotoblas_TSV110;
extern gotoblas_t  gotoblas_EMAG8180;
extern gotoblas_t  gotoblas_NEOVERSEN1;
#ifndef NO_SVE
extern gotoblas_t  gotoblas_NEOVERSEV1;
extern gotoblas_t  gotoblas_NEOVERSEN2;
extern gotoblas_t  gotoblas_ARMV8SVE;
#else
#define gotoblas_NEOVERSEV1 gotoblas_ARMV8
#define gotoblas_NEOVERSEN2 gotoblas_ARMV8
#define gotoblas_ARMV8SVE   gotoblas_ARMV8
#endif
extern gotoblas_t  gotoblas_THUNDERX3T110;
extern gotoblas_t  gotoblas_CORTEXA55;
#endif

extern void openblas_warning(int verbose, const char * msg);
#define FALLBACK_VERBOSE 1
#define NEOVERSEN1_FALLBACK "OpenBLAS : Your OS does not support SVE instructions. OpenBLAS is using Neoverse N1 kernels as a fallback, which may give poorer performance.\n"

#define NUM_CORETYPES   16

/*
 * In case asm/hwcap.h is outdated on the build system, make sure
 * that HWCAP_CPUID is defined 
 */
#ifndef HWCAP_CPUID
#define HWCAP_CPUID (1 << 11)
#endif
#ifndef HWCAP_SVE
#define HWCAP_SVE (1 << 22)
#endif

#define get_cpu_ftr(id, var) ({					\
		__asm__ __volatile__ ("mrs %0, "#id : "=r" (var));		\
	})

static char *corename[] = {
  "armv8",
  "cortexa53",
  "cortexa57",
  "cortexa72",
  "cortexa73",
  "falkor",
  "thunderx",
  "thunderx2t99",
  "tsv110",
  "emag8180",
  "neoversen1",
  "neoversev1",
  "neoversen2",
  "thunderx3t110",
  "cortexa55",
  "armv8sve",
  "unknown"
};

char *gotoblas_corename(void) {
  if (gotoblas == &gotoblas_ARMV8)        return corename[ 0];
  if (gotoblas == &gotoblas_CORTEXA53)    return corename[ 1];
  if (gotoblas == &gotoblas_CORTEXA57)    return corename[ 2];
  if (gotoblas == &gotoblas_CORTEXA72)    return corename[ 3];
  if (gotoblas == &gotoblas_CORTEXA73)    return corename[ 4];
  if (gotoblas == &gotoblas_FALKOR)       return corename[ 5];
  if (gotoblas == &gotoblas_THUNDERX)     return corename[ 6];
  if (gotoblas == &gotoblas_THUNDERX2T99) return corename[ 7];
  if (gotoblas == &gotoblas_TSV110)       return corename[ 8];
  if (gotoblas == &gotoblas_EMAG8180)     return corename[ 9];
  if (gotoblas == &gotoblas_NEOVERSEN1)   return corename[10];
  if (gotoblas == &gotoblas_NEOVERSEV1)   return corename[11];
  if (gotoblas == &gotoblas_NEOVERSEN2)   return corename[12];
  if (gotoblas == &gotoblas_THUNDERX3T110) return corename[13];
  if (gotoblas == &gotoblas_CORTEXA55)    return corename[14];
  if (gotoblas == &gotoblas_ARMV8SVE)     return corename[15];
  return corename[NUM_CORETYPES];
}

static gotoblas_t *force_coretype(char *coretype) {
  int i ;
  int found = -1;
  char message[128];

  for ( i=0 ; i < NUM_CORETYPES; i++)
  {
    if (!strncasecmp(coretype, corename[i], 20))
    {
        found = i;
        break;
    }
  }

  switch (found)
  {
    case  0: return (&gotoblas_ARMV8);
    case  1: return (&gotoblas_CORTEXA53);
    case  2: return (&gotoblas_CORTEXA57);
    case  3: return (&gotoblas_CORTEXA72);
    case  4: return (&gotoblas_CORTEXA73);
    case  5: return (&gotoblas_FALKOR);
    case  6: return (&gotoblas_THUNDERX);
    case  7: return (&gotoblas_THUNDERX2T99);
    case  8: return (&gotoblas_TSV110);
    case  9: return (&gotoblas_EMAG8180);
    case 10: return (&gotoblas_NEOVERSEN1);
    case 11: return (&gotoblas_NEOVERSEV1);
    case 12: return (&gotoblas_NEOVERSEN2);
    case 13: return (&gotoblas_THUNDERX3T110);
    case 14: return (&gotoblas_CORTEXA55);
    case 15: return (&gotoblas_ARMV8SVE);
  }
  snprintf(message, 128, "Core not found: %s\n", coretype);
  openblas_warning(1, message);
  return NULL;
}

static gotoblas_t *get_coretype(void) {
  int implementer, variant, part, arch, revision, midr_el1;
  char coremsg[128];

#if (!defined OS_LINUX && !defined OS_ANDROID)
  return NULL;
#else

  if (!(getauxval(AT_HWCAP) & HWCAP_CPUID)) {
#ifdef __linux
        FILE *infile;
        char buffer[512], *p, *cpu_part = NULL, *cpu_implementer = NULL;
        p = (char *) NULL ;
	infile = fopen("/sys/devices/system/cpu/cpu0/regs/identification/midr_el1","r");
	if (!infile) return NULL;
	(void)fgets(buffer, sizeof(buffer), infile);
	midr_el1=strtoul(buffer,NULL,16);
	fclose(infile);
#else
    snprintf(coremsg, 128, "Kernel lacks cpuid feature support. Auto detection of core type failed !!!\n");
    openblas_warning(1, coremsg);
    return NULL;
#endif
  } else {
    get_cpu_ftr(MIDR_EL1, midr_el1);
  }
  /*
   * MIDR_EL1
   *
   * 31          24 23     20 19          16 15          4 3        0
   * -----------------------------------------------------------------
   * | Implementer | Variant | Architecture | Part Number | Revision |
   * -----------------------------------------------------------------
   */
  implementer = (midr_el1 >> 24) & 0xFF;
  part        = (midr_el1 >> 4)  & 0xFFF;

  switch(implementer)
  {
    case 0x41: // ARM
      switch (part)
      {
        case 0xd03: // Cortex A53
          return &gotoblas_CORTEXA53;
        case 0xd07: // Cortex A57
          return &gotoblas_CORTEXA57;
        case 0xd08: // Cortex A72
          return &gotoblas_CORTEXA72;
        case 0xd09: // Cortex A73
          return &gotoblas_CORTEXA73;
        case 0xd0c: // Neoverse N1
          return &gotoblas_NEOVERSEN1;
#ifndef NO_SVE
        case 0xd49:
	  if (!(getauxval(AT_HWCAP) & HWCAP_SVE)) {
	    openblas_warning(FALLBACK_VERBOSE, NEOVERSEN1_FALLBACK);  
	    return &gotoblas_NEOVERSEN1;
	  } else
            return &gotoblas_NEOVERSEN2;
	case 0xd40:
	  if (!(getauxval(AT_HWCAP) & HWCAP_SVE)) {
		  openblas_warning(FALLBACK_VERBOSE, NEOVERSEN1_FALLBACK);
	    return &gotoblas_NEOVERSEN1;
      	  }else
	    return &gotoblas_NEOVERSEV1;
#endif
	case 0xd05: // Cortex A55
	  return &gotoblas_CORTEXA55;
      }
      break;
    case 0x42: // Broadcom
      switch (part)
      {
        case 0x516: // Vulcan
          return &gotoblas_THUNDERX2T99;
      }
      break;
    case 0x43: // Cavium
      switch (part)
      {
        case 0x0a1: // ThunderX
          return &gotoblas_THUNDERX;
        case 0x0af: // ThunderX2
          return &gotoblas_THUNDERX2T99;
        case 0x0b8: // ThunderX3
          return &gotoblas_THUNDERX3T110;
      }
      break;
    case 0x48: // HiSilicon
      switch (part)
      {
        case 0xd01: // tsv110
          return &gotoblas_TSV110;
      }
      break;
    case 0x50: // Ampere
      switch (part)
      {
        case 0x000: // Skylark/EMAG8180
          return &gotoblas_EMAG8180;
      }
      break;
    case 0x51: // Qualcomm
      switch (part)
      {
        case 0xc00: // Falkor
          return &gotoblas_FALKOR;
      }
      break;
    default:
      snprintf(coremsg, 128, "Unknown CPU model - implementer %x part %x\n",implementer,part);
      openblas_warning(1, coremsg);
  }
#ifndef NO_SVE
  if ((getauxval(AT_HWCAP) & HWCAP_SVE)) {
    return &gotoblas_ARMV8SVE;
  }
#endif

  return NULL;
#endif
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
    snprintf(coremsg, 128, "Falling back to generic ARMV8 core\n");
    openblas_warning(1, coremsg);
    gotoblas = &gotoblas_ARMV8;
  }

  if (gotoblas && gotoblas->init) {
    strncpy(coren, gotoblas_corename(), 20);
    sprintf(coremsg, "Core: %s\n", coren);
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
