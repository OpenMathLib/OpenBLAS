/*********************************************************************/
/* Copyright 2009, 2010 The University of Texas at Austin.           */
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

#ifdef ARCH_X86
#define EXTERN extern
#else
#define EXTERN
#endif

EXTERN gotoblas_t  gotoblas_KATMAI;
EXTERN gotoblas_t  gotoblas_COPPERMINE;
EXTERN gotoblas_t  gotoblas_NORTHWOOD;
EXTERN gotoblas_t  gotoblas_BANIAS;
EXTERN gotoblas_t  gotoblas_ATHLON;

extern gotoblas_t  gotoblas_PRESCOTT;
extern gotoblas_t  gotoblas_ATOM;
extern gotoblas_t  gotoblas_NANO;
extern gotoblas_t  gotoblas_CORE2;
extern gotoblas_t  gotoblas_PENRYN;
extern gotoblas_t  gotoblas_DUNNINGTON;
extern gotoblas_t  gotoblas_NEHALEM;
extern gotoblas_t  gotoblas_OPTERON;
extern gotoblas_t  gotoblas_OPTERON_SSE3;
extern gotoblas_t  gotoblas_BARCELONA;

#define VENDOR_INTEL      1
#define VENDOR_AMD        2
#define VENDOR_CENTAUR    3
#define VENDOR_UNKNOWN   99

#define BITMASK(a, b, c) ((((a) >> (b)) & (c)))

static int get_vendor(void){
  int eax, ebx, ecx, edx;
  char vendor[13];

  cpuid(0, &eax, &ebx, &ecx, &edx);
  
  *(int *)(&vendor[0]) = ebx;
  *(int *)(&vendor[4]) = edx;
  *(int *)(&vendor[8]) = ecx;
  vendor[12] = (char)0;

  if (!strcmp(vendor, "GenuineIntel")) return VENDOR_INTEL;
  if (!strcmp(vendor, "AuthenticAMD")) return VENDOR_AMD;
  if (!strcmp(vendor, "CentaurHauls")) return VENDOR_CENTAUR;

  if ((eax == 0) || ((eax & 0x500) != 0)) return VENDOR_INTEL;

  return VENDOR_UNKNOWN;
}

static gotoblas_t *get_coretype(void){

  int eax, ebx, ecx, edx;
  int family, exfamily, model, vendor, exmodel;

  cpuid(1, &eax, &ebx, &ecx, &edx);

  family   = BITMASK(eax,  8, 0x0f);
  exfamily = BITMASK(eax, 20, 0xff);
  model    = BITMASK(eax,  4, 0x0f);
  exmodel  = BITMASK(eax, 16, 0x0f);

  vendor = get_vendor();

  if (vendor == VENDOR_INTEL){
    switch (family) {
    case 0x6:
      switch (exmodel) {
      case 0:
	if (model <= 0x7) return &gotoblas_KATMAI;
	if ((model == 0x8) || (model == 0xa) || (model == 0xb)) return &gotoblas_COPPERMINE;
	if ((model == 0x9) || (model == 0xd)) return &gotoblas_BANIAS;
	if (model == 14) return &gotoblas_BANIAS;
	if (model == 15) return &gotoblas_CORE2;
	return NULL;

      case 1:
	if (model == 6) return &gotoblas_CORE2;
	if (model == 7) return &gotoblas_PENRYN;
	if (model == 13) return &gotoblas_DUNNINGTON;
	if ((model == 10) || (model == 11) || (model == 14) || (model == 15)) return &gotoblas_NEHALEM;
	if (model == 12) return &gotoblas_ATOM;
	return NULL;

	  case 2:
		  //Intel Core (Clarkdale) / Core (Arrandale)
		  // Pentium (Clarkdale) / Pentium Mobile (Arrandale)
		  // Xeon (Clarkdale), 32nm
		  if (model ==  5) return &gotoblas_NEHALEM;
		  
		  //Intel Xeon Processor 5600 (Westmere-EP)
		  if (model == 12) return &gotoblas_NEHALEM;
		  return NULL;
      }
      case 0xf:
      if (model <= 0x2) return &gotoblas_NORTHWOOD;
      return &gotoblas_PRESCOTT;
    }
  }

  if (vendor == VENDOR_AMD){
    if (family <= 0xe) return &gotoblas_ATHLON;
    if (family == 0xf){
      if ((exfamily == 0) || (exfamily == 2)) {
	if (ecx & (1 <<  0)) return &gotoblas_OPTERON_SSE3; 
	else return &gotoblas_OPTERON;
      }  else {
	return &gotoblas_BARCELONA;
      }
    }
  }

  if (vendor == VENDOR_CENTAUR) {
    switch (family) {
    case 0x6:
      return &gotoblas_NANO;
      break;
    }
  }
  
  return NULL;
}

static char *corename[] = {
    "Unknown",
    "Katmai",
    "Coppermine",
    "Northwood",
    "Prescott",
    "Banias",
    "Atom",
    "Core2",
    "Penryn",
    "Dunnington",
    "Nehalem",
    "Athlon",
    "Opteron",
    "Opteron(SSE3)",
    "Barcelona",
    "Nano",
};

char *gotoblas_corename(void) {

  if (gotoblas == &gotoblas_KATMAI)       return corename[ 1];
  if (gotoblas == &gotoblas_COPPERMINE)   return corename[ 2];
  if (gotoblas == &gotoblas_NORTHWOOD)    return corename[ 3];
  if (gotoblas == &gotoblas_PRESCOTT)     return corename[ 4];
  if (gotoblas == &gotoblas_BANIAS)       return corename[ 5];
  if (gotoblas == &gotoblas_ATOM)         return corename[ 6];
  if (gotoblas == &gotoblas_CORE2)        return corename[ 7];
  if (gotoblas == &gotoblas_PENRYN)       return corename[ 8];
  if (gotoblas == &gotoblas_DUNNINGTON)   return corename[ 9];
  if (gotoblas == &gotoblas_NEHALEM)      return corename[10];
  if (gotoblas == &gotoblas_ATHLON)       return corename[11];
  if (gotoblas == &gotoblas_OPTERON_SSE3) return corename[12]; 
  if (gotoblas == &gotoblas_OPTERON)      return corename[13];
  if (gotoblas == &gotoblas_BARCELONA)    return corename[14];
  if (gotoblas == &gotoblas_NANO)         return corename[15];
  
  return corename[0];
}

void gotoblas_dynamic_init(void) {
  
  if (gotoblas) return;

  gotoblas = get_coretype();
  
#ifdef ARCH_X86
  if (gotoblas == NULL) gotoblas = &gotoblas_KATMAI;
#else
  if (gotoblas == NULL) gotoblas = &gotoblas_PRESCOTT;
#endif
  
  if (gotoblas && gotoblas -> init) {
    gotoblas -> init();
  } else {
    fprintf(stderr, "GotoBLAS : Architecture Initialization failed. No initialization function found.\n");
    exit(1);
  }
  
}

void gotoblas_dynamic_quit(void) {
  
  gotoblas = NULL;

}
