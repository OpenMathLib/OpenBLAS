/*******************************************************************************
Copyright (c) 2022, The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*******************************************************************************/

#include "common.h"

extern gotoblas_t  gotoblas_LA464;
extern gotoblas_t  gotoblas_LA264;
extern gotoblas_t  gotoblas_LA64_GENERIC;

extern void openblas_warning(int verbose, const char * msg);

#define NUM_CORETYPES    3

static char *corename[] = {
  "la464",
  "la264",
  "la64_generic",
  "unknown"
};

char *gotoblas_corename(void) {
  if (gotoblas == &gotoblas_LA464)        return corename[0];
  if (gotoblas == &gotoblas_LA264)        return corename[1];
  if (gotoblas == &gotoblas_LA64_GENERIC) return corename[2];
  return corename[NUM_CORETYPES];
}

static gotoblas_t *force_coretype(char *coretype) {
  int i;
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
    case  0: return (&gotoblas_LA464);
    case  1: return (&gotoblas_LA264);
    case  2: return (&gotoblas_LA64_GENERIC);
  }
  snprintf(message, 128, "Core not found: %s\n", coretype);
  openblas_warning(1, message);
  return NULL;
}

#define LASX_MASK           1<<7
#define LSX_MASK            1<<6
#define LOONGARCH_CFG0      0x00
#define LOONGARCH_CFG2      0x02
#define PRID_SERIES_MASK    0xf000
#define PRID_SERIES_LA464   0xc000
#define PRID_SERIES_LA264   0xa000

static uint32_t get_prid() {
  uint32_t reg = 0;
  __asm__ volatile (
    "cpucfg %0, %1 \n\t"
    : "+&r"(reg)
    : "r"(LOONGARCH_CFG0)
  );
  return reg;
}

static int support_lasx() {
  uint32_t reg = 0;
  __asm__ volatile (
    "cpucfg %0, %1 \n\t"
    : "+&r"(reg)
    : "r"(LOONGARCH_CFG2)
  );

  if (reg & LASX_MASK)
    return 1;
  return 0;
}

static int support_lsx() {
  uint32_t reg = 0;
  __asm__ volatile (
    "cpucfg %0, %1 \n\t"
    : "+&r"(reg)
    : "r"(LOONGARCH_CFG2)
  );

  if (reg & LSX_MASK)
    return 1;
  return 0;
}

static gotoblas_t *get_coretype(void) {
#ifdef __linux
  uint32_t prid = get_prid();
  switch (prid & PRID_SERIES_MASK) {
    case (PRID_SERIES_LA464):
      if (support_lasx())
        return &gotoblas_LA464;
      else
        return &gotoblas_LA64_GENERIC;
    break;

    case (PRID_SERIES_LA264):
      if (support_lsx())
        return &gotoblas_LA264;
      else
        return &gotoblas_LA64_GENERIC;
    break;

    default:
      return &gotoblas_LA64_GENERIC;
  }
#endif
  return &gotoblas_LA64_GENERIC;
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
