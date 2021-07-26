/*****************************************************************************
Copyright (c) 2011-2020, The OpenBLAS Project
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
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

**********************************************************************************/

#include <stdint.h>

#define CPU_UNKNOWN     0
#define CPU_LOONGSON3R5 1

#define LOONGARCH_CFG2  0x02
#define LOONGARCH_LASX  1<<7

static char *cpuname[] = {
  "UNKNOWN",
  "LOONGSON3R5"
};

int detect(void) {
    uint32_t reg = 0;

    __asm__ volatile (
        "cpucfg %0, %1 \n\t"
        : "+&r"(reg)
        : "r"(LOONGARCH_CFG2)
    );

    if (reg & LOONGARCH_LASX)
        return CPU_LOONGSON3R5;
    else
        return CPU_UNKNOWN;
}

char *get_corename(void) {
  return cpuname[detect()];
}

void get_architecture(void) {
  printf("LOONGARCH64");
}

void get_subarchitecture(void) {
  if (detect() == CPU_LOONGSON3R5) {
    printf("LOONGSON3R5");
  } else {
    printf("UNKNOWN");
  }
}

void get_subdirname(void) {
  printf("loongarch64");
}

void get_cpuconfig(void) {
  if (detect() == CPU_LOONGSON3R5) {
    printf("#define LOONGSON3R5\n");
    printf("#define L1_DATA_SIZE 65536\n");
    printf("#define L1_DATA_LINESIZE 64\n");
    printf("#define L2_SIZE 1048576\n");
    printf("#define L2_LINESIZE 64\n");
    printf("#define DTB_DEFAULT_ENTRIES 64\n");
    printf("#define DTB_SIZE 4096\n");
    printf("#define L2_ASSOCIATIVE 16\n");
  } else {
    printf("#define LOONGSON3R5\n");
    printf("#define L1_DATA_SIZE 65536\n");
    printf("#define L1_DATA_LINESIZE 64\n");
    printf("#define L2_SIZE 1048576\n");
    printf("#define L2_LINESIZE 64\n");
    printf("#define DTB_DEFAULT_ENTRIES 64\n");
    printf("#define DTB_SIZE 4096\n");
    printf("#define L2_ASSOCIATIVE 16\n");
  }
}

void get_libname(void){
  if (detect() == CPU_LOONGSON3R5) {
    printf("loongson3r5\n");
  } else {
    printf("loongarch64\n");
  }
}
