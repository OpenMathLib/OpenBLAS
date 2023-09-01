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
#include <sys/auxv.h>

/*  If LASX extension instructions supported,
 *  using core LOONGSON3R5
 *  If only LSX extension instructions supported,
 *  using core LOONGSON2K1000
 *  If neither LASX nor LSX extension instructions supported,
 *  using core LOONGSONGENERIC (As far as I know, there is no such
 *  CPU yet)
 */

#define CPU_GENERIC        0
#define CPU_LOONGSON3R5    1
#define CPU_LOONGSON2K1000 2

#define LA_HWCAP_LSX    (1<<4)
#define LA_HWCAP_LASX   (1<<5)

static char *cpuname[] = {
  "LOONGSONGENERIC",
  "LOONGSON3R5",
  "LOONGSON2K1000"
};

static char *cpuname_lower[] = {
  "loongsongeneric",
  "loongson3r5",
  "loongson2k1000"
};

int detect(void) {
#ifdef __linux
  int flag  = (int)getauxval(AT_HWCAP);

  if (flag & LA_HWCAP_LASX)
    return CPU_LOONGSON3R5;
  else if (flag & LA_HWCAP_LSX)
    return CPU_LOONGSON2K1000;
  else
    return CPU_GENERIC;
#endif
  return CPU_GENERIC;
}

char *get_corename(void) {
  return cpuname[detect()];
}

void get_architecture(void) {
  printf("LOONGARCH64");
}

void get_subarchitecture(void) {
  int d = detect();
  printf("%s", cpuname[d]);
}

void get_subdirname(void) {
  printf("loongarch64");
}

void get_cpuconfig(void) {
  int d = detect();
  switch (d) {
    case CPU_LOONGSON3R5:
      printf("#define LOONGSON3R5\n");
      printf("#define L1_DATA_SIZE 65536\n");
      printf("#define L1_DATA_LINESIZE 64\n");
      printf("#define L2_SIZE 1048576\n");
      printf("#define L2_LINESIZE 64\n");
      printf("#define DTB_DEFAULT_ENTRIES 64\n");
      printf("#define DTB_SIZE 4096\n");
      printf("#define L2_ASSOCIATIVE 16\n");
    break;

    case CPU_LOONGSON2K1000:
      printf("#define LOONGSON2K1000\n");
      printf("#define L1_DATA_SIZE 65536\n");
      printf("#define L1_DATA_LINESIZE 64\n");
      printf("#define L2_SIZE 262144\n");
      printf("#define L2_LINESIZE 64\n");
      printf("#define DTB_DEFAULT_ENTRIES 64\n");
      printf("#define DTB_SIZE 4096\n");
      printf("#define L2_ASSOCIATIVE 16\n");
    break;

    default:
      printf("#define LOONGSONGENERIC\n");
      printf("#define L1_DATA_SIZE 65536\n");
      printf("#define L1_DATA_LINESIZE 64\n");
      printf("#define L2_SIZE 262144\n");
      printf("#define L2_LINESIZE 64\n");
      printf("#define DTB_DEFAULT_ENTRIES 64\n");
      printf("#define DTB_SIZE 4096\n");
      printf("#define L2_ASSOCIATIVE 16\n");
    break;
  }
}

void get_libname(void){
  int d = detect();
  printf("%s", cpuname_lower[d]);
}
