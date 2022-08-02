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
#include <math.h>

/*
 * <https://loongson.github.io/LoongArch-Documentation/
 *  Loongson-3A5000-usermanual-EN.html#la464-processor-core>
 */
#define CPU_GENERIC     0
#define CPU_LA464       1
#define CPU_LA264       2

#define LOONGARCH_CFG0      0x00
#define LOONGARCH_CFG2      0x02
#define LOONGARCH_CFG10     0x10
#define LOONGARCH_CFG11     0x11
#define LOONGARCH_CFG12     0x12
#define LOONGARCH_CFG13     0x13
#define LOONGARCH_CFG14     0x14
#define LASX_MASK           1<<7
#define LSX_MASK            1<<6
#define PRID_SERIES_MASK    0xf000
#define PRID_SERIES_LA464   0xc000
#define PRID_SERIES_LA264   0xa000

#define CACHE_INFO_L1_IU    0
#define CACHE_INFO_L1_D     1
#define CACHE_INFO_L2_IU    2
#define CACHE_INFO_L3_IU    3
#define L1_IU_PRESENT_MASK  0x0001
#define L1_IU_UNITY_MASK    0x0002
#define L1_D_PRESENT_MASK   0x0004
#define L2_IU_PRESENT_MASK  0x0008
#define L2_IU_UNITY_MASK    0x0010
#define L2_D_PRESENT_MASK   0x0080
#define L3_IU_PRESENT_MASK  0x0400
#define L3_IU_UNITY_MASK    0x0800
#define L3_D_PRESENT_MASK   0x4000
#define CACHE_WAY_MINUS_1_MASK      0x0000ffff
#define CACHE_INDEX_LOG2_MASK       0x00ff0000
#define CACHE_LINESIZE_LOG2_MASK    0x7f000000

typedef struct {
  int size;
  int associative;
  int linesize;
  int unify;
} cache_info_t;

static char *cpuname[] = {
  "LA64_GENERIC",
  "LA464",
  "LA264"
};

static char *cpuname_lower[] = {
  "la64_generic",
  "la464",
  "la264"
};

static void get_cacheinfo(int type, cache_info_t *cacheinfo) {
  cache_info_t cache_info;
  memset(&cache_info, 0, sizeof(cache_info));
  uint32_t reg_10 = 0;
  __asm__ volatile (
    "cpucfg %0, %1 \n\t"
    : "+&r"(reg_10)
    : "r"(LOONGARCH_CFG10)
  );

  switch (type) {
    case CACHE_INFO_L1_IU:
      if (reg_10 & L1_IU_PRESENT_MASK) {
        uint32_t reg_11 = 0;
        cache_info.unify = reg_10 & L1_IU_UNITY_MASK;
        __asm__ volatile (
          "cpucfg %0, %1 \n\t"
          : "+&r"(reg_11)
          : "r"(LOONGARCH_CFG11)
        );
        cache_info.associative  = (reg_11 & CACHE_WAY_MINUS_1_MASK) + 1;
        cache_info.linesize = pow(2, (reg_11 & CACHE_LINESIZE_LOG2_MASK) >> 24);
        cache_info.size = cache_info.associative * cache_info.linesize *
                          pow(2, (reg_11 & CACHE_INDEX_LOG2_MASK) >> 16);
      }
    break;

    case CACHE_INFO_L1_D:
      if (reg_10 & L1_D_PRESENT_MASK) {
        uint32_t reg_12 = 0;
        cache_info.unify = reg_10 & L1_IU_UNITY_MASK;
        __asm__ volatile (
          "cpucfg %0, %1 \n\t"
          : "+&r"(reg_12)
          : "r"(LOONGARCH_CFG12)
        );
        cache_info.associative  = (reg_12 & CACHE_WAY_MINUS_1_MASK) + 1;
        cache_info.linesize = pow(2, (reg_12 & CACHE_LINESIZE_LOG2_MASK) >> 24);
        cache_info.size = cache_info.associative * cache_info.linesize *
                          pow(2, (reg_12 & CACHE_INDEX_LOG2_MASK) >> 16);
      }
    break;

    case CACHE_INFO_L2_IU:
      if (reg_10 & L2_IU_PRESENT_MASK) {
        uint32_t reg_13 = 0;
        cache_info.unify = reg_10 & L2_IU_UNITY_MASK;
        __asm__ volatile (
          "cpucfg %0, %1 \n\t"
          : "+&r"(reg_13)
          : "r"(LOONGARCH_CFG13)
        );
        cache_info.associative  = (reg_13 & CACHE_WAY_MINUS_1_MASK) + 1;
        cache_info.linesize = pow(2, (reg_13 & CACHE_LINESIZE_LOG2_MASK) >> 24);
        cache_info.size = cache_info.associative * cache_info.linesize *
                          pow(2, (reg_13 & CACHE_INDEX_LOG2_MASK) >> 16);
      }
    break;

    case CACHE_INFO_L3_IU:
      if (reg_10 & L3_IU_PRESENT_MASK) {
        uint32_t reg_14 = 0;
        cache_info.unify = reg_10 & L3_IU_UNITY_MASK;
        __asm__ volatile (
          "cpucfg %0, %1 \n\t"
          : "+&r"(reg_14)
          : "r"(LOONGARCH_CFG14)
        );
        cache_info.associative  = (reg_14 & CACHE_WAY_MINUS_1_MASK) + 1;
        cache_info.linesize = pow(2, (reg_14 & CACHE_LINESIZE_LOG2_MASK) >> 24);
        cache_info.size = cache_info.associative * cache_info.linesize *
                          pow(2, (reg_14 & CACHE_INDEX_LOG2_MASK) >> 16);
      }
    break;

    default:
    break;
  }
  *cacheinfo = cache_info;
}

static void get_cpucount(uint32_t *count) {
#ifdef __linux
  uint32_t num = 0;
  FILE *f = fopen("/proc/cpuinfo", "r");
  if (!f) return;
  char buf[200];
  while (fgets(buf, sizeof(buf), f))
  {
    if (!strncmp("processor", buf, 9))
      num ++;
  }
  fclose(f);
  *count = num;
#endif
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

static uint32_t get_prid() {
  uint32_t reg = 0;
  __asm__ volatile (
    "cpucfg %0, %1 \n\t"
    : "+&r"(reg)
    : "r"(LOONGARCH_CFG0)
  );
  return reg;
}

int detect(void) {
#ifdef __linux
  uint32_t prid = get_prid();
  switch (prid & PRID_SERIES_MASK) {
    case (PRID_SERIES_LA464):
      if (support_lasx())
        return CPU_LA464;
      else
        return CPU_GENERIC;
    break;

    case (PRID_SERIES_LA264):
      if (support_lsx())
        return CPU_LA264;
      else
        return CPU_GENERIC;
    break;

    default:
      return CPU_GENERIC;
  }
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
  cache_info_t info;
  uint32_t num_cores = 0;

  printf("#define %s\n", cpuname[detect()]);

  get_cacheinfo(CACHE_INFO_L1_IU, &info);
  if (info.size > 0) {
      printf("#define L1_CODE_SIZE %d\n", info.size);
      printf("#define L1_CODE_ASSOCIATIVE %d\n", info.associative);
      printf("#define L1_CODE_LINESIZE %d\n", info.linesize);
  }
  get_cacheinfo(CACHE_INFO_L1_D, &info);
  if (info.size > 0) {
      printf("#define L1_DATA_SIZE %d\n", info.size);
      printf("#define L1_DATA_ASSOCIATIVE %d\n", info.associative);
      printf("#define L1_DATA_LINESIZE %d\n", info.linesize);
  }
  get_cacheinfo(CACHE_INFO_L2_IU, &info);
  if (info.size > 0) {
    if (info.unify) {
      printf("#define L2_SIZE %d\n", info.size);
      printf("#define L2_ASSOCIATIVE %d\n", info.associative);
      printf("#define L2_LINESIZE %d\n", info.linesize);
    } else {
      printf("#define L2_CODE_SIZE %d\n", info.size);
      printf("#define L2_CODE_ASSOCIATIVE %d\n", info.associative);
      printf("#define L2_CODE_LINESIZE %d\n", info.linesize);
    }
  }
  get_cacheinfo(CACHE_INFO_L3_IU, &info);
  if (info.size > 0) {
    if (info.unify) {
      printf("#define L3_SIZE %d\n", info.size);
      printf("#define L3_ASSOCIATIVE %d\n", info.associative);
      printf("#define L3_LINESIZE %d\n", info.linesize);
    } else {
      printf("#define L3_CODE_SIZE %d\n", info.size);
      printf("#define L3_CODE_ASSOCIATIVE %d\n", info.associative);
      printf("#define L3_CODE_LINESIZE %d\n", info.linesize);
    }
  }
  get_cpucount(&num_cores);
  if (num_cores)
    printf("#define NUM_CORES %d\n", num_cores);
  printf("#define DTB_DEFAULT_ENTRIES 64\n");
  printf("#define DTB_SIZE 4096\n");
}

void get_libname(void){
  int d = detect();
  printf("%s", cpuname_lower[d]);
}
