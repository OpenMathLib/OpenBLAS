/*****************************************************************************
Copyright (c) 2026, The OpenBLAS Project
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

#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include "openblas_utest.h"

#if defined(DYNAMIC_ARCH) && defined(__x86_64__) && defined(__linux__)

/* Builds that drop some cores map their names to another core's table,
   so the selected name is only compared where every core is built. */
#if !defined(DYNAMIC_LIST) && !defined(NO_AVX) && !defined(NO_AVX2) && !defined(NO_AVX512)
#define CHECK_SELECTED_NAME
#endif

static void host_cpuid(int leaf, int sub, int *eax, int *ebx, int *ecx, int *edx)
{
    __asm__ __volatile__("cpuid"
                         : "=a"(*eax), "=b"(*ebx), "=c"(*ecx), "=d"(*edx)
                         : "a"(leaf), "c"(sub));
}

/* The OS has enabled the register state bits in XCR0 selected by mask. */
static int host_os_saves(int mask)
{
    int eax, ecx, edx, ebx;

    host_cpuid(1, 0, &eax, &ebx, &ecx, &edx);
    if ((ecx & (1 << 27)) == 0)
        return 0;
    __asm__ __volatile__(".byte 0x0f, 0x01, 0xd0"
                         : "=a"(eax), "=d"(edx) : "c"(0));
    return (eax & mask) == mask;
}

/* EBX of CPUID leaf 7, or 0 where the cpu has no such leaf. */
static int host_leaf7_ebx(void)
{
    int eax, ebx, ecx, edx;

    host_cpuid(0, 0, &eax, &ebx, &ecx, &edx);
    if (eax < 7)
        return 0;
    host_cpuid(7, 0, &eax, &ebx, &ecx, &edx);
    return ebx;
}

/* Forcing a core runs its initialisation, which is compiled for that
   core and can use instructions the host lacks, so names the host
   cannot run are skipped. The setparam objects of the three AVX512
   cores contain BMI1 instructions (andn), the AVX2 and AVX ones none,
   so the BMI1 flag is only needed for the first group. The feature
   bits are read with cpuid here rather than with
   __builtin_cpu_supports, which older compilers reject for these
   names. */
static int host_can_run(const char *name)
{
    int ebx7 = host_leaf7_ebx();
    int avx, avx2;

    avx = host_os_saves(0x6);
    avx2 = avx && (ebx7 & (1 << 5)) != 0;
    if (!strcmp(name, "SkylakeX") || !strcmp(name, "Cooperlake") ||
        !strcmp(name, "SapphireRapids"))
        return avx2 && host_os_saves(0xe6) && (ebx7 & (1 << 16)) != 0 &&
               (ebx7 & (1U << 31)) != 0 && (ebx7 & (1 << 3)) != 0;
    if (!strcmp(name, "Haswell") || !strcmp(name, "Zen"))
        return avx2;
    if (!strcmp(name, "Sandybridge") || !strcmp(name, "Bulldozer") ||
        !strcmp(name, "Piledriver") || !strcmp(name, "Steamroller") ||
        !strcmp(name, "Excavator"))
        return avx;
    return 1;
}

/* OPENBLAS_CORETYPE is only read when the library is loaded, so every
   name is tried in a new process that runs this binary with a suite
   filter matching no test. */
CTEST(coretype, force_by_name)
{
    static const char *names[] = {
        "Prescott", "Core2", "Nehalem", "Barcelona", "Sandybridge",
        "Bulldozer", "Piledriver", "Steamroller", "Excavator", "Haswell",
        "Zen", "SkylakeX", "Cooperlake", "SapphireRapids"
    };
    char self[4096], cmd[4352], line[256], selected[64];
    ssize_t len;
    size_t i;
    int rejected, status;
    FILE *p;

    len = readlink("/proc/self/exe", self, sizeof(self) - 1);
    if (len <= 0) {
        fprintf(stderr, "cannot read /proc/self/exe, skipping\n");
        return;
    }
    self[len] = '\0';

    for (i = 0; i < sizeof(names) / sizeof(names[0]); i++) {
        if (!host_can_run(names[i])) {
            fprintf(stderr, "skipping %s, not supported by this cpu\n", names[i]);
            continue;
        }
        snprintf(cmd, sizeof(cmd),
                 "OPENBLAS_VERBOSE=2 OPENBLAS_CORETYPE=%s '%s' coretype_none 2>&1",
                 names[i], self);
        p = popen(cmd, "r");
        ASSERT_NOT_NULL(p);
        rejected = 0;
        selected[0] = '\0';
        while (fgets(line, sizeof(line), p) != NULL) {
            if (strstr(line, "Core not found") != NULL) rejected = 1;
            if (strncmp(line, "Core: ", 6) == 0)
                sscanf(line + 6, "%63s", selected);
        }
        status = pclose(p);
        if (rejected) fprintf(stderr, "OPENBLAS_CORETYPE=%s not accepted\n", names[i]);
        ASSERT_EQUAL(0, status);
        ASSERT_EQUAL(0, rejected);
        ASSERT_NOT_EQUAL(0, selected[0]);
#ifdef CHECK_SELECTED_NAME
        if (strcasecmp(selected, names[i]) != 0)
            fprintf(stderr, "OPENBLAS_CORETYPE=%s selected %s\n", names[i], selected);
        ASSERT_EQUAL(0, strcasecmp(selected, names[i]));
#endif
    }
}
#endif
