/*****************************************************************************
Copyright (c) 2020, The OpenBLAS Project
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

#include <sys/wait.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include "common.h"

extern gotoblas_t  gotoblas_LOONGSON3R3;
extern gotoblas_t  gotoblas_LOONGSON3R4;

extern void openblas_warning(int verbose, const char * msg);

#define NUM_CORETYPES    2

static char *corename[] = {
  "loongson3r3",
  "loongson3r4",
  "UNKNOWN"
};

char *gotoblas_corename(void) {
  if (gotoblas == &gotoblas_LOONGSON3R3)    return corename[0];
  if (gotoblas == &gotoblas_LOONGSON3R4)    return corename[1];
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
    case  0: return (&gotoblas_LOONGSON3R3);
    case  1: return (&gotoblas_LOONGSON3R4);
  }
  snprintf(message, 128, "Core not found: %s\n", coretype);
  openblas_warning(1, message);
  return NULL;
}

#define MMI_MASK    0x00000010
#define MSA_MASK    0x00000020

int fd[2];
int support_cpucfg;

static void handler(int signum)
{
    close(fd[1]);
    exit(1);
}

/* Brief :  Function to check if cpucfg supported on loongson
 * Return:  1   supported
 *          0   not supported
 */
static int cpucfg_test(void) {
    pid_t pid;
    int status = 0;

    support_cpucfg = 0;
    pipe(fd);
    pid = fork();
    if (pid == 0) { /* Subprocess */
        struct sigaction act;
        close(fd[0]);
        /* Set signal action for SIGILL. */
        act.sa_handler = handler;
        sigaction(SIGILL,&act,NULL);

        /* Execute cpucfg in subprocess. */
        __asm__ volatile(
            ".insn              \n\t"
            ".word (0xc8080118) \n\t"
            :::
        );
        support_cpucfg = 1;
        write(fd[1],&support_cpucfg,sizeof(support_cpucfg));
        close(fd[1]);
        exit(0);
    } else if (pid > 0){ /* Parent process*/
        close(fd[1]);
        if ((waitpid(pid,&status,0) <= 0) ||
            (read(fd[0],&support_cpucfg,sizeof(support_cpucfg)) <= 0))
            support_cpucfg = 0;
        close(fd[0]);
    } else {
        support_cpucfg = 0;
    }

    return support_cpucfg;
}

static gotoblas_t *get_coretype_from_cpucfg(void) {
    int flag = 0;
    __asm__ volatile(
        ".insn                     \n\t"
        "dli    $8,    0x01        \n\t"
        ".word (0xc9084918)        \n\t"
        "usw    $9,    0x00(%0)    \n\t"
        :
        : "r"(&flag)
        : "memory"
    );
    if (flag & MSA_MASK)
        return (&gotoblas_LOONGSON3R4);
    if (flag & MMI_MASK)
        return (&gotoblas_LOONGSON3R3);
    return NULL;
}

static gotoblas_t *get_coretype_from_cpuinfo(void) {
#ifdef linux
  FILE *infile;
  char buffer[512], *p;

  p = (char *)NULL;
  //Check model name for Loongson3
  infile = fopen("/proc/cpuinfo", "r");
  while (fgets(buffer, sizeof(buffer), infile)){
    if (!strncmp("model name", buffer, 10)){
      p = strchr(buffer, ':') + 2;
      break;
    }
  }
  fclose(infile);
  if(p != NULL){
   if (strstr(p, "Loongson-3A3000") || strstr(p, "Loongson-3B3000"))
     return (&gotoblas_LOONGSON3R3);
   else if(strstr(p, "Loongson-3A4000") || strstr(p, "Loongson-3B4000"))
     return (&gotoblas_LOONGSON3R4);
   else
     return NULL;
  }
#endif
    return NULL;
}

static gotoblas_t *get_coretype(void) {
    int ret = 0;

    ret = cpucfg_test();
    if (ret == 1)
        return get_coretype_from_cpucfg();
    else
        return get_coretype_from_cpuinfo();
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
    snprintf(coremsg, 128, "Falling back to loongson3r3 core\n");
    openblas_warning(1, coremsg);
    gotoblas = &gotoblas_LOONGSON3R3;
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
