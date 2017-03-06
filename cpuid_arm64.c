/**************************************************************************
  Copyright (c) 2013, The OpenBLAS Project
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
  *****************************************************************************/

#include <string.h>

#define CPU_UNKNOWN     	0
#define CPU_ARMV8       	1
#define CPU_CORTEXA57       	2
#define CPU_VULCAN       	3
#define CPU_THUNDERX    	4
#define CPU_THUNDERX2T99   	5

static char *cpuname[] = {
  "UNKNOWN",
  "ARMV8" ,
  "CORTEXA57",
  "VULCAN",
  "THUNDERX",
  "THUNDERX2T99"
};

static char *cpuname_lower[] = {
  "unknown",
  "armv8" ,
  "cortexa57",
  "vulcan",
  "thunderx",
  "thunderx2t99"
};

int get_feature(char *search)
{

#ifdef linux
	FILE *infile;
  	char buffer[2048], *p,*t;
  	p = (char *) NULL ;

  	infile = fopen("/proc/cpuinfo", "r");

	while (fgets(buffer, sizeof(buffer), infile))
	{

		if (!strncmp("Features", buffer, 8))
		{
			p = strchr(buffer, ':') + 2;
			break;
		}
	}

	fclose(infile);


	if( p == NULL ) return 0;

	t = strtok(p," ");
	while( t = strtok(NULL," "))
	{
		if (!strcmp(t, search))   { return(1); }
	}

#endif
	return(0);
}


int detect(void)
{

#ifdef linux

	FILE *infile;
	char buffer[512], *p, *cpu_part = NULL, *cpu_implementer = NULL;
	p = (char *) NULL ;

	infile = fopen("/proc/cpuinfo", "r");
	while (fgets(buffer, sizeof(buffer), infile)) {
		if ((cpu_part != NULL) && (cpu_implementer != NULL)) {
			break;
		}

		if ((cpu_part == NULL) && !strncmp("CPU part", buffer, 8)) {
			cpu_part = strchr(buffer, ':') + 2;
			cpu_part = strdup(cpu_part);
		} else if ((cpu_implementer == NULL) && !strncmp("CPU implementer", buffer, 15)) {
			cpu_implementer = strchr(buffer, ':') + 2;
			cpu_implementer = strdup(cpu_implementer);
		}
	}

	fclose(infile);
	if(cpu_part != NULL && cpu_implementer != NULL) {
		if (strstr(cpu_part, "0xd07") && strstr(cpu_implementer, "0x41"))
			return CPU_CORTEXA57;
		else if (strstr(cpu_part, "0x516") && strstr(cpu_implementer, "0x42"))
			return CPU_VULCAN;
		else if (strstr(cpu_part, "0x0a1") && strstr(cpu_implementer, "0x43"))
			return CPU_THUNDERX;
		else if (strstr(cpu_part, "0xFFF") && strstr(cpu_implementer, "0x43")) /* TODO */
			return CPU_THUNDERX2T99;
	}

	p = (char *) NULL ;
	infile = fopen("/proc/cpuinfo", "r");
	while (fgets(buffer, sizeof(buffer), infile))
	{

		if ((!strncmp("model name", buffer, 10)) || (!strncmp("Processor", buffer, 9)) ||
		    (!strncmp("CPU architecture", buffer, 16)))
		{
			p = strchr(buffer, ':') + 2;
			break;
      		}
  	}

  	fclose(infile);

  	if(p != NULL)
	{

		if (strstr(p, "AArch64"))
		{
			return CPU_ARMV8;

		}


	}
#endif

	return CPU_UNKNOWN;
}

char *get_corename(void)
{
	return cpuname[detect()];
}

void get_architecture(void)
{
	printf("ARM64");
}

void get_subarchitecture(void)
{
	int d = detect();
	printf("%s", cpuname[d]);
}

void get_subdirname(void)
{
	printf("arm64");
}

void get_cpuconfig(void)
{

	int d = detect();
	switch (d)
	{

		case CPU_ARMV8:
    			printf("#define ARMV8\n");
    			printf("#define L1_DATA_SIZE 32768\n");
    			printf("#define L1_DATA_LINESIZE 64\n");
    			printf("#define L2_SIZE 262144\n");
    			printf("#define L2_LINESIZE 64\n");
    			printf("#define DTB_DEFAULT_ENTRIES 64\n");
    			printf("#define DTB_SIZE 4096\n");
    			printf("#define L2_ASSOCIATIVE 4\n");
			break;

		case CPU_VULCAN:
			printf("#define VULCAN                        \n");
			printf("#define HAVE_VFP                      \n");
			printf("#define HAVE_VFPV3                    \n");
			printf("#define HAVE_NEON                     \n");
			printf("#define HAVE_VFPV4                    \n");
			printf("#define L1_CODE_SIZE         32768    \n");
			printf("#define L1_CODE_LINESIZE     64       \n");
			printf("#define L1_CODE_ASSOCIATIVE  8        \n");
			printf("#define L1_DATA_SIZE         32768    \n");
			printf("#define L1_DATA_LINESIZE     64       \n");
			printf("#define L1_DATA_ASSOCIATIVE  8        \n");
			printf("#define L2_SIZE              262144   \n");
			printf("#define L2_LINESIZE          64       \n");
			printf("#define L2_ASSOCIATIVE       8        \n");
			printf("#define L3_SIZE              33554432 \n");
			printf("#define L3_LINESIZE          64       \n");
			printf("#define L3_ASSOCIATIVE       32       \n");
			printf("#define DTB_DEFAULT_ENTRIES  64       \n");
			printf("#define DTB_SIZE             4096     \n");
			break;

		case CPU_CORTEXA57:
			printf("#define CORTEXA57\n");
			printf("#define HAVE_VFP\n");
			printf("#define HAVE_VFPV3\n");
			printf("#define HAVE_NEON\n");
			printf("#define HAVE_VFPV4\n");
			printf("#define L1_CODE_SIZE 49152\n");
			printf("#define L1_CODE_LINESIZE 64\n");
			printf("#define L1_CODE_ASSOCIATIVE 3\n");
			printf("#define L1_DATA_SIZE 32768\n");
			printf("#define L1_DATA_LINESIZE 64\n");
			printf("#define L1_DATA_ASSOCIATIVE 2\n");
			printf("#define L2_SIZE 2097152\n");
			printf("#define L2_LINESIZE 64\n");
			printf("#define L2_ASSOCIATIVE 16\n");
			printf("#define DTB_DEFAULT_ENTRIES 64\n");
			printf("#define DTB_SIZE 4096\n");
			break;

		case CPU_THUNDERX:
			printf("#define ARMV8\n");
			printf("#define THUNDERX\n");
			printf("#define L1_DATA_SIZE 32768\n");
			printf("#define L1_DATA_LINESIZE 128\n");
			printf("#define L2_SIZE 16777216\n");
			printf("#define L2_LINESIZE 128\n");
			printf("#define DTB_DEFAULT_ENTRIES 64\n");
			printf("#define DTB_SIZE 4096\n");
			printf("#define L2_ASSOCIATIVE 16\n");
			break;

		case CPU_THUNDERX2T99:
			printf("#define VULCAN                        \n");
			printf("#define HAVE_VFP                      \n");
			printf("#define HAVE_VFPV3                    \n");
			printf("#define HAVE_NEON                     \n");
			printf("#define HAVE_VFPV4                    \n");
			printf("#define L1_CODE_SIZE         32768    \n");
			printf("#define L1_CODE_LINESIZE     64       \n");
			printf("#define L1_CODE_ASSOCIATIVE  8        \n");
			printf("#define L1_DATA_SIZE         32768    \n");
			printf("#define L1_DATA_LINESIZE     64       \n");
			printf("#define L1_DATA_ASSOCIATIVE  8        \n");
			printf("#define L2_SIZE              262144   \n");
			printf("#define L2_LINESIZE          64       \n");
			printf("#define L2_ASSOCIATIVE       8        \n");
			printf("#define L3_SIZE              33554432 \n");
			printf("#define L3_LINESIZE          64       \n");
			printf("#define L3_ASSOCIATIVE       32       \n");
			printf("#define DTB_DEFAULT_ENTRIES  64       \n");
			printf("#define DTB_SIZE             4096     \n");
			break;
	}
}


void get_libname(void)
{
	int d = detect();
	printf("%s", cpuname_lower[d]);
}

void get_features(void)
{

#ifdef linux
	FILE *infile;
  	char buffer[2048], *p,*t;
  	p = (char *) NULL ;

  	infile = fopen("/proc/cpuinfo", "r");

	while (fgets(buffer, sizeof(buffer), infile))
	{

		if (!strncmp("Features", buffer, 8))
		{
			p = strchr(buffer, ':') + 2;
			break;
      		}
  	}

  	fclose(infile);


	if( p == NULL ) return;

	t = strtok(p," ");
	while( t = strtok(NULL," "))
	{
	}

#endif
	return;
}


