/*
 * hwanalysis.c
 *
 *  Created on: 31.01.2013
 *      Author: dfeld
 */

#include <stdio.h>

#include "hwanalysis.h"

#include "cpuinfo.h"
#include "intelspeccache.h"
#include "compilerinfo.h"

int estimate_sse_version(int* simd_vers)
{
int version_number=0;
int i;
int count=0;

for(i=0;i<7;i++)
  {
    if(simd_vers[i])
      {
      version_number=i;
      count++;
      }
  }

  //if no vector unit was detected
  if(count==0)
    {
    version_number=-1;
    }

return version_number;
}



int get_hardware_cache_infos(int option)
{
	// if the given option is not legal, this functions returns the negative
	// value -1, otherwise the wanted value
	/* List of Options:
         * CPU_ID CPU_NAME
         * SSE_ID SSE_NAME
         * L1CACHE_SIZE L1CACHE_ASSO L1CACHE_LINE
         * L2CACHE_SIZE L2CACHE_ASSO L2CACHE_LINE
         * L3CACHE_SIZE L3CACHE_ASSO L3CACHE_LINE
	 */

	int simd_vers[7]={0,0,0,0,0,0,0};
	 /*
	  * array to check available SIMD units
	  * [0] - MMX
	  * [1] - SSE1
	  * [2] - SSE2
	  * [3] - SSE3
	  * [4] - SSE4.1
	  * [5] - SSE4.2
	  * [6] - AVX
	  */

	 int cpu_fam=0;
	 /*
	  * integer to identify processor family
	  * 0 - unknown
	  * 1 - GenuineIntel
	  * 2 - AuthenticAMD
	  */

	int l1_cache[3]={0,0,0};
	/*
	 * [0] - cache-size in KByte
	 * [1] - cache-associativity in -way
	 * [2] - cache-line-size in byte
	 */

	int l2_cache[3]={0,0,0};
	/*
	 * [0] - cache-size in KByte
	 * [1] - cache-associativity in -way
	 * [2] - cache-line-size in byte
	 */

	int l3_cache[3]={0,0,0};
	/*
	 * [0] - cache-size in KByte
	 * [1] - cache-associativity in -way
	 * [2] - cache-line-size in byte
	 */

	get_cpu_info(simd_vers,l1_cache,l2_cache,l3_cache,&cpu_fam);

	get_compiler_info();

	if(option==CPU_ID)
		return cpu_fam;
	if(option==SSE_ID)
		return estimate_sse_version(&simd_vers[0]);
	if(option==L1CACHE_SIZE)
		return l1_cache[0];
	if(option==L2CACHE_SIZE)
		return l2_cache[0];
	if(option==L3CACHE_SIZE)
		return l3_cache[0];
	if(option==L1CACHE_ASSO)
		return l1_cache[1];
	if(option==L2CACHE_ASSO)
		return l2_cache[1];
	if(option==L3CACHE_ASSO)
		return l3_cache[1];
	if(option==L1CACHE_LINE)
		return l1_cache[2];
	if(option==L2CACHE_LINE)
		return l2_cache[2];
	if(option==L3CACHE_LINE)
		return l3_cache[2];
	else
		return -1;
}
