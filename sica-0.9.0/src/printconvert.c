/*
 * printconvert.c
 *
 *  Created on: 31.01.2013
 *      Author: dfeld
 */

#include <string.h>
#include <stdio.h>

#include "hwanalysis.h"

void print_all_cache_information()
{
    printf("[SICA] CACHE-INFO\n");
    printf("[SICA] CPU-NAME: \t\t%s\n", id2cpu(get_hardware_cache_infos(CPU_ID)));
    printf("[SICA] SSE-NAME: \t\t%s\t\t%i Bit register\n", id2sse(get_hardware_cache_infos(SSE_ID)), id2regsize(get_hardware_cache_infos(SSE_ID)));
    printf("[SICA] L1-CACHE: \t\t%i KByte\t%i-way associative\t%i Byte cachelinesize\n", get_hardware_cache_infos(L1CACHE_SIZE), get_hardware_cache_infos(L1CACHE_ASSO), get_hardware_cache_infos(L1CACHE_LINE));
    printf("[SICA] L2-CACHE: \t\t%i KByte\t%i-way associative\t%i Byte cachelinesize\n", get_hardware_cache_infos(L2CACHE_SIZE), get_hardware_cache_infos(L2CACHE_ASSO), get_hardware_cache_infos(L2CACHE_LINE));
    printf("[SICA] L3-CACHE: \t\t%i KByte\t%i-way associative\t%i Byte cachelinesize\n", get_hardware_cache_infos(L3CACHE_SIZE), get_hardware_cache_infos(L3CACHE_ASSO), get_hardware_cache_infos(L3CACHE_LINE));
}

void print_l1cache_hierarchie(int regsize, int l1cachesize)
{
    printf("[SICA] Optimization for the following elements:\n\n");
    printf("                %s CORE\n",id2cpu(get_hardware_cache_infos(CPU_ID)));
    printf("\n");
    printf("                --------- \n");
    printf("%s            | %i-Bit |\n",id2sse(get_hardware_cache_infos(SSE_ID)),regsize);
    printf("                --------- \n");
    printf("                    |\n");
    printf("               ------------ \n");
    printf("L1 Cache      |  %i KByte  |\n",l1cachesize);
    printf("               ------------ \n");
}

void print_addl2cache_hierarchie(int l2cachesize)
{
        printf("                    |\n");
        printf("             ----------------- \n");
        printf("L2 Cache    |    %i KByte    |\n",l2cachesize);
        printf("             ----------------- \n");
}

char* id2sse(int identifier)
{
char* ssestring;
ssestring = (char*)malloc(8*sizeof(char)); //8 characters for STRING

if(identifier==-1)
{
strcpy(&ssestring[0], "NO"); 
}
if(identifier==0)
{
strcpy(&ssestring[0], "MMX"); 
}
if(identifier==1)
{
strcpy(&ssestring[0], "SSE1"); 
}
if(identifier==2)
{
strcpy(&ssestring[0], "SSE2"); 
}
if(identifier==3)
{
strcpy(&ssestring[0], "SSE3"); 
}
if(identifier==4)
{
strcpy(&ssestring[0], "SSE4.1"); 
}
if(identifier==5)
{
strcpy(&ssestring[0], "SSE4.2"); 
}
if(identifier==6)
{
strcpy(&ssestring[0], "AVX"); 
}

return ssestring;

}


char* id2cpu(int identifier)
{
char* cpustring;
cpustring = (char*)malloc(8*sizeof(char)); //8 characters for STRING

if(identifier==0)
{
strcpy(&cpustring[0], "UNKNOWN"); 
}
if(identifier==1)
{
strcpy(&cpustring[0], "INTEL"); 
}
if(identifier==2)
{
strcpy(&cpustring[0], "AMD"); 
}

return cpustring;

}


int id2regsize(int identifier)
{
int regsize;

if(identifier==0)
{
regsize = 64;
}
if((identifier>0)&&(identifier<6))
{
regsize = 128;
}
if(identifier== 6)
{
regsize = 256; 
}

return regsize;

}

