/*
 * cpuinfo.c
 *
 *  Created on: 31.01.2013
 *      Author: dfeld
 *
 * @remark the 'cpuid makro' that is used is widely spread among the internet 
 */

#include "cpuinfo.h"

#include "intelspeccache.h"
#include "amdspeccache.h"

#include "mathhelp.h"

#include <string.h>

#define cpuid(func_a,func_c,ax,bx,cx,dx)\
	__asm__ __volatile__ ("cpuid":\
	"=a" (ax), "=b" (bx), "=c" (cx), "=d" (dx) : "a" (func_a), "c" (func_c));

void get_register_info(unsigned int function_a, unsigned int function_c, char *binary_a, char *binary_b, char *binary_c, char *binary_d, const int BINARY_SIZE)
{
unsigned int a,b,c,d;

cpuid(function_a,function_c,a,b,c,d);
dec2bin(a,binary_a,BINARY_SIZE);
dec2bin(b,binary_b,BINARY_SIZE);
dec2bin(c,binary_c,BINARY_SIZE);
dec2bin(d,binary_d,BINARY_SIZE);

}


void get_cpu_name(char * cpu_name, char *binary_b, char *binary_c, char *binary_d, const int BINARY_SIZE)
{

	int i,char_counter=0;


	//EBX-register read-out
	int counter=0;
	char next_byte=0;
	int pot;
	while(counter<(BINARY_SIZE-1))
	{
	pot = 1;

		for(i=0; i<(counter%8); i++)
		{
		pot *= 2;
		}

	if(binary_b[(BINARY_SIZE-1)-counter]=='1')
	{
	next_byte += pot;
	}

	if(((counter+1)%8)==0)
	{
	cpu_name[char_counter++]=next_byte;
	next_byte=0;
	}

	counter++;


	}



	//EDX-register read-out
	counter=0;
	next_byte=0;
	while(counter<(BINARY_SIZE-1))
	{
	pot = 1;

		for(i=0; i<(counter%8); i++)
		{
		pot *= 2;
		}

	if(binary_d[(BINARY_SIZE-1)-counter]=='1')
	{
	next_byte += pot;
	}

	if(((counter+1)%8)==0)
	{
	cpu_name[char_counter++]=next_byte;
	next_byte=0;
	}

	counter++;


	}


	//ECX-register read-out
	counter=0;
	next_byte=0;
	while(counter<(BINARY_SIZE-1))
	{
	pot = 1;

		for(i=0; i<(counter%8); i++)
		{
		pot *= 2;
		}

	if(binary_c[(BINARY_SIZE-1)-counter]=='1')
	{
	next_byte += pot;
	}

	if(((counter+1)%8)==0)
	{
	cpu_name[char_counter++]=next_byte;
	next_byte=0;
	}

	counter++;


	}

	cpu_name[char_counter]='\0';

}

//determine the L1 und L2 cache sizes
int get_intelamd_l2cache_size(char *binary_c, const int BINARY_SIZE)
{
	return bin2dec(binary_c, BINARY_SIZE, 16, 16);;
}

int get_intelamd_l2cache_line_size(char *binary_c, const int BINARY_SIZE)
{
	return bin2dec(binary_c, BINARY_SIZE, 0, 8);
}


int get_intelamd_cacheline_size(char *binary_b, const int BINARY_SIZE)
{
	return bin2dec(binary_b, BINARY_SIZE, 8, 8);
}




int get_intelamd_l2cache_assiociativity(char *binary_c, const int BINARY_SIZE)
{
	int l2_cache_assiociativity=bin2dec(binary_c, BINARY_SIZE, 12, 4);

	//data taken from AMD- and INTEL-datasheets | CPUID-amd p.25 and intel_cpuid.pdf p.38 from docs
	if(l2_cache_assiociativity==0)
		return 0;
	if(l2_cache_assiociativity==1)
		return 1;
	if(l2_cache_assiociativity==2)
		return 2;
	if(l2_cache_assiociativity==4)
		return 4;
	if(l2_cache_assiociativity==6)
		return 8;
	if(l2_cache_assiociativity==8)
		return 16;
	if(l2_cache_assiociativity==10)
		return 32;
	if(l2_cache_assiociativity==11)
		return 48;
	if(l2_cache_assiociativity==12)
		return 64;
	if(l2_cache_assiociativity==13)
		return 96;
	if(l2_cache_assiociativity==14)
		return 128;
	if(l2_cache_assiociativity==15)
		return 9999; //all-associative
	else
		return 666;	//error in determination
}


void get_cpu_info(int * simd_vers, int * l1_cache,int * l2_cache,int *l3_cache,int * cpu_fam)
{
 const int BINARY_SIZE=33;		//32 bits plus sign=33
 char binary_a[BINARY_SIZE];
 char binary_b[BINARY_SIZE];
 char binary_c[BINARY_SIZE];
 char binary_d[BINARY_SIZE];


//processor information

get_register_info(0x0,0x0, binary_a, binary_b, binary_c, binary_d, BINARY_SIZE);

char cpu_name[256];
get_cpu_name(cpu_name, binary_b, binary_c, binary_d, BINARY_SIZE);

if(cpu_name[0]=='G')
{
	if(cpu_name[7]=='I')
	{
		*cpu_fam=1;
	}
}
else
{
	if(cpu_name[0]=='A')
	{
		*cpu_fam=2;
	}
	else
	{
		*cpu_fam=0;
	}
}



//register information

get_register_info(0x1,0x0, binary_a, binary_b, binary_c, binary_d, BINARY_SIZE);

//Check for MMX-Support
if(binary_d[(BINARY_SIZE-1)-23]=='1') //+1 because of sign digit
{
	simd_vers[0]=1;
}

//Check for SSE1-Support
if(binary_d[(BINARY_SIZE-1)-25]=='1') //+1 because of sign digit
{
	simd_vers[1]=1;
}

//Check for SSE2-Support
if(binary_d[(BINARY_SIZE-1)-26]=='1') //+1 because of sign digit
{
	simd_vers[2]=1;
}

//Check for SSE3-Support
if(binary_c[(BINARY_SIZE-1)-0]=='1') //+1 because of sign digit
{
	simd_vers[3]=1;
}

//Check for SSE4.1-Support
if(binary_c[(BINARY_SIZE-1)-19]=='1') //+1 because of sign digit
{
	simd_vers[4]=1;
}

//Check for SSE4.2-Support
if(binary_c[(BINARY_SIZE-1)-20]=='1') //+1 because of sign digit
{
	simd_vers[5]=1;
}

//Check for AVX-Support
if(binary_c[(BINARY_SIZE-1)-28]=='1') //+1 because of sign digit
{
	simd_vers[6]=1;
}

get_register_info(0x1,0x0, binary_a, binary_b, binary_c, binary_d, BINARY_SIZE);

//Check whether cpu is INTEL
if(*cpu_fam==1)
{
		get_register_info(0x80000006,0x0, binary_a, binary_b, binary_c, binary_d, BINARY_SIZE);

		l2_cache[0]=get_intelamd_l2cache_size(binary_c, BINARY_SIZE);
		l2_cache[1]=get_intelamd_l2cache_assiociativity(binary_c, BINARY_SIZE);
		l2_cache[2]=get_intelamd_l2cache_line_size(binary_c, BINARY_SIZE);

		int execution_time=0;
		int go_on=1;       //the cpuid instruction may have to be executed several times to get all information, go_on check whether to do so or not 
				       //the number of necessary iterations is available in EAX

		int check_leaf4=0;     //variable to catch register code FF abfaengt to describe whether Leaf4 has to be used to read all cache information (for newer cpus)

		while(go_on)
		{
			execution_time++;
			get_register_info(0x2,0x0, binary_a, binary_b, binary_c, binary_d, BINARY_SIZE);
			go_on = get_intel_specific_cache_info(binary_a, binary_b, binary_c, binary_d, BINARY_SIZE, execution_time, l1_cache, l3_cache, &check_leaf4);
		}

		//CPUID Leaf 4 Test
		int cache_type=1;
		if(check_leaf4)
		{
		int counter=0;
		int cache_level, associativity, partitions, line_size, sets, cache_size;

		while(cache_type>0)//Todo: Check whether all informations are captured until the first 4 elements of a are zero
		{
			cache_type=0;
			cache_level=0;
			associativity=0;
			partitions=0;
			line_size=0;
			sets=0;
			cache_size=0;

			get_register_info(0x4,counter, binary_a, binary_b, binary_c, binary_d, BINARY_SIZE);

			cache_type=0;
			cache_type=bin2dec(binary_a, BINARY_SIZE, 0, 5);
			cache_level=bin2dec(binary_a, BINARY_SIZE, 5, 3);
			associativity=bin2dec(binary_b, BINARY_SIZE, 22, 10);
			if(binary_a[9]=='1')	//Fully assiociative
			{
				associativity=9999;//TODO: Fix because of calculated exhaustive cache size when fully ass. what to do?
			}
			partitions=bin2dec(binary_b, BINARY_SIZE, 12, 10);
			line_size=bin2dec(binary_b, BINARY_SIZE, 0, 12);
			sets=bin2dec(binary_c, BINARY_SIZE, 0, 32);

			//Calculate the resulting cache_size: (INTEL CPUID p. 35)
			cache_size=(associativity+1)*(partitions+1)*(line_size+1)*(sets+1);

			//SAVE the results to the arrays
			if(cache_level==1 && cache_type==1) //Cache-Type==1 to get the DATA-Cache and NOT THE INSTRUCTION CACHE (2)
			{
				l1_cache[0]=cache_size/1024; //divided by 1024 to get KByte info
				l1_cache[1]=associativity+1;
				l1_cache[2]=line_size+1;
			}
			if(cache_level==2)
			{
				l2_cache[0]=cache_size/1024;
				l2_cache[1]=associativity+1;
				l2_cache[2]=line_size+1;
			}
			if(cache_level==3)
			{
				l3_cache[0]=cache_size/1024;
				l3_cache[1]=associativity+1;
				l3_cache[2]=line_size+1;
			}

			counter++;
		}
		}
}
else
{
        //Check whether cpu is AMD
	if(*cpu_fam==2)
	{
	get_register_info(0x80000006,0x0, binary_a, binary_b, binary_c, binary_d, BINARY_SIZE);

	l2_cache[0]=get_intelamd_l2cache_size(binary_c, BINARY_SIZE);
	l2_cache[1]=get_intelamd_l2cache_assiociativity(binary_c, BINARY_SIZE);
	l2_cache[2]=get_intelamd_l2cache_line_size(binary_c, BINARY_SIZE);

	get_register_info(0x80000005,0x0, binary_a, binary_b, binary_c, binary_d, BINARY_SIZE);

	l1_cache[0]=get_amd_l1cache_size(binary_c, BINARY_SIZE);
	l1_cache[1]=get_amd_l1cache_assiociativity(binary_c, BINARY_SIZE);
	l1_cache[2]=get_amd_l1cache_line_size(binary_c, BINARY_SIZE);
	}
}

}

