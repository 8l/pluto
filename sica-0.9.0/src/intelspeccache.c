/*
 * intelspeccache.h
 *
 *  Created on: 31.01.2013
 *      Author: dfeld
 */

#include "mathhelp.h"
#include "intelspeccache.h"

int get_intel_specific_cache_info(char *binary_a, char *binary_b, char *binary_c, char *binary_d, int BINARY_SIZE, int execution_time, int *l1_cache, int *l3_cache, int * check_leaf4)
	{
	int k;
	int execution_num=0;

	//________________________________________________
	//EAX READ-OUT
	//the EAX-register has to be treated in a unique way, as it carries the information for our induction variables

	//check Bit 31, it has to be 0 if information concerning TLB and cache are available
	if(binary_a[(BINARY_SIZE-1)-31]=='0')
	{
		//EAX's first byte enumerates, how often CPUID has to be called to get all pages of information concerning cache and TLB
		execution_num=bin2dec(binary_a, BINARY_SIZE, 0, 8);

		//read out remaining EAX byte-wise (3 bytes)
		for(k=1; k<4; k++)
		{
			int new_byte=bin2dec(binary_a, BINARY_SIZE, 8*k, 8);
			byte2info(new_byte, check_leaf4, l1_cache, l3_cache);

		}

	}

	//________________________________________________
	//EBX READ-OUT
	//check Bit 31, it has to be 0 if information concerning TLB and cache are available
	if(binary_b[(BINARY_SIZE-1)-31]=='0')
	{
		//read-out EBX byte-wise (4 bytes)
		for(k=0; k<4; k++)
		{
			int new_byte=bin2dec(binary_b, BINARY_SIZE, 8*k, 8);
			byte2info(new_byte, check_leaf4, l1_cache, l3_cache);

		}

	}

	//________________________________________________
	//ECX READ-OUT
	//check Bit 31, it has to be 0 if information concerning TLB and cache are available
	if(binary_c[(BINARY_SIZE-1)-31]=='0')
	{
		//read-out ECX byte-wise (4 bytes)
		for(k=0; k<4; k++)
		{
			int new_byte=bin2dec(binary_c, BINARY_SIZE, 8*k, 8);
			byte2info(new_byte, check_leaf4, l1_cache, l3_cache);

		}

	}

	//________________________________________________
	//EDX READ-OUT
	//check Bit 31, it has to be 0 if information concerning TLB and cache are available
	if(binary_d[(BINARY_SIZE-1)-31]=='0')
	{
		//read-out EDX byte-wise (4 bytes)
		for(k=0; k<4; k++)
		{
			int new_byte=bin2dec(binary_d, BINARY_SIZE, 8*k, 8);
			byte2info(new_byte, check_leaf4, l1_cache, l3_cache);
		}

	}


	//if not all necessary CPUID calls were performed, return true
	if(execution_time<execution_num)
	{
		return 1;
	}
	else
	{
		return 0;
	}

}


void byte2info(int new_byte, int * check_leaf4, int* l1_cache, int* l3_cache)
{
if(new_byte>0)
{
// check whether Leaf4 is necessary for further cache information
if(new_byte==0xFF)
  {
  *check_leaf4=1;
  }
}

set_l1_cache_parameter(new_byte, l1_cache);
set_l3_cache_parameter(new_byte, l3_cache);
}

//set the L1 cache information (from INTEL data sheet)
void set_l1_cache_parameter(int new_byte, int *l1_cache)
{
	if(new_byte==0x0A)
	{
		l1_cache[0]=8;
		l1_cache[1]=2;
		l1_cache[2]=32;
	}
	if(new_byte==0x0C)
	{
		l1_cache[0]=16;
		l1_cache[1]=4;
		l1_cache[2]=32;
	}
	if(new_byte==0x0D)
	{
		l1_cache[0]=16;
		l1_cache[1]=4;
		l1_cache[2]=64;
	}
	if(new_byte==0x2C)
	{
		l1_cache[0]=32;
		l1_cache[1]=8;
		l1_cache[2]=64;
	}
	if(new_byte==0x60)
	{
		l1_cache[0]=16;
		l1_cache[1]=8;
		l1_cache[2]=64;
	}
	if(new_byte==0x66)
	{
		l1_cache[0]=8;
		l1_cache[1]=4;
		l1_cache[2]=64;
	}
	if(new_byte==0x67)
	{
		l1_cache[0]=16;
		l1_cache[1]=4;
		l1_cache[2]=64;
	}
	if(new_byte==0x68)
	{
		l1_cache[0]=32;
		l1_cache[1]=4;
		l1_cache[2]=64;
	}
}

//L2 cache information can be determined numerically, even for INTEL

//set the L3 cache information (from INTEL data sheet)
void set_l3_cache_parameter(int new_byte, int *l3_cache)
{
	if(new_byte==0x22)
	{
		l3_cache[0]=512;
		l3_cache[1]=4;
		l3_cache[2]=64;
	}
	if(new_byte==0x23)
	{
		l3_cache[0]=1024;
		l3_cache[1]=8;
		l3_cache[2]=64;
	}
	if(new_byte==0x25)
	{
		l3_cache[0]=2048;
		l3_cache[1]=8;
		l3_cache[2]=64;
	}
	if(new_byte==0x29)
	{
		l3_cache[0]=4096;
		l3_cache[1]=8;
		l3_cache[2]=64;
	}
	if(new_byte==0x46)
	{
		l3_cache[0]=4096;
		l3_cache[1]=4;
		l3_cache[2]=64;
	}
	if(new_byte==0x47)
	{
		l3_cache[0]=8192;
		l3_cache[1]=8;
		l3_cache[2]=64;
	}
	if(new_byte==0x49)
	{
		l3_cache[0]=4096;
		l3_cache[1]=16;
		l3_cache[2]=64;
	}
	if(new_byte==0x4A)
	{
		l3_cache[0]=6144;
		l3_cache[1]=12;
		l3_cache[2]=64;
	}
	if(new_byte==0x4B)
	{
		l3_cache[0]=8192;
		l3_cache[1]=16;
		l3_cache[2]=64;
	}
	if(new_byte==0x4C)
	{
		l3_cache[0]=12288;
		l3_cache[1]=12;
		l3_cache[2]=64;
	}
	if(new_byte==0x4D)
	{
		l3_cache[0]=16384;
		l3_cache[1]=16;
		l3_cache[2]=64;
	}
	if(new_byte==0xD0)
	{
		l3_cache[0]=512;
		l3_cache[1]=4;
		l3_cache[2]=64;
	}
	if(new_byte==0xD1)
	{
		l3_cache[0]=1024;
		l3_cache[1]=4;
		l3_cache[2]=64;
	}
	if(new_byte==0xD2)
	{
		l3_cache[0]=2048;
		l3_cache[1]=4;
		l3_cache[2]=64;
	}
	if(new_byte==0xD6)
	{
		l3_cache[0]=1024;
		l3_cache[1]=8;
		l3_cache[2]=64;
	}
	if(new_byte==0xD7)
	{
		l3_cache[0]=2048;
		l3_cache[1]=8;
		l3_cache[2]=64;
	}
	if(new_byte==0xD8)
	{
		l3_cache[0]=4096;
		l3_cache[1]=8;
		l3_cache[2]=64;
	}
	if(new_byte==0xDC)
	{
		l3_cache[0]=1536;
		l3_cache[1]=12;
		l3_cache[2]=64;
	}
	if(new_byte==0xDD)
	{
		l3_cache[0]=3072;
		l3_cache[1]=12;
		l3_cache[2]=64;
	}
	if(new_byte==0xDE)
	{
		l3_cache[0]=6188;
		l3_cache[1]=12;
		l3_cache[2]=64;
	}
	if(new_byte==0xE2)
	{
		l3_cache[0]=2048;
		l3_cache[1]=16;
		l3_cache[2]=64;
	}
	if(new_byte==0xE3)
	{
		l3_cache[0]=4096;
		l3_cache[1]=16;
		l3_cache[2]=64;
	}
	if(new_byte==0xE4)
	{
		l3_cache[0]=8192;
		l3_cache[1]=16;
		l3_cache[2]=64;
	}
	if(new_byte==0xEA)
	{
		l3_cache[0]=12288;
		l3_cache[1]=24;
		l3_cache[2]=64;
	}
	if(new_byte==0xEB)
	{
		l3_cache[0]=18432;
		l3_cache[1]=24;
		l3_cache[2]=64;
	}
	if(new_byte==0xEC)
	{
		l3_cache[0]=24576;
		l3_cache[1]=24;
		l3_cache[2]=64;
	}
}


//string get_cache_tbl_string(int new_byte)
//{
//	if(new_byte==0x00)
//			return "NULL";
//	if(new_byte==0x01)
//			return "Instruction TLB: 4-KBPages, 4-way set associative, 32 entries";
//	if(new_byte==0x02)
//			return "Instruction TLB: 4-MB Pages, fully associative, 2 entries";
//	if(new_byte==0x03)
//			return "Data TLB: 4-KB Pages, 4-way set associative, 64 entries";
//	if(new_byte==0x04)
//			return "Data TLB: 4-MB Pages, 4-way set associative, 8 entries";
//	if(new_byte==0x05)
//			return "Data TLB: 4-MB Pages, 4-way set associative, 32 entries";
//	if(new_byte==0x06)
//			return "1 -level instruction cache: 8-KB, 4-way set associative, 32-byte line size";
//	if(new_byte==0x08)
//			return "1 -level instruction cache: 16-KB, 4-way set associative, 32-byte line size";
//	if(new_byte==0x09)
//			return "1st-level Instruction Cache: 32-KB, 4-way set associative, 64-byte line size";
//	if(new_byte==0x0A)
//			return "1 -level data cache: 8-KB, 2-way set associative, 32-byte line size";
//	if(new_byte==0x0C)
//			return "1 -level data cache: 16-KB, 4-way set associative, 32-byte line size";
//	if(new_byte==0x0D)
//			return "1st-level Data Cache: 16-KB, 4-way set associative, 64-byte line size, ECC";
//	if(new_byte==0x21)
//			return "256-KB L2 (MLC), 8-way set associative, 64-byte line size";
//	if(new_byte==0x22)
//			return "3 -level cache: 512 KB, 4-way set associative, sectored cache, 64-byte line size";
//	if(new_byte==0x23)
//			return "3 -level cache: 1-MB, 8-way set associative, sectored cache, 64-byte line size";
//	if(new_byte==0x25)
//			return "3 -level cache: 2-MB, 8-way set associative, sectored cache, 64-byte line size";
//	if(new_byte==0x29)
//			return "3 -level cache: 4-MB, 8-way set associative, sectored cache, 64-byte line size";
//	if(new_byte==0x2C)
//			return "1 -level data cache: 32-KB, 8-way set associative, 64-byte line size";
//	if(new_byte==0x30)
//			return "1 -level instruction cache: 32-KB, 8-way set associative, 64-byte line size";
//	if(new_byte==0x39)
//			return "2 -level cache: 128-KB, 4-way set associative, sectored cache, 64-byte line size";
//	if(new_byte==0x3A)
//			return "2nd-level cache: 192-KB, 6-way set associative, sectored cache, 64-byte line size";
//	if(new_byte==0x3B)
//			return "2 -level cache: 128-KB, 2-way set associative, sectored cache, 64-byte line size";
//	if(new_byte==0x3C)
//			return "2 -level cache: 256-KB, 4-way set associative, sectored cache, 64-byte line size";
//	if(new_byte==0x3D)
//			return "2nd-level cache: 384-KB, 6-way set associative, sectored cache, 64-byte line size";
//	if(new_byte==0x3E)
//			return "2nd-level cache: 512-KB, 4-way set associative, sectored cache, 64-byte line size";
//	if(new_byte==0x40)
//			return "No 2 -level cache or, if processor contains a valid 2 -level cache, no 3 -level cache";
//	if(new_byte==0x41)
//			return "2 -level cache: 128-KB, 4-way set associative, 32-byte line size";
//	if(new_byte==0x42)
//			return "2 -level cache: 256-KB, 4-way set associative, 32-byte line size";
//	if(new_byte==0x43)
//			return "2 -level cache: 512-KB, 4-way set associative, 32-byte line size";
//	if(new_byte==0x44)
//			return "2 -level cache: 1-MB, 4-way set associative, 32-byte line size";
//	if(new_byte==0x45)
//			return "2 -level cache: 2-MB, 4-way set associative, 32-byte line size";
//	if(new_byte==0x46)
//			return "3rd-level cache: 4-MB, 4-way set associative, 64-byte line size";
//	if(new_byte==0x47)
//			return "3rd-level cache: 8-MB, 8-way set associative, 64-byte line size";
//	if(new_byte==0x48)
//			return "2nd-level cache: 3-MB, 12-way set associative, 64-byte line size, unified on-die";
//	if(new_byte==0x49)
//			return "3rd-level cache: 4-MB, 16-way set associative, 64-byte line size (Intel Xeon processor MP, Family 0Fh, Model 06h)  / 2nd-level cache: 4-MB, 16-way set associative, 64-byte line size";
//	if(new_byte==0x4A)
//			return "3rd-level cache: 6-MB, 12-way set associative, 64-byte line size";
//	if(new_byte==0x4B)
//			return "3rd-level cache: 8-MB, 16-way set associative, 64-byte line size";
//	if(new_byte==0x4C)
//			return "3rd-level cache: 12-MB, 12-way set associative, 64-byte line size";
//	if(new_byte==0x4D)
//			return "3rd-level cache: 16-MB, 16-way set associative, 64-byte line size";
//	if(new_byte==0x4E)
//			return "2nd-level cache: 6-MB, 24-way set associative, 64-byte line size";
//	if(new_byte==0x50)
//			return "Instruction TLB: 4-KB, 2-MB or 4-MB pages, fully associative, 64 entries";
//	if(new_byte==0x51)
//			return "Instruction TLB: 4-KB, 2-MB or 4-MB pages, fully associative, 128 entries";
//	if(new_byte==0x52)
//			return "Instruction TLB: 4-KB, 2-MB or 4-MB pages, fully associative, 256 entries";
//	if(new_byte==0x55)
//			return "Instruction TLB: 2-MB or 4-MB pages, fully associative, 7 entries";
//	if(new_byte==0x56)
//			return "L1 Data TLB: 4-MB pages, 4-way set associative, 16 entries";
//	if(new_byte==0x57)
//			return "L1 Data TLB: 4-KB pages, 4-way set associative, 16 entries";
//	if(new_byte==0x5A)
//			return "Data TLB0: 2-MB or 4-MB pages, 4-way associative, 32 entries";
//	if(new_byte==0x5B)
//			return "Data TLB: 4-KB or 4-MB pages, fully associative, 64 entries";
//	if(new_byte==0x5C)
//			return "Data TLB: 4-KB or 4-MB pages, fully associative, 128 entries";
//	if(new_byte==0x5D)
//			return "Data TLB: 4-KB or 4-MB pages, fully associative, 256 entries";
//	if(new_byte==0x60)
//			return "1 -level data cache: 16-KB, 8-way set associative, sectored cache, 64-byte line size";
//	if(new_byte==0x66)
//			return "1 -level data cache: 8-KB, 4-way set associative, sectored cache, 64-byte line size";
//	if(new_byte==0x67)
//			return "1 -level data cache: 16-KB, 4-way set associative, sectored cache, 64-byte line size";
//	if(new_byte==0x68)
//			return "1 -level data cache: 32-KB, 4 way set associative, sectored cache, 64-byte line size";
//	if(new_byte==0x70)
//			return "Trace cache: 12K-uops, 8-way set associative";
//	if(new_byte==0x71)
//			return "Trace cache: 16K-uops, 8-way set associative";
//	if(new_byte==0x72)
//			return "Trace cache: 32K-uops, 8-way set associative";
//	if(new_byte==0x73)
//			return "Trace cache: 64K-uops, 8-way set associative";
//	if(new_byte==0x78)
//			return "2 -level cache: 1-MB, 4-way set associative, 64-byte line size";
//	if(new_byte==0x79)
//			return "2 -level cache: 128-KB, 8-way set associative, sectored cache, 64-byte line size";
//	if(new_byte==0x7A)
//			return "2 -level cache: 256-KB, 8-way set associative, sectored cache, 64-byte line size";
//	if(new_byte==0x7B)
//			return "2 -level cache: 512-KB, 8-way set associative, sectored cache, 64-byte line size";
//	if(new_byte==0x7C)
//			return "2 -level cache: 1-MB, 8-way set associative, sectored cache, 64-byte line size";
//	if(new_byte==0x7D)
//			return "2 -level cache: 2-MB, 8-way set associative, 64-byte line size";
//	if(new_byte==0x7F)
//			return "2 -level cache: 512-KB, 2-way set associative, 64-byte line size";
//	if(new_byte==0x82)
//			return "2 -level cache: 256-KB, 8-way set associative, 32-byte line size";
//	if(new_byte==0x83)
//			return "2 -level cache: 512-KB, 8-way set associative, 32-byte line size";
//	if(new_byte==0x84)
//			return "2 -level cache: 1-MB, 8-way set associative, 32-byte line size";
//	if(new_byte==0x85)
//			return "2 -level cache: 2-MB, 8-way set associative, 32-byte line size";
//	if(new_byte==0x86)
//			return "2 -level cache: 512-KB, 4-way set associative, 64-byte line size";
//	if(new_byte==0x87)
//			return "2 -level cache: 1-MB, 8-way set associative, 64-byte line size";
//	if(new_byte==0xB0)
//			return "Instruction TLB: 4-KB Pages, 4-way set associative, 128 entries";
//	if(new_byte==0xB1)
//			return "Instruction TLB: 2-MB pages, 4-way, 8 entries or 4M pages, 4-way, 4 entries";
//	if(new_byte==0xB2)
//			return "Instruction TLB: 4-KB pages, 4-way set associative, 64 entries";
//	if(new_byte==0xB3)
//			return "Data TLB: 4-KB Pages, 4-way set associative, 128 entries";
//	if(new_byte==0xB4)
//			return "Data TLB: 4-KB Pages, 4-way set associative, 256 entries";
//	if(new_byte==0xCA)
//			return "Shared 2nd-level TLB: 4 KB pages, 4-way set associative, 512 entries";
//	if(new_byte==0xD0)
//			return "512KB L3 Cache, 4-way set associative, 64-byte line size";
//	if(new_byte==0xD1)
//			return "1-MB L3 Cache, 4-way set associative, 64-byte line size";
//	if(new_byte==0xD2)
//			return "2-MB L3 Cache, 4-way set associative, 64-byte line size";
//	if(new_byte==0xD6)
//			return "1-MB L3 Cache, 8-way set associative, 64-byte line size";
//	if(new_byte==0xD7)
//			return "2-MB L3 Cache, 8-way set associative, 64-byte line size";
//	if(new_byte==0xD8)
//			return "4-MB L3 Cache, 8-way set associative, 64-byte line size";
//	if(new_byte==0xDC)
//			return "1.5-MB L3 Cache, 12-way set associative, 64-byte line size";
//	if(new_byte==0xDD)
//			return "3-MB L3 Cache, 12-way set associative, 64-byte line size";
//	if(new_byte==0xDE)
//			return "6-MB L3 Cache, 12-way set associative, 64-byte line size";
//	if(new_byte==0xE2)
//			return "2-MB L3 Cache, 16-way set associative, 64-byte line size";
//	if(new_byte==0xE3)
//			return "4-MB L3 Cache, 16-way set associative, 64-byte line size";
//	if(new_byte==0xE4)
//			return "8-MB L3 Cache, 16-way set associative, 64-byte line size";
//	if(new_byte==0xEA)
//			return "12-MB L3 Cache, 24-way set associative, 64-byte line size";
//	if(new_byte==0xEB)
//			return "18-MB L3 Cache, 24-way set associative, 64-byte line size";
//	if(new_byte==0xEC)
//			return "24-MB L3 Cache, 24-way set associative, 64-byte line size";
//	if(new_byte==0xF0)
//			return "64-byte Prefetching";
//	if(new_byte==0xF1)
//			return "128-byte Prefetching";
//	if(new_byte==0xFF)
//		return "[INFO] CPUID Leaf 2 does not report cache descriptor information; use CPUID Leaf 4 to query cache parameters";
//	else
//		return "No information found";
//
//}

