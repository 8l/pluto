/*
 * amdspeccache.c
 *
 *  Created on: 31.01.2013
 *      Author: dfeld
 */

#include "amdspeccache.h"


int get_amd_l1cache_size(char *binary_c, const int BINARY_SIZE)
{
	int pot,i,j;

	int l1_cache_size=0;
	for(i=0; i<8; i++)
		{
			pot = 1;

				for(j=0; j<i; j++)
				{
				pot *= 2;
				}

			if(binary_c[(BINARY_SIZE-1)-(i+24)]=='1')
			{
			l1_cache_size += pot;
			}
		}

	return l1_cache_size;

}

int get_amd_l1cache_line_size(char *binary_c, const int BINARY_SIZE)
{
	int pot,i,j;

	int l1_cache_line_size=0;
	for(i=0; i<8; i++)
		{
			pot = 1;

				for(j=0; j<i; j++)
				{
				pot *= 2;
				}

			if(binary_c[(BINARY_SIZE-1)-i]=='1')
			{
			l1_cache_line_size += pot;
			}
		}

	return l1_cache_line_size;
}

int get_amd_l1cache_assiociativity(char *binary_c, const int BINARY_SIZE)
{
	int pot,i,j;

	int l1_cache_assiociativity=0;
	for(i=0; i<8; i++)
		{
			pot = 1;

				for(j=0; j<i; j++)
				{
				pot *= 2;
				}

			if(binary_c[(BINARY_SIZE-1)-(i+16)]=='1')
			{
			l1_cache_assiociativity += pot;
			}
		}

	if((l1_cache_assiociativity>0)&&(l1_cache_assiociativity<255))
		return l1_cache_assiociativity;
	if(l1_cache_assiociativity==255)
		return 9999;
	else
		return 666;


}

