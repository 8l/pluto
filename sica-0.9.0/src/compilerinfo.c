/*
 * compilerinfo.c
 *
 *  Created on: 31.01.2013
 *      Author: dfeld
 */

#include <ctype.h>
#include <stdio.h>

void get_compiler_info()
{
		#define GCC_VERSION (__GNUC__ * 10000 \
                               + __GNUC_MINOR__ * 100 \
                               + __GNUC_PATCHLEVEL__)

//	//TODO
//	printf("[CACHE] Compiled with:\t\t\tGCC %i.%i.%i (number: %i)\n",__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__,GCC_VERSION);


	char exe_str[]="gcc -dumpversion"; //insted of --version
	FILE *ptr;

	if((ptr = popen(exe_str, "r")) == NULL)
	perror("Couldn't open pipe");

	char act_char;
	int found_digits=0;
	int version=0;
	int power=0;
	int version_seperated[3];

	while((act_char=getc(ptr))&&(act_char != EOF)&&(found_digits<3))
	{
		if(isdigit(act_char))
		{
			switch(found_digits)
			{
			case 0:
				power=10000;
				break;
			case 1:
				power=100;
				break;
			case 2:
				power=1;
				break;
			default:
				power=0;
				break;
			}

			version+=power*(act_char-48);
			version_seperated[found_digits]=act_char-48;
			/*-48 to convert from char to int*/
			found_digits++;
		}
	}

//	//TODO
//	printf("[CACHE] Runtime available Version:\tGCC %i.%i.%i (number: %i)\n", version_seperated[0], version_seperated[1], version_seperated[2],version);
//
//	if(version>40000)
//	{
//		printf("[CACHE] Vectorization is possible through GCC\n");
//		if(version<40500)
//		{
//			printf("[CACHE] BUT - We recommend to install a newer Version (4.5 or higher) of GCC on your System\n");
//		}
//	}
//	else
//	{
//		printf("[CACHE] Vectorization is NOT possible through GCC - need at least gcc4.0.1\n");
//	}

	fclose(ptr);
}

