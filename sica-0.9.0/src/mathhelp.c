/**
 * @file mathhelp.h
 * @author dfeld
 * @remark the 'dec2bin' implementation is based on a code that is widely spread among the internet 
 * Created on: 31.01.2013
 */

#include "mathhelp.h"

      // accepts a decimal integer and returns a binary coded string of size BINARY_SIZE
      //
      void dec2bin(long long decimal, char *binary, const int BINARY_SIZE)
      {
      int k = 0, n = 0;
      int neg_flag = 0;
      int remain;
      char temp[BINARY_SIZE];  //temp-Array passend zum binary-Array

      // take care of negative input
      if (decimal < 0)
      {
    	  decimal = -decimal;	//sould be obsolete because of cerr
    	  neg_flag = 1;	
      }
      do
      {
      remain = decimal % 2;
      // whittle down the decimal number
      decimal = decimal / 2;
      // converts digit 0 or 1 to character '0' or '1'
      temp[k++] = remain + '0';
      }
	while (decimal > 0);

     	 if (neg_flag)
     	 temp[k++] = '-'; // add - sign
     	 else
     	 temp[k++] = '+'; // add + instead of space sign

     	 // reverse the spelling
     	 binary[n++] = temp[--k]; //add the sign digit

     	// 32-k zeros to get a standard 32-bit-string
     	 int i;
     	 for(i=0; i<((BINARY_SIZE-1)-k); i++)
     		 binary[n++]='0';


     	while (k >= 0)
     	binary[n++] = temp[--k];

      }

int bin2dec(char *binary_array, const int BINARY_SIZE, int start, int length)
{
	int pot,value=0,i,j;

	for(i=0; i<length; i++)
	{
		pot = 1;

			for(j=0; j<i; j++)
			{
			pot *= 2;
			}

		if(binary_array[(BINARY_SIZE-1)-(start+i)]=='1')
		{
		value += pot;
		}
	}
	return value;
}


