/*
 * cache_math_func.c
 *
 *  Created on: 25.02.2013
 *      Author: dfeld
 */

#include <stdio.h>

#include "cache_math_func.h"

////////////////////////////////////////////////////////////////////////////////
//                              ECHELON DETERMINANT                           //
////////////////////////////////////////////////////////////////////////////////

void cache_mult_matrices(float* mult1, float* mult2, float* res, int rows1, int columns1, int rows2, int columns2)
{
	if(columns1!=rows2)
	{
		printf("[ERROR] Dimensions in Matrix Multiplication do not fit!\n");
	}
	else
	{

	}

}

int cache_minusone_pow(int exponent)
{
	if((exponent%2)==0)
	{
		return 1;
	}
	else
	{
		return -1;
	}
}

void cache_print_matrix(float* matrix, int rows, int columns)
{
	int i,j;
	for(i=0; i<rows; i++)
	{
		for(j=0; j<columns; j++)
		{
			printf("%f \t",matrix[columns*i+j]);
		}
		printf("\n");
	}
}

void cache_cpy_matrix(float* matrix, float* temp_matrix, int rows, int columns)
{
	int i, j;
	for(i=0; i<rows; i++)
		{
			for(j=0; j<columns; j++)
			{
				temp_matrix[columns*i+j]=matrix[i*columns+j];
			}
		}
}

int cache_echelon_form(float* temp_matrix, int rows, int columns)
{
	float factor;
	int i,j,k;
	int swap_counter=0;
	int finished_row=0;
	/* generate Echelon Form without manipulating the determinant through
	 * multipliers (besides -1 for row-swapping)
	 */
	for(j=0; j<columns; j++)
	{
		for(i=j+1; i<rows; i++)
		{
			//zero pivot element
			if(temp_matrix[columns*j+j]==0)
			{
				int count=j;
				//find swappable row
				while((temp_matrix[columns*count+j]==0)&&(count<rows))
				{
					count++;
				}
				if(count>(rows-1))
				{
					finished_row=1;
				}

				if(finished_row==0)
				{
					//swap two rows
					int t;
					float help;
					for(t=j;t<columns;t++)
					{
						help=temp_matrix[columns*j+t];
						temp_matrix[columns*j+t]=temp_matrix[columns*count+t];
						temp_matrix[columns*count+t]=help;
					}
					swap_counter++;
				}
			}

			if(finished_row==0)
			{
				factor=(temp_matrix[columns*i+j]/temp_matrix[columns*j+j]);
				for(k=j;k<columns;k++)
				{
					temp_matrix[columns*i+k]=temp_matrix[columns*i+k]-temp_matrix[columns*j+k]*factor;
				}
			}
		}
	}
	return swap_counter;
}


float cache_echelon_determinant(float* temp_matrix, int N)
{
	float det=1.0;
	int i;
	//returns the number of row swaps
	int swap_counter=cache_echelon_form(temp_matrix, N, N);

	//Multiply diagonal of the echelon form
	for(i=0;i<N;i++)
	{
		det=det*temp_matrix[N*i+i];
	}

	/* Swap-counter counts the number of row-swaps done, so that (-1)^sc must be
	 * multiplied
	 */
	det=det*cache_minusone_pow(swap_counter);

	return det;
}


float cache_determinant(float* matrix, int N)
{
	float temp_matrix[N*N];
	//generates a temporary working copy for the transformation
	cache_cpy_matrix(matrix, temp_matrix, N, N);
	float det = cache_echelon_determinant(temp_matrix, N);
	return det;
}


////////////////////////////////////////////////////////////////////////////////
//                              INVERSE                                       //
////////////////////////////////////////////////////////////////////////////////
void cache_setup_identitymatrix(float* ident_matrix,int N)
{
	int i,j;
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			ident_matrix[i*N+j]=0;
		}
	}
	for(i=0;i<N;i++)
	{
		ident_matrix[i*N+i]=1;
	}
}

void cache_echelon_inverse(float* inverse_matrix, float* matrix,int N)
{
	float factor;
	int i,j,k;

	cache_setup_identitymatrix(inverse_matrix,N);

	//invert it!

	int finished_row=0;
	//generate Echelon Form without manipulating the determinant through
	//multipliers (besides -1 for row-swapping) as first inverting step
	for(j=0; j<N; j++)
	{
		for(i=j+1; i<N; i++)
		{
			//zero pivot element
			if(matrix[N*j+j]==0)
			{
				int count=j;
				//find addable row
				while((matrix[N*count+j]==0)&&(count<N))
				{
					count++;
				}
				if(count>(N-1))
				{
					finished_row=1;
				}

				if(finished_row==0)
				{
					//ADD TWO ROWS
					int t;

					for(t=0;t<N;t++)
					{

						matrix[N*j+t]+=matrix[N*count+t];
						inverse_matrix[N*j+t]+=inverse_matrix[N*count+t];

					}
				}
			}

			if(finished_row==0)
			{
				factor=(matrix[N*i+j]/matrix[N*j+j]);
				/*
				 * ToDo: Check whether it is correct to just add from J on,
				 * because there should be only zeros left of it, is it possible
				 * that there are non-zeros because of row adding (I don't think
				 * so, but am tired :-) )
				 *
				 * If there are problems, the for-loop should start with 0
				 * instead of j, no problem, just not necessary
				 */
				for(k=j;k<N;k++)
				{
					//manipulating the source matrix
					matrix[N*i+k]=matrix[N*i+k]-matrix[N*j+k]*factor;
				}
				for(k=0;k<N;k++)
				{
				//manipulating the inverse-matrix elements
				inverse_matrix[N*i+k]=inverse_matrix[N*i+k]-inverse_matrix[N*j+k]*factor;
				}
			}
		}
	}

	float quotient;
	//Normalize the diagonal
	for(i=0;i<N;i++)
	{
		quotient=matrix[N*i+i];
		for(j=0;j<N;j++)
		{
			matrix[N*i+j]=matrix[N*i+j]/quotient;
			inverse_matrix[N*i+j]=inverse_matrix[N*i+j]/quotient;
		}
	}

	//reverse echelon transformation for upper triangle
	for(j=N-1; j>0; j--)
	{
		for(i=j-1; i>=0; i--)
		{
			factor=matrix[N*i+j];
			for(k=0;k<N;k++)
			{
				matrix[N*i+k]=matrix[N*i+k]-matrix[N*j+k]*factor;
				//manipulating the inverse elements
				inverse_matrix[N*i+k]=inverse_matrix[N*i+k]-inverse_matrix[N*j+k]*factor;
			}
		}
	}
}


void cache_inverse(float* matrix, float* inverse_matrix, int N)
{
	float temp_matrix[N*N];//for calculating the determinant
	//generates a temporary working copy for the transformation
	cache_cpy_matrix(matrix, temp_matrix, N, N);

	float det = cache_echelon_determinant(temp_matrix, N);
	if(det==0)
	{
		printf("[ERROR]: Matrix not invertible!!!\n");
	}
	else
	{
		cache_echelon_inverse(inverse_matrix,matrix, N);
	}
}

void cache_vec_times_matrix(float* solution_vec, float* vector, float* matrix, int rows, int columns)
{
	int i,j;
	float accu;
	for(j=0;j<columns;j++)
	{
		accu=0;
		for(i=0;i<rows;i++)
		{
			accu+=vector[i]*matrix[i*columns+j];
		}
		solution_vec[j]=accu;
	}
}
