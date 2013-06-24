/*
 * sica_math_func.c
 *
 *  Created on: 25.02.2013
 *      Author: dfeld
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>

#include "sica_math_func.h"

////////////////////////////////////////////////////////////////////////////////
//                              ECHELON DETERMINANT                           //
////////////////////////////////////////////////////////////////////////////////

void sica_mult_matrices(float* mult1, float* mult2, float* res, int rows1, int columns1, int rows2, int columns2)
{
	if(columns1!=rows2)
	{
		printf("[ERROR] Dimensions in Matrix Multiplication do not fit!\n");
	}
	else
	{

	}

}

int sica_minusone_pow(int exponent)
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

void sica_print_matrix(float** matrix, int rows, int columns)
{
	int i,j;
	for(i=0; i<rows; i++)
	{
		for(j=0; j<columns; j++)
		{
			printf("%f \t",matrix[i][j]);
		}
		printf("\n");
	}
}

void sica_cpy_matrix(float** matrix, float** temp_matrix, int rows, int columns)
{
	int i, j;
	for(i=0; i<rows; i++)
		{
			for(j=0; j<columns; j++)
			{
				temp_matrix[i][j]=matrix[i][j];
			}
		}
}

void sica_cpy_matrixINT2FLOAT(int** matrix, float** temp_matrix, int rows, int columns)
{
	int i, j;
	for(i=0; i<rows; i++)
		{
			for(j=0; j<columns; j++)
			{
				temp_matrix[i][j]=(float)matrix[i][j];
			}
		}
}

void sica_cpy_matrixFLOAT2INT(float** matrix, int** temp_matrix, int rows, int columns)
{
	int i, j;
	for(i=0; i<rows; i++)
		{
			for(j=0; j<columns; j++)
			{
				temp_matrix[i][j]=(int)matrix[i][j];
			}
		}
}

int sica_echelon_form(float** det_matrix, int rows, int columns)
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
			if(det_matrix[j][j]==0)
			{
				int count=j;
				//find swappable row
				while((det_matrix[count][j]==0)&&(count<rows))
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
						help=det_matrix[j][t];
						det_matrix[j][t]=det_matrix[count][t];
						det_matrix[count][t]=help;
					}
					swap_counter++;
				}
			}

			if(finished_row==0)
			{
				factor=(det_matrix[i][j]/det_matrix[j][j]);
				//printf("Calculating with pivot in line %i: L%i=L%i-%f*L%i\n",j,i,i,factor,j);
				for(k=j;k<columns;k++)
				{
					det_matrix[i][k]=det_matrix[i][k]-det_matrix[j][k]*factor;
				}
			}
		}
	}
	return swap_counter;
}


float sica_echelon_determinant(float** matrix, int N)
{
	float det=1.0;
	int i;

	float** det_matrix;//for calculating the determinant
	 det_matrix=(float**)malloc(N*sizeof(float*));
	for(i=0; i<N; i++)    {
		 det_matrix[i]=(float*)malloc(N*sizeof(float));
	}
	//generates a temporary working copy for the transformation
	sica_cpy_matrix(matrix,  det_matrix, N, N);

	//returns the number of row swaps
	int swap_counter=sica_echelon_form(det_matrix, N, N);

	//Multiply diagonal of the echelon form
	for(i=0;i<N;i++)
	{
		det=det*det_matrix[i][i];
	}

	/* Swap-counter counts the number of row-swaps done, so that (-1)^sc must be
	 * multiplied
	 */
	det=det*sica_minusone_pow(swap_counter);

	return det;
}


////////////////////////////////////////////////////////////////////////////////
//                              INVERSE                                       //
////////////////////////////////////////////////////////////////////////////////
void sica_setup_identitymatrix(float** ident_matrix,int N)
{
	int i,j;
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{
			ident_matrix[i][j]=0.0;
		}
	}
	for(i=0;i<N;i++)
	{
		ident_matrix[i][i]=1.0;
	}
}

void sica_echelon_inverse(float** inverse_matrix, float** matrix,int N)
{
	float factor;
	int i,j,k;

	sica_setup_identitymatrix(inverse_matrix,N);

	//invert it!

	int finished_row=0;
	//generate Echelon Form without manipulating the determinant through
	//multipliers (besides -1 for row-swapping) as first inverting step
	for(j=0; j<N; j++)
	{
		for(i=j+1; i<N; i++)
		{
			//zero pivot element
			if(matrix[j][j]==0)
			{
				int count=j;
				//find addable row
				while((matrix[count][j]==0)&&(count<N))
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

						matrix[j][t]+=matrix[count][t];
						inverse_matrix[j][t]+=inverse_matrix[count][t];

					}
				}
			}

			if(finished_row==0)
			{
				factor=(matrix[i][j]/matrix[j][j]);
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
					matrix[i][k]=matrix[i][k]-matrix[j][k]*factor;
				}
				for(k=0;k<N;k++)
				{
				//manipulating the inverse-matrix elements
				inverse_matrix[i][k]=inverse_matrix[i][k]-inverse_matrix[j][k]*factor;
				}
			}
		}
	}

	float quotient;
	//Normalize the diagonal
	for(i=0;i<N;i++)
	{
		quotient=matrix[i][i];
		for(j=0;j<N;j++)
		{
			matrix[i][j]=matrix[i][j]/quotient;
			inverse_matrix[i][j]=inverse_matrix[i][j]/quotient;
		}
	}

	//reverse echelon transformation for upper triangle
	for(j=N-1; j>0; j--)
	{
		for(i=j-1; i>=0; i--)
		{
			factor=matrix[i][j];
			for(k=0;k<N;k++)
			{
				matrix[i][k]=matrix[i][k]-matrix[j][k]*factor;
				//manipulating the inverse elements
				inverse_matrix[i][k]=inverse_matrix[i][k]-inverse_matrix[j][k]*factor;
			}
		}
	}
}


void sica_inverse(int** matrix, int** matrix_inverted, int N)
{
	int i;

	//create temporaray float arrays for the calculation
	float** temp_matrix;//for calculating the determinant
	temp_matrix=(float**)malloc(N*sizeof(float*));
	for(i=0; i<N; i++)    {
		temp_matrix[i]=(float*)malloc(N*sizeof(float));
	}

	float** temp_matrix_inverted;//for calculating the determinant
	temp_matrix_inverted=(float**)malloc(N*sizeof(float*));
	for(i=0; i<N; i++)    {
		temp_matrix_inverted[i]=(float*)malloc(N*sizeof(float));
	}

	//generates a temporary working copy for the transformation
	sica_cpy_matrixINT2FLOAT(matrix, temp_matrix, N, N);
	//sica_cpy_matrixINT2FLOAT(matrix_inverted, temp_matrix_inverted, N, N);

	float det = sica_echelon_determinant(temp_matrix, N);
	if(det==0)
	{
		printf("[ERROR]: Transformation-Matrix not invertible! Please try without the --sica option. We would appreciate if you could send us the code that fails.\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		sica_echelon_inverse(temp_matrix_inverted,temp_matrix, N);

		sica_cpy_matrixFLOAT2INT(temp_matrix_inverted, matrix_inverted, N, N);
	}
}

