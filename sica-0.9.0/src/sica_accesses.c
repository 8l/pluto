/*
 * sica_accesses.c
 *
 *  Created on: 04.03.2013
 *      Author: dfeld
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>

#include "hwanalysis.h"

#include "pluto.h"

#include "sica_accesses.h"

#include "cache_math_func.h"

int** sica_access_matrix_malloc(int rows, int cols)
{
	int y;
	int** the_matrix=(int**)malloc(rows*sizeof(int*));

	for(y=0; y<rows; y++)    {
		the_matrix[y]=(int*)malloc(cols*sizeof(int));
	}

	return the_matrix;
}

SICAAccess* sica_accesses_malloc()    {
	SICAAccess* sicaacc;

	sicaacc = (SICAAccess*)malloc(sizeof(SICAAccess));
	sicaacc->nrows          = 0;
	sicaacc->ncols          = 0;

	sicaacc->next           = NULL;

	return sicaacc;
}

void sica_copy_access_matrix(int** matrix1, int** matrix2, int rows, int cols)    {
	int x,y;

   	for(y=0;y<rows;y++)    {
   		for(x=0;x<cols;x++)    {
   			matrix1[y][x]=matrix2[y][x];
   		}
	}
}

int sica_compare_access_matrices(int** matrix1, int** matrix2, int rows, int cols)    {
	int x,y;

	int ret_val=ACCESS_IS_IDENTICAL;

   	for(y=0;y<rows;y++)    {
   		for(x=0;x<cols;x++)    {
   			if(matrix1[y][x]!=matrix2[y][x])    {
   				ret_val=ACCESS_IS_NOT_IDENTICAL;
   			}
		}
   	}

   	//CHECK WHETHER THE DIFFERENCE IS ONLY IN THE VECTORIZED STRIDE
   	if(ret_val==ACCESS_IS_NOT_IDENTICAL)    {
   		ret_val=ACCESS_IS_IN_SAME_STRIDE;
   		//CHECK THE LAST LINE IF IT IS IDENTICAL BESIDES THE LAST (CONSTANT) ELEMENT
	   		for(x=0;x<cols-1;x++)    {
	   			if(matrix1[rows-1][x]!=matrix2[rows-1][x])    {
	   				ret_val=ACCESS_IS_NOT_IDENTICAL;
	   			}
			}

	   	//CHECK THE OTHER ROWS
   	   	for(y=0;y<rows-1;y++)    {
   	   		for(x=0;x<cols;x++)    {
   	   			if(matrix1[y][x]!=matrix2[y][x])    {
   	   				ret_val=ACCESS_IS_NOT_IDENTICAL;
   	   			}
   			}
   	   	}

   	}

   	return ret_val;
}


int sica_get_entry_sum(int** matrix, int rows, int cols)    {
	int x,y;

	int entry_sum=0;

   	for(y=0;y<rows;y++)    {
   		for(x=0;x<cols;x++)    {
   			entry_sum+=matrix[y][x];
		}
   	}

   	return entry_sum;
}


void sica_vec_times_matrix(int* solution_vec, int* vector, int** matrix, int rows, int cols)
{
	int i,j;
	int accu;
	for(j=0;j<cols;j++)
	{
		accu=0;
		for(i=0;i<rows;i++)
		{
			accu+=vector[i]*matrix[i][j];
		}
		solution_vec[j]=accu;
	}
}

int sica_get_array_id(Band* act_band, char* array_name)    {
	int a;
	for(a=0; a < act_band->sicadata->nb_arrays; a++)    {
		if(!strcmp(act_band->sicadata->id2arrayname[a], array_name))    {
			return a;
		}
    }

	return -1;
}

void sica_set_vectorized_bands(Band** bands, int nbands)    {
	int i;
    for (i=0; i<nbands; i++) {
    	if(bands[i]->sicadata->vecloop>-1)
    	{
        bands[i]->sicadata->isvec=1;
        //ToDo: Here seems to arise an error that the loop index is always one too big (not the difference of starting counting at 0 instead of 1)
        IF_DEBUG(printf("[SICA] BAND %i, Original loop that is vectorized: t%i, row in tile_sizes: %i\n",i,bands[i]->sicadata->vecloop+1,bands[i]->sicadata->vecrow+1););
    	}else{
        bands[i]->sicadata->isvec=0;
    	}
    }
}

void sica_get_trans_matrix(Band** bands, int nbands)    {
	int i,x,y,s;
	for (i=0; i<nbands; i++) {
      Band* act_band=bands[i];
      for(s=0; s<act_band->loop->nstmts;s++)    {
    	//calculatin the column offset -> TODO: SCALAR DIMENSIONS ARE MISSING, TAKE THE TRANSFORMATION FROM PROG???
        int firstD = act_band->loop->depth;
        int lastD = act_band->loop->depth+act_band->width-1;
        int widthD=lastD-firstD;

        int coloffset=0;
        act_band->sicadata->transwidth=widthD+1;

        if(options->tile)    {
        	coloffset=widthD+1;
        }
        if(options->l2tile)    {
        	coloffset=2*(widthD+1);
        }

        int rowoffset=coloffset;

    	IF_DEBUG2(printf("[SICA] column offset: %i\n", coloffset););
    	IF_DEBUG2(printf("[SICA] row offset: %i\n\n", rowoffset););

    	////malloc the sicadata->trans matrices and fill it
    	act_band->sicadata->trans=(int**)malloc(act_band->sicadata->transwidth*sizeof(int*));
    	for(x=0; x < act_band->sicadata->transwidth; x++)    {
    		act_band->sicadata->trans[x]=(int*)malloc(act_band->sicadata->transwidth*sizeof(int));
    	}

    	act_band->sicadata->trans_inverted=(int**)malloc(act_band->sicadata->transwidth*sizeof(int*));
    	for(x=0; x < act_band->sicadata->transwidth; x++)    {
    		act_band->sicadata->trans_inverted[x]=(int*)malloc(act_band->sicadata->transwidth*sizeof(int));
    	}

    	//fill it
    	int addrowoffset=0;
    	for(y=0; y < act_band->sicadata->transwidth; y++)    {
    		/* if the row is a scalar dimension or a zero row, skip it (so if it looks like this (0,0,...,0,?))
    		 * so if the sum of all elements but the last is 0
    		 */
    		int sum = 0;
    		for(x=0;x<act_band->loop->stmts[s]->trans->ncols-1;x++)    {
    			sum += act_band->loop->stmts[s]->trans->val[y][x];
    		}
    		if(sum==0)    {
    			addrowoffset++;
    		IF_DEBUG2(printf("[SICA] Skipping row %i\n",y););
    		}
	   	    for(x=0; x < act_band->sicadata->transwidth; x++)    {
	    		act_band->sicadata->trans[y][x]=act_band->loop->stmts[s]->trans->val[y+rowoffset+addrowoffset][x+coloffset];
	   	    }
   	    }

    	// [SICA] invert it
    	sica_inverse(act_band->sicadata->trans, act_band->sicadata->trans_inverted, act_band->sicadata->transwidth);
      }
    }
	// [SICA] STOP extract the transformation matrices before prevectorize
}
