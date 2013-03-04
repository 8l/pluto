/*
 * sica.h
 *
 *  Created on: 19.02.2013
 *      Author: dfeld
 */

#ifndef SICA_H
#define SICA_H

#define SICA_STRING_SIZE 32

//#include "sica_accesses.h"

/* [SICA] structure for SICA necessary data */
struct sica_data{
	int isvec;				//is the band vectorizable?
	
    int vecloop; 			//vectorized loop t%i
    int vecrow; 			//row in tile_sizes related to vectorized loop t%i
    
    int sical1size; 		//l1 tile size for vectorized loop
    int sical2size; 		//l2 tile size for outermost loop
    
    int *upperboundoffset; 	//offsetarray for each statement for retiling upper bound (this depends on transformation matrix),
							//as this is modified in the process, this backup is needed

    char **id2arrayname; 	//array that stores the names of all arrays to specify a unique id for each
    int nb_arrays;

    int **trans;			//pointer to transformation matrices for each statement in this band, trans[s]->mat[][]
    int **trans_inverted;	//inverted transformation matrix
    int transwidth;			//width of the quadratic transformation matrix TODO: is it always quadratic

//    SICAAccess** sica_accesses_on_array; //array that stores the relevant accesses for vectorization for each array in a linked list
										//usage: sica_accesses_on_array[arrayID]->...
};
typedef struct sica_data SICAData;

#endif /* SICA_H */
