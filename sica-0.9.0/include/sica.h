/*
 * sica.h
 *
 *  Created on: 19.02.2013
 *      Author: dfeld
 */

#ifndef SICA_H
#define SICA_H

#define ACCESS_IS_NOT_IDENTICAL 0
#define ACCESS_IS_IDENTICAL 1
#define ACCESS_IS_IN_SAME_STRIDE 2

#define SICA_STRING_SIZE 32

struct sica_matrix{
	int **val;
};
typedef struct sica_matrix SICAMatrix;

/* [SICA] structure for SICA necessary data */
struct sica_data{
	int isvec;					//is the band vectorizable?
	
    int vecloop; 				//vectorized loop t%i
    int vecrow; 				//row in tile_sizes related to vectorized loop t%i
    
    int sical1size; 			//l1 tile size for vectorized loop
    int sical2size; 			//l2 tile size for outermost loop
    
    int *upperboundoffset; 		//offsetarray for each statement for retiling upper bound (this depends on transformation matrix),
								//as this is modified in the process, this backup is needed

    char **id2arrayname; 		//array that stores the names of all arrays to specify a unique id for each
    int nb_arrays;

    SICAMatrix *trans;				//pointer to transformation matrices for each statement in this band, trans[s]->mat[][]
    SICAMatrix *trans_inverted;		//inverted transformation matrix
    int *transwidth;			//width of the quadratic transformation matrix TODO: is it always quadratic

    int *tilewidth;			//number of tiled dimensions sometimes != stms->num_tiled_loops (because of scalar dims);

    int vec_accesses;			//number of counted (different) accesses by the vectorized loop in this band
    int innermost_vec_accesses;	//number of counted (different) INNERMOST accesses by the vectorized loop in this band
    int bytes_per_vecit;		//number of bytes that have to be loaded per iteration of the vectorized loop within this band

    int largest_data_type;		/*this value is/might be interesting for mixed precision calculations because the different
								 *SIMD registers should be seperated in pieces of this datatype in such a case
								 * EXAMPLE: this calculation should split the float accesses also in portions of 4 to fit (instead of possible 8)
								 *
								 * float  |----|    |----|    |----|    |----|    |
								 *                         +
								 * double |---- ----|---- ----|---- ----|---- ----|
								 *                         =
								 * double |---- ----|---- ----|---- ----|---- ----|
								 */

    int *coloffset;

    int *scalar_dims;			//array containing the scalar dimension for each statement
};
typedef struct sica_data SICAData;

#endif /* SICA_H */
