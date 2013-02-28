/*
 * sica.h
 *
 *  Created on: 19.02.2013
 *      Author: dfeld
 */

#ifndef SICA_H
#define SICA_H

#define SICA_STRING_SIZE 32

#include "pluto.h"

/* [SICA] structure for SICA necessary data */
struct sica_data{
	int isvec;
	
    int vecloop; //vectorized loop t%i
    int vecrow; //row in tile_sizes related to vectorized loop t%i
    
    int sical1size; //l1 tile size for vectorized loop
    int sical2size; //l2 tile size for outermost loop
    
    int *upperboundoffset; //offsetarray for each statement for retiling upper bound (this depends on transformation matrix),
                          //as this is modified in the process, this backup is needed

    char **id2arrayname; //array that stores the names of all arrays to specify a unique id for each
    int nb_arrays;
};
typedef struct sica_data SICAData;

#endif /* SICA_H */
