/*
 * sica.h
 *
 *  Created on: 19.02.2013
 *      Author: dfeld
 */

#ifndef SICA_H
#define SICA_H

#include "pluto.h"

/* [SICA] structure for SICA necessary data */
struct sica_data{
	int isvec;
	
    int vecloop; //vectorized loop t%i
    int vecrow; //row in tile_sizes related to vectorized loop t%i
    
    int sical1size; //l1 tile size for vectorized loop
    int sical2size; //l2 tile size for outermost loop
};
typedef struct sica_data SICAData;

#endif /* SICA_H */
