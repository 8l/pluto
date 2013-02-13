/*
 * intelspeccache.h
 *
 *  Created on: 31.01.2013
 *      Author: dfeld
 */

#ifndef INTELSPECCACHE_H_
#define INTELSPECCACHE_H_

void set_l3_cache_parameter(int, int *);

void set_l1_cache_parameter(int , int *);

int get_intel_specific_cache_info(char * , char * , char * , char * , int , int , int * , int * , int * );

#endif /* INTELSPECCACHE_H_ */
