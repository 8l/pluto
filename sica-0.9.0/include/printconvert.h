/*
 * printconvert.h
 *
 *  Created on: 31.01.2013
 *      Author: dfeld
 */

#include <string.h>
//#include <cstdlib>

#ifndef PRINTCONVERT_H_
#define PRINTCONVERT_H_

void print_all_cache_information();

void print_l1cache_hierarchie(int , int );

void print_addl2cache_hierarchie(int );

char* id2sse(int );

char* id2cpu(int );

int id2regsize(int );

#endif /* PRINTCONVERT_H_ */
