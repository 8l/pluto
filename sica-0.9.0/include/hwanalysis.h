/*
 * hwanalysis.h
 *
 *  Created on: 31.01.2013
 *      Author: dfeld
 */

#ifndef HWANALYSIS_H_
#define HWANALYSIS_H_

#define CPU_ID         0
#define CPU_NAME       1
#define SSE_ID         2
#define SSE_NAME       3

#define L1CACHE_SIZE   4
#define L1CACHE_ASSO   5
#define L1CACHE_LINE   6

#define L2CACHE_SIZE   7
#define L2CACHE_ASSO   8
#define L2CACHE_LINE   9

#define L3CACHE_SIZE  10
#define L3CACHE_ASSO  11
#define L3CACHE_LINE  12

#include "printconvert.h"

int get_hardware_cache_infos(int );

#endif /* HWANALYSIS_H_ */

