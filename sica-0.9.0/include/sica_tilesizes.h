/*
 * sica_tilesizes.h
 *
 *  Created on: 07.03.2013
 *      Author: dfeld
 */


#ifndef SICA_TILESIZES_H_
#define SICA_TILESIZES_H_

/* [SICA] structure for hardware information */
struct sica_hwinfo{
	int l1cachesize;
	int l2cachesize;
	int regsize;
	float ratio;
};
typedef struct sica_hwinfo SICAHardware;

#endif /* SICA_TILESIZES_H_ */
