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

void sica_get_l1size(SICAData* sicadata, SICAHardware* sicahardware, int nstmts);

int sica_get_l2size(SICAHardware* sicahardware);

#endif /* SICA_TILESIZES_H_ */
