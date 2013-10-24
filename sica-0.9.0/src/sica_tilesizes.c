/*
 * sica_tilesizes.c
 *
 *  Created on: 07.03.2013
 *      Author: dfeld
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>

#include "pluto.h"

#include "hwanalysis.h"
#include "sica.h"
#include "sica_tilesizes.h"

void sica_get_l1size(SICAData* sicadata, SICAHardware* sicahardware, int nstmts)    {
int s=nstmts;//TEMP: RECENT STATEMENT
//	for(s=0; s<nstmts; s++)    {
	int l1tilesize;

	printf("[SICA] Largest datatype for this band: %i and statement %i\n", sicadata->largest_data_type[s],s);

	//get the theoretical cache fitting quantity
	l1tilesize=(int)(sicahardware->ratio*(float)((sicahardware->l1cachesize*1024)/sicadata->bytes_per_vecit[s]));

	//printf("DBG: sicadata->bytes_per_vecit[s]=%i\n", sicadata->bytes_per_vecit[s]);

	//get the regsize-floor related to the largest available datatype in this band
	int regsizeinelements=(sicahardware->regsize/8)/(sicadata->largest_data_type[s]);

	IF_DEBUG(printf("[SICA] %i elements fit to the register\n", regsizeinelements););

	l1tilesize=(int)(l1tilesize/regsizeinelements)*regsizeinelements;
//	if(sicadata->largest_data_type[s]<0)    {
//		//printf("There is not scalar dimension %i in this band\n", s);
//		l1tilesize=0;
//	}

	sicadata->sical1size[s]=l1tilesize;

//	}
}

int sica_get_l2size(SICAHardware* sicahardware)    {
	return sicahardware->l2cachesize/sicahardware->l1cachesize;
}

