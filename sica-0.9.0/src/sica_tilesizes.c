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

int sica_get_l1size(SICAData* sicadata, SICAHardware* sicahardware)    {

	int l1tilesize;

	//get the theoretical cache fitting quantity
	l1tilesize=(int)(sicahardware->ratio*(float)((sicahardware->l1cachesize*1024)/sicadata->bytes_per_vecit));

	//get the regsize-floor related to the largest available datatype in this band
	int regsizeinelements=(sicahardware->regsize/8)/(sicadata->largest_data_type);
	IF_DEBUG(printf("[SICA] Largest datatype for this band: %i\n", sicadata->largest_data_type););
	IF_DEBUG(printf("[SICA] %i elements fit to the register\n", regsizeinelements););
	l1tilesize=(int)(l1tilesize/regsizeinelements)*regsizeinelements;

	return l1tilesize;
}

int sica_get_l2size(SICAHardware* sicahardware)    {
	return sicahardware->l2cachesize/sicahardware->l1cachesize;
}

