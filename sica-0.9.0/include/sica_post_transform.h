/*
 * sica_post_transform.h
 *
 *  Created on: 19.02.2013
 *      Author: dfeld
 */

#ifndef SICA_POST_TRANSFORM_H
#define SICA_POST_TRANSFORM_H

#include "pluto.h"
#include "sica.h"

int sica_pre_vectorize_band(Band *band, int num_tiling_levels, PlutoProg *prog);
int sica_pre_vectorize(PlutoProg *prog);

/* [SICA] visible from orginal post_transform.c */
int get_num_spatial_accesses(Ploop *loop, PlutoProg *prog);
int get_num_invariant_accesses(Ploop *loop, PlutoProg *prog);
int get_num_accesses(Ploop *loop, PlutoProg *prog);

#endif /* SICA_POST_TRANSFORM_H */
