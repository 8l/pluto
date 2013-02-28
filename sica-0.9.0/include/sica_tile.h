/*
 * sica_tile.h
 *
 *  Created on: 19.02.2013
 *      Author: dfeld
 */

#ifndef SICA_TILE_H
#define SICA_TILE_H

void sica_tile_band(PlutoProg *prog, Band *band, int *tile_sizes);
void sica_tile(PlutoProg *prog, scoplib_scop_p);
void sica_tile_scattering_dims(PlutoProg *prog, Band **bands, int nbands, int l2);

#endif /* SICA_TILE_H */
