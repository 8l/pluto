/*
 * sica_retile.h
 *
 *  Created on: 19.02.2013
 *      Author: dfeld
 */

#ifndef SICA_RETILE_H
#define SICA_RETILE_H

void sica_retile_band(PlutoProg *prog, Band *band, int offset);
void sica_retile_scattering_dims(PlutoProg *prog, Band **bands, int nbands, int l2);

#endif /* SICA_RETILE_H */
