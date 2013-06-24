/*
 * sica_func.h
 *
 *  Created on: 06.03.2013
 *      Author: dfeld
 */


#ifndef SICA_FUNC_H_
#define SICA_FUNC_H_

#define SICA_DEFAULT_DATA_BYTES 4 //suggest a 4 byte data element for the vector registers as default if no one is specified

#include "pluto.h"

void sica_malloc_and_init_sicadata(Band **bands, int nbands);

void sica_free_sicadata(Band **bands, int nbands);

void sica_get_band_specific_tile_sizes(Band* act_band);

void sica_print_matrix_with_coloffset(int** matrix, int rows, int cols, int coloffset);

void sica_print_array_accesses_structures(Band* act_band, SICAAccess** sica_accesses_on_array);

void sica_print_fuse_structure(Band **bands, int nbands);

void sica_check_for_stmts_in_scalar_dim(int s, Band* act_band, SICAStmtList* sica_stmts_in_scalar_dim);

#endif /* SICA_FUNC_H_ */
