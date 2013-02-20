/*
 * sica_post_transform.c
 *
 *  Created on: 19.02.2013
 *      Author: dfeld
 */

#include <stdio.h>
#include <assert.h>

#include "pluto.h"
#include "sica_post_transform.h"
#include "post_transform.h"
#include "program.h"
#include "transforms.h"

/* Vectorize first loop in band that meets criteria */
int sica_pre_vectorize_band(Band *band, int num_tiling_levels, PlutoProg *prog)
{
    int num, l;

    /* Band has to be the innermost band as well */
    if (!pluto_is_band_innermost(band, num_tiling_levels)) return 0;

    Ploop **loops;

    loops = pluto_get_loops_under(band->loop->stmts, band->loop->nstmts, 
            band->loop->depth + num_tiling_levels*band->width, prog, &num);

    for (l=0; l<num; l++) {
        if (!pluto_loop_is_parallel(prog, loops[l])) continue;
        int s, t, a;
        a = get_num_accesses(loops[l], prog);
        s = get_num_spatial_accesses(loops[l], prog);
        t = get_num_invariant_accesses(loops[l], prog);
        /* Vectorize only if each access has either spatial or temporal
         * reuse */
        /* if accesses haven't been provided, a would be 0 */
        if (a >= 1 && a == s + t) break;
    }

    if (l < num) {
        /* [SICA] Store information on loop that is vectorized in this band*/
        band->sicadata->vecloop=loops[l]->depth;
        band->sicadata->vecrow=l;

        pluto_make_innermost(loops[l], prog);
        IF_DEBUG(printf("[Pluto] Loop to be vectorized: "););
        IF_DEBUG(pluto_loop_print(loops[l]););
        return 1;
    }

    return 0;
}


int sica_pre_vectorize(PlutoProg *prog)
{
    int nbands, i;
    Band **bands;
    bands = pluto_get_outermost_permutable_bands(prog, &nbands);
    int retval = 0;
    for (i=0; i<nbands; i++) {
        retval |= sica_pre_vectorize_band(bands[i], 0, prog); 
    }
    if (retval) pluto_transformations_pretty_print(prog);
    pluto_bands_free(bands, nbands);
    return 0;
}
