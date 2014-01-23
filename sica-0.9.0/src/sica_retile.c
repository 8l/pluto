/*
 * sica_retile.c
 *
 *  Created on: 19.02.2013
 *      Author: dfeld
 */

#include <stdio.h>
#include <assert.h>

#include "pluto.h"
#include "sica_post_transform.h"
#include "program.h"
#include "transforms.h"

#include "sica_retile.h"

/* Manipulates statement domain and transformation to tile scattering
 * dimensions from firstD to lastD */
void sica_retile_band(PlutoProg *prog, Band *band, int offset)
{
    IF_DEBUG(      printf("[SICA] Performing ");      );
	IF_DEBUG(if(band->sicadata->isvec){printf("VECTORIZED");}else{printf("NON-VECTORIZED");});
    IF_DEBUG(      printf(" sica_retile_band step\n");  );

    int j, s;
    int npar;

    npar = prog->npar;

    int firstD = band->loop->depth;
    int lastD = band->loop->depth+band->width-1;

    IF_DEBUG2(printf("\n\n[SICA] firstD=%i lastD=%i for the following transformation\n",firstD,lastD););

    for (s=0; s<band->loop->nstmts; s++) {

        int tile_sizes[prog->num_hyperplanes];

        /* [SICA] START SIMD-specific tile sizes */
        for (j=0; j<prog->num_hyperplanes; j++)   {
            tile_sizes[j] = 1;
        }
    	if(band->sicadata->isvec)    {/* only set real retile sizes if this loop is vectorized */
    		if((!offset)&&(options->l2tile))    {
            tile_sizes[0]=band->sicadata->sical2size;
    		} else {
            tile_sizes[band->sicadata->vecrow]=band->sicadata->sical1size[s];
    		}
    	}
        /* [SICA] STOP SIMD-specific tile sizes */



//    	int skipped_scalar_dims=0;

        int sica_dimensions2tile=band->sicadata->tilewidth[s];
    	if(options->l2tile)    {
        	sica_dimensions2tile=(band->sicadata->tilewidth[s])/2;
    	}

	/* [SICA] NEW APPROACH */
	Stmt *stmt = band->loop->stmts[s];

		int row_offset = -2*sica_dimensions2tile*offset;
		int col_offset = sica_dimensions2tile*offset;

		IF_DEBUG2(printf("[SICA] before re-tiling\n"); );
        IF_DEBUG2(pluto_constraints_print(stdout, stmt->domain);  );

int dim;
for(dim=0;dim<sica_dimensions2tile;dim++)    {
	IF_DEBUG2(printf("dim: %i, line: %i\n", dim, sica_dimensions2tile-dim););

	IF_DEBUG2(printf("(%i,%i) orig: %i, new: %i\n",(stmt->domain->nrows-1)-(2*dim)+row_offset,dim+col_offset, stmt->domain->val[(stmt->domain->nrows-1)-(2*dim)+row_offset][dim+col_offset],tile_sizes[sica_dimensions2tile-dim-1]););
	IF_DEBUG2(printf("(%i,%i) orig: %i, new: %i\n",(stmt->domain->nrows-1)-(2*dim)+row_offset,stmt->domain->ncols-1,stmt->domain->val[(stmt->domain->nrows-1)-(2*dim)+row_offset][stmt->domain->ncols-1], band->sicadata->upperboundoffset[s]+tile_sizes[sica_dimensions2tile-dim-1]-1););
	IF_DEBUG2(printf("(%i,%i) orig: %i, new: %i\n",(stmt->domain->nrows-1)-(2*dim)-1+row_offset,dim+col_offset,stmt->domain->val[(stmt->domain->nrows-1)-(2*dim)-1+row_offset][dim+col_offset],-tile_sizes[sica_dimensions2tile-dim-1]););
                /* upper bound of the tile-loop  X * tx + Y >= 0 (the relating outer tile-loop is stored in column depth+vecrow for l1 and depth+vecrow for l2) */
                /* X */
                IF_DEBUG2(printf("[SICA] Accessing line %i\n",(stmt->domain->nrows-1)-(2*sica_dimensions2tile)+row_offset); );
                stmt->domain->val[(stmt->domain->nrows-1)-(2*dim)+row_offset][dim+col_offset] = tile_sizes[sica_dimensions2tile-dim-1];
                /* Y */
                /* ToDo: is it correct that in stmt->trans->val there is NO OFFSET?!? LIKE THAT IT WORKS -> figure out the relation between trans element and domain value */
                stmt->domain->val[(stmt->domain->nrows-1)-(2*dim)+row_offset][stmt->domain->ncols-1] = band->sicadata->upperboundoffset[s]+tile_sizes[sica_dimensions2tile-dim-1]-1;
                //IF_DEBUG(printf("[SICA] Statement %i -> Backuped upper bound offset: %i, current offset would be %i\n", s, band->sicadata->upperboundoffset[s], -stmt->trans->val[(depth-firstD)+1+depth][stmt->dim+prog->npar]););

                /* lower bound of the tile-loop  X * tx >= 0 (the relating outer tile-loop is stored in column depth+vecrow for l1 and depth+vecrow for l2) */
                /* X */
                stmt->domain->val[(stmt->domain->nrows-1)-(2*dim)-1+row_offset][dim+col_offset] = -tile_sizes[sica_dimensions2tile-dim-1];
}
		IF_DEBUG2(printf("[SICA] after re-tiling\n"););
		IF_DEBUG2(pluto_constraints_print(stdout, stmt->domain););

        } /* all statements */

    /* Re-detect hyperplane types (H_SCALAR, H_LOOP) */
    pluto_detect_hyperplane_types(prog);

}

/* Tiles scattering functions for all bands; l2=1 => perform l2 tiling */
void sica_retile_scattering_dims(PlutoProg *prog, Band **bands, int nbands, int l2)
{
//	printf("[SICA] Performing sica_retile_scattering_dims step\n");

    int b=0;

    for (b=0; b<nbands; b++) {

        IF_DEBUG(printf("[SICA] SICA modification for band number %i\n",b); );

        if (l2) {
            sica_retile_band(prog, bands[b], 0);
        }else{
        	if(!options->l2tile)    {
            sica_retile_band(prog, bands[b], 0);
        	}else{
            sica_retile_band(prog, bands[b], 1);
        	}
        }
    } /* all bands */
    /* Re-detect hyperplane types (H_SCALAR, H_LOOP) */
    pluto_detect_hyperplane_types(prog);
    pluto_detect_hyperplane_types_stmtwise(prog);

    // print_hyperplane_properties(prog);
    // pluto_transformations_pretty_print(prog);
}

