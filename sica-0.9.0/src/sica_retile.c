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
void sica_retile_band(PlutoProg *prog, Band *band, int *tile_sizes, int offset)
{
    IF_DEBUG(      printf("[SICA] Performing ");      );
	IF_DEBUG(if(band->sicadata->isvec){printf("VECTORIZED");}else{printf("VECTORIZED");});
    IF_DEBUG(      printf(" sica_retile_band step\n");  );
	
    int s;
    int depth, npar;

    npar = prog->npar;

    int firstD = band->loop->depth;
    int lastD = band->loop->depth+band->width-1;
    
    int widthD=lastD-firstD;

    IF_DEBUG2(printf("\n\n[SICA] firstD=%i lastD=%i for the following transformation\n",firstD,lastD););
    
    for (s=0; s<band->loop->nstmts; s++) {
        for (depth=firstD; depth<=lastD; depth++)    {
            Stmt *stmt = band->loop->stmts[s];
            
            int hyp_type = (stmt->hyp_types[depth + depth - firstD] == H_SCALAR)? H_SCALAR: 
                H_TILE_SPACE_LOOP;

            // pluto_constraints_print(stdout, stmt->domain);
            if (hyp_type != H_SCALAR) {
                assert(tile_sizes[depth-firstD] >= 1);
                
                /* add offset on l1 tiling if l2 tiling will be performed */
                int row_offset = -2*(widthD+1)*offset;
                int col_offset = (widthD+1)*offset;
                                
                /* upper bound of the tile-loop  X * tx + Y >= 0 (the relating outer tile-loop is stored in column depth+vecrow for l1 and depth+vecrow for l2) */
                /* X */
                stmt->domain->val[(stmt->domain->nrows-1)-(2*(depth-firstD))+row_offset][depth-firstD+col_offset] = tile_sizes[widthD-(depth-firstD)];
                /* Y */
                /* ToDo: is it correct that in stmt->trans->val there is NO OFFSET?!? LIKE THAT IT WORKS -> figure out the relation between trans element and domeain value */
                stmt->domain->val[(stmt->domain->nrows-1)-(2*(depth-firstD))+row_offset][stmt->domain->ncols-1] = band->sicadata->upperboundoffset[s]+tile_sizes[widthD-(depth-firstD)]-1;
                //IF_DEBUG(printf("[SICA] Statement %i -> Backuped upper bound offset: %i, current offset would be %i\n", s, band->sicadata->upperboundoffset[s], -stmt->trans->val[(depth-firstD)+1+depth][stmt->dim+prog->npar]););

                /* lower bound of the tile-loop  X * tx >= 0 (the relating outer tile-loop is stored in column depth+vecrow for l1 and depth+vecrow for l2) */
                /* X */
                stmt->domain->val[(stmt->domain->nrows-1)-(2*(depth-firstD))-1+row_offset][depth-firstD+col_offset] = -tile_sizes[widthD-(depth-firstD)];
                
                //IF_DEBUG(printf("\n[SICA] MODIFIED: [%i,%i] with tile_sizes[%i]=%i\n", (stmt->domain->nrows-1)-(2*(depth-firstD)), depth-firstD, widthD-(depth-firstD), tile_sizes[widthD-(depth-firstD)]);          );
                //IF_DEBUG(printf("[SICA] MODIFIED: [%i,%i] with tile_sizes[%i]=%i\n", (stmt->domain->nrows-1)-(2*(depth-firstD))-1, depth-firstD, widthD-(depth-firstD), tile_sizes[widthD-(depth-firstD)]);          );
                //IF_DEBUG(printf("[SICA] MODIFIED: [%i,%i] with tile_sizes[%i]=%i\n", (stmt->domain->nrows-1)-(2*(depth-firstD))-1, stmt->domain->ncols-1, widthD-(depth-firstD), tile_sizes[widthD-(depth-firstD)]); );
                
                // printf("Stmt %d: depth: %d\n", stmt->id+1,depth);
                // pluto_matrix_print(stdout, stmt->trans);

            }
            
            IF_DEBUG2(printf("[SICA] after re-tiling\n");            );
            IF_DEBUG2(pluto_constraints_print(stdout, stmt->domain); );
        } /* all statements */
    } /* all scats to be tiled */

    /* Re-detect hyperplane types (H_SCALAR, H_LOOP) */
    pluto_detect_hyperplane_types(prog);
    

	///////////////////////////////////////////////////////////////////////
}

/* Tiles scattering functions for all bands; l2=1 => perform l2 tiling */
void sica_retile_scattering_dims(PlutoProg *prog, Band **bands, int nbands, int l2)
{
//	printf("[SICA] Performing sica_retile_scattering_dims step\n");
	
    int i, j, b;
    int tile_sizes[prog->num_hyperplanes];
    int l2_tile_size_ratios[prog->num_hyperplanes];

    Stmt **stmts = prog->stmts;

    for (b=0; b<nbands; b++) {
        /* [SICA] START SIMD-specific tile sizes */
        IF_DEBUG(printf("[SICA] SICA modification for band number %i\n",b); );
        for (j=0; j<prog->num_hyperplanes; j++)   {
            tile_sizes[j] = 1;
            l2_tile_size_ratios[j] = 1;
        }
    	if(bands[b]->sicadata->isvec)    {/* only set real retile sizes if this loop is vectorized */
            tile_sizes[bands[b]->sicadata->vecrow]=bands[b]->sicadata->sical1size;
            l2_tile_size_ratios[0]=bands[b]->sicadata->sical2size;
    	}
        /* [SICA] STOP SIMD-specific tile sizes */
    	
        if (l2) {
            sica_retile_band(prog, bands[b], l2_tile_size_ratios, 0);
        }else{
        	if(!options->l2tile)    {
            sica_retile_band(prog, bands[b], tile_sizes, 0);
        	}else{
            sica_retile_band(prog, bands[b], tile_sizes, 1);
        	}
        }
    } /* all bands */
    /* Re-detect hyperplane types (H_SCALAR, H_LOOP) */
    pluto_detect_hyperplane_types(prog);
    pluto_detect_hyperplane_types_stmtwise(prog);

    // print_hyperplane_properties(prog);
    // pluto_transformations_pretty_print(prog);
}

