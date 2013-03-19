/*
 * sica_tile.c
 *
 *  Created on: 19.02.2013
 *      Author: dfeld
 */

#include <stdio.h>
#include <assert.h>

#include <string.h>

#include "pluto.h"
#include "sica_post_transform.h"
#include "program.h"
#include "transforms.h"

#include "hwanalysis.h"

#include "sica_tile.h"
#include "sica_retile.h"
#include "sica_accesses.h"
#include "sica_func.h"
#include "sica_tilesizes.h"
#include "sica_math_func.h"

/* Manipulates statement domain and transformation to tile scattering 
 * dimensions from firstD to lastD */
void sica_tile_band(PlutoProg *prog, Band *band, int *tile_sizes)
{
    int i, j, s;
    int depth, npar;

    Stmt **stmts = prog->stmts;
    npar = prog->npar;

    int firstD = band->loop->depth;
    int lastD = band->loop->depth+band->width-1;

    for (depth=firstD; depth<=lastD; depth++)    {
        for (s=0; s<band->loop->nstmts; s++) {
            Stmt *stmt = band->loop->stmts[s];
            /* 1. Specify tiles in the original domain. 
             * NOTE: tile shape info comes in here */

            /* 1.1 Add additional dimensions */
            char iter[6];
            sprintf(iter, "zT%d", stmt->dim);

            int hyp_type = (stmt->hyp_types[depth + depth - firstD] == H_SCALAR)? H_SCALAR: 
                H_TILE_SPACE_LOOP;

            /* 1.2 Specify tile shapes in the original domain */
            // pluto_constraints_print(stdout, stmt->domain);
            int num_domain_supernodes = 0;
            if (hyp_type != H_SCALAR) {
                assert(tile_sizes[depth-firstD] >= 1);
		
		/* [SICA] count the really tiled dimensions */
		band->sicadata->tilewidth[s]++;

                /* Domain supernodes aren't added for scalar dimensions */
                // printf("S%d dim: %d %d\n", stmt->id+1, stmt->dim, depth-firstD);
                pluto_stmt_add_dim(stmt, num_domain_supernodes, depth, iter, hyp_type, prog);
                /* Add relation b/w tile space variable and intra-tile variables like
                 * 32*xt <= 2t+i <= 32xt + 31 */
                /* Lower bound */
                pluto_constraints_add_inequality(stmt->domain);

                for (j=num_domain_supernodes+1; j<stmt->dim+npar; j++) {
                    stmt->domain->val[stmt->domain->nrows-1][j] = 
                        stmt->trans->val[firstD+(depth-firstD)+1+(depth-firstD)][j];
                }

                stmt->domain->val[stmt->domain->nrows-1][num_domain_supernodes] = 
                    -tile_sizes[depth-firstD];

                stmt->domain->val[stmt->domain->nrows-1][stmt->domain->ncols-1] = 
                    stmt->trans->val[(depth-firstD)+1+depth][stmt->dim+prog->npar];

                PlutoConstraints *lb = pluto_constraints_select_row(stmt->domain, 
                        stmt->domain->nrows-1);
                pluto_update_deps(stmt, lb, prog);
                pluto_constraints_free(lb);

                /* Upper bound */
                pluto_constraints_add_inequality(stmt->domain);
                for (j=num_domain_supernodes+1; j<stmt->dim+npar; j++) {
                    stmt->domain->val[stmt->domain->nrows-1][j] = 
                        -stmt->trans->val[firstD+(depth-firstD)+1+(depth-firstD)][j];
                }

                stmt->domain->val[stmt->domain->nrows-1][num_domain_supernodes] 
                    = tile_sizes[depth-firstD];

                stmt->domain->val[stmt->domain->nrows-1][stmt->domain->ncols-1] = 
                    -stmt->trans->val[(depth-firstD)+1+depth][stmt->dim+prog->npar] 
                    +tile_sizes[depth-firstD]-1;
                
                /* [SICA] store the upper bound offset for the retiling (as this value can be changed in the intermediate steps) */
                /* [SICA] ToDo: Check wheather it is OK not to distinguish between L1 and L2 tiling but just taking the last (overwritten) value */ 
                band->sicadata->upperboundoffset[s]=-stmt->trans->val[(depth-firstD)+1+depth][stmt->dim+prog->npar];
                printf("Setting upper bound offset to %i\n", band->sicadata->upperboundoffset[s]);

                PlutoConstraints *ub = pluto_constraints_select_row(stmt->domain,
                        stmt->domain->nrows-1);
                pluto_update_deps(stmt, ub, prog);
                pluto_constraints_free(ub);

                num_domain_supernodes++;

                 IF_DEBUG2(printf("after adding tile constraints\n");     );
                 IF_DEBUG2(pluto_constraints_print(stdout, stmt->domain); );

                // printf("Stmt %d: depth: %d\n", stmt->id+1,depth);
                // pluto_matrix_print(stdout, stmt->trans);

            }else{
                /* Scattering function for tile space iterator is set the
                 * same as its associated domain iterator  
                 * Dimension is not a loop; tile it trivially
                 */
                pluto_stmt_add_hyperplane(stmt, H_SCALAR, depth);
                for (j=0; j<stmt->dim+npar+1; j++) {
                    stmt->trans->val[depth][j] = 
                        stmt->trans->val[firstD+(depth-firstD)+1+(depth-firstD)][j];
                }
            }
            stmt->num_tiled_loops++;
            stmt->first_tile_dim = firstD ;
            stmt->last_tile_dim = lastD;
        } /* all statements */
    } /* all scats to be tiled */

    int max = 0, curr;
    for (i=0; i<prog->nstmts; i++) {
        max = PLMAX(stmts[i]->trans->nrows, max);
    }
    for (i=0; i<prog->nstmts; i++) {
        curr = stmts[i]->trans->nrows;
        for (j=curr; j < max; j++) {
            pluto_sink_transformation(stmts[i], stmts[i]->trans->nrows, prog);
        }
    }

    // print_hyperplane_properties(prog);
    curr = prog->num_hyperplanes;
    for (depth=curr; depth<max; depth++)    {
        pluto_prog_add_hyperplane(prog, depth, H_UNKNOWN);
    }
    /* Re-detect hyperplane types (H_SCALAR, H_LOOP) */
    pluto_detect_hyperplane_types(prog);

    // print_hyperplane_properties(prog);
    // pluto_transformations_pretty_print(prog);
}


/* Updates the statement domains and transformations to represent the new
 * tiled code. A schedule of tiles is created for parallel execution if
 * --parallel is on 
 *
 *  Pre-vectorization is also done inside a tile
 *
 *  */
void sica_tile(PlutoProg *prog)
{

    IF_DEBUG(print_all_cache_information());

    SICAHardware* sica_hardware;
    sica_hardware=(SICAHardware*)malloc(sizeof(SICAHardware));
    sica_hardware->regsize = id2regsize(get_hardware_cache_infos(SSE_ID));
    sica_hardware->l1cachesize = get_hardware_cache_infos(L1CACHE_SIZE);
    sica_hardware->l2cachesize = get_hardware_cache_infos(L2CACHE_SIZE);
    sica_hardware->ratio=1.0;

    print_l1cache_hierarchy(sica_hardware->regsize, sica_hardware->l1cachesize);

    if(options->l2tile)
    {
        print_addl2cache_hierarchy(sica_hardware->l2cachesize);
    }
    printf("\n");


    int nbands, i;
    Band **bands;
    bands = pluto_get_outermost_permutable_bands(prog, &nbands);
    IF_DEBUG(printf("Outermost tilable bands\n"););
    IF_DEBUG(pluto_bands_print(bands, nbands););
    
    /* [SICA] allocate and initialize the SICAData memory on each band */
    sica_malloc_and_init_sicadata(bands, nbands);
    
    /* Now, we are ready to tile */
    if (options->lt >= 0 && options->ft >= 0)   {
        /* User option specified tiling */

        assert(options->ft <= prog->num_hyperplanes-1);
        assert(options->lt <= prog->num_hyperplanes-1);
        assert(options->ft <= options->lt);

        /* L1 tiling */
        sica_tile_scattering_dims(prog, bands, nbands, 0);

        if (options->l2tile)    {
            sica_tile_scattering_dims(prog, bands, nbands, 1);
        }
    }else{
        /* L1 tiling */
        sica_tile_scattering_dims(prog, bands, nbands, 0);
        if (options->l2tile)    {
            /* L2 tiling */
            sica_tile_scattering_dims(prog, bands, nbands, 1);
        }
    }

    if (options->intratileopt) {
        int retval = 0;
        for (i=0; i<nbands; i++) {
            retval |= pluto_intra_tile_optimize_band(bands[i], 1, prog); 
        }
        if (retval) pluto_detect_transformation_properties(prog);
        if (retval && !options->silent) {
            printf("After intra_tile_opt\n");
            pluto_transformations_pretty_print(prog);
        }
    }

    /* Detect properties again after tiling */
    pluto_detect_transformation_properties(prog);

    /* [SICA] get the transformation matrices BEFORE prevector transformation */
    pluto_transformations_pretty_print(prog);

    sica_get_trans_matrix(bands, nbands);
    sica_print_fuse_structure(bands, nbands);

    if (options->prevector) {
        int retval = 0;
        for (i=0; i<nbands; i++) {
            int num_tiling_levels = options->tile + options->l2tile;
            retval |= sica_pre_vectorize_band(bands[i], num_tiling_levels, prog); 
        }
        if (retval) pluto_detect_transformation_properties(prog);
        if (retval && !options->silent) {
            printf("After pre_vectorize:\n");
            pluto_transformations_pretty_print(prog);
        }
    }

    /* [SICA] set the bands either to be vectorized or not */
    sica_set_vectorized_bands(bands, nbands);

    /* [SICA] calculate the band specific tile quantities */
    for (i=0; i<nbands; i++) {

    	Band* act_band=bands[i];

		printf("[SICA] tiling band %i\n",i);

    	sica_get_band_specific_tile_sizes(act_band);

    	if(act_band->sicadata->isvec)    {
    	    //printf("VEC\tbands[%i], nstmts=%i\n", i, act_band->loop->nstmts);
    		if(act_band->sicadata->vec_accesses>0)    {
        		IF_DEBUG(printf("[SICA] vectorized accesses in this band: %i\n", act_band->sicadata->vec_accesses););
            	printf("[SICA] percentage of INNERMOST vectorized accesses: %.2f\n", 100.0*(float)act_band->sicadata->innermost_vec_accesses/(float)act_band->sicadata->vec_accesses);
            	IF_DEBUG(printf("[SICA] bytes to be loaded by the vectorized accesses: %i Bytes\n", act_band->sicadata->bytes_per_vecit););

    		act_band->sicadata->sical1size=sica_get_l1size(act_band->sicadata, sica_hardware); // [SICA] HERE A FUNCTION SHOULD BE CALLED THAT CALCULATES THE SICA SIZES FOR THAT BAND
    		act_band->sicadata->sical2size=sica_get_l2size(sica_hardware);  // [SICA] HERE A FUNCTION SHOULD BE CALLED THAT CALCULATES THE GLOBAL SIZE

    	    printf("[SICA] tile sizes for band %i -> Level-1: %i, Level-2: %i\n\n",i,act_band->sicadata->sical1size, act_band->sicadata->sical2size );

    		} else {
        		act_band->sicadata->sical1size=1;
        		act_band->sicadata->sical2size=1;

            	printf("[SICA] NO vectorized access\n");
        	    printf("[SICA] tile sizes for band %i -> Level-1: %i, Level-2: %i\n\n",i,act_band->sicadata->sical1size, act_band->sicadata->sical2size );
    		}
    	}
    }

    /* [SICA] Modify the tile sizes by SICA approach START */
    sica_retile_scattering_dims(prog, bands, nbands, 0); /* L1 tiling */
    if (options->l2tile)    {
        sica_retile_scattering_dims(prog, bands, nbands, 1); /* L2 tiling */
    }

    /* [SICA] Free the SICAData memory */

    if (options->parallel) {
        create_tile_schedule(prog, bands, nbands);
    }

    /* [SICA] free the SICAData memory on each band */
    //sica_free_sicadata(bands, nbands);

    pluto_bands_free(bands, nbands);
}




/* Tiles scattering functions for all bands; l2=1 => perform l2 tiling */
void sica_tile_scattering_dims(PlutoProg *prog, Band **bands, int nbands, int l2)
{
    int i, j, b;
    int depth;
    int tile_sizes[prog->num_hyperplanes];
    int l2_tile_size_ratios[prog->num_hyperplanes];

    Stmt **stmts = prog->stmts;

    /* [SICA] This is just the default tiling to enable potential vectorization in all dimensions*/
    /* [SICA] Nevertheless, the values can simply be the default ones for this (every positive value but '1' should work) */
    /* [SICA] As these values are overwritten anyway and to avoid any problems with possibly changed or disrupted DEFAULT_L1_TILE_SIZE it is hard-coded */
    for (j=0; j<prog->num_hyperplanes; j++)   {
        tile_sizes[j] = 32;//DEFAULT_L1_TILE_SIZE;
        l2_tile_size_ratios[j] = 8;
    }

    for (b=0; b<nbands; b++) {
        if (l2) {
            sica_tile_band(prog, bands[b], l2_tile_size_ratios);
        }else{
            sica_tile_band(prog, bands[b], tile_sizes);
        }
    } /* all bands */

    int max = 0, curr;
    for (i=0; i<prog->nstmts; i++) {
        max = PLMAX(stmts[i]->trans->nrows, max);
    }
    for (i=0; i<prog->nstmts; i++) {
        curr = stmts[i]->trans->nrows;
        for (j=curr; j < max; j++) {
            pluto_sink_transformation(stmts[i], stmts[i]->trans->nrows, prog);
        }
    }

    // print_hyperplane_properties(prog);
    curr = prog->num_hyperplanes;
    for (depth=curr; depth<max; depth++)    {
        pluto_prog_add_hyperplane(prog, depth, H_UNKNOWN);
    }
    /* Re-detect hyperplane types (H_SCALAR, H_LOOP) */
    pluto_detect_hyperplane_types(prog);
    pluto_detect_hyperplane_types_stmtwise(prog);

    // print_hyperplane_properties(prog);
    // pluto_transformations_pretty_print(prog);
}
