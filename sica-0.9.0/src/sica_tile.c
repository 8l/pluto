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

#include "sica_tile.h"
#include "sica_retile.h"

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
void sica_tile(PlutoProg *prog, scoplib_scop_p scop)
{
    int nbands, i, s;
    Band **bands;
    bands = pluto_get_outermost_permutable_bands(prog, &nbands);
    IF_DEBUG(printf("Outermost tilable bands\n"););
    IF_DEBUG(pluto_bands_print(bands, nbands););
    
    /* [SICA] allocate and initialize the SICAData memory on each band */
    for (i=0; i<nbands; i++) {
        bands[i]->sicadata=(SICAData*)malloc(sizeof(SICAData));
        bands[i]->sicadata->isvec=-1;
        bands[i]->sicadata->vecloop=-1;
        bands[i]->sicadata->vecrow=-1;
        bands[i]->sicadata->sical1size=-1;
        bands[i]->sicadata->sical2size=-1;
        bands[i]->sicadata->upperboundoffset=(int*)malloc((bands[i]->loop->nstmts)*sizeof(int));
        for(s=0;s<bands[i]->loop->nstmts;s++)    {
        	bands[i]->sicadata->upperboundoffset[s]=-1;
        }
        bands[i]->sicadata->nb_arrays=0;
    }
    

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
    
    /* [SICA] Print the sicadata */
    for (i=0; i<nbands; i++) {
    	if(bands[i]->sicadata->vecloop>-1)
    	{
        bands[i]->sicadata->isvec=1;
        //ToDo: Here seems to arise an error that the loop index is always one too big (not the difference of starting counting at 0 instead of 1)
        IF_DEBUG(printf("[SICA] BAND %i, Original loop that is vectorized: t%i, row in tile_sizes: %i\n",i,bands[i]->sicadata->vecloop+1,bands[i]->sicadata->vecrow+1););
    	}else{
        bands[i]->sicadata->isvec=0;
    	}
    }
    
    /* [SICA] calculate the band specific tile quantities */
    for (i=0; i<nbands; i++) {
	    int a,s,r,w,x,y,t;

    	Band* act_band=bands[i];

		///////////////////////////////////////////////
	    //TEST THE NEW SICA

		/* [SICA] get an array<->id relation
		 *
		 * act_band->loop->stmts[s]->reads[r]->sym_id is not set in PluTo
		 * array with nb_arrays elements storing the name of the array
		 */

        /* [SICA] get number of all reads and writes as an upper bound of the number of arrays */
		int max_nb_arrays=0;
	    for(s=0; s<act_band->loop->nstmts;s++)
	    {
	    	max_nb_arrays+=act_band->loop->stmts[s]->nreads; /* number of total reads in this statement */
	    	max_nb_arrays+=act_band->loop->stmts[s]->nwrites; /* number of total writes in this statement */
	    }
	    printf("[SICA] max_nb_arrays=%i\n",max_nb_arrays);

		//[SICA] malloc the id2arrayname
		act_band->sicadata->id2arrayname=(char**)malloc(max_nb_arrays*sizeof(char*));
		for(a=0; a < max_nb_arrays; a++)    {
			act_band->sicadata->id2arrayname[a]=(char*)malloc(SICA_STRING_SIZE*sizeof(char));
        }

		//[SICA] init the id2arrayname
		for(a=0; a < max_nb_arrays; a++)    {
			for(t=0; t < SICA_STRING_SIZE; t++)    {
				act_band->sicadata->id2arrayname[a][t]='\0';
			}
        }

		/* [SICA] get the different arrays */
	    for(s=0; s<act_band->loop->nstmts;s++)    {
	    	/* [SICA] reads */
		    for(r=0;r<act_band->loop->stmts[s]->nreads;r++)    {
		    	/* [SICA] go through the array */
				int arraynameisnew=1; //store whether the prospected array name is already available or not
		    	a=0;
				while(strcmp(act_band->sicadata->id2arrayname[a],"\0"))    { //while there is a non-empty string in this array position
					//if there is one, check if it is equal to the one analysed at the moment
					if(!strcmp(act_band->sicadata->id2arrayname[a], act_band->loop->stmts[s]->reads[r]->name))    { //so they are identical
						arraynameisnew=0;
						break;
					}

					a++;
				}
				if(arraynameisnew) {
					//printf("[SICA] Copying the new name '%s' to the array to position %i\n", act_band->loop->stmts[s]->reads[r]->name, a);
					act_band->sicadata->id2arrayname[a]=strcpy(act_band->sicadata->id2arrayname[a], act_band->loop->stmts[s]->reads[r]->name);
					bands[i]->sicadata->nb_arrays++;
					//printf("[SICA] COPY SUCCEEDED\n");
				}
		    }



		    /* [SICA] writes */
		    for(w=0;w<act_band->loop->stmts[s]->nwrites;w++)
		    {
		    	/* [SICA] go through the array */
				int arraynameisnew=1; //store whether the prospected array name is already available or not
		    	a=0;
				while(strcmp(act_band->sicadata->id2arrayname[a],"\0"))    { //while there is a non-empty string in this array position
					//if there is one, check if it is equal to the one analysed at the moment
					if(!strcmp(act_band->sicadata->id2arrayname[a], act_band->loop->stmts[s]->writes[w]->name))    { //so they are identical
						arraynameisnew=0;
						break;
					}

					a++;
				}
				if(arraynameisnew) {
					//printf("[SICA] Copying the new name '%s' to the array to position %i\n", act_band->loop->stmts[s]->reads[r]->name, a);
					act_band->sicadata->id2arrayname[a]=strcpy(act_band->sicadata->id2arrayname[a], act_band->loop->stmts[s]->writes[w]->name);
					bands[i]->sicadata->nb_arrays++;
					//printf("[SICA] COPY SUCCEEDED\n");
				}
		    }
	    }

	    // [SICA] print the id2arrayname
	    printf("[SICA] Detected %i different arrays\n", bands[i]->sicadata->nb_arrays);
		for(a=0; a < bands[i]->sicadata->nb_arrays; a++)    {
				printf("[SICA] Array-ID: %i, name: %s\n", a, act_band->sicadata->id2arrayname[a]);
        }



	    printf("\tbands[%i], width=%i\n", i, act_band->width);
	    printf("\tbands[%i], nstmts=%i\n", i, act_band->loop->nstmts);

	    for(s=0; s<act_band->loop->nstmts;s++)
	    {
		    printf("\t\t stmt=%i: nreads=%i\n", s, act_band->loop->stmts[s]->nreads);
		    for(r=0;r<act_band->loop->stmts[s]->nreads;r++)
		    {
		    	printf("\t\t\tMATRIX:\n");
		    	for(x=0;x<act_band->loop->stmts[s]->reads[r]->mat->alloc_nrows;x++)
		    	{
		    		printf("\t\t\t");
			    	for(y=0;y<act_band->loop->stmts[s]->reads[r]->mat->alloc_ncols;y++)
			    	{
			    		printf("%i ", act_band->loop->stmts[s]->reads[r]->mat->val[x][y]);
			    	}
		            printf("\n");
		    	}
		    	printf("\t\t\t read=%i, name=%s\n", r, act_band->loop->stmts[s]->reads[r]->name);
		    	printf("\t\t\t read=%i, type=%s\n", r, act_band->loop->stmts[s]->reads[r]->symbol->data_type);
		    	printf("\n");
		    	pluto_matrix_print(stdout, act_band->loop->stmts[s]->reads[r]->mat);

		    	//calculatin the column offset
		        int firstD = act_band->loop->depth;
		        int lastD = act_band->loop->depth+act_band->width-1;
		        int widthD=lastD-firstD;
		    	printf("\t\t\t column offset: %i\n", widthD+1);
		    	printf("\n");
		    }

		    printf("\t\t stmt=%i: nwrites=%i\n", s, act_band->loop->stmts[s]->nwrites);
		    for(w=0;w<act_band->loop->stmts[s]->nwrites;w++)
		    {
		    	//printf("\t\t\t ARRAY-Id: %i\n", act_band->loop->stmts[s]->writes[w]->mat[0]);
		    	printf("\t\t\t write=%i, name=%s\n", w, act_band->loop->stmts[s]->writes[w]->name);
		    	printf("\t\t\t write=%i, type=%s\n", w, act_band->loop->stmts[s]->writes[w]->symbol->data_type);
		    	printf("\n");
		    	pluto_matrix_print(stdout, act_band->loop->stmts[s]->writes[w]->mat);
		    }


		    //get the PluTo defined data types + (null)->default

		    //recognize scalar dimensions by finding 0-mat and set them to be scalar ("0")

	    }



		///////////////////////////////////////////////

    	if(bands[i]->sicadata->isvec)
    	{
    	    printf("VEC\tbands[%i], nstmts=%i\n", i, act_band->loop->nstmts);
    		act_band->sicadata->sical1size=16; // [SICA] HERE A FUNCTION SHOULD BE CALLED THAT CALCULATES THE SICA SIZES FOR THAT BAND
    		act_band->sicadata->sical2size=4;  // [SICA] HERE A FUNCTION SHOULD BE CALLED THAT CALCULATES THE GLOBAL SIZE

    	    // CACHE Getting first informations about the accesses for cache-tiling
    	    get_cache_access_amount_function(scop, prog, act_band->sicadata->vecrow,&act_band->sicadata->sical1size,&act_band->sicadata->sical2size); //TODO

    	    printf("[CACHE] Level-1: %i, Level-2: %i\n",act_band->sicadata->sical1size, act_band->sicadata->sical2size );

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
