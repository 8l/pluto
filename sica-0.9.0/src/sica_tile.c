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

#include "cache_math_func.h"

#include "sica_accesses.h"

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
    int nbands, i, a, s, r, w, x, y, t;
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
        bands[i]->sicadata->transwidth=0;
        bands[i]->sicadata->vec_accesses=0;
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

    /* [SICA] get the transformation matrices BEFORE prevector transformation */
    sica_get_trans_matrix(bands, nbands);

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
					act_band->sicadata->nb_arrays++;
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
					act_band->sicadata->nb_arrays++;
					//printf("[SICA] COPY SUCCEEDED\n");
				}
		    }
	    }

	    // [SICA] print the id2arrayname
	    printf("[SICA] Detected %i different arrays\n", act_band->sicadata->nb_arrays);
		for(a=0; a < act_band->sicadata->nb_arrays; a++)    {
				printf("[SICA] Array-ID: %i, name: %s\n", a, act_band->sicadata->id2arrayname[a]);
        }

		// [SICA] setup the empty array of sica_accesses structures
		SICAAccess** sica_accesses_on_array;
		sica_accesses_on_array=(SICAAccess**)malloc(act_band->sicadata->nb_arrays*sizeof(SICAAccess*));
		for(a=0; a<act_band->sicadata->nb_arrays; a++)    {
			sica_accesses_on_array[a]=sica_accesses_malloc();
		}


	    printf("\tbands[%i], width=%i\n", i, act_band->width);
	    printf("\tbands[%i], nstmts=%i\n", i, act_band->loop->nstmts);

	    for(s=0; s<act_band->loop->nstmts;s++)
	    {
	    	// [SICA] print the transformation matrix
	    	printf("T(S%i):\n",act_band->loop->stmts[s]->id);
			pluto_matrix_print(stdout, act_band->loop->stmts[s]->trans);

	    	// [SICA] Print the extracted transformation matrix
	    	printf("[SICA] Transformation matrix:\n");
	    	for(y=0; y < act_band->sicadata->transwidth; y++)    {
			    for(x=0; x < act_band->sicadata->transwidth; x++)    {
			    	printf("%i ",act_band->sicadata->trans[y][x]);
		 	   }
			    printf("\n");
	    	}

	    	// [SICA] Print the inverted transformation matrix
	    	printf("[SICA] Inverted Transformation matrix:\n");
	    	for(y=0; y < act_band->sicadata->transwidth; y++)    {
		   	    for(x=0; x < act_band->sicadata->transwidth; x++)    {
		   	    	printf("%i ",act_band->sicadata->trans_inverted[y][x]);
		   	    }
		   	    printf("\n");
	    	}

	    	//->START READ ANALYSIS
		    printf("[SICA] Analyse the READ accesses for vectorization\n");
	    	printf("stmt=%i: nreads=%i\n", s, act_band->loop->stmts[s]->nreads);
		    for(r=0;r<act_band->loop->stmts[s]->nreads;r++)
		    {
		    	printf("\t\t\t read=%i, id: %i, name=%s, type=%s\n", r, sica_get_array_id(act_band, act_band->loop->stmts[s]->reads[r]->name), act_band->loop->stmts[s]->reads[r]->name, act_band->loop->stmts[s]->reads[r]->symbol->data_type);
		    	pluto_matrix_print(stdout, act_band->loop->stmts[s]->reads[r]->mat);

		    	//get the access matrix offset caused by tiling dimensions
		    	int access_offset=act_band->sicadata->transwidth;
		    	if(options->l2tile)    {
		    		access_offset=2*act_band->sicadata->transwidth;
		    	}

		    	//get the access matrix
		    	printf("\t\t\t READ-ACCESS-MATRIX:\n");
		    	for(x=0;x<act_band->loop->stmts[s]->reads[r]->mat->alloc_nrows;x++)    {
		    		printf("\t\t\t ");
			    	for(y=0;y<act_band->sicadata->transwidth;y++)    {
			    		printf("%i ", act_band->loop->stmts[s]->reads[r]->mat->val[x][access_offset+y]);
			    	}
		            printf("\n");
		    	}
		    	printf("\n");

		    	//go through all dimensions of this access and transform them
		    	int** orig_access_mat;
		    	int** trans_access_mat;

		    	//malloc the transformed mat to fill it
	    		orig_access_mat=(int**)malloc(act_band->loop->stmts[s]->reads[r]->mat->alloc_nrows*sizeof(int*));
	    		trans_access_mat=(int**)malloc(act_band->loop->stmts[s]->reads[r]->mat->alloc_nrows*sizeof(int*));

		    	for(y=0; y<act_band->loop->stmts[s]->reads[r]->mat->alloc_nrows; y++)    {
		    		orig_access_mat[y]=(int*)malloc(act_band->loop->stmts[s]->reads[r]->mat->alloc_ncols*sizeof(int));
		    		trans_access_mat[y]=(int*)malloc(act_band->loop->stmts[s]->reads[r]->mat->alloc_ncols*sizeof(int));
		    	}

		    	int orig_access_iterators[act_band->sicadata->transwidth];
		    	int trans_access_iterators[act_band->sicadata->transwidth];

		    	//fill the two access and row arrays
		    	for(y=0;y<act_band->loop->stmts[s]->reads[r]->mat->alloc_nrows;y++)    {

			    	for(x=0;x<act_band->loop->stmts[s]->reads[r]->mat->alloc_ncols;x++)    {
			    		orig_access_mat[y][x] = act_band->loop->stmts[s]->reads[r]->mat->val[y][x];
			    		trans_access_mat[y][x] = act_band->loop->stmts[s]->reads[r]->mat->val[y][x]; //after transformation overwrite the trans part
			    	}

			    	for(x=0;x<act_band->sicadata->transwidth;x++)    {
			    		orig_access_iterators[x] = act_band->loop->stmts[s]->reads[r]->mat->val[y][access_offset+x];
			    		trans_access_iterators[x] = 0;
				    	}

		    		sica_vec_times_matrix(trans_access_iterators,orig_access_iterators,act_band->sicadata->trans_inverted,act_band->sicadata->transwidth,act_band->sicadata->transwidth);

			    	//print the two access arrays
		    		//printf("ORIG-ACCESS:\t");
		    		//for(y=0;y<act_band->sicadata->transwidth;y++)    {
				    //		printf("%i ", orig_access[y]);
				    //	}
		    		//printf("\n");
                    //
		    		//printf("TRANS-ACCESS:\t");
		    		//for(y=0;y<act_band->sicadata->transwidth;y++)    {
				    //		printf("%i ", trans_access[y]);
				    //	}
		    		///printf("\n");

			    	//overwrite the transformed parts in the trans_row
			    	for(x=0;x<act_band->sicadata->transwidth;x++)    {
			    		trans_access_mat[y][access_offset+x] = trans_access_iterators[x];
			    	}
		    	}

		    	int act_array_id=sica_get_array_id(act_band, act_band->loop->stmts[s]->reads[r]->name);
		    	printf("\t\t\t[SICA] Looking up the following array name: '%s' with id: '%i' and transformed access matrix:\n", act_band->loop->stmts[s]->reads[r]->name, act_array_id);
		    	//This is now the transformed access
		    	for(y=0;y<act_band->loop->stmts[s]->reads[r]->mat->alloc_nrows;y++)    {
		    		printf("\t\t\t");
		    		for(x=0;x<act_band->loop->stmts[s]->reads[r]->mat->alloc_ncols;x++)    {
							printf("%i ", trans_access_mat[y][x]);
						}
		    		printf("\n");
		    	}

		    	//Now we have all necessary data concerning this access, check if we have to count it for vectorized accesses!

		    	//check if the access is a just a scalar type, not a real array access (e.g. alpha or a[const][const])
		    	int entry_sum=sica_get_entry_sum(trans_access_mat, act_band->loop->stmts[s]->reads[r]->mat->alloc_nrows, act_band->loop->stmts[s]->reads[r]->mat->alloc_ncols);

		    	//IF it is a real array access...
		    	if(entry_sum)    {
		    		printf("\t\t\tThe Access on Array '%s' is an array access!\n", act_band->loop->stmts[s]->reads[r]->name);
			    	for(y=0;y<act_band->loop->stmts[s]->reads[r]->mat->alloc_nrows;y++)    {

			    		//...CHECK IF there is one or more dimensions that is accessed by the vectorized loop and...
						if(act_band->sicadata->isvec&&trans_access_mat[y][access_offset+act_band->sicadata->vecrow])    {
							printf("\t\t\t[SICA] VECTORIZATION: Array '%s' accesses dimension '%i' by a vectorized loop!\n", act_band->loop->stmts[s]->reads[r]->name, y);

							int sica_access_is_new=1;

							//...CHECK IF the access (matrix) is already recognized for this array (e.g. C[i][j]=C[i][j]+...). ...
							printf("[SICA] ACT-POINTER: %p\n",sica_accesses_on_array[act_array_id]);
							SICAAccess* act_access_temp=sica_accesses_on_array[act_array_id];

							while(act_access_temp->next)    {
								//CHECK IF THE act_access_temp->mat is equal to the actual one
								int check_comparison=0;
								check_comparison=sica_compare_access_matrices(trans_access_mat, act_access_temp->access_mat, act_band->loop->stmts[s]->reads[r]->mat->alloc_nrows, act_band->loop->stmts[s]->reads[r]->mat->alloc_ncols);

								//IF THEY ARE EQUAL
								if(check_comparison)    {
									printf("THIS ACCESS IS ALREADY RECOGNIZED AND THEREFORE NOT ADDED!\”");
									sica_access_is_new=0;
									break;
								}
								act_access_temp=act_access_temp->next;
							}

							//...IF NOT, add it to the linked list for this array and increase the counter ...
							if(sica_access_is_new)    {
								act_access_temp->nrows = act_band->loop->stmts[s]->reads[r]->mat->alloc_nrows;
								act_access_temp->ncols = act_band->loop->stmts[s]->reads[r]->mat->alloc_ncols;
								act_band->sicadata->vec_accesses++;
								printf("[SICA] adding the new access on pointer %p!\n",act_access_temp);
								act_access_temp->next=sica_accesses_malloc();
								printf("[SICA] after malloc             pointer %p!\n",act_access_temp);

								act_access_temp->access_mat=sica_access_matrix_malloc(act_band->loop->stmts[s]->reads[r]->mat->alloc_nrows, act_band->loop->stmts[s]->reads[r]->mat->alloc_ncols);
								sica_copy_access_matrix(act_access_temp->access_mat, trans_access_mat, act_band->loop->stmts[s]->reads[r]->mat->alloc_nrows, act_band->loop->stmts[s]->reads[r]->mat->alloc_ncols);
							}

						}
						if(!act_band->sicadata->isvec)    {
							printf("\t\t\t[SICA] THIS IS NOT A VECTORIZED BAND\n");
						}

			    	}
		    	} else {
		    		printf("\t\t\tThe Access on Array '%s' is NO array access!\n", act_band->loop->stmts[s]->reads[r]->name);
		    	}

		    }
		    //->STOP READ ANALYSIS


	    	//->START WRITE ANALYSIS
		    printf("\n[SICA] Analyse the WRITE accesses for vectorization\n");
		    printf("stmt=%i: nwrites=%i\n", s, act_band->loop->stmts[s]->nwrites);
		    for(w=0;w<act_band->loop->stmts[s]->nwrites;w++)
		    {
		    	printf("\t\t\tWRITE-MATRIX:\n");
		    	for(x=0;x<act_band->loop->stmts[s]->writes[w]->mat->alloc_nrows;x++)
		    	{
		    		printf("\t\t\t");
			    	for(y=0;y<act_band->loop->stmts[s]->writes[w]->mat->alloc_ncols;y++)
			    	{
			    		printf("%i ", act_band->loop->stmts[s]->writes[w]->mat->val[x][y]);
			    	}
		            printf("\n");
		    	}
		    	printf("\t\t\t write=%i, name=%s\n", w, act_band->loop->stmts[s]->writes[w]->name);
		    	printf("\t\t\t write=%i, type=%s\n", w, act_band->loop->stmts[s]->writes[w]->symbol->data_type);
		    	printf("\n");
		    	pluto_matrix_print(stdout, act_band->loop->stmts[s]->writes[w]->mat);
		    }
	    	//->STOP WRITE ANALYSIS


	    }

	    printf("[SICA] REPORT: There are %i accesses relevant for vectorization\n", act_band->sicadata->vec_accesses);

	    printf("[SICA] Print the accesses on arrays structure:\n");
		for(a=0; a<act_band->sicadata->nb_arrays; a++)    {
			printf("\tArray '%s': \n",act_band->sicadata->id2arrayname[a]);
			SICAAccess* act_access_temp=sica_accesses_on_array[a];
			while(act_access_temp->next)    {
				printf("\t\t%p:\n", act_access_temp);
				act_access_temp=act_access_temp->next;
			}
		}
	    //get the PluTo defined data types + (null)->default and calculate the tile quantities

	    //NOW LOKK UP HOW MANY ACCESSES OF WHICH TYPE ARE AVAILABLE



		///////////////////////////////////////////////

    	if(act_band->sicadata->isvec)
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
