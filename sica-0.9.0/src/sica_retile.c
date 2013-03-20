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
	IF_DEBUG(if(band->sicadata->isvec){printf("VECTORIZED");}else{printf("VECTORIZED");});
    IF_DEBUG(      printf(" sica_retile_band step\n");  );
	
    int j, s;
    int depth, npar;

    npar = prog->npar;

    int firstD = band->loop->depth;
    int lastD = band->loop->depth+band->width-1;
    
    int widthD=lastD-firstD;

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

            printf("[SICA] before re-tiling\n");            
           pluto_constraints_print(stdout, stmt->domain);  

int dim;
for(dim=0;dim<sica_dimensions2tile;dim++)    {
printf("dim: %i, line: %i\n", dim, sica_dimensions2tile-dim);

printf("(%i,%i) orig: %i, new: %i\n",(stmt->domain->nrows-1)-(2*dim)+row_offset,dim+col_offset, stmt->domain->val[(stmt->domain->nrows-1)-(2*dim)+row_offset][dim+col_offset],tile_sizes[sica_dimensions2tile-dim-1]);
printf("(%i,%i) orig: %i, new: %i\n",(stmt->domain->nrows-1)-(2*dim)+row_offset,stmt->domain->ncols-1,stmt->domain->val[(stmt->domain->nrows-1)-(2*dim)+row_offset][stmt->domain->ncols-1], band->sicadata->upperboundoffset[s]+tile_sizes[sica_dimensions2tile-dim-1]-1);
printf("(%i,%i) orig: %i, new: %i\n",(stmt->domain->nrows-1)-(2*dim)-1+row_offset,dim+col_offset,stmt->domain->val[(stmt->domain->nrows-1)-(2*dim)-1+row_offset][dim+col_offset],-tile_sizes[sica_dimensions2tile-dim-1]);
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
            printf("[SICA] after re-tiling\n");            
           pluto_constraints_print(stdout, stmt->domain);  


	/* [SICA] OLD APPROACH */
//        for (depth=firstD; depth<=lastD; depth++)    {
//           Stmt *stmt = band->loop->stmts[s];
//            
//            int hyp_type = (stmt->hyp_types[depth + depth - firstD] == H_SCALAR)? H_SCALAR: 
//                H_TILE_SPACE_LOOP;
//
//            // pluto_constraints_print(stdout, stmt->domain);
//            IF_DEBUG2(printf("[SICA] before re-tiling\n");            );
////            IF_DEBUG2(pluto_constraints_print(stdout, stmt->domain);  );
//           if (hyp_type != H_SCALAR) {
//               assert(tile_sizes[depth-firstD] >= 1);
//               /* add offset on l1 tiling if l2 tiling will be performed */
//               int row_offset = -2*sica_dimensions2tile*offset;//-2*(widthD+1)*offset;
//               int col_offset = sica_dimensions2tile*offset;//(widthD+1)*offset;
//
//                //printf("DBG: row_offset=%i\n", row_offset);
//                //printf("DBG: col_offset=%i\n", col_offset);
//
////                printf("IS: width=%i\n", widthD);
////                printf("S%i act_band->sicadata->transwidth[s]=%i\n", s+1, band->sicadata->transwidth[s]);
//
////                printf("%i-%i+%i+%i,%i+%i-%i)\n", (stmt->domain->nrows-1),(2*(depth-firstD)),row_offset,2*skipped_scalar_dims,depth-firstD,col_offset,skipped_scalar_dims);
//
////                printf("%i, %i\n", (stmt->domain->nrows-1)-(2*(depth-firstD))+row_offset+2*skipped_scalar_dims,depth-firstD+col_offset-skipped_scalar_dims);
//                /* [SICA] ERROR CHECK if there is a value that should be resetted and was zero before */
//
//		/* [SICA] checkif this dimension was simply tiled to one before or if any error occurs */
//		if((stmt->domain->val[(stmt->domain->nrows-1)-(2*(depth-firstD))+row_offset+2*skipped_scalar_dims][depth-firstD+col_offset-skipped_scalar_dims] == 1)&&(stmt->domain->val[(stmt->domain->nrows-1)-(2*(depth-firstD))-1+row_offset+2*skipped_scalar_dims][depth-firstD+col_offset-skipped_scalar_dims] == -1)&&(stmt->domain->val[(stmt->domain->nrows-1)-(2*(depth-firstD))+row_offset+2*skipped_scalar_dims][stmt->domain->ncols-1] == 0))    {
//			printf("[SICA] WARNING: There is a tile dimension that was previously tiled with tile-size '1'\n");
//		} else {
//
//               	if((stmt->domain->val[(stmt->domain->nrows-1)-(2*(depth-firstD))+row_offset+2*skipped_scalar_dims][depth-firstD+col_offset-skipped_scalar_dims])==0||(stmt->domain->val[(stmt->domain->nrows-1)-(2*(depth-firstD))+row_offset+2*skipped_scalar_dims][stmt->domain->ncols-1])==0||(stmt->domain->val[(stmt->domain->nrows-1)-(2*(depth-firstD))-1+row_offset+2*skipped_scalar_dims][depth-firstD+col_offset-skipped_scalar_dims])==0)    {
//                		printf("[SICA] ERROR in S%i: there is something going wrong in the retiling step with the original values!\n", s+1);
//                		printf("(%i,%i)=%i",(stmt->domain->nrows-1)-(2*(depth-firstD))+row_offset+2*skipped_scalar_dims,stmt->domain->ncols-1,stmt->domain->val[(stmt->domain->nrows-1)-(2*(depth-firstD))+row_offset+2*skipped_scalar_dims][stmt->domain->ncols-1]);
//               		printf("(nrows-1=%i)-(act_rel_row=%i)+(row_offset=%i)+2*(scalar_skipped=%i)\n",stmt->domain->nrows-1,(2*(depth-firstD)),row_offset,skipped_scalar_dims);
//                		printf("pos=%i + col_offset=%i - skipped_scalar=%i\n",depth-firstD,col_offset,skipped_scalar_dims);
//                	}
//	
//                	/* [SICA] ERROR CHECK if there is a value that should be taken for reset and was zero before */
//                	if(tile_sizes[widthD-(depth-firstD)]<=0)    {
//                		printf("[SICA] ERROR in S%i: coordinates (%i,%i)\n", s+1, (stmt->domain->nrows-1)-(2*(depth-firstD))+row_offset+2*skipped_scalar_dims, depth-firstD+col_offset-skipped_scalar_dims);
//                		printf("[SICA] ERROR in S%i: there is something going wrong in the retiling step with the new tile values!\n", s+1);
//                		printf("[SICA] ERROR in S%i: tile_sizes[%i]=%i, band->sicadata->upperboundoffset[%i]=%i\n", s+1, widthD-(depth-firstD),tile_sizes[widthD-(depth-firstD)],s,band->sicadata->upperboundoffset[s]);
//                	}
//		}
//                                
//                /* upper bound of the tile-loop  X * tx + Y >= 0 (the relating outer tile-loop is stored in column depth+vecrow for l1 and depth+vecrow for l2) */
//                /* X */
//                IF_DEBUG2(printf("[SICA] Accessing line %i\n",(stmt->domain->nrows-1)-(2*(depth-firstD))+row_offset); );
//                stmt->domain->val[(stmt->domain->nrows-1)-(2*(depth-firstD))+row_offset+2*skipped_scalar_dims][depth-firstD+col_offset-skipped_scalar_dims] = tile_sizes[widthD-(depth-firstD)];
//                /* Y */
//                /* ToDo: is it correct that in stmt->trans->val there is NO OFFSET?!? LIKE THAT IT WORKS -> figure out the relation between trans element and domain value */
//                stmt->domain->val[(stmt->domain->nrows-1)-(2*(depth-firstD))+row_offset+2*skipped_scalar_dims][stmt->domain->ncols-1] = band->sicadata->upperboundoffset[s]+tile_sizes[widthD-(depth-firstD)]-1;
//                //IF_DEBUG(printf("[SICA] Statement %i -> Backuped upper bound offset: %i, current offset would be %i\n", s, band->sicadata->upperboundoffset[s], -stmt->trans->val[(depth-firstD)+1+depth][stmt->dim+prog->npar]););
//
//                /* lower bound of the tile-loop  X * tx >= 0 (the relating outer tile-loop is stored in column depth+vecrow for l1 and depth+vecrow for l2) */
//                /* X */
//                stmt->domain->val[(stmt->domain->nrows-1)-(2*(depth-firstD))-1+row_offset+2*skipped_scalar_dims][depth-firstD+col_offset-skipped_scalar_dims] = -tile_sizes[widthD-(depth-firstD)];
//               
//                //IF_DEBUG(printf("\n[SICA] MODIFIED: [%i,%i] with tile_sizes[%i]=%i\n", (stmt->domain->nrows-1)-(2*(depth-firstD)), depth-firstD, widthD-(depth-firstD), tile_sizes[widthD-(depth-firstD)]);          );
//                //IF_DEBUG(printf("[SICA] MODIFIED: [%i,%i] with tile_sizes[%i]=%i\n", (stmt->domain->nrows-1)-(2*(depth-firstD))-1, depth-firstD, widthD-(depth-firstD), tile_sizes[widthD-(depth-firstD)]);          );
//                //IF_DEBUG(printf("[SICA] MODIFIED: [%i,%i] with tile_sizes[%i]=%i\n", (stmt->domain->nrows-1)-(2*(depth-firstD))-1, stmt->domain->ncols-1, widthD-(depth-firstD), tile_sizes[widthD-(depth-firstD)]); );
//                
//                // printf("Stmt %d: depth: %d\n", stmt->id+1,depth);
//                // pluto_matrix_print(stdout, stmt->trans);
//
//            } else {
//            	IF_DEBUG2(printf("[SICA] SKIPPED LINE %i\n", depth + depth - firstD); );
//           	//skipped_scalar_dims++;
//            }
//            IF_DEBUG2(printf("[SICA] SKIPPED_SCALAR_DIMS=%i\n",skipped_scalar_dims); );
            IF_DEBUG2(printf("[SICA] after re-tiling\n");            );
            IF_DEBUG2(pluto_constraints_print(stdout, stmt->domain); );
        } /* all statements */
//    } /* all scats to be tiled */

    /* Re-detect hyperplane types (H_SCALAR, H_LOOP) */
    pluto_detect_hyperplane_types(prog);
    

	///////////////////////////////////////////////////////////////////////
}

/* Tiles scattering functions for all bands; l2=1 => perform l2 tiling */
void sica_retile_scattering_dims(PlutoProg *prog, Band **bands, int nbands, int l2)
{
//	printf("[SICA] Performing sica_retile_scattering_dims step\n");
	
    int i, j, b;

    Stmt **stmts = prog->stmts;

    IF_DEBUG(printf("[SICA] SICA modification for band number %i\n",b); );

    for (b=0; b<nbands; b++) {
    	
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

