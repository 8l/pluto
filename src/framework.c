/*
 * PLUTO: An automatic parallelizer and locality optimizer
 * 
 * Copyright (C) 2007--2008 Uday Kumar Bondhugula
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * A copy of the GNU General Public Licence can be found in the file
 * `LICENSE' in the top-level directory of this distribution. 
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "math_support.h"
#include "constraints.h"
#include "pluto.h"

#include <isl/constraint.h>
#include <isl/mat.h>
#include <isl/set.h>
#include "candl/candl.h"

static void eliminate_farkas_multipliers(PlutoConstraints *farkas_cst, int num_elim);

/**
 *
 * Each constraint row is represented as follows
 *
 *      [comm. vol bound | mapping coeff.s for S1, S2,... |constant]
 * Size:[    npar+1      | (nvar+1)*nstmts                | 1      ]
 *
 * npar - number of parameters in whole program
 * nvar - number of parameters in whole program
 *
 */
#if 0
static PlutoConstraints *get_permutability_constraints_uniform_dep (Dep *dep)
{
    int cst_offset;
    int j, dest_stmt;
    PlutoConstraints *cst;

    /* constant dependences */
    /* uniform self-edge, no need to apply farkas */
    dest_stmt = dep->dest;


    cst_offset = npar+1+dest_stmt*(nvar+1);

    cst = constraints_alloc(2, CST_WIDTH);
    cst->ncols = CST_WIDTH;

    if (!IS_RAR(dep->type)) {
        cst->nrows = 2;
        /* Tiling legality constraint */
        for (j=0; j<nvar; j++)  {
            cst->val[0][cst_offset+j] = -dep->h->val[j][nvar+npar];
        }
        /* Translation coefficient */
        cst->val[0][cst_offset+nvar]=0;

        /* Add bounding function */
        for (j=0; j<npar; j++)  {
            cst->val[1][j] = 0;
        }
        cst->val[1][npar] = 1;
        for (j=cst_offset; j<cst_offset+nvar; j++)  {
            cst->val[1][j] = -cst->val[0][j];
        }
        cst->val[1][cst_offset+nvar]=0;
    }else{
        /* Add bounding function */
        for (j=0; j<npar; j++)  {
            cst->val[0][j] = 0;
        }
        cst->val[0][npar] = 1;
        for (j=cst_offset; j<cst_offset+nvar; j++)  {
            cst->val[0][j] = dep->h->val[j-cst_offset][nvar+npar];
        }
        cst->val[0][cst_offset+nvar]=0;
        cst->nrows=1;
    }

    return cst;
}
#endif


/* Builds legality constraints for a non-uniform dependence */
static PlutoConstraints *get_permutability_constraints_nonuniform_dep(Dep *dep, const PlutoProg *prog)
{
    PlutoConstraints *farkas_cst, *comm_farkas_cst, *cst;
    int src_stmt, dest_stmt, j, k;
    int src_offset, dest_offset;

    int nvar = prog->nvar;
    int npar = prog->npar;
    Stmt *stmts = prog->stmts;
    int nstmts = prog->nstmts;

    dest_stmt = dep->dest;
    src_stmt = dep->src;

    /* Non-uniform dependence - farkas lemma comes in */
    /* Apply farkas lemma, eliminate farkas multipliers using
     * fourier-motzkin 
     * 
     * -- farkas_cst format for legality --
     * [ mapping coeff for src | ... for dest |farkas multipliers|constant]
     * SIZE: [nvar+1 | nvar+1 | dep.dpolytope->nrows+1 | 1]
     *
     * -- farkas_cst format for bounding function --
     * [bounding func | mapping coeff for src | ... for dest |farkas multipliers|constant]
     * SIZE: [npar+1| nvar+1 | nvar+1 | dep.dpolytope->nrows+1 | 1]
     *
     */
    if (src_stmt != dest_stmt)  {
        /* Inter-statement non-uniform dep */
        farkas_cst = pluto_constraints_alloc(MAX_FARKAS_CST, 2*nvar+2+dep->dpolytope->nrows+2);
        comm_farkas_cst = pluto_constraints_alloc(MAX_FARKAS_CST, npar+1+2*nvar+2+dep->dpolytope->nrows+2);

        farkas_cst->nrows = (2*nvar+npar+1)+1+dep->dpolytope->nrows+1;
        farkas_cst->ncols = 2*(nvar+1)+dep->dpolytope->nrows+2;

        comm_farkas_cst->nrows = (2*nvar+npar+1)+1+dep->dpolytope->nrows+1;
        comm_farkas_cst->ncols = npar+1+2*(nvar+1)+dep->dpolytope->nrows+2;
    }else{
        /* Intra-statement non-uniform dependence */
        farkas_cst = pluto_constraints_alloc(MAX_FARKAS_CST, nvar+1+dep->dpolytope->nrows+2);
        comm_farkas_cst = pluto_constraints_alloc(MAX_FARKAS_CST, npar+1+nvar+1+dep->dpolytope->nrows+2);

        farkas_cst->nrows = (2*nvar+npar+1)+1+dep->dpolytope->nrows+1;
        farkas_cst->ncols = (nvar+1)+dep->dpolytope->nrows+2;

        comm_farkas_cst->nrows = (2*nvar+npar+1)+1+dep->dpolytope->nrows+1;
        comm_farkas_cst->ncols = npar+1+(nvar+1)+dep->dpolytope->nrows+2;
    }


    /* Initialize all to zero */
    for (j=0; j<farkas_cst->nrows; j++)  {
        for (k=0; k<farkas_cst->ncols; k++)  {
            farkas_cst->val[j][k] = 0;
        }
    }

    for (j=0; j<comm_farkas_cst->nrows; j++)  {
        for (k=0; k<comm_farkas_cst->ncols; k++)  {
            comm_farkas_cst->val[j][k] = 0;
        }
    }

    if (src_stmt != dest_stmt)  {

        /* Add tiling legality constraints */
        for (j=0; j<2*nvar+npar+1; j++)  {
            if (j < nvar)   {
                /* src stmt coeff */
                farkas_cst->val[j][j] = -1;
            }else if (j < 2*nvar)   {
                /* dest stmt coeff */
                farkas_cst->val[j][j+1] = 1;
            }else if (j < 2*nvar+npar)  {
                /* Do nothing - all coeff multipliers stay zero */
                /* since structure parameters not in our affine mappings */
            }else{
                /* j = 2*nvar+npar */
                /* Translation coefficients in the affine mappings */
                farkas_cst->val[j][nvar] = -1;
                farkas_cst->val[j][2*nvar+1] = 1;
                /* \lambda_0 */
                farkas_cst->val[j][farkas_cst->ncols-2] = -1;
            } 

            /* Set coeff's for farkas multipliers (all except \lambda_0) */
            for (k=2*nvar+2; k<2*nvar+2+dep->dpolytope->nrows; k++)  {
                /* Note that dep polytope is dpolytope->nrows x (2*nvar+npar+1) */
                farkas_cst->val[j][k] = -dep->dpolytope->val[k-2*nvar-2][j];
            }
            farkas_cst->val[j][farkas_cst->ncols-1] = 0;
        }

        /* Since the above are equalities - add sigma negative */
        for (k=0; k<farkas_cst->ncols; k++)    {
            farkas_cst->val[2*nvar+npar+1][k] = 0;
            for (j=0; j<2*nvar+npar+1; j++)  {
                farkas_cst->val[2*nvar+npar+1][k] -= farkas_cst->val[j][k];
            }
        }

        /* All Farkas multipliers are non-negative */
        for (j=0; j<dep->dpolytope->nrows+1; j++)  {
            for (k=0; k<dep->dpolytope->nrows+1; k++)  {
                farkas_cst->val[2*nvar+npar+2+j][2*nvar+2+k] = ((j==k)?1:0);
            }
        }

        /* Bounding function constraints */
        for (k=0; k<npar; k++)  {
            comm_farkas_cst->val[2*nvar+k][k] = 1;
        }

        comm_farkas_cst->val[2*nvar+npar][npar] = 1;

        for (j=0; j<2*nvar+npar+1; j++)
            for (k=0; k<farkas_cst->ncols-dep->dpolytope->nrows-2; k++)
                comm_farkas_cst->val[j][npar+1+k] = -farkas_cst->val[j][k];

        for (j=0; j<2*nvar+npar+1; j++)
            for (k=farkas_cst->ncols-dep->dpolytope->nrows-2; k<farkas_cst->ncols; k++)
                comm_farkas_cst->val[j][npar+1+k] = farkas_cst->val[j][k];

        /* Add opp inequality since the above were equalities */
        for (k=0; k<comm_farkas_cst->ncols; k++)    {
            comm_farkas_cst->val[2*nvar+npar+1][k] = 0;
            for (j=0; j<2*nvar+npar+1; j++) {
                comm_farkas_cst->val[2*nvar+npar+1][k] -= comm_farkas_cst->val[j][k];
            }
        }

        for (j=2*nvar+npar+2; j<farkas_cst->nrows; j++)
            for (k=0; k<farkas_cst->ncols; k++)
                comm_farkas_cst->val[j][npar+1+k] = farkas_cst->val[j][k];
        
        eliminate_farkas_multipliers(farkas_cst, farkas_cst->ncols-2*nvar-3);
        eliminate_farkas_multipliers(comm_farkas_cst, comm_farkas_cst->ncols-npar-1-2*nvar-3);

        /* constraints_print(stdout, farkas_cst); */

    }else{
        /* Source stmt == Dest stmt */

        for (j=0; j<2*nvar+npar+1; j++)  {
            if (j < nvar)   {
                /* src stmt coeff */
                farkas_cst->val[j][j] = -1;
            }else if (j < 2*nvar)   {
                /* dest stmt coeff */
                farkas_cst->val[j][j-nvar] = 1;
            }else if (j < 2*nvar+npar)  {
                /* Do nothing - all coeff multipliers stay zero */
                /* NOTE: structure parameters not in our affine mappings */
            }else{
                /* Translation coefficient gets subtracted out */
                farkas_cst->val[j][nvar] = 0;
                farkas_cst->val[j][farkas_cst->ncols-2] = -1;
            } 

            /* Set coeff's for farkas multipliers */
            for (k=nvar+1; k<nvar+1+dep->dpolytope->nrows; k++)  {
                farkas_cst->val[j][k] = -dep->dpolytope->val[k-nvar-1][j];
            }
            farkas_cst->val[j][farkas_cst->ncols-1] = 0;
        }

        /* Since the above are equalities - add sigma negative */
        for (k=0; k<farkas_cst->ncols; k++)    {
            farkas_cst->val[2*nvar+npar+1][k] = 0;
            for (j=0; j<2*nvar+npar+1; j++)  {
                farkas_cst->val[2*nvar+npar+1][k] -= farkas_cst->val[j][k];
            }
        }

        /* All farkas multipliers are positive */
        for (j=0; j<dep->dpolytope->nrows+1; j++)  {
            for (k=0; k<dep->dpolytope->nrows+1; k++)  {
                farkas_cst->val[2*nvar+npar+2+j][nvar+1+k] = ((j==k)?1:0);
            }
        }

        /* Bounding function constraints */
        for (k=0; k<npar; k++)  {
            comm_farkas_cst->val[2*nvar+k][k] = 1;
        }

        comm_farkas_cst->val[2*nvar+npar][npar] = 1;

        for (j=0; j<2*nvar+npar+1; j++)
            for (k=0; k<farkas_cst->ncols-dep->dpolytope->nrows-2; k++)
                comm_farkas_cst->val[j][npar+1+k] = -farkas_cst->val[j][k];

        for (j=0; j<2*nvar+npar+1; j++)
            for (k=farkas_cst->ncols-dep->dpolytope->nrows-2; k<farkas_cst->ncols; k++)
                comm_farkas_cst->val[j][npar+1+k] = farkas_cst->val[j][k];

        /* Add opp inequality since the above were equalities */
        for (k=0; k<comm_farkas_cst->ncols; k++)    {
            comm_farkas_cst->val[2*nvar+npar+1][k] = 0;
            for (j=0; j<2*nvar+npar+1; j++) {
                comm_farkas_cst->val[2*nvar+npar+1][k] -= comm_farkas_cst->val[j][k];
            }
        }

        for (j=2*nvar+npar+2; j<farkas_cst->nrows; j++)
            for (k=0; k<farkas_cst->ncols; k++)
                comm_farkas_cst->val[j][npar+1+k] = farkas_cst->val[j][k];

        eliminate_farkas_multipliers(farkas_cst, farkas_cst->ncols-nvar-2);
        eliminate_farkas_multipliers(comm_farkas_cst, comm_farkas_cst->ncols-npar-1-nvar-2);

        /* constraints_print(stdout, farkas_cst); */
    }

    /* Aggregate permutability and bounding function constraints together in
     * global format format */

    /* Initialize everything to zero */
    cst = pluto_constraints_alloc(farkas_cst->nrows + comm_farkas_cst->nrows, CST_WIDTH);
    cst->ncols = CST_WIDTH;

    for (k=0; k<farkas_cst->nrows+comm_farkas_cst->nrows; k++)   {
        for (j=0; j<cst->ncols; j++)  {
            cst->val[k][j] = 0;
        }
    }

    src_offset = npar+1+src_stmt*(nvar+1);
    dest_offset = npar+1+dest_stmt*(nvar+1);

    /* Permutability constraints */
    if (!IS_RAR(dep->type)) {
        for (k=0; k<farkas_cst->nrows; k++)   {
            for (j=0; j<nvar+1; j++)  {
                cst->val[cst->nrows+k][src_offset+j] = farkas_cst->val[k][j];
                if (src_stmt != dest_stmt) {
                    cst->val[cst->nrows+k][dest_offset+j] = farkas_cst->val[k][nvar+1+j];
                }
            }
            /* constant part */
            if (src_stmt == dest_stmt)  {
                cst->val[cst->nrows+k][cst->ncols-1] = farkas_cst->val[k][nvar+1];
            }else{
                cst->val[cst->nrows+k][cst->ncols-1] = farkas_cst->val[k][2*nvar+2];
            }
        }
        cst->nrows = farkas_cst->nrows;
    }

    if (!options->nobound)   {
        /* Add bounding constraints */
        src_offset = npar+1+src_stmt*(nvar+1);
        dest_offset = npar+1+dest_stmt*(nvar+1);

        for (k=0; k<comm_farkas_cst->nrows; k++)   {
            for (j=0; j<npar+1; j++)  {
                cst->val[cst->nrows+k][j] = comm_farkas_cst->val[k][j];
            }
            for (j=0; j<nvar+1; j++)  {
                cst->val[cst->nrows+k][src_offset+j] = comm_farkas_cst->val[k][npar+1+j];
                if (src_stmt != dest_stmt) cst->val[cst->nrows+k][dest_offset+j] = comm_farkas_cst->val[k][npar+1+nvar+1+j];
            }
            /* constant part */
            if (src_stmt == dest_stmt)  {
                cst->val[cst->nrows+k][cst->ncols-1] = comm_farkas_cst->val[k][npar+1+nvar+1];
            }else{
                cst->val[cst->nrows+k][cst->ncols-1] = comm_farkas_cst->val[k][npar+1+2*nvar+2];
            }
        }
        cst->nrows += comm_farkas_cst->nrows;
    }


    /* Coefficients of those variables that don't appear in the outer loop
     * are useless */
    for (k=0; k<nvar; k++)    {
        if (!stmts[src_stmt].is_outer_loop[k])  {
            for (j=0; j < cst->nrows; j++)   {
                cst->val[j][src_offset+k] = 0;
            }
        }
        if (src_stmt != dest_offset && !stmts[dest_stmt].is_outer_loop[k])  {
            for (j=0; j < farkas_cst->nrows+comm_farkas_cst->nrows; j++)   {
                cst->val[j][dest_offset+k] = 0;
            }
        }
    }

    pluto_constraints_free(farkas_cst);
    pluto_constraints_free(comm_farkas_cst);

    return cst;
}


PlutoConstraints *get_permutability_constraints(Dep *deps, int ndeps, 
        const PlutoProg *prog)
{
    int i, dest_stmt, src_stmt;
    Dep *dep;
    static PlutoConstraints *globcst = NULL;
    static PlutoConstraints **depcst = NULL;

    int nstmts = prog->nstmts;
    int nvar = prog->nvar;
    int npar = prog->npar;

    if (!depcst)   {
        depcst = (PlutoConstraints **) malloc(ndeps*sizeof(PlutoConstraints *));
        for (i=0; i<ndeps; i++) {
            depcst[i] = NULL;
        }
    }

    int total_cst_rows = 0;

// #pragma omp parallel for private(i,dep,dest_stmt,src_stmt) reduction(+:total_cst)
    for (i=0; i<ndeps; i++) {
        dep = &deps[i];

        dest_stmt = dep->dest;
        src_stmt = dep->src;

        if (options->rar == 0 && IS_RAR(dep->type))  {
            continue;
        }

        if (!depcst[i]) {
            /* First time, get the constraints */

            // Candl doesn't separate out uniform depedences and
            // h-transformation
            // if (src_stmt == dest_stmt && IS_UNIFORM(deps[i].type)) {
                /* Uniform self-edge */
                // depcst[i] = get_permutability_constraints_uniform_dep(dep);
            // }else{
                /* Non-uniform dependences */
            depcst[i] = get_permutability_constraints_nonuniform_dep(dep, prog);
            // }

            IF_DEBUG(fprintf(stdout, "After dep: %d; num_constraints: %d\n", i+1, depcst[i]->nrows));
            total_cst_rows += depcst[i]->nrows;
        }
    }

    if (!globcst) globcst = pluto_constraints_alloc(total_cst_rows, CST_WIDTH);
    globcst->ncols = CST_WIDTH;
    globcst->nrows = 0;

    for (i=0; i<ndeps; i++) {
        dep = &deps[i];

        if (options->rar == 0 && IS_RAR(dep->type))  {
            continue;
        }

        /* Note that dependences would be marked satisfied (in
         * pluto_auto_transform) only after all possible independent solutions 
         * are found to the formulation
         */ 
        if (dep_is_satisfied(dep)) continue;

        /* Subsequent calls can just use the old ones */
        pluto_constraints_add(globcst, depcst[i]);

        IF_DEBUG(fprintf(stdout, "After dep: %d; num_constraints: %d\n", i+1, globcst->nrows));
        if (globcst->nrows >= 0.7*MAX_CONSTRAINTS)  {
            IF_DEBUG(fprintf(stdout, "After dep: %d; num_constraints_simplified: %d\n", i+1, globcst->nrows));
        }
        pluto_constraints_simplify(globcst);
        IF_DEBUG2(pluto_constraints_print(stdout, globcst));
    }

    pluto_constraints_simplify(globcst);

    IF_DEBUG(fprintf(stdout, "After all dependences: num constraints: %d\n", globcst->nrows));

    return globcst;
}


/* PlutoConstraints to avoid trivial solutions (all zeros) */
PlutoConstraints *get_non_trivial_sol_constraints(const PlutoProg *prog)
{
    PlutoConstraints *nzcst;
    int i, j, stmt_offset;

    Stmt *stmts = prog->stmts;
    int nstmts = prog->nstmts;
    int nvar = prog->nvar;
    int npar = prog->npar;

    nzcst = pluto_constraints_alloc(nstmts, CST_WIDTH);
    nzcst->ncols = CST_WIDTH;

    for (i=0; i<nstmts; i++) {
        /* Don't add the constraint if enough solutions have been found */
        if (stmts[i].num_ind_sols >= stmts[i].dim)   {
            IF_DEBUG2(fprintf(stdout, "non-zero cst: skipping stmt %d\n", i));
            continue;
        }
        stmt_offset = npar+1+i*(nvar+1);
        for (j=0; j<nvar; j++)  {
            if (stmts[i].is_outer_loop[j] == 1)
                nzcst->val[nzcst->nrows][stmt_offset+j] = 1;
        }
        nzcst->val[nzcst->nrows][CST_WIDTH-1] = -1;
        nzcst->nrows++;
    }

    return nzcst;
}


/*
 * Eliminates the last num_elim variables from farkas_cst -- these are the
 * farkas multipliers
 */
static void eliminate_farkas_multipliers(PlutoConstraints *farkas_cst, int num_elim)
{
    int i;
    int best_elim;

    /* printf("To start with: %d constraints, %d to be eliminated out of %d\n", 
            farkas_cst->nrows, num_elim, farkas_cst->ncols-1); */

    for (i=0; i<num_elim; i++)  {
        best_elim = best_elim_candidate(farkas_cst, num_elim-i);
        fourier_motzkin_eliminate(farkas_cst, best_elim);
        /* printf("After elimination of %d variable: %d constraints\n", 
                num_elim-i, farkas_cst->nrows); */
        /* constraints_print(stdout, farkas_cst); */
    }

}


/*
 * Construct a PlutoMatrix with the same content as the given isl_mat.
 */
static PlutoMatrix *pluto_matrix_from_isl_mat(__isl_keep isl_mat *mat)
{
    int i, j;
    int rows, cols;
    isl_int v;
    PlutoMatrix *pluto;

    rows = isl_mat_rows(mat);
    cols = isl_mat_cols(mat);
    pluto = pluto_matrix_alloc(rows, cols);

    isl_int_init(v);

    for (i = 0; i < rows; ++i)
       for (j = 0; j < cols; ++j) {
           isl_mat_get_element(mat, i, j, &v);
           pluto->val[i][j] = isl_int_get_si(v);
       }

    isl_int_clear(v);

    return pluto;
}


/*
 * Construct a non-parametric basic set from the constraints in cst.
 */
static __isl_give isl_basic_set *isl_basic_set_from_pluto_constraints(
       isl_ctx *ctx, const PlutoConstraints *cst)
{
    int i, j;
    isl_int v;
    isl_dim *dim;
    isl_constraint *c;
    isl_basic_set *bset;

    isl_int_init(v);

    dim = isl_dim_set_alloc(ctx, 0, cst->ncols - 1);
    bset = isl_basic_set_universe(isl_dim_copy(dim));

    for (i = 0; i < cst->nrows; ++i) {
       if (cst->is_eq[i])
           c = isl_equality_alloc(isl_dim_copy(dim));
       else
           c = isl_inequality_alloc(isl_dim_copy(dim));

       isl_int_set_si(v, cst->val[i][cst->ncols - 1]);
       isl_constraint_set_constant(c, v);

       for (j = 0; j < cst->ncols - 1; ++j) {
           isl_int_set_si(v, cst->val[i][j]);
           isl_constraint_set_coefficient(c, isl_dim_set, j, v);
       }

       bset = isl_basic_set_add_constraint(bset, c);
    }

    isl_dim_free(dim);

    isl_int_clear(v);

    return bset;
}


/*
 * Negate the single constraint in cst.
 */
static void negate_constraint(PlutoConstraints *cst)
{
    int i;

    for (i = 0; i < cst->ncols; ++i)
       cst->val[0][i] = -cst->val[0][i];
}


/*
 * Returns linear independence constraints for a single statement.
 *
 * In particular, if H contains the first rows of an affine transformation,
 * then return a constraint on the coefficients of the next row that
 * ensures that this next row is linearly independent of the first rows.
 * Furthermore, the constraint is constructed in such a way that it allows
 * for a solution when combined with the other constraints on the coefficients
 * (currcst), provided any such constraint can be constructed.
 *
 * We do this by computing a basis for the null space of H and returning
 * a constraint that enforces the sum of these linear expressions
 * over the coefficients to be strictly greater than zero.
 * In this sum, some of the linear expressions may be negated to ensure
 * that a solution exists.
 *
 * The return value is a list of constraints, the first *orthonum corresponding
 * to the linear expressions that form a basis of the null space
 * and the final constraint the actual linear independence constraint.
 *
 * If the null space is 0-dimensional, *orthonum is zero and the return
 * value is NULL
 */
PlutoConstraints **get_stmt_ortho_constraints(Stmt *stmt, const PlutoProg *prog,
        HyperplaneProperties *hProps, const PlutoConstraints *currcst,
       int *orthonum)
{
    int i, j, k, p, q;
    PlutoConstraints **orthcst;
    isl_ctx *ctx;
    isl_int v;
    isl_mat *h;
    isl_basic_set *isl_currcst;

    int nvar = prog->nvar;
    int npar = prog->npar;
    int nstmts = prog->nstmts;

    if (stmt->num_ind_sols >= stmt->dim) {
        *orthonum = 0;
        return NULL;
    }

    /* Get rid of the variables that don't appear in the domain of this
     * statement and also beta rows
     */
    for (i = 0, p = 0; i < nvar; i++)
        if (stmt->is_outer_loop[i])
           p++;
    for (j = 0, q = 0; j < stmt->trans->nrows; j++)
       if (hProps[j].type != H_SCALAR)
           q++;

    if (q == 0) {
        /* no need to add any orthogonality constraints */
        *orthonum = 0;
        return NULL;
    }

    ctx = isl_ctx_alloc();
    assert(ctx);
    isl_int_init(v);

    h = isl_mat_alloc(ctx, q, p);

    p=0; 
    q=0;
    for (i=0; i<nvar; i++) {
        if (stmt->is_outer_loop[i])    {
            q=0;
            for (j=0; j<stmt->trans->nrows; j++) {
                /* Skip rows of h that are zero */
                if (hProps[j].type != H_SCALAR)   {
                   isl_int_set_si(v, stmt->trans->val[j][i]);
                   h = isl_mat_set_element(h, q, p, v);
                    q++;
                }
            }
            p++;
        }
    }

    h = isl_mat_right_kernel(h);

    PlutoMatrix *ortho = pluto_matrix_from_isl_mat(h);

    isl_mat_free(h);

    orthcst = (PlutoConstraints **) malloc(nvar*sizeof(PlutoConstraints *)); 

    for (i=0; i<nvar; i++)  {
        orthcst[i] = pluto_constraints_alloc(1, CST_WIDTH);
        orthcst[i]->ncols = CST_WIDTH;
    }

    /* Positive orthant only */
    /* An optimized version where the constraints are added as
     * c_1 >= 0, c_2 >= 0, ..., c_n >= 0, c_1+c_2+..+c_n >= 1
     *
     * basically only look in the orthogonal space where everything is
     * non-negative
     *
     * All of these constraints are added later to 
     * the global constraint matrix
     */

    /* Normalize ortho first */
    for (j=0; j<ortho->ncols; j++)    {
        if (ortho->val[0][j] == 0) continue;
        int colgcd = abs(ortho->val[0][j]);
        for (i=1; i<ortho->nrows; i++)    {
            if (ortho->val[i][j] == 0)  break;
            colgcd = gcd(colgcd,abs(ortho->val[i][j]));
        }
        if (i == ortho->nrows)   {
            if (colgcd > 1)    {
                for (k=0; k<ortho->nrows; k++)    {
                    ortho->val[k][j] /= colgcd;
                }
            }
        }
    }
    // pluto_matrix_print(stdout, ortho); 

    isl_currcst = isl_basic_set_from_pluto_constraints(ctx, currcst);

    assert(p == ortho->nrows);
    p=0;
    for (i=0; i<ortho->ncols; i++) {
        isl_basic_set *orthcst_i;

        j=0;
        for (q=0; q<nvar; q++) {
            if (stmt->is_outer_loop[q])    {
                orthcst[p]->val[0][npar+1+(stmt->id)*(nvar+1)+q] = ortho->val[j][i];
                j++;
            }
        }
        orthcst[p]->nrows = 1;
        orthcst[p]->val[0][CST_WIDTH-1] = -1;
        orthcst_i = isl_basic_set_from_pluto_constraints(ctx, orthcst[p]);
        orthcst[p]->val[0][CST_WIDTH-1] = 0;
        orthcst_i = isl_basic_set_intersect(orthcst_i,
                isl_basic_set_copy(isl_currcst));
        if (isl_basic_set_is_empty(orthcst_i))
            negate_constraint(orthcst[p]);
        isl_basic_set_free(orthcst_i);
        p++;
        assert(p<=nvar-1);
    }

    // pluto_matrix_print(stdout, stmt->trans);

    if (p > 0)  {
        /* Sum of all of the above is the last constraint */
        for(j=0; j<CST_WIDTH; j++)  {
            for (i=0; i<p; i++) {
                orthcst[p]->val[0][j] += orthcst[i]->val[0][j];
            }
        }
        orthcst[p]->nrows = 1;
        orthcst[p]->val[0][CST_WIDTH-1] = -1;
        p++;
    }

#if 0
    /* Since each of the ortho constraints is tried and the
     * best of the solutions will be kept; give all constraints for the 
     * statement
     * */

    p=0;
    for (i=0; i<ncols; i++) {
        for (j=0; j<ncols; j++) {
            if (ortho->val[j][i] != 0) break;
        }
        /* Ignore all zero cols */
        if (j==ncols) continue;

        /* We have a non-zero col */
        j=0;
        for (q=0; q<nvar; q++) {
            if (stmt->is_outer_loop[q])    {
                orthcst[p]->val[0][npar+1+(stmt->id)*(nvar+1)+q] = ortho->val[j][i];
                j++;
            }
        }
        orthcst[p]->nrows = 1;
        orthcst[p]->val[0][CST_WIDTH-1] = -1;
        p++;
    }
#endif

    *orthonum = p;

    /* Free the unnecessary ones */
    for (i=p; i<nvar; i++)    {
        pluto_constraints_free(orthcst[i]);
    }

    /* printf("Ortho constraints: %d\n", *orthonum); */
    // for (i=0; i<*orthonum; i++) {
        // IF_DEBUG2(constraints_print(stdout, orthcst[i]));
    // }

    pluto_matrix_free(ortho);
    isl_int_clear(v);
    isl_basic_set_free(isl_currcst);
    isl_ctx_free(ctx);

    return orthcst;
}


/*
 * Check whether the dependence is carried at level 'level'
 * (works whether the dep is const or non-const, inter-stmt or
 * self edge
 */
bool dep_satisfaction_test(Dep *dep, PlutoProg *prog, int level)
{
    static PlutoConstraints *cst = NULL;
    int i, j, src, dest, *sol;

    int nvar = prog->nvar;
    int npar = prog->npar;

    Stmt *stmts = prog->stmts;

    src = dep->src;
    dest = dep->dest;

    assert(level < stmts[src].trans->nrows);
    assert(level < stmts[dest].trans->nrows);

    if (!cst || cst->alloc_nrows < 1+dep->dpolytope->nrows)   {
        if (cst) pluto_constraints_free(cst);
        /* rougly allocate twice to prevent frequent increase */
        cst = pluto_constraints_alloc(2*(1+dep->dpolytope->nrows), 2*nvar+npar+1);
    }
    cst->ncols = 2*nvar+npar+1;

    /*
     * constraint format 
     * \phi(src) - \phi (dest) >= 0
     * (reverse of satisfaction)
     */

    for (j=0; j<nvar; j++)    {
        cst->val[0][j] = stmts[src].trans->val[level][j];
    }
    for (j=nvar; j<2*nvar; j++)    {
        cst->val[0][j] = -stmts[dest].trans->val[level][j-nvar];
    }
    cst->val[0][2*nvar+npar] = 
        stmts[src].trans->val[level][nvar] - stmts[dest].trans->val[level][nvar];

    for (i=0; i<dep->dpolytope->nrows; i++)  {
        for (j=0; j<2*nvar+npar+1; j++)  {
            cst->val[1+i][j] = dep->dpolytope->val[i][j];
        }
    }

    cst->nrows = 1+dep->dpolytope->nrows;

    /* if no solution exists, the dependence is carried, i.e., no points
     * satisfy \geq 0 */ 
    sol = pluto_constraints_solve(cst);

    bool retval = (sol)? false:true;
    free(sol);

    return retval;
}


int get_dep_direction(const Dep *dep, const PlutoProg *prog, int level)
{
    static PlutoConstraints *cst = NULL;
    int i, j, src, dest;

    int nvar = prog->nvar;
    int npar = prog->npar;
    Stmt *stmts = prog->stmts;

    src = dep->src;
    dest = dep->dest;

    assert(level < stmts[src].trans->nrows);
    assert(level < stmts[dest].trans->nrows);

    if (!cst || cst->alloc_nrows < 2+dep->dpolytope->nrows)   {
        if (cst) pluto_constraints_free(cst);
        /* Rougly allocate twice to prevent frequent increase */
        cst = pluto_constraints_alloc(2*(2+dep->dpolytope->nrows), 2*nvar+npar+1);
    }
    cst->ncols = 2*nvar+npar+1;

    /*
     * Check for zero
     *
     * To test \phi (dest) - \phi(src) = 0, we try 
     *
     * \phi(dest) - \phi(src) >= 1
     */

    for (j=0; j<nvar; j++)    {
        cst->val[0][j] = -stmts[src].trans->val[level][j];
    }
    for (j=nvar; j<2*nvar; j++)    {
        cst->val[0][j] = stmts[dest].trans->val[level][j-nvar];
    }
    cst->val[0][2*nvar+npar] = 
        -stmts[src].trans->val[level][nvar] + stmts[dest].trans->val[level][nvar]-1;

    for (i=0; i<dep->dpolytope->nrows; i++)  {
        for (j=0; j<2*nvar+npar+1; j++)  {
            cst->val[1+i][j] = dep->dpolytope->val[i][j];
        }
    }

    cst->nrows = 1+dep->dpolytope->nrows;

    int *sol = pluto_constraints_solve(cst);

    if (!sol)   {
        for (j=0; j<nvar; j++)    {
            cst->val[0][j] = stmts[src].trans->val[level][j];
        }
        for (j=nvar; j<2*nvar; j++)    {
            cst->val[0][j] = -stmts[dest].trans->val[level][j-nvar];
        }
        cst->val[0][2*nvar+npar] = 
            stmts[src].trans->val[level][nvar] - stmts[dest].trans->val[level][nvar]-1;

        for (i=0; i<dep->dpolytope->nrows; i++)  {
            for (j=0; j<2*nvar+npar+1; j++)  {
                cst->val[1+i][j] = dep->dpolytope->val[i][j];
            }
        }

        cst->nrows = 1+dep->dpolytope->nrows;

        sol = pluto_constraints_solve(cst);

        /* If no solution exists, all points satisfy \phi (dest) - \phi (src) = 0 */
        if (!sol)   {
            return DEP_ZERO;
        }
    }


    /*
     * Check for PLUS
     * Constraint format 
     * \phi(dest) - \phi (src) <= -1
     * (reverse of plus)
     */

    for (j=0; j<nvar; j++)    {
        cst->val[0][j] = stmts[src].trans->val[level][j];
    }
    for (j=nvar; j<2*nvar; j++)    {
        cst->val[0][j] = -stmts[dest].trans->val[level][j-nvar];
    }
    cst->val[0][2*nvar+npar] = 
        stmts[src].trans->val[level][nvar] - stmts[dest].trans->val[level][nvar] -1;

    for (i=0; i<dep->dpolytope->nrows; i++)  {
        for (j=0; j<2*nvar+npar+1; j++)  {
            cst->val[1+i][j] = dep->dpolytope->val[i][j];
        }
    }

    cst->nrows = 1+dep->dpolytope->nrows;

    free(sol);
    sol = pluto_constraints_solve(cst);

    if (!sol)   {
        return DEP_PLUS;
    }

    /*
     * Check for MINUS
     *
     * Constraint format 
     * \phi(dest) - \phi (src) >= 1
     * reverse of minus, we alraedy know that it's not zero
     */

    for (j=0; j<nvar; j++)    {
        cst->val[0][j] = -stmts[src].trans->val[level][j];
    }
    for (j=nvar; j<2*nvar; j++)    {
        cst->val[0][j] = stmts[dest].trans->val[level][j-nvar];
    }
    cst->val[0][2*nvar+npar] = 
        -stmts[src].trans->val[level][nvar] + stmts[dest].trans->val[level][nvar] -1;

    for (i=0; i<dep->dpolytope->nrows; i++)  {
        for (j=0; j<2*nvar+npar+1; j++)  {
            cst->val[1+i][j] = dep->dpolytope->val[i][j];
        }
    }

    cst->nrows = 1+dep->dpolytope->nrows;

    free(sol);
    sol = pluto_constraints_solve(cst);

    if (!sol)   {   
        return DEP_MINUS;
    }

    free(sol);

    /* Neither ZERO, nor PLUS, nor MINUS, has to be STAR */
    return DEP_STAR;
}
