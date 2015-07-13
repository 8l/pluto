/*
 * PLUTO: An automatic parallelier and locality optimizer
 * 
 * Copyright (C) 2007-2008 Uday Bondhugula
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * A copy of the GNU General Public Licence can be found in the 
 * top-level directory of this program (`COPYING') 
 *
 */
#ifndef _MATH_SUPPORT_H
#define _MATH_SUPPORT_H

#include <stdio.h>

#include "isl/mat.h"
#include "gmp.h"

#define PLMAX(a,b) ((a>=b)?(a):(b))
#define PLMIN(a,b) ((a<=b)?(a):(b))
#define PLABS(a) ((a>=0)?(a):(-a))

#define int64 long long int
#define LONG_LONG_INT_MAX 0x7FFFFFFFFFFFFFFFL

/* A matrix */
struct plutoMatrix{
    /* The values */
    int64 **val;

    int nrows;
    int ncols;

    /* Pre-allocated number of rows */
    int alloc_nrows;
    int alloc_ncols;
};
typedef struct plutoMatrix PlutoMatrix;

void pluto_matrix_print(FILE *, const PlutoMatrix *);
void pluto_matrix_read(FILE *, const PlutoMatrix *);
PlutoMatrix *pluto_matrix_alloc(int nrows, int ncols);
void pluto_matrix_free(PlutoMatrix *mat);
PlutoMatrix *pluto_matrix_dup(const PlutoMatrix *src);
PlutoMatrix *pluto_matrix_identity(int size);
void pluto_matrix_set(PlutoMatrix *mat, int val);
PlutoMatrix *pluto_matrix_input(FILE *fp);

PlutoMatrix *pluto_matrix_inverse(PlutoMatrix *mat);
PlutoMatrix *pluto_matrix_product(const PlutoMatrix *mat1, 
        const PlutoMatrix *mat2);
int pluto_matrix_get_rank(const PlutoMatrix *mat);

void pluto_matrix_add_row(PlutoMatrix *mat, int pos);
void pluto_matrix_add_col(PlutoMatrix *mat, int pos);
void pluto_matrix_remove_row(PlutoMatrix *mat, int pos);
void pluto_matrix_remove_col(PlutoMatrix *, int);
void pluto_matrix_zero_row(PlutoMatrix *mat, int pos);
void pluto_matrix_zero_col(PlutoMatrix *mat, int pos);

void pluto_matrix_move_col(PlutoMatrix *mat, int r1, int r2);
void pluto_matrix_interchange_cols(PlutoMatrix *mat, int c1, int c2);
void pluto_matrix_interchange_rows(PlutoMatrix *mat, int r1, int r2);

void pluto_matrix_normalize_row(PlutoMatrix *mat, int pos);
void pluto_matrix_negate_row(PlutoMatrix *mat, int pos);
void pluto_matrix_add(PlutoMatrix *mat1, const PlutoMatrix *mat2);
void gaussian_eliminate(PlutoMatrix *mat, int start, int end);

int64 lcm(int64 a, int64 b);
int64 gcd(int64 a, int64 b);
int64 *min_lexical(int64 *a, int64 *b, int64 num);

char *concat(const char *prefix, const char *suffix);
void pluto_affine_function_print(FILE *fp, int64 *func, int ndims, char **vars);

void pluto_matrix_reverse_rows(PlutoMatrix *mat);
void pluto_matrix_negate(PlutoMatrix *mat);

int pluto_vector_is_parallel(PlutoMatrix *mat1, int r1, PlutoMatrix *mat2, int r2);
int pluto_vector_is_normal(PlutoMatrix *mat1, int r1, PlutoMatrix *mat2, int r2);

PlutoMatrix *pluto_matrix_from_isl_mat(__isl_keep isl_mat *mat);

long long isl_val_get_num_ll(__isl_keep isl_val *v);
void mpz_set_sll(mpz_t n, long long sll);

#endif
