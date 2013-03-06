/*
 * cache_math_func.c
 *
 *  Created on: 25.02.2013
 *      Author: dfeld
 */

#ifndef CACHE_MATH_FUNC_H_
#define CACHE_MATH_FUNC_H_

void cache_mult_matrices(float* mult1, float* mult2, float* res, int rows1, int columns1, int rows2, int columns2);

int cache_minusone_pow(int exponent);

void cache_print_matrix(float* matrix, int rows, int columns);

void cache_cpy_matrix(float* matrix, float* temp_matrix, int rows, int columns);

int cache_echelon_form(float* temp_matrix, int rows, int columns);

float cache_echelon_determinant(float* temp_matrix, int N);

float cache_determinant(float* matrix, int N);

void cache_setup_identitymatrix(float* ident_matrix,int N);

void cache_echelon_inverse(float* inverse_matrix, float* matrix,int N);

void cache_inverse(float* matrix, float* inverse_matrix, int N);

void cache_vec_times_matrix(float* solution_vec, float* vector, float* matrix, int rows, int columns);

#endif /* CACHE_MATH_FUNC_H_ */
