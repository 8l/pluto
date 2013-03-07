/*
 * sica_math_func.c
 *
 *  Created on: 25.02.2013
 *      Author: dfeld
 */

#ifndef SICA_MATH_FUNC_H_
#define SICA_MATH_FUNC_H_

void sica_mult_matrices(float* mult1, float* mult2, float* res, int rows1, int columns1, int rows2, int columns2);

int sica_minusone_pow(int exponent);

void sica_print_matrix(float** matrix, int rows, int columns);

void sica_cpy_matrix(float** matrix, float** temp_matrix, int rows, int columns);

int sica_echelon_form(float** temp_matrix, int rows, int columns);

float sica_echelon_determinant(float** temp_matrix, int N);

float sica_determinant(float** matrix, int N);

void sica_setup_identitymatrix(float** ident_matrix,int N);

void sica_echelon_inverse(float** inverse_matrix, float** matrix,int N);

void sica_inverse(int** matrix, int** inverse_matrix, int N);

#endif /* SICA_MATH_FUNC_H_ */
