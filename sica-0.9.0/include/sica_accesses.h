/*
 * sica_accesses.h
 *
 *  Created on: 04.03.2013
 *      Author: dfeld
 */


#ifndef SICA_ACCESSES_H_
#define SICA_ACCESSES_H_

#include "pluto.h"

struct sica_stmt_list{
  int stmt_nb;

  struct sica_stmt_list* next; /**< Next statement in the linked list */
};
typedef struct sica_stmt_list SICAStmtList;

struct sica_access_matrices{
  int** access_mat; //array that stores the transformated accesses that are accessed by the vectorized dimension
  int nrows;
  int ncols;

  struct sica_access_matrices* next; /**< Next statement in the linked list */
};
typedef struct sica_access_matrices SICAAccess;

SICAAccess* sica_accesses_malloc();

int** sica_access_matrix_malloc(int rows, int cols);

void sica_copy_access_matrix(int** matrix1, int** matrix2, int rows, int cols);

int sica_compare_access_matrices(int** matrix1, int** matrix2, int nrows, int ncols);

int sica_check_for_entry(int** matrix, int rows, int cols);

void sica_vec_times_matrix(int* solution_vec, int* vector, int** matrix, int rows, int columns);

int sica_get_array_id(Band* act_band, char* array_name);

void sica_set_vectorized_bands(Band** bands, int nbands);

void sica_get_trans_matrix(Band** bands, int nbands);

#endif /* SICA_ACCESSES_H_ */
