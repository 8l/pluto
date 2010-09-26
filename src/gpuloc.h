#ifndef _GPULOC_H
#define _GPULOC_H

#include <barvinok/barvinok.h>
#include <cloog/isl/cloog.h>

#include "pluto.h"

struct array_info {
    /* Name of the array. */
    char *name;
    /* Number of indices. */
    int dim;
    /* Number of elements accessed per iteration of host loops. */
    isl_pw_qpolynomial *transfer_size;
    /* Amount of shared memory required per block. */
    isl_pw_qpolynomial_fold *shared_size;
};

struct localizer_info {
    /* Number of rows in the original schedule computed by Pluto. */
    int len;
    /* First tile dimension. */
    int first;
    /* Number of tile dimensions. */
    int tile_len;
    /* Number of rows in the tiled schedule. */
    int tiled_len;
    /* Number of loops on the host. */
    int host_len;
    /* Number of loops up to sequential loop on GPU. */
    int gpu_len;
    /* Global, tiled, schedule. */
    isl_union_map *sched;
    PlutoProg *prog;
    CloogState *state;
    int indent;
    FILE *input;
    FILE *host_c;
    FILE *kernel_h;
    FILE *kernel_c;
    FILE *dst;

    const char *type;

    /* Identifier of current kernel. */
    int kernel_id;

    /* Number of arrays. */
    int n_array;
    struct array_info *array;
    /* Total number of array elements accessed per iteration of host loops. */
    isl_pw_qpolynomial *transfer_size;
    /* Maximal number of array elements accessed
     * over all iterations of host loops.
     */
    isl_pw_qpolynomial_fold *max_transfer_size;
    /* Maximal total amount of shared memory required. */
    isl_pw_qpolynomial_fold *max_shared_size;

    int upper;
    int n;
};

int gpuloc(PlutoProg *prog, PlutoOptions *options, const char *input);

void gpu_print_host_user(struct localizer_info *loc,
    struct clast_user_stmt *u);

#endif
