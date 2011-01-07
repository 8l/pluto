#include <assert.h>

#include <isl/polynomial.h>
#include <isl/union_set.h>
#include <cloog/isl/cloog.h>

#include "cuda.h"
#include "cuda_common.h"
#include "gpucode.h"
#include "program.h"
#include "schedule.h"

struct cuda_array_bound {
	isl_int size;
	isl_qpolynomial *lb;
};

struct cuda_array_info {
	isl_dim *dim;
	/* Name of the array. */
	char *name;
	/* Number of indices. */
	unsigned n_index;
	/* For each index, a bound on the array in that direction. */
	isl_pw_qpolynomial_fold **bound;
	/* For each index, bound[i] specialized to the current kernel. */
	isl_pw_qpolynomial_fold **local_bound;
	/* Is the array used in the current kernel? */
	int active;
	/* For each index, size and offset of piece in shared memory. */
	struct cuda_array_bound *shared_bound;
};

struct cuda_gen {
	struct cuda_info cuda;
	struct gpucode_info code;
	struct gpucode_info kernel_code;

	isl_ctx *ctx;
	PlutoProg *prog;
	PlutoOptions *options;
	CloogState *state;

	int n_array;
	struct cuda_array_info *array;

	/* Identifier of current kernel. */
	int kernel_id;

	/* First tile dimension. */
	int tile_first;
	/* Number of tile dimensions. */
	int tile_len;

	/* Number of dimensions determining shared memory. */
	int shared_len;

	/* Number of rows in the untiled schedule. */
	int untiled_len;
	/* Number of rows in the tiled schedule. */
	int tiled_len;
	/* Number of rows in schedule after tiling/wrapping over threads. */
	int thread_tiled_len;

	/* Global untiled schedule. */
	isl_union_map *sched;
	/* Local tiled schedule. */
	isl_union_map *tiled_sched;

	int n_grid;
	int n_block;
	int grid_dim[2];
	int block_dim[3];
	int *tile_size;
};

/* Compute bounds on the host arrays based on the accessed elements.
 */
static int extract_array_info(__isl_take isl_set *array, void *user)
{
	int i;
	struct cuda_gen *gen = (struct cuda_gen *)user;
	const char *name;
	int n_index;
	isl_pw_qpolynomial_fold **bounds;
	isl_pw_qpolynomial_fold **local_bounds;
	struct cuda_array_bound *shared_bound;

	n_index = isl_set_dim(array, isl_dim_set);
	name = isl_set_get_tuple_name(array);
	bounds = isl_alloc_array(isl_set_get_ctx(array),
				 isl_pw_qpolynomial_fold *, n_index);
	assert(bounds);
	local_bounds = isl_calloc_array(isl_set_get_ctx(array),
				 isl_pw_qpolynomial_fold *, n_index);
	assert(local_bounds);
	shared_bound = isl_alloc_array(isl_set_get_ctx(array),
				 struct cuda_array_bound, n_index);
	assert(shared_bound);
	gen->array[gen->n_array].dim = isl_set_get_dim(array);
	gen->array[gen->n_array].name = strdup(name);
	gen->array[gen->n_array].n_index = n_index;
	gen->array[gen->n_array].bound = bounds;
	gen->array[gen->n_array].local_bound = local_bounds;
	gen->array[gen->n_array].shared_bound = shared_bound;

	for (i = 0; i < n_index; ++i) {
		isl_int_init(shared_bound[i].size);
		shared_bound[i].lb = NULL;
	}

	for (i = 0; i < n_index; ++i) {
		isl_dim *dim;
		isl_qpolynomial *qp, *one;
		isl_pw_qpolynomial *pwqp;
		isl_pw_qpolynomial_fold *pwf;

		dim = isl_set_get_dim(array);
		one = isl_qpolynomial_one(isl_dim_copy(dim));
		qp = isl_qpolynomial_var(dim, isl_dim_set, i);
		qp = isl_qpolynomial_add(qp, one);
		pwqp = isl_pw_qpolynomial_alloc(isl_set_copy(array), qp);
		pwf = isl_pw_qpolynomial_bound(pwqp, isl_fold_max, NULL);

		bounds[i] = pwf;
	}

	gen->n_array++;

	isl_set_free(array);
	return 0;
}

static void collect_array_info(struct cuda_gen *gen, PlutoProg *prog)
{
	isl_union_set *arrays;

	arrays = isl_union_map_range(isl_union_map_copy(prog->read));
	arrays = isl_union_set_union(arrays,
			isl_union_map_range(isl_union_map_copy(prog->write)));
	arrays = isl_union_set_coalesce(arrays);

	gen->n_array = isl_union_set_n_set(arrays);
	gen->array = isl_alloc_array(prog->ctx,
				     struct cuda_array_info, gen->n_array);
	assert(gen->array);
	gen->n_array = 0;
	isl_union_set_foreach_set(arrays, &extract_array_info, gen);
	isl_union_set_free(arrays);
}

static void free_array_info(struct cuda_gen *gen)
{
	int i, j;

	for (i = 0; i < gen->n_array; ++i) {
		free(gen->array[i].name);
		for (j = 0; j < gen->array[i].n_index; ++j) {
			isl_pw_qpolynomial_fold_free(gen->array[i].bound[j]);
			isl_int_clear(gen->array[i].shared_bound[j].size);
		}
		isl_dim_free(gen->array[i].dim);
		free(gen->array[i].bound);
		free(gen->array[i].shared_bound);
	}
	free(gen->array);
}

static void declare_device_arrays(struct cuda_gen *gen)
{
	int i;

	for (i = 0; i < gen->n_array; ++i)
		fprintf(gen->cuda.host_c, "%s *dev_%s;\n",
			gen->options->type, gen->array[i].name);
}

static void print_array_size(struct cuda_gen *gen, FILE *out,
	struct cuda_array_info *array)
{
	int i;
	isl_printer *prn;

	prn = isl_printer_to_file(gen->ctx, out);
	prn = isl_printer_set_output_format(prn, ISL_FORMAT_C);
	for (i = 0; i < array->n_index; ++i) {
		prn = isl_printer_print_str(prn, "(");
		prn = isl_printer_print_pw_qpolynomial_fold(prn,
							    array->bound[i]);
		prn = isl_printer_print_str(prn, ") * ");
	}
	prn = isl_printer_print_str(prn, "sizeof(");
	prn = isl_printer_print_str(prn, gen->options->type);
	prn = isl_printer_print_str(prn, ")");
	isl_printer_free(prn);
}

static void allocate_device_arrays(struct cuda_gen *gen)
{
	int i;

	for (i = 0; i < gen->n_array; ++i) {
		fprintf(gen->cuda.host_c, "cudaMalloc(&dev_%s, ", 
			gen->array[i].name);
		print_array_size(gen, gen->cuda.host_c, &gen->array[i]);
		fprintf(gen->cuda.host_c, ");\n"); 
	}
}

static void copy_arrays_to_device(struct cuda_gen *gen)
{
	int i;

	for (i = 0; i < gen->n_array; ++i) {
		fprintf(gen->cuda.host_c, "assert(sizeof(%s) == ",
			gen->array[i].name);
		print_array_size(gen, gen->cuda.host_c, &gen->array[i]);
		fprintf(gen->cuda.host_c, ");\n"); 
		fprintf(gen->cuda.host_c, "cudaMemcpy(dev_%s, %s, ",
			gen->array[i].name, gen->array[i].name);
		print_array_size(gen, gen->cuda.host_c, &gen->array[i]);
		fprintf(gen->cuda.host_c, ", cudaMemcpyHostToDevice);\n");
	}
}

static void copy_arrays_from_device(struct cuda_gen *gen)
{
	int i;

	for (i = 0; i < gen->n_array; ++i) {
		fprintf(gen->cuda.host_c, "cudaMemcpy(%s, dev_%s, ",
			gen->array[i].name, gen->array[i].name);
		print_array_size(gen, gen->cuda.host_c, &gen->array[i]);
		fprintf(gen->cuda.host_c, ", cudaMemcpyDeviceToHost);\n");
	}
}

static void select_tile_dimensions(struct cuda_gen *gen, PlutoProg *prog)
{
	int last;

	getOutermostTilableBand(prog, &gen->tile_first, &last);
	gen->tile_len = last - gen->tile_first + 1;

	gen->untiled_len = prog->num_hyperplanes;
}

/* Check if the first loop in the select tilable band is parallel.
 * If not, apply a wavefront transformation and drop the corresponding
 * loop from the tilable band.
 */
static void compute_global_schedule(struct cuda_gen *gen, PlutoProg *prog,
	PlutoOptions *options)
{
	isl_union_map *schedule;

	schedule = extract_schedule(prog);
	if (prog->hProps[gen->tile_first].dep_prop != PARALLEL) {
		int i;
		isl_dim *dim;
		isl_map *front;

		dim = isl_union_map_get_dim(schedule);
		front = wavefront(dim, gen->untiled_len,
				  gen->tile_first, gen->tile_len);
		schedule = isl_union_map_apply_range(schedule,
						 isl_union_map_from_map(front));
		gen->tile_first += 1;
		gen->tile_len -= 1;
		for (i = 0; i < gen->tile_len; ++i)
			prog->hProps[gen->tile_first + i].dep_prop = PARALLEL;
	}

	gen->sched = schedule;
}

static void read_sizes_from_file(struct cuda_gen *gen, const char *filename,
	int *sizes, int len)
{
	int i;
	FILE *file;

	file = fopen(filename, "r");
	if (!file)
		return;

	for (i = 0; i < len; ++i)
		if (fscanf(file, "%d", &sizes[i]) < 1)
			break;

	fclose(file);
}

/* Read user specified sizes from "tile.sizes", "block.sizes" and "grid.sizes"
 * after filling in some potentially useful defaults.
 */
static void read_sizes(struct cuda_gen *gen)
{
	int n;
	HyperplaneProperties *hProps = gen->prog->hProps;

	gen->tile_size = isl_alloc_array(gen->ctx, int, gen->tile_len);
	assert(gen->tile_size);
	for (n = 0; n < gen->tile_len; ++n)
		gen->tile_size[n] = gen->options->tile_size;
	read_sizes_from_file(gen, "tile.sizes", gen->tile_size, gen->tile_len);

	for (n = 0; n < gen->tile_len; ++n)
		if (hProps[gen->tile_first + n].dep_prop != PARALLEL)
			break;

	gen->n_block = (n <= 3) ? n : 3;
	switch (gen->n_block) {
	case 1:
		gen->block_dim[0] = 512;
		break;
	case 2:
		gen->block_dim[0] = 32;
		gen->block_dim[1] = 16;
		break;
	default:
		gen->block_dim[0] = 32;
		gen->block_dim[1] = 4;
		gen->block_dim[2] = 4;
		break;
	}
	read_sizes_from_file(gen, "block.sizes", gen->block_dim, gen->n_block);

	gen->n_grid = (n <= 2) ? n : 2;
	switch (gen->n_grid) {
	case 1:
		gen->grid_dim[0] = 65536;
		break;
	default:
		gen->grid_dim[0] = 256;
		gen->grid_dim[1] = 256;
		break;
	}
	read_sizes_from_file(gen, "grid.sizes", gen->grid_dim, gen->n_grid);
}

static void clear_cuda_gen(struct cuda_gen *gen)
{
	free_array_info(gen);
	free(gen->tile_size);
	isl_union_map_free(gen->sched);
}

static void declare_host_iterators(struct cuda_gen *gen)
{
	int i;

	if (gen->tile_first == 0)
		return;

	fprintf(gen->cuda.host_c, "int ");
	for (i = 0; i < gen->tile_first; ++i) {
		if (i)
			fprintf(gen->cuda.host_c, ", ");
		fprintf(gen->cuda.host_c, "h%d", i);
	}
	fprintf(gen->cuda.host_c, ";\n");
}

static void print_list(FILE *out, int len, int *list)
{
	int i;

	for (i = 0; i < len; ++i) {
		if (i)
			fprintf(out, ", ");
		fprintf(out, "%d", list[i]);
	}
}

static void print_kernel_launch(struct cuda_gen *gen,
	__isl_keep isl_union_set *arrays)
{
	int i;
	int first = 1;
	unsigned nparam;
	isl_dim *dim;

	print_indent(gen->code.dst, gen->code.indent);
	fprintf(gen->code.dst, "kernel%d <<<k%d_dimGrid, k%d_dimBlock>>> (",
		gen->kernel_id, gen->kernel_id, gen->kernel_id);
	fprintf(gen->cuda.kernel_c, "__global__ void kernel%d(",
		gen->kernel_id);
	fprintf(gen->cuda.kernel_h, "__global__ void kernel%d(",
		gen->kernel_id);

	for (i = 0; i < gen->n_array; ++i) {
		isl_dim *dim;
		isl_set *arr;
		int empty;

		dim = isl_dim_copy(gen->array[i].dim);
		arr = isl_union_set_extract_set(arrays, dim);
		empty = isl_set_fast_is_empty(arr);
		isl_set_free(arr);
		if (empty)
			continue;

		if (!first) {
			fprintf(gen->code.dst, ", ");
			fprintf(gen->cuda.kernel_c, ", ");
			fprintf(gen->cuda.kernel_h, ", ");
		}

		fprintf(gen->code.dst, "dev_%s", gen->array[i].name);
		fprintf(gen->cuda.kernel_c, "%s *%s",
			gen->options->type, gen->array[i].name);
		fprintf(gen->cuda.kernel_h, "%s *%s",
			gen->options->type, gen->array[i].name);

		first = 0;
	}

	dim = isl_union_set_get_dim(arrays);
	nparam = isl_dim_size(dim, isl_dim_param);
	for (i = 0; i < nparam; ++i) {
		const char *name = isl_dim_get_name(dim, isl_dim_param, i);
		if (!first) {
			fprintf(gen->code.dst, ", ");
			fprintf(gen->cuda.kernel_c, ", ");
			fprintf(gen->cuda.kernel_h, ", ");
		}
		fprintf(gen->code.dst, "%s", name);
		fprintf(gen->cuda.kernel_c, "int %s", name);
		fprintf(gen->cuda.kernel_h, "int %s", name);
		first = 0;
	}
	isl_dim_free(dim);

	for (i = 0; i < gen->tile_first; ++i) {
		if (!first) {
			fprintf(gen->code.dst, ", ");
			fprintf(gen->cuda.kernel_c, ", ");
			fprintf(gen->cuda.kernel_h, ", ");
		}
		fprintf(gen->code.dst, "h%d", i);
		fprintf(gen->cuda.kernel_c, "int h%d", i);
		fprintf(gen->cuda.kernel_h, "int h%d", i);
		first = 0;
	}

	fprintf(gen->code.dst, ");\n");
	fprintf(gen->cuda.kernel_c, ")\n");
	fprintf(gen->cuda.kernel_h, ");\n");
}

/* Construct a map from a domain of dimensionality "len"
 * to a domain of dimensionality "len" + "tile_len" that tiles
 * the "tile_len" coordinates starting at "first".
 * In particular, [s_i] -> [s_i / tile_size[i], s_i % tile_size[i]].
 * "dim" prescribes the parameters.
 */
static __isl_give isl_map *tile(__isl_take isl_dim *dim, int len,
        int first, int tile_len, int *tile_size)
{
	int i;
	isl_int v;
	isl_basic_map *bmap;
	isl_constraint *c;

	isl_int_init(v);

	dim = isl_dim_add(dim, isl_dim_in, len);
	dim = isl_dim_add(dim, isl_dim_out, len + tile_len);
	bmap = isl_basic_map_universe(isl_dim_copy(dim));

	for (i = 0; i < len - tile_len; ++i) {
		int j = i < first ? i : i + tile_len;
		int k = i < first ? i : i + 2 * tile_len;

		c = isl_equality_alloc(isl_dim_copy(dim));
		isl_int_set_si(v, -1);
		isl_constraint_set_coefficient(c, isl_dim_in, j, v);
		isl_int_set_si(v, 1);
		isl_constraint_set_coefficient(c, isl_dim_out, k, v);
		bmap = isl_basic_map_add_constraint(bmap, c);
	}

	for (i = 0; i < tile_len; ++i) {
		c = isl_equality_alloc(isl_dim_copy(dim));
		isl_int_set_si(v, -1);
		isl_constraint_set_coefficient(c, isl_dim_in, first + i, v);
		isl_int_set_si(v, tile_size[i]);
		isl_constraint_set_coefficient(c, isl_dim_out, first + i, v);
		isl_int_set_si(v, 1);
		isl_constraint_set_coefficient(c, isl_dim_out,
						first + i + tile_len, v);
		bmap = isl_basic_map_add_constraint(bmap, c);

		c = isl_inequality_alloc(isl_dim_copy(dim));
		isl_int_set_si(v, 1);
		isl_constraint_set_coefficient(c, isl_dim_out,
						first + i + tile_len, v);
		bmap = isl_basic_map_add_constraint(bmap, c);
	
		c = isl_inequality_alloc(isl_dim_copy(dim));
		isl_int_set_si(v, -1);
		isl_constraint_set_coefficient(c, isl_dim_out,
						first + i + tile_len, v);
		isl_int_set_si(v, tile_size[i] - 1);
		isl_constraint_set_constant(c, v);
		bmap = isl_basic_map_add_constraint(bmap, c);
	}

	isl_dim_free(dim);
	isl_int_clear(v);

	return isl_map_from_basic_map(bmap);
}

/* Construct a map from a domain of dimensionality "len"
 * to a domain of dimensionality "len" + "wrap_len" that "wraps"
 * the "wrap_len" coordinates starting at "first" according to "wrap_size".
 * In particular, [s_i] -> [s_i, s_i % wrap_size[i]].
 * To do so, we need extra variables corresponding to [s_i / wrap_size[i]],
 * that are projected out at the end.
 * "dim" prescribes the parameters.
 */
static __isl_give isl_map *wrap(__isl_take isl_dim *dim, int len,
        int first, int wrap_len, int *wrap_size)
{
	int i;
	isl_basic_map *bmap;
	isl_constraint *c;

	dim = isl_dim_add(dim, isl_dim_in, len);
	dim = isl_dim_add(dim, isl_dim_out, len + 2 * wrap_len);
	bmap = isl_basic_map_universe(isl_dim_copy(dim));

	for (i = 0; i < len; ++i) {
		int k = i < first + wrap_len ? i : i + 2 * wrap_len;

		c = isl_equality_alloc(isl_dim_copy(dim));
		isl_constraint_set_coefficient_si(c, isl_dim_in, i, -1);
		isl_constraint_set_coefficient_si(c, isl_dim_out, k, 1);
		bmap = isl_basic_map_add_constraint(bmap, c);
	}

	for (i = 0; i < wrap_len; ++i) {
		c = isl_equality_alloc(isl_dim_copy(dim));
		isl_constraint_set_coefficient_si(c, isl_dim_out,
						    first + i, -1);
		isl_constraint_set_coefficient_si(c, isl_dim_out,
						    first + wrap_len + i, 1);
		isl_constraint_set_coefficient_si(c, isl_dim_out,
				    first + 2 * wrap_len + i, wrap_size[i]);
		bmap = isl_basic_map_add_constraint(bmap, c);

		c = isl_inequality_alloc(isl_dim_copy(dim));
		isl_constraint_set_coefficient_si(c, isl_dim_out,
						    first + wrap_len + i, 1);
		bmap = isl_basic_map_add_constraint(bmap, c);
	
		c = isl_inequality_alloc(isl_dim_copy(dim));
		isl_constraint_set_coefficient_si(c, isl_dim_out,
						    first + wrap_len + i, -1);
		isl_constraint_set_constant_si(c, wrap_size[i] - 1);
		bmap = isl_basic_map_add_constraint(bmap, c);
	}

	isl_dim_free(dim);

	bmap = isl_basic_map_project_out(bmap, isl_dim_out,
				first + 2 * wrap_len, wrap_len);

	return isl_map_from_basic_map(bmap);
}

/* Equate the "n" dimensions of "set" starting at "first" to
 * freshly created parameters named prefix%d.
 */
static __isl_give isl_set *parametrize(__isl_take isl_set *set,
	int first, int n, const char *prefix)
{
	int i;
	unsigned nparam;
	isl_int v;
	isl_dim *dim;
	isl_basic_set *bset;
	isl_constraint *c;
	char name[20];

	nparam = isl_set_dim(set, isl_dim_param);
	set = isl_set_add_dims(set, isl_dim_param, n);

	for (i = 0; i < n; ++i) {
		snprintf(name, sizeof(name), "%s%d", prefix, i);
		set = isl_set_set_dim_name(set, isl_dim_param,
					    nparam + i, name);
	}

	dim = isl_set_get_dim(set);
	bset = isl_basic_set_universe(isl_dim_copy(dim));

	isl_int_init(v);

	for (i = 0; i < n; ++i) {
		c = isl_equality_alloc(isl_dim_copy(dim));
		isl_int_set_si(v, -1);
		isl_constraint_set_coefficient(c, isl_dim_param, nparam + i, v);
		isl_int_set_si(v, 1);
		isl_constraint_set_coefficient(c, isl_dim_set, first + i, v);
		bset = isl_basic_set_add_constraint(bset, c);
	}

	isl_int_clear(v);
	isl_dim_free(dim);

	return isl_set_intersect(set, isl_set_from_basic_set(bset));
}

static __isl_give isl_set *parametrization(__isl_take isl_dim *dim,
	int len, int first, int n, const char *prefix)
{
	isl_set *set;

	dim = isl_dim_add(dim, isl_dim_set, len);
	set = isl_set_universe(dim);

	return parametrize(set, first, n, prefix);
}

/* Tile the B loops over the tile sizes and then tile/wrap
 * the T1 loops over the blocks.
 */
static __isl_give isl_union_map *tile_schedule(struct cuda_gen *gen,
	__isl_take isl_union_map *sched)
{
	isl_dim *dim;
	isl_map *tiling, *block_tiling;

	dim = isl_union_map_get_dim(sched);
	tiling = tile(isl_dim_copy(dim), gen->untiled_len,
		      gen->tile_first, gen->tile_len, gen->tile_size);

	if (gen->options->cuda_wrap)
		block_tiling = wrap(dim, gen->untiled_len + gen->tile_len,
				gen->tile_first, gen->n_grid, gen->grid_dim);
	else
		block_tiling = tile(dim, gen->untiled_len + gen->tile_len,
				gen->tile_first, gen->n_grid, gen->grid_dim);

	gen->tiled_len = gen->untiled_len + gen->tile_len + gen->n_grid;

	tiling = isl_map_apply_range(tiling, block_tiling);

	sched = isl_union_map_apply_range(sched,
					     isl_union_map_from_map(tiling));

	gen->shared_len = gen->tile_first + gen->tile_len + gen->n_grid;

	return sched;
}

static __isl_give isl_union_map *parametrize_tiled_schedule(
	struct cuda_gen *gen, __isl_take isl_union_map *sched)
{
	isl_dim *dim;
	isl_set *par;

	dim = isl_union_map_get_dim(sched);
	par = parametrization(isl_dim_copy(dim),
		gen->tiled_len, 0, gen->tile_first, "h");
	sched = isl_union_map_intersect_range(sched,
						isl_union_set_from_set(par));

	par = parametrization(dim, gen->tiled_len,
		gen->tile_first + gen->n_grid, gen->n_grid, "b");
	sched = isl_union_map_intersect_range(sched,
						isl_union_set_from_set(par));

	return sched;
}

/* Tile/wrap the P1 loops over the threads.
 */
static __isl_give isl_union_map *thread_tile_schedule(struct cuda_gen *gen,
	__isl_take isl_union_map *sched)
{
	isl_dim *dim;
	isl_map *tiling;
	isl_set *par;

	dim = isl_union_map_get_dim(sched);

	if (gen->options->cuda_wrap)
		tiling = wrap(isl_dim_copy(dim), gen->tiled_len,
				gen->shared_len, gen->n_block, gen->block_dim);
	else
		tiling = tile(isl_dim_copy(dim), gen->tiled_len,
				gen->shared_len, gen->n_block, gen->block_dim);
	gen->thread_tiled_len = gen->tiled_len + gen->n_block;

	sched = isl_union_map_apply_range(sched,
					     isl_union_map_from_map(tiling));

	par = parametrization(dim, gen->thread_tiled_len,
		gen->tile_first + gen->tile_len + gen->n_grid + gen->n_block,
		gen->n_block, "t");
	sched = isl_union_map_intersect_range(sched,
						isl_union_set_from_set(par));

	gen->shared_len = gen->tile_first + gen->tile_len + gen->n_grid;

	return sched;
}

/* If the user asked for it, scale the shared memory tile loops
 * (T1P and T2) of "sched" by gen->tile_size[i].
 * If we are not performing "wrapping", then additionally scale the T1P
 * loops by gen->grid_dim[i].
 */
static __isl_give isl_union_map *scale_tile_loops(struct cuda_gen *gen,
	__isl_take isl_union_map *sched)
{
	int i;
	isl_dim *dim;
	isl_basic_map *scale;
	isl_constraint *c;

	if (!gen->options->cuda_scale_tile_loops)
		return sched;

	dim = isl_union_map_get_dim(sched);
	dim = isl_dim_add(dim, isl_dim_in, gen->tiled_len);
	dim = isl_dim_add(dim, isl_dim_out, gen->tiled_len);
	scale = isl_basic_map_universe(isl_dim_copy(dim));

	for (i = 0; i < gen->tiled_len; ++i) {
		int f = 1;

		if (i >= gen->tile_first && i < gen->tile_first + gen->n_grid) {
			f = gen->tile_size[i - gen->tile_first];
			if (!gen->options->cuda_wrap)
				f *= gen->grid_dim[i - gen->tile_first];
		} else if (i >= gen->tile_first + gen->n_grid &&
			   i < gen->tile_first + gen->n_grid + gen->tile_len) {
			f = gen->tile_size[i - (gen->tile_first + gen->n_grid)];
		}

		c = isl_equality_alloc(isl_dim_copy(dim));
		isl_constraint_set_coefficient_si(c, isl_dim_in, i, f);
		isl_constraint_set_coefficient_si(c, isl_dim_out, i, -1);
		scale = isl_basic_map_add_constraint(scale, c);
	}

	isl_dim_free(dim);

	sched = isl_union_map_apply_range(sched,
		isl_union_map_from_map(isl_map_from_basic_map(scale)));

	return sched;
}

/* If we are not performing "wrapping" and if the user asked for it,
 * scale the thread tile loops (P1T) of "sched" by gen->block_dim[i].
 */
static __isl_give isl_union_map *scale_thread_tile_loops(struct cuda_gen *gen,
	__isl_take isl_union_map *sched)
{
	int i;
	isl_dim *dim;
	isl_basic_map *scale;
	isl_constraint *c;

	if (gen->options->cuda_wrap)
		return sched;
	if (!gen->options->cuda_scale_tile_loops)
		return sched;

	dim = isl_union_map_get_dim(sched);
	dim = isl_dim_add(dim, isl_dim_in, gen->thread_tiled_len);
	dim = isl_dim_add(dim, isl_dim_out, gen->thread_tiled_len);
	scale = isl_basic_map_universe(isl_dim_copy(dim));

	for (i = 0; i < gen->thread_tiled_len; ++i) {
		int f = 1;

		if (i >= gen->shared_len &&
		    i < gen->shared_len + gen->n_block)
			f = gen->block_dim[i - gen->shared_len];

		c = isl_equality_alloc(isl_dim_copy(dim));
		isl_constraint_set_coefficient_si(c, isl_dim_in, i, f);
		isl_constraint_set_coefficient_si(c, isl_dim_out, i, -1);
		scale = isl_basic_map_add_constraint(scale, c);
	}

	isl_dim_free(dim);

	sched = isl_union_map_apply_range(sched,
		isl_union_map_from_map(isl_map_from_basic_map(scale)));

	return sched;
}

/* If we are not performing "wrapping" and if the user asked for it,
 * scale the first "n_tile" loops of "sched" by gen->block_dim[i].
 */
static __isl_give isl_union_map *scale_access_tile_loops(struct cuda_gen *gen,
	__isl_take isl_union_map *sched, int len, int n_tile)
{
	int i;
	isl_dim *dim;
	isl_basic_map *scale;
	isl_constraint *c;

	if (gen->options->cuda_wrap)
		return sched;
	if (!gen->options->cuda_scale_tile_loops)
		return sched;

	dim = isl_union_map_get_dim(sched);
	dim = isl_dim_add(dim, isl_dim_in, len);
	dim = isl_dim_add(dim, isl_dim_out, len);
	scale = isl_basic_map_universe(isl_dim_copy(dim));

	for (i = 0; i < len; ++i) {
		int f = 1;

		if (i < n_tile)
			f = gen->block_dim[i];

		c = isl_equality_alloc(isl_dim_copy(dim));
		isl_constraint_set_coefficient_si(c, isl_dim_in, i, f);
		isl_constraint_set_coefficient_si(c, isl_dim_out, i, -1);
		scale = isl_basic_map_add_constraint(scale, c);
	}

	isl_dim_free(dim);

	sched = isl_union_map_apply_range(sched,
		isl_union_map_from_map(isl_map_from_basic_map(scale)));

	return sched;
}

static void print_cloog_shared_body(struct cuda_gen *gen,
	__isl_keep isl_set *context, __isl_keep isl_union_map *sched,
	int len)
{
	int i;
	CloogOptions *options;
	CloogDomain *cloog_context;
	CloogUnionDomain *ud;
	CloogInput *input;
	struct clast_stmt *stmt;
	char name[20];

	sched = isl_union_map_copy(sched);
	sched = isl_union_map_align_params(sched, isl_set_get_dim(context));

	options = cloog_options_malloc(gen->state);
	options->language = LANGUAGE_C;
	options->strides = 1;
	options->sh = 1;

	ud = cloog_union_domain_from_isl_union_map(sched);
	for (i = 0; i < len; ++i) {
		snprintf(name, sizeof(name), "c%d", i);
		ud = cloog_union_domain_set_name(ud, CLOOG_SCAT, i, name);
	}
	cloog_context = cloog_domain_from_isl_set(isl_set_copy(context));
	input = cloog_input_alloc(cloog_context, ud);

	stmt = cloog_clast_create_from_input(input, options);
	clast_pprint(gen->cuda.kernel_c, stmt,
			gen->kernel_code.indent, options);

	cloog_clast_free(stmt);
	cloog_options_free(options);
}

/* Add "len" parameters p[i] called prefix%d,
 * with bounds to 0 <= p[i] < size[i].
 */
__isl_give isl_set *add_bounded_parameters(__isl_take isl_set *set,
	int len, int *size, const char *prefix)
{
	int i;
	unsigned nparam;
	isl_int v;
	isl_dim *dim;
	isl_basic_set *bset;
	isl_constraint *c;
	char name[20];

	nparam = isl_set_dim(set, isl_dim_param);
	set = isl_set_add_dims(set, isl_dim_param, len);

	for (i = 0; i < len; ++i) {
		snprintf(name, sizeof(name), "%s%d", prefix, i);
		set = isl_set_set_dim_name(set, isl_dim_param,
					    nparam + i, name);
	}

	dim = isl_set_get_dim(set);
	bset = isl_basic_set_universe(isl_dim_copy(dim));

	isl_int_init(v);

	for (i = 0; i < len; ++i) {
		c = isl_inequality_alloc(isl_dim_copy(dim));
		isl_int_set_si(v, 1);
		isl_constraint_set_coefficient(c, isl_dim_param, nparam + i, v);
		bset = isl_basic_set_add_constraint(bset, c);
	
		c = isl_inequality_alloc(isl_dim_copy(dim));
		isl_int_set_si(v, -1);
		isl_constraint_set_coefficient(c, isl_dim_param, nparam + i, v);
		isl_int_set_si(v, size[i] - 1);
		isl_constraint_set_constant(c, v);
		bset = isl_basic_set_add_constraint(bset, c);
	}

	isl_int_clear(v);
	isl_dim_free(dim);

	return isl_set_intersect(set, isl_set_from_basic_set(bset));
}

static void print_shared_body(struct cuda_gen *gen,
	__isl_keep isl_set *shared_domain, __isl_keep isl_union_map *sched,
	int len)
{
	isl_set *context;

	context = isl_set_copy(shared_domain);
	context = parametrize(context, 0, gen->shared_len, "g");
	context = isl_set_project_out(context, isl_dim_set, 0, gen->shared_len);
	context = add_bounded_parameters(context,
					gen->n_block, gen->block_dim, "t");

	print_cloog_shared_body(gen, context, sched, len);

	isl_set_free(context);
}

/* Construct a schedule for iterating over all elements in the given
 * piece of an array.  We essentially build an identity mapping,
 * except that we move the iteration over the final array index
 * to the first dimension so that it gets tiled/wrapped over threadIdx.x,
 * which should improve coalescing.
 * We subsequently also perform the tiling/wrapping over the threads.
 */
static __isl_give isl_union_map *access_schedule(struct cuda_gen *gen,
	__isl_take isl_set *access)
{
	int i;
	isl_dim *dim;
	isl_basic_map *bmap;
	isl_map *sched;
	isl_union_map *usched;
	isl_map *tiling;
	isl_set *par;
	unsigned nvar = isl_set_dim(access, isl_dim_set);
	isl_constraint *c;
	int n_tile;

	dim = isl_set_get_dim(access);
	dim = isl_dim_from_domain(dim);
	dim = isl_dim_add(dim, isl_dim_out, nvar);
	bmap = isl_basic_map_universe(isl_dim_copy(dim));

	for (i = 0; i < nvar; ++i) {
		int j = (i + 1 == nvar) ? 0 : i + 1;

		c = isl_equality_alloc(isl_dim_copy(dim));
		isl_constraint_set_coefficient_si(c, isl_dim_in, i, 1);
		isl_constraint_set_coefficient_si(c, isl_dim_out, j, -1);
		bmap = isl_basic_map_add_constraint(bmap, c);
	}

	isl_dim_free(dim);
	sched = isl_map_from_basic_map(bmap);
	sched = isl_map_intersect_domain(sched, access);

	n_tile = gen->n_block;
	if (n_tile > nvar)
		n_tile = nvar;

	dim = isl_map_get_dim(sched);
	dim = isl_dim_drop(dim, isl_dim_in, 0, nvar);
	dim = isl_dim_drop(dim, isl_dim_out, 0, nvar);
	if (gen->options->cuda_wrap)
		tiling = wrap(isl_dim_copy(dim), nvar, 0,
				n_tile, gen->block_dim);
	else
		tiling = tile(isl_dim_copy(dim), nvar, 0,
				n_tile, gen->block_dim);
	sched = isl_map_apply_range(sched, tiling);

	par = parametrization(dim, nvar + n_tile, n_tile, n_tile, "t");
	usched = isl_union_map_from_map(sched);
	usched = isl_union_map_intersect_range(usched,
						isl_union_set_from_set(par));

	usched = scale_access_tile_loops(gen, usched, nvar + n_tile, n_tile);

	return usched;
}

static void print_shared_access(struct cuda_gen *gen,
	__isl_keep isl_set *shared_domain, __isl_take isl_set *access,
	const char *type)
{
	const char *array_name;
	char *name;
	isl_ctx *ctx;
	isl_union_map *sched;
	unsigned nvar = isl_set_dim(access, isl_dim_set);
	int n_tile;

	if (isl_set_fast_is_empty(access)) {
		isl_set_free(access);
		return;
	}

	ctx = isl_set_get_ctx(access);
	array_name = isl_set_get_tuple_name(access);
	name = isl_alloc_array(ctx, char,
		    strlen(type) + sizeof("_shared_") + strlen(array_name));
	sprintf(name, "%s_shared_%s", type, array_name);
	access = isl_set_set_tuple_name(access, name);
	free(name);

	sched = access_schedule(gen, access);

	n_tile = gen->n_block;
	if (n_tile > nvar)
		n_tile = nvar;

	print_shared_body(gen, shared_domain, sched, nvar + n_tile);

	isl_union_map_free(sched);
}

/* Print code for reading into or writing from shared memory.
 */
static void print_shared_accesses(struct cuda_gen *gen,
	__isl_keep isl_set *shared_domain, __isl_keep isl_union_map *access,
	const char *type)
{
	int i;
	isl_dim *dim;
	isl_map *proj;
	isl_set *par;
	isl_union_set *access_set;

	access = isl_union_map_copy(access);
	access = isl_union_map_apply_domain(access,
					isl_union_map_copy(gen->tiled_sched));
	dim = isl_union_map_get_dim(access);
	proj = projection(dim, gen->tiled_len, gen->shared_len);
	access = isl_union_map_apply_domain(access,
			isl_union_map_from_map(proj));
	access = isl_union_map_intersect_domain(access,
			isl_union_set_from_set(isl_set_copy(shared_domain)));

	dim = isl_union_map_get_dim(access);
	par = parametrization(dim, gen->shared_len, 0, gen->shared_len, "g");
	access = isl_union_map_intersect_domain(access,
						isl_union_set_from_set(par));

	access_set = isl_union_map_range(access);
	access_set = isl_union_set_coalesce(access_set);

	for (i = 0; i < gen->n_array; ++i) {
		isl_dim *dim;
		isl_set *set;

		dim = isl_union_set_get_dim(access_set);
		dim = isl_dim_add(dim, isl_dim_set, gen->array[i].n_index);
		dim = isl_dim_set_tuple_name(dim, isl_dim_set,
						gen->array[i].name);
		set = isl_union_set_extract_set(access_set, dim);

		print_shared_access(gen, shared_domain, set, type);
	}

	isl_union_set_free(access_set);
}

/* This function is called for each leaf in the clast of the kernel code.
 * We first specialize the schedule to the site of the leaf and
 * print code for reading into shared memory, performing the actual
 * computations and writing from shard memory, with the required
 * synchronizations.
 */
static void print_kernel_user(struct gpucode_info *code,
	struct clast_user_stmt *u)
{
	struct cuda_gen *gen = code->user;
	isl_dim *dim;
	isl_set *par;
	isl_set *shared_domain;
	isl_union_map *local_sched;

	shared_domain = extract_entire_host_domain(u);

	local_sched = isl_union_map_intersect_range(
		    isl_union_map_copy(gen->tiled_sched),
		    isl_union_set_from_set(extend(isl_set_copy(shared_domain),
						  gen->tiled_len)));
	dim = isl_union_map_get_dim(local_sched);
	par = parametrization(dim, gen->tiled_len, 0, gen->shared_len, "g");
	local_sched = isl_union_map_intersect_range(local_sched,
						isl_union_set_from_set(par));

	local_sched = thread_tile_schedule(gen, local_sched);
	local_sched = scale_thread_tile_loops(gen, local_sched);

	print_shared_accesses(gen, shared_domain, gen->prog->read, "read");

	print_indent(gen->cuda.kernel_c, gen->kernel_code.indent);
	fprintf(gen->cuda.kernel_c, "__syncthreads();\n");

	print_shared_body(gen, shared_domain, local_sched,
			    gen->thread_tiled_len);

	print_indent(gen->cuda.kernel_c, gen->kernel_code.indent);
	fprintf(gen->cuda.kernel_c, "__syncthreads();\n");

	print_shared_accesses(gen, shared_domain, gen->prog->write, "write");

	print_indent(gen->cuda.kernel_c, gen->kernel_code.indent);
	fprintf(gen->cuda.kernel_c, "__syncthreads();\n");

	isl_union_map_free(local_sched);
	isl_set_free(shared_domain);
}

/* Use CLooG to generate code for the outer gen->shared_first loops
 * of the local schedule "sched".
 * The pretty printing of this code is handled by gpu_print_host_stmt,
 * which calls print_kernel_user for each iteration of the shared tile loops.
 */
static void print_cloog_kernel_body(struct cuda_gen *gen,
	__isl_keep isl_set *context, __isl_keep isl_union_map *sched)
{
	int i;
	CloogOptions *options;
	CloogDomain *cloog_context;
	CloogUnionDomain *ud;
	CloogInput *input;
	struct clast_stmt *stmt;
	char name[20];

	sched = isl_union_map_copy(sched);
	sched = isl_union_map_align_params(sched, isl_set_get_dim(context));

	options = cloog_options_malloc(gen->state);
	options->language = LANGUAGE_C;
	options->strides = 1;
	options->sh = 1;
	options->stop = gen->shared_len;
	options->f = gen->tiled_len;
	options->l = gen->tiled_len;
	options->save_domains = 1;
	options->noscalars = 1;

	ud = cloog_union_domain_from_isl_union_map(sched);
	for (i = 0; i < gen->shared_len; ++i) {
		snprintf(name, sizeof(name), "g%d", i);
		ud = cloog_union_domain_set_name(ud, CLOOG_SCAT, i, name);
	}
	cloog_context = cloog_domain_from_isl_set(isl_set_copy(context));
	input = cloog_input_alloc(cloog_context, ud);

	stmt = cloog_clast_create_from_input(input, options);

	gen->kernel_code.indent = 4;
	gen->kernel_code.dst = gen->cuda.kernel_c;
	gen->kernel_code.print_user_stmt = &print_kernel_user;
	gen->kernel_code.user = gen;
	gpu_print_host_stmt(&gen->kernel_code, stmt);

	cloog_clast_free(stmt);
	cloog_options_free(options);
}

static void print_kernel_iterators(struct cuda_gen *gen)
{
	int i;
	const char *block_dims[] = { "blockIdx.x", "blockIdx.y" };
	const char *thread_dims[] = { "threadIdx.x", "threadIdx.y",
					"threadIdx.z" };

	print_indent(gen->cuda.kernel_c, 4);
	fprintf(gen->cuda.kernel_c, "int ");
	for (i = 0; i < gen->tiled_len + gen->n_block; ++i) {
		if (i)
			fprintf(gen->cuda.kernel_c, ", ");
		if (i < gen->shared_len)
			fprintf(gen->cuda.kernel_c, "g%d", i);
		else
			fprintf(gen->cuda.kernel_c, "c%d", i);
	}
	fprintf(gen->cuda.kernel_c, ";\n");

	if (gen->n_grid > 0) {
		print_indent(gen->cuda.kernel_c, 4);
		fprintf(gen->cuda.kernel_c, "int ");
		for (i = 0; i < gen->n_grid; ++i) {
			if (i)
				fprintf(gen->cuda.kernel_c, ", ");
			fprintf(gen->cuda.kernel_c, "b%d = %s",
							i, block_dims[i]);
		}
		fprintf(gen->cuda.kernel_c, ";\n");
	}

	if (gen->n_block > 0) {
		print_indent(gen->cuda.kernel_c, 4);
		fprintf(gen->cuda.kernel_c, "int ");
		for (i = 0; i < gen->n_block; ++i) {
			if (i)
				fprintf(gen->cuda.kernel_c, ", ");
			fprintf(gen->cuda.kernel_c, "t%d = %s",
							i, thread_dims[i]);
		}
		fprintf(gen->cuda.kernel_c, ";\n");
	}
}

static void print_shared_arrays(struct cuda_gen *gen)
{
	int i, j;

	for (i = 0; i < gen->n_array; ++i) {
		if (!gen->array[i].active)
			continue;
		print_indent(gen->cuda.kernel_c, 4);
		fprintf(gen->cuda.kernel_c, "__shared__ %s shared_%s",
			gen->options->type, gen->array[i].name);
		for (j = 0; j < gen->array[i].n_index; ++j) {
			fprintf(gen->cuda.kernel_c, "[");
			isl_int_print(gen->cuda.kernel_c,
				gen->array[i].shared_bound[j].size, 0);
			fprintf(gen->cuda.kernel_c, "]");
		}
		fprintf(gen->cuda.kernel_c, ";\n");
	}
}

static void print_kernel_body(struct cuda_gen *gen,
	__isl_keep isl_set *host_domain, __isl_keep isl_union_map *sched)
{
	isl_set *context;

	context = isl_set_copy(host_domain);
	context = parametrize(context, 0, gen->tile_first, "h");
	context = isl_set_project_out(context, isl_dim_set, 0, gen->tile_first);
	context = add_bounded_parameters(context,
					gen->n_grid, gen->grid_dim, "b");

	print_kernel_iterators(gen);
	print_shared_arrays(gen);

	fprintf(gen->cuda.kernel_c, "\n");

	print_cloog_kernel_body(gen, context, sched);

	isl_set_free(context);
}

struct cuda_size_info {
	isl_basic_set *bset;
	struct cuda_array_bound *bound;
	int pos;
};

/* Given a constraint from the basic set describing the bounds on
 * an array index, check if it is a lower bound, say i >= b(x), and,
 * if so, check whether the expression "i - b(x) + 1" has a constant
 * upper bound.  If so, and if this bound is smaller than any bound
 * derived from earlier constraints, set the size to this bound on
 * the expression and the lower bound to b(x).
 * Note that we currently only consider lower bounds in which
 * the coefficient of the array index is one.  In the general case,
 * m i >= b(x), we would have to work with ceil(b(x)/m).
 */
static int compute_size_in_direction(__isl_take isl_constraint *c, void *user)
{
	struct cuda_size_info *size = user;
	unsigned nparam;
	isl_int v;

	nparam = isl_basic_set_dim(size->bset, isl_dim_param);

	isl_int_init(v);

	isl_constraint_get_coefficient(c, isl_dim_set, size->pos, &v);

	if (isl_int_is_one(v)) {
		isl_int d;
		isl_dim *dim;
		isl_set *set;
		isl_qpolynomial *lb, *qp, *var, *one;
		isl_pw_qpolynomial *pwqp;
		isl_pw_qpolynomial_fold *pwf;

		lb = isl_qpolynomial_from_constraint(c, isl_dim_set, size->pos);
		dim = isl_qpolynomial_get_dim(lb);
		one = isl_qpolynomial_one(isl_dim_copy(dim));
		var = isl_qpolynomial_var(dim, isl_dim_set, size->pos);
		qp = isl_qpolynomial_sub(var, isl_qpolynomial_copy(lb));
		qp = isl_qpolynomial_add(qp, one);
		set = isl_set_from_basic_set(isl_basic_set_copy(size->bset));
		pwqp = isl_pw_qpolynomial_alloc(set, qp);
		pwqp = isl_pw_qpolynomial_move_dims(pwqp, isl_dim_set, 0,
						    isl_dim_param, 0, nparam);
		pwf = isl_pw_qpolynomial_bound(pwqp, isl_fold_max, NULL);
		dim = isl_pw_qpolynomial_fold_get_dim(pwf);
		qp = isl_pw_qpolynomial_fold_eval(pwf, isl_point_zero(dim));
		if (!isl_qpolynomial_is_infty(qp)) {
			isl_int_init(d);
			isl_qpolynomial_is_cst(qp, &v, &d);
			isl_int_fdiv_q(v, v, d);
			if (isl_int_is_neg(size->bound->size) ||
			    isl_int_lt(v, size->bound->size)) {
				isl_int_set(size->bound->size, v);
				isl_qpolynomial_free(size->bound->lb);
				size->bound->lb = isl_qpolynomial_copy(lb);
			}
			isl_int_clear(d);
		}
		isl_qpolynomial_free(lb);
		isl_qpolynomial_free(qp);
	} else
		isl_constraint_free(c);

	isl_int_clear(v);

	return 0;
}

/* Given a basic map "bounds" that maps parameters and input dimensions
 * to a single output dimension, look for an expression in the parameters
 * and input dimensions such that the range of the output dimension shifted
 * by this expression is a constant.
 *
 * In particular, we currently only consider lower bounds on the output
 * dimension as candidate expressions.
 */
static void compute_array_dim_size(struct cuda_gen *gen,
	struct cuda_array_bound *bound, __isl_take isl_basic_map *bounds)
{
	struct cuda_size_info size;

	isl_int_set_si(bound->size, -1);
	bound->lb = NULL;

	size.bound = bound;
	size.pos = isl_basic_map_dim(bounds, isl_dim_in);
	size.bset = isl_basic_map_wrap(bounds);
	size.bset = isl_basic_set_flatten(size.bset);
	isl_basic_set_foreach_constraint(size.bset, &compute_size_in_direction,
					&size);
	isl_basic_set_free(size.bset);

	assert(isl_int_is_nonneg(bound->size));
}

/* Compute the size of the shared array corresonding to the given array,
 * based on the given accesses from the current kernel,
 * as well as the offset of the shared piece in the original array.
 *
 * We project the accesses on each index in turn and look for a parametric
 * offset such that the size is constant.
 */
static void compute_array_shared_size(struct cuda_gen *gen,
	struct cuda_array_info *array, __isl_take isl_map *access)
{
	int i;

	if (isl_map_fast_is_empty(access)) {
		isl_map_free(access);
		array->active = 0;
		return;
	}

	array->active = 1;

	for (i = 0; i < array->n_index; ++i) {
		isl_map *access_i;
		isl_basic_map *bounds;

		access_i = isl_map_copy(access);
		access_i = isl_map_project_out(access_i, isl_dim_out, 0, i);
		access_i = isl_map_project_out(access_i, isl_dim_out,
					    i + 1, array->n_index - (i + 1));
		bounds = isl_map_simple_hull(access_i);
		compute_array_dim_size(gen, &array->shared_bound[i], bounds);
	}

	isl_map_free(access);
}

/* Compute the sizes of all shared arrays for the current kernel,
 * as well as the offsets of the shared pieces in the original arrays.
 */
static void compute_shared_size(struct cuda_gen *gen)
{
	int i;
	isl_dim *dim;
	isl_map *proj;
	isl_set *par;
	isl_union_map *access;

	access = isl_union_map_union(isl_union_map_copy(gen->prog->read),
				     isl_union_map_copy(gen->prog->write));
	access = isl_union_map_apply_domain(access,
					isl_union_map_copy(gen->tiled_sched));

	dim = isl_union_map_get_dim(access);
	proj = projection(dim, gen->tiled_len, gen->shared_len);
	access = isl_union_map_apply_domain(access,
					isl_union_map_from_map(proj));

	dim = isl_union_map_get_dim(access);
	par = parametrization(dim, gen->shared_len, 0, gen->shared_len, "g");
	access = isl_union_map_intersect_domain(access,
						isl_union_set_from_set(par));

	for (i = 0; i < gen->n_array; ++i) {
		isl_dim *dim;
		isl_map *acc;

		dim = isl_union_map_get_dim(access);
		dim = isl_dim_add(dim, isl_dim_out, gen->array[i].n_index);
		dim = isl_dim_set_tuple_name(dim, isl_dim_out,
						gen->array[i].name);
		dim = isl_dim_add(dim, isl_dim_in, gen->shared_len);
		acc = isl_union_map_extract_map(access, dim);

		compute_array_shared_size(gen, &gen->array[i], acc);
	}

	isl_union_map_free(access);
}

/* Free all array information that is local to the current kernel.
 */
static void free_local_array_info(struct cuda_gen *gen)
{
	int i, j;

	for (i = 0; i < gen->n_array; ++i) {
		if (!gen->array[i].active)
			continue;
		for (j = 0; j < gen->array[i].n_index; ++j) {
			isl_qpolynomial_free(gen->array[i].shared_bound[j].lb);
			gen->array[i].shared_bound[j].lb = NULL;
			isl_pw_qpolynomial_fold_free(gen->array[i].local_bound[j]);
			gen->array[i].local_bound[j] = NULL;
		}
	}
}

static void print_iterator_list(FILE *out, int len, const char *prefix,
	int parens)
{
	int i;

	fprintf(out, "(");
	for (i = 0; i < len; ++i) {
		if (i)
			fprintf(out, ", ");
		if (parens)
			fprintf(out, "(%s%d)", prefix, i);
		else
			fprintf(out, "%s%d", prefix, i);
	}
	fprintf(out, ")");
}

/* Print an access to the element in the global memory copy of the
 * given array that corresponds to element [a0][a1]... of the original array.
 * The copy in global memory has been linearized, so we need to take
 * the array size into acount.
 */
static void print_global_index(isl_ctx *ctx, FILE *out,
	struct cuda_array_info *array)
{
	int i;
	isl_printer *prn;

	fprintf(out, "%s[", array->name);
	for (i = 0; i + 1 < array->n_index; ++i)
		fprintf(out, "(");
	for (i = 0; i < array->n_index; ++i) {
		if (i) {
			prn = isl_printer_to_file(ctx, out);
			prn = isl_printer_set_output_format(prn, ISL_FORMAT_C);
			prn = isl_printer_print_str(prn, ") * (");
			prn = isl_printer_print_pw_qpolynomial_fold(prn,
							array->local_bound[i]);
			prn = isl_printer_print_str(prn, ") + ");
			isl_printer_free(prn);
		}
		fprintf(out, "a%d", i);
	}
	fprintf(out, "]");
}

/* Print an access to the element in the shared memory copy of the
 * given array that corresponds to element [a0][a1]... of the original array.
 * Since the array in shared memory is just a shifted copy of part
 * of the original array, we simply need to subtract the lower bound,
 * which was computed in compute_array_shared_size.
 */
static void print_local_index(FILE *out, struct cuda_array_info *array)
{
	int i;
	isl_ctx *ctx;
	isl_printer *prn;

	ctx = isl_dim_get_ctx(array->dim);
	fprintf(out, "shared_%s", array->name);
	for (i = 0; i < array->n_index; ++i) {
		fprintf(out, "[a%d - (", i);
		prn = isl_printer_to_file(ctx, out);
		prn = isl_printer_set_output_format(prn, ISL_FORMAT_C);
		prn = isl_printer_print_qpolynomial(prn,
						array->shared_bound[i].lb);
		isl_printer_free(prn);
		fprintf(out, ")]");
	}
}

/* Print '#define's for copying data from global memory to shared
 * memory and back for the given array.
 */
static void print_array_copy_defines(struct cuda_gen *gen,
	struct cuda_array_info *array)
{
	int i;
	const char *type[] = { "read", "write" };

	for (i = 0; i < 2; ++i) {
		fprintf(gen->cuda.kernel_c, "#define %s_shared_%s",
			type[i], array->name);
		print_iterator_list(gen->cuda.kernel_c, array->n_index, "a", 0);
		fprintf(gen->cuda.kernel_c, " %s_shared_%s_",
			type[i], array->name);
		print_iterator_list(gen->cuda.kernel_c, array->n_index, "a", 1);
		fprintf(gen->cuda.kernel_c, "\n");

		fprintf(gen->cuda.kernel_c, "#define %s_shared_%s_",
			type[i], array->name);
		print_iterator_list(gen->cuda.kernel_c, array->n_index, "a", 0);
		if (i) {
			fprintf(gen->cuda.kernel_c, " ");
			print_global_index(gen->ctx, gen->cuda.kernel_c, array);
			fprintf(gen->cuda.kernel_c, " = ");
			print_local_index(gen->cuda.kernel_c, array);
		} else {
			fprintf(gen->cuda.kernel_c, " ");
			print_local_index(gen->cuda.kernel_c, array);
			fprintf(gen->cuda.kernel_c, " = ");
			print_global_index(gen->ctx, gen->cuda.kernel_c, array);
		}
		fprintf(gen->cuda.kernel_c, "\n");
	}
}

static void print_copy_defines(struct cuda_gen *gen)
{
	int i;

	for (i = 0; i < gen->n_array; ++i) {
		if (!gen->array[i].active)
			continue;
		print_array_copy_defines(gen, &gen->array[i]);
	}
}

/* This function is called for each access to an array in some statement
 * in the original code.
 * Replace that access by an access to shared memory.
 * Since the array in shared memory is just
 * a shifted copy of part of the original array, we simply need
 * to subtract the lower bound, which was computed
 * in compute_array_shared_size.
 */
static void print_access(__isl_take isl_map *access, Stmt *stmt, void *user)
{
	struct cuda_gen *gen = user;
	int i;
	const char *name;
	unsigned n_index;
	struct cuda_array_info *array = NULL;
	isl_printer *prn;
	isl_basic_map *aff;

	name = isl_map_get_tuple_name(access, isl_dim_out);
	fprintf(gen->cuda.kernel_c, "shared_%s", name);

	for (i = 0; i < gen->n_array; ++i) {
		if (strcmp(name, gen->array[i].name))
			continue;
		array = &gen->array[i];
	}
	assert(array);

	n_index = isl_map_dim(access, isl_dim_out);
	aff = isl_map_affine_hull(access);

	for (i = 0; i < n_index; ++i) {
		isl_constraint *c;
		isl_qpolynomial *qp;
		int ok;

		ok = isl_basic_map_has_defining_equality(aff,
							isl_dim_out, i, &c);
		assert(ok);
		qp = isl_qpolynomial_from_constraint(c, isl_dim_out, i);

		fprintf(gen->cuda.kernel_c, "[(");
		prn = isl_printer_to_file(gen->ctx, gen->cuda.kernel_c);
		prn = isl_printer_print_qpolynomial(prn, qp);
		prn = isl_printer_print_str(prn, ") - (");
		prn = isl_printer_set_output_format(prn, ISL_FORMAT_C);
		prn = isl_printer_print_qpolynomial(prn,
						array->shared_bound[i].lb);
		isl_printer_free(prn);
		fprintf(gen->cuda.kernel_c, ")]");
		isl_qpolynomial_free(qp);
	}

	isl_basic_map_free(aff);
}

static void print_name_list(FILE *out, int len, char **list)
{
	int i;

	fprintf(out, "(");
	for (i = 0; i < len; ++i) {
		if (i)
			fprintf(out, ", ");
		fprintf(out, "%s", list[i]);
	}
	fprintf(out, ")");
}

static void print_statement_defines(struct cuda_gen *gen)
{
	int i;
	char name[20];
	isl_union_set *domain;

	domain = isl_union_map_domain(isl_union_map_copy(gen->tiled_sched));

	for (i = 0; i < gen->prog->nstmts; ++i) {
		isl_dim *dim;
		isl_set *domain_i;
		int empty;
		Stmt *stmt = &gen->prog->stmts[i];

		snprintf(name, sizeof(name), "S_%d", i);

		dim = isl_union_set_get_dim(domain);
		dim = isl_dim_add(dim, isl_dim_set, stmt->dim);
		dim = isl_dim_set_tuple_name(dim, isl_dim_set, name);
		domain_i = isl_union_set_extract_set(domain, dim);
		empty = isl_set_fast_is_empty(domain_i);
		isl_set_free(domain_i);
		if (empty)
			continue;

		fprintf(gen->cuda.kernel_c, "#define S_%d", i);
		print_iterator_list(gen->cuda.kernel_c, stmt->dim, "i", 0);
		fprintf(gen->cuda.kernel_c, " S_%d_", i);
		print_iterator_list(gen->cuda.kernel_c, stmt->dim, "i", 1);
		fprintf(gen->cuda.kernel_c, "\n");

		fprintf(gen->cuda.kernel_c, "#define S_%d_", i);
		print_name_list(gen->cuda.kernel_c, stmt->dim, stmt->iterators);
		fprintf(gen->cuda.kernel_c, " ");
		print_stmt_body(gen->cuda.kernel_c, gen->ctx, gen->prog,
				stmt, name, &print_access, gen);
		fprintf(gen->cuda.kernel_c, "\n");
	}
	isl_union_set_free(domain);
}

/* The sizes of the arrays on the host that have been computed by
 * extract_array_info may depend on the parameters.  Use the extra
 * constraints on the parameters that are valid at "host_domain"
 * to simplify these expressions.
 */
static void localize_bounds(struct cuda_gen *gen,
	__isl_keep isl_set *host_domain)
{
	int i, j;
	isl_set *context;
	unsigned nvar;

	context = isl_set_copy(host_domain);
	nvar = isl_set_dim(host_domain, isl_dim_set);
	context = isl_set_project_out(host_domain, isl_dim_set, 0, nvar);

	for (i = 0; i < gen->n_array; ++i) {
		struct cuda_array_info *array = &gen->array[i];

		if (!array->active)
			continue;

		for (j = 0; j < array->n_index; ++j) {
			isl_pw_qpolynomial_fold *pwf;

			pwf = isl_pw_qpolynomial_fold_copy(array->bound[j]);
			pwf = isl_pw_qpolynomial_fold_gist(pwf,
							isl_set_copy(context));
			array->local_bound[j] = pwf;
		}
	}
	isl_set_free(context);
}

/* This function is called for each leaf in the clast of the host code.
 * We first specialize the schedule to the site of the leaf, compute
 * the size of shared memory and then print the body of host code
 * and the associated kernel (through a call to print_kernel_body).
 */
static void print_host_user(struct gpucode_info *code,
	struct clast_user_stmt *u)
{
	struct cuda_gen *gen = code->user;
	isl_set *host_domain;
	isl_union_map *access;
	isl_union_map *local_sched;
	isl_union_set *arrays;

	host_domain = extract_entire_host_domain(u);

	local_sched = isl_union_map_intersect_range(
		    isl_union_map_copy(gen->sched),
		    isl_union_set_from_set(extend(isl_set_copy(host_domain),
						  gen->untiled_len)));
	access = isl_union_map_union(isl_union_map_copy(gen->prog->read),
				     isl_union_map_copy(gen->prog->write));
	access = isl_union_map_apply_domain(access,
					    isl_union_map_copy(local_sched));
	arrays = isl_union_map_range(access);

	print_indent(code->dst, code->indent);
	fprintf(code->dst, "dim3 k%d_dimBlock(", gen->kernel_id);
	print_list(code->dst, gen->n_block, gen->block_dim);
	fprintf(code->dst, ");\n");

	print_indent(code->dst, code->indent);
	fprintf(code->dst, "dim3 k%d_dimGrid(", gen->kernel_id);
	print_list(code->dst, gen->n_grid, gen->grid_dim);
	fprintf(code->dst, ");\n");

	gen->tiled_sched = tile_schedule(gen, local_sched);
	gen->tiled_sched = parametrize_tiled_schedule(gen, gen->tiled_sched);
	gen->tiled_sched = scale_tile_loops(gen, gen->tiled_sched);

	compute_shared_size(gen);
	localize_bounds(gen, host_domain);

	print_copy_defines(gen);
	print_statement_defines(gen);
	print_kernel_launch(gen, arrays);

	fprintf(gen->cuda.kernel_c, "{\n");

	print_kernel_body(gen, host_domain, gen->tiled_sched);

	fprintf(gen->cuda.kernel_c, "}\n");

	free_local_array_info(gen);
	isl_union_map_free(gen->tiled_sched);
	isl_union_set_free(arrays);
	isl_set_free(host_domain);

	gen->kernel_id++;
}

/* Use CLooG to generate code for the outer gen->tile_first loops
 * of the global schedule in gen->sched.
 * The pretty printing of this code is handled by gpu_print_host_stmt,
 * which calls print_host_user for each kernel invocation location.
 */
static void print_cloog_host_code(struct cuda_gen *gen)
{
	int i;
	isl_set *context;
	isl_union_map *sched;
	CloogOptions *options;
	CloogDomain *cloog_context;
	CloogUnionDomain *ud;
	CloogInput *input;
	struct clast_stmt *stmt;
	char name[20];

	options = cloog_options_malloc(gen->state);
	options->language = LANGUAGE_C;
	options->otl = 0;
	options->strides = 1;
	options->stop = gen->tile_first;
	options->f = gen->untiled_len;
	options->l = gen->untiled_len;
	options->save_domains = 1;
	options->noscalars = 1;

	sched = isl_union_map_copy(gen->sched);
	ud = cloog_union_domain_from_isl_union_map(sched);
	for (i = 0; i < options->stop; ++i) {
		snprintf(name, sizeof(name), "h%d", i);
		ud = cloog_union_domain_set_name(ud, CLOOG_SCAT, i, name);
	}
	context = isl_set_copy(gen->prog->context);
	cloog_context = cloog_domain_from_isl_set(context);
	input = cloog_input_alloc(cloog_context, ud);

	stmt = cloog_clast_create_from_input(input, options);

	gen->code.indent = 0;
	gen->code.dst = gen->cuda.host_c;
	gen->code.print_user_stmt = &print_host_user;
	gen->code.user = gen;
	gpu_print_host_stmt(&gen->code, stmt);

	cloog_clast_free(stmt);
	cloog_options_free(options);
}

static void print_host_code(struct cuda_gen *gen)
{
	fprintf(gen->cuda.host_c, "{\n");
	print_cloog_macros(gen->cuda.host_c);
	print_cloog_macros(gen->cuda.kernel_c);

	declare_device_arrays(gen);
	declare_host_iterators(gen);

	allocate_device_arrays(gen);
	copy_arrays_to_device(gen);

	gen->kernel_id = 0;
	print_cloog_host_code(gen);

	copy_arrays_from_device(gen);

	fprintf(gen->cuda.host_c, "}\n");
}

/* Replace the scop in the "input" file by equivalent code
 * that uses the GPU.  "prog" is assumed to correspond to this scop
 * and the statements are assumed to have been assigned a valid "trans"
 * by Pluto that ensures that all dependences are non-negative in the
 * coordinate directions.
 *
 * We first select the outermost band of tilable dimensions.
 * If the first of these if not parallel, we apply a wavefront
 * transformation to make sure that the remaining loops are parallel.
 * In this case, the dimension corresponding to the wavefront is no longer
 * considered to be part of the tilable band.
 * We then have three block of dimension
 *
 *	H		B			G
 *
 * The tilable band "B" is first tiled according to "tile.sizes", resulting
 * in
 *
 *	H	T		P		G
 *
 * For each iteration of the T loop and for each array, we compute
 * the array elements accessed by that iteration, construct a rectangular
 * box around it and shift it to the origin.  The result is used
 * as shared memory for the array.
 *
 * We then split off at most 2 parallel loops from the T loops and
 * at most 3 parallel loops from the P loops
 *
 *	H	T1	T2	P1	P2	G
 *
 * The T1/P1 loops are then tiled or "wrapped" over the blocks/threads,
 * according to "grid.sizes"/"block.sizes".
 *
 *	H	T1P T1P	T2	P1T P1P	P2	G
 *
 * Finally, the T1P and P1P iterators are equated to the block and
 * thread dimensions respectively and so are effectively removed.
 * The H loops are run on the host.  The T1P, T2, P1T, P2 and G loops
 * are run on the GPU.
 *
 * Code is generated in three stages.  We first generate code for the
 * host (the H loops), with iterators h%d.  Then, for each leaf node
 * of the resulting AST, we generate code for the shared loops (up to
 * and including T2), with iterators g%d and after equating the H loops
 * to h%d parameters and the T1P loops to the block dimensions.
 * Finally, we generate code for the remaining loops in a similar fashion.
 *
 * The function frees "prog".
 */
int cuda(PlutoProg *prog, PlutoOptions *options, const char *input)
{
	struct cuda_gen gen;

	gen.ctx = prog->ctx;
	gen.prog = prog;
	gen.options = options;
	gen.state = cloog_isl_state_malloc(prog->ctx);

	cuda_open_files(&gen.cuda, input);

	collect_array_info(&gen, prog);

	select_tile_dimensions(&gen, prog);
	compute_global_schedule(&gen, prog, options);
	read_sizes(&gen);

	print_host_code(&gen);

	cloog_state_free(gen.state);
	clear_cuda_gen(&gen);
	pluto_prog_free(prog);

	cuda_close_files(&gen.cuda);

	return 0;
}
