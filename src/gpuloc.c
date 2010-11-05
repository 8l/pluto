#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#include <isl_constraint.h>
#include <isl_union_map.h>
#include <isl_union_set.h>
#include <cloog/isl/cloog.h>

#include "gpuloc.h"
#include "gpucode.h"
#include "program.h"

static void print_cloog_macros(FILE *dst)
{
    fprintf(dst,
        "#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))\n");
    fprintf(dst,
        "#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))\n");
    fprintf(dst, "#define max(x,y)    ((x) > (y) ? (x) : (y))\n");
    fprintf(dst, "#define min(x,y)    ((x) < (y) ? (x) : (y))\n");
}

/* Construct equality constraints equating each of the "n"
 * variables the affine expression in the corresponding row of "matrix".
 */
static __isl_give isl_mat *extract_equalities(isl_ctx *ctx,
        PlutoMatrix *matrix, int first, int n, int npar)
{
    int i, j;
    int n_col;
    isl_int v;
    isl_mat *eq;

    n_col = matrix->ncols;

    isl_int_init(v);
    eq = isl_mat_alloc(ctx, n, n + n_col + npar);

    for (i = 0; i < n; ++i) {
        isl_int_set_si(v, 0);
        for (j = 0; j < n; ++j)
            eq = isl_mat_set_element(eq, i, j, v);
        isl_int_set_si(v, -1);
        eq = isl_mat_set_element(eq, i, i, v);
        for (j = 0; j < n_col; ++j) {
            isl_int_set_si(v, matrix->val[first + i][j]);
            eq = isl_mat_set_element(eq, i, n + j, v);
        }
        isl_int_set_si(v, 0);
        for (j = 0; j < npar; ++j)
            eq = isl_mat_set_element(eq, i, n + n_col + j, v);
    }

    isl_int_clear(v);

    return eq;
}

/* Construct an isl_map performing the same operation
 * as multiplication by the given matrix.
 */
static __isl_give isl_map *pluto_matrix_to_isl_map(PlutoMatrix *matrix,
        __isl_take isl_dim *dim)
{
    int n_row, n_col, npar;
    isl_ctx *ctx;
    isl_mat *eq, *ineq;
    isl_basic_map *bmap;

    ctx = isl_dim_get_ctx(dim);
    n_row = matrix->nrows;
    n_col = matrix->ncols;

    npar = isl_dim_size(dim, isl_dim_param);
    ineq = isl_mat_alloc(ctx, 0, n_row + n_col + npar);
    eq = extract_equalities(ctx, matrix, 0, n_row, npar);

    bmap = isl_basic_map_from_constraint_matrices(dim, eq, ineq,
            isl_dim_out, isl_dim_in, isl_dim_div, isl_dim_cst, isl_dim_param);
    return isl_map_from_basic_map(bmap);
}

/* Extract the schedule mapping each iteration domain to a common
 * iteration space from the "isl_domain" and "trans" of each
 * statement in "prog".
 */
static __isl_give isl_union_map *extract_schedule(PlutoProg *prog)
{
    int i;
    isl_dim *dim;
    isl_union_map *schedule;

    dim = isl_dim_set_alloc(prog->ctx, prog->npar, 0);
    dim = set_dim_names(dim, isl_dim_param, prog->params);
    schedule = isl_union_map_empty(dim);

    for (i = 0; i <prog->nstmts; ++i) {
        Stmt *stmt = &prog->stmts[i];
        isl_map *schedule_i;

        dim = isl_dim_from_domain(isl_set_get_dim(stmt->isl_domain));
        dim = isl_dim_add(dim, isl_dim_out, stmt->trans->nrows);
        schedule_i = pluto_matrix_to_isl_map(stmt->trans, dim);
        schedule_i = isl_map_intersect_domain(schedule_i,
                                        isl_set_copy(stmt->isl_domain));

        schedule = isl_union_map_union(schedule,
                                        isl_union_map_from_map(schedule_i));
    }

    return schedule;
}

/* Construct a map from a domain of dimensionality "len"
 * to a domain of dimensionality "len" + "tile_len" that tiles
 * the "tile_len" coordinates starting at "first".
 * "dim" prescribes the parameters.
 */
static __isl_give isl_map *tile(__isl_take isl_dim *dim, int len,
        int first, int tile_len, int tile_size)
{
    int i;
    isl_int v;
    isl_basic_map *bmap;
    isl_constraint *c;

    isl_int_init(v);

    dim = isl_dim_add(dim, isl_dim_in, len);
    dim = isl_dim_add(dim, isl_dim_out, len + tile_len);
    bmap = isl_basic_map_universe(isl_dim_copy(dim));

    for (i = 0; i < len; ++i) {
        int j = i < first ? i : i + tile_len;

        c = isl_equality_alloc(isl_dim_copy(dim));
        isl_int_set_si(v, -1);
        isl_constraint_set_coefficient(c, isl_dim_in, i, v);
        isl_int_set_si(v, 1);
        isl_constraint_set_coefficient(c, isl_dim_out, j, v);
        bmap = isl_basic_map_add_constraint(bmap, c);
    }

    for (i = 0; i < tile_len; ++i) {
        c = isl_inequality_alloc(isl_dim_copy(dim));
        isl_int_set_si(v, 1);
        isl_constraint_set_coefficient(c, isl_dim_in, first + i, v);
        isl_int_set_si(v, -tile_size);
        isl_constraint_set_coefficient(c, isl_dim_out, first + i, v);
        bmap = isl_basic_map_add_constraint(bmap, c);

        c = isl_inequality_alloc(isl_dim_copy(dim));
        isl_int_set_si(v, -1);
        isl_constraint_set_coefficient(c, isl_dim_in, first + i, v);
        isl_int_set_si(v, tile_size);
        isl_constraint_set_coefficient(c, isl_dim_out, first + i, v);
        isl_int_set_si(v, tile_size - 1);
        isl_constraint_set_constant(c, v);
        bmap = isl_basic_map_add_constraint(bmap, c);
    }

    isl_dim_free(dim);
    isl_int_clear(v);

    return isl_map_from_basic_map(bmap);
}

/* Construct a map from a len-dimensional domain to
 * a (len-n)-dimensional domain that projects out the n coordinates
 * starting at first.
 * "dim" prescribes the parameters.
 */
static __isl_give isl_map *project_out(__isl_take isl_dim *dim,
    int len, int first, int n)
{
    int i, j;
    isl_constraint *c;
    isl_basic_map *bmap;
    isl_int v;

    isl_int_init(v);

    dim = isl_dim_add(dim, isl_dim_in, len);
    dim = isl_dim_add(dim, isl_dim_out, len - n);
    bmap = isl_basic_map_universe(isl_dim_copy(dim));

    for (i = 0, j = 0; i < len; ++i) {
        if (i >= first && i < first + n)
            continue;
        c = isl_equality_alloc(isl_dim_copy(dim));
        isl_int_set_si(v, -1);
        isl_constraint_set_coefficient(c, isl_dim_in, i, v);
        isl_int_set_si(v, 1);
        isl_constraint_set_coefficient(c, isl_dim_out, j, v);
        bmap = isl_basic_map_add_constraint(bmap, c);
        ++j;
    }
    isl_dim_free(dim);

    isl_int_clear(v);

    return isl_map_from_basic_map(bmap);
}

/* Construct a projection that maps a src_len dimensional domain
 * to its first dst_len coordinates.
 * "dim" prescribes the parameters.
 */
static __isl_give isl_map *projection(__isl_take isl_dim *dim,
    int src_len, int dst_len)
{
    return project_out(dim, src_len, dst_len, src_len - dst_len);
}

/* Construct a map that maps a domain of dimensionality "len"
 * to another domain of the same dimensionality such that
 * coordinate "first" of the range is equal to the sum of the "wave_len"
 * coordinates starting at "first" in the domain.
 * The remaining coordinates in the range are equal to the corresponding ones
 * in the domain.
 * "dim" prescribes the parameters.
 */
static __isl_give isl_map *wavefront(__isl_take isl_dim *dim, int len,
        int first, int wave_len)
{
    int i;
    isl_int v;
    isl_basic_map *bmap;
    isl_constraint *c;

    isl_int_init(v);

    dim = isl_dim_add(dim, isl_dim_in, len);
    dim = isl_dim_add(dim, isl_dim_out, len);
    bmap = isl_basic_map_universe(isl_dim_copy(dim));

    for (i = 0; i < len; ++i) {
        if (i == first)
            continue;

        c = isl_equality_alloc(isl_dim_copy(dim));
        isl_int_set_si(v, -1);
        isl_constraint_set_coefficient(c, isl_dim_in, i, v);
        isl_int_set_si(v, 1);
        isl_constraint_set_coefficient(c, isl_dim_out, i, v);
        bmap = isl_basic_map_add_constraint(bmap, c);
    }

    c = isl_equality_alloc(isl_dim_copy(dim));
    isl_int_set_si(v, -1);
    for (i = 0; i < wave_len; ++i)
        isl_constraint_set_coefficient(c, isl_dim_in, first + i, v);
    isl_int_set_si(v, 1);
    isl_constraint_set_coefficient(c, isl_dim_out, first, v);
    bmap = isl_basic_map_add_constraint(bmap, c);

    isl_dim_free(dim);
    isl_int_clear(v);

    return isl_map_from_basic_map(bmap);
}

static __isl_give isl_set *extract_host_domain(struct clast_user_stmt *u)
{
    return isl_set_from_cloog_domain(cloog_domain_copy(u->domain));
}

/* Extract the set of scattering dimension values for which the given
 * sequence of user statements is executed.
 * In principle, this set should be the same for each of the user
 * statements in the sequence, but we compute the union just to be safe.
 */
static __isl_give isl_set *extract_entire_host_domain(struct clast_user_stmt *u)
{
    struct clast_stmt *s;
    isl_set *host_domain = NULL;

    for (s = &u->stmt; s; s = s->next) {
        isl_set *set_i;

        assert(CLAST_STMT_IS_A(s, stmt_user));
        u = (struct clast_user_stmt *) s;

        set_i = extract_host_domain(u);

        if (!host_domain)
            host_domain = set_i;
        else
            host_domain = isl_set_union(host_domain, set_i);
        assert(host_domain);
    }

    return isl_set_coalesce(host_domain);
}

/* Extend "set" with unconstrained coordinates to a total length of "dst_len".
 */
static __isl_give isl_set *extend(__isl_take isl_set *set, int dst_len)
{
    int n_set;
    isl_dim *dim;
    isl_map *map;

    dim = isl_set_get_dim(set);
    n_set = isl_dim_size(dim, isl_dim_set);
    dim = isl_dim_drop(dim, isl_dim_set, 0, n_set);
    map = projection(dim, dst_len, n_set);
    map = isl_map_reverse(map);

    return isl_set_apply(set, map);
}

/* Extend the domain of "umap" with unconstrained coordinates
 * from an original length of "src_len" to a length of "dst_len".
 */
static __isl_give isl_union_map *extend_domain(__isl_take isl_union_map *umap,
    int src_len, int dst_len)
{
    isl_dim *dim;
    isl_map *proj;

    dim = isl_union_map_get_dim(umap);
    proj = project_out(dim, dst_len, src_len, dst_len - src_len);
    proj = isl_map_reverse(proj);
    return isl_union_map_apply_domain(umap, isl_union_map_from_map(proj));
}

/* Project "dom" onto its first "host_len" coordinates and the
 * coordinate "proj_dim" and simplify the result in the context of "host".
 */
static __isl_give isl_set *bounds_on(__isl_keep isl_union_set *dom,
    int len, int proj_dim, __isl_keep isl_set *host, int host_len)
{
    isl_dim *dim;
    isl_map *proj1, *proj2;
    isl_union_set *uhost;
    isl_set *bounds;

    dim = isl_union_set_get_dim(dom);
    proj1 = project_out(isl_dim_copy(dim), len, proj_dim + 1,
                        len - (proj_dim + 1));
    proj2 = project_out(isl_dim_copy(dim), proj_dim + 1,
                        host_len, proj_dim - host_len);
    proj1 = isl_map_apply_range(proj1, proj2);

    dom = isl_union_set_copy(dom);
    dom = isl_union_set_apply(dom, isl_union_map_from_map(proj1));

    proj1 = projection(isl_dim_copy(dim), host_len + 1, host_len);
    proj1 = isl_map_reverse(proj1);

    uhost = isl_union_set_from_set(isl_set_copy(host));
    uhost = isl_union_set_apply(uhost, isl_union_map_from_map(proj1));

    dom = isl_union_set_gist(dom, uhost);
    dom = isl_union_set_coalesce(dom);

    dim = isl_dim_add(dim, isl_dim_set, host_len + 1);
    bounds = isl_union_set_extract_set(dom, dim);
    isl_union_set_free(dom);

    return bounds;
}

/* Construct a set of the given number of coordinates ("len")
 * where the first loc->first + loc->tile_len coordinates are equated
 * to extra h%d and b%d parameters.
 */
static __isl_give isl_set *parametrization(struct localizer_info *loc,
    __isl_take isl_dim *dim, int len)
{
    int i;
    int param_len = loc->first + loc->tile_len;
    int offset;
    unsigned nparam;
    char name[20];
    isl_basic_set *bset;
    isl_constraint *c;
    isl_int v;

    if (len < param_len)
        param_len = len;
    nparam = isl_dim_size(dim, isl_dim_param);
    dim = isl_dim_add(dim, isl_dim_param, param_len);
    dim = isl_dim_add(dim, isl_dim_set, len);
    offset = nparam;
    for (i = 0; i < loc->first + 1 && offset + i < nparam + param_len; ++i) {
        snprintf(name, sizeof(name), "h%d", i);
        dim = isl_dim_set_name(dim, isl_dim_param, offset + i, name);
    }
    offset += loc->first + 1;
    for (i = 0; i < loc->tile_len - 1 && offset + i < nparam + param_len; ++i) {
        snprintf(name, sizeof(name), "b%d", i);
        dim = isl_dim_set_name(dim, isl_dim_param, offset + i, name);
    }

    isl_int_init(v);

    bset = isl_basic_set_universe(isl_dim_copy(dim));
    for (i = 0; i < param_len; ++i) {
        c = isl_equality_alloc(isl_dim_copy(dim));
        isl_int_set_si(v, -1);
        isl_constraint_set_coefficient(c, isl_dim_param, nparam + i, v);
        isl_int_set_si(v, 1);
        isl_constraint_set_coefficient(c, isl_dim_set, i, v);
        bset = isl_basic_set_add_constraint(bset, c);
    }
    isl_dim_free(dim);

    isl_int_clear(v);

    return isl_set_from_basic_set(bset);
}

/* Replace the first loc->first + loc_tile coordinates of the range
 * by parameters namd h%d and b%d.
 * "len" is the total number of coordinates in the range.
 */
static __isl_give isl_union_map *parameterize_range(struct localizer_info *loc,
    __isl_take isl_union_map *umap, int len)
{
    int remove;
    isl_dim *dim;
    isl_set *p;
    isl_map *proj;

    dim = isl_union_map_get_dim(umap);
    p = parametrization(loc, dim, len);
    umap = isl_union_map_intersect_range(umap, isl_union_set_from_set(p));

    remove = loc->first + loc->tile_len;
    if (len < remove)
        remove = len;
    dim = isl_union_map_get_dim(umap);
    proj = project_out(dim, len, 0, remove);
    umap = isl_union_map_apply_range(umap, isl_union_map_from_map(proj));

    return umap;
}

/* Relate points in the scheduled domain (according to "schedule")
 * to array elements accessed in those points.
 */
static __isl_give isl_union_map *compute_scheduled_accesses(
    struct localizer_info *loc, __isl_keep isl_union_map *schedule)
{
    isl_union_map *access;

    access = isl_union_map_union(isl_union_map_copy(loc->prog->read),
                                 isl_union_map_copy(loc->prog->write));
    access = isl_union_map_coalesce(access);
    access = isl_union_map_apply_domain(access, isl_union_map_copy(schedule));

    return access;
}

/* Project domain of "umap" of dimension "orig_len" to its first "target_len"
 * coordinates.
 */
static __isl_give isl_union_map *project_domain(__isl_take isl_union_map *umap,
    int orig_len, int target_len)
{
    isl_dim *dim;
    isl_map *proj;

    dim = isl_union_map_get_dim(umap);
    proj = projection(dim, orig_len, target_len);
    umap = isl_union_map_apply_domain(umap, isl_union_map_from_map(proj));
    umap = isl_union_map_coalesce(umap);

    return umap;
}

/* Construct an isl_map between two domains of dimensionality "len"
 * that increments the final coordinate (leaving all other coordinates
 * untouched).
 */
static __isl_give isl_map *compute_next(__isl_take isl_dim *dim, int len)
{
    int i;
    isl_constraint *c;
    isl_basic_map *bmap;
    isl_int v;

    isl_int_init(v);

    dim = isl_dim_add(dim, isl_dim_in, len);
    dim = isl_dim_add(dim, isl_dim_out, len);
    bmap = isl_basic_map_universe(isl_dim_copy(dim));

    for (i = 0; i < len - 1; ++i) {
        c = isl_equality_alloc(isl_dim_copy(dim));
        isl_int_set_si(v, -1);
        isl_constraint_set_coefficient(c, isl_dim_in, i, v);
        isl_int_set_si(v, 1);
        isl_constraint_set_coefficient(c, isl_dim_out, i, v);
        bmap = isl_basic_map_add_constraint(bmap, c);
    }

    c = isl_equality_alloc(isl_dim_copy(dim));
    isl_int_set_si(v, -1);
    isl_constraint_set_coefficient(c, isl_dim_in, len - 1, v);
    isl_constraint_set_constant(c, v);
    isl_int_set_si(v, 1);
    isl_constraint_set_coefficient(c, isl_dim_out, i, v);
    bmap = isl_basic_map_add_constraint(bmap, c);

    isl_dim_free(dim);

    isl_int_clear(v);

    return isl_map_from_basic_map(bmap);
}

/* Compute the set of elements accessed by a given iteration
 * that were not accessed by the previous iteration of the inner loop.
 */
static __isl_give isl_union_map *compute_copy_in(
    struct localizer_info *loc, __isl_take isl_union_map *access)
{
    isl_dim *dim;
    isl_map *next;
    isl_union_map *prev_access;
    isl_union_map *copy_in;

    dim = isl_union_map_get_dim(access);
    next = compute_next(dim, loc->gpu_len);
    prev_access = isl_union_map_apply_domain(isl_union_map_copy(access),
                                             isl_union_map_from_map(next));
    copy_in = isl_union_map_subtract(access, prev_access);
    copy_in = isl_union_map_compute_divs(copy_in);
    copy_in = isl_union_map_coalesce(copy_in);
    return copy_in;
}

/* Compute the set of elements accessed by a given iteration
 * that are not accessed by the next iteration of the inner loop.
 */
static __isl_give isl_union_map *compute_copy_out(
    struct localizer_info *loc, __isl_take isl_union_map *access)
{
    isl_dim *dim;
    isl_map *next;
    isl_union_map *next_access;
    isl_union_map *copy_out;

    dim = isl_union_map_get_dim(access);
    next = compute_next(dim, loc->gpu_len);
    next_access = isl_union_map_apply_range(isl_union_map_from_map(next),
                                             isl_union_map_copy(access));
    copy_out = isl_union_map_subtract(access, next_access);
    copy_out = isl_union_map_compute_divs(copy_out);
    copy_out = isl_union_map_coalesce(copy_out);
    return copy_out;
}

/* Compute the set of elements accessed by a given iteration
 * that are also accessed by the next iteration of the inner loop.
 */
static __isl_give isl_union_map *compute_reorganize(
    struct localizer_info *loc, __isl_take isl_union_map *access)
{
    isl_dim *dim;
    isl_map *next;
    isl_union_map *next_access;
    isl_union_map *reorganize;

    dim = isl_union_map_get_dim(access);
    next = compute_next(dim, loc->gpu_len);
    next_access = isl_union_map_apply_range(isl_union_map_from_map(next),
                                             isl_union_map_copy(access));
    reorganize = isl_union_map_intersect(access, next_access);
    reorganize = isl_union_map_compute_divs(reorganize);
    reorganize = isl_union_map_coalesce(reorganize);
    return reorganize;
}

/* Print (un)bundling code for the elements of an array identified by "set".
 */
static int print_bundle(struct localizer_info *loc, __isl_take isl_set *set,
    __isl_keep isl_set *context, int unbundle)
{
    int i;
    int nparam;
    CloogOptions *options;
    CloogDomain *domain;
    CloogDomain *cloog_context;
    CloogUnionDomain *ud;
    CloogInput *input;
    struct clast_stmt *stmt;
    const char *array_name;
    char *name;

    array_name = isl_set_get_tuple_name(set);
    name = malloc(sizeof("unbundle_") + strlen(array_name));
    sprintf(name, "%sbundle_%s", unbundle ? "un" : "", array_name);

    options = cloog_options_malloc(loc->state);
    options->language = LANGUAGE_C;
    options->strides = 1;

    nparam = isl_set_dim(set, isl_dim_param);
    ud = cloog_union_domain_alloc(nparam);
    set = isl_set_set_tuple_name(set, NULL);
    domain = cloog_domain_from_isl_set(isl_set_copy(set));
    ud = cloog_union_domain_add_domain(ud, name, domain, NULL, NULL);
    for (i = 0; i < nparam; ++i) {
        const char *name;
        name = isl_set_get_dim_name(set, isl_dim_param, i);
        ud = cloog_union_domain_set_name(ud, CLOOG_PARAM, i, name);
    }
    cloog_context = cloog_domain_from_isl_set(isl_set_copy(context));
    input = cloog_input_alloc(cloog_context, ud);

    stmt = cloog_clast_create_from_input(input, options);
    clast_pprint(loc->dst, stmt, loc->indent, options);

    free(name);
    isl_set_free(set);
    cloog_clast_free(stmt);
    cloog_options_free(options);

    return 0;
}

/* Print (un)bundling code for each array.
 */
static void bundle(struct localizer_info *loc, __isl_keep isl_union_map *access,
    __isl_keep isl_set *host_domain, int unbundle)
{
    int i, j;
    isl_set *context;
    int host_len = isl_set_dim(host_domain, isl_dim_set);
    int nparam = isl_set_dim(host_domain, isl_dim_param);
    char name[20];

    context = isl_set_copy(host_domain);
    context = isl_set_move_dims(context, isl_dim_param, nparam,
                                isl_dim_set, 0, host_len);
    for (i = 0; i < host_len; ++i) {
        snprintf(name, sizeof(name), "h%d", i);
        context = isl_set_set_dim_name(context, isl_dim_param, nparam + i, name);
    }

    print_indent(loc->dst, loc->indent);
    fprintf(loc->dst, "p = transfer_buffer;\n");
    for (i = 0; i < loc->n_array; ++i) {
        isl_dim *dim;
        isl_map *map;
        dim = isl_union_map_get_dim(access);
        dim = isl_dim_add(dim, isl_dim_in, loc->host_len);
        dim = isl_dim_add(dim, isl_dim_out, loc->array[i].dim);
        dim = isl_dim_set_tuple_name(dim, isl_dim_out, loc->array[i].name);
        map = isl_union_map_extract_map(access, dim);
        map = isl_map_move_dims(map, isl_dim_param, nparam,
                                isl_dim_in, 0, host_len);
        for (j = 0; j < host_len; ++j) {
            snprintf(name, sizeof(name), "h%d", j);
            map = isl_map_set_dim_name(map, isl_dim_param, nparam + j, name);
        }
        print_bundle(loc, isl_map_range(map), context, unbundle);
    }

    isl_set_free(context);
}

/* Increment loc->n if "c" is a lower bound (or an upper bound
 * if loc->upper is set.
 */
static int count_bound(__isl_take isl_constraint *c, void *user)
{
    struct localizer_info *loc = (struct localizer_info *)user;
    int dim = isl_constraint_dim(c, isl_dim_set);
    isl_int v;
    int sgn;

    isl_int_init(v);

    isl_constraint_get_coefficient(c, isl_dim_set, dim - 1, &v);
    sgn = isl_int_sgn(v);
    assert(sgn != 0);
    if ((loc->upper && sgn < 0) || (!loc->upper && sgn > 0) ||
            isl_constraint_is_equality(c))
        loc->n++;

    isl_int_clear(v);

    isl_constraint_free(c);
    return 0;
}

/* Print the given constraint as a bound on the final coordinate.
 */
static int print_bound(__isl_take isl_constraint *c, void *user)
{
    int i;
    struct localizer_info *loc = (struct localizer_info *)user;
    int nparam = isl_constraint_dim(c, isl_dim_param);
    int dim = isl_constraint_dim(c, isl_dim_set);
    isl_int v, t;
    int sgn;
    int first;
    const char *name;

    isl_int_init(v);

    isl_constraint_get_coefficient(c, isl_dim_set, dim - 1, &v);
    sgn = isl_int_sgn(v);
    assert(sgn != 0);
    if ((loc->upper && sgn < 0) || (!loc->upper && sgn > 0) ||
            isl_constraint_is_equality(c)) {
        isl_int_init(t);
        sgn = -sgn;
        isl_int_abs(v, v);
        if (loc->n)
            fprintf(loc->dst, ", ");
        if (isl_int_cmp_si(v, 1) > 0)
            fprintf(loc->dst, loc->upper ? "floord(" : "ceild(");
        first = 1;
        isl_constraint_get_constant(c, &t);
        if (!isl_int_is_zero(t)) {
            if (sgn < 0)
                isl_int_neg(t, t);
            isl_int_print(loc->dst, t, 0);
            first = 0;
        }
        for (i = 0; i < nparam; ++i) {
            isl_constraint_get_coefficient(c, isl_dim_param, i, &t);
            if (isl_int_is_zero(t))
                continue;
            if (sgn < 0)
                isl_int_neg(t, t);
            if (!first && isl_int_is_pos(t))
                fprintf(loc->dst, "+");
            isl_int_print(loc->dst, t, 0);
            name = isl_constraint_get_dim_name(c, isl_dim_param, i);
            fprintf(loc->dst, "*%s", name);
            first = 0;
        }
        for (i = 0; i < dim - 1; ++i) {
            isl_constraint_get_coefficient(c, isl_dim_set, i, &t);
            if (isl_int_is_zero(t))
                continue;
            if (sgn < 0)
                isl_int_neg(t, t);
            if (!first && isl_int_is_pos(t))
                fprintf(loc->dst, "+");
            isl_int_print(loc->dst, t, 0);
            fprintf(loc->dst, "*h%d", i);
            first = 0;
        }
        if (first)
            fprintf(loc->dst, "0");
        if (isl_int_cmp_si(v, 1) > 0) {
            fprintf(loc->dst, ", ");
            isl_int_print(loc->dst, v, 0);
            fprintf(loc->dst, ")");
        }
        if (loc->n)
            fprintf(loc->dst, ")");
        isl_int_clear(t);
        loc->n++;
    }

    isl_int_clear(v);

    isl_constraint_free(c);
    return 0;
}

/* Print the lower (or upper) bound on the final coordinate of "bset".
 */
static int print_basic_set_bound(__isl_take isl_basic_set *bset, void *user)
{
    int i;
    struct localizer_info *loc = (struct localizer_info *)user;

    loc->n = 0;
    isl_basic_set_foreach_constraint(bset, &count_bound, user);

    for (i = 0; i < loc->n - 1; ++i) {
        fprintf(loc->dst, loc->upper ? "min(" : "max(");
    }

    loc->n = 0;
    isl_basic_set_foreach_constraint(bset, &print_bound, user);

    isl_basic_set_free(bset);
    return 0;
}

/* Print the lower (or upper) bound on the final coordinate of "bound".
 */
static void print_bounds(struct localizer_info *loc,
    __isl_keep isl_set *bound, int upper)
{
    loc->upper = upper;

    assert(isl_set_n_basic_set(bound) == 1);
    isl_set_foreach_basic_set(bound, &print_basic_set_bound, loc);
}

static int extract_array_info(__isl_take isl_map *access, void *user)
{
    struct localizer_info *loc = (struct localizer_info *)user;
    const char *name;

    name = isl_map_get_tuple_name(access, isl_dim_out);
    loc->array[loc->n_array].name = strdup(name);
    loc->array[loc->n_array].dim = isl_map_dim(access, isl_dim_out);
    loc->array[loc->n_array].transfer_size = isl_map_card(access);

    loc->n_array++;
    return 0;
}

static void collect_array_info(struct localizer_info *loc)
{
    isl_union_map *access;

    access = compute_scheduled_accesses(loc, loc->sched);
    access = project_domain(access, loc->tiled_len, loc->host_len);
    loc->n_array = isl_union_map_n_map(access);
    loc->array = isl_alloc_array(ctx, struct array_info, loc->n_array);
    assert(loc->array);
    loc->n_array = 0;
    isl_union_map_foreach_map(access, &extract_array_info, loc);
    isl_union_map_free(access);
}

static void free_array_info(struct localizer_info *loc)
{
    int i;

    for (i = 0; i < loc->n_array; ++i) {
        free(loc->array[i].name);
        isl_pw_qpolynomial_free(loc->array[i].transfer_size);
    }
    free(loc->array);
}

/* Compute a piecewise qpolynomial that maps an element of the given
 * relation to the number of other image elements associated to the
 * same domain element that are lexicographically smaller than
 * the range element.
 */
static isl_union_pw_qpolynomial *compute_ranking(__isl_keep isl_union_map *umap)
{
    isl_union_set *range;
    isl_union_map *lex_gt;
    isl_union_pw_qpolynomial *rank;
    isl_union_map *lex_smaller_image, *domain_map, *in_image;

    range = isl_union_map_range(isl_union_map_copy(umap));
    lex_gt = isl_union_set_lex_gt_union_set(isl_union_set_copy(range),
                                            isl_union_set_copy(range));
    isl_union_set_free(range);
    lex_smaller_image = isl_union_map_range_map(isl_union_map_copy(umap));
    lex_smaller_image = isl_union_map_apply_range(lex_smaller_image, lex_gt);
    domain_map = isl_union_map_domain_map(isl_union_map_copy(umap));
    in_image = isl_union_map_apply_range(domain_map, isl_union_map_copy(umap));
    lex_smaller_image = isl_union_map_intersect(lex_smaller_image, in_image);
    rank = isl_union_map_card(lex_smaller_image);
    return rank;
}

/* Turn the first loc->first + loc->tile_len coordinates of "block_domain"
 * into parameters.
 */
static __isl_give isl_set *parameterize_block_domain(struct localizer_info *loc,
    __isl_take isl_set *block_domain)
{
    int i;
    int host_len = loc->first + 1;
    int block_len = loc->tile_len - 1;
    int nparam = isl_set_dim(block_domain, isl_dim_param);
    char name[20];

    block_domain = isl_set_move_dims(block_domain, isl_dim_param, nparam,
                                     isl_dim_set, 0, host_len + block_len);
    for (i = 0; i < host_len; ++i) {
        snprintf(name, sizeof(name), "h%d", i);
        block_domain = isl_set_set_dim_name(block_domain, isl_dim_param,
                                            nparam + i, name);
    }
    for (i = 0; i < block_len; ++i) {
        snprintf(name, sizeof(name), "b%d", i);
        block_domain = isl_set_set_dim_name(block_domain, isl_dim_param,
                                            nparam + host_len + i, name);
    }

    return block_domain;
}

struct sched_data {
    isl_union_map *umap;
    int max_len;
    const char *type;
    int pos;
};

/* Extend the range of "map" with zero coordinates up to a size of "len".
 */
static __isl_give isl_map *extend_range(__isl_take isl_map *map, int len)
{
    int i;
    int n_out = isl_map_dim(map, isl_dim_out);

    map = isl_map_add_dims(map, isl_dim_out, len - n_out);
    for (i = 0; i < len - n_out; ++i)
        map = isl_map_fix_si(map, isl_dim_out, n_out + i, 0);

    return map;
}

/* Given an access relation of the form
 *
 *      [h b s] -> A[i]
 *
 * construct a schedule that schedules each at
 *
 *      h b s "data->pos" i 0
 *
 * The domain of the schedule is
 *
 *      "data->type"_A[s i]
 *
 * h and b are available as parameters in the corresponding statements.
 *
 * We essentially construct an identity map on a flattened version
 * of the input relation and then add and remove the appropriate dimensions
 * from domain and range.
 */
static int create_sched(__isl_take isl_map *access, void *user)
{
    isl_ctx *ctx;
    struct sched_data *data = (struct sched_data *)user;
    const char *array_name;
    char *name;
    isl_dim *dim;
    isl_map *proj;
    isl_map *sched;
    isl_set *set;
    int sched_len, arr_len;

    ctx = isl_map_get_ctx(access);

    array_name = isl_map_get_tuple_name(access, isl_dim_out);
    name = isl_alloc_array(ctx, char,
                        strlen(data->type) + sizeof("_") + strlen(array_name));
    sprintf(name, "%s_%s", data->type, array_name);

    sched_len = isl_map_dim(access, isl_dim_in);
    arr_len = isl_map_dim(access, isl_dim_out);
    set = isl_map_wrap(access);
    set = isl_set_flatten(set);
    dim = isl_set_get_dim(set);
    dim = isl_dim_drop(dim, isl_dim_set, 0, sched_len + arr_len);
    sched = project_out(isl_dim_copy(dim), sched_len + arr_len, 0, sched_len - 1);
    sched = isl_map_reverse(sched);
    sched = isl_map_intersect_range(sched, set);

    proj = project_out(dim, sched_len + arr_len + 1, sched_len, 1);
    proj = isl_map_reverse(proj);
    proj = isl_map_fix_si(proj, isl_dim_out, sched_len, data->pos);
    sched = isl_map_apply_range(sched, proj);
    sched = extend_range(sched, data->max_len);
    sched = isl_map_set_tuple_name(sched, isl_dim_in, name);

    data->umap = isl_union_map_add_map(data->umap, sched);

    free(name);

    return 0;
}

/* Construct a schedule for each of the accesses in "copy"
 * according to create_sched.
 */
static __isl_give isl_union_map *add_copy_sched(
    __isl_take isl_union_map *sched, int max_len,
    __isl_keep isl_union_map *copy, const char *type, int pos)
{
    isl_dim *dim;
    struct sched_data sd;

    dim = isl_union_map_get_dim(sched);
    sd.max_len = max_len;
    sd.umap = isl_union_map_empty(dim);
    sd.type = type;
    sd.pos = pos;
    isl_union_map_foreach_map(copy, &create_sched, &sd);

    return isl_union_map_union(sched, sd.umap);
}

/* For each iteration i of "seq_domain", add a synchronization
 * of "type" at
 *
 *      i "pos" 0
 */
static __isl_give isl_union_map *add_sync_sched(
    __isl_take isl_union_map *sched, int max_len,
    __isl_keep isl_set *seq_domain, const char *type, int pos)
{
    int len = isl_set_dim(seq_domain, isl_dim_set);
    isl_map *map;

    seq_domain = isl_set_copy(seq_domain);
    seq_domain = isl_set_add_dims(seq_domain, isl_dim_set, 1);
    seq_domain = isl_set_fix_si(seq_domain, isl_dim_set, len, pos);

    map = isl_map_from_range(seq_domain);
    map = isl_map_set_tuple_name(map, isl_dim_in, type);

    map = extend_range(map, max_len);

    return isl_union_map_union(sched, isl_union_map_from_map(map));
}

/* Replace the thread loops t in "schedule" by
 *
 *      t = 512 t' + threadIdx_{x,y,z}
 */
static __isl_give isl_union_map *insert_threading(struct localizer_info *loc,
    __isl_take isl_union_map *schedule, int max_len)
{
    int i;
    isl_dim *dim;
    isl_basic_map *bmap;
    isl_union_map *umap;
    int nparam;
    isl_int v;
    isl_constraint *c;

    isl_int_init(v);

    dim = isl_union_map_get_dim(schedule);
    nparam = isl_dim_size(dim, isl_dim_param);

    assert(loc->tile_len - 1 <= 3);
    dim = isl_dim_add(dim, isl_dim_param, loc->tile_len - 1);
    for (i = 0; i < loc->tile_len - 1; ++i) {
        const char *dims[] = { "threadIdx_x", "threadIdx_y", "threadIdx_z" };
        dim = isl_dim_set_name(dim, isl_dim_param, nparam + i, dims[i]);
    }

    dim = isl_dim_add(dim, isl_dim_in, max_len);
    dim = isl_dim_add(dim, isl_dim_out, max_len);
    bmap = isl_basic_map_universe(isl_dim_copy(dim));
    for (i = 0; i < max_len; ++i) {
        if (i < loc->gpu_len + 1 || i >= loc->gpu_len + loc->tile_len) {
            c = isl_equality_alloc(isl_dim_copy(dim));
            isl_int_set_si(v, -1);
            isl_constraint_set_coefficient(c, isl_dim_in, i, v);
            isl_int_set_si(v, 1);
            isl_constraint_set_coefficient(c, isl_dim_out, i, v);
            bmap = isl_basic_map_add_constraint(bmap, c);
        } else {
            c = isl_equality_alloc(isl_dim_copy(dim));
            isl_int_set_si(v, -1);
            isl_constraint_set_coefficient(c, isl_dim_in, i, v);
            isl_int_set_si(v, 512);
            isl_constraint_set_coefficient(c, isl_dim_out, i, v);
            isl_int_set_si(v, 1);
            isl_constraint_set_coefficient(c, isl_dim_param,
                                           nparam + i - (loc->gpu_len + 1), v);
            bmap = isl_basic_map_add_constraint(bmap, c);
        }
    }

    isl_dim_free(dim);
    isl_int_clear(v);

    umap = isl_union_map_from_map(isl_map_from_basic_map(bmap));

    return isl_union_map_apply_range(schedule, umap);
}

/* Insert a threadIdx_x parameter with constraints 0 <= threadIdx_x < 512.
 */
__isl_give isl_set *context_insert_threading(struct localizer_info *loc,
    __isl_take isl_set *context)
{
    int i;
    int nparam;
    isl_dim *dim;
    isl_basic_set *bset;
    isl_int v;
    isl_constraint *c;

    isl_int_init(v);

    nparam = isl_set_dim(context, isl_dim_param);
    context = isl_set_add_dims(context, isl_dim_param, loc->tile_len - 1);
    for (i = 0; i < loc->tile_len - 1; ++i) {
        const char *dims[] = { "threadIdx_x", "threadIdx_y", "threadIdx_z" };
        context = isl_set_set_dim_name(context, isl_dim_param,
                                            nparam + i, dims[i]);
    }

    dim = isl_set_get_dim(context);
    bset = isl_basic_set_universe(isl_dim_copy(dim));

    for (i = 0; i < loc->tile_len - 1; ++i) {
        c = isl_inequality_alloc(isl_dim_copy(dim));
        isl_int_set_si(v, 1);
        isl_constraint_set_coefficient(c, isl_dim_param, nparam + i, v);
        bset = isl_basic_set_add_constraint(bset, c);

        c = isl_inequality_alloc(isl_dim_copy(dim));
        isl_int_set_si(v, -1);
        isl_constraint_set_coefficient(c, isl_dim_param, nparam + i, v);
        isl_int_set_si(v, 511);
        isl_constraint_set_constant(c, v);
        bset = isl_basic_set_add_constraint(bset, c);
    }

    isl_dim_free(dim);
    isl_int_clear(v);

    context = isl_set_intersect(context, isl_set_from_basic_set(bset));

    return context;
}

/* Generate and print code for the body of the kernel.
 * We first insert an extra coordinate after the sequential loop on the GPU
 * to impose a relative order within this sequential loop between the
 * original statements, copying and synchronization.
 * In particular, we have
 *
 *      -2      copy_in
 *      -1      sync
 *       0      original statements
 *       1      sync
 *       2      copy_out
 *       3      reorganize
 *       4      sync
 *       5      swap buffers
 *
 * The copy statements are included before the parallel loop is
 * wrapped over the threads.
 * The synchronization statements are included afterward.
 */
static void print_kernel_body(struct localizer_info *loc,
    __isl_keep isl_set *host_domain, __isl_keep isl_set *block_domain,
    __isl_take isl_union_map *copy_in, __isl_take isl_union_map *copy_out,
    __isl_take isl_union_map *reorganize)
{
    int i;
    isl_dim *dim;
    isl_map *proj;
    isl_union_set *sched_domain;
    isl_set *seq_domain;
    CloogOptions *options;
    CloogDomain *cloog_context;
    CloogUnionDomain *ud;
    CloogInput *input;
    struct clast_stmt *stmt;
    int max_len, offset;
    char name[20];
    isl_union_map *sched;

    sched = isl_union_map_copy(loc->sched);

    sched_domain = isl_union_map_range(isl_union_map_copy(sched));
    dim = isl_union_map_get_dim(sched);
    dim = isl_dim_add(dim, isl_dim_set, loc->tiled_len);
    seq_domain = isl_union_set_extract_set(sched_domain, dim);
    isl_union_set_free(sched_domain);
    seq_domain = isl_set_project_out(seq_domain, isl_dim_set,
                                loc->gpu_len, loc->tiled_len - loc->gpu_len);

    dim = isl_union_map_get_dim(sched);
    proj = project_out(dim, loc->tiled_len + 1, loc->gpu_len, 1);
    proj = isl_map_reverse(proj);
    proj = isl_map_fix_si(proj, isl_dim_out, loc->gpu_len, 0);
    sched = isl_union_map_apply_range(sched, isl_union_map_from_map(proj));

    max_len = loc->tiled_len + 1;
    dim = isl_union_map_get_dim(sched);
    proj = project_out(dim, max_len, 0, 0);
    for (i = 0; i < loc->n_array; ++i)
        if (loc->gpu_len + 1 + loc->array[i].dim > max_len)
            max_len = loc->gpu_len + 1 + loc->array[i].dim;
    proj = extend_range(proj, max_len);
    sched = isl_union_map_apply_range(sched, isl_union_map_from_map(proj));

    sched = add_copy_sched(sched, max_len, copy_in, "copyin", -2);
    sched = add_copy_sched(sched, max_len, copy_out, "copyout", 2);
    sched = add_copy_sched(sched, max_len, reorganize, "reorganize", 3);

    sched = insert_threading(loc, sched, max_len);

    sched = add_sync_sched(sched, max_len, seq_domain, "sync1", -1);
    sched = add_sync_sched(sched, max_len, seq_domain, "sync2", 1);
    sched = add_sync_sched(sched, max_len, seq_domain, "sync3", 4);
    sched = add_sync_sched(sched, max_len, seq_domain, "swap_buffers", 5);

    sched = parameterize_range(loc, sched, max_len);

    isl_set_free(seq_domain);

    block_domain = parameterize_block_domain(loc, isl_set_copy(block_domain));
    block_domain = context_insert_threading(loc, block_domain);

    max_len -= loc->gpu_len - 1;

    fprintf(loc->kernel_c, "    int s");
    for (i = 0; i < loc->tile_len - 1; ++i)
        fprintf(loc->kernel_c, ", t%d", i);
    for (i = 0; i < max_len - (loc->tile_len + 1); ++i)
        fprintf(loc->kernel_c, ", l%d", i);
    fprintf(loc->kernel_c, ";\n");
    fprintf(loc->kernel_c, "\n");

    options = cloog_options_malloc(loc->state);
    options->language = LANGUAGE_C;
    options->strides = 1;
    options->f = 2 + (loc->tile_len - 1);
    options->l = max_len;

    sched = isl_union_map_align_params(sched, isl_set_get_dim(block_domain));
    ud = cloog_union_domain_from_isl_union_map(sched);

    offset = 0;
    ud = cloog_union_domain_set_name(ud, CLOOG_SCAT, offset, "s");
    offset++;
    ud = cloog_union_domain_set_name(ud, CLOOG_SCAT, offset, "dummy");
    offset++;
    for (i = 0; i < loc->tile_len - 1; ++i) {
        snprintf(name, sizeof(name), "t%d", i);
        ud = cloog_union_domain_set_name(ud, CLOOG_SCAT, offset + i, name);
    }
    offset += loc->tile_len - 1;
    for (i = 0; i < max_len - offset; ++i) {
        snprintf(name, sizeof(name), "l%d", i);
        ud = cloog_union_domain_set_name(ud, CLOOG_SCAT, offset + i, name);
    }

    cloog_context = cloog_domain_from_isl_set(block_domain);
    input = cloog_input_alloc(cloog_context, ud);

    stmt = cloog_clast_create_from_input(input, options);
    clast_pprint(loc->kernel_c, stmt, 4, options);

    cloog_clast_free(stmt);
    cloog_options_free(options);

    isl_union_map_free(copy_in);
    isl_union_map_free(copy_out);
    isl_union_map_free(reorganize);
}

/* Print the kernel invocation, the kernel header and the kernel body
 * (through print_kernel_body).
 * host_domain is the set of iterations on the host for which the kernel
 * is invoked.
 * bounds are the bounds on the block ids.
 * block_domain is the set of host iterations and block ids for which
 * the kernel is run.
 * The domains of copy_in, copy_out and reorganize contain
 * the host iterations, the block ids and the sequential loop on the gpu.
 */
static void print_kernel(struct localizer_info *loc,
    __isl_keep isl_set *host_domain, __isl_keep isl_set *block_domain,
    __isl_keep isl_set **bounds,
    __isl_take isl_union_map *copy_in, __isl_take isl_union_map *copy_out,
    __isl_take isl_union_map *reorganize)
{
    int i, j, k;
    char name[20];
    unsigned nparam;
    FILE *host_file;

    /* for now */
    assert(loc->tile_len == 2);
    print_indent(loc->dst, loc->indent);
    fprintf(loc->dst, "dim3 k%d_dimBlock(512);\n", loc->kernel_id);
    print_indent(loc->dst, loc->indent);
    fprintf(loc->dst, "dim3 k%d_dimGrid((", loc->kernel_id);
    print_bounds(loc, bounds[0], 1);
    fprintf(loc->dst, ")-(");
    print_bounds(loc, bounds[0], 0);
    fprintf(loc->dst, ")+1);\n");

    print_indent(loc->dst, loc->indent);
    fprintf(loc->dst, "kernel%d <<<k%d_dimGrid, k%d_dimBlock>>> (dev_buffer",
            loc->kernel_id, loc->kernel_id, loc->kernel_id);
    fprintf(loc->kernel_c,
            "__global__ void kernel%d(%s *dev_buffer", loc->kernel_id, loc->type);
    fprintf(loc->kernel_h,
            "__global__ void kernel%d(%s *dev_buffer", loc->kernel_id, loc->type);
    nparam = isl_set_dim(host_domain, isl_dim_param);
    for (i = 0; i < nparam; ++i) {
        const char *name = isl_set_get_dim_name(host_domain, isl_dim_param, i);
        fprintf(loc->dst, ", %s", name);
        fprintf(loc->kernel_c, ", int %s", name);
        fprintf(loc->kernel_h, ", int %s", name);
    }
    for (i = 0; i < loc->host_len; ++i) {
        fprintf(loc->dst, ", h%d", i);
        fprintf(loc->kernel_c, ", int h%d", i);
        fprintf(loc->kernel_h, ", int h%d", i);
    }
    fprintf(loc->dst, ");\n");
    fprintf(loc->kernel_c, ")\n{\n");
    fprintf(loc->kernel_h, ");\n");
    host_file = loc->dst;
    loc->dst = loc->kernel_c;
    assert(loc->tile_len - 1 <= 2);
    for (i = 0; i < loc->tile_len - 1; ++i) {
        const char *dims[] = { "x", "y" };
        fprintf(loc->kernel_c, "    int b%d = ", i);
        print_bounds(loc, bounds[0], 0);
        fprintf(loc->kernel_c, " + blockIdx.%s;\n", dims[i]);
    }
    fprintf(loc->kernel_c, "    __shared__ %s L[2][", loc->type);
    isl_pw_qpolynomial_fold_print(loc->max_shared_size, loc->kernel_c,
                                  ISL_FORMAT_C);
    fprintf(loc->kernel_c, "];\n");
    for (i = 0; i < loc->n_array; ++i) {
        for (j = 0; j < 2; ++j) {
            fprintf(loc->kernel_c, "    %s *L%s_%s = L[%d]",
                    loc->type, j ? "2" : "", loc->array[i].name, j);
            for (k = 0; k < i; ++k) {
                fprintf(loc->kernel_c, " + (");
                isl_pw_qpolynomial_fold_print(loc->array[k].shared_size,
                                         loc->kernel_c, ISL_FORMAT_C);
                fprintf(loc->kernel_c, ")");
            }
            fprintf(loc->kernel_c, ";\n");
        }
        fprintf(loc->kernel_c, "    %s *buffer_%s = dev_buffer",
                loc->type, loc->array[i].name);
        for (k = 0; k < i; ++k) {
            isl_pw_qpolynomial *transfer_size;
            transfer_size = isl_pw_qpolynomial_copy(loc->array[k].transfer_size);
            transfer_size = isl_pw_qpolynomial_gist(transfer_size,
                                                    isl_set_copy(host_domain));
            for (j = 0; j < loc->first + 1; ++j) {
                snprintf(name, sizeof(name), "h%d", j);
                transfer_size = isl_pw_qpolynomial_set_dim_name(
                    transfer_size, isl_dim_set, j, name);
            }

            fprintf(loc->kernel_c, " + (");
            isl_pw_qpolynomial_print(transfer_size,
                                     loc->kernel_c, ISL_FORMAT_C);
            fprintf(loc->kernel_c, ")");

            isl_pw_qpolynomial_free(transfer_size);
        }
        fprintf(loc->kernel_c, ";\n");
    }
    loc->dst = host_file;

    print_kernel_body(loc, host_domain, block_domain,
                      copy_in, copy_out, reorganize);
    fprintf(loc->kernel_c, "}\n");
}

/* Extract the access in stmt->text starting at position identifier
 * and of length identifier_len, and rewrite the index expression A[i]
 * to L_A[rho(i)].
 * 
 * The access in C notation is first copied to "buffer" (which
 * has been allocated by the caller and should be of sufficient size)
 * and slightly modified to a map in isl notation.
 * This string is then parsed by isl.
 * The result is a mapping from the original iteration domain to array space,
 * D -> A.
 * rho maps a pair of scattered dimension and array index to a number,
 * [S -> A] -> N.
 * time_loop_proj maps the iteration domain to the same number of scattering
 * dimensions as used in rho, D -> S.
 *
 * We first construct a map [D -> A] -> [S -> A] from time_loop_proj
 * and an identity relation on A.
 * From the access relation, we construct a mapping D -> [D -> A].
 * Combining the mappings D -> [D -> A], [D -> A] -> [S -> A] and
 * [S -> A] -> N, we obtain a mapping D -> N.
 * Note that each of these mapping is single-valued.
 */
static void extract_access(struct localizer_info *loc, Stmt *stmt,
    char *buffer, int identifier, int identifier_len,
    int access, int access_len, char *name,
    __isl_keep isl_union_pw_qpolynomial *rho,
    __isl_keep isl_map *time_loop_proj)
{
    int i;
    int pos = 0;
    isl_ctx *ctx;
    isl_dim *dim;
    isl_map *map, *id;
    isl_union_map *umap;
    isl_pw_qpolynomial *pwqp;

    ctx = isl_union_map_get_ctx(loc->sched);
    pos += sprintf(buffer, "[");
    for (i = 0; i < loc->prog->npar; ++i) {
        if (i)
            pos += sprintf(buffer + pos, ",");
        pos += sprintf(buffer + pos, "%s", loc->prog->params[i]);
    }
    pos += sprintf(buffer + pos, "] -> { %s[", name);
    for (i = 0; i < stmt->dim; ++i) {
        if (i)
            pos += sprintf(buffer + pos, ",");
        pos += sprintf(buffer + pos, "%s", stmt->iterators[i]);
    }
    pos += sprintf(buffer + pos, "] -> ");
    memcpy(buffer + pos, stmt->text + identifier, identifier_len);
    pos += identifier_len;
    pos += sprintf(buffer + pos, "[");
    for (i = 0; i < access_len; ++i) {
        if (stmt->text[access + i] == ']') {
            buffer[pos++] = ',';
            ++i;
        } else
            buffer[pos++] = stmt->text[access + i];
    }
    pos += sprintf(buffer + pos, "] }");
    map = isl_map_read_from_str(ctx, buffer, -1);

    dim = isl_map_get_dim(map);
    dim = isl_dim_range(dim);
    id = isl_map_identity(dim);
    umap = isl_union_map_product(
                isl_union_map_from_map(isl_map_copy(time_loop_proj)),
                isl_union_map_from_map(id));

    map = isl_map_domain_map(map);
    map = isl_map_reverse(map);
    umap = isl_union_map_apply_range(isl_union_map_from_map(map), umap);

    rho = isl_union_map_apply_union_pw_qpolynomial(umap,
                                            isl_union_pw_qpolynomial_copy(rho));

    dim = isl_union_pw_qpolynomial_get_dim(rho);
    dim = isl_dim_add(dim, isl_dim_set, stmt->dim);
    dim = isl_dim_set_tuple_name(dim, isl_dim_set, name);
    pwqp = isl_union_pw_qpolynomial_extract_pw_qpolynomial(rho, dim);
    pwqp = isl_pw_qpolynomial_coalesce(pwqp);

    fprintf(loc->kernel_c, "L_");
    fwrite(stmt->text + identifier, 1, identifier_len, loc->kernel_c);
    fprintf(loc->kernel_c, "[");
    isl_pw_qpolynomial_print(pwqp, loc->kernel_c, ISL_FORMAT_C);
    fprintf(loc->kernel_c, "]");
    isl_pw_qpolynomial_free(pwqp);

    isl_union_pw_qpolynomial_free(rho);
}

/* Print stmt->text to loc->kernel_c, replacing each access A[i] by L_A[rho(i)].
 */
static void convert_accesses(struct localizer_info *loc, Stmt *stmt, char *name,
    __isl_keep isl_union_pw_qpolynomial *rho,
    __isl_keep isl_map *time_loop_proj)
{
    int i, j;
    isl_ctx *ctx;
    size_t text_len = strlen(stmt->text);
    size_t len = 50;
    char *buffer;
    int printed = 0;
    int identifier = -1;
    int end = -1;

    ctx = isl_map_get_ctx(time_loop_proj);
    for (i = 0; i < loc->prog->npar; ++i)
        len += strlen(loc->prog->params[i]);
    for (i = 0; i < stmt->dim; ++i)
        len += strlen(stmt->iterators[i]);
    buffer = isl_alloc_array(ctx, char, len);
    assert(buffer);

    for (i = 0; i < text_len; ++i) {
        if (identifier < 0 && isalpha(stmt->text[i])) {
            fwrite(stmt->text + printed, 1, i - printed, loc->kernel_c);
            identifier = i;
            end = -1;
        } else if (identifier >= 0 && end < 0 && isalnum(stmt->text[i]))
            continue;
        else if (identifier >= 0 && end < 0 && isblank(stmt->text[i]))
            end = i;
        else if (identifier >= 0 && end >= 0 && isblank(stmt->text[i]))
            continue;
        else if (identifier >= 0 && stmt->text[i] == '[') {
            if (end < 0)
                end = i;
            for (j = i + 1; j < text_len; ++j)
                if (stmt->text[j] == ']' && stmt->text[j + 1] != '[')
                    break;
            extract_access(loc, stmt, buffer, identifier, end - identifier,
                            i + 1, j - i - 1, name, rho, time_loop_proj);
            i = j;
            end = identifier = -1;
            printed = i + 1;
        } else {
            end = identifier = -1;
        }
    }
    fwrite(stmt->text + printed, 1, text_len - printed, loc->kernel_c);

    free(buffer);
}

static void print_transfer_size_init(struct localizer_info *loc,
    __isl_keep isl_pw_qpolynomial *size)
{
    int i;
    char name[20];

    for (i = 0; i < loc->host_len; ++i) {
        snprintf(name, sizeof(name), "h%d", i);
        size = isl_pw_qpolynomial_set_dim_name(size, isl_dim_set, i, name);
    }
    print_indent(loc->dst, loc->indent);
    fprintf(loc->dst, "transfer_size = (");
    isl_pw_qpolynomial_print(size, loc->dst, ISL_FORMAT_C);
    fprintf(loc->dst, ") * sizeof(%s);\n", loc->type);
}

/* Project schedule onto first loc->gpu_len coordinates.
 */
static __isl_give isl_union_map *compute_time_loop_proj(
    struct localizer_info *loc, __isl_keep isl_union_map *sched)
{
    isl_dim *dim;
    isl_map *proj;

    sched = isl_union_map_copy(sched);
    dim = isl_union_map_get_dim(sched);
    proj = project_out(dim, loc->tiled_len, loc->gpu_len,
                        loc->tiled_len - loc->gpu_len);
    sched = isl_union_map_apply_range(sched, isl_union_map_from_map(proj));

    return sched;
}

static void print_statement_define(struct localizer_info *loc,
    __isl_keep isl_union_map *sched, __isl_keep isl_union_pw_qpolynomial *rho)
{
    int i, j;
    char name[20];
    isl_union_map *time_loop_proj;

    time_loop_proj = compute_time_loop_proj(loc, sched);

    for (i = 0; i < loc->prog->nstmts; ++i) {
        Stmt *stmt = &loc->prog->stmts[i];
        isl_dim *dim;
        isl_map *proj_i;

        snprintf(name, sizeof(name), "S_%d", i);

        fprintf(loc->kernel_c, "#define %s(", name);
        for (j = 0; j < stmt->dim; ++j) {
            if (j)
                fprintf(loc->kernel_c, ",");
            fprintf(loc->kernel_c, "%s", stmt->iterators[j]);
        }
        fprintf(loc->kernel_c, ") %s_(", name);
        for (j = 0; j < stmt->dim; ++j) {
            if (j)
                fprintf(loc->kernel_c, ",");
            fprintf(loc->kernel_c, "(%s)", stmt->iterators[j]);
        }
        fprintf(loc->kernel_c, ")\n");

        fprintf(loc->kernel_c, "#define %s_(", name);
        for (j = 0; j < stmt->dim; ++j) {
            if (j)
                fprintf(loc->kernel_c, ",");
            fprintf(loc->kernel_c, "%s", stmt->iterators[j]);
        }
        fprintf(loc->kernel_c, ") ");

        dim = isl_union_map_get_dim(time_loop_proj);
        dim = isl_dim_add(dim, isl_dim_in, stmt->dim);
        dim = isl_dim_set_tuple_name(dim, isl_dim_in, name);
        dim = isl_dim_add(dim, isl_dim_out, loc->gpu_len);
        proj_i = isl_union_map_extract_map(time_loop_proj, dim);

        convert_accesses(loc, stmt, name, rho, proj_i);

        isl_map_free(proj_i);

        fprintf(loc->kernel_c, "\n");
    }

    isl_union_map_free(time_loop_proj);
}

/* Print a #define for copy or reorganize statement "stmt_name".
 * In particular for each array A, print
 *
 *  #define stmt_name_A(s,i) lhs_array[lhs(s,i)] = rhs_array[rhs(s,i)]
 *
 * where s is the sequential loopon the GPU and i is the array index.
 * The host iterators h%d and the block ids b%d are treated as parameters.
 */
static void print_copy_define(struct localizer_info *loc,
    __isl_keep isl_union_map *map, const char *stmt_name,
    const char *lhs_array, __isl_keep isl_union_pw_qpolynomial *lhs,
    const char *rhs_array, __isl_keep isl_union_pw_qpolynomial *rhs)
{
    int i, j;
    isl_union_set *dom;
    char name[20];
    int block_len = loc->first + loc->tile_len;

    dom = isl_union_map_wrap(isl_union_map_copy(map));
    lhs = isl_union_pw_qpolynomial_copy(lhs);
    lhs = isl_union_pw_qpolynomial_gist(lhs, isl_union_set_copy(dom));
    lhs = isl_union_pw_qpolynomial_coalesce(lhs);
    rhs = isl_union_pw_qpolynomial_copy(rhs);
    rhs = isl_union_pw_qpolynomial_gist(rhs, dom);
    rhs = isl_union_pw_qpolynomial_coalesce(rhs);
    for (i = 0; i < loc->n_array; ++i) {
        isl_dim *dim;
        isl_pw_qpolynomial *rhs_i;
        isl_pw_qpolynomial *lhs_i;

        dim = isl_union_pw_qpolynomial_get_dim(rhs);
        dim = isl_dim_add(dim, isl_dim_in, loc->gpu_len);
        dim = isl_dim_add(dim, isl_dim_out, loc->array[i].dim);
        dim = isl_dim_set_tuple_name(dim, isl_dim_out, loc->array[i].name);
        dim = isl_dim_wrap(dim);
        rhs_i = isl_union_pw_qpolynomial_extract_pw_qpolynomial(rhs, dim);

        dim = isl_union_pw_qpolynomial_get_dim(lhs);
        dim = isl_dim_add(dim, isl_dim_in, loc->gpu_len);
        dim = isl_dim_add(dim, isl_dim_out, loc->array[i].dim);
        dim = isl_dim_set_tuple_name(dim, isl_dim_out, loc->array[i].name);
        dim = isl_dim_wrap(dim);
        lhs_i = isl_union_pw_qpolynomial_extract_pw_qpolynomial(lhs, dim);

        fprintf(loc->kernel_c, "#define %s_%s(", stmt_name, loc->array[i].name);
        for (j = 0; j < loc->host_len; ++j) {
            snprintf(name, sizeof(name), "h%d", j);
            rhs_i = isl_pw_qpolynomial_set_dim_name(rhs_i,
                                              isl_dim_set, j, name);
            lhs_i = isl_pw_qpolynomial_set_dim_name(lhs_i,
                                              isl_dim_set, j, name);
        }
        for (j = 0; j < loc->tile_len - 1; ++j) {
            snprintf(name, sizeof(name), "b%d", j);
            rhs_i = isl_pw_qpolynomial_set_dim_name(rhs_i,
                                          isl_dim_set, loc->host_len + j, name);
            lhs_i = isl_pw_qpolynomial_set_dim_name(lhs_i,
                                          isl_dim_set, loc->host_len + j, name);
        }
        for (j = 0; j < 1 + loc->array[i].dim; ++j) {
            snprintf(name, sizeof(name), "i%d", j);
            rhs_i = isl_pw_qpolynomial_set_dim_name(rhs_i,
                                              isl_dim_set, block_len + j, name);
            lhs_i = isl_pw_qpolynomial_set_dim_name(lhs_i,
                                              isl_dim_set, block_len + j, name);
            if (j)
                fprintf(loc->kernel_c, ",");
            fprintf(loc->kernel_c, "i%d", j);
        }
        fprintf(loc->kernel_c, ") %s_%s_(", stmt_name, loc->array[i].name);
        for (j = 0; j < 1 + loc->array[i].dim; ++j) {
            if (j)
                fprintf(loc->kernel_c, ",");
            fprintf(loc->kernel_c, "(i%d)", j);
        }
        fprintf(loc->kernel_c, ")\n");
        fprintf(loc->kernel_c, "#define %s_%s_(",
                stmt_name, loc->array[i].name);
        for (j = 0; j < 1 + loc->array[i].dim; ++j) {
            if (j)
                fprintf(loc->kernel_c, ",");
            fprintf(loc->kernel_c, "i%d", j);
        }
        fprintf(loc->kernel_c, ") %s_%s[", lhs_array, loc->array[i].name);
        isl_pw_qpolynomial_print(lhs_i, loc->kernel_c, ISL_FORMAT_C);
        fprintf(loc->kernel_c, "] = %s_%s[", rhs_array, loc->array[i].name);
        isl_pw_qpolynomial_print(rhs_i, loc->kernel_c, ISL_FORMAT_C);
        fprintf(loc->kernel_c, "]\n");
        isl_pw_qpolynomial_free(lhs_i);
        isl_pw_qpolynomial_free(rhs_i);
    }
    isl_union_pw_qpolynomial_free(rhs);
    isl_union_pw_qpolynomial_free(lhs);
}

/* For each array A, print
 *
 *  #define reorganize_A(s,i) L2[rho(s+1,i)] = L[rho(s,i)]
 *
 * where s is the sequential loopon the GPU and i is the array index.
 * The host iterators h%d and the block ids b%d are treated as parameters.
 */
void print_reorganize_define(struct localizer_info *loc,
    __isl_keep isl_union_map *reorganize,
    __isl_keep isl_union_pw_qpolynomial *rho)
{
    isl_dim *dim;
    isl_map *next;
    isl_union_set *array;
    isl_union_map *id;
    isl_union_pw_qpolynomial *rho_next;

    dim = isl_union_map_get_dim(reorganize);
    next = compute_next(dim, loc->gpu_len);
    array = isl_union_pw_qpolynomial_domain(isl_union_pw_qpolynomial_copy(rho));
    array = isl_union_map_range(isl_union_set_unwrap(array));
    id = isl_union_set_identity(array);
    id = isl_union_map_product(isl_union_map_from_map(next), id);
    rho_next = isl_union_pw_qpolynomial_copy(rho);
    rho_next = isl_union_map_apply_union_pw_qpolynomial(id, rho_next);

    print_copy_define(loc, reorganize, "reorganize", "L2", rho_next, "L", rho);

    isl_union_pw_qpolynomial_free(rho_next);
}

static void print_swap_buffers_define(struct localizer_info *loc)
{
    int i;
    fprintf(loc->kernel_c, "#define swap_buffers() do { %s *t; ", loc->type);
    for (i = 0; i < loc->n_array; ++i) {
        fprintf(loc->kernel_c, "t = L_%s; ", loc->array[i].name);
        fprintf(loc->kernel_c, "L_%s = L2_%s; ",
                loc->array[i].name, loc->array[i].name);
        fprintf(loc->kernel_c, "L2_%s = t; ", loc->array[i].name);
    }
    fprintf(loc->kernel_c, "} while (0)\n");
}

static void print_defines(struct localizer_info *loc)
{
    fprintf(loc->kernel_c, "#define threadIdx_x ((int)threadIdx.x)\n");
    fprintf(loc->kernel_c, "#define threadIdx_y ((int)threadIdx.y)\n");
    fprintf(loc->kernel_c, "#define threadIdx_z ((int)threadIdx.z)\n");
    fprintf(loc->kernel_c, "#define sync1 __syncthreads\n");
    fprintf(loc->kernel_c, "#define sync2 __syncthreads\n");
    fprintf(loc->kernel_c, "#define sync3 __syncthreads\n");
}

/* Compute the number of shared memory locations needed per array
 * and per block and compute the maximal total number of shared memory
 * locations over the whole program run.
 *
 * For each array, we first compute the number of memory locations
 * needed for each iteration of the sequential loop in the block
 * and then compute an upper bound over all iterations of this loop.
 * The sum is kept in total_size and a bound on total_size
 * is stored in loc->max_shared_size.
 */
static void compute_shared_sizes(struct localizer_info *loc,
    __isl_keep isl_union_map *gpu_access, __isl_keep isl_set *block_domain)
{
    int i, j;
    char name[20];
    isl_dim *dim;
    isl_map *flatten_map;
    isl_pw_qpolynomial_fold *total_size;

    dim = isl_union_map_get_dim(gpu_access);
    dim = isl_dim_add(dim, isl_dim_in, loc->gpu_len - 1);
    dim = isl_dim_add(dim, isl_dim_out, 1);
    dim = isl_dim_wrap(dim);
    flatten_map = isl_set_flatten_map(isl_set_universe(dim));

    dim = isl_union_map_get_dim(gpu_access);
    dim = isl_dim_add(dim, isl_dim_set, loc->gpu_len - 1);
    total_size = isl_pw_qpolynomial_fold_zero(dim, isl_fold_max);

    for (i = 0; i < loc->n_array; ++i) {
        int offset;
        isl_map *map_i;
        isl_pw_qpolynomial *size;
        isl_dim *dim = isl_union_map_get_dim(gpu_access);
        isl_pw_qpolynomial_fold *shared_size;

        dim = isl_dim_add(dim, isl_dim_in, loc->gpu_len);
        dim = isl_dim_add(dim, isl_dim_out, loc->array[i].dim);
        dim = isl_dim_set_tuple_name(dim, isl_dim_out, loc->array[i].name);
        map_i = isl_union_map_extract_map(gpu_access, dim);
        size = isl_map_card(map_i);
        size = isl_map_apply_pw_qpolynomial(isl_map_copy(flatten_map), size);
        shared_size = isl_pw_qpolynomial_bound(size, isl_fold_max, NULL);
        for (j = 0; j < loc->first + 1; ++j) {
            snprintf(name, sizeof(name), "h%d", j);
            shared_size = isl_pw_qpolynomial_fold_set_dim_name(shared_size,
                                                isl_dim_set, j, name);
        }
        offset = loc->first + 1;
        for (j = 0; j < loc->tile_len - 1; ++j) {
            snprintf(name, sizeof(name), "b%d", j);
            shared_size = isl_pw_qpolynomial_fold_set_dim_name(shared_size,
                                                 isl_dim_set, offset + j, name);
        }
        loc->array[i].shared_size = isl_pw_qpolynomial_fold_copy(shared_size);
        loc->array[i].shared_size =
            isl_pw_qpolynomial_fold_gist(loc->array[i].shared_size,
                                         isl_set_copy(block_domain));
        total_size = isl_pw_qpolynomial_fold_add(total_size, shared_size);
    }

    loc->max_shared_size = isl_pw_qpolynomial_fold_bound(total_size, NULL);

    isl_map_free(flatten_map);
}

/* Concatenate host_domain with bounds.
 */
static __isl_give isl_set *compute_block_domain(struct localizer_info *loc,
    __isl_keep isl_set *host_domain, __isl_keep isl_set **bounds)
{
    isl_dim *dim;
    isl_map *proj;
    isl_set *block_domain;

    dim = isl_union_map_get_dim(loc->sched);
    proj = projection(dim, loc->first + 2, loc->first + 1);
    proj = isl_map_reverse(proj);
    block_domain = isl_set_copy(host_domain);
    block_domain = isl_set_apply(block_domain, proj);
    assert(loc->tile_len == 2);
    block_domain = isl_set_intersect(block_domain, isl_set_copy(bounds[0]));

    return block_domain;
}

/* This function is called for each leaf in the clast of the host code.
 * We first specialize the schedule to the site of the leaf, derive
 * some characteristics and then print the body of host code
 * and the associated kernel (through a call to print_kernel).
 */
void gpu_print_host_user(struct localizer_info *loc, struct clast_user_stmt *u)
{
    int i;
    isl_union_map *local_sched;
    isl_set *host_domain = NULL;
    isl_set *block_domain;
    isl_union_set *local_domain;
    isl_set *bounds[loc->tile_len - 1];
    isl_union_map *copy_in, *copy_out, *reorganize;
    isl_union_map *local_access, *host_access, *gpu_access;
    isl_union_map *lifted_host_access;
    isl_pw_qpolynomial *size;
    isl_union_pw_qpolynomial *sigma, *rho;

    host_domain = extract_entire_host_domain(u);

    size = isl_pw_qpolynomial_copy(loc->transfer_size);
    size = isl_pw_qpolynomial_gist(size, isl_set_copy(host_domain));

    local_sched = isl_union_map_intersect_range(isl_union_map_copy(loc->sched),
                isl_union_set_from_set(extend(isl_set_copy(host_domain),
                                              loc->tiled_len)));
    local_access = compute_scheduled_accesses(loc, local_sched);
    local_domain = isl_union_map_range(isl_union_map_copy(local_sched));

    for (i = 0; i < loc->tile_len - 1; ++i) {
        bounds[i] = bounds_on(local_domain, loc->tiled_len,
                        loc->first + 1 + i, host_domain, loc->host_len);
    }
    block_domain = compute_block_domain(loc, host_domain, bounds);

    gpu_access = project_domain(isl_union_map_copy(local_access),
                                    loc->tiled_len, loc->gpu_len);
    host_access = project_domain(local_access, loc->tiled_len, loc->host_len);
    lifted_host_access = extend_domain(isl_union_map_copy(host_access),
                                       loc->host_len, loc->gpu_len);

    compute_shared_sizes(loc, gpu_access, block_domain);

    copy_in = compute_copy_in(loc, isl_union_map_copy(gpu_access));
    copy_out = compute_copy_out(loc, isl_union_map_copy(gpu_access));
    reorganize = compute_reorganize(loc, isl_union_map_copy(gpu_access));

    rho = compute_ranking(gpu_access);
    sigma = compute_ranking(lifted_host_access);

    print_transfer_size_init(loc, size);

    print_cloog_macros(loc->kernel_c);
    print_statement_define(loc, local_sched, rho);

    print_copy_define(loc, copy_in, "copyin", "L", rho, "buffer", sigma);
    print_copy_define(loc, copy_out, "copyout", "buffer", sigma, "L", rho);
    print_reorganize_define(loc, reorganize, rho);

    print_swap_buffers_define(loc);
    print_defines(loc);

    bundle(loc, host_access, host_domain, 0);
    print_indent(loc->dst, loc->indent);
    fprintf(loc->dst,
        "cudaMemcpy(dev_buffer, transfer_buffer, transfer_size, "
        "cudaMemcpyHostToDevice);\n");

    print_kernel(loc, host_domain, block_domain, bounds,
                 copy_in, copy_out, reorganize);

    print_indent(loc->dst, loc->indent);
    fprintf(loc->dst,
        "cudaMemcpy(transfer_buffer, dev_buffer, transfer_size, "
        "cudaMemcpyDeviceToHost);\n");
    bundle(loc, host_access, host_domain, 1);

    for (i = 0; i < loc->n_array; ++i)
        isl_pw_qpolynomial_fold_free(loc->array[i].shared_size);
    isl_union_map_free(gpu_access);
    isl_union_map_free(host_access);
    isl_union_map_free(lifted_host_access);
    isl_set_free(bounds[0]);
    isl_set_free(host_domain);
    isl_set_free(block_domain);
    isl_union_set_free(local_domain);
    isl_union_map_free(local_sched);
    isl_union_pw_qpolynomial_free(sigma);
    isl_union_pw_qpolynomial_free(rho);
    isl_pw_qpolynomial_free(size);
    isl_pw_qpolynomial_fold_free(loc->max_shared_size);

    loc->kernel_id++;
}

static void print_bundle_array_access(FILE *dst, struct array_info *array)
{
    int j;

    fprintf(dst, "%s", array->name);
    for (j = 0; j < array->dim; ++j)
        fprintf(dst, "[i%d]", j);
}

static void print_bundle_defines(FILE *dst, struct array_info *array,
    int unbundle)
{
    int j;

    fprintf(dst, "#define %sbundle_%s(", unbundle ? "un" : "", array->name);
    for (j = 0; j < array->dim; ++j) {
        if (j)
            fprintf(dst, ",");
        fprintf(dst, "i%d", j);
    }
    fprintf(dst, ") ");
    if (unbundle) {
        print_bundle_array_access(dst, array);
        fprintf(dst, " = *p++");
    } else {
        fprintf(dst, "*p++ = ");
        print_bundle_array_access(dst, array);
    }
    fprintf(dst, "\n");
}

/* Print the code that will run on the host.
 * In particular, use CLooG to generate code for the outer "loc->host_len" loops
 * of the global schedule in "loc->sched".
 * The pretty printing of this code is handled by gpu_print_host_stmt,
 * which calls gpu_print_host_user for each kernel invocation location.
 */
static void print_host_code(struct localizer_info *loc)
{
    int i;
    CloogOptions *options;
    CloogDomain *cloog_context;
    CloogUnionDomain *ud;
    CloogInput *input;
    struct clast_stmt *stmt;
    char name[20];

    options = cloog_options_malloc(loc->state);
    options->language = LANGUAGE_C;
    options->otl = 0;
    options->strides = 1;
    options->stop = loc->host_len;
    options->f = loc->gpu_len;
    options->l = loc->tiled_len;;
    options->save_domains = 1;
    options->noscalars = 1;

    ud = cloog_union_domain_from_isl_union_map(isl_union_map_copy(loc->sched));
    for (i = 0; i < options->stop; ++i) {
        snprintf(name, sizeof(name), "h%d", i);
        ud = cloog_union_domain_set_name(ud, CLOOG_SCAT, i, name);
    }
    cloog_context = cloog_domain_from_isl_set(isl_set_copy(loc->prog->context));
    input = cloog_input_alloc(cloog_context, ud);

    stmt = cloog_clast_create_from_input(input, options);
    loc->dst = loc->host_c;
    fprintf(loc->dst, "{\n");
    print_cloog_macros(loc->dst);
    for (i = 0; i < loc->n_array; ++i) {
        print_bundle_defines(loc->dst, &loc->array[i], 0);
        print_bundle_defines(loc->dst, &loc->array[i], 1);
    }
    fprintf(loc->dst, "%s *p;\n", loc->type);
    fprintf(loc->dst, "%s transfer_buffer[", loc->type);
    isl_pw_qpolynomial_fold_print(loc->max_transfer_size, loc->dst, ISL_FORMAT_C);
    fprintf(loc->dst, "];\n");
    fprintf(loc->dst, "%s *dev_buffer;\n", loc->type);
    fprintf(loc->dst, "size_t transfer_size;\n");
    fprintf(loc->dst, "int ");
    for (i = 0; i < options->stop; ++i) {
        if (i)
            fprintf(loc->dst, ", ");
        fprintf(loc->dst, "h%d", i);
    }
    fprintf(loc->dst, ";\n");
    fprintf(loc->dst, "cudaMalloc(&dev_buffer, (");
    isl_pw_qpolynomial_fold_print(loc->max_transfer_size, loc->dst, ISL_FORMAT_C);
    fprintf(loc->dst, ") * sizeof(%s));\n", loc->type);
    loc->indent = 0;
    gpu_print_host_stmt(loc, stmt);
    fprintf(loc->dst, "cudaFree(dev_buffer);\n");
    fprintf(loc->dst, "}\n");

    cloog_clast_free(stmt);
    cloog_options_free(options);
}

static char *skip_spaces(char *s)
{
    while (isspace(*s))
        ++s;
    return s;
}

static int is_begin_scop(char *line)
{
    line = skip_spaces(line);
    if (*line != '#')
        return 0;
    line = skip_spaces(line + 1);
    if (strncmp(line, "pragma", sizeof("pragma") - 1))
        return 0;
    line = skip_spaces(line + sizeof("pragma") - 1);
    if (strncmp(line, "scop", sizeof("scop") - 1))
        return 0;
    return 1;
}

static int is_end_scop(char *line)
{
    line = skip_spaces(line);
    if (*line != '#')
        return 0;
    line = skip_spaces(line + 1);
    if (strncmp(line, "pragma", sizeof("pragma") - 1))
        return 0;
    line = skip_spaces(line + sizeof("pragma") - 1);
    if (strncmp(line, "endscop", sizeof("endscop") - 1))
        return 0;
    return 1;
}

/* Open the "input" file for reading and open the host .cu file
 * and the kernel .hu and .cu files for writing.
 * Add the necessary includes and copy all code from the input
 * file up to the openscop pragma to the host .cu file.
 */
static void open_files(struct localizer_info *loc, const char *input)
{
    char name[PATH_MAX];
    const char *base;
    const char *ext;
    int len;
    char line[1024];

    base = strrchr(input, '/');
    if (base)
        base++;
    else
        base = input;
    ext = strrchr(base, '.');
    len = ext ? ext - base : strlen(base);

    memcpy(name, base, len);
    strcpy(name + len, "_host.cu");
    loc->host_c = fopen(name, "w");

    strcpy(name + len, "_kernel.cu");
    loc->kernel_c = fopen(name, "w");

    strcpy(name + len, "_kernel.hu");
    loc->kernel_h = fopen(name, "w");
    fprintf(loc->host_c, "#include \"%s\"\n", name);
    fprintf(loc->kernel_c, "#include \"%s\"\n", name);
    fprintf(loc->kernel_h, "#include \"cuda.h\"\n\n");

    loc->input = fopen(input, "r");
    while (fgets(line, sizeof(line), loc->input)) {
        fprintf(loc->host_c, "%s", line);
        if (is_begin_scop(line))
            break;
    }
}

/* Copy all code starting at the endscop pragma from the input
 * file to the host .cu file and close all input and output files.
 */
static void close_files(struct localizer_info *loc)
{
    char line[1024];

    while (fgets(line, sizeof(line), loc->input)) {
        if (is_end_scop(line)) {
            fprintf(loc->host_c, "%s", line);
            break;
        }
    }
    while (fgets(line, sizeof(line), loc->input)) {
        fprintf(loc->host_c, "%s", line);
    }

    fclose(loc->input);
    fclose(loc->kernel_c);
    fclose(loc->kernel_h);
    fclose(loc->host_c);
}

static void select_tile_dimensions(struct localizer_info *loc,
    PlutoProg *prog)
{
    int last;

    getOutermostTilableBand(prog, &loc->first, &last);
    /* need at least two loops */
    assert(last > loc->first);
    /* only tile two loops for now */
    last = loc->first + 1;

    loc->tile_len = last - loc->first + 1;
    loc->len = prog->num_hyperplanes;

    loc->host_len = loc->first + 1;
    loc->gpu_len = loc->first + loc->tile_len + 1;
    loc->tiled_len = loc->len + loc->tile_len;
}

/* Extract the schedule computed by Pluto and apply skewing and tiling.
 */
static __isl_give isl_union_map *compute_global_schedule(
    struct localizer_info *loc, PlutoProg *prog, PlutoOptions *options)
{
    isl_dim *dim;
    isl_map *tiling, *front;
    isl_union_map *schedule;

    schedule = extract_schedule(prog);

    dim = isl_union_map_get_dim(schedule);

    tiling = tile(isl_dim_copy(dim), loc->len, loc->first, loc->tile_len,
                    options->tile_size);

    front = wavefront(isl_dim_copy(dim), loc->len, loc->first, loc->tile_len);
    tiling = isl_map_apply_range(front, tiling);

    front = wavefront(dim, loc->len + loc->tile_len, loc->first, loc->tile_len);
    tiling = isl_map_apply_range(tiling, front);

    schedule = isl_union_map_apply_range(schedule,
                                         isl_union_map_from_map(tiling));

    return schedule;
}

/* Compute the maximal total transfer size over all iterations of
 * the host loops, or at least an upper bound on this maximum.
 */
static void compute_max_transfer_size(struct localizer_info *loc)
{
    int i;

    assert(loc->n_array > 0);
    loc->transfer_size = isl_pw_qpolynomial_copy(loc->array[0].transfer_size);
    for (i = 1; i < loc->n_array; ++i)
        loc->transfer_size = isl_pw_qpolynomial_add(loc->transfer_size,
                    isl_pw_qpolynomial_copy(loc->array[i].transfer_size));
    loc->max_transfer_size = isl_pw_qpolynomial_bound(
            isl_pw_qpolynomial_copy(loc->transfer_size), isl_fold_max, NULL);
}

/* Replace the scop in the "input" file by code that equivalent code
 * that uses the GPU.  "prog" is assumed to correspond to this scop
 * and the statements are assumed to have been assigned a valid "trans"
 * by Pluto that ensures that all dependences are non-negative in the
 * coordinate directions.
 *
 * We first select a sequence of two loops to tile.
 * Then we skew these two loop such that the outer of the two
 * carries all dependences carried by the original two loops.
 * In particular, we apply a skew S := {[i,j] -> [i+j,j]}.
 * The two loops are then tiled with tile size options->tile_size.
 * Finally, the tile loops are skewed again to ensure that the outer
 * tile dimension carries all dependences.
 *
 * In summary, let P be the transformation computed by Pluto, then
 * the final schedule is
 *
 * (id flat_cross ((S flat_cross id) after T after S) flat_cross id) after P
 *
 * with T the tiling and id identity mappings of the appropriate dimensions.
 *
 * The first dimensions, up to and including the sequential loop over
 * the tiles are run on the host and are named h%d.
 * The parallel loop over the tiles is spread over the blocks and is
 * called b0.
 * The sequential loop inside the tiles is run on the GPU and is called s.
 * The parallel loop inside the tiles is wrapped over the threads.
 * That is, it is subjected to
 *
 *      [threadIdx_x] -> { [t] -> [t'] : t = 512 t' + threadIdx_x and
 *                         0 <= threadIdx_x < 512 }
 *
 * Any remaining loops inside the wrapped parallel loop are called l%d.
 *
 * Code is first generated for the host.  The loop structure is obtained
 * by letting CLooG generate code for the first dimensions of the global
 * schedule.  The bodies of the sequential loops over the tiles are then
 * set to copying data to the GPU, kernel invocation and copying data back
 * from the GPU.
 * For each kernel invocation, the global schedule is specialized for
 * the invocation site.  Extra statements are added to the schedule
 * for copying code in and out of shared memory and for synchronization
 * and code is generated for each kernel by calling CLooG.
 *
 * The function frees "prog" and "options".
 */
int gpuloc(PlutoProg *prog, PlutoOptions *options, const char *input)
{
    struct localizer_info loc;

    loc.state = cloog_isl_state_malloc(prog->ctx);
    loc.type = options->type ? options->type : "float";
    loc.prog = prog;

    open_files(&loc, input);
    select_tile_dimensions(&loc, prog);
    loc.sched = compute_global_schedule(&loc, prog, options);
    collect_array_info(&loc);
    compute_max_transfer_size(&loc);

    loc.kernel_id = 0;
    print_host_code(&loc);

    isl_pw_qpolynomial_free(loc.transfer_size);
    isl_pw_qpolynomial_fold_free(loc.max_transfer_size);
    free_array_info(&loc);
    isl_union_map_free(loc.sched);
    cloog_state_free(loc.state);
    pluto_options_free(options);
    pluto_prog_free(prog);

    close_files(&loc);

    return 0;
}
