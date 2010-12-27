#include <assert.h>
#include <ctype.h>

#include <isl/constraint.h>

#include "schedule.h"
#include "program.h"

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
__isl_give isl_union_map *extract_schedule(PlutoProg *prog)
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

/* Construct a map that maps a domain of dimensionality "len"
 * to another domain of the same dimensionality such that
 * coordinate "first" of the range is equal to the sum of the "wave_len"
 * coordinates starting at "first" in the domain.
 * The remaining coordinates in the range are equal to the corresponding ones
 * in the domain.
 * "dim" prescribes the parameters.
 */
__isl_give isl_map *wavefront(__isl_take isl_dim *dim, int len,
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

/* Construct a map from a len-dimensional domain to
 * a (len-n)-dimensional domain that projects out the n coordinates
 * starting at first.
 * "dim" prescribes the parameters.
 */
__isl_give isl_map *project_out(__isl_take isl_dim *dim,
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
__isl_give isl_map *projection(__isl_take isl_dim *dim,
    int src_len, int dst_len)
{
    return project_out(dim, src_len, dst_len, src_len - dst_len);
}

/* Extend "set" with unconstrained coordinates to a total length of "dst_len".
 */
__isl_give isl_set *extend(__isl_take isl_set *set, int dst_len)
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

/* Extract the access in stmt->text starting at position identifier
 * and of length identifier_len, and rewrite the index expression A[i]
 * by calling "print_access".
 * 
 * The access in C notation is first copied to "buffer" (which
 * has been allocated by the caller and should be of sufficient size)
 * and slightly modified to a map in isl notation.
 * This string is then parsed by isl and passed on to "print_access".
 */
static void extract_access(isl_ctx *ctx, PlutoProg *prog, Stmt *stmt,
    char *buffer, int identifier, int identifier_len,
    int access, int access_len, const char *name,
    void (*print_access)(__isl_take isl_map *access, Stmt *stmt, void *user),
    void *user)
{
    int i;
    int pos = 0;
    isl_map *map;

    pos += sprintf(buffer, "[");
    for (i = 0; i < prog->npar; ++i) {
        if (i)
            pos += sprintf(buffer + pos, ",");
        pos += sprintf(buffer + pos, "%s", prog->params[i]);
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

    print_access(map, stmt, user);
}

/* Print stmt->text to out, calling "print_access" to rewrite each access.
 * "name" is the name of the statement and is of the form S_%d.
 */
void print_stmt_body(FILE *out, isl_ctx *ctx, PlutoProg *prog, 
    Stmt *stmt, const char *name,
    void (*print_access)(__isl_take isl_map *access, Stmt *stmt, void *user),
    void *user)
{
    int i, j;
    size_t text_len = strlen(stmt->text);
    size_t len = 50;
    char *buffer;
    int printed = 0;
    int identifier = -1;
    int end = -1;

    for (i = 0; i < prog->npar; ++i)
        len += strlen(prog->params[i]);
    for (i = 0; i < stmt->dim; ++i)
        len += strlen(stmt->iterators[i]);
    buffer = isl_alloc_array(ctx, char, len);
    assert(buffer);

    for (i = 0; i < text_len; ++i) {
        if (identifier < 0 && isalpha(stmt->text[i])) {
            fwrite(stmt->text + printed, 1, i - printed, out);
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
            extract_access(ctx, prog, stmt, buffer,
			    identifier, end - identifier,
                            i + 1, j - i - 1, name, print_access, user);
            i = j;
            end = identifier = -1;
            printed = i + 1;
        } else {
            end = identifier = -1;
        }
    }
    fwrite(stmt->text + printed, 1, text_len - printed, out);

    free(buffer);
}
