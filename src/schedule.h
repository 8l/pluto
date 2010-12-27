#ifndef _SCHEDULE_H
#define _SCHEDULE_H

#include "pluto.h"

__isl_give isl_union_map *extract_schedule(PlutoProg *prog);
__isl_give isl_map *wavefront(__isl_take isl_dim *dim, int len,
        int first, int wave_len);
__isl_give isl_map *project_out(__isl_take isl_dim *dim,
	int len, int first, int n);
__isl_give isl_map *projection(__isl_take isl_dim *dim,
	int src_len, int dst_len);
__isl_give isl_set *extend(__isl_take isl_set *set, int dst_len);

void print_stmt_body(FILE *out, isl_ctx *ctx, PlutoProg *prog, 
    Stmt *stmt, const char *name,
    void (*print_access)(__isl_take isl_map *access, Stmt *stmt, void *user),
    void *user);

#endif
