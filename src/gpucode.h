#ifndef _GPUCODE_H
#define _GPUCODE_H

#include <cloog/isl/cloog.h>

#include "gpuloc.h"
#include "pluto.h"

void print_cloog_macros(FILE *dst);
void print_indent(FILE *dst, int indent);
void gpu_print_host_stmt(struct localizer_info *loc, struct clast_stmt *s);

#endif
