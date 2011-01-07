/*
 * PLUTO: An automatic parallelizer and locality optimizer
 * 
 * Copyright (C) 2007--2008 Uday Kumar Bondhugula
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * A copy of the GNU General Public Licence can be found in the file
 * `LICENSE' in the top-level directory of this distribution. 
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <getopt.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <isl/options.h>

#include "pluto.h"

#include "clan/clan.h"
#include "candl/candl.h"

#include "math_support.h"
#include "post_transform.h"
#include "ddg.h"
#include "program.h"
#include "gpuloc.h"
#include "cuda.h"

PlutoOptions *options;

static struct isl_arg_choice fuse_choice[] = {
    {"no",    NO_FUSE},
    {"max",   MAXIMAL_FUSE},
    {"smart", SMART_FUSE},
    {0}
};

static struct isl_arg_choice dep_choice[] = {
    {"clan",  DEP_CLAN},
    {"isl",   DEP_ISL},
    {0}
};

static void print_version()
{
    printf("PLUTO 0.6.0 - An automatic parallelizer and locality optimizer\n\
Copyright (C) 2007--2008  Uday Kumar Bondhugula\n\
This is free software; see the source for copying conditions.  There is NO\n\
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n");
}

struct isl_arg options_arg[] = {
ISL_ARG_CHILD(PlutoOptions, isl, "isl", isl_options_arg, "isl options")
ISL_ARG_BOOL(PlutoOptions, tile, 0, "tile", 0, "Tile for locality")
ISL_ARG_BOOL(PlutoOptions, parallel, 0, "parallel", 0,
    "Automatically parallelize using OpenMP pragmas")
ISL_ARG_ALIAS("parallelize")
ISL_ARG_BOOL(PlutoOptions, l2tile, 0, "l2tile", 0,
    "Tile a second time (typically for L2 cache)")
ISL_ARG_BOOL(PlutoOptions, multipipe, 0, "multipipe", 0,
    "Extract two degrees of pipelined parallelism if possible; "
    "by default one degree is extracted (if it exists)")
ISL_ARG_BOOL(PlutoOptions, rar, 0, "rar", 0, "Consider RAR dependences too")
ISL_ARG_BOOL(PlutoOptions, unroll, 0, "unroll", 0, "Unroll-jam")
ISL_ARG_INT(PlutoOptions, ufactor, 0, "ufactor", "factor", 8,
    "Unroll-jam factor")
ISL_ARG_BOOL(PlutoOptions, prevector, 0, "prevector", 1,
    "Make code amenable to compiler auto-vectorization (with ICC)")
/* Default context is no context */
ISL_ARG_INT(PlutoOptions, context, 0, "context", "context", -1,
    "Lower bound on parameters")
ISL_ARG_CHOICE(PlutoOptions, dep, 0, "dep", dep_choice, DEP_CLAN,
    "Dependence tester")
ISL_ARG_BOOL(PlutoOptions, lastwriter, 0, "lastwriter", 0,
    "Work with refined dependences "
    "(last conflicting access is computed for RAW/WAW)")
ISL_ARG_BOOL(PlutoOptions, bee, 0, "bee", 0, "Generate pragmas for Bee+Cl@k")
ISL_ARG_PHANTOM_BOOL('i', "indent", NULL, "Indent generated code")
ISL_ARG_BOOL(PlutoOptions, silent, 'q', "silent", 0,
    "Silent mode; no output as long as everything goes fine")
ISL_ARG_CHOICE(PlutoOptions, fuse, 0, "fuse", fuse_choice, SMART_FUSE,
    "Fusion heuristic; no: do not fuse across SCCs of data dependence graph; "
    "max: maximal fusion; smart: heuristic (in between no and max)")
ISL_ARG_BOOL(PlutoOptions, debug, 0, "debug", 0, "Verbose output")
ISL_ARG_BOOL(PlutoOptions, moredebug, 0, "moredebug", 0, "More verbose output")
ISL_ARG_BOOL(PlutoOptions, gpuloc, 0, "gpuloc", 0, "Localize for GPU")
ISL_ARG_BOOL(PlutoOptions, cuda, 0, "cuda", 0, "Generate CUDA code")
ISL_ARG_BOOL(PlutoOptions, cuda_scale_tile_loops, 0,
	"cuda-scale-tile-loops", 1, NULL)
ISL_ARG_BOOL(PlutoOptions, cuda_wrap, 0, "cuda-wrap", 1, NULL)
ISL_ARG_STR(PlutoOptions, type, 't', "type", "type", "float",
    "Element type of arrays")
ISL_ARG_INT(PlutoOptions, tile_size, 'S', "tile-size", "size",
    DEFAULT_L1_TILE_SIZE, NULL)
/* Override for first and last levels to tile */
ISL_ARG_INT_F(PlutoOptions, ft, 0, "ft", NULL, -1, NULL, ISL_ARG_HIDDEN)
ISL_ARG_INT_F(PlutoOptions, lt, 0, "lt", NULL, -1, NULL, ISL_ARG_HIDDEN)
/* Override for first and last cloog options */
ISL_ARG_INT_F(PlutoOptions, cloogf, 0, "cloogf", NULL, -1, NULL, ISL_ARG_HIDDEN)
ISL_ARG_INT_F(PlutoOptions, cloogl, 0, "cloogl", NULL, -1, NULL, ISL_ARG_HIDDEN)
/* Experimental */
ISL_ARG_BOOL_F(PlutoOptions, polyunroll, 0, "polyunroll", 0, NULL,
    ISL_ARG_HIDDEN)
ISL_ARG_BOOL_F(PlutoOptions, bound, 0, "bound", 1, NULL, ISL_ARG_HIDDEN)
ISL_ARG_BOOL_F(PlutoOptions, scalpriv, 0, "scalpriv", 0, NULL, ISL_ARG_HIDDEN)
ISL_ARG_VERSION(&print_version)
ISL_ARG_END
};

ISL_ARG_DEF(options, PlutoOptions, options_arg)

struct plutoArg {
    PlutoOptions *options;
    char *srcFileName;
};
typedef struct plutoArg PlutoArg;

struct isl_arg arg_arg[] = {
ISL_ARG_CHILD(PlutoArg, options, NULL, options_arg, NULL)
ISL_ARG_ARG(PlutoArg, srcFileName, "input", NULL)
ISL_ARG_FOOTER("To report bugs, please send an email to "
    "<pluto-development@googlegroups.com>")
ISL_ARG_END
};

ISL_ARG_DEF(arg, PlutoArg, arg_arg)

int main(int argc, char *argv[])
{
    int i;

    FILE *src_fp;

    char outFileName[256] = "";

    char cloogFileName[256];
    FILE *cloogfp, *outfp;

    PlutoArg *arg;

    arg = arg_new_with_defaults();
    argc = arg_parse(arg, argc, argv, ISL_ARG_ALL);
    options = arg->options;

    src_fp  = fopen(arg->srcFileName, "r");

    if (!src_fp)   {
        fprintf(stderr, "pluto: error opening source file: '%s'\n",
                        arg->srcFileName);
        return 5;
    }

    /* Extract polyhedral representation from input program */
    scoplib_scop_p scop;

    clan_options_p clanOptions = clan_options_malloc();

    scop = clan_scop_extract(src_fp, clanOptions);

    if (!scop || !scop->statement)   {
        fprintf(stderr, "Error extracting polyhedra from source file: \'%s'\n",
                arg->srcFileName);
        return 1;
    }

    /* IF_DEBUG(clan_scop_print_dot_scop(stdout, scop, clanOptions)); */

    /* Convert clan scop to Pluto program */
    PlutoProg *prog = scop_to_pluto_prog(scop, options);

    clan_options_free(clanOptions);

    /* Backup irregular program portion in .scop. */
    char* irroption = scoplib_scop_tag_content(scop, "<irregular>",
                                            "</irregular>");

    scoplib_scop_free(scop);

    IF_DEBUG2(deps_print(stdout, prog->deps, prog->ndeps));
    IF_DEBUG2(stmts_print(stdout, prog->stmts, prog->nstmts));

    /* Create the data dependence graph */
    prog->ddg = ddg_create(prog);
    ddg_compute_scc(prog);

    int dim_sum=0;
    for (i=0; i<prog->nstmts; i++) {
        dim_sum += prog->stmts[i].dim;
    }

    /* Make options consistent */
    if (options->multipipe == 1 && options->parallel == 0)    {
        fprintf(stdout, "Warning: multipipe needs parallel to be on; turning on parallel\n");
        options->parallel = 1;
    }

    /* Disable pre-vectorization if tile is not on */
    if (options->tile == 0 && options->prevector == 1) {
        /* If code will not be tiled, pre-vectorization does not make
         * sense */
        if (!options->silent)   {
            fprintf(stdout, "[Pluto] Warning: pre-vectorization does not fit (--tile is off)\n");
        }
        options->prevector = 0;
    }

    if (!options->silent)   {
        fprintf(stdout, "[Pluto] Number of statements: %d\n", prog->nstmts);
        fprintf(stdout, "[Pluto] Total number of loops: %d\n", dim_sum);
        fprintf(stdout, "[Pluto] Number of deps: %d\n", prog->ndeps);
        fprintf(stdout, "[Pluto] Maximum domain dimensionality: %d\n", prog->nvar);
        fprintf(stdout, "[Pluto] Number of parameters: %d\n", prog->npar);
    }

    /* Auto transformation */
    pluto_auto_transform(prog);

    if (!options->silent)   {
        fprintf(stdout, "[Pluto] Affine transformations [<iter coeff's> <const>]\n\n");
    }

    Stmt *stmts = prog->stmts;
    int nstmts = prog->nstmts;

    /* Print out the transformations */
    if (!options->silent)   {
        for (i=0; i<nstmts; i++) {
            fprintf(stdout, "T(S%d): ", i+1);
            int level;
            printf("(");
            for (level=0; level<prog->num_hyperplanes; level++) {
                if (level > 0) printf(", ");
                pretty_print_affine_function(stdout, &stmts[i], level);
            }
            printf(")\n");

            pluto_matrix_print(stdout, stmts[i].trans);
        }

        print_hyperplane_properties(prog->hProps, prog->num_hyperplanes);
    }

    if (options->cuda) {
        int r = cuda(prog, options, arg->srcFileName);
        arg_free(arg);
        return r;
    }

    if (options->gpuloc) {
        int r = gpuloc(prog, options, arg->srcFileName);
        arg_free(arg);
        return r;
    }

    if (options->tile)   {
        pluto_tile(prog);
    }

    if (options->parallel)   {
        int outermostBandStart, outermostBandEnd;
        getOutermostTilableBand(prog, &outermostBandStart, &outermostBandEnd);

        /* Obtain pipelined parallelization by skewing the tile space */
        bool retval = create_tile_schedule(prog, outermostBandStart, outermostBandEnd);

        /* Even if the user hasn't supplied --tile and there is only pipelined
         * parallelism, we will warn the user, but anyway do fine-grained 
         * parallelization
         */
        if (retval && options->tile == 0)   {
            printf("WARNING: --tile is not used and there is pipelined parallelism\n");
            printf("\t This leads to finer grained parallelism; add --tile to the list\n");
            printf("\t of cmd-line options for a better coarse-grained parallelized code.\n");
        }
    }

    if (options->prevector) {
        pre_vectorize(prog);
    }else{
        /* Create an empty .vectorize file */
        fopen(".vectorize", "w");
    }

    if (options->tile && !options->silent)  {
        fprintf(stdout, "[Pluto] After tiling:\n");
        print_hyperplane_properties(prog->hProps, prog->num_hyperplanes);
    }

    if (options->parallel)  {
        /* Generate meta info for insertion of OpenMP pragmas */
        generate_openmp_pragmas(prog);
    }


    if (options->unroll || options->polyunroll)    {
        /* Will generate a .unroll file */
        /* plann needs a .params */
        FILE *paramsFP = fopen(".params", "w");
        if (paramsFP)   {
            int i;
            for (i=0; i<prog->npar; i++)  {
                fprintf(paramsFP, "%s\n", prog->params[i]);
            }
            fclose(paramsFP);
        }
        detect_unrollable_loops(prog);
    }else{
        /* Create an empty .unroll file */
        fopen(".unroll", "w");
    }

    if (options->polyunroll)    {
        /* Experimental */
        for (i=0; i<prog->num_hyperplanes; i++)   {
            if (prog->hProps[i].unroll)  {
                unroll_phis(prog, i, options->ufactor);
            }
        }
    }

    /* The .cloog file name */
    strcpy(cloogFileName, arg->srcFileName);
    cloogFileName[strlen(arg->srcFileName)-2] = '\0';

    if (options->parallel && options->multipipe)   {
        strcat(cloogFileName, ".par2d.cloog");
    }else if (options->parallel)   {
        strcat(cloogFileName, ".par.cloog");
    }else if (options->tile)  {
        strcat(cloogFileName, ".tiled.cloog");
    }else{
        strcat(cloogFileName, ".opt.cloog");
    }

    cloogfp = fopen(cloogFileName, "w+");

    /* Remove .c extension and append a new one */
    strcpy(outFileName, arg->srcFileName);
    outFileName[strlen(arg->srcFileName)-2] = '\0';
    strcat(outFileName, ".pluto.c");

    outfp = fopen(outFileName, "w");

    if (!cloogfp)   {
        fprintf(stderr, "Can't open .cloog file: %s\n", cloogFileName);
        return 2;
    }


    /* Generate the .cloog file */
    IF_DEBUG(printf("[Pluto] Generating Cloog file\n"));
    print_cloog_file(cloogfp, prog);
    /* Add the <irregular> tag from clan, if any */
    if (irroption != NULL) {
        fprintf(cloogfp, "<irregular>\n%s\n</irregular>\n\n", irroption);
        free(irroption);
    }
    rewind(cloogfp);

    if (!outfp) {
        fprintf(stderr, "Can't open file %s for writing\n", outFileName);
        return 1;
    }

    /* Generate code using Cloog and add necessary stuff before/after code */
    pluto_codegen(cloogfp, outfp, prog);

    fclose(cloogfp);

    arg_free(arg);

    pluto_prog_free(prog);

    return 0;
}
