/*
 * PLUTO: An automatic parallelizer and locality optimizer
 * 
 * Copyright (C) 2007-2012 Uday Bondhugula
 *
 * This file is part of Pluto.
 *
 * Pluto is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
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
#include <libgen.h>

#include <unistd.h>
#include <sys/time.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include "osl/scop.h"
#include "osl/generic.h"
#include "osl/extensions/irregular.h"

#include "pluto.h"
#include "transforms.h"
#include "math_support.h"
#include "post_transform.h"
#include "program.h"
#include "version.h"

#include "clan/clan.h"
#include "candl/candl.h"
#include "candl/scop.h"

#include "pet.h"

PlutoOptions *options;

void usage_message(void)
{
    fprintf(stdout, "Usage: polycc <input.c> [options] [-o output]\n");
    fprintf(stdout, "\nOptions:\n");
    fprintf(stdout, "       --pet                     Use libpet for polyhedral extraction instead of clan [default - clan]\n");
    fprintf(stdout, "\n");
    fprintf(stdout, "       --tile                    Tile for locality\n");
    fprintf(stdout, "       --intratileopt            Optimize intra-tile execution order for locality\n");
    fprintf(stdout, "       --l2tile                  Tile a second time (typically for L2 cache) - disabled by default \n");
    fprintf(stdout, "       --parallel                Automatically parallelize (generate OpenMP pragmas)\n");
    fprintf(stdout, "    or --parallelize\n");
    fprintf(stdout, "       --lbtile                  Enables full-dimensional concurrent start\n");
    fprintf(stdout, "    or --diamond-tile\n");
    fprintf(stdout, "       --partlbtile              Enables one-dimensional concurrent start\n");
    fprintf(stdout, "       --[no]prevector           Transform for and mark loops for (icc) vectorization (enabled by default)\n");
    fprintf(stdout, "       --innerpar                Choose pure inner parallelism over pipelined/wavefront parallelism\n");
    fprintf(stdout, "       --multipipe               Extract two or more degrees of pipelined parallelism if possible;\n");
    fprintf(stdout, "                                     by default one degree is extracted (if it exists)\n");
    fprintf(stdout, "       --variables_not_global    Variables not declared globally (if so, macros provide variable declarations)\n");
    fprintf(stdout, "\n   Runtime               Options related to compilation for a runtime\n");
    fprintf(stdout, "       --dynschedule             Dynamically schedule tasks on processors using Synthesized Runtime Interface\n");
    fprintf(stdout, "                                     (for shared and distributed memory)\n");
    fprintf(stdout, "       --dynschedule_graph       Dynamically schedule tasks on processors using Intel TBB Flow Graph\n");
    fprintf(stdout, "                                     (only for shared-memory)\n");
    fprintf(stdout, "       --dataflow                Alias to --dynschedule\n");
    fprintf(stdout, "\n   Architecture          Options related to compilation for a specific architecture\n");
#ifdef PLUTO_OPENCL
    fprintf(stdout, "       --opencl                  Generate OpenCL code \n");
#endif
    fprintf(stdout, "       --distmem                 Parallelize for distributed-memory clusters (generate MPI)\n");
    fprintf(stdout, "       --mpiomp                  Parallelize for shared-memory along with distributed-memory (generate MPI+OpenMP)\n");
    fprintf(stdout, "       --data_dist               Performs data tiling and with distmem option distributes the data across multiple compute nodes at the granularity of data tiles\n");
    fprintf(stdout, "       --data_tile_opt           Trys to hoist mod and divide operations out of innermost loop, to be used with data_dist option\n");
    fprintf(stdout, "\n   Analysis              Options related to analyzing generated code (for Runtime or Architecture)\n");
    fprintf(stdout, "       --timereport              Generate code to report communication volume (for distributed-memory only)\n");
    fprintf(stdout, "                                     and analysis of time (for distributed-memory or runtime)\n");
    fprintf(stdout, "\n   Communication code    Options related to communication code generation for distributed-memory\n");
    fprintf(stdout, "       --commopt                 Generate communication code using Flow-Out (FO) scheme (enabled by default)\n");
    fprintf(stdout, "       --commopt_foifi           Generate communication code using Flow-Out Intersection Flow-In (FOIFI) scheme\n");
    fprintf(stdout, "       --commopt_fop             Generate communication code using Flow-Out Partitioning (FOP) scheme (multicast pack by default)\n");
    fprintf(stdout, "       --fop_unicast_runtime     FOP: generate code to choose between unicast and multicast pack at runtime\n\n");
    fprintf(stdout, "       --rar                  Consider RAR dependences too (disabled by default)\n");
    fprintf(stdout, "       --[no]unroll           Unroll-jam (disabled by default)\n");
    fprintf(stdout, "       --ufactor=<factor>     Unroll-jam factor (default is 8)\n");
    fprintf(stdout, "       --[no]prevector        Make code amenable to compiler auto-vectorization (with ICC) - enabled by default\n");
    fprintf(stdout, "       --context=<context>    Parameters are at least as much as <context>\n");
    fprintf(stdout, "       --forceparallel=<depth>  Depth (1-indexed) to force parallel\n");
    fprintf(stdout, "       --[no]isldep              Use ISL-based dependence tester (enabled by default)\n");
    fprintf(stdout, "       --islsolve [default]      Use ISL as ILP solver (default)\n");
    fprintf(stdout, "       --pipsolve                Use PIP as ILP solver\n");
    fprintf(stdout, "       --readscop             Read input from a .scop file\n");
    fprintf(stdout, "       --[no]lastwriter          Work with refined dependences (last conflicting access is computed for RAW/WAW)\n");
    fprintf(stdout, "                                 (enabled by default with --distmem; disabled otherwise)\n");
    fprintf(stdout, "       --bee                     Generate pragmas for Bee+Cl@k\n\n");
    fprintf(stdout, "       --indent  | -i            Indent generated code (disabled by default)\n");
    fprintf(stdout, "       --silent  | -q            Silent mode; no output as long as everything goes fine (disabled by default)\n");
    fprintf(stdout, "       --help    | -h            Print this help menu\n");
    fprintf(stdout, "       --version | -v            Display version number\n");
    fprintf(stdout, "\n   Fusion                Options to control fusion heuristic\n");
    fprintf(stdout, "       --nofuse                  Do not fuse across SCCs of data dependence graph\n");
    fprintf(stdout, "       --maxfuse                 Maximal fusion\n");
    fprintf(stdout, "       --smartfuse [default]     Heuristic (in between nofuse and maxfuse)\n");
    fprintf(stdout, "\n   Index Set Splitting        \n");
    fprintf(stdout, "       --iss                  \n");
    fprintf(stdout, "\n   Code generation       Options to control Cloog code generation\n");
    fprintf(stdout, "       --nocloogbacktrack        Do not call Cloog with backtrack (default - backtrack)\n");
    fprintf(stdout, "       --cloogsh                 Ask Cloog to use simple convex hull (default - off)\n");
    fprintf(stdout, "       --codegen-context=<value> Parameters are at least as much as <value>\n");
    fprintf(stdout, "\n   Debugging\n");
    fprintf(stdout, "       --debug                   Verbose/debug output\n");
    fprintf(stdout, "       --moredebug               More verbose/debug output\n");
    fprintf(stdout, "\nTo report bugs, please email <pluto-development@googlegroups.com>\n\n");
}

static double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

int main(int argc, char *argv[])
{
    int i;

    double t_start, t_c, t_d, t_t, t_all, t_start_all;

    t_c = 0.0;

    t_start_all = rtclock();

    FILE *src_fp;

    int option;
    int option_index = 0;

    int nolastwriter = 0;

    char *srcFileName;

    FILE *cloogfp, *outfp, *dynschedfp, *sigmafp, *headerfp;

    struct pet_scop *pscop;

    dynschedfp = NULL;
    sigmafp = NULL;
    headerfp = NULL;

    if (argc <= 1)  {
        usage_message();
        return 1;
    }

    options = pluto_options_alloc();

    const struct option pluto_options[] =
    {
        {"tile", no_argument, &options->tile, 1},
        {"notile", no_argument, &options->tile, 0},
        {"intratileopt", no_argument, &options->intratileopt, 1},
        {"nointratileopt", no_argument, &options->intratileopt, 0},
        {"lbtile", no_argument, &options->lbtile, 1},
        {"pet", no_argument, &options->pet, 1},
        {"diamond-tile", no_argument, &options->lbtile, 1},
        {"part-diamond-tile", no_argument, &options->partlbtile, 1},
        {"partlbtile", no_argument, &options->partlbtile, 1},
        {"dynschedule", no_argument, &options->dynschedule, 1},
        {"dataflow", no_argument, &options->dynschedule, 1},
        {"dynschedule_graph", no_argument, &options->dynschedule_graph, 1},
        {"dynschedule_graph_old", no_argument, &options->dynschedule_graph_old, 1},
        {"dyn_trans_deps_tasks", no_argument, &options->dyn_trans_deps_tasks, 1},
        {"debug", no_argument, &options->debug, 1},
        {"moredebug", no_argument, &options->moredebug, 1},
        {"rar", no_argument, &options->rar, 1},
        {"identity", no_argument, &options->identity, 1},
        {"nofuse", no_argument, &options->fuse, NO_FUSE},
        {"maxfuse", no_argument, &options->fuse, MAXIMAL_FUSE},
        {"smartfuse", no_argument, &options->fuse, SMART_FUSE},
        {"parallel", no_argument, &options->parallel, 1},
        {"parallelize", no_argument, &options->parallel, 1},
        {"innerpar", no_argument, &options->innerpar, 1},
        {"iss", no_argument, &options->iss, 1},
        {"distmem", no_argument, &options->distmem, 1},
        {"no_multi_level_distribution", no_argument, &options->multi_level_distribution, 0},
        {"multi_level_distribution", no_argument, &options->multi_level_distribution, 1},
        {"commopt", no_argument, &options->commopt, 1},
        {"commopt_fop", no_argument, &options->commopt_fop, 1},
        {"fop_unicast_runtime", no_argument, &options->fop_unicast_runtime, 1},
        {"commopt_foifi", no_argument, &options->commopt_foifi, 1},
        {"nocommopt", no_argument, &options->commopt, 0},
        {"timereport", no_argument, &options->timereport, 1},
        {"variables_not_global", no_argument, &options->variables_not_global, 1},
        {"mpiomp", no_argument, &options->mpiomp, 1},
        {"blockcyclic", no_argument, &options->blockcyclic, 1},
        {"unroll", no_argument, &options->unroll, 1},
        {"nounroll", no_argument, &options->unroll, 0},
        {"polyunroll", no_argument, &options->polyunroll, 1},
        {"bee", no_argument, &options->bee, 1},
        {"ufactor", required_argument, 0, 'u'},
        {"prevector", no_argument, &options->prevector, 1},
        {"noprevector", no_argument, &options->prevector, 0},
        {"codegen-context", required_argument, 0, 'c'},
        {"coeff-bound", required_argument, 0, 'C'},
        {"cloogf", required_argument, 0, 'F'},
        {"cloogl", required_argument, 0, 'L'},
        {"cloogsh", no_argument, &options->cloogsh, 1},
        {"nocloogbacktrack", no_argument, &options->cloogbacktrack, 0},
        {"cyclesize", required_argument, 0, 'S'},
        {"forceparallel", required_argument, 0, 'p'},
        {"ft", required_argument, 0, 'f'},
        {"lt", required_argument, 0, 'l'},
        {"multipipe", no_argument, &options->multipipe, 1},
        {"l2tile", no_argument, &options->l2tile, 1},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        {"indent", no_argument, 0, 'i'},
        {"silent", no_argument, &options->silent, 1},
        {"lastwriter", no_argument, &options->lastwriter, 1},
        {"nolastwriter", no_argument, &nolastwriter, 1},
        {"nodepbound", no_argument, &options->nodepbound, 1},
        {"scalpriv", no_argument, &options->scalpriv, 1},
        {"isldep", no_argument, &options->isldep, 1},
        {"candldep", no_argument, &options->candldep, 1},
        {"isldepaccesswise", no_argument, &options->isldepaccesswise, 1},
        {"isldepstmtwise", no_argument, &options->isldepaccesswise, 0},
        {"noisldepcoalesce", no_argument, &options->isldepcoalesce, 0},
        {"readscop", no_argument, &options->readscop, 1},
        {"pipsolve", no_argument, &options->pipsolve, 1},
        {"islsolve", no_argument, &options->islsolve, 1},
        {"fusesends", no_argument, &options->fusesends, 1},
        {"data_dist", no_argument, &options->data_dist, 1},
        {"verify_output", no_argument, &options->verify_output, 1},
        {"data_tile_opt", no_argument, &options->data_tile_opt, 1},
        {"nodata_tile_opt", no_argument, &options->data_tile_opt,0},
        {"identity_data_dist", no_argument, &options->identity_data_dist,1},
        {"global_opt", no_argument, &options->global_opt,1},
        {"noglobal_opt", no_argument, &options->global_opt,0},
        {"compute_pi", no_argument, &options->compute_pi,1},
        {"donot_compute_pi", no_argument, &options->compute_pi,0},
        {"num_tiles_per_dim", required_argument, 0, 'T'},
        {"num_parts", required_argument, 0, 'U'},
        {"time", no_argument, &options->time, 1},
        {0, 0, 0, 0}
    };


    /* Read command-line options */
    while (1) {
        option = getopt_long(argc, argv, "bhiqvf:l:F:L:c:o:", pluto_options,
                &option_index);

        if (option == -1)   {
            break;
        }

        switch (option) {
            case 0:
                break;
            case 'F':
                options->cloogf = atoi(optarg);
                break;
            case 'L':
                options->cloogl = atoi(optarg);
                break;
            case 'S':
                options->cyclesize = atoi(optarg);
                break;
            case 'T':
                options->num_tiles_per_dim = atoi(optarg);
                break;
            case 'U':
                options->num_inital_partitions = atoi(optarg);
                break;
            case 'b':
                options->bee = 1;
                break;
            case 'c':
                options->codegen_context = atoi(optarg);
                break;
            case 'C':
                options->coeff_bound = atoi(optarg);
                if (options->coeff_bound <= 0) {
                    printf("ERROR: coeff-bound should be at least 1\n");
                    return 2;
                }
                break;
            case 'f':
                options->ft = atoi(optarg);
                break;
            case 'g':
                break;
            case 'h':
                usage_message();
                return 2;
            case 'i':
                /* Handled in polycc */
                break;
            case 'l':
                options->lt = atoi(optarg);
                break;
            case 'm':
                break;
            case 'n':
                break;
            case 'o':
                options->out_file = strdup(optarg);
                break;
            case 'p':
                options->forceparallel = atoi(optarg);
                break;
            case 'q':
                options->silent = 1;
                break;
            case 's':
                break;
            case 'u':
                options->ufactor = atoi(optarg);
                break;
            case 'v':
                printf("PLUTO %s - An automatic parallelizer and locality optimizer\n\
Copyright (C) 2007--2008  Uday Kumar Bondhugula\n\
This is free software; see the source for copying conditions.  There is NO\n\
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n", PLUTO_VERSION);
                pluto_options_free(options);
                return 3;
            default:
                usage_message();
                pluto_options_free(options);
                return 4;
        }
    }

    if (optind <= argc-1)   {
        srcFileName = alloca(strlen(argv[optind])+1);
        strcpy(srcFileName, argv[optind]);
    }else{
        /* No non-option argument was specified */
        usage_message();
        pluto_options_free(options);
        return 5;
    }

    if (options->fusesends && options->mpiomp) {
        fprintf(stderr, "[pluto] Error: fusesends should not be used with mpiomp\n");
        return 7;
    }

    if (!options->isldepaccesswise && options->distmem) {
        fprintf(stderr, "[pluto] Error: --isldepstmtwise can't be used for distmem parallelization\n");
        return 7;
    }

    /* Make options consistent */
    if (options->isldep && options->candldep) {
        printf("[pluto] Error: only one of isldep and candldep should be specified)\n");
        pluto_options_free(options);
        usage_message();
        return 1;
    }

    /* isldep is the default */
    if (!options->isldep && !options->candldep) {
        options->isldep = 1;
    }

    if (options->lastwriter && options->candldep) {
        printf("[pluto] ERROR: --lastwriter is only supported with --isldep\n");
        pluto_options_free(options);
        usage_message();
        return 1;
    }

    /* lastwriter is default under distmem */
    if (options->distmem && nolastwriter == 0)    {
        printf("[pluto] Turning on lastwriter\n");
        options->lastwriter = 1;
    }

    if (options->lastwriter && nolastwriter) {
        printf("[pluto] WARNING: both --lastwriter, --nolastwriter are on\n");
        printf("[pluto] disabling --lastwriter\n");
        options->lastwriter = 0;
    }

    if (options->identity == 1) {
        options->partlbtile = 0;
        options->lbtile = 0;
    }

    if (options->dynschedule && options->dynschedule_graph) {
        options->dynschedule_graph = 0;
    }

    if (options->dynschedule && options->dynschedule_graph_old) {
        options->dynschedule_graph_old = 0;
    }

    if (options->dynschedule_graph && options->dynschedule_graph_old) {
        options->dynschedule_graph_old = 0;
    }

    if (options->dynschedule_graph || options->dynschedule_graph_old) {
        assert(options->distmem == 0);
    }

    if (options->distmem == 1 && options->dynschedule == 0 && options->parallel == 0)    {
        options->parallel = 1;
    }

    if (options->partlbtile == 1 && options->lbtile == 0)    {
        options->lbtile = 1;
    }

    if (options->lbtile == 1 && options->tile == 0)    {
        options->tile = 1;
    }

    if (options->multipipe == 1 && options->parallel == 0)    {
        fprintf(stdout, "Warning: multipipe needs parallel to be on; turning on parallel\n");
        options->parallel = 1;
    }

    if ((options->dynschedule || options->dynschedule_graph || options->dynschedule_graph_old) && !options->tile)  {
        fprintf(stderr, "[Pluto] WARNING: dynschedule needs tile to be on; turning on tile\n");
        options->tile = 1;
    }

    if (options->distmem == 1 && options->multi_level_distribution == 1) {
        options->multipipe = 1; // even for dynschedule
    }else{
        if (options->multipipe == 1 && options->parallel == 0)    {
            fprintf(stdout, "Warning: multipipe needs parallel to be on; turning on parallel\n");
            options->parallel = 1;
        }

        if ((options->dynschedule || options->dynschedule_graph || options->dynschedule_graph_old) && (options->multipipe))  {
            fprintf(stderr, "[Pluto] WARNING: --multipipe option not needed with --dynschedule; turning off multipipe\n");
            options->multipipe = 0;
        }
    }

    if ((options->dynschedule || options->dynschedule_graph || options->dynschedule_graph_old) && (options->parallel))  {
        fprintf(stderr, "[Pluto] WARNING: --parallel option not needed with --dynschedule; turning off parallel\n");
        options->parallel = 0;
    }

    //reset commopt and commopt_foifi when commopt_fop is selected, default of commopt is 1
    if(options->commopt && options->commopt_fop) {
        options->commopt = 0;
    }
    if(options->commopt_foifi && options->commopt_fop) {
        options->commopt_foifi = 0;
    }

    //reset commopt when commopt_foifi is selected, default of commopt is 1
    if(options->commopt && options->commopt_foifi) {
        options->commopt = 0;
    }

    if(options->data_dist){
    	options->verify_output = 1;
    }


    /* Extract polyhedral representation from osl scop */
    PlutoProg *prog = NULL; 

    osl_scop_p scop = NULL;
    char *irroption = NULL;

    /* Extract polyhedral representation from input program */
    if (options->pet) {
        isl_ctx *pctx = isl_ctx_alloc();
        pscop = pet_scop_extract_from_C_source(pctx, srcFileName, NULL);

        if (!pscop) {
            fprintf(stdout, "[pluto] No SCoPs extracted or error extracting SCoPs  using pet\n");
            pluto_options_free(options);
            isl_ctx_free(pctx);
            return 12;
        }
        t_start = rtclock();
        prog = pet_to_pluto_prog(pscop, pctx, options);
        t_d = rtclock() - t_start;

        pet_scop_free(pscop);
        isl_ctx_free(pctx);

        FILE *srcfp = fopen(".srcfilename", "w");
        if (srcfp)    {
            fprintf(srcfp, "%s\n", srcFileName);
            fclose(srcfp);
        }
    }else{
        /* Extract polyhedral representation from clan scop */
        if(!strcmp(srcFileName, "stdin")){  //read from stdin
            src_fp = stdin;
            osl_interface_p registry = osl_interface_get_default_registry();
            t_start = rtclock();
            scop = osl_scop_pread(src_fp, registry, PLUTO_OSL_PRECISION);
            t_d = rtclock() - t_start;
        }else{  // read from regular file

            src_fp  = fopen(srcFileName, "r");

            if (!src_fp)   {
                fprintf(stderr, "pluto: error opening source file: '%s'\n", srcFileName);
                pluto_options_free(options);
                return 6;
            }

            /* Extract polyhedral representation from input program */
            clan_options_p clanOptions = clan_options_malloc();

            if (options->readscop){
                osl_interface_p registry = osl_interface_get_default_registry();
                t_start = rtclock();
                scop = osl_scop_pread(src_fp, registry, PLUTO_OSL_PRECISION);
                t_d = rtclock() - t_start;
            }else{
                t_start = rtclock();
                scop = clan_scop_extract(src_fp, clanOptions);
                t_d = rtclock() - t_start;
            }

            if (!scop || !scop->statement)   {
                fprintf(stderr, "Error extracting polyhedra from source file: \'%s'\n",
                        srcFileName);
                pluto_options_free(options);
                return 8;
            }
            FILE *srcfp = fopen(".srcfilename", "w");
            if (srcfp)    {
                fprintf(srcfp, "%s\n", srcFileName);
                fclose(srcfp);
            }

            clan_options_free(clanOptions);

            /* IF_DEBUG(clan_scop_print_dot_scop(stdout, scop, clanOptions)); */
        }

        /* Convert clan scop to Pluto program */
        prog = scop_to_pluto_prog(scop, options);

        /* Backup irregular program portion in .scop. */
        osl_irregular_p irreg_ext = NULL;
        irreg_ext = osl_generic_lookup(scop->extension, OSL_URI_IRREGULAR);
        if (irreg_ext!=NULL)
            irroption = osl_irregular_sprint(irreg_ext);  //TODO: test it
        osl_irregular_free(irreg_ext);
    }
    IF_MORE_DEBUG(pluto_prog_print(stdout, prog));

    int dim_sum=0;
    for (i=0; i<prog->nstmts; i++) {
        dim_sum += prog->stmts[i]->dim;
    }

    
   
    if (!options->silent)   {
        fprintf(stdout, "[pluto] Number of statements: %d\n", prog->nstmts);
        fprintf(stdout, "[pluto] Total number of loops: %d\n", dim_sum);
        fprintf(stdout, "[pluto] Number of deps: %d\n", prog->ndeps);
        fprintf(stdout, "[pluto] Maximum domain dimensionality: %d\n", prog->nvar);
        fprintf(stdout, "[pluto] Number of parameters: %d\n", prog->npar);
    }

    if (options->iss) {
        // PlutoConstraints *dom = pluto_constraints_read(stdin);
        // printf("Input set\n");
        // pluto_constraints_compact_print(stdout, dom);
        // PlutoConstraints **doms = malloc(1*sizeof(PlutoConstraints *));
        // doms[0] = dom;
        // pluto_find_iss(doms, 1, 1, NULL);
        // PlutoMatrix *mat = pluto_matrix_input(stdin);
        // pluto_constraints_print(stdout, dom);
        // pluto_matrix_print(stdout, mat);
        // PlutoConstraints *farkas = farkas_affine(dom, mat);
        //pluto_constraints_pretty_print(stdout, farkas);
        // pluto_constraints_free(dom);
        // pluto_options_free(options);
        pluto_iss_dep(prog);
    }

    t_start = rtclock();
    /* Auto transformation */
    if (!options->identity) {
        pluto_auto_transform(prog);
    }
    t_t = rtclock() - t_start;

    if (options->identity_data_dist && options->data_dist) {
    	pluto_data_dist_identity_trans(prog);
    }

    pluto_detect_transformation_properties(prog);

    if (!options->silent)   {
        fprintf(stdout, "[pluto] Affine transformations [<iter coeff's> <param> <const>]\n\n");
        /* Print out transformations */
        pluto_transformations_pretty_print(prog);
        pluto_print_hyperplane_properties(prog);
    }

    if (options->tile)   {
        pluto_tile(prog);
    }else{
        if (options->intratileopt) {
            pluto_intra_tile_optimize(prog, 0);
        }
    }

    if (options->parallel && !options->tile && !options->identity)   {
        /* Obtain wavefront/pipelined parallelization by skewing if
         * necessary */
        int nbands;
        Band **bands;
        bands = pluto_get_outermost_permutable_bands(prog, &nbands);
        bool retval = pluto_create_tile_schedule(prog, bands, nbands);
        pluto_bands_free(bands, nbands);

        /* If the user hasn't supplied --tile and there is only pipelined
         * parallelism, we will warn the user */
        if (retval)   {
            printf("[pluto] WARNING: pipelined parallelism exists and --tile is not used.\n");
            printf("use --tile for better parallelization \n");
            IF_DEBUG(fprintf(stdout, "[pluto] After skewing:\n"););
            IF_DEBUG(pluto_transformations_pretty_print(prog););
            IF_DEBUG(pluto_print_hyperplane_properties(prog););
        }
    }

    if (options->unroll || options->polyunroll)    {
        /* Will generate a .unroll file */
        /* plann/plorc needs a .params */
        FILE *paramsFP = fopen(".params", "w");
        if (paramsFP)   {
            int i;
            for (i=0; i<prog->npar; i++)  {
                fprintf(paramsFP, "%s\n", prog->params[i]);
            }
            fclose(paramsFP);
        }
        pluto_detect_mark_unrollable_loops(prog);
    }

    if (options->polyunroll)    {
        /* Experimental */
        for (i=0; i<prog->num_hyperplanes; i++)   {
            if (prog->hProps[i].unroll)  {
                unroll_phis(prog, i, options->ufactor);
            }
        }
    }

    if(!options->pet && !strcmp(srcFileName, "stdin")){
        //input stdin == output stdout
        pluto_populate_scop(scop, prog, options);
        osl_scop_print(stdout, scop);
    }else{  // do the usual Pluto stuff

        /* NO MORE TRANSFORMATIONS BEYOND THIS POINT */
        /* Since meta info about loops
         * is printed to be processed by scripts - if transformations are
         * performed, changed loop order/iterator names will be missed  */
        gen_unroll_file(prog);

        char *outFileName;
        char *headerFileName; // used only by options->distmem, options->dynschedule, options->dynschedule_graph
        char *cloogFileName;
        char *dynschedFileName;
        char *sigmaFileName = NULL; // used only by options->distmem, options->dynschedule, options->dynschedule_graph
        char *piFileName = NULL; // used only by options->distmem
        char *bname, *basec;
        if (options->out_file == NULL)  {
            /* Get basename, remove .c extension and append a new one */
            basec = strdup(srcFileName);
            bname = basename(basec);

            if (strlen(bname) >= 2 && !strcmp(bname+strlen(bname)-2, ".c")) {
                outFileName = malloc(strlen(bname)-2+strlen(".pluto.c")+1);
                strncpy(outFileName, bname, strlen(bname)-2);
                outFileName[strlen(bname)-2] = '\0';
            }else{
                outFileName = malloc(strlen(bname)+strlen(".pluto.c")+1);
                strcpy(outFileName, bname);
            }
            strcat(outFileName, ".pluto.c");
        }else{
            basec = strdup(options->out_file);
            bname = basename(basec);

            outFileName = malloc(strlen(options->out_file)+1);
            strcpy(outFileName, options->out_file);
        }

        if (strlen(bname) >= 2 && !strcmp(bname+strlen(bname)-2, ".c")) {
            headerFileName = malloc(strlen(bname)-2+strlen(".h")+1);
            cloogFileName = malloc(strlen(bname)-2+strlen(".pluto.cloog")+1);
            dynschedFileName = malloc(strlen(bname)-2+strlen(".pluto.append.c")+1);
            strncpy(headerFileName, bname, strlen(bname)-2);
            strncpy(cloogFileName, bname, strlen(bname)-2);
            strncpy(dynschedFileName, bname, strlen(bname)-2);
            headerFileName[strlen(bname)-2] = '\0';
            cloogFileName[strlen(bname)-2] = '\0';
            dynschedFileName[strlen(bname)-2] = '\0';
        }else{
            headerFileName = malloc(strlen(bname)+strlen(".h")+1);
            cloogFileName = malloc(strlen(bname)+strlen(".pluto.cloog")+1);
            dynschedFileName = malloc(strlen(bname)+strlen(".pluto.append.c")+1);
            strcpy(headerFileName, bname);
            strcpy(cloogFileName, bname);
            strcpy(dynschedFileName, bname);
        }
        strcat(headerFileName, ".h");
        strcat(cloogFileName, ".pluto.cloog");
        strcat(dynschedFileName, ".pluto.append.c");
        free(basec);

        cloogfp = fopen(cloogFileName, "w+");
        if (!cloogfp)   {
            fprintf(stderr, "[Pluto] Can't open .cloog file: '%s'\n", cloogFileName);
            free(cloogFileName);
            pluto_options_free(options);
            pluto_prog_free(prog);
            return 9;
        }
        free(cloogFileName);

        outfp = fopen(outFileName, "w");
        if (!outfp) {
            fprintf(stderr, "[Pluto] Can't open file '%s' for writing\n", outFileName);
            free(outFileName);
            pluto_options_free(options);
            pluto_prog_free(prog);
            fclose(cloogfp);
            return 10;
        }

        int retval = 1;

#ifdef PLUTO_OPENCL
    if (options->distmem || options->opencl || options->dynschedule || options->dynschedule_graph || options->data_dist)  
#else
        if (options->distmem || options->dynschedule || options->dynschedule_graph || options->data_dist)  
#endif
        {
            sigmaFileName = malloc(strlen("sigma_")+strlen(outFileName)+1);
            strcpy(sigmaFileName, "sigma_");
            strcat(sigmaFileName, outFileName);

            sigmafp = fopen(sigmaFileName, "w");
            if (!sigmafp) {
                fprintf(stderr, "[Pluto] Can't open file: '%s'\n", sigmaFileName);
                free(sigmaFileName);
                pluto_options_free(options);
                pluto_prog_free(prog);
                return 12;
            }
            fprintf(sigmafp, "#include \"%s\"\n", headerFileName);
            if (options->distmem)  {
                fprintf(sigmafp, "#include \"polyrt.h\"\n");
            if(options->data_dist)
            fprintf(sigmafp, "#include \"buffer_manager.h\"\n");

            }

            headerfp = fopen(headerFileName, "w");
            if (!headerfp) {
                fprintf(stderr, "[Pluto] Can't open file: '%s'\n", headerFileName);
                free(headerFileName);
                pluto_options_free(options);
                pluto_prog_free(prog);
                return 13;
            }
            fprintf(headerfp, "#include \"polyrt.h\"\n");
            if(options->data_dist)
            fprintf(headerfp, "#include \"buffer_manager.h\"\n");

            piFileName = malloc(strlen("pi_")+strlen(outFileName)+1);
            strcpy(piFileName, "pi_");
            strcat(piFileName, outFileName);

            FILE *pifp = fopen(piFileName, "w");
            free(piFileName);
            if (!pifp) {
                fprintf(stderr, "[Pluto] Can't open file: '%s'\n", piFileName);
                pluto_options_free(options);
                pluto_prog_free(prog);
                return 14;
            }

#ifdef PLUTO_OPENCL
            if (options->distmem || options->opencl)  
#else
            if (options->distmem)  
#endif
            {
                retval = pluto_distmem_parallelize(prog, sigmafp, headerfp, pifp);

            }else{ // shared-memory
                if (options->dynschedule || options->dynschedule_graph) {
                    retval = pluto_dynschedule_parallelize(prog, sigmafp, headerfp, pifp);
                }
                else if (options->data_dist) {
					fprintf(outfp, "#include \"%s\"\n", headerFileName);
                    retval = pluto_shared_memory_data_dist(prog, headerfp, outfp);
                }
            }


            fclose(pifp);
            IF_DEBUG(pluto_transformations_pretty_print(prog));
            IF_DEBUG(pluto_print_hyperplane_properties(prog));
        }else if (options->dynschedule_graph_old) {
            dynschedfp = fopen(dynschedFileName, "w");
            if (!dynschedfp) {
                fprintf(stderr, "[Pluto] Can't open file %s for writing\n", dynschedFileName);
                free(dynschedFileName);
                pluto_options_free(options);
                pluto_prog_free(prog);
                fclose(cloogfp);
                fclose(outfp);
                return 11;
            }
        }

        if (!options->pet) pluto_detect_scalar_dimensions(prog);
        if (options->moredebug) {
            printf("After scalar dimension detection (final transformations)\n");
            pluto_transformations_pretty_print(prog);
        }

        /* Generate .cloog file */
        pluto_gen_cloog_file(cloogfp, prog);
        /* Add the <irregular> tag from clan, if any */
        if(!options->pet) {
            if (irroption) {
                fprintf(cloogfp, "<irregular>\n%s\n</irregular>\n\n", irroption);
                free(irroption);
            }
        }
        rewind(cloogfp);

        /* Very important: Dont change the order of calls to print_dynsched_file
         * between pluto_gen_cloog_file() and pluto_*_codegen()
         */
        if (options->dynschedule_graph_old) {
            print_dynsched_file(srcFileName, cloogfp, dynschedfp, prog);
            rewind(cloogfp);
        }

        /* Generate code using Cloog and add necessary stuff before/after code */
        if (options->distmem && !retval){
            fprintf(outfp, "#include \"%s\"\n", headerFileName);
            t_start = rtclock();
            pluto_distmem_codegen(prog, cloogfp, sigmafp, outfp, headerfp);
            t_c = rtclock() - t_start;
        }else if (options->dynschedule && !retval){ // shared-memory
            fprintf(outfp, "#include \"%s\"\n", headerFileName);
            pluto_dynschedule_codegen(prog, sigmafp, outfp, headerfp);
        }else if (options->dynschedule_graph && !retval){ // shared-memory
            fprintf(outfp, "#include \"%s\"\n", headerFileName);
            pluto_dynschedule_graph_codegen(prog, sigmafp, outfp, headerfp);
        }
        else if(!options->distmem && !options->dynschedule && options->data_dist){


//        	if(!headerfp){
//				headerfp = fopen(headerFileName, "w");
//				if (!headerfp) {
//					fprintf(stderr, "[Pluto] Can't open file: '%s'\n", headerFileName);
//					free(headerFileName);
//					pluto_options_free(options);
//					pluto_prog_free(prog);
//					return 13;
//				}
//				fprintf(headerfp, "#include \"polyrt.h\"\n");
//				if(options->data_dist)
//				fprintf(headerfp, "#include \"buffer_manager.h\"\n");
//        	}
//
//			pluto_shared_memory_data_dist(prog, headerfp);

//            fprintf(outfp, "#include \"%s\"\n", headerFileName);
// t_start = rtclock();
//            pluto_sharedmem_data_dist_codegen(prog, cloogfp, outfp, headerfp);
            // t_c = rtclock() - t_start;

        }else{
            if (options->distmem) { // no parallel loops to distribute
                // ensure only one processor will output the data
                fprintf(outfp, "#include <mpi.h>\n\n");
                fprintf(outfp, "#define MPI \n\n");
                fprintf(outfp, "\n##ifndef GLOBAL_MY_RANK\n\tint my_rank;\n##endif\n");
                fprintf(outfp, "\tMPI_Init(NULL, NULL);\n");
                fprintf(outfp, "\tMPI_Comm_rank(MPI_COMM_WORLD, &my_rank);\n");
                fprintf(outfp, "\tMPI_Finalize();\n");
            }
            t_start = rtclock();
            pluto_multicore_codegen(cloogfp, (options->dynschedule_graph_old) ? dynschedfp : outfp, prog);
            t_c = rtclock() - t_start;
        }

#ifdef PLUTO_OPENCL
        if (options->distmem || options->opencl || options->dynschedule || options->dynschedule_graph)  
#else
        if (options->distmem || options->dynschedule || options->dynschedule_graph)  
#endif
        {
            fclose(sigmafp);
            fclose(headerfp);
        }

        FILE *tmpfp = fopen(".outfilename", "w");
        if (tmpfp)    {
            fprintf(tmpfp, "%s\n", outFileName);
            fclose(tmpfp);
            PLUTO_MESSAGE(printf( "[Pluto] Output written to %s\n", outFileName););
        }
        free(outFileName);

        if (options->distmem || options->dynschedule || options->dynschedule_graph || options->data_dist) {
            tmpfp = fopen(".sigmafilename", "w");
            if (tmpfp) {
                fprintf(tmpfp, "%s\n", sigmaFileName);
                fclose(tmpfp);
            }
            free(sigmaFileName);

            tmpfp = fopen(".headerfilename", "w");
            if (tmpfp) {
                fprintf(tmpfp, "%s\n", headerFileName);
                fclose(tmpfp);
            }
        }

        if (options->dynschedule_graph_old) {
            tmpfp = fopen(".appendfilename", "w");
            if (tmpfp) {
                fprintf(tmpfp, "%s\n", dynschedFileName);
                fclose(tmpfp);
            }
        }
        free(dynschedFileName);

        /* create main file for dynschedule_graph_old */
        if (options->dynschedule_graph_old) {
            fprintf(outfp,"\n\n#include \"scheduler.h\"\n");
            fprintf(outfp,"\n\n\tgenerate_dag();\n");
            fprintf(outfp,"\t##ifdef __STATIC_SCHEDULE__\n");
            fprintf(outfp,"\t\tinit_schedule(atoi(argv[1]));\n");
            fprintf(outfp,"\t\tschedule_dag_using_mcp();\n");
            fprintf(outfp,"\t##endif\n");
            fprintf(outfp,"\tdag_execute();\n");
            fprintf(outfp,"\t##ifdef __STATIC_SCHEDULE__\n");
            fprintf(outfp,"\t\tfree_schedule();\n");
            fprintf(outfp,"\t##endif\n");
            fprintf(outfp,"\tfree_dag();\n");
            fclose(dynschedfp);
        }
        fclose(cloogfp);
        fclose(outfp);
        free(headerFileName);
    }


    t_all = rtclock() - t_start_all;

    if (options->time && !options->silent) {
        printf("\n[pluto] Timing statistics\n[pluto] SCoP extraction + dependence analysis time: %0.6lfs\n", t_d);
        printf("[pluto] Auto-transformation time: %0.6lfs\n", t_t);
        printf("[pluto] Code generation time: %0.6lfs\n", t_c);
        printf("[pluto] Other/Misc time: %0.6lfs\n", t_all-t_c-t_t-t_d);
        printf("[pluto] Total time: %0.6lfs\n", t_all);
        printf("[pluto] All times: %0.6lf %0.6lf %.6lf %.6lf\n", t_d, t_t, t_c,
             t_all-t_c-t_t-t_d);
    }

    pluto_prog_free(prog);
    pluto_options_free(options);

    osl_scop_free(scop);

    return 0;
}
