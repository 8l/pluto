# PLuTo
#
# Copyright (C) Uday Kumar Bondhugula 2007
#

SUBDIRS = 

MAINTAINERCLEANFILES = Makefile.in

bin_PROGRAMS = pluto
pluto_SOURCES = \
	constraints.c \
	math_support.c \
	ddg.c \
	pluto.c \
	post_transform.c \
	framework.c \
	gpucode.c \
	gpucode.h \
	gpuloc.c \
	gpuloc.h \
	cuda_common.h \
	cuda_common.c \
	main.c \
	constraints.h \
	math_support.h \
	pluto.h \
	program.c \
	program.h \
	schedule.c \
	schedule.h \
	tile.c
pluto_CPPFLAGS = -DSCOPLIB_INT_T_IS_LONGLONG -DCLOOG_INT_GMP \
   -I../piplib-1.4.0/include \
   -I../clan-0.6.0/include \
   -I../clan-0.6.0/scoplib/scoplib-0.2.0/include \
   -I../candl-0.4.0/include \
   -I../isl/include -I../cloog-isl/include \
   -I../barvinok/polylib/include -I../barvinok
pluto_LDADD = ../candl-0.4.0/source/libcandl.la \
   ../clan-0.6.0/source/libclan.la \
   ../piplib-1.4.0/libpiplib64.la \
   ../cloog-isl/libcloog-isl.la ../barvinok/libbarvinok.la ../isl/libisl.la \
   -lm -lstdc++
INCLUDES = 
