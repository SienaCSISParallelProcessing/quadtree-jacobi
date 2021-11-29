#
# Makefile for Jacobi iteration on an adaptive quadtree, pthreads version
#
# Jim Teresco, CSIS-335, Siena College, Fall 2021
#
CFILES = quadtree.c quadtree_solver.c com_abort.c
OFILES = $(CFILES:.c=.o)
CC=gcc
CFLAGS=-g

quadtree_solver:	$(OFILES)
	$(CC) $(CFLAGS) -o quadtree_solver $(OFILES) -lm -lpthread

quadtree.o:	quadtree.h macros.h
quadtree_solver.o:	quadtree.h

clean::
	/bin/rm -f quadtree_solver $(OFILES)
