#
# Makefile for Jacobi iteration on an adaptive quadtree, serial version
#
# Jim Teresco, CS 341, Mount Holyoke College
#
CFILES = quadtree.c quadtree_solver.c com_abort.c
OFILES = $(CFILES:.c=.o)
CC=gcc

quadtree_solver:	$(OFILES)
	$(CC) $(CFLAGS) -o quadtree_solver $(OFILES) -lm

quadtree.o:	quadtree.h macros.h
quadtree_solver.o:	quadtree.h

clean::
	/bin/rm -f quadtree_solver $(OFILES)
