/*
  quadtree header file for Jacobi on quadtree project

  Jim Teresco

  Initial implementation: Sun Mar  5 19:11:28 EST 2006

*/

#ifndef __QUADTREE_H__
#define __QUADTREE_H__

#include <stdio.h>

/* the recursive quadtree structure */
typedef struct quadtree {
  /* solution values (previous and current) */
  double value;
  double previous;

  /* bounding box */
  double top, left, bottom, right;

  /* parent pointer */
  struct quadtree *parent;

  /* child pointers */
  /* children[0] is northwest quadrant */
  /* children[1] is northeast quadrant */
  /* children[2] is southwest quadrant */
  /* children[3] is southeast quadrant */
  /* children[3] set to non-zero when children[0] is NULL indicates a mark
     for refinement */
  struct quadtree *children[4];

} quadtree_type;

/* callback function type for quadtree visitor function */
typedef void QUADRANT_VISITOR_FN(void *, struct quadtree *);

extern struct quadtree *new_quadtree(double, double, double, double, double, 
				     struct quadtree *);
extern void free_quadtree(struct quadtree *);
extern void refine_leaf_quadrant(struct quadtree *);
extern int is_leaf_quadrant(struct quadtree *);
extern void visit_all_leaf_quadrants(struct quadtree *, QUADRANT_VISITOR_FN *,
				     void *);
extern double quadrant_centroid_x(struct quadtree *);
extern double quadrant_centroid_y(struct quadtree *);
extern int quadrant_contains(struct quadtree *, double x, double y);
extern struct quadtree *child_quadrant_containing(struct quadtree *, 
						  double x, double y);
extern void print_quadrant(FILE *fp, struct quadtree *);
extern void print_leaf_quadrants(FILE *fp, struct quadtree *);
extern int quadrant_level(struct quadtree *);

extern double quadrant_value(struct quadtree *);
extern double quadrant_previous(struct quadtree *);
extern double quadrant_side_value(struct quadtree *, char dir);
extern double quadrant_side_previous(struct quadtree *, char dir);
extern void quadrant_set_value(struct quadtree *, double);
extern void quadrant_set_previous(struct quadtree *, double);

extern struct quadtree *neighbor_quadrant(struct quadtree *, char dir);

extern void mark_quadrant_for_refinement(struct quadtree *);
extern int is_quadrant_marked_for_refinement(struct quadtree *);

extern int num_leaf_quadrants(struct quadtree *);

#endif
