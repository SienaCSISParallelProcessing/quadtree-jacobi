/*
  Jacobi solver on quadtree project main program

  Jim Teresco

  Initial implementation:
  Sun Mar  5 19:40:41 EST 2006

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "quadtree.h"
#include "macros.h"

/********************************/
/* problem definition functions */
/********************************/
static double initial_cond(double x, double y) {

  return 0;
}

/* the y=1 boundary condition */
static double bc_top(double x) {

  if (x < 0.5) return 1-2*x;
  return 0;
}

/* the y=0 boundary condition */
static double bc_bottom(double x) {

  return 0;
}

/* the x=0 boundary condition */
static double bc_left(double y) {

  if (y < 0.5) return 0;
  return 2*(y-0.5);
}

/* the x=1 boundary condition */
static double bc_right(double y) {

  return 0;
}

/* use this to apply special conditions like a point heat source in
   the interior of the domain.  Return 1 if the leaf quadrant gets its
   values modified (modify both value and previous for simplicity), 0
   otherwise.
*/
static int apply_other_bc(struct quadtree *leaf) {

  if (quadrant_contains(leaf, 0.8, 0.2)) {
    quadrant_set_value(leaf, 3);
    quadrant_set_previous(leaf, 3);
    return 1;
  }

  return 0;
}

/* function to set initial conditions on a given quadrant (callback) */
static void set_init_cond(void *data, struct quadtree *leaf) {

  quadrant_set_value(leaf, initial_cond(quadrant_centroid_x(leaf),
					quadrant_centroid_y(leaf)));
}

/* function to refine all quadrants in the form of a visitor callback */
static void do_refine(void *data, struct quadtree *leaf) {

  refine_leaf_quadrant(leaf);
}

/* function to mark leaf quadrants if they have neighbors whose values
   differ by more than the refinement tolerance - form of a visitor
   callback */
static void check_if_refinement_needed(void *data, struct quadtree *leaf) {

  struct quadtree *nbr;
  double *tol = (double *)data;
  int level = quadrant_level(leaf);
  char *dirs = "nsew";
  int dir;

  /* check each neighbor */
  for (dir = 0; dir < 4; dir++) {
    nbr = neighbor_quadrant(leaf, dirs[dir]);
    /* only compare neighbor values if the neighbor is at least as refined */
    if (nbr && level <= quadrant_level(nbr) && 
	fabs(quadrant_value(leaf) - quadrant_value(nbr)) > *tol) {
      mark_quadrant_for_refinement(leaf);
      return;
    }
  }
}

/* function to do a refinement on any quadrant that has been marked
   for refinement */
static void refine_if_marked(void *data, struct quadtree *leaf) {

  int *any_ref = (int *)data;

  if (is_quadrant_marked_for_refinement(leaf)) {
    refine_leaf_quadrant(leaf);
    *any_ref += 1;
  }
}

/* check error and refine if necessary, return 0 if no refinement
   was needed, number of refinements done otherwise */
int calc_error_and_refine(struct quadtree *root, double tol) {
  int any_refs = 0;

  /* mark quadrants for refinement */
  visit_all_leaf_quadrants(root, check_if_refinement_needed, &tol);

  /* do the refinements */
  visit_all_leaf_quadrants(root, refine_if_marked, &any_refs);

  return any_refs;
}

/* do a Jacobi iteration on one cell computing "previous" from "value" */
void do_jacobi_iter_phase1(void *data, struct quadtree *leaf) {

  /* data is NULL for this one */
  double north, south, east, west;
  struct quadtree *nbr;
  
  /* apply special boundary conditions if necessary */
  if (apply_other_bc(leaf)) return;

  /* look north */
  nbr = neighbor_quadrant(leaf, 'n');
  if (nbr) {
    north = quadrant_side_value(nbr, 's');
  }
  else {
    north = bc_top(quadrant_centroid_x(leaf));
  }
  
  /* look south */
  nbr = neighbor_quadrant(leaf, 's');
  if (nbr) {
    south = quadrant_side_value(nbr, 'n');
  }
  else {
    south = bc_bottom(quadrant_centroid_x(leaf));
  }
  
  /* look west */
  nbr = neighbor_quadrant(leaf, 'w');
  if (nbr) {
    west = quadrant_side_value(nbr, 'e');
  }
  else {
    west = bc_left(quadrant_centroid_y(leaf));
  }
  
  /* look east */
  nbr = neighbor_quadrant(leaf, 'e');
  if (nbr) {
    east = quadrant_side_value(nbr, 'w');
  }
  else {
    east = bc_right(quadrant_centroid_y(leaf));
  }

  quadrant_set_previous(leaf,(north + east + south + west)*0.25);
}

/* do a Jacobi iteration on one cell computing "value" from "previous" */
void do_jacobi_iter_phase2(void *data, struct quadtree *leaf) {

  /* data is a pointer to the max error (abs(previous-value)) for
     this iteration so far */
  double *maxerr = (double *)data;
  double north, south, east, west;
  struct quadtree *nbr;
  
  /* apply special boundary conditions if necessary */
  if (apply_other_bc(leaf)) return;

  /* look north */
  nbr = neighbor_quadrant(leaf, 'n');
  if (nbr) {
    north = quadrant_side_previous(nbr, 's');
  }
  else {
    north = bc_top(quadrant_centroid_x(leaf));
  }
  
  /* look south */
  nbr = neighbor_quadrant(leaf, 's');
  if (nbr) {
    south = quadrant_side_previous(nbr, 'n');
  }
  else {
    south = bc_bottom(quadrant_centroid_x(leaf));
  }
  
  /* look west */
  nbr = neighbor_quadrant(leaf, 'w');
  if (nbr) {
    west = quadrant_side_previous(nbr, 'e');
  }
  else {
    west = bc_left(quadrant_centroid_y(leaf));
  }
  
  /* look east */
  nbr = neighbor_quadrant(leaf, 'e');
  if (nbr) {
    east = quadrant_side_previous(nbr, 'w');
  }
  else {
    east = bc_right(quadrant_centroid_y(leaf));
  }

  quadrant_set_value(leaf, (north + east + south + west)*0.25);

  if (*maxerr < fabs(quadrant_value(leaf) - quadrant_previous(leaf))) {
    *maxerr = fabs(quadrant_value(leaf) - quadrant_previous(leaf));
  }
}


/* main driver program */
int main(int argc, char *argv[]) {

  struct quadtree *root;
  int i;
  int init_ref;
  double jac_tol;
  int jac_maxiter;
  double ref_tol;
  int ref_maxiter;
  int jac_iter_num;
  double max_jac_diff;
  FILE *fp;
  int ref_count;
  int output_level;
  char filename[256];

  /* parse command-line parameters */
  if (argc != 7) {
    fprintf(stderr, "Usage: %s: init_ref jac_tol jac_maxiter ref_tol ref_maxiteriter output_level\n", argv[0]);
    return 1;
  }

  /* first parameter is number of initial refinement levels */
  init_ref = atoi(argv[1]);

  /* second parameter is Jacobi iteration tolerance */
  jac_tol = atof(argv[2]);

  /* third parameter is maximum number of Jacobi iterations */
  jac_maxiter = atoi(argv[3]);

  /* fourth parameter is refinement tolerance */
  ref_tol = atof(argv[4]);

  /* fifth parameter is maximum number of refinements */
  ref_maxiter = atoi(argv[5]);

  /* sixth parameter is an output level */
  /* 0 = no solution output */
  /* 1 = final solution only */
  /* 2 = solution after each jacobi convergence */
  /* 3 = solution after each iteration */
  output_level = atoi(argv[6]);

  printf("Starting with parameters\n");
  printf("Initial number of refinement levels: %d\n", init_ref);
  printf("Jacobi tolerance: %10.8f\n", jac_tol);
  printf("Maximum number of Jacobi iterations: %d\n", jac_maxiter);
  printf("Refinement tolerance: %10.8f\n", ref_tol);
  printf("Maximum number of refinements: %d\n", ref_maxiter);
  printf("Output level: %d\n", output_level);

  /* create the root quadtree */
  root = new_quadtree(1,0,0,1,0,NULL);

  /* refine down the requested number of levels globally to start up */
  for (i=0; i<init_ref; i++) {
    visit_all_leaf_quadrants(root, do_refine, NULL);
  }

  /* insert initial conditions */
  visit_all_leaf_quadrants(root, set_init_cond, NULL);

  /* print the initial solution */
  if (output_level >= 3) {
    GRACEFULLY_OPEN(fp, "solution.0.0.dat", "w");
    print_leaf_quadrants(fp, root);
    GRACEFULLY_CLOSE(fp, "solution.0.0.dat");
  }

  /* do the computation */
  
  /* outer loop controls the adaptive refinement -- keep going until
     we compute a solution to acceptable accuracy. */
  ref_count = -1;
  do {
    printf("Compute on quadtree with %d leaves...\n", num_leaf_quadrants(root));
    ref_count++;

    /* Jacobi iteration until the solution settles down or we 
       hit the max number of iterations */
    jac_iter_num = 0;
    max_jac_diff = 0;

    do {
      visit_all_leaf_quadrants(root, do_jacobi_iter_phase1, NULL);

      /* are we printing every iteration? */
      if (output_level >= 3) {
	sprintf(filename, "solution.%d.%d.dat", ref_count, jac_iter_num+1);
	GRACEFULLY_OPEN(fp, filename, "w");
	print_leaf_quadrants(fp, root);
	GRACEFULLY_CLOSE(fp, filename);
      }

      max_jac_diff = 0;
      visit_all_leaf_quadrants(root, do_jacobi_iter_phase2, &max_jac_diff);

      jac_iter_num += 2;  /* we always do 2 iters */

      /* are we printing every iteration? */
      if (output_level >= 3) {
	sprintf(filename, "solution.%d.%d.dat", ref_count, jac_iter_num);
	GRACEFULLY_OPEN(fp, filename, "w");
	print_leaf_quadrants(fp, root);
	GRACEFULLY_CLOSE(fp, filename);
      }

      printf("Completed Jacobi iteration %d, max_jac_diff=%f\n",
	     jac_iter_num, max_jac_diff);
    } while ((jac_iter_num < jac_maxiter) && (max_jac_diff > jac_tol));

    /* are we printing every solution on each grid only? */
      if (output_level == 2) {
	sprintf(filename, "solution.%d.%d.dat", ref_count, jac_iter_num);
	GRACEFULLY_OPEN(fp, filename, "w");
	print_leaf_quadrants(fp, root);
	GRACEFULLY_CLOSE(fp, filename);
      }


  } while ((ref_count <= ref_maxiter) && calc_error_and_refine(root,ref_tol));

  /* are we printing the solution at the end only? */
  if (output_level == 1) {
    GRACEFULLY_OPEN(fp, "solution.dat", "w");
    print_leaf_quadrants(fp, root);
    GRACEFULLY_CLOSE(fp, "solution.dat");
  }

  /* clean up */
  free_quadtree(root);

  return 0;
}
