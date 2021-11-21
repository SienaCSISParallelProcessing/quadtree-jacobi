/*
  quadtree functions for Jacobi on quadtree project

  Jim Teresco

  Initial implementation: Sun Mar  5 19:50:13 EST 2006

*/

#include "quadtree.h"
#include "macros.h"

/* quadtree functions */

/* create a new node */
struct quadtree *new_quadtree(double top, double left, double bottom, 
			      double right, double value,
			      struct quadtree *parent) {

  struct quadtree *treenode;
  SAFE_MALLOC(treenode,struct quadtree *,sizeof(struct quadtree));
  treenode->value = value;
  treenode->previous = value;
  treenode->top = top;
  treenode->left = left;
  treenode->bottom = bottom;
  treenode->right = right;
  treenode->parent = parent;
  treenode->children[0] = NULL;
  treenode->children[1] = NULL;
  treenode->children[2] = NULL;
  treenode->children[3] = NULL;
  return treenode;
}

/* free a node and all of its children */
void free_quadtree(struct quadtree *root) {
  int i;

  for (i=0; i<4; i++) {
    if (root->children[i] != NULL) free_quadtree(root->children[i]);
  }

  free(root);

}

/* refine a leaf quadrant, copying down solution values */
void refine_leaf_quadrant(struct quadtree *leaf) {

  /* assert that this is a leaf -- all children must be null */
  ASSERT(is_leaf_quadrant(leaf));
  
  /* create new northwest child */
  leaf->children[0] = new_quadtree(leaf->top, leaf->left,
				   (leaf->top+leaf->bottom)/2,
				   (leaf->left+leaf->right)/2,
				   leaf->value, leaf);
  /* create new northeast child */
  leaf->children[1] = new_quadtree(leaf->top, (leaf->left+leaf->right)/2,
				   (leaf->top+leaf->bottom)/2, leaf->right,
				   leaf->value, leaf);
  /* create new southwest child */
  leaf->children[2] = new_quadtree((leaf->top+leaf->bottom)/2, leaf->left,
				   leaf->bottom,(leaf->left+leaf->right)/2,
				   leaf->value, leaf);
  /* create new southeast child */
  leaf->children[3] = new_quadtree((leaf->top+leaf->bottom)/2,
				   (leaf->left+leaf->right)/2,
				   leaf->bottom, leaf->right,
				   leaf->value, leaf);
}

/* is this a leaf quadrant? */
int is_leaf_quadrant(struct quadtree *node) {

  return (node->children[0] == NULL);
}

/* visit all leaf quadrants */
void visit_all_leaf_quadrants(struct quadtree *root, 
			      QUADRANT_VISITOR_FN *func,
			      void *data) {
  int i;

  if (is_leaf_quadrant(root)) {
    (*func)(data, root);
  }
  else {
    for (i=0; i<4; i++) {
      visit_all_leaf_quadrants(root->children[i], func, data); 
    }
  }
}

/* compute the centroid x coordinate of a quadrant */
double quadrant_centroid_x(struct quadtree *node) {

  return (node->left+node->right)/2;
}

/* compute the centroid y coordinate of a quadrant */
double quadrant_centroid_y(struct quadtree *node) {

  return (node->top+node->bottom)/2;
}

/* is a given point inside this quadrant? */
int quadrant_contains(struct quadtree *node, double x, double y) {

  return ((x <= node->right) && (x >= node->left) &&
	  (y <= node->top) && (y >= node->bottom));
}

/* which child of this quadrant contains the given point? */
struct quadtree *child_quadrant_containing(struct quadtree *node, 
					   double x, double y) {

  ASSERT(quadrant_contains(node, x, y));

  if (x < quadrant_centroid_x(node)) {
    if (y < quadrant_centroid_y(node)) return node->children[2];
    return node->children[0];
  }
  /* else x > quadrant_centroid_x(node) */
  if (y < quadrant_centroid_y(node)) return node->children[3];
  return node->children[1];
}

static void print_quadrant_voidfp(void *voidfp, struct quadtree *node) {

  print_quadrant((FILE *)voidfp, node);
}

/* print the coordinates and values of a quadrant */
void print_quadrant(FILE *fp, struct quadtree *node) {

  fprintf(fp, "%12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n", 
	  node->top, node->left, node->bottom, node->right, 
	  quadrant_centroid_x(node), quadrant_centroid_y(node), 
	  node->value, node->previous);
}

/* print out the coordinates and solution values of leaf quadrants */
void print_leaf_quadrants(FILE *fp, struct quadtree *root) {

  visit_all_leaf_quadrants(root, print_quadrant_voidfp, fp);
}

/* compute the level of the given quadrant where the parent is level 0 */
int quadrant_level(struct quadtree *node) {

  int level = 0;
  while (node->parent != NULL) {
    level++;
    node = node->parent;
  }

  return level;
}

/* get the value of a quadrant */
double quadrant_value(struct quadtree *node) {

  return node->value;
}

/* get the previous value of a quadrant */
double quadrant_previous(struct quadtree *node) {

  return node->previous;
}

/* get the value of a quadrant as an average of its children
   bounding the given direction */
double quadrant_side_value(struct quadtree *node, char dir) {

  // a leaf?  easy!
  if (is_leaf_quadrant(node)) return node->value;

  // average the leaf quadrant values from the two children
  // that bound on the dir side
  double children_vals = 0.0;
  switch (dir) {
  case 'n':
    children_vals = quadrant_side_value(node->children[0], 'n') +
      quadrant_side_value(node->children[1], 'n');
    break;
  case 's':
    children_vals = quadrant_side_value(node->children[2], 's') +
      quadrant_side_value(node->children[3], 's');
    break;
  case 'e':
    children_vals = quadrant_side_value(node->children[1], 'e') +
      quadrant_side_value(node->children[3], 'e');
    break;
  case 'w':
    children_vals = quadrant_side_value(node->children[0], 'w') +
      quadrant_side_value(node->children[2], 'w');
    break;
  }
  return children_vals / 2.0;
}

/* get the previous of a quadrant as an average of its children
   bounding the given direction */
double quadrant_side_previous(struct quadtree *node, char dir) {

  // a leaf?  easy!
  if (is_leaf_quadrant(node)) return node->previous;

  // average the leaf quadrant previouss from the two children
  // that bound on the dir side
  double children_vals = 0.0;
  switch (dir) {
  case 'n':
    children_vals = quadrant_side_previous(node->children[0], 'n') +
      quadrant_side_previous(node->children[1], 'n');
    break;
  case 's':
    children_vals = quadrant_side_previous(node->children[2], 's') +
      quadrant_side_previous(node->children[3], 's');
    break;
  case 'e':
    children_vals = quadrant_side_previous(node->children[1], 'e') +
      quadrant_side_previous(node->children[3], 'e');
    break;
  case 'w':
    children_vals = quadrant_side_previous(node->children[0], 'w') +
      quadrant_side_previous(node->children[2], 'w');
    break;
  }
  return children_vals / 2.0;
}

/* set the value of a quadrant */
void quadrant_set_value(struct quadtree *node, double value) {

  node->value = value;
}

/* set the previous value of a quadrant */
void quadrant_set_previous(struct quadtree *node, double value) {

  node->previous = value;
}

/* find neighbor quadrant in the given direction (specified by
   characters 'n', 's', 'e', 'w').  Returns NULL if the neighbor would
   be outside of the root quadrant.  Returns a quadrant at most at the
   same level. */
struct quadtree *neighbor_quadrant(struct quadtree *node, char dir) {

  /* stop if we search down to this level */
  int maxlevel;
  /* coordinates that would be in a neighbor, if it exists */
  double search_x, search_y;
  struct quadtree *answer;

  maxlevel = quadrant_level(node);

  search_x = quadrant_centroid_x(node);
  search_y = quadrant_centroid_y(node);
  
  switch (dir) {
  case 'n':
    /* find +y neighbor */
    search_y = node->top + (node->top - search_y);
    break;
  case 's':
    /* find -y neighbor */
    search_y = node->bottom - (search_y - node->bottom);
    break;
  case 'e':
    /* find +x neighbor */
    search_x = node->right + (node->right - search_x);
    break;
  case 'w':
    /* find -x neighbor */
    search_x = node->left - (search_x - node->left);
    break;
  default:
    FAIL;
  }
  
  /* now we know a point in the quadrant we want to find, go looking.
     Look up the tree until we find a quadrant that contains it (if we
     get past the root, we're outside the whole thing) and then look down
     in that quadrant up to the maxlevel */
  answer = node;
  while (answer && !quadrant_contains(answer, search_x, search_y)) {
    answer = answer->parent;
  }

  /* we went off the top */
  if (!answer) return NULL;

  /* see if we can refine */
  while ((quadrant_level(answer) < maxlevel) && !is_leaf_quadrant(answer)) {
    answer = child_quadrant_containing(answer, search_x, search_y);
  }

  return answer;
}

/* mark the quadrant for refinement by setting children[3] to 1 */
void mark_quadrant_for_refinement(struct quadtree *leaf) {

  ASSERT(is_leaf_quadrant(leaf));

  leaf->children[3] = (struct quadtree *)1;
}

/* check if a quadrant is marked for refinement */
int is_quadrant_marked_for_refinement(struct quadtree *node) {

  return (node->children[3] == (struct quadtree *)1);
}

/* count the leaf quadrants */
int num_leaf_quadrants(struct quadtree *node) {

  if (is_leaf_quadrant(node)) return 1;

  return num_leaf_quadrants(node->children[0]) +
    num_leaf_quadrants(node->children[1]) +
    num_leaf_quadrants(node->children[2]) +
    num_leaf_quadrants(node->children[3]);
}

static void count_visits(void *data, struct quadtree *leaf) {

  int *count = (int *)data;

  (*count)++;
}

/* count the leaf quadrants visited by the visitor (testing only) */
int num_leaf_quadrants_by_visitor(struct quadtree *node) {

  int count = 0;

  visit_all_leaf_quadrants(node, count_visits, &count);

  return count;
}
