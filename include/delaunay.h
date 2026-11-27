#ifndef VOR3D_DELAUNAY_H
#define VOR3D_DELAUNAY_H

#include <stdbool.h>
#include "points.h"
#include "trimesh.h"

typedef struct random_context s_random_context;

/* Tetrahedral mesh. Results own their memory: free with free_tetmesh.
 * Every point is a vertex of some tet. Each tet is positively oriented: exact
 * orient3d sign, the convention of signed_volume_tetra (whose float result can
 * round a near-flat tet to <= 0). neighbors[4*t+k] is the tet sharing the face
 * opposite vertex k of tet t, or -1 on the mesh boundary. */
typedef struct tetmesh {
    s_points  points;     /* compacted; relative input order kept */
    double   *weights;    /* aligned with points; NULL if unweighted */
    int       N_tets;
    int      *tets;       /* flat, 4 per tet, indices into points */
    int      *neighbors;  /* flat, 4 per tet, tet indices or -1 */
} s_tetmesh;

/* Delaunay triangulation of the convex hull of pts, or the regular (weighted
 * Delaunay) triangulation if weights != NULL (power distance |x-p|^2 - w).
 * TOL_dup: near-duplicate dropping tolerance.
 * rng: breaks point-location walk ties (NULL = fixed order).
 * out_kept_idx (nullable, length pts->N): index of input point i in t.points,
 *   or -1 if dropped (near-duplicate, or redundant weighted point).
 * Failure (fewer than 4 points, degenerate input with no tetrahedron,
 *   allocation failure): zeroed result, out_kept_idx filled with -1. */
s_tetmesh delaunay_3d(const s_points *pts, const double *weights, double TOL_dup,
                      s_random_context *rng, int *out_kept_idx);

/* Constrained tetrahedralization of the interior of a closed, outward-oriented
 * trimesh: interior tets only; neighbors == -1 on the surface. points: the
 * mesh points (relative order kept), then any Steiner points.
 * Failure: zeroed result. */
s_tetmesh tetrahedralize_trimesh(const s_trimesh *mesh, double EPS_DEG, double TOL);

/* False iff t is the zeroed failure result (tets == NULL). */
bool tetmesh_is_valid(const s_tetmesh *t);
void free_tetmesh(s_tetmesh *t);

#endif
