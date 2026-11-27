#ifndef VOR3D_DELAUNAY_H
#define VOR3D_DELAUNAY_H

#include <stdbool.h>
#include "points.h"
#include "trimesh.h"
#include "vor3d_ctx.h"

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
 * ctx (nullable): tol_dup, randd01/rctx (walk order; the result does not depend on it).
 * out_kept_idx (nullable, length pts->N): index of input point i in t.points,
 *   or -1 if dropped (near-duplicate, or redundant weighted point).
 * Failure (fewer than 4 points, degenerate input with no tetrahedron,
 *   allocation failure): zeroed result, out_kept_idx filled with -1. */
s_tetmesh delaunay_3d(const s_points *pts, const double *weights,
                      const s_vor3d_ctx *ctx, int *out_kept_idx);

/* Constrained tetrahedralization of the interior of a closed, outward-oriented
 * trimesh: interior tets only; neighbors == -1 on the surface. points: the
 * mesh points (relative order kept), then any Steiner points.
 * ctx (nullable): tol_dup.
 * Failure: zeroed result. */
s_tetmesh tetrahedralize_trimesh(const s_trimesh *mesh, const s_vor3d_ctx *ctx);

/* False iff t is the zeroed failure result (tets == NULL). */
bool tetmesh_is_valid(const s_tetmesh *t);
void free_tetmesh(s_tetmesh *t);

#endif
