#ifndef VOR3D_ALPHA_H
#define VOR3D_ALPHA_H

#include <stdbool.h>
#include "points.h"
#include "trimesh.h"
#include "delaunay.h"   /* s_tetmesh */
#include "vor3d_ctx.h"

/* Alpha complex K of a (weighted) point set.
 * Power distance |x-p|^2 - w (w = 0 when unweighted); alpha is a SQUARED
 * (power) radius.  Results own their memory: free with free_alpha_complex. */
typedef struct alpha_complex {
    s_points  points;    /* vertices of K, compacted; every point is in K */
    double   *weights;   /* aligned with points; NULL if unweighted */
    int       N_tets, N_faces, N_edges;
    int      *tets;      /* flat, 4 per tet, indices into points */
    int      *faces;     /* flat, 3 per face */
    int      *edges;     /* flat, 2 per edge */
    double    alpha;
} s_alpha_complex;

/* K = every simplex (tet, face, edge, vertex) of the (weighted) Delaunay
 * triangulation whose orthoball has squared radius <= alpha and is
 * empty/attached-free (Edelsbrunner).  Lists contain ALL simplices of K of that
 * dimension (not only maximal ones), so K is closed: faces of listed tets are
 * listed, etc.  Isolated vertices are points that appear in no edge.
 *
 * Order: within every simplex the vertex indices are ascending, and each list
 * is sorted lexicographically (tet orientation is NOT encoded; use
 * signed_volume_tetra if needed).  points keep the relative input order.
 *
 * weights: NULL = unweighted.  alpha: any finite value (weighted, alpha <= 0 is
 *   legitimate; alpha = 0 is the union of balls with r_i^2 = w_i).
 * ctx (nullable): tol_dup, randd01/rctx (DT walk order), forwarded to the DT build.
 * out_kept_idx (nullable, length pts->N): index of input point i in ac.points,
 *   or -1 if not in K (duplicate, redundant weighted point, or outside K).
 *
 * Failure (fewer than 4 points, NaN alpha, degenerate input -- coplanar / all
 *   duplicates, i.e. no real tetrahedron -- or allocation failure): returns a
 *   zeroed struct with ac.alpha = NAN (alpha_complex_is_valid false), and
 *   out_kept_idx filled with -1.  Empty K is NOT a failure (alpha too small):
 *   it returns points.N = N_* = 0 with ac.alpha == alpha.
 * free_alpha_complex is safe on both. */
s_alpha_complex alpha_complex_3d(const s_points *pts, const double *weights,
                                 double alpha, const s_vor3d_ctx *ctx,
                                 int *out_kept_idx);
bool alpha_complex_is_valid(const s_alpha_complex *ac);
void free_alpha_complex(s_alpha_complex *ac);


/* Alpha shape: boundary surface of the REGULARIZED alpha shape.
 * Results own their memory: free with free_alpha_shape. */
typedef struct alpha_shape {
    s_trimesh surface;        /* boundary of the REGULARIZED shape, outward, edge-manifold */
    double    alpha;          /* alpha actually used (input, or auto) */
    double    volume;         /* enclosed volume == sum of in-tet volumes */
    int       N_components;   /* connected components of the in-set */
    int       N_promoted;     /* tets changed by manifold repair */
    int       N_pts_dropped;  /* input points touching no in-tet (incl. duplicates
                                 and redundant weighted points) */
    s_tetmesh tets;           /* the in-tets (want_tets), zeroed otherwise: they
                                 tetrahedralize the shape; points start with
                                 surface.points (same indices), then the interior
                                 vertices; neighbors == -1 across the surface */
} s_alpha_shape;

/* Boundary surface of the REGULARIZED alpha shape of pts.
 *
 * A tetrahedron of the (weighted) Delaunay triangulation is "in" iff its
 * squared (power) radius is <= alpha; the returned trimesh is the boundary of
 * the union of the in-tets ("regularized" = tets only, no dangling lower-dim
 * simplices), oriented outward, so that its enclosed volume equals the total
 * volume of the in-tets.
 *
 * Parameters:
 *   weights : NULL = unweighted; else power distance |x-p|^2 - w.
 *   alpha   : SQUARED (power) radius threshold.  alpha -> +inf approaches the
 *             convex hull; smaller alpha carves concavities and interior
 *             cavities.  Pass NAN to auto-select (step 5 below).  Any other
 *             value is used as given (weighted, alpha <= 0 is legitimate).
 *   ctx     : nullable; tol_dup, randd01/rctx, forwarded to the DT build.
 *   out_kept_idx (nullable, length pts->N): index of input point i in
 *             as.surface.points, or -1 if it is not a surface vertex.
 *   want_tets: also return the in-tets in as.tets (e.g. to split the shape
 *             into convex pieces without re-tetrahedralizing the surface).
 *
 * Algorithm:
 *   1. (Weighted) Delaunay triangulation of the (deduplicated) points, with
 *      the big-tetrahedron sentinels kept, so every real face has a neighbor.
 *      Redundant weighted points (empty power cell) are dropped here.
 *   2. Marking: each real tet is flagged "in" iff its squared (power) radius
 *      is <= alpha -- the same tet test as alpha_complex_3d, so before repair
 *      the in-set is exactly the tets of the alpha complex.  Sentinel-touching
 *      tets are always "out".
 *   3. Manifold repair.  The boundary of the in-set union is always closed but
 *      need not be a 2-manifold: it can pinch along an edge.  Consider the ring
 *      of tets sharing one Delaunay edge (a full cycle, since the big tetra is
 *      kept).  Walking the ring, the surface crosses that edge with exactly two
 *      faces per maximal run of in-tets.  One in-run => 2 faces => manifold;
 *      two or more in-runs => 4, 6, ... faces => a non-manifold "hourglass"
 *      pinch where separate sheets meet along the edge.  This happens wherever
 *      the in-set is locally thin (near-tangent blobs, one-tet-thick bridges).
 *      The repair DILATES to leave exactly one in-run per edge: it fills the
 *      out-"gaps" (maximal runs of out-tets) between in-runs by promoting those
 *      out-tets to in, keeping only the largest gap open so the shape is
 *      distorted as little as possible.  A gap containing a sentinel (or a tet
 *      locked by an earlier step) cannot be filled and is the one left open; if
 *      two or more such gaps remain (a rare hull pinch) it falls back to
 *      DEMOTING all but the largest in-run and locking those tets so they are
 *      never re-added.  Dilation is preferred over carving because removing
 *      tets could expose input points.  Filling one edge perturbs neighboring
 *      rings, so the whole in-set is swept until a pass changes nothing;
 *      termination is guaranteed (the in-count only grows via promotions and
 *      the locked-count only grows via demotions, both bounded).  The number of
 *      tets changed is reported in N_promoted.
 *   4. Boundary extraction: every face separating an in-tet from an out-tet is
 *      emitted, wound so its normal points away from the in-tet (outward), and
 *      the referenced vertices are compacted into a fresh trimesh.
 *   5. Auto-alpha (when the caller passes alpha = NAN): the smallest alpha,
 *      taken from the sorted real-tet squared power radii, such that every
 *      point that is a vertex of some real tet (i.e. every non-duplicate,
 *      non-redundant point) is also a vertex of an in-tet.  Weighted, the
 *      result may be <= 0; that is valid.  This is the TIGHTEST wrap that
 *      still contains all points: it guarantees coverage, validity and
 *      connectivity, but the surface is deliberately jagged (it under-fills
 *      volume and may carry small handles).  For a smoother result pass a
 *      larger explicit alpha.  The value used is reported in the alpha field.
 *
 * Manifold guarantee / limitation: the output is EDGE-manifold (every edge has
 * exactly two incident faces).  Rare vertex-only pinches (two shells touching
 * at a single vertex, not along an edge) are NOT split and remain in the
 * output; the underlying trimesh validity check is edge-based and accepts them.
 *
 * PERFORMANCE -- input order matters: points are inserted into the DT in the
 * given order, and a spatially structured order over degenerate positions
 * (e.g. a regular lattice listed in raster order) makes every insertion land
 * on the coplanar/cospherical frontier of the previous ones, forcing the
 * exact-arithmetic predicate stage throughout (measured ~10x slower on a 74k
 * lattice).  Shuffling the input is the CALLER's responsibility; jittering
 * lattice positions by a tiny fraction of the spacing helps further.
 *
 * Failure (fewer than 4 points, degenerate input -- coplanar / all
 * duplicates --, empty in-set at this alpha, allocation): returns a zeroed
 * struct with alpha = NAN (alpha_shape_is_valid false); out_kept_idx is
 * filled with -1.  free_alpha_shape is safe on both. */
s_alpha_shape alpha_shape_3d(const s_points *pts, const double *weights,
                             double alpha, const s_vor3d_ctx *ctx,
                             int *out_kept_idx, bool want_tets);
bool alpha_shape_is_valid(const s_alpha_shape *as);
void free_alpha_shape(s_alpha_shape *as);

#endif
