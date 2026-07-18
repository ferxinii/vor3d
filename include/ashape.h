#ifndef VOR3D_ASHAPE
#define VOR3D_ASHAPE

#include "points.h"
#include "trimesh.h"

/* Forward decl of the PRNG context (definition in src/random.h); see scplx.h. */
typedef struct random_context s_random_context;

/* Diagnostics filled by alpha_shape_3d (all optional; pass NULL to ignore). */
typedef struct ashape_info {
    double alpha_used;    /* the alpha actually applied (== input, or auto) */
    int    N_tets_in;     /* tets in the regularized complex (after repair) */
    int    N_promoted;    /* tets force-added by manifold repair (Phase 3) */
    int    N_components;   /* connected components of the in-set */
    int    N_pts_dropped; /* kept DT points not on/inside the surface */
    double volume;        /* total volume of the in-tets (== enclosed volume) */
} s_ashape_info;

/* Boundary surface of the REGULARIZED alpha shape of pts (unweighted).
 *
 * A tetrahedron of the Delaunay triangulation is "in" iff its squared
 * circumradius is <= alpha; the returned trimesh is the boundary of the union
 * of the in-tets ("regularized" = tets only, no dangling lower-dim simplices),
 * oriented outward, so that its enclosed volume equals the total volume of the
 * in-tets.
 *
 * Parameters:
 *   alpha   : SQUARED radius threshold.  alpha -> +inf approaches the convex
 *             hull; smaller alpha carves concavities and interior cavities.
 *             Pass alpha <= 0 to auto-select (see step 5 below).
 *   TOL_dup : near-duplicate dropping tolerance, forwarded to the DT builder.
 *   info    : optional diagnostics (may be NULL).
 *
 * Algorithm:
 *   1. Delaunay triangulation of the (deduplicated) points, with the
 *      big-tetrahedron sentinels kept, so every real face has a neighbor.
 *   2. Marking: each real tet is flagged "in" iff its squared circumradius is
 *      <= alpha.  Sentinel-touching tets are always "out".
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
 *      tets changed is reported in info->N_promoted.
 *   4. Boundary extraction: every face separating an in-tet from an out-tet is
 *      emitted, wound so its normal points away from the in-tet (outward), and
 *      the referenced vertices are compacted into a fresh trimesh.
 *   5. Auto-alpha (when the caller passes alpha <= 0): the smallest alpha, taken
 *      from the sorted real-tet circumradii, such that every point that is a
 *      vertex of some real tet is also a vertex of an in-tet.  This is the
 *      TIGHTEST wrap that still contains all points: it guarantees coverage,
 *      validity and connectivity, but the surface is deliberately jagged (it
 *      under-fills volume and may carry small handles).  For a smoother result
 *      pass a larger explicit alpha.  The value used is reported in
 *      info->alpha_used.
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
 * Returns trimesh_NAN if the points are degenerate (fewer than 4 after
 * deduplication, or coplanar so that no real tetrahedron exists), if the
 * in-set is empty at this alpha, or on allocation failure. */
s_trimesh alpha_shape_3d(const s_points *pts, double alpha, double TOL_dup,
                         s_ashape_info *info, s_random_context *rng);

#endif
