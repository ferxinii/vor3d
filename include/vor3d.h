#ifndef VOR3D_MAIN_H
#define VOR3D_MAIN_H

#include <stdbool.h>
#include "vdiagram.h"
#include "trimesh.h"   /* s_random_context is forward-declared in scplx.h */


/* ---- Domain ----
 *
 * The caller picks the domain; the library picks the clip kernel. A domain is
 * built once and can be reused for many diagrams (it is never copied into a
 * diagram). All builders return a zeroed domain on error (check with
 * vor3d_domain_is_valid).
 *
 * vor3d_domain_from_trimesh: closed, outward-oriented mesh (the s_trimesh
 *   constructors' invariant). The library tests convexity exactly
 *   (trimesh_is_convex: one connected component, every edge locally convex by
 *   exact orient3d). Convex -> convex kernel, no CDT is built. Non-convex ->
 *   the interior CDT is built. A near-convex mesh is non-convex; to use its hull
 *   instead, build the domain from the hull (vor3d_domain_from_points on the
 *   mesh points) -- the library never swaps the domain silently.
 *   keep_exterior (non-convex only): 0 stores interior tets only (lean);
 *   1 also stores the exterior pockets of the hull, flagged, which makes
 *   vor3d_domain_contains a single exact walk. Same diagrams either way.
 * vor3d_domain_from_convh / _from_points: always convex (convex hull of the
 *   points, built with EPS_DEG). */
s_vor3d_domain vor3d_domain_from_trimesh(const s_trimesh *mesh, double EPS_DEG,
                                         double TOL, int keep_exterior);
s_vor3d_domain vor3d_domain_from_convh(const s_convh *convh);
s_vor3d_domain vor3d_domain_from_points(const s_points *points, double EPS_DEG);
bool vor3d_domain_is_valid(const s_vor3d_domain *d);
void free_vor3d_domain(s_vor3d_domain *d);

/* True iff p is inside or on the domain (on = within TOL of the hull boundary
 * for the half-space test; TOL = 0 makes that test exact).
 *   convex: half-space test against the hull faces.
 *   non-convex: bbox + hull reject, then a point-location walk on the CDT.
 *     Exact when the domain was built with keep_exterior = 1. With the lean CDT
 *     (keep_exterior = 0) a walk that lands in a tet is still conclusive, but a
 *     walk that leaves the complex through a concavity is not, and the test
 *     falls back to the surface winding number (robust; can be wrong only for
 *     points within rounding noise of the surface). */
bool vor3d_domain_contains(const s_vor3d_domain *d, s_point p, double EPS_DEG, double TOL);


/* ---- Diagram ----
 *
 * Voronoi (weights == NULL) or power/Laguerre (weights index-aligned with
 * seeds->p, distance^2 - weight) diagram of EXACTLY the given seeds, clipped to
 * the domain. Seeds are not filtered by position: a seed outside the domain
 * still owns the part of the domain closest to it (possibly nothing -> an
 * empty cell). To keep only interior seeds, filter with vor3d_domain_contains
 * first.
 *
 * out_kept_idx (nullable, length seeds->N): out_kept_idx[i] is the index of
 *   input seed i in vd.seeds / vd.vcells, or -1 if it has no cell slot
 *   (duplicate seed, or power-redundant weighted seed). Empty cells KEEP their
 *   slot (N_pieces = 0, volume = 0). On failure every entry is -1.
 * want_surface: also build each cell's boundary trimesh(es) (vcell.surface).
 * merge_orphans (non-convex only; convex cells are always connected):
 *   reassign each disconnected cell component ("orphan", cut off across a
 *   concavity) to a neighbouring seed so every seed owns one connected region
 *   (see MERGING_ORPHANS.md). Needs the surface machinery, so it forces
 *   want_surface on.
 * rng: caller-owned PRNG breaking point-location walk ties; NULL = fixed
 *   deterministic order. Must outlive the call.
 *
 * vd.vol_rel_err is always computed; the diagram is never rejected for it.
 * Returns a zeroed s_vdiagram (vcells == NULL) on failure. */
s_vdiagram vor3d_diagram(const s_points *seeds, const double *weights,
                         const s_vor3d_domain *domain,
                         double EPS_DEG, double TOL, s_random_context *rng,
                         int *out_kept_idx,
                         bool want_surface, bool merge_orphans);

#endif
