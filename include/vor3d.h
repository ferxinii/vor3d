#ifndef VOR3D_MAIN_H
#define VOR3D_MAIN_H

#include <stdbool.h>
#include "points.h"
#include "convh.h"
#include "trimesh.h"
typedef struct random_context s_random_context;


/* ---- Cells and diagram ---- */

/* A Voronoi cell is represented as an array of convex pieces: 1 piece in a
 * convex domain, >= 1 in a non-convex one. An EMPTY cell (seed outside the
 * domain whose cell misses it, or a weighted seed dominated everywhere in the
 * domain) keeps its slot with N_pieces = 0, pieces = NULL, volume = 0 and no
 * surface. */
typedef struct vcell {
    int        seed_id;
    s_convh   *pieces;
    int        N_pieces;
    double     volume;
    s_trimesh *surface;     /* array of the cell's boundary meshes, one per connected
                               component (a non-convex cell straddling a domain fold can
                               pinch into >1 lobe). NULL/0 if extraction failed. */
    int        N_surface;   /* number of surface components */
} s_vcell;

typedef struct vdiagram {
    s_points  seeds;      /* kept seeds, compacted (see out_kept_idx in vor3d_diagram) */
    s_vcell  *vcells;     /* length seeds.N, index-aligned with seeds */
    /* Combinatorial cell adjacency (Delaunay dual of the seed DT): cells i and j
     * are adjacent iff they share a Voronoi/power facet. CSR layout: cell i's
     * neighbours are adj_neighbor[adj_offset[i] .. adj_offset[i+1]-1]. */
    int      *adj_offset;     /* size seeds.N+1 */
    int      *adj_neighbor;   /* size adj_offset[seeds.N] */
    double    domain_volume;
    double    vol_rel_err;    /* (domain_volume - sum of cell volumes) / domain_volume */
    bool      is_convex;      /* kernel that built it */
} s_vdiagram;

/* The region a diagram is clipped to. Built by the vor3d_domain_from_*
 * functions below, which also decide is_convex; the convex kernel is used iff
 * is_convex, the non-convex (CDT) kernel otherwise. */
typedef struct vor3d_domain {
    bool       is_convex;
    s_trimesh  surface;   /* boundary mesh (owns memory), always present */
    s_convh    convh;     /* convex: the domain itself; non-convex: its convex hull */
    s_point    convh_vertex_avg;   /* average of convh's vertices (strictly interior) */
    s_point    bbox_min, bbox_max; /* axis-aligned bounding box of convh */
    struct simplical_complex *cdt;  /* non-convex only (private): CDT of the
                             interior, plus flagged exterior pockets if built with
                             keep_exterior; NULL when is_convex */
    bool       cdt_has_exterior;  /* cdt was built with keep_exterior = 1 */
    double     volume;    /* convex: volume of convh; non-convex: sum of interior tets */
} s_vor3d_domain;


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

/* False iff vd is the zeroed failure result (vcells == NULL). */
bool vdiagram_is_valid(const s_vdiagram *vd);
/* Index of the cell containing x (any piece, boundary included), or -1. */
int find_inside_which_vcell(const s_vdiagram *vd, s_point x, double EPS_DEG, double TOL);
void free_vdiagram(s_vdiagram *vd);
void free_vcell(s_vcell *c);


/* ---- Plotting (TEMPORARY: to be removed or refactored onto hlibs/gplot.h, see TODO.txt) ---- */
void plot_vcell(const s_vor3d_domain *domain, const s_vcell *vcell, char *f_name, const s_point ranges[2]);
void plot_all_vcells(const s_vdiagram *vdiagram, char *f_name, const s_point ranges[2], char *view_command);
void plot_vdiagram_differentviews(const s_vdiagram *vdiagram, char *f_name, const s_point ranges[2]);

#endif
