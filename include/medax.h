#ifndef VOR3D_MEDAX_H
#define VOR3D_MEDAX_H

#include "trimesh.h"
#include "points.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Paper defaults (Dey-Zhao): theta0 = pi/8, rho0 = 8.
 * RHO_SQRT3 calibrates the paper's rho0 (tuned against umbrella-triangle
 * circumradii ~ spacing/sqrt(3)) to our R = r_sample (the spacing itself):
 * effective ratio threshold is rho0 / sqrt(3). Defined here, in one place. */
#define MEDAX_THETA0_DEFAULT (M_PI / 8.0)
#define MEDAX_RHO0_DEFAULT   8.0
#define MEDAX_RHO_EFF(rho0)  ((rho0) / 1.7320508075688772)

/* Point-sample approximation of the (inner) medial axis, optionally with its
 * connectivity ("the medial graph").
 *
 * verts.p[i] is a Voronoi vertex of the surface-sample Voronoi diagram that
 * lies on the approximate medial axis; radius[i] is the radius of the
 * corresponding (empty) medial ball, i.e. the distance from verts.p[i] to the
 * surface samples defining it. Together they sample the Medial Axis Transform;
 * this point set IS the medial axis.
 *
 * graph (the "medial graph"): NULL unless build_medial_graph was requested. It
 * is the 1-skeleton of the Voronoi diagram restricted to the inner-medial tets
 * -- a node per medial point (same indexing as verts), an edge per pair of
 * face-adjacent tets (dual Voronoi edge), weight = |verts.p[i] - verts.p[j]|
 * (arclength proxy). Every edge is a certified-inside corridor: the two tets
 * share 3 samples, so their circumspheres intersect, so the medial balls
 * overlap (|c_i - c_j| <= r_i + r_j) and the segment between the points is
 * covered by the union of the two (interior) balls.
 *
 * It is deliberately plain DATA in CSR form: node i's neighbors are
 * adj[adj_off[i] .. adj_off[i+1]); adj_w[k] is edge adj[k]'s weight; undirected,
 * so each edge is stored in both endpoints' slices and adj_off[N] == 2*n_edges.
 * vor3d only produces it; traversal/analysis (Dijkstra, centering, corridors)
 * belongs to the consumer -- see MEDIAL_AXIS_BT.md. */
typedef struct medial_graph {
    int      N;         /* number of nodes; == verts.N */
    int     *adj_off;   /* size N+1 */
    int     *adj;       /* size adj_off[N] (= 2 * n_edges) */
    double  *adj_w;     /* size adj_off[N], parallel to adj */
} s_medial_graph;

typedef struct medax {
    s_points        verts;
    double         *radius;   /* size verts.N */
    s_medial_graph *graph;    /* the medial graph; NULL unless requested */
} s_medax;

#define medax_NAN (s_medax){ points_NAN, NULL, NULL }

/* Extract the inner approximate medial axis of a closed trimesh.
 * r_sample: target surface sampling spacing. Extra sample points are scattered
 *   on the faces (area-proportional, stratified) at this spacing; <= 0 means
 *   "mesh vertices only", sensible only for meshes already dense w.r.t. local
 *   thickness. Rule of thumb: r_sample <= 1/4 of the local half-thickness.
 * theta0/rho0: filter parameters (pass the *_DEFAULT macros normally).
 * randd01/rctx: RNG callback for the scatter (same pattern as vor3d_in_bp_pds);
 *   may be NULL when r_sample <= 0.
 * TOL: duplicate-point tolerance for the internal DT (as in dt_builder).
 * EPS_DEG: degeneracy epsilon (circumcenters, normals).
 * build_medial_graph: if true, also build ma.graph (the medial graph); if
 *   false, ma.graph is NULL and only the point set (verts, radius) is produced.
 * Returns a zero struct (medax_NAN) on error. */
s_medax medax_from_trimesh(const s_trimesh *mesh, double r_sample,
                           double theta0, double rho0,
                           double (*randd01)(void*), void *rctx,
                           double TOL, double EPS_DEG, bool build_medial_graph);

void free_medax(s_medax *ma);

/* Which notion of "center" medax_center computes (the node mass / objective). The
 * first three are geodesic Frechet means over the medial graph -- argmin_i
 * sum_j w_j * d(i,j)^2, restricted to the largest connected component (d = graph
 * shortest-path distance) -- differing only in the node mass w_j. They all return
 * a medial node, hence interior with clearance ma->radius[node]. Cost is
 * dominated by computing w_j: CONTRIB builds a power (Laguerre) volume partition
 * over all balls (volume_contribution_spheres) and is by far the most expensive;
 * R3 and UNIFORM are O(N) and skip it, so they run several times faster (the
 * partition, not the descent, is the bottleneck -- see MEDIAL_AXIS_BT.md). */
typedef enum {
    /* w_j = GOVERNED VOLUME of ball j: its per-ball contribution
     * vol(B_j INTERSECT L_j) to the union of all medial balls (the Laguerre/power
     * partition), so sum_j w_j == vol(Omega) and overlaps are counted once. The
     * canonical, most accurate center; also the slowest. If the union
     * decomposition degenerates, falls back internally to R3. */
    MEDAX_CENTER_CONTRIB   = 0,
    /* w_j = r_j^3 (ball volume, overlaps double-counted). A cheap O(N) proxy for
     * CONTRIB; the center typically lands within a few percent of it. */
    MEDAX_CENTER_R3        = 1,
    /* w_j = 1 (unweighted). The pure geodesic center of the medial graph itself,
     * independent of ball sizes -- the "shape midpoint". O(N). */
    MEDAX_CENTER_UNIFORM   = 2,
    /* Not a Frechet mean: the node of LARGEST medial radius in the largest
     * component (center of the deepest inscribed ball). Needs no graph traversal
     * beyond component labeling; effectively free. out_energy receives that
     * radius (the clearance), not an energy. */
    MEDAX_CENTER_MAXRADIUS = 3,
} e_medax_center;

/* Canonical center of the medial axis; see e_medax_center for the mode/objective.
 * The result is a medial node (interior, clearance ma->radius[node]). Requires
 * ma->graph (build_medial_graph = true). Returns the node index into
 * ma->verts / ma->radius, or -1 on error.
 *
 * mode: which center to compute (see e_medax_center). MEDAX_CENTER_MAXRADIUS
 *   ignores the descent knobs below.
 *
 * Every knob takes -1 (or <=0) to use its default; out_energy may be NULL:
 *   exact_below (<0 -> 2500): largest-component size below which the exact
 *     all-pairs optimum is computed (guaranteed global, O(N^2)); at/above it,
 *     warm-started multi-start gradient descent (O(k) Dijkstra). 0 -> always
 *     descent; INT_MAX -> always exact.
 *   K (<=0 -> 3): max descent seeds (deep local-radius maxima after geodesic
 *     non-max suppression).
 *   beta (<=0 -> 2.0): NMS seed separation = beta * r_seed (larger -> fewer, more
 *     spread seeds; ~1-2 gives one seed per thick region).
 *   polish_hops (<0 -> 8): when a descent stalls at a 1-hop local minimum, search
 *     this many hops out for a lower-energy node and, if found, restart the
 *     descent there (escapes flat-basin ties; the H-ball is ~O(H^2) on the
 *     sheet-like graph, evaluated only at convergence). 0 disables it.
 *   out_energy: if non-NULL, receives sum_j w_j d(center,j)^2 for the Frechet
 *     modes, or the node clearance for MEDAX_CENTER_MAXRADIUS. */
int medax_center(const s_medax *ma, e_medax_center mode, int exact_below, int K,
                 double beta, int polish_hops, double *out_energy);

/* Diagnostic: expose the descent seeds and where each converges, for a Frechet
 * mode. Fills seeds[] with the (<=K) seed node indices (deep local-radius maxima
 * after NMS) and results[] with each seed's gradient-descent result node index;
 * both are caller arrays of size >= K. mode selects the node mass as in
 * medax_center (MEDAX_CENTER_MAXRADIUS is not a Frechet mode and is treated as
 * MEDAX_CENTER_CONTRIB here). Knobs default as in medax_center. Returns the seed
 * count, or -1 on error; medax_center keeps the lowest-energy result. */
int medax_center_seeds(const s_medax *ma, e_medax_center mode, int K, double beta,
                       int polish_hops, int *seeds, int *results);

/* Step 2 introspection (debug): return just the surface sample set that
 * medax_from_trimesh would build internally -- mesh vertices followed by the
 * area-proportional stratified face scatter at spacing r_sample. Lets the
 * scatter be visualized (e.g. medax_write_vtk on the result, radius NULL)
 * before any filtering exists. Caller owns the returned points (free_points).
 * Returns points_NAN on error. */
s_points medax_debug_sample_surface(const s_trimesh *mesh, double r_sample,
                                    double (*randd01)(void*), void *rctx,
                                    double EPS_DEG);

/* Debug/visualization: write as VTK polydata (POINTS + VERTICES cells, LINES
 * cells for the graph edges when present, POINT_DATA scalar "radius").
 * Returns 0 on success, non-zero on error. */
int medax_write_vtk(const s_medax *ma, const char *path);

#endif
