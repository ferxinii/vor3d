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

/* Canonical center of the medial axis: the volume-weighted geodesic Frechet mean
 * over the medial graph -- argmin_i sum_j w_j * d(i,j)^2 -- restricted to the
 * largest connected component (d = graph shortest-path distance). The node mass
 * w_j is the GOVERNED VOLUME of ball j: its per-ball contribution
 * vol(B_j INTERSECT L_j) to the union of all medial balls (the Laguerre/power
 * partition, via volume_contribution_spheres), so sum_j w_j == vol(Omega) and
 * ball overlaps are counted once (MEDIAL_AXIS_BT.md sec.3). If the union
 * decomposition degenerates, falls back to the ball-volume proxy w_j = r_j^3.
 * The result is a medial node, hence interior with clearance ma->radius[node].
 * Requires ma->graph (build_medial_graph = true). Returns the node index into
 * ma->verts / ma->radius, or -1 on error.
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
 *   out_energy: if non-NULL, receives sum_j w_j d(center,j)^2. */
int medax_center(const s_medax *ma, int exact_below, int K, double beta,
                 int polish_hops, double *out_energy);

/* Diagnostic: expose the descent seeds and where each converges. Fills seeds[]
 * with the (<=K) seed node indices (deep local-radius maxima after NMS) and
 * results[] with each seed's gradient-descent result node index; both are caller
 * arrays of size >= K. Knobs default as in medax_center (K, beta, polish_hops
 * take -1/<=0 for defaults). Returns the seed count, or -1 on error;
 * medax_center keeps the lowest-energy result. */
int medax_center_seeds(const s_medax *ma, int K, double beta, int polish_hops,
                       int *seeds, int *results);

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
