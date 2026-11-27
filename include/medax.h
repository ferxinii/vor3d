#ifndef VOR3D_MEDAX_H
#define VOR3D_MEDAX_H

#include <stdbool.h>
#include "trimesh.h"
#include "points.h"
#include "vor3d_opts.h"

/* Point-sample approximation of the (inner) medial axis, optionally with its
 * connectivity ("the medial graph").
 *
 * verts.p[i] is a Voronoi vertex of the surface-sample Voronoi diagram that
 * lies on the approximate medial axis; radius[i] is the radius of the
 * corresponding (empty) medial ball, i.e. the distance from verts.p[i] to the
 * surface samples defining it. Together they sample the Medial Axis Transform;
 * this point set IS the medial axis.
 *
 * graph (the "medial graph"): NULL unless want_graph was requested. It
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

/* Extract the inner approximate medial axis of a closed trimesh (Dey-Zhao).
 * The surface is sampled, the samples are Delaunay-triangulated, and the
 * Voronoi vertices dual to "medial" Delaunay edges are kept.
 * r_sample: target surface sampling spacing. Extra sample points are scattered
 *   on the faces (area-proportional, stratified) at this spacing; <= 0 means
 *   "mesh vertices only", sensible only for meshes already dense w.r.t. local
 *   thickness. Rule of thumb: r_sample <= 1/4 of the local half-thickness.
 * theta0, rho0: the Dey-Zhao edge filters; <= 0 selects the paper defaults
 *   (theta0 = pi/8, rho0 = 8), which is what you want normally. An edge is
 *   medial if it makes an angle > theta0 with every surface plane around its
 *   sample (angle test), or if it is long w.r.t. the local sampling, ratio
 *   > rho0 (ratio test). Smaller values keep more of the axis, including
 *   spurious "hairs" from sampling noise; larger values give a cleaner axis
 *   but may lose thin parts.
 * opts (nullable): randd01/rctx drive the scatter and are required when
 *   r_sample > 0; tol_dup merges samples in the internal DT.
 * want_graph: also build ma.graph (the medial graph); otherwise ma.graph is
 *   NULL and only the point set (verts, radius) is produced.
 * On failure returns an invalid medax (medax_is_valid false). Zero medial
 * vertices is NOT a failure. free_medax is safe on both. */
s_medax medax_from_trimesh(const s_trimesh *mesh, double r_sample,
                           double theta0, double rho0,
                           const s_vor3d_opts *opts, bool want_graph);

bool medax_is_valid(const s_medax *ma);
void free_medax(s_medax *ma);

/* Debug/visualization: write as VTK polydata (POINTS + VERTICES cells, LINES
 * cells for the graph edges when present, POINT_DATA scalar "radius").
 * Returns 0 on success, non-zero on error. */
int medax_write_vtk(const s_medax *ma, const char *path);

#endif
