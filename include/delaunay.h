#ifndef VOR3D_DELAUNAY
#define VOR3D_DELAUNAY

#include <stdbool.h>
#include "scplx.h"
#include "trimesh.h"

/* Test Delaunayness of scplx */

typedef enum delaunay_test_type {
    DELAUNAY_TEST_STRICT,
    DELAUNAY_TEST_NONSTRICT
} e_delaunay_test_type;

bool are_locally_delaunay(const s_scplx *scplx, const s_ncell *ncell, int id_opposite, 
                         e_delaunay_test_type type);
bool is_delaunay_3d(const s_scplx *scplx, e_delaunay_test_type type);



/* Build a 3D Delaunay triangulation from points (copy stored in the returned s_scplx).
 * out_Nreal: pass non-NULL only when the first *out_Nreal entries of points are "real"
 *   seeds and the rest are mirror images.  On input *out_Nreal is that count; on return
 *   it holds how many of those real seeds survived deduplication and compaction inside
 *   the DT -- which may be less than the original value when seeds lie closer than
 *   TOL_duplicates to each other.  The caller must use this updated value as Nreal when
 *   calling voronoi_from_delaunay_3d, because the compacted dt->points array packs
 *   surviving real seeds first; if any were dropped, mirrored seeds slide into positions
 *   < original_Nreal and would be mistaken for real seeds.  Pass NULL to ignore. */
s_scplx construct_dt_3d(const s_points *points, const double *weights,
                        bool keep_big_tetra, double TOL_duplicates, int *out_Nreal);




/* Incremental DT builder - allows inserting mirror points into an existing DT
 * without rebuilding from scratch.  Typical use:
 *   s_dt_builder b = dt_builder_begin(seeds, NULL, TOL);
 *   // ... inspect b.dt for neighbor info, compute mirrors ...
 *   dt_builder_extend(&b, mirrors, TOL);
 *   s_scplx dt = dt_builder_end(&b, false, &Nreal, NULL, 0);
 * The big-tetrahedron sentinels remain at indices [0..3] until dt_builder_end. */
typedef struct {
    s_scplx dt;        /* big tetra present: sentinels [0..3], seeds [4..] */
    bool   *_ignored;  /* per-point deduplication flags; size = dt.points.N */
    void   *_stack;    /* opaque: heap-allocated flip stack */
    int     _N_original; /* number of original seeds (before any extend calls) */
    s_point _bbox_min;          /* accumulated bbox min of all real points inserted so far */
    s_point _bbox_max;          /* accumulated bbox max of all real points inserted so far */
    s_point _sentinel_center;   /* centroid of the four sentinel vertices (fixed after begin) */
    double  _sentinel_inradius; /* current inscribed-sphere radius of the big tetrahedron */
    int    *_l2g_ids;  /* owned local->global id map for exact-local builds; NULL otherwise */
} s_dt_builder;

/* Phase 1: build DT of seeds, keep big tetra in place.
 * Returns zero-initialised s_dt_builder (._stack == NULL) on error. */
/* bb_min_hint / bb_max_hint: optional AABB hint (e.g. bounding-polytope AABB) that
 * expands the big-tetrahedron to safely contain future mirror points even when few
 * seeds are provided.  Pass NULL for both to use the seed bounding box only. */
s_dt_builder dt_builder_begin(const s_points *seeds, const double *weights, double TOL_dup,
                               const s_point *bb_min_hint, const s_point *bb_max_hint);

/* Exact-id build (CDT): all builder predicates run on the cdt_predicates
 * registry (registry index == scplx point index), no weights.  Requires
 * cdt_predicates_init() to have been called; clears and repopulates the
 * registry for the seeds.  Add Steiners with dt_builder_extend_lnc. */
s_dt_builder dt_builder_begin_exact(const s_points *seeds, double TOL_dup,
                                    const s_point *bb_min_hint, const s_point *bb_max_hint);

/* Exact-id build of a LOCAL cavity DT (Phase B): predicates run on the SAME
 * cdt_predicates registry as the global CDT DT, without clearing it.  seeds are
 * the cavity vertex coords; l2g_real[i] is the registry id of seed i.  The 4
 * local big-tetra sentinels are registered into scratch registry slots
 * [scratch_base .. scratch_base+3] (reserve them past the global point count).
 * Vertex ids stay local (sentinels 0..3, seeds 4+i); the seam translates.
 * bb_min_hint/bb_max_hint (NULL = seed bbox) override the big-tetra bounding
 * box: cavity fills pass an INFLATED box so the finite sentinels behave like
 * points at infinity -- near-degenerate slab cavities have sliver
 * circumspheres of radius ~L^2/h that would otherwise swallow bbox-scale
 * sentinels, leaving the local DT with no real tets at all. */
s_dt_builder dt_builder_begin_exact_local(const s_points *seeds, double TOL_dup,
                                          const int *l2g_real, int n, int scratch_base,
                                          const s_point *bb_min_hint,
                                          const s_point *bb_max_hint);

/* Phase 2: insert additional points (e.g. mirrors) into the existing DT.
 * Returns false on error. */
bool dt_builder_extend(s_dt_builder *b, const s_points *new_points, double TOL_dup);

/* Insert one implicit Steiner s*V[v1] + (1-s)*V[v2] (cdt_point_set_lnc
 * convention) into an exact-id builder; registers the exact LNC and stores the
 * rounded double in points.p.  Returns false on error. */
bool dt_builder_extend_lnc(s_dt_builder *b, int v1, int v2, double s, double TOL_dup);

/* Phase 3: remove big tetra cells, compact points, free builder state.
 * If out_Nreal != NULL, sets *out_Nreal = count of surviving original seeds.
 * If kept_idx != NULL, it is composed in place over its first N_kept_idx
 * entries: any entry already holding a valid index into the pre-compaction
 * "seeds passed to dt_builder_begin" range is rewritten to that seed's index
 * in the compacted s_scplx.points (and hence, via voronoi_from_delaunay_3d,
 * in the resulting s_vdiagram's seeds/vcells); entries already < 0 are left
 * untouched. Lets a caller-built index array (e.g. one produced by an
 * earlier filtering stage) be updated to the final mapping with no extra
 * allocation -- this reuses the remap table dt_builder_end already builds
 * for compaction. Pass kept_idx = NULL, N_kept_idx = 0 to skip. */
s_scplx dt_builder_end(s_dt_builder *b, bool keep_big_tetra, int *out_Nreal,
                       int *kept_idx, int N_kept_idx);




int N_ncells_not_big_tetra(s_scplx *scplx);  /* Only makes sense if scplx has kept big_tetra */


void extract_alpha_complex(s_scplx *scplx, bool has_big_tetra, double alpha, s_dynarray *buff_ncellPTR,
                           s_hash_table *out_faces, s_hash_table *out_edges,
                           bool *out_vertices);



/* Build a Constrained Delaunay Tetrahedralization of the interior of a closed
 * trimesh, using only the mesh boundary vertices (plus any Steiner midpoints
 * inserted on edges that are missing from the initial DT).
 * Returns an s_scplx containing only the interior tetrahedra (centroid inside
 * the trimesh). The point set may be a superset of mesh->points if Steiner
 * points were added.  Returns an empty s_scplx ({0}) on error. */
s_scplx tetrahedralize_interior_trimesh(const s_trimesh *mesh,
                                        double EPS_DEG, double TOL);

/* Like tetrahedralize_interior_trimesh, but returns the FULL tetrahedralization
 * of the convex hull of the mesh vertices (interior + exterior "pocket" tets;
 * only the ghost/sentinel tets are removed). Every returned tet carries
 * nc->interior = true iff it lies inside the trimesh. Because the complex is
 * convex, point location terminates for any query: a point inside the hull
 * lands in a tet (read its .interior flag); a point outside the hull is blocked
 * at a hull boundary face (opposite == NULL) and is therefore outside the mesh.
 * Intended as a robust, exact inside/outside oracle.
 * Returns an empty s_scplx ({0}) on error. */
s_scplx tetrahedralize_domain_flagged(const s_trimesh *mesh,
                                       double EPS_DEG, double TOL);



/* Insert a single point into an existing (sentinel-stripped) DT via Bowyer-Watson.
 * Returns the new point index on success, -1 on failure or near-duplicate. */
int scplx_insert_point(s_scplx *dt, s_point p, double TOL);


#endif
