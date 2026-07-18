#ifndef VOR3D_SIMPLICAL_COMPLEX_H
#define VOR3D_SIMPLICAL_COMPLEX_H

#include <stdio.h>
#include "points.h"

/* Forward decl of the PRNG context (definition in src/random.h). A complex may
 * carry a caller-owned RNG used to break point-location walk ties deterministically
 * (replaces the former process-global libc rand()); NULL => a fixed-order fallback,
 * which is still deterministic (loop-detection + brute-force backstop stay sound). */
typedef struct random_context s_random_context;

typedef struct simplical_complex {  // May in stack
    s_points points;
    double *weights;  // size point.N, May be NULL if unweighted
    int N_ncells;
    struct ncell *head;  // Linked list of ncells
    int mark_stamp;
    int walk_stamp;      // in_ncell_walk loop-detection stamp; independent of mark_stamp
                         // so the walk never collides with CDT/other mark_token users.
    int walk_hint_vid;   // vertex id near the last successful walk's result; the next
                         // walk starts at point2tet[walk_hint_vid] instead of a random
                         // tet (O(1) hops for spatially coherent insertions/queries).
                         // Any start tet is valid, so a stale hint (post-compaction id
                         // reuse) is harmless; out-of-range/NULL falls back to random.
    struct ncell **point2tet;  // [v] = one tet containing vertex v; NULL if unused
    int exact_ids;  // 0 = coordinate predicates (default); 1 = cdt_predicates by vertex id.
                    // Set only for CDT builds (see dt_predseam.h); must be zero for weighted.
    const int *l2g_ids;  // exact mode only: translate this scplx's local vertex id ->
                         // cdt_predicates registry id.  NULL = identity (the global CDT DT).
                         // Non-NULL for Phase B local cavity DTs (their ids are local).
    struct ncell_pool *pool;  // slab allocator for this complex's ncells (scplx.c).
                              // Lazily created by malloc_ncell, so (s_scplx){0} stays a
                              // valid empty complex and by-value transfers keep working.
                              // Cells live in stable blocks: pointers never move, and
                              // free_complex releases whole blocks (no per-cell free).
                              // Compile with -DNCELL_NO_POOL for plain per-cell malloc
                              // (ASAN-friendly: use-after-free of a flipped tet stays
                              // detectable instead of silently reading a recycled cell).
    s_random_context *rng;    // caller-owned PRNG for walk tie-breaking; NULL = fixed
                              // deterministic order (see the forward decl above). Set by
                              // the build entry points; propagates through by-value copies.
} s_scplx;


typedef struct ncell {  // Must live in heap
    int vertex_id[4];
    struct ncell *opposite[4];
    struct ncell *next;  // Linked list of cells
    struct ncell *prev;
    int mark_token;  // Used to mark particular ncells: mark_token == mark_stamp
    int mark_token2; // Private traversal mark for ncells_incident_face, so its
                     // flood-fill never clobbers a caller's mark_token dedup.
    int walk_token;  // Private to in_ncell_walk: == scplx.walk_stamp iff this cell
                     // was already visited in the current walk (loop detection).
    bool mask_alpha;  // if it belongs to the alpha complex (for a given alpha)
    bool in_stack;
    bool interior;   // CDT domain classification: true = tet inside the trimesh.
                     // Set only by tetrahedralize_domain_flagged (keep_exterior).
} s_ncell;


typedef struct dynarray s_dynarray;
typedef struct hash_table s_hash_table;

void free_complex(s_scplx *scplx);
void assert_point2tet(const s_scplx *scplx);
void print_ncell(const s_ncell *ncell);
void print_scomplex(const s_scplx *scplx);


int test_point_in_ncell(const s_scplx *scplx, const s_ncell *ncell, s_point query);
s_ncell *bruteforce_find_ncell_containing(const s_scplx *scplx, s_point p);
/* Locate the ncell containing the query by straight-line walk. Returns NULL when
 * the walk needs to step across a boundary (opposite==NULL) face -- i.e. the
 * query is outside the complex; this is exact for a convex complex (the only
 * kind on which walking is sound). */
s_ncell *in_ncell_walk(const s_scplx *scplx, s_point p);
s_ncell *in_ncell_walk_id(const s_scplx *scplx, int point_id);  /* walk by registered point id (seam-routed) */
void plot_dt_3d_differentviews(const s_scplx *scplx, char *f_name, s_point ranges[2]);


/* HELPERS */
s_ncell *malloc_ncell(s_scplx *scplx);            /* zeroed cell from scplx's pool */
void free_ncell(s_scplx *scplx, s_ncell *ncell);  /* return cell to scplx's pool */

int ncells_incident_face(s_scplx *scplx, s_ncell *ncell, int dim_face,
                         const int *v_localid, s_dynarray *out);
int vertex_neighbors(s_scplx *scplx, int v_global,
                     s_ncell *start_cell, int v_local,
                     int skip_below, s_dynarray *out_ids,  /* ids >= skip_below are ignored */
                     s_dynarray *scratch_cells);
void face_localid_of_adjacent_ncell(const s_ncell *ncell, int dim_face,
                                    const int v_localid[3-dim_face], int id_adjacent,
                                    int out_v_localid[3-dim_face]);
s_ncell *next_ncell_ridge_cycle(const s_ncell *ncell, int v_localid_main, int v_localid_2, 
                                int *new_v_localid_main, int *new_v_localid_2);
void walk_ridge_cycle_and_check_ncells(const s_ncell *start, const int omitted[2],
                                       bool (*check)(void *ctx, const s_ncell *), void *ctx);

double power_distance_point_vertex(const s_scplx *scplx, int vid, s_point p);

void extract_vertices_ncell(const s_scplx *scplx, const s_ncell *ncell, s_point out[4]);
void extract_weights_ncell(const s_scplx *scplx, const s_ncell *ncell, double out[4]);
void extract_vertices_and_weights_ncell(const s_scplx *scplx, const s_ncell *ncell,
                                        s_point p_out[4], double w_out[4]);
void extract_ids_face(const s_ncell *ncell, int dim_face, const int *v_localid, int *out);
void extract_vertices_face(const s_scplx *scplx, const s_ncell *ncell, int dim_face,
                           const int v_localid[3-dim_face], s_point out[dim_face+1]);
void extract_weights_face(const s_scplx *scplx, const s_ncell *ncell, int dim_face,
                          const int v_localid[3-dim_face], double out[dim_face+1]);
void extract_vertices_and_weights_face(const s_scplx *scplx, const s_ncell *ncell, 
                                       int dim_face, const int v_localid[3-dim_face],
                                       s_point p_out[dim_face+1], double w_out[dim_face+1]);


#endif
