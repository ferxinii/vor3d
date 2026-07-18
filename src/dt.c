/* Algorithm from article: 
 * Hugo Ledoux. Computing the 3D Voronoi Diagram Robustly: An Easy Explanation. 
 * But considering the possibility of weights (regular triangulation) */
#include "delaunay.h"
#include "scplx.h"
#include "points.h"
#include "gtests.h"
#include "dt_predseam.h"  /* Phase 1: id-based predicate seam (dead exact branch until Phase 2) */
#include "voronoi_predicates.h"  /* orient3d_dd for the sentinel regularity reductions */
#include "dynarray.h"     /* star enumeration in the edge-unsplit preconditions */
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>


/* orient3d of four s_point (Shewchuk convention: sign det3(a-d, b-d, c-d)). */
static int or3_pts(s_point a, s_point b, s_point c, s_point d)
{
    s_point plane[3] = {a, b, c};
    return test_orientation(plane, d);
}

/* sign det3(a-b, c-d, e-f), robust (voronoi_predicates.cpp). */
static int dd_pts(s_point a, s_point b, s_point c, s_point d, s_point e, s_point f)
{
    return orient3d_dd(a.x,a.y,a.z, b.x,b.y,b.z, c.x,c.y,c.z,
                       d.x,d.y,d.z, e.x,e.y,e.z, f.x,f.y,f.z);
}

/* --- DT_TRACE=1: verbose insertion/flip tracing for tiny debug inputs ------ */
static int dt_trace_on(void)
{
    static int cached = -1;
    if (cached < 0) cached = getenv("DT_TRACE") != NULL;
    return cached;
}

static int tet_orient_sign(const s_scplx *s, const s_ncell *nc)
{
    return or3_pts(s->points.p[nc->vertex_id[0]], s->points.p[nc->vertex_id[1]],
                   s->points.p[nc->vertex_id[2]], s->points.p[nc->vertex_id[3]]);
}

static void trace_dump_tets(const s_scplx *s, const char *tag)
{
    if (!dt_trace_on()) return;
    fprintf(stderr, "[trace] ---- tets %s ----\n", tag);
    for (const s_ncell *nc = s->head; nc; nc = nc->next) {
        int o = tet_orient_sign(s, nc);
        fprintf(stderr, "[trace]   {%d,%d,%d,%d} orient=%+d%s\n",
                nc->vertex_id[0], nc->vertex_id[1], nc->vertex_id[2],
                nc->vertex_id[3], o, o == 0 ? "  <-- FLAT" : "");
    }
    /* Full non-regularity sweep, sentinel tets included (the returned complex
     * drops them, so the test-side brute force cannot see these). */
    int nbad = 0;
    for (const s_ncell *nc = s->head; nc; nc = nc->next) {
        for (int i = 0; i < 4; i++) {
            if (!nc->opposite[i]) continue;
            if (!are_locally_delaunay(s, nc, i, DELAUNAY_TEST_NONSTRICT)) {
                nbad++;
                fprintf(stderr, "[trace]   NON-REGULAR face {%d,%d,%d} of {%d,%d,%d,%d} (opp vtx %d)\n",
                        nc->vertex_id[(i+1)%4], nc->vertex_id[(i+2)%4],
                        nc->vertex_id[(i+3)%4],
                        nc->vertex_id[0], nc->vertex_id[1], nc->vertex_id[2],
                        nc->vertex_id[3], i);
            }
        }
    }
    fprintf(stderr, "[trace]   => %d non-regular half-faces %s\n", nbad, tag);
}

/* Phase 6: local regularity when ANY of the five points (four sphere-defining
 * vertices def[0..3], or the query o) is a sentinel (id < 4).
 *
 * Derivation (PLAN_DE_PREDICATES.md; validated in exact rational arithmetic by
 * tests/validate_sentinel_reductions.py, 2400/2400):
 * Regularity is in1 = sign(D5(def, o)) * sign(orient3d(def)), where D5 is the
 * 5x5 lifted in-sphere determinant, rows [x, y, z, x^2+y^2+z^2 - w, 1].
 * Sentinels are points at infinity receding at the SAME RATE (the big-tetra
 * corners are equidistant from its centre): def[i] = C + R*u_i, R -> inf.
 * Expanding D5 by row-multilinearity, the lift column can be taken by at most
 * ONE sentinel row, so with >= 2 sentinels the leading power of R mixes one
 * term per sentinel; with equal norms |u_i| the mix collapses to determinants
 * of sentinel coordinate DIFFERENCES (s-t), and the finite weights drop out.
 * (This is why the earlier "strike one sentinel, keep the rest as ordinary
 * points" attempt was wrong: it kept only one of the leading terms.)
 *
 * With F[] = finite definers and S[] = sentinel definers, BOTH in tet order
 * (this makes all row-permutation parities cancel between D5 and orient3d --
 * validated), and dd(a,b,c,d,e,f) = sign det3(a-b, c-d, e-f):
 *
 *   o finite, ns=1 (C):  or3(F0,F1,F2,o)        * or3(F0,F1,F2,S0)
 *   o finite, ns=2 (D):  dd(F0,o, F1,o, S0,S1)  * or3(F0,F1,S0,S1)
 *   o finite, ns=3 (E):  dd(F0,o, S1,S0, S2,S0) * or3(F0,S0,S1,S2)
 *   o finite, ns=4 (X):  +1  (the initial big tetra conflicts with everything)
 *   o sent,   ns=0 (B):  -1  (a point at infinity is outside any finite sphere)
 *   o sent,   ns=1 (M1): dd(F1,F0, F2,F0, S0,o) * or3(F0,F1,F2,S0)
 *   o sent,   ns=2 (M2): -dd(F0,F1, S1,S0, o,S0) * or3(F0,F1,S0,S1)
 *   o sent,   ns=3 (M3): -or3(S1,S2,o,S0)        * or3(F0,S0,S1,S2)
 *
 * Geometrically: the degenerate orthosphere is a halfspace bounded by the
 * plane through the finite definers parallel to the sentinel chord (D: s-t)
 * or to the sentinel face (E: plane(s,t,u)).  Sentinel coordinates enter only
 * robust orientation predicates (size-stable, Finding G); the orthosphere
 * itself never sees a sentinel.  Assumes the big tetra is present (ids 0..3
 * are sentinels), which holds during construction and for complexes kept with
 * keep_big_tetra.
 *
 * Returns test_orthosphere convention: +1 inside/conflict, -1 outside.
 * Resolution ladder on ties: L1 symbolic leading term (below); L2 exact finite
 * orthosphere (orthosphere_tiebreak); L3 weight-SoS (insphere_sos).  A 0 can
 * only escape if the five points are exactly coplanar, which the flat-tet
 * guard upstream prevents. */
/* Ladder level L2: finite-orthosphere evaluation of the same 5 points, used to
 * break ties of the symbolic reduction (leading R-term exactly zero, e.g. the
 * two finite definers of Case D collinear with the query).  On a tie the true
 * sign lives in the subleading R-orders; the robust finite determinant computes
 * ALL orders exactly at the actual sentinel radius.  Finding G blow-ups are not
 * a concern here: they affect configurations whose LEADING term is nonzero, and
 * those never reach this fallback. */
static int orthosphere_tiebreak(const s_scplx *scplx, const int def[4], int o_id)
{
    s_point c[4]; double w[4];
    for (int i = 0; i < 4; i++) {
        c[i] = scplx->points.p[def[i]];
        w[i] = scplx->weights ? scplx->weights[def[i]] : 0.0;
    }
    return test_orthosphere(4, c, w, scplx->points.p[o_id],
                            scplx->weights ? scplx->weights[o_id] : 0.0);
}

/* Ladder level L3 (Phase 6f): deterministic weight-SoS when the exact finite
 * determinant D5 of the 5 points (rows def[0..3], o) is itself zero, i.e. the
 * actual configuration is genuinely degenerate (weighted-cospherical).
 *
 * Perturbation model: w_i -> w_i + eps^(2^id_i), eps -> 0+.  Distinct global
 * ids give distinct eps-powers; the SMALLEST id dominates (sentinels 0..3
 * resolve first, in fixed order).  D5 is linear in each w_j, which appears
 * only in entry (row j, lift column) with coefficient -1, so
 *     dD5/dw_j = (-1)^rr_j * Omega_j
 * with rr_j = j's 0-based row index and Omega_j = orient3d of the OTHER FOUR
 * points in row order.  Perturbing two weights strikes the lift column twice
 * (proportional rows), so ALL cross-terms vanish: the perturbed sign is
 * decided by the smallest id whose Omega is nonzero.  If every Omega is zero
 * the five points are coplanar -- pre-guarded by the flat-tet check upstream.
 * The predicate is in1 = sign(perturbed D5) * sign(orient3d(def)) as in the
 * master identity.  Validated in exact rational arithmetic against the
 * eps-perturbed determinant, incl. chain depth >= 2 and two-sided facet
 * consistency: tests/validate_sos.py (941/941). */
static int insphere_sos(const s_scplx *scplx, const int def[4], int o_id)
{
    int ids[5] = { def[0], def[1], def[2], def[3], o_id };
    s_point P[5];
    for (int i = 0; i < 4; i++) P[i] = scplx->points.p[def[i]];
    P[4] = scplx->points.p[o_id];

    int base = or3_pts(P[0], P[1], P[2], P[3]);
    if (base == 0) return 0;   /* flat definers: pre-guarded upstream */

    int order[5] = { 0, 1, 2, 3, 4 };
    for (int i = 1; i < 5; i++) {   /* insertion sort by ascending id */
        int k = order[i], j = i - 1;
        while (j >= 0 && ids[order[j]] > ids[k]) { order[j+1] = order[j]; j--; }
        order[j+1] = k;
    }

    for (int n = 0; n < 5; n++) {
        int rr = order[n];
        s_point q[4]; int m = 0;
        for (int i = 0; i < 5; i++) if (i != rr) q[m++] = P[i];
        int om = or3_pts(q[0], q[1], q[2], q[3]);
        if (om != 0) return ((rr % 2 == 0) ? 1 : -1) * om * base;
    }
    return 0;   /* all five coplanar: unreachable past the flat guard */
}

static int insphere_sentinel(const s_scplx *scplx, const int def[4], int o_id)
{
    s_point F[4], S[4];   /* Case B has 4 finite definers; Case X has 4 sentinels */
    int nf = 0, ns = 0;
    for (int i = 0; i < 4; i++) {
        if (def[i] < 4) S[ns++] = scplx->points.p[def[i]];
        else            F[nf++] = scplx->points.p[def[i]];
    }
    s_point o = scplx->points.p[o_id];

    int r = 0;
    if (o_id >= 4) {   /* query finite */
        switch (ns) {
            case 1: r = or3_pts(F[0], F[1], F[2], o)
                      * or3_pts(F[0], F[1], F[2], S[0]);         break;
            case 2: r = dd_pts(F[0],o, F[1],o, S[0],S[1])
                      * or3_pts(F[0], F[1], S[0], S[1]);         break;
            case 3: r = dd_pts(F[0],o, S[1],S[0], S[2],S[0])
                      * or3_pts(F[0], S[0], S[1], S[2]);         break;
            case 4: return 1;
            default: assert(0 && "insphere_sentinel: bad pattern"); return 0;
        }
    } else {           /* query is a sentinel */
        switch (ns) {
            case 0: return -1;
            case 1: r = dd_pts(F[1],F[0], F[2],F[0], S[0],o)
                      * or3_pts(F[0], F[1], F[2], S[0]);         break;
            case 2: r = -dd_pts(F[0],F[1], S[1],S[0], o,S[0])
                      * or3_pts(F[0], F[1], S[0], S[1]);         break;
            case 3: r = -or3_pts(S[1], S[2], o, S[0])
                      * or3_pts(F[0], S[0], S[1], S[2]);         break;
            default: assert(0 && "insphere_sentinel: bad pattern"); return 0;
        }
    }
    if (r != 0) return r;                              /* L1: symbolic */
    r = orthosphere_tiebreak(scplx, def, o_id);        /* L2: exact finite */
    if (r != 0) return r;
    return insphere_sos(scplx, def, o_id);             /* L3: weight-SoS */
}

bool are_locally_delaunay(const s_scplx *scplx, const s_ncell *ncell, int id_opposite,
                          e_delaunay_test_type type)
{
    const int *v1 = ncell->vertex_id;
    const int *v2 = ncell->opposite[id_opposite]->vertex_id;

    if (dtp_orient(scplx, v1[0], v1[1], v1[2], v1[3]) == 0 ||
        dtp_orient(scplx, v2[0], v2[1], v2[2], v2[3]) == 0) {
        return false;  /* If any is flat, return false */
    }

    int opp_face_localid;
    face_localid_of_adjacent_ncell(ncell, 2, &id_opposite, id_opposite, &opp_face_localid);
    int opp_face_vertex_id = (ncell->opposite[id_opposite])->vertex_id[opp_face_localid];

    int in1;
    if (scplx->weights) {
        /* Weighted (regular triangulation) path.  Phase 6: keep the orthosphere
         * OFF the sentinels -- when the sphere-defining tet touches a sentinel,
         * reduce to a robust ORIENTATION of the finite points (Finding G).  This
         * is what makes the RT tile the convex hull regardless of sentinel size. */
        int nsent = 0;
        for (int i = 0; i < 4; i++) if (v1[i] < 4) nsent++;
        bool o_sent = (opp_face_vertex_id < 4);

        if (nsent == 0 && !o_sent) {
            /* Case A: all five points finite -> genuine finite orthosphere;
             * exact tie (weighted-cospherical) resolved by weight-SoS. */
            s_point coords1[4]; extract_vertices_ncell(scplx, ncell, coords1);
            double weights1[4]; extract_weights_ncell(scplx, ncell, weights1);
            in1 = test_orthosphere(4, coords1, weights1,
                    scplx->points.p[opp_face_vertex_id], scplx->weights[opp_face_vertex_id]);
            if (in1 == 0)
                in1 = insphere_sos(scplx, v1, opp_face_vertex_id);
        } else {
            /* Any sentinel among the 5 points: symbolic reduction to robust
             * orientation predicates (Cases B/C/D/E/M1-M3/X, see
             * insphere_sentinel). */
            in1 = insphere_sentinel(scplx, v1, opp_face_vertex_id);
        }
    } else if (!scplx->exact_ids) {
        /* Unweighted (non-CDT) path.  Historically this used the finite
         * insphere everywhere and tolerated two flaws on degenerate input
         * (diagnosed 2026-07-13 via the intermittent example_basic "invalid
         * volume" failures, seed-deterministic repro in tests/repro_basic.c):
         *   (1) sentinel-touching tets judged with the FINITE insphere are
         *       left in unfixable strict conflict near hull slivers (Finding
         *       G), so the complex is non-Delaunay BEFORE later insertions --
         *       voiding the star-insertion correctness lemma; and
         *   (2) exact cosphericity ties (mirror-symmetric Voronoi point sets)
         *       frozen as "regular" block the repair flips.
         * Fix: same treatment as the weighted path -- symbolic sentinel
         * reductions when any of the 5 points is a sentinel (weights read as
         * 0), and the deterministic weight-SoS chain on exact ties. */
        int nsent = 0;
        for (int i = 0; i < 4; i++) if (v1[i] < 4) nsent++;

        if (nsent > 0 || opp_face_vertex_id < 4) {
            in1 = insphere_sentinel(scplx, v1, opp_face_vertex_id);
        } else {
            in1 = dtp_insphere(scplx, v1[0], v1[1], v1[2], v1[3], opp_face_vertex_id);
            if (in1 == 0)
                in1 = insphere_sos(scplx, v1, opp_face_vertex_id);
        }
    } else {
        /* CDT/exact-ids path: untouched seam. */
        in1 = dtp_insphere(scplx, v1[0], v1[1], v1[2], v1[3], opp_face_vertex_id);
    }

    switch (type) {
        case DELAUNAY_TEST_STRICT:
            if (in1 == -1) return true;
            else return false;
        case DELAUNAY_TEST_NONSTRICT:
            if (in1 != 1) return true;
            else return false;
    }
}

bool is_delaunay_3d(const s_scplx *scplx, e_delaunay_test_type type)
{
    s_ncell *current = scplx->head;
    while (current) {
        const int *v = current->vertex_id;
        if (dtp_orient(scplx, v[0], v[1], v[2], v[3]) == 0) {
            fprintf(stderr, "Flat tetra!! vids: %d %d %d %d\n",
                current->vertex_id[0], current->vertex_id[1],
                current->vertex_id[2], current->vertex_id[3]);
            return false;
        }
        for (int ii=0; ii<4; ii++) {
            if (current->opposite[ii] &&
                !are_locally_delaunay(scplx, current, ii, type)) {
                    fprintf(stderr, "Insphere failure: ncell vids %d %d %d %d, opposite vids %d %d %d %d\n",
                        current->vertex_id[0], current->vertex_id[1],
                        current->vertex_id[2], current->vertex_id[3],
                        current->opposite[ii]->vertex_id[0],
                        current->opposite[ii]->vertex_id[1],
                        current->opposite[ii]->vertex_id[2],
                        current->opposite[ii]->vertex_id[3]);
                return false;
            }
        }
        current = current->next;
    }
    return true;
}


static bool p_locally_redundant_in_ncell(const s_scplx *scplx, const s_ncell *nc, int p_id)
{   /* In unweighted triangulation, all points are NON-redundant */
    if (!scplx->weights) return false;

    /* Sentinel-touching container (PLAN_REGULAR_DT 6b, corrected against the
     * references): route through insphere_sentinel's reduced predicates, same
     * +1 == in-conflict convention, so the criterion below (redundant iff NOT
     * in conflict) is unchanged.  A p STRICTLY inside any sentinel-touching
     * tet is strictly outside the finite hull, hence extremal, hence owns an
     * unbounded power region and can never be redundant -- and for those the
     * L1 symbolic term (== alpha-mol's reduced dets: delcx.h ninf==1 ~941,
     * ninf==2 ~1099, ninf==3 ~1216) already answers in-conflict.  The ladder
     * only decides the degenerate ON-hull ties (p exactly on a hull edge or
     * face of the current point set), where true redundancy DOES depend on
     * weights: L2 evaluates the finite orthosphere with the real weights, L3
     * weight-SoS.  alpha-mol resolves those ties with weight-blind coordinate
     * SoS and tolerates the resulting flat tets; we forbid flat tets, so the
     * weight-aware drop at insertion is the consistent choice.  This replaces
     * the old blanket any-sentinel early-out (only approximate on ties). */
    int nsent = 0;
    for (int i = 0; i < 4; i++)
        if (nc->vertex_id[i] < 4) nsent++;
    if (nsent > 0)
        return insphere_sentinel(scplx, nc->vertex_id, p_id) != 1;

    s_point v[4]; double w[4];
    extract_vertices_and_weights_ncell(scplx, nc, v, w);
    /* p is redundant iff its lifted point lies STRICTLY ABOVE the lifted facet
     * of the containing tet, i.e. pi(p^, orthosphere(nc)) > 0 (Vigo et al.
     * Sec. 4).  test_orthosphere returns -1 exactly in that case (+1 == in
     * conflict, must insert).  The exact tie pi == 0 (p on the orthosphere,
     * e.g. grid data with cospherical equal-weight points) is resolved by the
     * same weight-SoS chain as regularity: p on the orthosphere is NOT above
     * the lifted hull, and its power cell can be full-dimensional, so the old
     * "tie -> drop" heuristic silently lost real points; consistent SoS keeps
     * or drops per the perturbed configuration, matching the flips. */
    int s = test_orthosphere(4, v, w, scplx->points.p[p_id],
                             scplx->weights[p_id]);
    if (s == 0)
        s = insphere_sos(scplx, nc->vertex_id, p_id);
    return (s != 1);   /* s == 0 only for a flat containing tet: keep the
                        * conservative drop in that (corrupt) case */
}



// ----------------------------------- HELPERS ---------------------------------------------
//
static int id_where_equal_int(const int *arr, int N, int entry) 
{
    for (int ii=0; ii<N; ii++) if (arr[ii] == entry) return ii;

    fprintf(stderr, "id_where_equal_int: Could not find id=%d in [", entry);
    for (int jj = 0; jj < N; jj++) fprintf(stderr, "%d%s", arr[jj], jj<N-1?",":"");
    fprintf(stderr, "]\n");
    assert(1==0);
    return -10;
}

static int inarray(const int *arr1, int N, int a)
{
    for (int ii=0; ii<N; ii++) if (arr1[ii] == a) return 1;
    return 0;
}

// ----------------------------------- STACK ---------------------------------------------

#define INIT_STACK_SIZE 100

typedef struct dstack {
    s_ncell **entry;
    int size;
    int capacity;
} s_dstack;

static s_dstack stack_create(void)
{
    s_dstack stack;
    stack.size = 0;
    stack.capacity = INIT_STACK_SIZE;
    stack.entry = malloc(stack.capacity * sizeof(s_ncell *));
    return stack;
}

static void stack_free(s_dstack *stack)
{
    free(stack->entry);
}

static int stack_push(s_dstack *stack, s_ncell *ncell)
{
    if (ncell->in_stack) return 1;  /* Avoid duplicates */

    if (stack->size == stack->capacity) { /* Expand if needed */
        stack->capacity *= 2;
        s_ncell **tmp = realloc(stack->entry, stack->capacity * sizeof(s_ncell *));
        if (!tmp) { fprintf(stderr, "delaunay.c, stack_push: realloc failed increasing stack\n"); return 0; }
        stack->entry = tmp;
    }

    stack->entry[stack->size++] = ncell;
    ncell->in_stack = true;
    return 1;
}

static s_ncell *stack_pop(s_dstack *stack)
{   
    while (stack->size > 0) {
        s_ncell *ncell = stack->entry[--stack->size];
        if (ncell->in_stack) {
            ncell->in_stack = false;
            return ncell;
        }
    }
    return NULL;
}

static void stack_remove_ncell(s_dstack *stack, s_ncell *ncell) {
    int newSize = 0;
    for (int ii = 0; ii < stack->size; ii++)
        if (stack->entry[ii] != ncell)  /* Only copy entries that are not the target. */
            stack->entry[newSize++] = stack->entry[ii];
    stack->size = newSize;
    ncell->in_stack = false;
}



// ----------------------------------- FLIPS ---------------------------------------------

static inline void vertex_ids_ncell(s_ncell *ncell, int out[4])
{
    for (int ii=0; ii<4; ii++) 
        out[ii] = ncell->vertex_id[ii];
}

static inline void opposite_pointers_ncell(s_ncell *ncell, s_ncell *out[4])
{
    for (int ii=0; ii<4; ii++)
        out[ii] = ncell->opposite[ii];
}

static inline void set_ncell_vids(s_ncell *ncell, int v1, int v2, int v3, int v4)
{
    ncell->vertex_id[0] = v1;   ncell->vertex_id[1] = v2;
    ncell->vertex_id[2] = v3;   ncell->vertex_id[3] = v4;
}

static inline void set_ncell_opposite_pointers(s_ncell *ncell, s_ncell *o1, s_ncell *o2, s_ncell *o3, s_ncell *o4)
{
    ncell->opposite[0] = o1;    ncell->opposite[1] = o2;
    ncell->opposite[2] = o3;    ncell->opposite[3] = o4;
}

static inline void p2t_update(s_scplx *sc, s_ncell *nc)
{
    if (!sc->point2tet) return;
    for (int k = 0; k < 4; k++)
        sc->point2tet[nc->vertex_id[k]] = nc;
}

static int flip14(s_scplx *scplx, s_ncell *nc1, int point_id, s_dstack *stack)
{   
    scplx->N_ncells += 3;
    s_ncell *nc2 = malloc_ncell(scplx), *nc3 = malloc_ncell(scplx), *nc4 = malloc_ncell(scplx);
    if (!nc2 || !nc3 || !nc4) return 0;

    s_ncell *opp_aux[4]; opposite_pointers_ncell(nc1, opp_aux);
    int v_aux[4]; vertex_ids_ncell(nc1, v_aux);
    
    /* Update linked-list of ncells */
    /* ---- NC1 ---------------------------- */
    /* ---- NC1 --- NC2 --- NC3 --- NC4 ---- */
    nc4->next = nc1->next; if (nc4->next) (nc4->next)->prev = nc4;  /* Last migh be null! */
    nc3->next = nc4; nc4->prev = nc3;
    nc2->next = nc3; nc3->prev = nc2;
    nc1->next = nc2; nc2->prev = nc1;

    /* Update nc1 */
    nc1->vertex_id[0] = point_id;
    nc1->opposite[1] = nc2;  nc1->opposite[2] = nc3;  nc1->opposite[3] = nc4;

    /* scplx nc2 */
    set_ncell_vids(nc2,    v_aux[0], point_id, v_aux[2], v_aux[3]);
    set_ncell_opposite_pointers(nc2,    nc1, opp_aux[1], nc3, nc4);
    if (nc2->opposite[1]) {
        int opp_id; face_localid_of_adjacent_ncell(nc2, 2, &(int){1}, 1, &opp_id); 
        nc2->opposite[1]->opposite[opp_id] = nc2;
    }

    /* scplx nc3 */
    set_ncell_vids(nc3,    v_aux[0], v_aux[1], point_id, v_aux[3]);
    set_ncell_opposite_pointers(nc3,    nc1, nc2, opp_aux[2], nc4);
    if (nc3->opposite[2]) {
        int opp_id; face_localid_of_adjacent_ncell(nc3, 2, &(int){2}, 2, &opp_id);
        nc3->opposite[2]->opposite[opp_id] = nc3;
    }

    /* scplx nc4 */
    set_ncell_vids(nc4,    v_aux[0], v_aux[1], v_aux[2], point_id);
    set_ncell_opposite_pointers(nc4,    nc1, nc2, nc3, opp_aux[3]);
    if (nc4->opposite[3]) {
        int opp_id; face_localid_of_adjacent_ncell(nc4, 2, &(int){3}, 3, &opp_id);
        nc4->opposite[3]->opposite[opp_id] = nc4;
    }

    p2t_update(scplx, nc1);
    p2t_update(scplx, nc2);
    p2t_update(scplx, nc3);
    p2t_update(scplx, nc4);

    /* Push to stack */
    if (stack) {
        if (!stack_push(stack, nc1)) return 0;
        if (!stack_push(stack, nc2)) return 0;
        if (!stack_push(stack, nc3)) return 0;
        if (!stack_push(stack, nc4)) return 0;
    }
    return 1;
}


static inline void map_vid_lid(s_ncell *ncell, int v1, int *l1, int v2, int *l2, int v3, int *l3, int v4, int *l4)
{   /* Map vertex ids to localids */
    *l1 = id_where_equal_int(ncell->vertex_id, 4, v1);
    *l2 = id_where_equal_int(ncell->vertex_id, 4, v2);
    *l3 = id_where_equal_int(ncell->vertex_id, 4, v3);
    *l4 = id_where_equal_int(ncell->vertex_id, 4, v4);
}   

static int flip23(s_scplx *scplx, s_dstack *stack, s_ncell *nc1, int opp_cell_id, int opp_face_localid, s_ncell *OUT_PTRS[3])
{   
    scplx->N_ncells += 1;
    s_ncell *nc2 = nc1->opposite[opp_cell_id], *nc3 = malloc_ncell(scplx);
    if (!nc2) return 0;

    s_ncell *nc1_opp_old[4]; opposite_pointers_ncell(nc1, nc1_opp_old);
    s_ncell *nc2_opp_old[4]; opposite_pointers_ncell(nc2, nc2_opp_old);

    /* Update linked-list of ncells */
    /* ---- NC1 ----- NC2 ---------------- */
    /* ---- NC1 ----- NC2 ----- NC3 ------ */
    nc3->next = nc2->next; if (nc3->next) nc3->next->prev = nc3;
    nc2->next = nc3; nc3->prev = nc2;

    /* scplx important indices */
    int face_vid[3]; extract_ids_face(nc1, 2, &opp_cell_id, face_vid);
    int a = face_vid[0];
    int b = face_vid[1];
    int c = face_vid[2];
    int d = nc2->vertex_id[opp_face_localid];
    int p = nc1->vertex_id[opp_cell_id];

    int nc1_id_a, nc1_id_b, nc1_id_c, nc1_id_p;
    map_vid_lid(nc1, a, &nc1_id_a, b, &nc1_id_b, c, &nc1_id_c, p, &nc1_id_p);
    int nc2_id_a, nc2_id_b, nc2_id_c, nc2_id_d;
    map_vid_lid(nc2, a, &nc2_id_a, b, &nc2_id_b, c, &nc2_id_c, d, &nc2_id_d);

     /* Update nc1 */
    nc1->vertex_id[nc1_id_c] = d;
    nc1->opposite[nc1_id_a] = nc2;
    nc1->opposite[nc1_id_b] = nc3;
    nc1->opposite[nc1_id_p] = nc2_opp_old[nc2_id_c];
    if (nc2_opp_old[nc2_id_c]) {
        int opp_aux; face_localid_of_adjacent_ncell(nc1, 2, &nc1_id_p, nc1_id_p, &opp_aux);
        nc2_opp_old[nc2_id_c]->opposite[opp_aux] = nc1;
    }

    /* Update nc2 */
    nc2->vertex_id[nc2_id_a] = p;
    nc2->opposite[nc2_id_b] = nc3;
    nc2->opposite[nc2_id_c] = nc1;
    nc2->opposite[nc2_id_d] = nc1_opp_old[nc1_id_a];
    if (nc1_opp_old[nc1_id_a]) {
        int opp_aux; face_localid_of_adjacent_ncell(nc2, 2, &nc2_id_d, nc2_id_d, &opp_aux);
        nc1_opp_old[nc1_id_a]->opposite[opp_aux] = nc2;
    }

    /* Build nc3 */
    set_ncell_vids(nc3, p, c, d, a);
    set_ncell_opposite_pointers(nc3,  nc2_opp_old[nc2_id_b], nc1, nc1_opp_old[nc1_id_b], nc2);
    if (nc1_opp_old[nc1_id_b]) {
        int opp_aux; face_localid_of_adjacent_ncell(nc3, 2, &(int){2}, 2, &opp_aux);
        nc1_opp_old[nc1_id_b]->opposite[opp_aux] = nc3;
    }
    if (nc2_opp_old[nc2_id_b]) {
        int opp_aux; face_localid_of_adjacent_ncell(nc3, 2, &(int){0}, 0, &opp_aux);
        nc2_opp_old[nc2_id_b]->opposite[opp_aux] = nc3;
    }

    p2t_update(scplx, nc1);
    p2t_update(scplx, nc2);
    p2t_update(scplx, nc3);

    if (stack) {
        if (!stack_push(stack, nc1)) return 0;
        if (!stack_push(stack, nc2)) return 0;
        if (!stack_push(stack, nc3)) return 0;
    }
    if (OUT_PTRS) { OUT_PTRS[0] = nc1; OUT_PTRS[1] = nc2; OUT_PTRS[2] = nc3; }
    return 1;
}

static int can_perform_flip32(const s_ncell *ncell, int opp_cell_id, int *ridge_id_2)
{   /* Checks if tetra abpd exists */
    int face_vid[3]; extract_ids_face(ncell, 2, &opp_cell_id, face_vid);
    int opp_face_lid; face_localid_of_adjacent_ncell(ncell, 2, &opp_cell_id, opp_cell_id, &opp_face_lid);
    
    s_ncell *opp_ncell = ncell->opposite[opp_cell_id];
    for (int ii=0; ii<4; ii++) {
        s_ncell *opp_opp = opp_ncell->opposite[ii];
        if (opp_opp && opp_opp != ncell &&
            inarray(opp_opp->vertex_id, 4, ncell->vertex_id[opp_cell_id]) &&
            inarray(opp_opp->vertex_id, 4, opp_ncell->vertex_id[opp_face_lid])) {
                *ridge_id_2 = id_where_equal_int(ncell->vertex_id, 4, opp_ncell->vertex_id[ii]);
                /* ES/Vigo flippability: a 3->2 flip is valid only if the reflex
                 * edge (ridge) has degree EXACTLY 3 -- i.e. the ring
                 * ncell->n2->n3 closes back to ncell.  Finding a third tet
                 * sharing p and d is necessary but NOT sufficient; on a degree>3
                 * ridge the facet is not flippable yet and flip32's surgery would
                 * corrupt the complex (duplicate tets).  Defer so the cascade
                 * flips other facets first (which reduces this ridge's degree). */
                int m2, s2, m3, s3, m4, s4;
                s_ncell *n2 = next_ncell_ridge_cycle(ncell, opp_cell_id, *ridge_id_2, &m2, &s2);
                if (!n2) { if (dt_trace_on()) fprintf(stderr, "[trace]       f32: ring open at n2\n"); return 0; }
                s_ncell *n3 = next_ncell_ridge_cycle(n2, m2, s2, &m3, &s3);
                if (!n3) { if (dt_trace_on()) fprintf(stderr, "[trace]       f32: ring open at n3\n"); return 0; }
                s_ncell *n4 = next_ncell_ridge_cycle(n3, m3, s3, &m4, &s4);
                if (n4 != ncell) {
                    if (dt_trace_on()) fprintf(stderr, "[trace]       f32: ridge degree > 3\n");
                    return 0;
                }
                return 1;
        }
    }
    if (dt_trace_on()) fprintf(stderr, "[trace]       f32: no third tet with p+d found\n");
    return 0;
}


static int flip32(s_scplx *scplx, s_dstack *stack, s_ncell *nc1, int opp_cell_id, int ridge_id_2, int opp_face_localid, s_ncell *OUT_PTRS[2])
{
    scplx->N_ncells -= 1;
    s_ncell *nc2, *nc3;
    {
        int v2_main, v2_2, v3_main, v3_2;
        nc2 = next_ncell_ridge_cycle(nc1, opp_cell_id, ridge_id_2, &v2_main, &v2_2);
        nc3 = next_ncell_ridge_cycle(nc2, v2_main, v2_2, &v3_main, &v3_2);
    }
    s_ncell *nc1_opp_old[4]; opposite_pointers_ncell(nc1, nc1_opp_old);
    s_ncell *nc2_opp_old[4]; opposite_pointers_ncell(nc2, nc2_opp_old);
    s_ncell *nc3_opp_old[4]; opposite_pointers_ncell(nc3, nc3_opp_old);

    /* Remove nc3 from linked list */
    s_ncell *nc3_next = nc3->next;
    if (nc3->next) nc3->next->prev = nc3->prev;
    if (nc3->prev) nc3->prev->next = nc3_next;
    else scplx->head = nc3->next;

    /* scplx important indices */
    int lid_ridge[2] = {opp_cell_id, ridge_id_2};
    int vid_ridge[2]; extract_ids_face(nc1, 1, lid_ridge, vid_ridge);
    int p = nc1->vertex_id[opp_cell_id];
    int a = nc1->vertex_id[ridge_id_2];
    int b = vid_ridge[0];
    int c = vid_ridge[1];
    int d = nc2->vertex_id[opp_face_localid];

    int nc1_id_a, nc1_id_b, nc1_id_c, nc1_id_p;
    map_vid_lid(nc1, a, &nc1_id_a, b, &nc1_id_b, c, &nc1_id_c, p, &nc1_id_p);
    int nc2_id_a, nc2_id_b, nc2_id_c, nc2_id_d;
    map_vid_lid(nc2, a, &nc2_id_a, b, &nc2_id_b, c, &nc2_id_c, d, &nc2_id_d);
    int nc3_id_b, nc3_id_c, nc3_id_d, nc3_id_p;
    map_vid_lid(nc3, b, &nc3_id_b, c, &nc3_id_c, d, &nc3_id_d, p, &nc3_id_p);

     /* Update nc1 */
    nc1->vertex_id[nc1_id_c] = d;
    nc1->opposite[nc1_id_b] = nc2;
    nc1->opposite[nc1_id_p] = nc2_opp_old[nc2_id_c];
    nc1->opposite[nc1_id_a] = nc3_opp_old[nc3_id_c];
    if (nc2_opp_old[nc2_id_c]) {
        int opp_aux; face_localid_of_adjacent_ncell(nc1, 2, &nc1_id_p, nc1_id_p, &opp_aux);
        nc2_opp_old[nc2_id_c]->opposite[opp_aux] = nc1;
    }
    if (nc3_opp_old[nc3_id_c]) {
        int opp_aux; face_localid_of_adjacent_ncell(nc1, 2, &nc1_id_a, nc1_id_a, &opp_aux);
        nc3_opp_old[nc3_id_c]->opposite[opp_aux] = nc1;
    }

    /* Update nc2 */
    nc2->vertex_id[nc2_id_b] = p; 
    nc2->opposite[nc2_id_c] = nc1;
    nc2->opposite[nc2_id_d] = nc1_opp_old[nc1_id_b];
    nc2->opposite[nc2_id_a] = nc3_opp_old[nc3_id_b];
    if (nc1_opp_old[nc1_id_b]) {
        int opp_aux; face_localid_of_adjacent_ncell(nc2, 2, &nc2_id_d, nc2_id_d, &opp_aux);
        nc1_opp_old[nc1_id_b]->opposite[opp_aux] = nc2;
    }
    if (nc3_opp_old[nc3_id_b]) {
        int opp_aux; face_localid_of_adjacent_ncell(nc2, 2, &nc2_id_a, nc2_id_a, &opp_aux);
        nc3_opp_old[nc3_id_b]->opposite[opp_aux] = nc2;
    }

    if (stack) stack_remove_ncell(stack, nc3);
    free_ncell(scplx, nc3);
    p2t_update(scplx, nc1);
    p2t_update(scplx, nc2);
    if (stack) {
        if (!stack_push(stack, nc1)) return 0;
        if (!stack_push(stack, nc2)) return 0;
    }
    if (OUT_PTRS) { OUT_PTRS[0] = nc1; OUT_PTRS[1] = nc2; }
    return 1;
}


static int can_perform_flip44(const s_scplx *scplx, const s_ncell *ncell, int opp_cell_id, int *ridge_id_2)
{   /* In general, in config44 no need for coplanarity. But since this flip is only done on a degenerate case,
       it makes the test simpler. ridge_id_2 is the other vertex NOT belonging to the ridge. */
    /* Determine ridge along which we are in config44 */
    int opp_face_localid; face_localid_of_adjacent_ncell(ncell, 2, &opp_cell_id, opp_cell_id, &opp_face_localid);
    int face_vid[3]; extract_ids_face(ncell, 2, &opp_cell_id, face_vid);
    int pid = ncell->vertex_id[opp_cell_id];
    int did = ncell->opposite[opp_cell_id]->vertex_id[opp_face_localid];

    int o0 = dtp_orient(scplx, pid, face_vid[0], face_vid[1], did);
    int o1 = dtp_orient(scplx, pid, face_vid[1], face_vid[2], did);
    int o2 = dtp_orient(scplx, pid, face_vid[2], face_vid[0], did);
    /* Expected: 1 or 2 edges coplanar with pd (ridge through edge or vertex). */
    assert((o0==0)+(o1==0)+(o2==0) >= 1 && (o0==0)+(o1==0)+(o2==0) <= 2);

    /* NOTE (2026-07-18): TWO zero orients <=> a shared-face vertex r is
     * COLLINEAR with p and d, strictly between them.  On the WEIGHTED path
     * that sub-case is dispatched to the atomic edge-unsplit BEFORE flip44 is
     * consulted (flip_tetrahedra CASE_FLAT; PLAN_UNSPLIT_EDGE.md) -- the old
     * behavior of running flip44 there built tets around edge (p,d) THROUGH r
     * (transient flat tets) and deadlocked in 26/120 fixture-D orders.  The
     * unweighted path still reaches this function with 2 zeros and keeps the
     * legacy behavior (the trigger never pops there in practice). */

    int AUX_ridge_id_2[2];  int k=0;
    if (o0 == 0) { AUX_ridge_id_2[k++] = id_where_equal_int(ncell->vertex_id, 4, face_vid[2]); }
    if (o1 == 0) { AUX_ridge_id_2[k++] = id_where_equal_int(ncell->vertex_id, 4, face_vid[0]); }
    if (o2 == 0) { AUX_ridge_id_2[k++] = id_where_equal_int(ncell->vertex_id, 4, face_vid[1]); }

    for (int i=0; i<k; i++) {
        *ridge_id_2 = AUX_ridge_id_2[i];

        int nc2_id1, nc2_id2;
        s_ncell *nc2 = next_ncell_ridge_cycle(ncell, opp_cell_id, *ridge_id_2, &nc2_id1, &nc2_id2);
        if (!nc2) continue;
        int nc3_id1, nc3_id2;
        s_ncell *nc3 = next_ncell_ridge_cycle(nc2, nc2_id1, nc2_id2, &nc3_id1, &nc3_id2);
        if (!nc3) continue;
        int nc4_id1, nc4_id2;
        s_ncell *nc4 = next_ncell_ridge_cycle(nc3, nc3_id1, nc3_id2, &nc4_id1, &nc4_id2);
        if (!nc4) continue;

        if (nc4->opposite[nc4_id1] == ncell) {
            /* Duplicate-avoidance guard (local, O(4)): flip44's internal flip23
             * would replace {ncell, its neighbour} with three tets
             * {p,d}+(each edge of the shared face {a,c,b}).  Each such tet shares
             * a face CONTAINING p with ncell, so it already exists iff one of
             * ncell's neighbours across a face containing p (any face except the
             * shared one, opp_cell_id) also contains d (= did).  If so the flip
             * would duplicate that cell and corrupt the complex (same risk
             * CASE_P_IN_EDGE guards for flip23).  Defer: try the other candidate
             * ridge, else fall through to return 0 so the BW loop flips
             * elsewhere first. */
            bool would_dup = false;
            for (int f = 0; f < 4 && !would_dup; f++) {
                if (f == opp_cell_id) continue;         /* the shared face (opposite p) */
                const s_ncell *nb = ncell->opposite[f];
                if (nb && inarray(nb->vertex_id, 4, did)) would_dup = true;
            }
            if (would_dup) {
                if (dt_trace_on()) fprintf(stderr, "[trace]       f44: candidate ridge would duplicate a tet\n");
                continue;
            }
            return 1;
        }
        if (dt_trace_on()) fprintf(stderr, "[trace]       f44: candidate ridge ring did not close in 4\n");
    }
    return 0;
}


static int flip_tetrahedra(s_scplx *scplx, s_dstack *stack, s_ncell *ncell, int opp_cell_id, bool *ignored);

static int flip44(s_scplx *scplx, s_dstack *stack, s_ncell *ncell, int id_ridge_1, int id_ridge_2, s_ncell *OUT_PTRS[4], bool *ignored)
{
    /* 1) flip23 */
    /* Store some indices before flip */
    int opp_face_lid; face_localid_of_adjacent_ncell(ncell, 2, &id_ridge_1, id_ridge_1, &opp_face_lid);
    int opp_face_vid = ncell->opposite[id_ridge_1]->vertex_id[opp_face_lid];
    int lid_ridge[2] = {id_ridge_1, id_ridge_2};
    int vid_ridge[2]; extract_ids_face(ncell, 1, lid_ridge, vid_ridge);
    int p = ncell->vertex_id[id_ridge_1];
    int a = vid_ridge[0];
    int c = vid_ridge[1];
    int d = opp_face_vid;

    s_ncell *FLIP23_PTRS[3];
    if (!flip23(scplx, NULL, ncell, id_ridge_1, opp_face_lid, FLIP23_PTRS)) return 0;

    /* Find which ncell added shares ridge. Currently, no way to predict this. */
    s_ncell *nc5 = NULL;
    for (int ii=0; ii<3; ii++) {
        if (inarray(FLIP23_PTRS[ii]->vertex_id, 4, a) && inarray(FLIP23_PTRS[ii]->vertex_id, 4, c) &&
            inarray(FLIP23_PTRS[ii]->vertex_id, 4, d) && inarray(FLIP23_PTRS[ii]->vertex_id, 4, p)) {
            nc5 = FLIP23_PTRS[ii];
            if (stack) {
                if (!stack_push(stack, FLIP23_PTRS[(ii+1)%3])) return 0;
                if (!stack_push(stack, FLIP23_PTRS[(ii+2)%3])) return 0;
            }
            if (OUT_PTRS) { OUT_PTRS[0] = FLIP23_PTRS[(ii+1)%3]; OUT_PTRS[1] = FLIP23_PTRS[(ii+2)%3]; }
            break;
        }
    }
    assert(nc5 != NULL && "Could not perform flip44...");

    /* 2) flip32 */
    int nc5_p = id_where_equal_int(nc5->vertex_id, 4, ncell->vertex_id[id_ridge_1]);
    s_ncell *nc3 = nc5->opposite[id_where_equal_int(nc5->vertex_id, 4, ncell->vertex_id[id_ridge_1])];
    int nc3_id1 = id_where_equal_int(nc3->vertex_id, 4, opp_face_vid);
    int nc3_id2; face_localid_of_adjacent_ncell(nc5, 2, &nc5_p, nc5_p, &nc3_id2);
    int nc3_opp_face_lid; face_localid_of_adjacent_ncell(nc3, 2, &nc3_id1, nc3_id1, &nc3_opp_face_lid);

    s_ncell *FLIP32_PTRS[2];
    flip32(scplx, NULL, nc3, nc3_id1, nc3_id2, nc3_opp_face_lid, FLIP32_PTRS);
    if(stack) {
        if (!stack_push(stack, FLIP32_PTRS[0])) return 0;
        if (!stack_push(stack, FLIP32_PTRS[1])) return 0;

        /* The Bowyer-Watson loop pops each cell and checks exactly one face:
         * opposite[local_p], i.e. the face opposite the newly inserted point p.
         * This is correct for flip14/flip23/flip32, where every face shared between
         * two cells that both contain p is guaranteed locally Delaunay by construction
         * -- those are "interior" star faces that never need re-examination.
         *
         * flip44 breaks that invariant.  It calls flip23 then flip32 on a
         * sub-configuration determined by the CASE_FLAT ridge, not the outer
         * boundary of the star.  The two cells produced by flip32 share a face
         * that CONTAINS p (so it is opposite some old vertex, not opposite p).
         * That face was never an outer boundary candidate at any point in the
         * sequence, so no flip ever pushed it for checking.  When the BW loop
         * pops these two cells it checks their outer faces (opposite p), which
         * are fine -- and the interior shared face goes unexamined.
         *
         * Fix: check that shared face explicitly here, right after flip32 returns.
         * Use NONSTRICT so that coplanar inputs (where the tested point lies exactly
         * on the circumsphere) are accepted without flipping -- otherwise the two
         * valid diagonalisations of a degenerate square cycle forever. */
        for (int fi = 0; fi < 4; fi++) {
            if (FLIP32_PTRS[0]->opposite[fi] == FLIP32_PTRS[1]) {
                if (!are_locally_delaunay(scplx, FLIP32_PTRS[0], fi, DELAUNAY_TEST_NONSTRICT)) {
                    stack_remove_ncell(stack, FLIP32_PTRS[0]);
                    stack_remove_ncell(stack, FLIP32_PTRS[1]);
                    if (flip_tetrahedra(scplx, stack, FLIP32_PTRS[0], fi, ignored) == -1) return -1;
                }
                break;
            }
        }
    }
    if (OUT_PTRS) { OUT_PTRS[2] = FLIP32_PTRS[0]; OUT_PTRS[3] = FLIP32_PTRS[1]; }
    return 1;
}


/* ==========================================================================
 * Edge-unsplit (2k -> k): atomic removal of a vertex r lying exactly ON the
 * segment between two vertices p and d.  See PLAN_UNSPLIT_EDGE.md.
 *
 * Precondition (validated by can_perform_unsplit_edge): r's star is exactly
 * the "double fan"
 *     P_i = {r, p, a_i, a_i+1}   (ring around edge (r,p), i = 0..k-1 cyclic)
 *     D_i = {r, d, a_i, a_i+1}   (ring around edge (r,d))
 * with one common link ring a_0..a_k-1 and P_i glued to D_i across the
 * interior face {r, a_i, a_i+1}.  The operation replaces the 2k star tets by
 *     U_i = {p, d, a_i, a_i+1}
 * removing r.  Each U_i is the exact union of P_i and D_i (r on segment
 * (p,d)), so volume and the star's outer boundary are preserved.
 * ========================================================================== */

#define UNSPLIT_MAX_RING 64

typedef struct {
    int k;
    int r_vid, p_vid, d_vid;
    s_ncell *P[UNSPLIT_MAX_RING];   /* fan around (r,p); P[i] spans (a_i, a_i+1) */
    s_ncell *D[UNSPLIT_MAX_RING];   /* fan around (r,d); D[i] = P[i]'s fan mate  */
    int a[UNSPLIT_MAX_RING + 1];    /* link ring; a[k] == a[0] (walk closure)    */
} s_unsplit_cfg;

/* Walk the closed ring of tets around the edge of `start` NOT containing the
 * vertices at local slots lid_main/lid_sec, collecting cells and the cyclic
 * off-ridge vertex sequence (aseq[0] = vids[lid_main], aseq[1] = vids[lid_sec],
 * then one new vertex per step; aseq[k] closes back to aseq[0]).  Returns the
 * ring size k, or 0 on an open ring, cap overflow, or inconsistency. */
static int unsplit_collect_ring(s_ncell *start, int lid_main, int lid_sec,
                                s_ncell *cells[UNSPLIT_MAX_RING],
                                int aseq[UNSPLIT_MAX_RING + 1])
{
    s_ncell *cur = start;
    int main_l = lid_main, sec_l = lid_sec;
    cells[0] = start;
    aseq[0] = start->vertex_id[lid_main];
    aseq[1] = start->vertex_id[lid_sec];
    int k = 1;
    while (1) {
        int nmain, nsec;
        s_ncell *next = next_ncell_ridge_cycle(cur, main_l, sec_l, &nmain, &nsec);
        if (!next) return 0;                   /* open ring / broken adjacency */
        if (next == start) break;
        if (k >= UNSPLIT_MAX_RING) return 0;
        cells[k] = next;
        aseq[k + 1] = next->vertex_id[nsec];
        k++;
        cur = next; main_l = nmain; sec_l = nsec;
    }
    if (aseq[k] != aseq[0]) return 0;          /* walk closure consistency */
    return k;
}

static int can_perform_unsplit_edge(s_scplx *scplx, s_ncell *ncell, int opp_cell_id,
                                    int r_vid, int d_vid, s_unsplit_cfg *cfg)
{
    /* P0: never remove a sentinel (p or d MAY be sentinels; only r goes). */
    if (r_vid < 4) return 0;

    cfg->r_vid = r_vid;
    cfg->p_vid = ncell->vertex_id[opp_cell_id];
    cfg->d_vid = d_vid;

    s_ncell *dstart = ncell->opposite[opp_cell_id];
    if (!dstart) return 0;

    /* ncell = {r, p, x, y}: the two ring vertices are the shared-face
     * vertices other than r. */
    int lid_x = -1, lid_y = -1;
    for (int i = 0; i < 4; i++) {
        int v = ncell->vertex_id[i];
        if (v == cfg->p_vid || v == r_vid) continue;
        if (lid_x < 0) lid_x = i; else lid_y = i;
    }
    if (lid_x < 0 || lid_y < 0) return 0;

    /* P1: ring around edge (r,p). */
    int k = unsplit_collect_ring(ncell, lid_x, lid_y, cfg->P, cfg->a);
    if (k < 3) {
        if (dt_trace_on()) fprintf(stderr, "[trace]       unsplit: (r,p) ring open/degenerate (k=%d)\n", k);
        return 0;
    }
    cfg->k = k;
    for (int i = 0; i < k; i++)
        if (!inarray(cfg->P[i]->vertex_id, 4, r_vid) ||
            !inarray(cfg->P[i]->vertex_id, 4, cfg->p_vid)) return 0;

    /* Link ring sanity: distinct vertices, none of them r, p or d. */
    for (int i = 0; i < k; i++) {
        if (cfg->a[i] == r_vid || cfg->a[i] == cfg->p_vid || cfg->a[i] == d_vid)
            return 0;
        for (int j = i + 1; j < k; j++)
            if (cfg->a[i] == cfg->a[j]) return 0;
    }

    /* P2: ring around edge (r,d), started at the fan mate with the SAME
     * (x,y) orientation, so a matching double fan yields the identical
     * vertex sequence (both walks begin by crossing the face omitting x). */
    int dlid_x = -1, dlid_y = -1;
    for (int i = 0; i < 4; i++) {
        if (dstart->vertex_id[i] == ncell->vertex_id[lid_x]) dlid_x = i;
        if (dstart->vertex_id[i] == ncell->vertex_id[lid_y]) dlid_y = i;
    }
    if (dlid_x < 0 || dlid_y < 0) return 0;

    int dseq[UNSPLIT_MAX_RING + 1];
    int kd = unsplit_collect_ring(dstart, dlid_x, dlid_y, cfg->D, dseq);
    if (kd != k) {
        if (dt_trace_on()) fprintf(stderr, "[trace]       unsplit: ring mismatch k_p=%d k_d=%d\n", k, kd);
        return 0;
    }
    for (int i = 0; i < k; i++)
        if (!inarray(cfg->D[i]->vertex_id, 4, r_vid) ||
            !inarray(cfg->D[i]->vertex_id, 4, d_vid)) return 0;

    /* P3: identical link sequences ... */
    for (int i = 0; i <= k; i++)
        if (dseq[i] != cfg->a[i]) {
            if (dt_trace_on()) fprintf(stderr, "[trace]       unsplit: link sequences differ\n");
            return 0;
        }
    /* ... and P_i really is glued to D_i across {r, a_i, a_i+1} (the face of
     * P_i omitting p). */
    for (int i = 0; i < k; i++) {
        int lid_p = id_where_equal_int(cfg->P[i]->vertex_id, 4, cfg->p_vid);
        if (cfg->P[i]->opposite[lid_p] != cfg->D[i]) {
            if (dt_trace_on()) fprintf(stderr, "[trace]       unsplit: fan mates not glued pairwise\n");
            return 0;
        }
    }

    /* P4: star completeness -- r's FULL star must be exactly {P_i} u {D_i}.
     * Any other tet containing r falsifies the double-fan hypothesis. */
    {
        int lid_r = id_where_equal_int(ncell->vertex_id, 4, r_vid);
        int omit[3], m = 0;
        for (int j = 0; j < 4; j++) if (j != lid_r) omit[m++] = j;
        s_dynarray star = dynarray_initialize(sizeof(s_ncell *), 64);
        if (!star.items) return 0;
        int ok = ncells_incident_face(scplx, ncell, 0, omit, &star) &&
                 star.N == (unsigned)(2 * k);
        if (ok) {
            s_ncell **sc = star.items;
            for (unsigned s = 0; s < star.N && ok; s++) {
                bool found = false;
                for (int i = 0; i < k && !found; i++)
                    found = (sc[s] == cfg->P[i] || sc[s] == cfg->D[i]);
                ok = found;
            }
        }
        dynarray_free(&star);
        if (!ok) {
            if (dt_trace_on()) fprintf(stderr, "[trace]       unsplit: star of %d is not the double fan\n", r_vid);
            return 0;
        }
    }

    /* P5: no replacement tet may be degenerate. */
    for (int i = 0; i < k; i++)
        if (dtp_orient(scplx, cfg->p_vid, d_vid, cfg->a[i], cfg->a[i + 1]) == 0) {
            if (dt_trace_on()) fprintf(stderr, "[trace]       unsplit: replacement tet %d would be flat\n", i);
            return 0;
        }

    return 1;
}

/* Locate the current cell containing face {v0,v1,v2}; out_lid = local id of
 * its fourth vertex.  Looked up fresh through point2tet (falling back to a
 * full scan), so callers may hold NO cell pointers across cascade flips.
 * Returns 0 if the face no longer exists. */
static int unsplit_find_face(s_scplx *scplx, int v0, int v1, int v2,
                             s_ncell **out_cell, int *out_lid)
{
    s_ncell *found = NULL;
    s_ncell *seed = (scplx->point2tet && v0 >= 0 && v0 < scplx->points.N)
                        ? scplx->point2tet[v0] : NULL;
    if (seed && inarray(seed->vertex_id, 4, v0)) {
        int lid_v0 = id_where_equal_int(seed->vertex_id, 4, v0);
        int omit[3], m = 0;
        for (int j = 0; j < 4; j++) if (j != lid_v0) omit[m++] = j;
        s_dynarray star = dynarray_initialize(sizeof(s_ncell *), 64);
        if (!star.items) return 0;
        if (ncells_incident_face(scplx, seed, 0, omit, &star)) {
            s_ncell **sc = star.items;
            for (unsigned s = 0; s < star.N; s++)
                if (inarray(sc[s]->vertex_id, 4, v1) &&
                    inarray(sc[s]->vertex_id, 4, v2)) { found = sc[s]; break; }
        }
        dynarray_free(&star);
    } else {
        for (s_ncell *nc = scplx->head; nc; nc = nc->next)
            if (inarray(nc->vertex_id, 4, v0) && inarray(nc->vertex_id, 4, v1) &&
                inarray(nc->vertex_id, 4, v2)) { found = nc; break; }
    }
    if (!found) return 0;
    for (int j = 0; j < 4; j++) {
        int v = found->vertex_id[j];
        if (v != v0 && v != v1 && v != v2) { *out_cell = found; *out_lid = j; return 1; }
    }
    return 0;
}

static int unsplit_edge(s_scplx *scplx, s_dstack *stack, const s_unsplit_cfg *cfg,
                        bool *ignored)
{   /* 0 ERROR, 1 OK (flip41 convention) */
    int k = cfg->k, r = cfg->r_vid, p = cfg->p_vid, d = cfg->d_vid;

    if (dt_trace_on())
        fprintf(stderr, "[trace]     -> unsplit_edge: remove %d on segment (%d,%d), ring k=%d\n",
                r, p, d, k);

    /* S1: capture ALL outer neighbours before any rewrite (flip41 discipline).
     * X_i / Y_i contain no r, so they are never a P_j or D_j. */
    s_ncell *X[UNSPLIT_MAX_RING], *Y[UNSPLIT_MAX_RING];
    for (int i = 0; i < k; i++) {
        int lid_r_P = id_where_equal_int(cfg->P[i]->vertex_id, 4, r);
        int lid_r_D = id_where_equal_int(cfg->D[i]->vertex_id, 4, r);
        Y[i] = cfg->P[i]->opposite[lid_r_P];   /* across {p, a_i, a_i+1} */
        X[i] = cfg->D[i]->opposite[lid_r_D];   /* across {d, a_i, a_i+1} */
    }

    /* S2: reuse P_i as U_i = {p, d, a_i, a_i+1} (slots 0:p 1:d 2:a_i 3:a_i+1). */
    for (int i = 0; i < k; i++)
        set_ncell_vids(cfg->P[i], p, d, cfg->a[i], cfg->a[i + 1]);

    /* S3: gluing.  Face omitting slot 0 (p) = {d,a_i,a_i+1} -> X_i;
     * omitting 1 (d) = {p,a_i,a_i+1} -> Y_i; omitting 2 (a_i) = {p,d,a_i+1}
     * -> U_i+1; omitting 3 (a_i+1) = {p,d,a_i} -> U_i-1. */
    for (int i = 0; i < k; i++)
        set_ncell_opposite_pointers(cfg->P[i], X[i], Y[i],
                                    cfg->P[(i + 1) % k], cfg->P[(i + k - 1) % k]);

    /* S4: back-pointers (after every U's vids are final). */
    for (int i = 0; i < k; i++) {
        s_ncell *U = cfg->P[i];
        if (X[i]) {
            int lid = 0, aux;
            face_localid_of_adjacent_ncell(U, 2, &lid, lid, &aux);
            X[i]->opposite[aux] = U;
        }
        if (Y[i]) {
            int lid = 1, aux;
            face_localid_of_adjacent_ncell(U, 2, &lid, lid, &aux);
            Y[i]->opposite[aux] = U;
        }
    }

    /* S5: unlink and free the D fan. */
    for (int i = 0; i < k; i++) {
        s_ncell *Dc = cfg->D[i];
        if (Dc->next) Dc->next->prev = Dc->prev;
        if (Dc->prev) Dc->prev->next = Dc->next;
        else scplx->head = Dc->next;
        if (stack) stack_remove_ncell(stack, Dc);
        free_ncell(scplx, Dc);
    }

    /* S6: bookkeeping. */
    scplx->N_ncells -= k;
    ignored[r] = true;
    if (scplx->point2tet) scplx->point2tet[r] = NULL;
    for (int i = 0; i < k; i++) p2t_update(scplx, cfg->P[i]);

    /* S7: stack. */
    if (stack)
        for (int i = 0; i < k; i++)
            if (!stack_push(stack, cfg->P[i])) return 0;

    /* Post-surgery re-checks: the BW pop invariant does not cover these faces
     * (and unsplit can fire from nested contexts where p is not the current
     * insertion point).  Faces are recorded BY VERTEX IDS and re-located
     * through point2tet before each check, so cascade flips triggered here can
     * never leave a dangling pointer in the worklist. */
    int nfaces = 0;
    int (*faces)[3] = malloc(sizeof(int[3]) * (size_t)(3 * k));
    if (!faces) return 0;
    for (int i = 0; i < k; i++) {
        faces[nfaces][0] = d; faces[nfaces][1] = cfg->a[i]; faces[nfaces][2] = cfg->a[i + 1]; nfaces++;
        faces[nfaces][0] = p; faces[nfaces][1] = cfg->a[i]; faces[nfaces][2] = cfg->a[i + 1]; nfaces++;
        faces[nfaces][0] = p; faces[nfaces][1] = d;         faces[nfaces][2] = cfg->a[i];     nfaces++;
    }
    for (int f = 0; f < nfaces; f++) {
        s_ncell *cell; int lid;
        if (!unsplit_find_face(scplx, faces[f][0], faces[f][1], faces[f][2], &cell, &lid))
            continue;                      /* face consumed by a cascade flip */
        if (!cell->opposite[lid]) continue;
        if (!are_locally_delaunay(scplx, cell, lid, DELAUNAY_TEST_NONSTRICT)) {
            if (flip_tetrahedra(scplx, stack, cell, lid, ignored) == -1) {
                free(faces);
                return 0;
            }
        }
    }
    free(faces);
    return 1;
}


static int can_perform_flip41(const s_ncell *ncell, int opp_cell_id, int *redundant_localid)
{
    int face_vid[3]; extract_ids_face(ncell, 2, &opp_cell_id, face_vid);
    int opp_face_lid; face_localid_of_adjacent_ncell(ncell, 2, &opp_cell_id, opp_cell_id, &opp_face_lid);
    s_ncell *opp_ncell = ncell->opposite[opp_cell_id];
    int p_vid = ncell->vertex_id[opp_cell_id];

    /* For each edge of shared face, check if degree 3 tetrahedron exists */
    int degree3_edges[3][2], n_degree3 = 0;
    for (int i=0; i<3; i++) {
        int va = face_vid[i], vb = face_vid[(i+1)%3];
        for (int j=0; j<4; j++) {
            s_ncell *nb = opp_ncell->opposite[j];
            if (nb && nb != ncell &&
                inarray(nb->vertex_id, 4, p_vid) &&
                inarray(nb->vertex_id, 4, va)    &&
                inarray(nb->vertex_id, 4, vb)) {
                degree3_edges[n_degree3][0] = va;
                degree3_edges[n_degree3][1] = vb;
                n_degree3++;
                break;
            }
        }
    }

    if (n_degree3 != 2) {
        if (dt_trace_on()) fprintf(stderr, "[trace]       f41: n_degree3=%d (need 2)\n", n_degree3);
        return 0;
    }

    /* Redundant vertex is common to both degree-3 edges */
    int redundant_vid = -1;
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            if (degree3_edges[0][i] == degree3_edges[1][j])
                redundant_vid = degree3_edges[0][i];

    assert(redundant_vid != -1);
    /* Finding D: never flip away a big-tetra sentinel (id < 4).  Near the hull
     * the degree-3 pattern can land on a sentinel; removing it corrupts the hull
     * topology.  Reject the flip41 -- flip_tetrahedra falls through to flip32. */
    if (redundant_vid < 4) {
        if (dt_trace_on()) fprintf(stderr, "[trace]       f41: redundant vid %d is a sentinel\n", redundant_vid);
        return 0;
    }
    *redundant_localid = id_where_equal_int(ncell->vertex_id, 4, redundant_vid);

    /* The two degree-3 edges are NECESSARY but not SUFFICIENT for a 4->1 flip.
     * A hull-adjacent weighted configuration can pass the edge test yet have the
     * redundant vertex's star NOT be a clean 4-tet diamond (the 4 incident tets
     * span 5+ vertices).  flip41's pointer surgery assumes a diamond, so a false
     * positive corrupts the complex (stack overflow in face_localid_of_adjacent
     * _ncell writing multiple ids into a single-int slot).  Build the star and
     * require exactly 4 incident tets whose vertices, minus r, are exactly 4
     * distinct.  If not, reject -- flip_tetrahedra falls through to flip32. */
    const s_ncell *star[4]; star[0] = ncell; int k = 1;
    for (int i = 0; i < 4; i++) {
        if (i == *redundant_localid) continue;
        const s_ncell *nb = ncell->opposite[i];
        if (nb && inarray(nb->vertex_id, 4, redundant_vid) && k < 4) star[k++] = nb;
    }
    if (k != 4) {
        if (dt_trace_on()) fprintf(stderr, "[trace]       f41: star of %d has %d local tets (need 4)\n", redundant_vid, k);
        return 0;
    }
    int diamond[4], nd = 0;
    for (int i = 0; i < 4; i++)
        for (int m = 0; m < 4; m++) {
            int v = star[i]->vertex_id[m];
            if (v == redundant_vid) continue;
            if (!inarray(diamond, nd, v)) {
                if (nd >= 4) {
                    if (dt_trace_on()) fprintf(stderr, "[trace]       f41: star of %d not a diamond (5+ link vertices)\n", redundant_vid);
                    return 0;   /* 5th distinct vertex -> not a diamond */
                }
                diamond[nd++] = v;
            }
        }
    if (nd != 4) {
        if (dt_trace_on()) fprintf(stderr, "[trace]       f41: link has %d vertices (need 4)\n", nd);
        return 0;
    }
    return 1;
}


static int flip41(s_scplx *scplx, s_dstack *stack, s_ncell *ncell, int r_localid, bool *ignored)
{
    int r_vid = ncell->vertex_id[r_localid];
    scplx->N_ncells -= 3;
    ignored[r_vid] = true;

    /* Determine the 4 ncells in the star of redundant_vid */
    s_ncell *star[4]; 
    star[0] = ncell;
    int k = 1;
    for (int i = 0; i < 4; i++) {
        if (i == r_localid) continue;
        s_ncell *nb = ncell->opposite[i];
        if (nb && inarray(nb->vertex_id, 4, r_vid))
            star[k++] = nb;
    }
    assert(k == 4 && "flip41: redundant vertex does not have exactly 4 tetrahedra in its star");

    /* Precompute old opposite pointers before modifying anything */
    s_ncell *nc1_opp_old[4]; opposite_pointers_ncell(star[0], nc1_opp_old);
    s_ncell *nc2_opp_old[4]; opposite_pointers_ncell(star[1], nc2_opp_old);
    s_ncell *nc3_opp_old[4]; opposite_pointers_ncell(star[2], nc3_opp_old);
    s_ncell *nc4_opp_old[4]; opposite_pointers_ncell(star[3], nc4_opp_old);

    /* Local id of redundant vertex in each ncell */
    int nc1_r = id_where_equal_int(star[0]->vertex_id, 4, r_vid);
    int nc2_r = id_where_equal_int(star[1]->vertex_id, 4, r_vid);
    int nc3_r = id_where_equal_int(star[2]->vertex_id, 4, r_vid);
    int nc4_r = id_where_equal_int(star[3]->vertex_id, 4, r_vid);

    /* 4 vertices that will be kept */
    int a = ncell->vertex_id[(r_localid+1)%4];  
    int b = ncell->vertex_id[(r_localid+2)%4];  
    int c = ncell->vertex_id[(r_localid+3)%4];
    int d = -1;
    for (int i=0; i<4; i++) {
        int aux = star[1]->vertex_id[i];
        if (aux != a && aux != b && aux != c && aux != r_vid) d = aux;
    }
    assert(d != -1);

    /* Correct opposite pointer correspondence.  out[j] is the external neighbour
     * across the face opposite abcd[j], i.e. the old outer-neighbour of whichever
     * star tet is MISSING abcd[j].  IMPORTANT: compute this BEFORE overwriting
     * star[0]'s vertices -- star[0] is the tet missing d, so it must still read
     * as {r,a,b,c} here, otherwise its external neighbour (out[d]) is never
     * reconnected and is left as a one-directional (orphan) opposite pointer. */
    int abcd[4] = {a, b, c, d};
    s_ncell *out_raw[4]  = {nc1_opp_old[nc1_r], nc2_opp_old[nc2_r],
                            nc3_opp_old[nc3_r], nc4_opp_old[nc4_r]};
    s_ncell *out[4]      = {NULL, NULL, NULL, NULL};
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            if (!inarray(star[i]->vertex_id, 4, abcd[j]))
                { out[j] = out_raw[i]; break; }

    /* Now overwrite star[0] to be the surviving tetrahedron {a,b,c,d}. */
    set_ncell_vids(star[0], a, b, c, d);

    set_ncell_opposite_pointers(star[0], out[0], out[1], out[2], out[3]);

    /* Update outer neighbors. local ids now correctly match out[] */
    for (int j = 0; j < 4; j++) {
        if (out[j]) {
            int opp_aux; face_localid_of_adjacent_ncell(star[0], 2, &j, j, &opp_aux);
            out[j]->opposite[opp_aux] = star[0];
        }
    }

    /* Remove nc2, nc3, nc4 from linked list */
    if (star[1]->next) star[1]->next->prev = star[1]->prev;
    if (star[1]->prev) star[1]->prev->next = star[1]->next;
    else scplx->head = star[1]->next;

    if (star[2]->next) star[2]->next->prev = star[2]->prev;
    if (star[2]->prev) star[2]->prev->next = star[2]->next;
    else scplx->head = star[2]->next;

    if (star[3]->next) star[3]->next->prev = star[3]->prev;
    if (star[3]->prev) star[3]->prev->next = star[3]->next;
    else scplx->head = star[3]->next;

    if (stack) {
        stack_remove_ncell(stack, star[1]);
        stack_remove_ncell(stack, star[2]);
        stack_remove_ncell(stack, star[3]);
        if (!stack_push(stack, star[0])) return 0;
    }

    free_ncell(scplx, star[1]);
    free_ncell(scplx, star[2]);
    free_ncell(scplx, star[3]);
    if (scplx->point2tet) scplx->point2tet[r_vid] = NULL;
    p2t_update(scplx, star[0]);

    /* Finding E: same invariant hole as flip44 (see the comment in flip44).  The
     * BW loop only re-checks the face opposite the inserted point, but flip41
     * rewired all four faces of the surviving tet star[0]; a face that changed
     * its incident-cell pair may be left non-regular and never re-examined.  Under
     * exact general position ES theory says they stay regular, but the weighted
     * float path has no SoS, so re-check star[0]'s faces explicitly (NONSTRICT to
     * avoid degenerate flip cycles).  flip_tetrahedra may free/replace star[0], so
     * break after the first flip -- the flipped cells are on the stack. */
    if (stack) {
        for (int fi = 0; fi < 4; fi++) {
            if (!star[0]->opposite[fi]) continue;
            if (!are_locally_delaunay(scplx, star[0], fi, DELAUNAY_TEST_NONSTRICT)) {
                if (flip_tetrahedra(scplx, stack, star[0], fi, ignored) == -1) return 0;
                break;
            }
        }
    }
    return 1;
}



// --------------------------------- MAIN ALGORITHM -------------------------------------------

static void regular_tetrahedron(s_point centre, double inradius, s_point out[4])
{   /* Direct 3D construction. */
    double R = 3.0 * inradius;  /* circumradius */

    /* One vertex on top */
    out[0].x = centre.x;
    out[0].y = centre.y;
    out[0].z = centre.z + R;

    /* Three vertices on bottom circle at height centre.z - inradius */
    double r_circle = sqrt(R*R - inradius*inradius);  /* = R*sqrt(8/9) = 2*sqrt(2)*inradius */
    for (int ii = 0; ii < 3; ii++) {
        double angle = 2.0 * M_PI * ii / 3.0;
        out[ii+1].x = centre.x + r_circle * cos(angle);
        out[ii+1].y = centre.y + r_circle * sin(angle);
        out[ii+1].z = centre.z - inradius;
    }
}

static int initialize_scplx(const s_points *points, const double *weights, s_scplx *out,
                             const s_point *bb_min_hint, const s_point *bb_max_hint)
{
    /* scplx->points is extended for the extra nodes of big_ncell, put at the beginning. */
    s_points scplx_points = { .N = points->N + 4,
                              .p = malloc(sizeof(s_point) * (points->N + 4)) };
    if (!scplx_points.p) return 0;
    for (int ii=0; ii<points->N; ii++)
        scplx_points.p[ii+4] = points->p[ii];

    /* Initialize the regular tetrahedron big enough to contain all points. */
    s_point bmin, bmax;
    if (points->N > 0) {
        bounding_box_points(points, &bmin, &bmax);
    } else if (bb_min_hint && bb_max_hint) {
        bmin = *bb_min_hint; bmax = *bb_max_hint;
    } else {
        bmin = (s_point){.x=0,.y=0,.z=0}; bmax = (s_point){.x=0,.y=0,.z=0};
    }
    /* Expand bounding box with optional hint (needed when seeds are few/coincident). */
    if (bb_min_hint) {
        bmin.x = fmin(bmin.x, bb_min_hint->x);
        bmin.y = fmin(bmin.y, bb_min_hint->y);
        bmin.z = fmin(bmin.z, bb_min_hint->z);
    }
    if (bb_max_hint) {
        bmax.x = fmax(bmax.x, bb_max_hint->x);
        bmax.y = fmax(bmax.y, bb_max_hint->y);
        bmax.z = fmax(bmax.z, bb_max_hint->z);
    }
    double bdiag = distance(bmin, bmax);
    s_point centre = scale_point(sum_points(bmax, bmin), 0.5);
    /* Add max radius so big tet contains all balls, not just centers. */
    double maxr = 0.0;
    if (weights) {
        for (int ii = 0; ii < points->N; ii++)
            maxr = fmax(maxr, sqrt(weights[ii]));
    }
    /* inradius min: bdiag/2 + 2*maxr, but the bigger it makes less errors. */
    double inradius = (2*bdiag + 8*maxr);    // TODO TEST WHY TIGHTENING RESULTS IN FAILS
    // double inradius = 1.1 * (bdiag/2 + 2*maxr);  
    /* Deterministic, unperturbed regular tetrahedron.  The sentinel regularity
     * reductions (insphere_sentinel) assume equal-norm sentinel directions
     * about the big-tetra centre, which the exact construction satisfies to
     * ulp; the old rand() jitter broke that premise at 1e-6 relative scale AND
     * made hull-adjacent topology nondeterministic.  Degeneracies between
     * sentinels and data are handled by the exact predicates' tie paths
     * (orthosphere_tiebreak now; SoS in Phase 6f), not by randomness. */
    regular_tetrahedron(centre, inradius, scplx_points.p);

    out->points.N = points->N + 4;
    out->points = scplx_points;
    s_ncell *big_ncell = malloc_ncell(out);
    if (!big_ncell) { free_points(&scplx_points); return 0; }
    for (int ii=0; ii<4; ii++) {
        big_ncell->vertex_id[ii] = ii;
        big_ncell->opposite[ii] = NULL;
    }
    out->head = big_ncell;
    out->N_ncells = 1;
    if (weights) {
        out->weights = malloc(sizeof(double) * scplx_points.N);
        if (!out->weights) {free_ncell(out, big_ncell); free_points(&scplx_points); return 0; }
        for (int ii=0; ii<4; ii++) out->weights[ii] = 0;
        for (int ii=4, jj=0; ii<points->N+4; ii++) out->weights[ii] = weights[jj++];
    } else out->weights = NULL;

    out->point2tet = calloc((size_t)scplx_points.N, sizeof(s_ncell *));
    if (!out->point2tet) { free_ncell(out, big_ncell); free_points(&scplx_points); return 0; }
    for (int k = 0; k < 4; k++)
        out->point2tet[k] = big_ncell;

    return 1;
}


typedef enum type_union_tetra {
    CASE_CONVEX,
    CASE_NON_CONVEX,
    CASE_FLAT,
    CASE_P_IN_EDGE,
    CASE_ERROR
} e_type_union_tetra;

/* Exact id-only mirror of determine_case (CDT_REFACTOR.md Appendix A).  Goes
 * through the dtp_* seam so ids are translated (local cavity DTs) and uses only
 * RELATIVE orient comparisons, so the genericPoint parity flip in dtp_orient is
 * irrelevant (the CASE result is invariant under a global sign flip). */
static e_type_union_tetra determine_case_exact(const s_scplx *s,
                                               int a, int b, int c, int p, int d)
{
    int op = dtp_orient(s, a, b, c, p);
    int tp = dtp_point_in_triangle(s, p, a, b, c);   /* -1 OUT, 0 BND, +1 IN */
    if (tp == 0) return CASE_P_IN_EDGE;

    int o2 = dtp_orient(s, a, b, c, d);
    int itype;  /* 0 EMPTY, 1 DEGENERATE, 2 NONDEGENERATE (of segment pd vs face) */
    if (op == 0 && o2 == 0) {
        /* both endpoints in f's plane: does CLOSED pd meet CLOSED abc? */
        int meets = (dtp_point_in_triangle(s, p, a, b, c) >= 0) ||
                    (dtp_point_in_triangle(s, d, a, b, c) >= 0) ||
                    dtp_segments_cross(s, p, d, a, b) ||
                    dtp_segments_cross(s, p, d, b, c) ||
                    dtp_segments_cross(s, p, d, c, a);
        itype = meets ? 1 : 0;
    } else if (op == 0 || o2 == 0) {
        int te = dtp_point_in_triangle(s, (op == 0) ? p : d, a, b, c);
        itype = (te >= 0) ? 1 : 0;   /* IN or BOUNDARY */
    } else if (op == o2) {
        itype = 0;                   /* same side, no crossing */
    } else {
        int s1 = dtp_orient(s, a, b, p, d);
        int s2 = dtp_orient(s, b, c, p, d);
        int s3 = dtp_orient(s, c, a, p, d);
        if ((s1==0 && s2==s3) || (s2==0 && s1==s3) || (s3==0 && s1==s2) ||
            (s1==0 && s2==0) || (s1==0 && s3==0) || (s2==0 && s3==0))
            itype = 1;
        else if (s1==s2 && s2==s3)
            itype = 2;
        else
            itype = 0;
    }

    switch (itype) {
        case 1: return (tp == 1) ? CASE_CONVEX : CASE_FLAT;
        case 0: return CASE_NON_CONVEX;
        case 2: return CASE_CONVEX;
        default: return CASE_ERROR;
    }
}

static e_type_union_tetra determine_case(const s_scplx *scplx,
                                         int fa, int fb, int fc, int pid, int did)
{   /* two ncells sharing face (fa,fb,fc), with opposite vertices pid and did. */
    if (scplx->exact_ids)
        return determine_case_exact(scplx, fa, fb, fc, pid, did);

    /* Default: the coordinate (already-robust, EPS=TOL=0) implementation. */
    const s_point *P = scplx->points.p;
    s_point vertices_face[3] = { P[fa], P[fb], P[fc] };
    s_point p = P[pid], d = P[did];

    e_geom_test test_p_in_face = test_point_in_triangle_3D(vertices_face, p, 0, 0);
    if (test_p_in_face == TEST_BOUNDARY) return CASE_P_IN_EDGE;

    /* Checks intersection of segment pd with face */
    switch (test_segment_triangle_intersect_3D((s_point[2]){p,d}, vertices_face, 0, 0)) {
        case INTERSECT_DEGENERATE:
            /* either pd has endpoint inside abc, or pd intersects edge/vertex */
            if (test_p_in_face == TEST_IN) return CASE_CONVEX;
            else return CASE_FLAT;
        case INTERSECT_EMPTY: return CASE_NON_CONVEX;
        case INTERSECT_NONDEGENERATE: return CASE_CONVEX;
        case INTERSECT_ERROR:
            fprintf(stderr, "ERROR IN DETERMINE CASE! Degenerate triangle?\n");
            return CASE_ERROR;
        default: return CASE_ERROR;
    }
}

static int flip_tetrahedra(s_scplx *scplx, s_dstack *stack, s_ncell *ncell, int opp_cell_id, bool *ignored)
{   /* -1 ERROR, 0 NOT FLIPPED, 1 FLIPPED */
    if (!ncell->opposite[opp_cell_id]) return 0;

    int face_ids[3];
    extract_ids_face(ncell, 2, &opp_cell_id, face_ids);

    int opp_face_localid;
    face_localid_of_adjacent_ncell(ncell, 2, &opp_cell_id, opp_cell_id, &opp_face_localid);
    int opp_face_vertex_id = (ncell->opposite[opp_cell_id])->vertex_id[opp_face_localid];

    int pid = ncell->vertex_id[opp_cell_id];

    e_type_union_tetra fcase = determine_case(scplx, face_ids[0], face_ids[1], face_ids[2], pid, opp_face_vertex_id);
    if (dt_trace_on())
        fprintf(stderr, "[trace]   flip? tet {%d,%d,%d,%d} face {%d,%d,%d} p=%d opp=%d case=%d\n",
                ncell->vertex_id[0], ncell->vertex_id[1], ncell->vertex_id[2],
                ncell->vertex_id[3], face_ids[0], face_ids[1], face_ids[2],
                pid, opp_face_vertex_id, (int)fcase);
    switch (fcase) {
        int ridge_id_2;
        int redundant_localid;
        case CASE_ERROR:
            fprintf(stderr, "  flip_tetrahedra: tet vids=[%d,%d,%d,%d]  opp_cell_id=%d  opp_face_vid=%d\n",
                    ncell->vertex_id[0], ncell->vertex_id[1],
                    ncell->vertex_id[2], ncell->vertex_id[3],
                    opp_cell_id, opp_face_vertex_id);
            return 0;
        case CASE_CONVEX:
            if (dt_trace_on()) fprintf(stderr, "[trace]     -> flip23\n");
            if (!flip23(scplx, stack, ncell, opp_cell_id, opp_face_localid, NULL)) return -1;
            return 1;
        case CASE_NON_CONVEX:
            if (scplx->weights &&
                can_perform_flip41(ncell, opp_cell_id, &redundant_localid)) {
                    if (dt_trace_on())
                        fprintf(stderr, "[trace]     -> flip41 (remove vid %d)\n",
                                ncell->vertex_id[redundant_localid]);
                    if (!flip41(scplx, stack, ncell, redundant_localid, ignored))
                        return -1;
                    else return 1;
            } else if (can_perform_flip32(ncell, opp_cell_id, &ridge_id_2)) {
                if (dt_trace_on()) fprintf(stderr, "[trace]     -> flip32\n");
                if (!flip32(scplx, stack, ncell, opp_cell_id, ridge_id_2, opp_face_localid, NULL))
                    return -1;
                else return 1;
            } else {
                if (dt_trace_on())
                    fprintf(stderr, "[trace]     -> BLOCKED (non-convex: no flip41/flip32)\n");
                return 0;
            }
        case CASE_FLAT:
            /* Through-vertex sub-case (PLAN_UNSPLIT_EDGE.md): 2 of the 3
             * (p, face-edge, d) orients zero <=> a shared-face vertex r is
             * collinear with, and strictly between, p and d.  In the weighted
             * RT the popped non-regular verdict certifies r locally redundant
             * on the segment: route to the atomic edge-unsplit; flip44's
             * through-vertex gamble is never taken (see NOTE in
             * can_perform_flip44).  Defer while r's star is not yet the double
             * fan.  Unweighted keeps legacy behavior (on-edge vertices are
             * legitimate there and this trigger never pops). */
            if (scplx->weights) {
                int z0 = (dtp_orient(scplx, pid, face_ids[0], face_ids[1], opp_face_vertex_id) == 0);
                int z1 = (dtp_orient(scplx, pid, face_ids[1], face_ids[2], opp_face_vertex_id) == 0);
                int z2 = (dtp_orient(scplx, pid, face_ids[2], face_ids[0], opp_face_vertex_id) == 0);
                if (z0 + z1 + z2 == 2) {
                    int r_vid = (z0 && z1) ? face_ids[1]
                              : (z1 && z2) ? face_ids[2] : face_ids[0];
                    s_unsplit_cfg cfg;
                    if (can_perform_unsplit_edge(scplx, ncell, opp_cell_id, r_vid,
                                                 opp_face_vertex_id, &cfg)) {
                        if (!unsplit_edge(scplx, stack, &cfg, ignored)) return -1;
                        return 1;
                    }
                    if (dt_trace_on())
                        fprintf(stderr, "[trace]     -> DEFERRED (unsplit of %d not applicable yet)\n", r_vid);
                    return 0;
                }
            }
            if (can_perform_flip44(scplx, ncell, opp_cell_id, &ridge_id_2)) {
                if (dt_trace_on()) fprintf(stderr, "[trace]     -> flip44\n");
                if (flip44(scplx, stack, ncell, opp_cell_id, ridge_id_2, NULL, ignored) == -1)
                    return -1;
                else return 1;
            } else {
                if (dt_trace_on())
                    fprintf(stderr, "[trace]     -> BLOCKED (flat: no flip44)\n");
                return 0;
            }
        case CASE_P_IN_EDGE: {
            s_ncell *nbr = ncell->opposite[opp_cell_id];
            /* Ring-of-3 check: if any face of nbr (other than the shared face) borders an
             * existing star tet, then flip23 would create a duplicate tet.  Use flip32 on
             * the ring {ncell, nbr, existing_star} around the collinear edge instead. */
            int ring3_ridge_id2 = -1;
            for (int _fi = 0; _fi < 4; _fi++) {
                if (_fi == opp_face_localid) continue;
                s_ncell *_ext = nbr->opposite[_fi];
                if (_ext && _ext != ncell &&
                    inarray(_ext->vertex_id, 4, ncell->vertex_id[opp_cell_id])) {
                    ring3_ridge_id2 = id_where_equal_int(ncell->vertex_id, 4,
                                                         nbr->vertex_id[_fi]);
                    break;
                }
            }
            if (ring3_ridge_id2 >= 0) {
                if (dt_trace_on()) fprintf(stderr, "[trace]     -> flip32 (p-in-edge ring3)\n");
                if (!flip32(scplx, stack, ncell, opp_cell_id, ring3_ridge_id2, opp_face_localid, NULL))
                    return -1;
                return 1;
            }
            if (dt_trace_on()) fprintf(stderr, "[trace]     -> flip23 (p-in-edge)\n");
            if (!flip23(scplx, stack, ncell, opp_cell_id, opp_face_localid,
                         NULL))
                return -1;
            return 1;
        }
    }
    return 0;
}

static bool point_close_to_ncell_vertex(s_scplx *scplx, s_ncell *ncell, s_point point, double TOL)
{
    const double TOL2 = TOL*TOL;
    if (distance_squared(scplx->points.p[ncell->vertex_id[0]], point) <= TOL2 || 
        distance_squared(scplx->points.p[ncell->vertex_id[1]], point) <= TOL2 ||
        distance_squared(scplx->points.p[ncell->vertex_id[2]], point) <= TOL2 ||
        distance_squared(scplx->points.p[ncell->vertex_id[3]], point) <= TOL2)
        return true;
    else return false;
}

static int insert_one_point(s_scplx *scplx, int point_id, s_dstack *stack, double TOL_dup, bool *ignored)
{   /* -1: ERROR, 0: not inserted, 1: inserted */
    s_point point = scplx->points.p[point_id];
    s_ncell *container_ncell = in_ncell_walk_id(scplx, point_id);

    if (scplx->exact_ids && getenv("CDT_CHECK_WALK")) {
        /* Exact point-in-tet: for each face, the query must be on the same side
         * as the opposite vertex (or on the face) -- robust to tet orientation. */
        const int *V = container_ncell->vertex_id;
        bool inside = true;
        for (int i = 0; i < 4 && inside; i++) {
            int f[3], k = 0;
            for (int j = 0; j < 4; j++) if (j != i) f[k++] = V[j];
            int oo = cdt_orient3d(f[0], f[1], f[2], V[i]);
            int oq = cdt_orient3d(f[0], f[1], f[2], point_id);
            if (oo != 0 && oq != 0 && oo != oq) inside = false;
        }
        if (!inside)
            fprintf(stderr, "WALK BUG: container [%d,%d,%d,%d] does NOT contain point %d (exact)\n",
                    V[0], V[1], V[2], V[3], point_id);
    }

    if (dt_trace_on())
        fprintf(stderr, "[trace] INSERT %d (%.3g,%.3g,%.3g) container {%d,%d,%d,%d}\n",
                point_id, point.x, point.y, point.z,
                container_ncell->vertex_id[0], container_ncell->vertex_id[1],
                container_ncell->vertex_id[2], container_ncell->vertex_id[3]);

    if (point_close_to_ncell_vertex(scplx, container_ncell, point, TOL_dup) ||
        p_locally_redundant_in_ncell(scplx, container_ncell, point_id)) {
        if (dt_trace_on()) fprintf(stderr, "[trace]   -> dropped (dup/redundant)\n");
        ignored[point_id] = true;
        return 0;
    }

    if (!flip14(scplx, container_ncell, point_id, stack)) return -1;

    while (stack->size > 0) {
        s_ncell *current = stack_pop(stack);
        if (!inarray(current->vertex_id, 4, point_id)) continue;
        int opp_cell_id = id_where_equal_int(current->vertex_id, 4, point_id);
        if (!current->opposite[opp_cell_id]) continue;

        /* Zero-volume tets (point on a DT edge) have insphere=0 so NONSTRICT skips them.
         * Force a flip to let the collinear-edge cascade complete. */
        {
            bool on_boundary;
            if (scplx->exact_ids) {
                int f_ids[3]; extract_ids_face(current, 2, &opp_cell_id, f_ids);
                on_boundary = (dtp_point_in_triangle(scplx, point_id, f_ids[0], f_ids[1], f_ids[2]) == 0);
            } else {
                s_point face_pts[3];
                extract_vertices_face(scplx, current, 2, &opp_cell_id, face_pts);
                on_boundary = (test_point_in_triangle_3D(face_pts, point, 0, 0) == TEST_BOUNDARY);
            }
            if (on_boundary) {
                if (dt_trace_on())
                    fprintf(stderr, "[trace]   pop {%d,%d,%d,%d}: p ON opposite face -> forced flip\n",
                            current->vertex_id[0], current->vertex_id[1],
                            current->vertex_id[2], current->vertex_id[3]);
                if (flip_tetrahedra(scplx, stack, current, opp_cell_id, ignored) == -1) return -1;
                continue;
            }
        }
        if (!are_locally_delaunay(scplx, current, opp_cell_id, DELAUNAY_TEST_NONSTRICT)) {
            if (dt_trace_on())
                fprintf(stderr, "[trace]   pop {%d,%d,%d,%d}: NOT locally regular -> flip\n",
                        current->vertex_id[0], current->vertex_id[1],
                        current->vertex_id[2], current->vertex_id[3]);
            if (flip_tetrahedra(scplx, stack, current, opp_cell_id, ignored) == -1) return -1;
        }
    }

    return 1;
}

int scplx_insert_point(s_scplx *dt, s_point p, double TOL)
{
    int new_id = dt->points.N;
    s_point *tmp = realloc(dt->points.p, (size_t)(new_id + 1) * sizeof(s_point));
    if (!tmp) return -1;
    dt->points.p = tmp;
    dt->points.p[new_id] = p;
    dt->points.N = new_id + 1;

    bool *ignored = calloc((size_t)(new_id + 1), sizeof(bool));
    if (!ignored) { dt->points.N--; return -1; }

    s_dstack stack = stack_create();
    if (!stack.entry) { free(ignored); dt->points.N--; return -1; }

    int res = insert_one_point(dt, new_id, &stack, TOL, ignored);
    stack_free(&stack);
    free(ignored);

    if (res <= 0) { dt->points.N--; return -1; }
    return new_id;
}

/*
 * kept_idx, if non-NULL, is composed in place: for each i in [0, N_kept_idx),
 * if kept_idx[i] is a valid index into the pre-compaction "real seed" range
 * (i.e. was produced by an earlier filtering stage as an index there), it is
 * rewritten to that seed's final compacted index, or -1 if the seed was
 * dropped here as a near-duplicate. Reuses the remap table already built
 * for compaction below -- no extra allocation needed.
 */
static void remove_ignored_points(s_scplx *scplx, bool *ignored, bool keep_big_tetra,
                                  int *kept_idx, int N_kept_idx)
{
    if (!keep_big_tetra) {
        /* Mark first 4 indices as ignored */
        for (int i = 0; i < 4; i++) ignored[i] = true;
        
        /* Remove any ncell referencing big tetra vertices */
        s_ncell *current = scplx->head;
        while (current) {
            s_ncell *next = current->next;
            for (int ii=0; ii<4; ii++) if (current->vertex_id[ii] < 4) { 
                if (current->next) (current->next)->prev = current->prev;
                if (current->prev) (current->prev)->next = next;
                else scplx->head = current->next;

                /* Update opposite's opposite to NULL */
                for (int jj=0; jj<4; jj++) if (current->opposite[jj]) {
                    for (int kk=0; kk<4; kk++) if (current->opposite[jj]->opposite[kk] == current) {
                        current->opposite[jj]->opposite[kk] = NULL;
                        break;
                    }
                }

                free_ncell(scplx, current);
                scplx->N_ncells--;
                break;
            }
            current = next;
        }
    }

    /* Compact points / weights, and build remap table */
    int *remap = malloc(sizeof(int) * scplx->points.N);
    int k = 0;
    for (int i = 0; i < scplx->points.N; i++) {
        if (ignored[i]) { remap[i] = -1; continue; }
        scplx->points.p[k] = scplx->points.p[i];
        if (scplx->weights) scplx->weights[k] = scplx->weights[i];
        remap[i] = k++;
    }
    scplx->points.N = k;
    scplx->points.p = realloc(scplx->points.p, sizeof(s_point) * k);
    if (scplx->weights) scplx->weights = realloc(scplx->weights, sizeof(double) * k);
    if (scplx->point2tet) {
        s_ncell **tmp = realloc(scplx->point2tet, (size_t)k * sizeof(s_ncell *));
        if (tmp) scplx->point2tet = tmp;
        memset(scplx->point2tet, 0, (size_t)k * sizeof(s_ncell *));
    }

    /* Update vertex ids and rebuild point2tet with compacted indices */
    for (s_ncell *c = scplx->head; c; c = c->next) {
        for (int i = 0; i < 4; i++) {
            assert(remap[c->vertex_id[i]] != -1 && "Ignored point still referenced by tetra.");
            c->vertex_id[i] = remap[c->vertex_id[i]];
        }
        if (scplx->point2tet)
            for (int i = 0; i < 4; i++)
                scplx->point2tet[c->vertex_id[i]] = c;
    }

    if (kept_idx)
        for (int i = 0; i < N_kept_idx; i++)
            if (kept_idx[i] >= 0)
                kept_idx[i] = remap[4 + kept_idx[i]];

    free(remap);
}




/* BUILDER */
static s_dt_builder dt_builder_begin_impl(const s_points *seeds, const double *weights,
                              double TOL_dup, const s_point *bb_min_hint,
                              const s_point *bb_max_hint, bool exact,
                              const int *l2g_real, int scratch_base,
                              s_random_context *rng)
{
    s_dt_builder b = {0};
    assert(!(exact && weights) && "exact-id mode does not support weights");

    bool *ignored = calloc(seeds->N + 4, sizeof(bool));
    if (!ignored) return b;

    s_dstack *stack = malloc(sizeof(s_dstack));
    if (!stack) { free(ignored); return b; }
    *stack = stack_create();
    if (!stack->entry) { free(ignored); free(stack); return b; }

    if (!initialize_scplx(seeds, weights, &b.dt, bb_min_hint, bb_max_hint)) {
        free(ignored); stack_free(stack); free(stack); return b;
    }
    /* Attach the caller's PRNG BEFORE the seed-insertion loop below so every
     * point-location walk during construction breaks ties deterministically. */
    b.dt.rng = rng;

    if (exact && !l2g_real) {
        /* Exact GLOBAL DT: OWN the registry -- register every point (sentinels
         * 0..3 + seeds) as an explicit entry.  Registry index == scplx point
         * index.  initialize_scplx already filled points.p[] completely. */
        b.dt.exact_ids = 1;
        cdt_predicates_clear();
        for (int i = 0; i < b.dt.points.N; i++) {
            s_point *p = &b.dt.points.p[i];
            cdt_point_set_explicit(i, p->x, p->y, p->z);
        }
    } else if (exact) {
        /* Exact LOCAL cavity DT: SHARE the global registry (do not clear it).
         * Vertex ids stay local; the seam translates via l2g_ids.  Real seeds
         * map to their existing global ids; the 4 local sentinels are registered
         * into scratch registry slots [scratch_base..scratch_base+3]. */
        b.dt.exact_ids = 1;
        int local_N = b.dt.points.N;            /* 4 + n */
        int *l2g = malloc((size_t)local_N * sizeof(int));
        if (!l2g) { free(ignored); stack_free(stack); free(stack);
                    free_complex(&b.dt); return (s_dt_builder){0}; }
        for (int j = 0; j < 4; j++)        l2g[j] = scratch_base + j;
        for (int i = 4; i < local_N; i++)  l2g[i] = l2g_real[i - 4];
        b.dt.l2g_ids = l2g;
        b._l2g_ids   = l2g;
        for (int j = 0; j < 4; j++)
            cdt_point_set_explicit(scratch_base + j, b.dt.points.p[j].x,
                                   b.dt.points.p[j].y, b.dt.points.p[j].z);
    }

    for (int ii = 4; ii < b.dt.points.N; ii++) {
        if (insert_one_point(&b.dt, ii, stack, TOL_dup, ignored) == -1) {
            free(ignored); stack_free(stack); free(stack);
            free(b._l2g_ids); free_complex(&b.dt);
            return (s_dt_builder){0};
        }
        if (dt_trace_on()) {
            char tag[64]; snprintf(tag, sizeof tag, "after insert %d", ii);
            trace_dump_tets(&b.dt, tag);
        }
    }

    b._ignored = ignored;
    b._stack = stack;
    b._N_original = seeds->N;

    /* Accumulated bbox: same union as initialize_scplx used for the big tetra. */
    if (seeds->N > 0) {
        bounding_box_points(seeds, &b._bbox_min, &b._bbox_max);
    } else {
        b._bbox_min = b._bbox_max = (s_point){.x=0,.y=0,.z=0};
    }
    if (bb_min_hint) {
        b._bbox_min.x = fmin(b._bbox_min.x, bb_min_hint->x);
        b._bbox_min.y = fmin(b._bbox_min.y, bb_min_hint->y);
        b._bbox_min.z = fmin(b._bbox_min.z, bb_min_hint->z);
    }
    if (bb_max_hint) {
        b._bbox_max.x = fmax(b._bbox_max.x, bb_max_hint->x);
        b._bbox_max.y = fmax(b._bbox_max.y, bb_max_hint->y);
        b._bbox_max.z = fmax(b._bbox_max.z, bb_max_hint->z);
    }

    /* Sentinel centroid and inradius from the built vertices.
     * The big tet is nearly regular (perturbation 1e-6), so inradius ~ circumradius/3. */
    b._sentinel_center = (s_point){.x=0,.y=0,.z=0};
    for (int i = 0; i < 4; i++) {
        b._sentinel_center.x += b.dt.points.p[i].x * 0.25;
        b._sentinel_center.y += b.dt.points.p[i].y * 0.25;
        b._sentinel_center.z += b.dt.points.p[i].z * 0.25;
    }
    double max_R = 0;
    for (int i = 0; i < 4; i++) {
        double dx = b.dt.points.p[i].x - b._sentinel_center.x;
        double dy = b.dt.points.p[i].y - b._sentinel_center.y;
        double dz = b.dt.points.p[i].z - b._sentinel_center.z;
        double R = sqrt(dx*dx + dy*dy + dz*dz);
        if (R > max_R) max_R = R;
    }
    b._sentinel_inradius = max_R / 3.0;

    return b;
}

s_dt_builder dt_builder_begin(const s_points *seeds, const double *weights, double TOL_dup,
                              const s_point *bb_min_hint, const s_point *bb_max_hint,
                              s_random_context *rng)
{
    return dt_builder_begin_impl(seeds, weights, TOL_dup, bb_min_hint, bb_max_hint,
                                 false, NULL, 0, rng);
}

/* Exact-id build: predicates run on the cdt_predicates registry (no weights).
 * The registry must be init'd (cdt_predicates_init) beforehand; this clears and
 * repopulates it for the seeds.  Steiners are then added via
 * dt_builder_extend_lnc, real points via dt_builder_extend. */
s_dt_builder dt_builder_begin_exact(const s_points *seeds, double TOL_dup,
                                    const s_point *bb_min_hint, const s_point *bb_max_hint,
                                    s_random_context *rng)
{
    return dt_builder_begin_impl(seeds, NULL, TOL_dup, bb_min_hint, bb_max_hint,
                                 true, NULL, 0, rng);
}

s_dt_builder dt_builder_begin_exact_local(const s_points *seeds, double TOL_dup,
                                          const int *l2g_real, int n, int scratch_base,
                                          const s_point *bb_min_hint,
                                          const s_point *bb_max_hint,
                                          s_random_context *rng)
{
    (void)n;  /* n == seeds->N; the map is sized from points.N inside the impl */
    return dt_builder_begin_impl(seeds, NULL, TOL_dup, bb_min_hint, bb_max_hint,
                                 true, l2g_real, scratch_base, rng);
}


bool dt_builder_extend(s_dt_builder *b, const s_points *new_points, double TOL_dup)
{
    /* Expand accumulated bbox to include new points, then grow the big tetra if needed
     * so every new point falls within its inscribed sphere (same 2xbdiag formula as
     * initialize_scplx, applied to the now-larger bounding box). */
    for (int i = 0; i < new_points->N; i++) {
        s_point p = new_points->p[i];
        b->_bbox_min.x = fmin(b->_bbox_min.x, p.x);
        b->_bbox_min.y = fmin(b->_bbox_min.y, p.y);
        b->_bbox_min.z = fmin(b->_bbox_min.z, p.z);
        b->_bbox_max.x = fmax(b->_bbox_max.x, p.x);
        b->_bbox_max.y = fmax(b->_bbox_max.y, p.y);
        b->_bbox_max.z = fmax(b->_bbox_max.z, p.z);
    }
    double needed = 2.0 * distance(b->_bbox_min, b->_bbox_max);
    if (needed > b->_sentinel_inradius) {
        double factor = needed / b->_sentinel_inradius;
        for (int i = 0; i < 4; i++) {
            b->dt.points.p[i].x = b->_sentinel_center.x + factor * (b->dt.points.p[i].x - b->_sentinel_center.x);
            b->dt.points.p[i].y = b->_sentinel_center.y + factor * (b->dt.points.p[i].y - b->_sentinel_center.y);
            b->dt.points.p[i].z = b->_sentinel_center.z + factor * (b->dt.points.p[i].z - b->_sentinel_center.z);
        }
        b->_sentinel_inradius = needed;
        /* Exact mode: sentinels moved -> refresh their registry coords. */
        if (b->dt.exact_ids)
            for (int i = 0; i < 4; i++)
                cdt_point_set_explicit(i, b->dt.points.p[i].x, b->dt.points.p[i].y, b->dt.points.p[i].z);
    }

    int N_old = b->dt.points.N;
    int N_add = new_points->N;
    int N_new = N_old + N_add;

    s_point *tmp_p = realloc(b->dt.points.p, N_new * sizeof(s_point));
    if (!tmp_p) return false;
    b->dt.points.p = tmp_p;
    memcpy(&b->dt.points.p[N_old], new_points->p, N_add * sizeof(s_point));
    b->dt.points.N = N_new;

    bool *tmp_ig = realloc(b->_ignored, N_new * sizeof(bool));
    if (!tmp_ig) return false;
    b->_ignored = tmp_ig;
    memset(&b->_ignored[N_old], 0, N_add * sizeof(bool));

    if (b->dt.weights) {
        double *tmp_w = realloc(b->dt.weights, N_new * sizeof(double));
        if (!tmp_w) return false;
        b->dt.weights = tmp_w;
        for (int i = N_old; i < N_new; i++) b->dt.weights[i] = 0.0;
    }

    if (b->dt.point2tet) {
        s_ncell **tmp_p2t = realloc(b->dt.point2tet, N_new * sizeof(s_ncell *));
        if (!tmp_p2t) return false;
        b->dt.point2tet = tmp_p2t;
        for (int i = N_old; i < N_new; i++) b->dt.point2tet[i] = NULL;
    }

    /* Exact mode: new explicit points must be in the registry before insertion. */
    if (b->dt.exact_ids)
        for (int i = N_old; i < N_new; i++)
            cdt_point_set_explicit(i, b->dt.points.p[i].x, b->dt.points.p[i].y, b->dt.points.p[i].z);

    s_dstack *stack = (s_dstack *)b->_stack;
    for (int i = N_old; i < N_new; i++)
        if (insert_one_point(&b->dt, i, stack, TOL_dup, b->_ignored) == -1)
            return false;

    return true;
}

/* Insert one implicit Steiner point s*V[v1] + (1-s)*V[v2] (cdt_point_set_lnc
 * convention) into an exact-id builder.  Registers the exact LNC BEFORE
 * insertion (so the insertion predicates see the true implicit position) and
 * stores the rounded double in points.p (the output representation, consumed by
 * float clients).  v1, v2 must be registered explicit points.  Returns like
 * dt_builder_extend. */
bool dt_builder_extend_lnc(s_dt_builder *b, int v1, int v2, double s, double TOL_dup)
{
    assert(b->dt.exact_ids && "dt_builder_extend_lnc requires an exact-id builder");

    s_point A = b->dt.points.p[v1];
    s_point B = b->dt.points.p[v2];
    s_point M = { .x = s*A.x + (1.0-s)*B.x,
                  .y = s*A.y + (1.0-s)*B.y,
                  .z = s*A.z + (1.0-s)*B.z };

    /* Steiners are on a segment between existing points, so essentially
     * interior -- but the rounded M can land a hair outside the accumulated
     * bbox, so grow the big tetra if the safety margin is exceeded (same rule
     * as dt_builder_extend).  extend_lnc is always exact, so refresh the moved
     * sentinels in the registry. */
    b->_bbox_min.x = fmin(b->_bbox_min.x, M.x); b->_bbox_max.x = fmax(b->_bbox_max.x, M.x);
    b->_bbox_min.y = fmin(b->_bbox_min.y, M.y); b->_bbox_max.y = fmax(b->_bbox_max.y, M.y);
    b->_bbox_min.z = fmin(b->_bbox_min.z, M.z); b->_bbox_max.z = fmax(b->_bbox_max.z, M.z);
    double needed = 2.0 * distance(b->_bbox_min, b->_bbox_max);
    if (needed > b->_sentinel_inradius) {
        double factor = needed / b->_sentinel_inradius;
        for (int i = 0; i < 4; i++) {
            b->dt.points.p[i].x = b->_sentinel_center.x + factor * (b->dt.points.p[i].x - b->_sentinel_center.x);
            b->dt.points.p[i].y = b->_sentinel_center.y + factor * (b->dt.points.p[i].y - b->_sentinel_center.y);
            b->dt.points.p[i].z = b->_sentinel_center.z + factor * (b->dt.points.p[i].z - b->_sentinel_center.z);
        }
        b->_sentinel_inradius = needed;
        for (int i = 0; i < 4; i++)
            cdt_point_set_explicit(i, b->dt.points.p[i].x, b->dt.points.p[i].y, b->dt.points.p[i].z);
    }

    int new_id = b->dt.points.N;
    int N_new = new_id + 1;

    s_point *tmp_p = realloc(b->dt.points.p, N_new * sizeof(s_point));
    if (!tmp_p) return false;
    b->dt.points.p = tmp_p;
    b->dt.points.p[new_id] = M;
    b->dt.points.N = N_new;

    bool *tmp_ig = realloc(b->_ignored, N_new * sizeof(bool));
    if (!tmp_ig) return false;
    b->_ignored = tmp_ig;
    b->_ignored[new_id] = false;

    if (b->dt.point2tet) {
        s_ncell **tmp_p2t = realloc(b->dt.point2tet, N_new * sizeof(s_ncell *));
        if (!tmp_p2t) return false;
        b->dt.point2tet = tmp_p2t;
        b->dt.point2tet[new_id] = NULL;
    }

    cdt_point_set_lnc(new_id, v1, v2, s);   /* before insertion */

    s_dstack *stack = (s_dstack *)b->_stack;
    if (insert_one_point(&b->dt, new_id, stack, TOL_dup, b->_ignored) == -1)
        return false;

    return true;
}


s_scplx dt_builder_end(s_dt_builder *b, bool keep_big_tetra, int *out_Nreal,
                       int *kept_idx, int N_kept_idx)
{
    // if (!is_delaunay_3d(&b->dt, DELAUNAY_TEST_NONSTRICT)) {
    //     fprintf(stderr, "WARNING: DT is not Delaunay!\n");
    //     write_points_to_csv("error.csv", "w", &b->dt.points);
    //     // s_ncell *c = b->dt.head;
    //     // while (c) {
    //     //     s_point v[4]; extract_vertices_ncell(&b->dt, c, v);
    //     //     if (test_orientation(v, v[3]) == 0) {
    //     //         fprintf(stderr, "Flat tetra coords:\n");
    //     //         fprintf(stderr, "  v0(%d): %.6f %.6f %.6f\n", c->vertex_id[0], v[0].x, v[0].y, v[0].z);
    //     //         fprintf(stderr, "  v1(%d): %.6f %.6f %.6f\n", c->vertex_id[1], v[1].x, v[1].y, v[1].z);
    //     //         fprintf(stderr, "  v2(%d): %.6f %.6f %.6f\n", c->vertex_id[2], v[2].x, v[2].y, v[2].z);
    //     //         fprintf(stderr, "  v3(%d): %.6f %.6f %.6f\n", c->vertex_id[3], v[3].x, v[3].y, v[3].z);
    //     //     }
    //     //     c = c->next;
    //     // }
    //     exit(1);
    // }

    /* No truly-flat (zero-volume) tet is ever allowed in the result -- in exact
     * mode SoS guarantees it by construction; in float mode a flat tet means a
     * genuine bug (degenerate input that slipped through).  Hard-fail either
     * way.  Ghost/big-tetra tets (a vertex < 4) are skipped: they are the
     * bounding sentinels, discarded downstream. */
    for (s_ncell *c = b->dt.head; c; c = c->next) {
        const int *v = c->vertex_id;
        if (v[0] < 4 || v[1] < 4 || v[2] < 4 || v[3] < 4) continue;
        if (dtp_orient(&b->dt, v[0], v[1], v[2], v[3]) == 0) {
            fprintf(stderr, "FATAL: zero-volume (flat) tetra in DT: %d %d %d %d "
                    "(exact_ids=%d) -- flat tets are never allowed.\n",
                    v[0], v[1], v[2], v[3], b->dt.exact_ids);
            abort();
        }
    }

    if (out_Nreal != NULL) {
        int surviving = 0;
        for (int i = 4; i < 4 + b->_N_original; i++)
            if (!b->_ignored[i]) surviving++;
        *out_Nreal = surviving;
    }

    remove_ignored_points(&b->dt, b->_ignored, keep_big_tetra, kept_idx, N_kept_idx);

    s_dstack *stack = (s_dstack *)b->_stack;
    stack_free(stack);
    free(stack);
    free(b->_ignored);
    b->_stack = NULL;
    b->_ignored = NULL;

    /* Exact-local map is consumed by the build + flat scan above; the returned
     * complex is used only topologically (face_in_dt, commit), so drop it. */
    free(b->_l2g_ids);
    b->_l2g_ids = NULL;
    b->dt.l2g_ids = NULL;

    // assert_point2tet(&b->dt);
    return b->dt;
}


s_scplx construct_dt_3d(const s_points *points, const double *weights,
                        bool keep_big_tetra, double TOL_duplicates, int *out_Nreal,
                        s_random_context *rng)
{
    s_dt_builder b = dt_builder_begin(points, weights, TOL_duplicates, NULL, NULL, rng);
    if (!b._stack) {
        fprintf(stderr, "construct_dt_3d: Error.\n");
        return (s_scplx){0};
    }
    /* Preserve the old behaviour: *out_Nreal on entry specifies how many of the
     * first points->N entries are "real" seeds (the rest are mirrors added by
     * extend_sites_mirroring).  Store that count so dt_builder_end uses it. */
    if (out_Nreal != NULL && *out_Nreal <= points->N)
        b._N_original = *out_Nreal;
    return dt_builder_end(&b, keep_big_tetra, out_Nreal, NULL, 0);
}




int N_ncells_not_big_tetra(s_scplx *scplx)
{
    int N_remove = 0;

    s_ncell *current = scplx->head;
    while (current) {
        for (int ii=0; ii<4; ii++) if (current->vertex_id[ii] < 4) {
            N_remove++;
            break;
        }
        current = current->next;
    }
    return scplx->N_ncells - N_remove;
}

