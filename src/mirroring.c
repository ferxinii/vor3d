
#include "points.h"
#include "dynarray.h"
#include "voronoi_predicates.h"
#include "robust_predicates.h"
#include "vdiagram.h"
#include "scplx.h"
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <stdio.h>

/* Identifies a constraint.  TRI* idx is unused (set to -1).
 * CSEED idx is the constraint's position in the best_idx[] array. */
typedef enum { TRI0, TRI1, TRI2, CSEED } e_ctype;
typedef struct { e_ctype type; int idx; } s_cid;

static int lp_D(const s_point face[3], s_point s,
         const s_points *seeds, const int *best_idx,
         s_cid i, s_cid j)
{
    s_point A = face[0], B = face[1], C = face[2];

    /* CSEED x CSEED */
    if (i.type == CSEED && j.type == CSEED) {
        s_point ti = seeds->p[best_idx[i.idx]];
        s_point tj = seeds->p[best_idx[j.idx]];
        return lp_det2(A.x,A.y,A.z, B.x,B.y,B.z, C.x,C.y,C.z,
                       s.x,s.y,s.z,
                       ti.x,ti.y,ti.z, tj.x,tj.y,tj.z);
    }

    /* CSEED x TRI */
    if (i.type == CSEED && j.type != CSEED) {
        s_point ti = seeds->p[best_idx[i.idx]];
        if (j.type == TRI0) return -lp_D_T0_S(A.x,A.y,A.z, C.x,C.y,C.z, s.x,s.y,s.z, ti.x,ti.y,ti.z);
        if (j.type == TRI1) return -lp_D_T1_S(A.x,A.y,A.z, B.x,B.y,B.z, s.x,s.y,s.z, ti.x,ti.y,ti.z);
        if (j.type == TRI2) return -lp_D_T2_S(B.x,B.y,B.z, C.x,C.y,C.z, s.x,s.y,s.z, ti.x,ti.y,ti.z);
    }

    /* TRI x CSEED */
    if (i.type != CSEED && j.type == CSEED) {
        s_point tj = seeds->p[best_idx[j.idx]];
        if (i.type == TRI0) return lp_D_T0_S(A.x,A.y,A.z, C.x,C.y,C.z, s.x,s.y,s.z, tj.x,tj.y,tj.z);
        if (i.type == TRI1) return lp_D_T1_S(A.x,A.y,A.z, B.x,B.y,B.z, s.x,s.y,s.z, tj.x,tj.y,tj.z);
        if (i.type == TRI2) return lp_D_T2_S(B.x,B.y,B.z, C.x,C.y,C.z, s.x,s.y,s.z, tj.x,tj.y,tj.z);
    }

    /* TRI x TRI */
    if (i.type == j.type) return 0;  // Repeated arguments == 0
    else if (i.type == TRI0 && j.type == TRI1) return 1;
    else if (i.type == TRI1 && j.type == TRI0) return -1;
    else if (i.type == TRI1 && j.type == TRI2) return 1;
    else if (i.type == TRI2 && j.type == TRI1) return -1;
    else if (i.type == TRI2 && j.type == TRI0) return 1;
    else if (i.type == TRI0 && j.type == TRI2) return -1;

    assert(0 && "lp_D: case not implemented");
    return 0;
}


static int lp_feasible(const s_point face[3], s_point s,
                const s_points *seeds, const int *best_idx,
                s_cid L, s_cid M, s_cid N)
{
    s_point A = face[0], B = face[1], C = face[2];

    /* duplicate CSEED constraints -> degenerate intersection, return 0 */
    if (L.type == CSEED && M.type == CSEED && L.idx == M.idx) return 0;
    if (L.type == CSEED && N.type == CSEED && L.idx == N.idx) return 0;
    if (M.type == CSEED && N.type == CSEED && M.idx == N.idx) return 0;

    /* swap so TRI always precedes CSEED in (L,M); sign is unchanged by this swap */
    if (L.type == CSEED && M.type != CSEED) { s_cid tmp = L; L = M; M = tmp; }

    /* TRI x TRI x CSEED */
    if (L.type != CSEED && M.type != CSEED && N.type == CSEED) {
        s_point ti = seeds->p[best_idx[N.idx]];
        if (L.type == TRI0 && M.type == TRI1)
            return lp_feasible_T0_T1_S(A.x,A.y,A.z, s.x,s.y,s.z, ti.x,ti.y,ti.z);
        if (L.type == TRI1 && M.type == TRI2)
            return lp_feasible_T1_T2_S(B.x,B.y,B.z, s.x,s.y,s.z, ti.x,ti.y,ti.z);
        if (L.type == TRI0 && M.type == TRI2)
            return lp_feasible_T0_T2_S(C.x,C.y,C.z, s.x,s.y,s.z, ti.x,ti.y,ti.z);
    }

    /* TRI x CSEED x CSEED */
    if (L.type != CSEED && M.type == CSEED && N.type == CSEED) {
        s_point tj = seeds->p[best_idx[M.idx]];
        s_point tl = seeds->p[best_idx[N.idx]];
        if (L.type == TRI0)
            return lp_feasible_T0_S_S(A.x,A.y,A.z, B.x,B.y,B.z, C.x,C.y,C.z,
                                      s.x,s.y,s.z, tj.x,tj.y,tj.z, tl.x,tl.y,tl.z);
        if (L.type == TRI1)
            return lp_feasible_T1_S_S(A.x,A.y,A.z, B.x,B.y,B.z, C.x,C.y,C.z,
                                      s.x,s.y,s.z, tj.x,tj.y,tj.z, tl.x,tl.y,tl.z);
        if (L.type == TRI2)
            return lp_feasible_T2_S_S(A.x,A.y,A.z, B.x,B.y,B.z, C.x,C.y,C.z,
                                      s.x,s.y,s.z, tj.x,tj.y,tj.z, tl.x,tl.y,tl.z);
    }

    /* TRI x CSEED x TRI */
    if (L.type != CSEED && M.type == CSEED && N.type != CSEED) {
        if (N.type == L.type) return 0;
        s_point tj = seeds->p[best_idx[M.idx]];
        if (L.type == TRI0 && N.type == TRI1)
            return lp_feasible_T0_S_T1(A.x,A.y,A.z, C.x,C.y,C.z, s.x,s.y,s.z, tj.x,tj.y,tj.z);
        if (L.type == TRI0 && N.type == TRI2)
            return lp_feasible_T0_S_T2(A.x,A.y,A.z, C.x,C.y,C.z, s.x,s.y,s.z, tj.x,tj.y,tj.z);
        if (L.type == TRI1 && N.type == TRI0)
            return lp_feasible_T1_S_T0(A.x,A.y,A.z, B.x,B.y,B.z, s.x,s.y,s.z, tj.x,tj.y,tj.z);
        if (L.type == TRI1 && N.type == TRI2)
            return lp_feasible_T1_S_T2(A.x,A.y,A.z, B.x,B.y,B.z, s.x,s.y,s.z, tj.x,tj.y,tj.z);
        if (L.type == TRI2 && N.type == TRI0)
            return lp_feasible_T2_S_T0(B.x,B.y,B.z, C.x,C.y,C.z, s.x,s.y,s.z, tj.x,tj.y,tj.z);
        if (L.type == TRI2 && N.type == TRI1)
            return lp_feasible_T2_S_T1(B.x,B.y,B.z, C.x,C.y,C.z, s.x,s.y,s.z, tj.x,tj.y,tj.z);
    }

    /* CSEED x CSEED x TRI */
    if (L.type == CSEED && M.type == CSEED && N.type != CSEED) {
        s_point tj = seeds->p[best_idx[L.idx]];
        s_point tk = seeds->p[best_idx[M.idx]];
        if (N.type == TRI0)
            return lp_feasible_S_S_T0(A.x,A.y,A.z, B.x,B.y,B.z, C.x,C.y,C.z,
                                      s.x,s.y,s.z, tj.x,tj.y,tj.z, tk.x,tk.y,tk.z);
        if (N.type == TRI1)
            return lp_feasible_S_S_T1(A.x,A.y,A.z, B.x,B.y,B.z, C.x,C.y,C.z,
                                      s.x,s.y,s.z, tj.x,tj.y,tj.z, tk.x,tk.y,tk.z);
        if (N.type == TRI2)
            return lp_feasible_S_S_T2(A.x,A.y,A.z, B.x,B.y,B.z, C.x,C.y,C.z,
                                      s.x,s.y,s.z, tj.x,tj.y,tj.z, tk.x,tk.y,tk.z);
    }

    /* TRI x TRI x TRI: intersection is a vertex, always inside unless any two match */
    if (L.type != CSEED && M.type != CSEED && N.type != CSEED)
        return (L.type == M.type || N.type == L.type || N.type == M.type) ? 0 : 1;

    /* CSEED x CSEED x CSEED */
    if (L.type == CSEED && M.type == CSEED && N.type == CSEED) {
        s_point tL = seeds->p[best_idx[L.idx]];
        s_point tM = seeds->p[best_idx[M.idx]];
        s_point tN = seeds->p[best_idx[N.idx]];
        return lp_det2(A.x,A.y,A.z, B.x,B.y,B.z, C.x,C.y,C.z,
                       s.x,s.y,s.z, tL.x,tL.y,tL.z, tM.x,tM.y,tM.z)
             * lp_det3(A.x,A.y,A.z, B.x,B.y,B.z, C.x,C.y,C.z,
                       s.x,s.y,s.z, tL.x,tL.y,tL.z, tM.x,tM.y,tM.z,
                       tN.x,tN.y,tN.z);
    }

    fprintf(stderr, "lp_feasible: unhandled constraint combination\n");
    return -1;
}


/* Update one 1D-LP bound. Returns false if the interval becomes empty or
 * the line is globally infeasible. When concurrent=true, a zero feasibility
 * result while tightening a non-null slot signals a concurrent triangle. */
static bool lp_process_cj(const s_point face[3], s_point s,
                           const s_points *seeds, const int *seeds_idx,
                           s_cid ci, s_cid cj,
                           s_cid *lo, bool *lo_null,
                           s_cid *hi, bool *hi_null)
{
    int d = lp_D(face, s, seeds, seeds_idx, ci, cj);
    if (d == 0) {
        /* C_j parallel to C_i: use whichever bound is non-null as reference
         * (same sign for any point on C_i when C_j is parallel). */
        s_cid ref = *lo_null ? *hi : *lo;
        int v = lp_feasible(face, s, seeds, seeds_idx, ci, ref, cj);
        if (v < 0) return false;
    } else if (d > 0) {
        if (*lo_null) {
            *lo = cj; *lo_null = false;
        } else {
            int v = lp_feasible(face, s, seeds, seeds_idx, ci, cj, *lo);
            if (v > 0) *lo = cj;
        }
    } else {
        if (*hi_null) {
            *hi = cj; *hi_null = false;
        } else {
            int v = lp_feasible(face, s, seeds, seeds_idx, ci, cj, *hi);
            if (v > 0) *hi = cj;
        }
    }
    if (!*lo_null && !*hi_null) {
        int empty = lp_feasible(face, s, seeds, seeds_idx, ci, *lo, *hi);
        if (empty < 0) return false;
    }
    return true;
}


static bool lp_should_mirror(const s_point face[3], const s_points *seeds,
                        int N, int seeds_idx[N], int id,
                        int (*randint)(void *rctx, int), void *rctx)
{
    s_point s = seeds->p[seeds_idx[id]];
    s_cid T0 = {TRI0, -1}, T1 = {TRI1, -1}, T2 = {TRI2, -1};

    /* Step 0: triangle degeneracy. Collinear vertices give orient3d = 0 for
     * every fourth point; at most two axis-offset candidates fit any plane,
     * so all three returning 0 is the exact collinearity condition. */
    {
        double ax = face[0].x, ay = face[0].y, az = face[0].z;
        double bx = face[1].x, by = face[1].y, bz = face[1].z;
        double cx = face[2].x, cy = face[2].y, cz = face[2].z;
        if (orient3d(ax,ay,az, bx,by,bz, cx,cy,cz, ax+1,ay,az) == 0 &&
            orient3d(ax,ay,az, bx,by,bz, cx,cy,cz, ax,ay+1,az) == 0 &&
            orient3d(ax,ay,az, bx,by,bz, cx,cy,cz, ax,ay,az+1) == 0)
            return false;
    }

    /* Step 1: initial candidate = T0 ∩ T1 (a triangle vertex). */
    s_cid A = T0, B = T1;

    /* Fisher-Yates shuffle. */
    int order[N];
    for (int i = 0; i < N; i++) order[i] = i;
    for (int i = N-1; i > 0; i--) {
        int j = randint(rctx, i+1);
        int tmp = order[i]; order[i] = order[j]; order[j] = tmp;
    }

    /* gate: seed 123 (x≈0.833) on the x=1 face (all verts have x=1) */
    for (int oi = 0; oi < N; oi++) {
        int i = order[oi];
        if (i == id) continue;
        s_cid ci = {CSEED, i};

        int outer_feas = lp_feasible(face, s, seeds, seeds_idx, A, B, ci);
        if (outer_feas >= 0) continue;

        /* 1D LP on C_i's line. lo/hi track interval endpoints; *_null = true
         * means that slot is unbounded (no constraint has closed it yet). */
        s_cid lo = T0, hi = T0; /* dummy init to silence uninit warnings */
        bool lo_null, hi_null;

        /* Step 3a: seed the interval from C_T0. */
        int d_i_T0 = lp_D(face, s, seeds, seeds_idx, ci, T0);
        if (d_i_T0 != 0) {
            /* Generic case: one slot from C_T0, other stays open. */
            if (d_i_T0 > 0) { lo = T0; lo_null = false; hi_null = true; }
            else             { hi = T0; hi_null = false; lo_null = true; }
            if (!lp_process_cj(face, s, seeds, seeds_idx, ci, T1,
                                &lo, &lo_null, &hi, &hi_null))
                return false;
        } else {
            /* Degenerate case: C_T0 is parallel to C_i.
             * Check whole line feasibility w.r.t. C_T0, then seed from C_T1. */
            int par_check = lp_feasible(face, s, seeds, seeds_idx, ci, T1, T0);
            if (par_check < 0) return false;
            int d_i_T1 = lp_D(face, s, seeds, seeds_idx, ci, T1);
            if (d_i_T1 > 0) { lo = T1; lo_null = false; hi_null = true; }
            else             { hi = T1; hi_null = false; lo_null = true; }
        }

        /* Step 3c: C_T2, then all SEED constraints before C_i in random order. */
        if (!lp_process_cj(face, s, seeds, seeds_idx, ci, T2,
                            &lo, &lo_null, &hi, &hi_null))
            return false;

        for (int oj = 0; oj < oi; oj++) {
            int jj = order[oj];
            if (jj == id) continue;
            s_cid cj = {CSEED, jj};
            if (!lp_process_cj(face, s, seeds, seeds_idx, ci, cj,
                                &lo, &lo_null, &hi, &hi_null))
                return false;
        }

        /* Step 4: new candidate is C_i ∩ lo. */
        A = ci;
        B = lo;
    }

    return true;
}


int extract_mirrored_points(const s_bpoly *bp, double EPS_degenerate,
                            s_scplx *dt, int N_seeds,
                            int (*randint)(void *rctx, int), void *rctx,
                            s_dynarray *out_points)
{   /* Use DT neighbor relationships as LP constraint set instead of all seeds.
     * LP with fewer constraints has a larger feasible region, so we may generate
     * extra mirrors (false positives) but never miss a needed one (no false negatives).
     *
     * Loop order: seed-outer, face-inner — neighbors are collected once per seed and
     * reused across all faces, rather than recomputed for every (face, seed) pair. */
    out_points->N = 0;

    /* seeds_idx holds [vi, neighbor_0, ..., neighbor_n] for the current seed.
     * Reused across seeds; grows automatically via dynarray. */
    s_dynarray seeds_idx  = dynarray_initialize(sizeof(int), 32);
    s_dynarray scratch_cells = dynarray_initialize(sizeof(s_ncell *), 64);
    if (!seeds_idx.items || !scratch_cells.items) goto error;

    for (int i = 0; i < N_seeds; i++) {
        int vi = 4 + i;
        if (!dt->point2tet[vi]) continue;  /* deduplicated seed */

        int vi_lid = -1;
        for (int k = 0; k < 4; k++)
            if (dt->point2tet[vi]->vertex_id[k] == vi) { vi_lid = k; break; }

        /* Build seeds_idx = [vi, neighbors...] for this seed. */
        seeds_idx.N = 0;
        if (!dynarray_push(&seeds_idx, &vi)) goto error;
        if (vertex_neighbors(dt, vi, dt->point2tet[vi], vi_lid,
                             4, &seeds_idx, &scratch_cells) < 0) goto error;

        /* seeds_idx.items is now a contiguous int[] usable by lp_should_mirror.
         * id=0 selects vi as the seed under test; all others are constraints. */
        int  n_total   = (int)seeds_idx.N;
        int *idx_array = (int *)seeds_idx.items;

        s_point p = dt->points.p[vi];

        for (int ff = 0; ff < bp->convh.Nf; ff++) {
            s_point face[3];
            convh_get_face(&bp->convh, ff, face);
            s_point n = bp->convh.fnormals[ff];
            double n2 = dot_prod(n, n);
            if (n2 < EPS_degenerate * EPS_degenerate) continue;

            if (!lp_should_mirror(face, &dt->points, n_total, idx_array, 0, randint, rctx))
                continue;

            double f = dot_prod(n, face[0]) - dot_prod(n, p);
            s_point p_mirror = { .x = p.x + 2*(f/n2)*n.x,
                                 .y = p.y + 2*(f/n2)*n.y,
                                 .z = p.z + 2*(f/n2)*n.z };
            if (!dynarray_push(out_points, &p_mirror)) goto error;
        }
    }

    dynarray_free(&seeds_idx);
    dynarray_free(&scratch_cells);
    return 1;

    error:
        dynarray_free(&seeds_idx);
        dynarray_free(&scratch_cells);
        fprintf(stderr, "Error in extend_sites_mirroring_dt!\n");
        return 0;
}

