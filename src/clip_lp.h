#ifndef CLIP_LP_H
#define CLIP_LP_H

/* -----------------------------------------------------------------------------
 * clip_lp.h -- 3D LP-feasibility dispatch layer for the non-convex clip.
 *
 * Header-only. Turns tagged constraint IDs (tet face T, or seed bisector S)
 * into exact sign queries by dispatching to the lp3_* predicates in
 * voronoi_predicates. This is the 3D analogue of the lp_D / lp_feasible
 * dispatchers in mirroring.c, and follows lp3_predicates.tex.
 *
 * A vertex of (cell ^ tet) is a triple of constraint boundaries; each is either
 *   T : a tetrahedron face  (n_T, defined by 3 of the tet's 4 vertices)
 *   S : a seed bisector      s->t  (favouring the cell seed s over neighbour t)
 * The feasibility of that vertex against a fourth constraint is
 *   feasible = sign(Delta) * sign(Gamma).
 *
 * Geometry sources:
 *   tet[4]     -- the four vertices of the clipping tetrahedron (T geometry)
 *   s          -- the cell's own seed (common endpoint of every S bisector)
 *   seeds->p[] -- neighbour seeds; a CT_SEED cid names its partner via .idx
 *   sigma_f[4] -- per-tet face interior-orientation signs (precomputed once by
 *                 lp3_tet_face_orient); consumed only by T-query cases.
 * --------------------------------------------------------------------------- */

#include "points.h"
#include "predicates_lp.h"
#include "predicates_lp_power.h" /* power_lp2_/power_lp3_ weighted dispatch */
#include "predicates_basic.h"   /* orient3d */
#include "predicates_power.h"   /* power_incircle3d: weighted co-orthocircularity */
#include "gtests.h"              /* test_insphere, test_orthosphere */

/* Constraint identifier. CT_TET: idx = face 0..3 (face f = the 3 verts != v[f]).
 * CT_SEED: idx = partner seed index into seeds->p[]. */
typedef enum { CT_TET, CT_SEED } e_ctype3;
typedef struct { e_ctype3 type; int idx; } s_cid3;

/* ---- tet combinatorics -------------------------------------------------- */

/* The three vertex indices of face f (the ones != f), ascending. */
static inline void lp3_face_idx(int f, int out[3])
{
    int k = 0;
    for (int i = 0; i < 4; i++) if (i != f) out[k++] = i;
}

/* The two shared-edge vertex indices of faces f and g (the ones != f and != g). */
static inline void lp3_shared_edge(int f, int g, int edge[2])
{
    int k = 0;
    for (int i = 0; i < 4; i++) if (i != f && i != g) edge[k++] = i;
}

/* The single vertex common to faces f,g,h == the remaining index (0+1+2+3=6). */
static inline int lp3_common_vertex(int f, int g, int h)
{
    return 6 - f - g - h;
}

/* Per-tet face interior-orientation signs. sigma_f[f] = orient3d(face f, v[f]):
 * the side of the opposite (interior) vertex relative to the raw face triple.
 * Multiplied into the T-query predicates so "interior" reads as feasible > 0. */
static inline void lp3_tet_face_orient(const s_point tet[4], int sigma_f[4])
{
    for (int f = 0; f < 4; f++) {
        int v[3]; lp3_face_idx(f, v);
        s_point A = tet[v[0]], B = tet[v[1]], C = tet[v[2]], W = tet[f];
        sigma_f[f] = orient3d(A.x,A.y,A.z, B.x,B.y,B.z, C.x,C.y,C.z, W.x,W.y,W.z);
    }
}

/* ---- exact concyclicity (cocircular 4 points) -------------------------- */

/* Do a,b,c,d lie on a common circle?  Concyclic = coplanar AND on a common
 * circle.  Gates the coplanarity precondition of incircle3d (valid only for
 * coplanar input): first an exact orient3d, then incircle3d.  a,b,c must be a
 * real triangle; a collinear triple has no circumcircle (incircle3d returns its
 * sentinel 2) and this reports 0.  Exact -- returns 1 iff concyclic, else 0. */
static inline int points_concyclic(s_point a, s_point b, s_point c, s_point d)
{
    if (test_orientation((s_point[]){ a, b, c }, d) != 0) return 0;   /* not coplanar */
    return incircle3d(a.x, a.y, a.z, b.x, b.y, b.z,
                      c.x, c.y, c.z, d.x, d.y, d.z) == 0;
}

/* Weighted (power) analogue: are a,b,c,d co-orthocircular?  Same coplanarity
 * gate (unweighted -- coplanarity is purely geometric); power_incircle3d ==
 * 0 iff d is on the weighted circle through (a,b,c).  wa,wb,wc,wd index-align
 * with a,b,c,d respectively. */
static inline int power_points_concyclic(s_point a, double wa, s_point b, double wb,
                                         s_point c, double wc, s_point d, double wd)
{
    if (test_orientation((s_point[]){ a, b, c }, d) != 0) return 0;   /* not coplanar */
    return power_incircle3d(a.x, a.y, a.z, wa, b.x, b.y, b.z, wb,
                            c.x, c.y, c.z, wc, d.x, d.y, d.z, wd) == 0;
}

/* ---- helpers to pull the S/T pieces out of a constraint triple ---------- */

/* Split (a,b,c) into the tet-face ids (T) and seed cids (S). */
static inline void lp3_split3(s_cid3 a, s_cid3 b, s_cid3 c,
                              int Tf[3], int *nT, s_cid3 Sc[3], int *nS)
{
    s_cid3 in[3] = { a, b, c };
    *nT = 0; *nS = 0;
    for (int i = 0; i < 3; i++) {
        if (in[i].type == CT_TET) Tf[(*nT)++] = in[i].idx;
        else                      Sc[(*nS)++] = in[i];
    }
}

/* ---- lp3_D : slope factor sign det[n1;n2;n3] ---------------------------- */
/* Dispatch by the type multiset of the triple. Used for slope decisions in a
 * Seidel-style solver; the enumerate-and-hull clip does not call it, but it is
 * provided for completeness and matches the lp3_D row of the catalogue.       */
static inline int lp3_D(const s_point tet[4], s_point s, const s_points *seeds,
                        s_cid3 a, s_cid3 b, s_cid3 c)
{
    int Tf[3], nT, nS; s_cid3 Sc[3];
    lp3_split3(a, b, c, Tf, &nT, Sc, &nS);

    if (nT == 0) {                                   /* SSS: orient3d(s,t1,t2,t3) */
        s_point t1 = seeds->p[Sc[0].idx], t2 = seeds->p[Sc[1].idx], t3 = seeds->p[Sc[2].idx];
        return orient3d(t1.x,t1.y,t1.z, t2.x,t2.y,t2.z, t3.x,t3.y,t3.z, s.x,s.y,s.z);
    }
    if (nT == 1) {                                   /* TSS: lp3_D_TSS */
        int fv[3]; lp3_face_idx(Tf[0], fv);
        s_point A = tet[fv[0]], Q = tet[fv[1]], R = tet[fv[2]];
        s_point t1 = seeds->p[Sc[0].idx], t2 = seeds->p[Sc[1].idx];
        return lp3_D_TSS(A.x,A.y,A.z, Q.x,Q.y,Q.z, R.x,R.y,R.z,
                         s.x,s.y,s.z, t1.x,t1.y,t1.z, t2.x,t2.y,t2.z);
    }
    if (nT == 2) {                                   /* TTS: (t-s).(B-A) on shared edge */
        int e[2]; lp3_shared_edge(Tf[0], Tf[1], e);
        s_point A = tet[e[0]], B = tet[e[1]], t = seeds->p[Sc[0].idx];
        /* sign((t-s).(B-A)) via lp2_D_T1_S(A,B,s,t); per-edge orientation sign is
         * applied by the caller when this feeds a 1D-edge lo/hi decision. */
        return lp2_D_T1_S(A.x,A.y,A.z, B.x,B.y,B.z, s.x,s.y,s.z, t.x,t.y,t.z);
    }
    /* nT == 3, TTT: determinant of three face normals -- a fixed per-tet
     * orientation sign. Not used by the enumerate-and-hull clip. */
    return orient3d(tet[0].x,tet[0].y,tet[0].z, tet[1].x,tet[1].y,tet[1].z,
                    tet[2].x,tet[2].y,tet[2].z, tet[3].x,tet[3].y,tet[3].z);
}

/* ---- lp3_feasible : feasibility of vertex (a^b^c) against query d -------- */
/* Returns >0 (strictly inside d), 0 (on d's boundary / degenerate), <0 (outside).
 * For T queries the raw predicate is multiplied by sigma_f[query face] so that
 * the tet interior is the feasible side.
 *
 * seed_weights: NULL-able, index-aligned with seeds->p[] (NULL = unweighted,
 * matching the dt_builder_begin idiom). s_w is the cell's own seed s's weight
 * (meaningless/unused when seed_weights is NULL). When non-NULL, every
 * S-involving branch dispatches to the weighted power_lp2_/power_lp3_
 * predicate instead of its unweighted counterpart; the two purely
 * combinatorial branches (TTS vs T, TTT vs T) are weight-independent and
 * untouched -- see the Phase 2 predicate table. */
static inline int lp3_feasible(const s_point tet[4], const int sigma_f[4], s_point s,
                               const s_points *seeds,
                               const double *seed_weights, double s_w,
                               s_cid3 a, s_cid3 b, s_cid3 c, s_cid3 d)
{
    int Tf[3], nT, nS; s_cid3 Sc[3];
    lp3_split3(a, b, c, Tf, &nT, Sc, &nS);

    if (nT == 0) {                                   /* ---- SSS vertex ---- */
        s_point t1 = seeds->p[Sc[0].idx], t2 = seeds->p[Sc[1].idx], t3 = seeds->p[Sc[2].idx];
        if (d.type == CT_SEED) {                     /* SSS vs S == insphere */
            /* circumcenter of (s,t1,t2,t3); feasible iff u outside circumsphere.
             * (In a valid DT this never rejects.) Sign validated by
             * tests/clip_lp_validate.c. */
            s_point sph[4] = { s, t1, t2, t3 };
            s_point u = seeds->p[d.idx];
            if (seed_weights) {
                double sph_w[4] = { s_w, seed_weights[Sc[0].idx], seed_weights[Sc[1].idx],
                                     seed_weights[Sc[2].idx] };
                double w_u = seed_weights[d.idx];
                return -test_orthosphere(4, sph, sph_w, u, w_u);
            }
            return -test_insphere(sph, u);
        } else {                                     /* SSS vs T */
            int fv[3]; lp3_face_idx(d.idx, fv);
            s_point A = tet[fv[0]], Q = tet[fv[1]], R = tet[fv[2]];
            if (seed_weights) {
                double wt1 = seed_weights[Sc[0].idx], wt2 = seed_weights[Sc[1].idx],
                       wt3 = seed_weights[Sc[2].idx];
                return sigma_f[d.idx] *
                       power_lp3_feasible_SSS_T(s.x,s.y,s.z, t1.x,t1.y,t1.z, t2.x,t2.y,t2.z,
                                                t3.x,t3.y,t3.z, A.x,A.y,A.z, Q.x,Q.y,Q.z, R.x,R.y,R.z,
                                                s_w, wt1, wt2, wt3);
            }
            return sigma_f[d.idx] *
                   lp3_feasible_SSS_T(s.x,s.y,s.z, t1.x,t1.y,t1.z, t2.x,t2.y,t2.z,
                                      t3.x,t3.y,t3.z, A.x,A.y,A.z, Q.x,Q.y,Q.z, R.x,R.y,R.z);
        }
    }

    if (nT == 1) {                                   /* ---- TSS vertex ---- */
        int f = Tf[0];
        s_point t1 = seeds->p[Sc[0].idx], t2 = seeds->p[Sc[1].idx];
        if (d.type == CT_SEED) {                     /* TSS vs S */
            int fv[3]; lp3_face_idx(f, fv);
            s_point A = tet[fv[0]], Q = tet[fv[1]], R = tet[fv[2]], u = seeds->p[d.idx];
            if (seed_weights) {
                double wt1 = seed_weights[Sc[0].idx], wt2 = seed_weights[Sc[1].idx],
                       wu = seed_weights[d.idx];
                return power_lp3_feasible_TSS_S(A.x,A.y,A.z, Q.x,Q.y,Q.z, R.x,R.y,R.z,
                                                s.x,s.y,s.z, t1.x,t1.y,t1.z, t2.x,t2.y,t2.z,
                                                u.x,u.y,u.z, s_w, wt1, wt2, wu);
            }
            return lp3_feasible_TSS_S(A.x,A.y,A.z, Q.x,Q.y,Q.z, R.x,R.y,R.z,
                                      s.x,s.y,s.z, t1.x,t1.y,t1.z, t2.x,t2.y,t2.z,
                                      u.x,u.y,u.z);
        } else {                                     /* TSS vs T (query face d) */
            int g = d.idx;
            int e[2]; lp3_shared_edge(f, g, e);
            /* T1 (definer) = face f, apex = tet[g];  T2 (query) = face g, apex = tet[f]. */
            s_point A = tet[e[0]], E = tet[e[1]], C = tet[g], D = tet[f];
            /* Orientation of the query face T2 must match the predicate's own
             * n_T2 = (E-A)x(D-A) vertex order, so compute it inline rather than
             * from the ascending-order sigma_f[] (which can differ by a swap). */
            int sq = orient3d(A.x,A.y,A.z, E.x,E.y,E.z, D.x,D.y,D.z,
                              tet[g].x,tet[g].y,tet[g].z);
            if (seed_weights) {
                double wt1 = seed_weights[Sc[0].idx], wt2 = seed_weights[Sc[1].idx];
                return sq *
                       power_lp3_feasible_TSS_T(A.x,A.y,A.z, E.x,E.y,E.z, C.x,C.y,C.z, D.x,D.y,D.z,
                                                s.x,s.y,s.z, t1.x,t1.y,t1.z, t2.x,t2.y,t2.z,
                                                s_w, wt1, wt2);
            }
            return sq *
                   lp3_feasible_TSS_T(A.x,A.y,A.z, E.x,E.y,E.z, C.x,C.y,C.z, D.x,D.y,D.z,
                                      s.x,s.y,s.z, t1.x,t1.y,t1.z, t2.x,t2.y,t2.z);
        }
    }

    if (nT == 2) {                                   /* ---- TTS vertex (edge) ---- */
        int e[2]; lp3_shared_edge(Tf[0], Tf[1], e);
        s_point A = tet[e[0]], B = tet[e[1]];
        if (d.type == CT_SEED) {                     /* TTS vs S */
            s_point t = seeds->p[Sc[0].idx], u = seeds->p[d.idx];
            if (seed_weights) {
                double wt = seed_weights[Sc[0].idx], wu = seed_weights[d.idx];
                return power_lp3_feasible_TTS_S(A.x,A.y,A.z, B.x,B.y,B.z, s.x,s.y,s.z,
                                                t.x,t.y,t.z, u.x,u.y,u.z, s_w, wt, wu);
            }
            return lp3_feasible_TTS_S(A.x,A.y,A.z, B.x,B.y,B.z, s.x,s.y,s.z,
                                      t.x,t.y,t.z, u.x,u.y,u.z);
        } else {                                     /* TTS vs T : automatic */
            /* Crossing lies on tet edge (faces Tf[0],Tf[1]); it is on a query
             * face iff that face is one of the two, else strictly interior.
             * Purely combinatorial (no bisector value involved) -> weight-
             * independent, unchanged whether or not seed_weights is given. */
            return (d.idx == Tf[0] || d.idx == Tf[1]) ? 0 : 1;
        }
    }

    /* nT == 3, TTT vertex == the tet vertex common to faces Tf[0..2]. */
    {
        int V = lp3_common_vertex(Tf[0], Tf[1], Tf[2]);
        if (d.type == CT_SEED) {                     /* TTT vs S : bisector value at V */
            s_point Vp = tet[V], t = seeds->p[d.idx];
            /* lp2_feasible_T0_T1_S(V,s,t) = sign(|V-t|^2 - |V-s|^2) = feasibility
             * (>0 when V is strictly closer to s). */
            if (seed_weights) {
                double wt = seed_weights[d.idx];
                return power_lp2_feasible_T0_T1_S(Vp.x,Vp.y,Vp.z, s.x,s.y,s.z, t.x,t.y,t.z,
                                                  s_w, wt);
            }
            return lp2_feasible_T0_T1_S(Vp.x,Vp.y,Vp.z, s.x,s.y,s.z, t.x,t.y,t.z);
        } else {                                     /* TTT vs T : combinatorial */
            return (d.idx == V) ? 1 : 0;             /* strictly inside only the opposite face,
                                                       * weight-independent */
        }
    }
}

/* ---- cheap frontier feasibility filter --------------------------------------
 * Exact one-sided test used to prune the flood frontier before a full clip.
 *
 * cell_i = intersection over neighbours t of the half-space {x closer to s than
 * t}. If cell_i ^ tet has a point p, then for every bisector s->t that p lies on
 * the s-side, so the s-side half-space meets the tet; a half-space meets a convex
 * polytope iff it contains a vertex of it (the linear bisector value is minimised
 * at a vertex). Hence:
 *     cell_i ^ tet != empty  ==>  every bisector s->t has some tet vertex on the
 *                                 s-side.
 * The contrapositive is exact: if any bisector has all four tet vertices strictly
 * on the t-side, cell_i ^ tet is empty -- safe to reject (frontier). Passing is
 * only necessary, so a survivor still gets the full clip. Uses only the validated
 * lp2_feasible_T0_T1_S (bisector value sign), no lp3_D and no orientation data.
 * Returns 1 (maybe overlaps -> clip), 0 (definitely disjoint -> frontier).
 *
 * seed_weights: NULL-able, index-aligned with seeds->p[] (NULL = unweighted).
 * s_w is the cell seed s's own weight (unused when seed_weights is NULL).
 * The "minimized at a vertex" soundness argument above is stated for the
 * bisector's own favoring test, whatever that is; a power bisector is still
 * an affine half-space (a_S untouched, only its offset shifts -- Phase 2),
 * so the identical argument holds verbatim once every "s-side" check below
 * uses the WEIGHTED favoring test instead of the unweighted one. Without
 * this, the filter would test cell membership against the wrong (unweighted)
 * half-space when the cell is actually weighted, which could wrongly reject
 * (as "frontier") a cell that genuinely overlaps the tet under the true
 * weighted metric -- a real correctness gap, not just a missed
 * optimization. */
static inline int lp3_cell_maybe_hits_tet(const s_point tet[4], s_point s,
                                          const s_points *seeds,
                                          const double *seed_weights, double s_w,
                                          const int *nbr, int NB)
{
    for (int j = 0; j < NB; j++) {
        s_point t = seeds->p[nbr[j]];
        int on_s_side = 0;
        for (int k = 0; k < 4; k++) {
            int fe = seed_weights
                ? power_lp2_feasible_T0_T1_S(tet[k].x,tet[k].y,tet[k].z,
                                             s.x,s.y,s.z, t.x,t.y,t.z, s_w, seed_weights[nbr[j]])
                : lp2_feasible_T0_T1_S(tet[k].x,tet[k].y,tet[k].z,
                                      s.x,s.y,s.z, t.x,t.y,t.z);
            if (fe >= 0) { on_s_side = 1; break; }
        }
        if (!on_s_side) return 0;               /* this bisector cuts the whole tet away */
    }
    return 1;
}

/* Same one-sided filter against a triangle (a boundary face). Used to classify
 * which cells reach the domain surface (Part II interior-cell shortcut): passing
 * is necessary for cell_i to touch the triangle, so a cell that truly touches the
 * boundary is never rejected -> it can never be mis-classified as interior. */
static inline int lp3_cell_maybe_hits_tri(const s_point tri[3], s_point s,
                                          const s_points *seeds,
                                          const double *seed_weights, double s_w,
                                          const int *nbr, int NB)
{
    for (int j = 0; j < NB; j++) {
        s_point t = seeds->p[nbr[j]];
        int on_s_side = 0;
        for (int k = 0; k < 3; k++) {
            int fe = seed_weights
                ? power_lp2_feasible_T0_T1_S(tri[k].x,tri[k].y,tri[k].z,
                                             s.x,s.y,s.z, t.x,t.y,t.z, s_w, seed_weights[nbr[j]])
                : lp2_feasible_T0_T1_S(tri[k].x,tri[k].y,tri[k].z,
                                      s.x,s.y,s.z, t.x,t.y,t.z);
            if (fe >= 0) { on_s_side = 1; break; }
        }
        if (!on_s_side) return 0;
    }
    return 1;
}

#endif /* CLIP_LP_H */
