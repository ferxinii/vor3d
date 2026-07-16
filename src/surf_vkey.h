/*
 * Combinatorial vertex keys + per-cell surface assembly for non-convex
 * Voronoi-cell boundary extraction.  See SURFACE_EXTRACTION_PLAN.md.
 *
 * A surface vertex is identified by the exact combination of constraint planes
 * that meet there, using GLOBAL ids (CDT vertex ids, real-seed ids) that are
 * identical no matter which tet or which cell computes the point.  Two floats
 * for the "same" seam vertex, computed independently in two adjacent pieces,
 * therefore carry the same key and stitch exactly -- no coordinate tolerance.
 *
 * Header-only (static inline), same style as dynarray.h / hash.h.
 */
#ifndef VOR3D_SURF_VKEY_H
#define VOR3D_SURF_VKEY_H

#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include "points.h"
#include "trimesh.h"
#include "dynarray.h"
#include "hash.h"
#include "UF.h"         /* Union-Find for connected-component splitting */
#include "gtests.h"     /* exact orient2d + test_insphere */
#include "voronoi_predicates.h"   /* exact lp3_* feasibility for edge ordering */
#include "robust_predicates.h"    /* orient3d */

/* ---- key type ---------------------------------------------------------- */

typedef enum { VK_TV = 0, VK_VB, VK_VV, VK_EF, VK_EE } e_vkind;

/* Fixed-size POD, no padding (7 x int32 = 28 bytes): safe to memcmp / hash by
 * bytes.  Unused fields are -1.  Fields are canonicalised (sorted) so equal
 * points always produce byte-identical keys. */
typedef struct vkey {
    int32_t tag;
    int32_t a, b, c, d, e, f;
} s_vkey;

/* Per-hull-vertex label: identity key + plane membership (for face classify).
 * fmask bit i set => vertex lies on local tet face i (0..3).
 * bis[k] = partner seed id of a bisector {s, bis[k]} the vertex lies on (-1 unused).
 *
 * Degenerate (cospherical) seed sets make several candidate constructions land
 * on the SAME point (e.g. a grid-cube corner is the circumcentre of every tet
 * triangulating the cube: up to ~6 VV duals, possibly plus EF/VB crossings).
 * When the hull dedupes them into one vertex the label must not lose plane
 * evidence, so it keeps the UNION of bisector partners (SV_NBIS slots) and all
 * distinct merged keys: `key` is the canonical (lexicographically smallest)
 * one -- deterministic, so every piece welds the vertex identically -- and
 * alt[0..nalt) are the others, consulted for complete plane membership. */
#define SV_NBIS 12
#define SV_NALT 8
typedef struct vlabel {
    s_vkey  key;                 /* canonical identity (min over merged keys) */
    uint8_t fmask;
    int8_t  nalt;
    int32_t bis[SV_NBIS];
    s_vkey  alt[SV_NALT];        /* other merged constructions of this point */
} s_vlabel;

/* One triangle corner emitted to a cell's surface accumulator. */
typedef struct surf_vtx {
    s_vlabel lb;
    s_point  p;
} s_surf_vtx;

/* Canonical id of the supporting plane of a piece face, using GLOBAL ids so the
 * two tets sharing an interior face (and the several tet-pieces spanning one
 * bisector facet) land on the SAME plane. tag: PK_BIS -> a,b = sorted seed pair;
 * PK_FACE -> a,b,c = sorted CDT face-triangle vertex ids (d unused). */
typedef enum { PK_BIS = 0, PK_FACE } e_pkind;
typedef struct pkey { int32_t tag, a, b, c; } s_pkey;   /* 16 bytes, no padding */

/* One piece face emitted to the accumulator: its supporting plane + 3 corners. */
typedef struct surf_tri {
    s_pkey     plane;
    s_surf_vtx v[3];
} s_surf_tri;


/* ---- small canonicalisation helpers ------------------------------------ */

static inline void vk_sort2(int32_t *x, int32_t *y)
{
    if (*x > *y) { int32_t t = *x; *x = *y; *y = t; }
}
static inline void vk_sort3(int32_t *x, int32_t *y, int32_t *z)
{
    vk_sort2(x, y); vk_sort2(y, z); vk_sort2(x, y);
}
static inline void vk_sort4(int32_t *w, int32_t *x, int32_t *y, int32_t *z)
{
    vk_sort2(w, x); vk_sort2(y, z); vk_sort2(w, y); vk_sort2(x, z); vk_sort2(x, y);
}

static inline s_pkey pkey_bis(int s, int t)
{
    int32_t a = s, b = t; vk_sort2(&a, &b);
    return (s_pkey){ PK_BIS, a, b, -1 };
}
static inline s_pkey pkey_face(int v0, int v1, int v2)
{
    int32_t a = v0, b = v1, c = v2; vk_sort3(&a, &b, &c);
    return (s_pkey){ PK_FACE, a, b, c };
}


/* ---- key constructors -------------------------------------------------- */

/* (b) tet vertex: a single global CDT vertex id. */
static inline s_vkey vkey_TV(int v)
{
    return (s_vkey){ VK_TV, v, -1, -1, -1, -1, -1 };
}

/* (d) tet-edge ^ bisector: CDT edge {e0,e1} (global) x seed pair {s,t}. */
static inline s_vkey vkey_VB(int e0, int e1, int s, int t)
{
    int32_t a = e0, b = e1, c = s, d = t;
    vk_sort2(&a, &b); vk_sort2(&c, &d);
    return (s_vkey){ VK_VB, a, b, c, d, -1, -1 };
}

/* (a) Voronoi vertex: circumcentre of DT tet {s,j1,j2,j3} (global seed ids). */
static inline s_vkey vkey_VV(int s, int j1, int j2, int j3)
{
    int32_t a = s, b = j1, c = j2, d = j3;
    vk_sort4(&a, &b, &c, &d);
    return (s_vkey){ VK_VV, a, b, c, d, -1, -1 };
}

/* (c) Voronoi-edge {s,j1,j2} ^ tet-face triangle {t0,t1,t2} (all global). */
static inline s_vkey vkey_EF(int s, int j1, int j2, int t0, int t1, int t2)
{
    int32_t a = s, b = j1, c = j2, d = t0, e = t1, f = t2;
    vk_sort3(&a, &b, &c); vk_sort3(&d, &e, &f);
    return (s_vkey){ VK_EF, a, b, c, d, e, f };
}

/* (c') Voronoi-edge {s,j1,j2} ^ tet-EDGE {e0,e1} (all global): the crossing
 * lies exactly ON a CDT edge (a face-feasibility test returned 0), so the
 * face-based EF key would differ between the two faces sharing the edge.  A
 * line ^ line crossing is unique, so this form is identical in every tet and
 * face around that CDT edge. */
static inline s_vkey vkey_EE(int s, int j1, int j2, int e0, int e1)
{
    int32_t a = s, b = j1, c = j2, d = e0, e = e1;
    vk_sort3(&a, &b, &c); vk_sort2(&d, &e);
    return (s_vkey){ VK_EE, a, b, c, d, e, -1 };
}


/* ---- hash.h callbacks -------------------------------------------------- */

static inline size_t vkey_hash(const void *k)
{
    const unsigned char *p = (const unsigned char *)k;
    size_t h = 1469598103934665603ULL;              /* FNV-1a 64 */
    for (size_t i = 0; i < sizeof(s_vkey); i++) { h ^= p[i]; h *= 1099511628211ULL; }
    return h;
}
static inline bool vkey_equals(const void *a, const void *b)
{
    return memcmp(a, b, sizeof(s_vkey)) == 0;
}


/* ---- label merge (two candidates deduped to one hull vertex) ----------- */

/* Total order on vkeys (lexicographic on the canonical index fields): picks the
 * canonical key of a merged (degenerate) vertex, and is the symbolic tiebreak
 * for the exact edge ordering.  Identical everywhere, so every piece/facet
 * makes the same choice. */
static inline int sv_vkey_cmp(s_vkey A, s_vkey B)
{
    if (A.tag != B.tag) return A.tag < B.tag ? -1 : 1;
    if (A.a != B.a) return A.a < B.a ? -1 : 1;
    if (A.b != B.b) return A.b < B.b ? -1 : 1;
    if (A.c != B.c) return A.c < B.c ? -1 : 1;
    if (A.d != B.d) return A.d < B.d ? -1 : 1;
    if (A.e != B.e) return A.e < B.e ? -1 : 1;
    if (A.f != B.f) return A.f < B.f ? -1 : 1;
    return 0;
}

/* Label for a fresh candidate: canonical key `k`, tet-face mask, and up to
 * three bisector partners (-1 = unused). Remaining slots are cleared. */
static inline s_vlabel vlabel_make3(s_vkey k, uint8_t fmask, int t1, int t2, int t3)
{
    s_vlabel L;
    L.key = k; L.fmask = fmask; L.nalt = 0;
    L.bis[0] = t1; L.bis[1] = t2; L.bis[2] = t3;
    for (int i = 3; i < SV_NBIS; i++) L.bis[i] = -1;
    memset(L.alt, 0, sizeof(L.alt));
    return L;
}

static inline void vlabel_add_bis(s_vlabel *dst, int t)
{
    if (t < 0) return;
    for (int k = 0; k < SV_NBIS; k++) if (dst->bis[k] == t) return;   /* already present */
    for (int k = 0; k < SV_NBIS; k++) if (dst->bis[k] < 0) { dst->bis[k] = t; return; }
    fprintf(stderr, "vlabel_add_bis: >%d bisector partners at one vertex (dropped %d)\n",
            SV_NBIS, t);
}

static inline void vlabel_add_alt(s_vlabel *dst, s_vkey k)
{
    if (sv_vkey_cmp(k, dst->key) == 0) return;
    for (int i = 0; i < dst->nalt; i++) if (sv_vkey_cmp(k, dst->alt[i]) == 0) return;
    if (dst->nalt < SV_NALT) { dst->alt[dst->nalt++] = k; return; }
    fprintf(stderr, "vlabel_add_alt: >%d coincident constructions at one vertex\n", SV_NALT);
}

/* Coincident candidates deduped to one hull vertex: union the plane evidence
 * and keep the smallest key as the canonical identity (order-independent, so
 * every piece containing the point picks the same one). */
static inline void vlabel_merge(s_vlabel *dst, const s_vlabel *src)
{
    dst->fmask |= src->fmask;
    for (int k = 0; k < SV_NBIS; k++) vlabel_add_bis(dst, src->bis[k]);
    if (sv_vkey_cmp(src->key, dst->key) < 0) {
        vlabel_add_alt(dst, dst->key);
        dst->key = src->key;
    } else {
        vlabel_add_alt(dst, src->key);
    }
    for (int i = 0; i < src->nalt; i++) vlabel_add_alt(dst, src->alt[i]);
}

/* Smallest partner seed id present in all three vertices' bisector lists (-1 if
 * none): the three vertices of a hull face all lie on the bisector
 * {s, that_partner}. Smallest => deterministic across the pieces tiling a
 * facet (a non-degenerate triangle determines its plane, so distinct common
 * partners can only name the same plane; the tie must still break identically). */
static inline int vlabel_common_bis(const s_vlabel *a, const s_vlabel *b,
                                    const s_vlabel *c)
{
    int best = -1;
    for (int i = 0; i < SV_NBIS; i++) {
        int t = a->bis[i];
        if (t < 0) continue;
        int inb = 0, inc = 0;
        for (int k = 0; k < SV_NBIS; k++) {
            if (b->bis[k] == t) inb = 1;
            if (c->bis[k] == t) inc = 1;
        }
        if (inb && inc && (best < 0 || t < best)) best = t;
    }
    return best;
}


/* ---- per-cell assembly: per-plane edge XOR + combinatorial subdivision -----
 *
 * Every piece face carries its supporting plane (s_pkey). For each plane we XOR
 * the DIRECTED edges of its triangles (an edge and its reverse cancel):
 *   - interior triangulation diagonals cancel within a piece;
 *   - an interior tet face is contributed by both adjacent tets with opposite
 *     orientation -> its whole boundary cancels -> the plane vanishes;
 *   - a boundary tet face or a bisector facet (possibly spanning several tet
 *     pieces) leaves its outer boundary loop.
 * Before XOR, every triangle edge is SUBDIVIDED at the vertices lying on it,
 * detected COMBINATORIALLY: an edge (u,v) on facet P and neighbour plane Q has
 * supporting line P^Q, and w lies on it iff on_plane(w,P) && on_plane(w,Q) (its
 * vkey places it on both). Both facets sharing the line subdivide identically,
 * so no T-junction can open a seam -- the identification is exact, never a
 * coordinate test. (The position of w along the edge -- which sub-edges to emit
 * -- currently uses a float axis order; SURFACE_EXACT_PLAN.md Phase B replaces it
 * with the exact `ord` predicate composed from lp3_feasible_*.) The surviving
 * edges per plane are walked into closed loops and ear-clipped (exact orient2d). */

/* Axis-aligned projection to 2D: keep the two coordinates != drop, in an order
 * that preserves orientation for the +drop normal. Coordinates are exact (a
 * subset of the input doubles) -> orient2d on them is exact. */
static inline s_point2d surf_proj2d(s_point p, int drop)
{
    switch (drop) {
        case 0:  return (s_point2d){{{ p.y, p.z }}};
        case 1:  return (s_point2d){{{ p.z, p.x }}};
        default: return (s_point2d){{{ p.x, p.y }}};
    }
}

/* ---- exact clockwise angular order (for face tracing) --------------------
 * Rank of direction V->X in a CLOCKWISE sweep starting just after direction
 * V->U: right half-plane first (0), exactly opposite (1), left half-plane (2),
 * exactly the incoming direction (U-turn) last (3).  Exact orient2d on the
 * projected doubles; the collinear same/opposite split compares signs on the
 * dominant axis of U-V (exact float compares). */
static inline int sv_dir_rank(s_point2d V, s_point2d U, s_point2d X)
{
    int o = test_orientation_2d((s_point2d[]){ V, U }, X);
    if (o < 0) return 0;
    if (o > 0) return 2;
    double du = U.x - V.x, dx = X.x - V.x;
    if (fabs(U.y - V.y) > fabs(du)) { du = U.y - V.y; dx = X.y - V.y; }
    return ((du > 0) == (dx > 0)) ? 3 : 1;
}

/* Does direction V->A come strictly before V->B in that clockwise sweep?
 * Same-rank candidates order by a second exact orient2d; exactly-parallel
 * same-direction candidates (post-subdivision leftovers) order closer-first. */
static inline int sv_cw_before(s_point2d V, s_point2d U, s_point2d A, s_point2d B)
{
    int ra = sv_dir_rank(V, U, A), rb = sv_dir_rank(V, U, B);
    if (ra != rb) return ra < rb;
    int o = test_orientation_2d((s_point2d[]){ V, A }, B);
    if (o != 0) return o < 0;
    return (fabs(A.x - V.x) + fabs(A.y - V.y)) < (fabs(B.x - V.x) + fabs(B.y - V.y));
}

/* Is X strictly inside the circumcircle of (A,B,C)?  Exact, orientation-agnostic.
 * Returns 0 for a collinear (degenerate) triangle: it has no circumcircle. */
static inline int sv_in_circumcircle(s_point2d A, s_point2d B, s_point2d C, s_point2d X)
{
    int o = test_orientation_2d((s_point2d[]){ A, B }, C);
    if (o == 0) return 0;
    s_point2d tri[3];
    if (o > 0) { tri[0] = A; tri[1] = B; tri[2] = C; }   /* test_incircle wants CCW */
    else       { tri[0] = A; tri[1] = C; tri[2] = B; }
    return test_incircle(tri, X) > 0;
}

/* One directed triangle edge, for pairing triangles across shared edges. */
typedef struct { int lo, hi, t, c; } s_sv_tedge;  /* c = corner OPPOSITE the edge */

static inline int sv_tedge_cmp(const void *a, const void *b)
{
    const s_sv_tedge *x = (const s_sv_tedge *)a, *y = (const s_sv_tedge *)b;
    if (x->lo != y->lo) return x->lo - y->lo;
    return x->hi - y->hi;
}

/* Lawson-flip the INTERIOR diagonals of a triangulated simple polygon until no
 * edge is locally non-Delaunay.  tris holds ntri triangles as local ring indices
 * into P, all wound consistently; both are modified in place.
 *
 * Boundary edges are shared by only one triangle here, so they are never
 * flipped: the polygon's region, its boundary and its vertex set are all
 * preserved (Steiner-free -- neighbouring facets keep conforming).  What changes
 * is only WHICH diagonals are used.
 *
 * This matters because ear clipping alone produces arbitrarily bad triangles: on
 * a Voronoi cut face the clip stalls (collinear ring vertices are never ears)
 * and the fan fallback below then fans the rest from a single hub -- measured on
 * one cell: one vertex in 111 of 151 triangles, needles down to 0.0016 deg, and
 * triangles whose apex sits ~1e-15 from the opposite edge.  Segment recovery in
 * the CDT cannot converge on that (it would need Steiner points spaced ~1e-15
 * apart).  Flipping to Delaunay gives the best-conditioned triangulation of this
 * vertex set.  A flip is taken when the edge is non-Delaunay OR when it borders
 * a zero-area triangle, so degenerate slivers are actively removed.
 * Requires a strictly convex quad, so the flip is always geometrically valid. */
static inline void sv_flip_to_delaunay(const s_point2d *P, int *tris, int ntri)
{
    if (ntri < 2) return;

    s_sv_tedge *E    = (s_sv_tedge *)malloc(sizeof(s_sv_tedge) * (size_t)(3 * ntri));
    char       *hot  = (char *)malloc((size_t)ntri);
    if (!E || !hot) { free(E); free(hot); return; }   /* keep the ear-clip result */

    const int max_rounds = 4 * ntri + 8;
    for (int round = 0; round < max_rounds; round++) {
        int ne = 0;
        for (int t = 0; t < ntri; t++)
            for (int c = 0; c < 3; c++) {
                int u = tris[3*t + (c+1)%3], v = tris[3*t + (c+2)%3];
                E[ne].lo = u < v ? u : v;
                E[ne].hi = u < v ? v : u;
                E[ne].t  = t;
                E[ne].c  = c;
                ne++;
            }
        qsort(E, (size_t)ne, sizeof(s_sv_tedge), sv_tedge_cmp);
        memset(hot, 0, (size_t)ntri);

        int nflips = 0;
        for (int i = 0; i + 1 < ne; i++) {
            if (E[i].lo != E[i+1].lo || E[i].hi != E[i+1].hi) continue;  /* boundary edge */

            int t1 = E[i].t,   c1 = E[i].c;
            int t2 = E[i+1].t, c2 = E[i+1].c;
            if (hot[t1] || hot[t2]) continue;      /* already re-shaped this round */

            /* t1 carries the edge DIRECTED u1 -> v1; t2 carries it reversed. */
            int u1 = tris[3*t1 + (c1+1)%3];
            int v1 = tris[3*t1 + (c1+2)%3];
            int w1 = tris[3*t1 + c1];
            int w2 = tris[3*t2 + c2];

            /* The flip replaces diagonal (u1,v1) with (w1,w2); it is valid only
             * if the quad is STRICTLY convex, i.e. u1 and v1 lie on opposite
             * sides of the new diagonal. */
            int a1 = test_orientation_2d((s_point2d[]){ P[w1], P[w2] }, P[u1]);
            int a2 = test_orientation_2d((s_point2d[]){ P[w1], P[w2] }, P[v1]);
            if (a1 == 0 || a2 == 0 || (a1 > 0) == (a2 > 0)) continue;

            int d1 = test_orientation_2d((s_point2d[]){ P[u1], P[v1] }, P[w1]);
            int d2 = test_orientation_2d((s_point2d[]){ P[u1], P[v1] }, P[w2]);
            int degenerate = (d1 == 0 || d2 == 0);   /* zero-area triangle: flip it away */

            if (!degenerate &&
                !sv_in_circumcircle(P[u1], P[v1], P[w1], P[w2])) continue;  /* Delaunay */

            tris[3*t1 + 0] = u1; tris[3*t1 + 1] = w2; tris[3*t1 + 2] = w1;
            tris[3*t2 + 0] = v1; tris[3*t2 + 1] = w1; tris[3*t2 + 2] = w2;
            hot[t1] = hot[t2] = 1;
            nflips++;
        }
        if (nflips == 0) break;
    }

    free(E); free(hot);
}

/* Triangulate a simple loop (vertex ids into `coords`, length n) into `faces_out`
 * (int[3] per triangle): ear-clip for the topology, then Lawson-flip the interior
 * diagonals to Delaunay for the shape.  Exact orient2d/incircle decisions.  Falls
 * back to a fan if no ear is found (still Steiner-free -> still closed; the flip
 * pass then repairs the fan).  Returns 1 OK, -1 alloc error. */
static inline int surf_earclip(const s_point *coords, const int *loop, int n,
                               s_dynarray *faces_out)
{
    if (n < 3) return 1;

    /* normal (FP) only picks the drop axis; the projection itself is exact. */
    s_point nrm = {{{ 0, 0, 0 }}};
    for (int i = 1; i + 1 < n; i++) {
        s_point a = coords[loop[0]], b = coords[loop[i]], c = coords[loop[i+1]];
        s_point cr = cross_prod(subtract_points(b, a), subtract_points(c, a));
        if (norm(cr) > norm(nrm)) nrm = cr;
    }
    int drop = coord_with_largest_component_3D(nrm);

    s_point2d *P    = (s_point2d *)malloc(sizeof(s_point2d) * (size_t)n);
    int       *V    = (int *)malloc(sizeof(int) * (size_t)n);   /* active ring, indices 0..n-1 */
    int       *tris = (int *)malloc(sizeof(int) * 3 * (size_t)(n - 2));  /* local indices */
    if (!P || !V || !tris) { free(P); free(V); free(tris); return -1; }
    for (int i = 0; i < n; i++) { P[i] = surf_proj2d(coords[loop[i]], drop); V[i] = i; }
    int m = n, nt = 0, rc = 1;

    /* Polygon orientation from its lexicographically-lowest vertex (always a
     * strictly convex corner): exact via orient2d. */
    int e0 = 0;
    for (int i = 1; i < n; i++)
        if (P[i].x < P[e0].x || (P[i].x == P[e0].x && P[i].y < P[e0].y)) e0 = i;
    int poly_or = test_orientation_2d((s_point2d[]){ P[(e0+n-1)%n], P[e0] }, P[(e0+1)%n]);
    if (poly_or == 0) {              /* corner degenerate (e.g. a bridge twin): shoelace */
        double s2 = 0;
        for (int i = 0; i < n; i++) {
            s_point2d u = P[i], v = P[(i+1)%n];
            s2 += u.x * v.y - v.x * u.y;
        }
        poly_or = (s2 > 0) - (s2 < 0);
        if (poly_or == 0) poly_or = 1;                    /* degenerate: fan still closes */
    }

    /* Clip only STRICTLY-convex ears. A flat (collinear) vertex is never clipped
     * as an ear -- doing so would create a chord skipping it, along the straight
     * cell edge, which collides with the identical chord from the neighbouring
     * facet (both facets share that edge). Leaving it in the ring preserves its
     * two boundary edges (matching the neighbour) and connects it via interior
     * diagonals instead. */
    int guard = 0, guard_max = 3 * n + 4;
    while (m > 3 && guard++ < guard_max) {
        /* Drop zero-length ring edges: clipping around a keyhole bridge
         * eventually makes the two copies of a bridge endpoint adjacent. */
        for (int i = 0; i < m && m > 3; ) {
            int a = V[i], b = V[(i+1)%m];
            if (P[a].x == P[b].x && P[a].y == P[b].y) {
                for (int k = (i+1)%m; k < m - 1; k++) V[k] = V[k+1];
                m--; i = 0;                   /* removal can create a new pair */
            } else i++;
        }
        if (m == 3) break;
        int clipped = 0;
        for (int i = 0; i < m; i++) {
            int a = V[(i+m-1)%m], b = V[i], c = V[(i+1)%m];
            int turn = test_orientation_2d((s_point2d[]){ P[a], P[b] }, P[c]);
            if (turn != poly_or) continue;                /* reflex/collinear: not an ear */
            int ok = 1;                                   /* no other vertex inside a,b,c */
            for (int j = 0; j < m && ok; j++) {
                int vj = V[j];
                if (vj == a || vj == b || vj == c) continue;
                if ((P[vj].x == P[a].x && P[vj].y == P[a].y) ||
                    (P[vj].x == P[b].x && P[vj].y == P[b].y) ||
                    (P[vj].x == P[c].x && P[vj].y == P[c].y)) continue;  /* corner's bridge twin */
                int o1 = test_orientation_2d((s_point2d[]){ P[a], P[b] }, P[vj]);
                int o2 = test_orientation_2d((s_point2d[]){ P[b], P[c] }, P[vj]);
                int o3 = test_orientation_2d((s_point2d[]){ P[c], P[a] }, P[vj]);
                if ((o1 == poly_or || o1 == 0) &&
                    (o2 == poly_or || o2 == 0) &&
                    (o3 == poly_or || o3 == 0)) ok = 0;   /* inside/on -> blocks ear */
            }
            if (!ok) continue;
            tris[3*nt + 0] = a; tris[3*nt + 1] = b; tris[3*nt + 2] = c;
            nt++;
            for (int k = i; k < m - 1; k++) V[k] = V[k+1];
            m--; clipped = 1; break;
        }
        if (!clipped) break;                              /* no ear -> fan the rest */
    }
    if (m == 3) {
        tris[3*nt + 0] = V[0]; tris[3*nt + 1] = V[1]; tris[3*nt + 2] = V[2];
        nt++;
    } else {
        for (int i = 0; i < m; i++)          /* fanning across a bridge is unsafe: report */
            for (int j = i + 1; j < m; j++)
                if (P[V[i]].x == P[V[j]].x && P[V[i]].y == P[V[j]].y) {
                    fprintf(stderr, "surf_earclip: fan fallback on a bridged ring "
                                    "(facet with hole; report)\n");
                    j = m; i = m;
                }
        for (int i = 1; i + 1 < m; i++) {                 /* fan fallback (Steiner-free) */
            tris[3*nt + 0] = V[0]; tris[3*nt + 1] = V[i]; tris[3*nt + 2] = V[i+1];
            nt++;
        }
    }

    /* Shape pass: same vertices, same boundary, better diagonals. */
    sv_flip_to_delaunay(P, tris, nt);

    for (int t = 0; t < nt; t++) {
        int tri[3] = { loop[tris[3*t]], loop[tris[3*t + 1]], loop[tris[3*t + 2]] };
        if (tri[0] == tri[1] || tri[1] == tri[2] || tri[2] == tri[0]) continue; /* bridge-twin sliver */
        if (!dynarray_push(faces_out, tri)) { rc = -1; break; }
    }

    free(P); free(V); free(tris);
    return rc;
}

/* ---- polygon-with-holes assembly (annular cut facets) ----------------------
 * A facet pierced by a tunnel of the cell traces into more than one loop on
 * its plane: the outer boundary plus one loop per hole, holes wound OPPOSITE
 * to the outer (the XOR hands every loop with the facet region on the same
 * side).  Ear-clipping each loop independently double-covers every hole (a
 * full outer disk plus a reversed hole disk: two coincident opposite sheets
 * that break any downstream CDT), so nested loops are first merged into ONE
 * weakly-simple ring per outer loop with keyhole bridges:
 *   ring = ... w, h, hole..., h, w, ...
 * (w outer / h hole vertex, mutually visible, each appearing twice; the hole
 * is spliced in its stored winding, which is already the reverse of the
 * outer's).  surf_earclip handles the duplicated bridge vertices. */

/* w strictly inside segment (u,v); all three collinear (checked by caller). */
static inline int sv_between2d(s_point2d u, s_point2d v, s_point2d w)
{
    if (fabs(v.x - u.x) >= fabs(v.y - u.y))
        return (u.x < w.x && w.x < v.x) || (v.x < w.x && w.x < u.x);
    return (u.y < w.y && w.y < v.y) || (v.y < w.y && w.y < u.y);
}

/* Does segment (p,q) cross edge (a,b) properly, or touch its interior / have
 * an endpoint inside it?  Exact orient2d; endpoint-to-endpoint contact is NOT
 * a hit (the caller skips edges incident to the bridge endpoints by id). */
static inline int sv_seg_hits_edge(s_point2d p, s_point2d q, s_point2d a, s_point2d b)
{
    int o1 = test_orientation_2d((s_point2d[]){ p, q }, a);
    int o2 = test_orientation_2d((s_point2d[]){ p, q }, b);
    int o3 = test_orientation_2d((s_point2d[]){ a, b }, p);
    int o4 = test_orientation_2d((s_point2d[]){ a, b }, q);
    if (o1*o2 < 0 && o3*o4 < 0) return 1;             /* proper crossing */
    if (o1 == 0 && sv_between2d(p, q, a)) return 1;   /* grazing contacts */
    if (o2 == 0 && sv_between2d(p, q, b)) return 1;
    if (o3 == 0 && sv_between2d(a, b, p)) return 1;
    if (o4 == 0 && sv_between2d(a, b, q)) return 1;
    return 0;
}

/* Signed double-area of a loop in the drop-projection (float shoelace: only
 * the SIGN is consumed, and real facet loops have macroscopic area). */
static inline double sv_loop_area2(const s_point *coords, const int *loop, int n, int drop)
{
    double s2 = 0;
    s_point2d u = surf_proj2d(coords[loop[n-1]], drop);
    for (int i = 0; i < n; i++) {
        s_point2d v = surf_proj2d(coords[loop[i]], drop);
        s2 += u.x * v.y - v.x * u.y;
        u = v;
    }
    return s2;
}

/* Is p strictly inside the loop?  Crossing parity of the +x ray, with the
 * side test exact (orient2d); p exactly on the boundary counts as outside. */
static inline int sv_point_in_loop(s_point2d p, const s_point *coords,
                                   const int *loop, int n, int drop)
{
    int cross = 0;
    for (int i = 0; i < n; i++) {
        s_point2d a = surf_proj2d(coords[loop[i]], drop);
        s_point2d b = surf_proj2d(coords[loop[(i+1)%n]], drop);
        if ((a.y > p.y) == (b.y > p.y)) continue;        /* no y-straddle */
        int o = test_orientation_2d((s_point2d[]){ a, b }, p);
        if (o == 0) return 0;                            /* on the boundary */
        if ((b.y > a.y) ? (o > 0) : (o < 0)) cross ^= 1; /* hit right of p */
    }
    return cross;
}

/* Candidate keyhole bridge, sorted nearest-first (deterministic tiebreak). */
typedef struct sv_bpair { double d2; int hi, wi; } s_sv_bpair;
static inline int sv_bpair_cmp(const void *A, const void *B)
{
    const s_sv_bpair *a = (const s_sv_bpair *)A, *b = (const s_sv_bpair *)B;
    if (a->d2 != b->d2) return a->d2 < b->d2 ? -1 : 1;
    if (a->hi != b->hi) return a->hi - b->hi;
    return a->wi - b->wi;
}

/* Bridge admissibility: (p,q) = projected (hole vtx hid, ring vtx wid).  It
 * must hit no blocking edge (every loop edge of the plane + earlier bridges,
 * minus edges incident to hid/wid) and pass through no other loop vertex. */
static inline int sv_bridge_ok(int hid, int wid, s_point2d p, s_point2d q,
                               const int *eu, const int *ev, int ne,
                               const int *plb, int total,
                               const s_point *coords, int drop)
{
    for (int e = 0; e < ne; e++) {
        if (eu[e] == hid || ev[e] == hid || eu[e] == wid || ev[e] == wid) continue;
        if (sv_seg_hits_edge(p, q, surf_proj2d(coords[eu[e]], drop),
                                   surf_proj2d(coords[ev[e]], drop))) return 0;
    }
    for (int j = 0; j < total; j++) {                /* no vertex ON the bridge */
        int x = plb[j];
        if (x == hid || x == wid) continue;
        s_point2d X = surf_proj2d(coords[x], drop);
        if ((X.x == p.x && X.y == p.y) || (X.x == q.x && X.y == q.y)) continue;
        if (test_orientation_2d((s_point2d[]){ p, q }, X) == 0 &&
            sv_between2d(p, q, X)) return 0;
    }
    return 1;
}

static inline int sv_id_count(const int *arr, int n, int id)
{
    int c = 0;
    for (int i = 0; i < n; i++) if (arr[i] == id) c++;
    return c;
}

/* Append one final (ear-clippable) loop to the global loop buffers. */
static inline int sv_append_loop(const int *loop, int n, s_dynarray *lbuf_da,
                                 s_dynarray *loff_da, int *nloops)
{
    for (int j = 0; j < n; j++) if (!dynarray_push(lbuf_da, &loop[j])) return -1;
    int off = (int)lbuf_da->N;
    if (!dynarray_push(loff_da, &off)) return -1;
    (*nloops)++;
    return 0;
}

/* Classify one plane's traced loops by nesting and emit ear-clippable rings:
 * loops wound WITH the net region orientation are outers, opposite ones are
 * holes; each hole is keyhole-merged into its smallest containing outer.  A
 * hole with no containing outer or no admissible bridge is emitted verbatim
 * (the pre-fix behaviour: it double-covers) with a warning.  plb/plo hold
 * the plane's loops (plo[0] = 0, loop i = plb[plo[i]..plo[i+1])), pn >= 1.
 * Returns 0 OK, -1 alloc error. */
static inline int sv_emit_plane_loops(const s_point *coords, int drop,
                                      const int *plb, const int *plo, int pn,
                                      s_dynarray *lbuf_da, s_dynarray *loff_da,
                                      int *nloops, int cell_seed)
{
    if (pn <= 0) return 0;
    if (pn == 1)                                     /* the common case */
        return sv_append_loop(plb, plo[1], lbuf_da, loff_da, nloops);

    int total = plo[pn];
    double *a2 = (double *)malloc(sizeof(double) * pn);
    int *parent = (int *)malloc(sizeof(int) * pn);
    if (!a2 || !parent) { free(a2); free(parent); return -1; }
    double asum = 0;
    for (int i = 0; i < pn; i++) {
        a2[i] = sv_loop_area2(coords, plb + plo[i], plo[i+1] - plo[i], drop);
        asum += a2[i];
    }
    int S = (asum > 0) - (asum < 0);
    int nholes = 0;
    for (int i = 0; i < pn; i++) {
        parent[i] = -1;
        if (S && a2[i] * S < 0) nholes++;
    }
    if (nholes == 0) {                               /* disjoint facets only */
        int rc = 0;
        for (int i = 0; i < pn && rc == 0; i++)
            rc = sv_append_loop(plb + plo[i], plo[i+1] - plo[i], lbuf_da, loff_da, nloops);
        free(a2); free(parent);
        return rc;
    }

    /* parent of each hole = smallest-|area| outer loop containing it */
    for (int i = 0; i < pn; i++) {
        if (a2[i] * S >= 0) continue;
        s_point2d p = surf_proj2d(coords[plb[plo[i]]], drop);
        for (int o = 0; o < pn; o++) {
            if (o == i || a2[o] * S <= 0) continue;
            if (!sv_point_in_loop(p, coords, plb + plo[o], plo[o+1] - plo[o], drop)) continue;
            if (parent[i] < 0 || fabs(a2[o]) < fabs(a2[parent[i]])) parent[i] = o;
        }
        if (parent[i] < 0)
            fprintf(stderr, "sv_emit_plane_loops: cell %d hole loop without a "
                            "containing outer loop (emitted uncut)\n", cell_seed);
    }

    int *eu = (int *)malloc(sizeof(int) * (total + nholes));   /* blocking edges */
    int *ev = (int *)malloc(sizeof(int) * (total + nholes));
    int *ring = (int *)malloc(sizeof(int) * (total + 2 * nholes));
    if (!eu || !ev || !ring) {
        free(a2); free(parent); free(eu); free(ev); free(ring); return -1;
    }
    int ne = 0;
    for (int i = 0; i < pn; i++)
        for (int j = plo[i]; j < plo[i+1]; j++) {
            eu[ne] = plb[j];
            ev[ne] = plb[(j + 1 < plo[i+1]) ? j + 1 : plo[i]];
            ne++;
        }

    int rc = 0;
    for (int o = 0; o < pn && rc == 0; o++) {        /* holes ride with their outer */
        if (a2[o] * S < 0) continue;
        int rn = plo[o+1] - plo[o];
        memcpy(ring, plb + plo[o], sizeof(int) * rn);
        for (int h = 0; h < pn && rc == 0; h++) {
            if (parent[h] != o) continue;
            const int *hl = plb + plo[h];
            int hn = plo[h+1] - plo[h];
            /* nearest admissible (hole vtx, ring vtx) pair; previously used
             * bridge endpoints (duplicated in ring) and pinch-repeated
             * vertices are skipped so the splice position is unambiguous */
            s_sv_bpair *pair = (s_sv_bpair *)malloc(sizeof(s_sv_bpair) * (size_t)hn * rn);
            if (!pair) { rc = -1; break; }
            int np = 0;
            for (int i = 0; i < hn; i++) {
                if (sv_id_count(hl, hn, hl[i]) > 1) continue;
                s_point2d p = surf_proj2d(coords[hl[i]], drop);
                for (int j = 0; j < rn; j++) {
                    if (ring[j] == hl[i] || sv_id_count(ring, rn, ring[j]) > 1) continue;
                    s_point2d q = surf_proj2d(coords[ring[j]], drop);
                    if (p.x == q.x && p.y == q.y) continue;
                    double dx = q.x - p.x, dy = q.y - p.y;
                    pair[np++] = (s_sv_bpair){ dx*dx + dy*dy, i, j };
                }
            }
            qsort(pair, np, sizeof(s_sv_bpair), sv_bpair_cmp);
            int found = -1;
            for (int k = 0; k < np; k++) {
                s_point2d p = surf_proj2d(coords[hl[pair[k].hi]], drop);
                s_point2d q = surf_proj2d(coords[ring[pair[k].wi]], drop);
                if (sv_bridge_ok(hl[pair[k].hi], ring[pair[k].wi], p, q,
                                 eu, ev, ne, plb, total, coords, drop)) { found = k; break; }
            }
            if (found < 0) {
                fprintf(stderr, "sv_emit_plane_loops: cell %d no admissible bridge "
                                "for a hole loop (emitted uncut)\n", cell_seed);
                rc = sv_append_loop(hl, hn, lbuf_da, loff_da, nloops);
                free(pair);
                continue;
            }
            int hpos = pair[found].hi, iw = pair[found].wi;
            int hid = hl[hpos], wid = ring[iw];
            free(pair);
            /* splice: ... w, h, hole in stored (reverse) winding, h, w, ... */
            memmove(ring + iw + hn + 3, ring + iw + 1, sizeof(int) * (rn - iw - 1));
            for (int k = 0; k < hn; k++) ring[iw + 1 + k] = hl[(hpos + k) % hn];
            ring[iw + 1 + hn] = hid;
            ring[iw + 2 + hn] = wid;
            rn += hn + 2;
            eu[ne] = hid; ev[ne] = wid; ne++;        /* later bridges must clear it */
        }
        if (rc == 0) rc = sv_append_loop(ring, rn, lbuf_da, loff_da, nloops);
    }
    for (int h = 0; h < pn && rc == 0; h++)          /* orphan holes: verbatim */
        if (a2[h] * S < 0 && parent[h] < 0)
            rc = sv_append_loop(plb + plo[h], plo[h+1] - plo[h], lbuf_da, loff_da, nloops);

    free(a2); free(parent); free(eu); free(ev); free(ring);
    return rc;
}

/* Undirected edge -> incident-triangle bucket, for splitting the assembled
 * surface into connected components at non-manifold (pinch) edges. */
typedef struct e2tkey { int32_t lo, hi; } s_e2tkey;
#define SV_EMAXT 12
typedef struct e2tval { int32_t lo, hi, cnt; int32_t t[SV_EMAXT]; } s_e2tval;

static inline size_t e2t_hash(const void *k)
{
    const unsigned char *p = (const unsigned char *)k;
    size_t h = 1469598103934665603ULL;
    for (size_t i = 0; i < sizeof(s_e2tkey); i++) { h ^= p[i]; h *= 1099511628211ULL; }
    return h;
}
static inline bool e2t_equals(const void *a, const void *b)
{
    return memcmp(a, b, sizeof(s_e2tkey)) == 0;
}

/* Edge XOR bucket: key is (plane, lo, hi); value carries the same fields plus a
 * signed net count (+1 for lo->hi, -1 for hi->lo). Survivors have net != 0. */
typedef struct edgekey { s_pkey plane; int32_t lo, hi; } s_edgekey;   /* 24B, no pad */
typedef struct edgerec { int32_t net; s_pkey plane; int32_t lo, hi; } s_edgerec;

static inline size_t edgekey_hash(const void *k)
{
    const unsigned char *p = (const unsigned char *)k;
    size_t h = 1469598103934665603ULL;
    for (size_t i = 0; i < sizeof(s_edgekey); i++) { h ^= p[i]; h *= 1099511628211ULL; }
    return h;
}
static inline bool edgekey_equals(const void *a, const void *b)
{
    return memcmp(a, b, sizeof(s_edgekey)) == 0;
}
static inline int edgerec_cmp_plane(const void *A, const void *B)
{
    const s_edgerec *a = (const s_edgerec *)A, *b = (const s_edgerec *)B;
    return memcmp(&a->plane, &b->plane, sizeof(s_pkey));
}

/* ---- combinatorial facet tracker helpers ------------------------------- */

static inline size_t pkey_hash(const void *k)
{
    const unsigned char *p = (const unsigned char *)k;
    size_t h = 1469598103934665603ULL;
    for (size_t i = 0; i < sizeof(s_pkey); i++) { h ^= p[i]; h *= 1099511628211ULL; }
    return h;
}
static inline bool pkey_equals(const void *a, const void *b)
{
    return memcmp(a, b, sizeof(s_pkey)) == 0;
}

static inline int vk_in3(int x, int a, int b, int c) { return x == a || x == b || x == c; }
static inline int vk_in4(int x, int a, int b, int c, int d) { return x==a||x==b||x==c||x==d; }

/* Does the surface vertex `vk` lie on facet plane P?  Purely combinatorial: the
 * vkey encodes the constraint planes meeting at the vertex; match them against
 * P's identity (seed pair for a bisector, CDT-vertex triple for a face). This
 * membership is COMPLETE -- it catches a vertex sitting on a straight facet edge
 * even when that facet's own triangulation skipped it (the T-junction case). */
static inline int on_plane(s_vkey vk, s_pkey P)
{
    if (P.tag == PK_BIS) {                      /* bisector between seeds {P.a,P.b} */
        switch (vk.tag) {
            case VK_VB: return vk.c == P.a && vk.d == P.b;              /* its bisector {c,d} */
            case VK_VV: return vk_in4(P.a, vk.a,vk.b,vk.c,vk.d)
                            && vk_in4(P.b, vk.a,vk.b,vk.c,vk.d);        /* equidistant to both */
            case VK_EF:
            case VK_EE: return vk_in3(P.a, vk.a,vk.b,vk.c)
                            && vk_in3(P.b, vk.a,vk.b,vk.c);             /* Voronoi-edge seeds */
            default:    return 0;                                      /* VK_TV */
        }
    }
    switch (vk.tag) {                           /* PK_FACE: boundary face {P.a,P.b,P.c} */
        case VK_TV: return vk.a == P.a || vk.a == P.b || vk.a == P.c;  /* a face vertex */
        case VK_VB: return vk_in3(vk.a, P.a,P.b,P.c)
                        && vk_in3(vk.b, P.a,P.b,P.c);                   /* an edge of the face */
        case VK_EF: return vk.d == P.a && vk.e == P.b && vk.f == P.c;  /* its face */
        case VK_EE: return vk_in3(vk.d, P.a,P.b,P.c)
                        && vk_in3(vk.e, P.a,P.b,P.c);                   /* on an edge of the face */
        default:    return 0;                                         /* VK_VV */
    }
}

/* Complete plane membership for a (possibly degenerate, merged) vertex: the
 * canonical key alone under-reports for a coincident cluster, so also consult
 * the alternate merged keys and, for a bisector plane of this cell, the
 * unioned partner list. */
static inline int lb_on_plane(const s_vlabel *L, s_pkey P, int cell_seed)
{
    if (on_plane(L->key, P)) return 1;
    for (int i = 0; i < L->nalt; i++) if (on_plane(L->alt[i], P)) return 1;
    if (P.tag == PK_BIS && (P.a == cell_seed || P.b == cell_seed)) {
        int t = (P.a == cell_seed) ? P.b : P.a;
        for (int k = 0; k < SV_NBIS; k++) if (L->bis[k] == t) return 1;
    }
    return 0;
}

/* ---- exact edge ordering (Phase B) ------------------------------------------
 * ord_line(a,b, ...) returns sign(pos_a - pos_b) along the facet-edge line P^Q,
 * = -sign(Delta_b) * feas(P,Q,R_a; R_b)  (lp3_predicates.tex Thm 8.1).  feas is
 * an exact lp3_* predicate; the slope sign Delta_b = det[n_Rb; n_P; n_Q] is a
 * plain float determinant (nonzero and well-separated for any real vertex, so its
 * sign is reliable -- only feas needs exactness, to return 0 on coincidence). */

static inline s_point sv_plane_nrm(s_pkey P, const s_points *sd, const s_points *cv)
{
    if (P.tag == PK_BIS) return subtract_points(sd->p[P.b], sd->p[P.a]);
    return cross_prod(subtract_points(cv->p[P.b], cv->p[P.a]),
                      subtract_points(cv->p[P.c], cv->p[P.a]));
}
static inline int sv_det3s(s_point a, s_point b, s_point c)
{
    double d = a.x*(b.y*c.z-b.z*c.y) - a.y*(b.x*c.z-b.z*c.x) + a.z*(b.x*c.y-b.y*c.x);
    return (d > 0) - (d < 0);
}
/* the OTHER boundary face carrying CDT edge {e0,e1}, besides facet notF (-1 none). */
static inline int sv_other_face(int e0, int e1, s_pkey notF, const s_pkey *fac, int nf)
{
    for (int f = 0; f < nf; f++) {
        if (fac[f].tag != PK_FACE || pkey_equals(&fac[f], &notF)) continue;
        if (vk_in3(e0, fac[f].a, fac[f].b, fac[f].c) && vk_in3(e1, fac[f].a, fac[f].b, fac[f].c))
            return f;
    }
    return -1;
}
/* partner seeds of X other than the cell seed cs (into out[]); returns count. */
static inline int sv_partners(s_vkey X, int cs, int out[3])
{
    int n = 0, sd[4], k = 0;
    if (X.tag == VK_VV)      { sd[0]=X.a; sd[1]=X.b; sd[2]=X.c; sd[3]=X.d; k=4; }
    else if (X.tag == VK_EF) { sd[0]=X.a; sd[1]=X.b; sd[2]=X.c; k=3; }
    else if (X.tag == VK_VB) { sd[0]=X.c; sd[1]=X.d; k=2; }
    for (int i = 0; i < k; i++) if (sd[i] != cs) out[n++] = sd[i];
    return n;
}

/* Exact position sign along the line P^Q; MAY return 0 (coincident on the line).
 * The 0 result is the exact coincidence certificate used by the merge (step 1.6). */
static inline int ord_line_raw(s_vkey a, s_vkey b, s_point ca, s_point cb,
                               s_pkey P, s_pkey Q, int cs,
                               const s_points *sd, const s_points *cv,
                               const s_pkey *facet, int nfac)
{
    s_point s  = sd->p[cs];
    s_point nP = sv_plane_nrm(P, sd, cv), nQ = sv_plane_nrm(Q, sd, cv);
    int feas = 2, db = 0;                        /* feas sentinel 2 = unhandled */

    if (P.tag == PK_BIS && Q.tag == PK_BIS) {    /* L1: Voronoi edge; verts VV/EF */
        int t1 = (P.a==cs)?P.b:P.a, t2 = (Q.a==cs)?Q.b:Q.a;
        s_point T1 = sd->p[t1], T2 = sd->p[t2];
        if (b.tag == VK_VV) {                    /* R_b = bisector B(cs,t3b) */
            int pb[3], nb = sv_partners(b, cs, pb), t3b = -1;
            for (int i=0;i<nb;i++) if (pb[i]!=t1 && pb[i]!=t2) t3b = pb[i];
            if (t3b < 0) goto fallback;
            db = sv_det3s(subtract_points(sd->p[t3b], s), nP, nQ);
            s_point u = sd->p[t3b];
            if (a.tag == VK_VV) {                /* SSS vs S = insphere */
                int pa[3], na = sv_partners(a, cs, pa), t3a = -1;
                for (int i=0;i<na;i++) if (pa[i]!=t1 && pa[i]!=t2) t3a = pa[i];
                if (t3a < 0) goto fallback;
                s_point sph[4] = { s, T1, T2, sd->p[t3a] };
                feas = -test_insphere(sph, u);
            } else if (a.tag == VK_EF) {         /* TSS vs S */
                s_point A=cv->p[a.d], Qf=cv->p[a.e], R=cv->p[a.f];
                feas = lp3_feasible_TSS_S(A.x,A.y,A.z, Qf.x,Qf.y,Qf.z, R.x,R.y,R.z,
                                          s.x,s.y,s.z, T1.x,T1.y,T1.z, T2.x,T2.y,T2.z,
                                          u.x,u.y,u.z);
            }
        } else if (b.tag == VK_EF) {             /* R_b = b's face */
            s_point A2=cv->p[b.d], Q2=cv->p[b.e], R2=cv->p[b.f];
            db = sv_det3s(cross_prod(subtract_points(Q2,A2),subtract_points(R2,A2)), nP, nQ);
            if (a.tag == VK_VV) {                /* SSS vs T */
                int pa[3], na = sv_partners(a, cs, pa), t3a = -1;
                for (int i=0;i<na;i++) if (pa[i]!=t1 && pa[i]!=t2) t3a = pa[i];
                if (t3a < 0) goto fallback;
                s_point U = sd->p[t3a];
                feas = lp3_feasible_SSS_T(s.x,s.y,s.z, T1.x,T1.y,T1.z, T2.x,T2.y,T2.z,
                                          U.x,U.y,U.z, A2.x,A2.y,A2.z, Q2.x,Q2.y,Q2.z, R2.x,R2.y,R2.z);
            } else if (a.tag == VK_EF) {         /* TSS vs T (general query face) */
                s_point A1=cv->p[a.d], Q1=cv->p[a.e], R1=cv->p[a.f];
                feas = lp3_feasible_TSS_T_gen(A1.x,A1.y,A1.z, Q1.x,Q1.y,Q1.z, R1.x,R1.y,R1.z,
                                              A2.x,A2.y,A2.z, Q2.x,Q2.y,Q2.z, R2.x,R2.y,R2.z,
                                              s.x,s.y,s.z, T1.x,T1.y,T1.z, T2.x,T2.y,T2.z);
            }
        }
    }
    else if (P.tag != Q.tag) {                   /* L2: bisector ^ face; verts EF/VB */
        s_pkey BP = (P.tag==PK_BIS)?P:Q, FP = (P.tag==PK_FACE)?P:Q;
        int t = (BP.a==cs)?BP.b:BP.a;
        s_point T = sd->p[t];
        if (b.tag == VK_EF) {                    /* R_b = other bisector B(cs,t') */
            int pb[3], nb = sv_partners(b, cs, pb), tp = -1;
            for (int i=0;i<nb;i++) if (pb[i]!=t) tp = pb[i];
            if (tp < 0) goto fallback;
            db = sv_det3s(subtract_points(sd->p[tp], s), nP, nQ);
            s_point u = sd->p[tp];
            if (a.tag == VK_EF) {                /* TSS vs S */
                int pa[3]; sv_partners(a, cs, pa); int ta = (pa[0]!=t)?pa[0]:pa[1];
                s_point A=cv->p[a.d],Qf=cv->p[a.e],R=cv->p[a.f], Ta=sd->p[ta];
                feas = lp3_feasible_TSS_S(A.x,A.y,A.z,Qf.x,Qf.y,Qf.z,R.x,R.y,R.z,
                                          s.x,s.y,s.z, T.x,T.y,T.z, Ta.x,Ta.y,Ta.z, u.x,u.y,u.z);
            } else if (a.tag == VK_VB) {         /* TTS vs S */
                s_point A=cv->p[a.a], B=cv->p[a.b];
                feas = lp3_feasible_TTS_S(A.x,A.y,A.z,B.x,B.y,B.z, s.x,s.y,s.z,
                                          T.x,T.y,T.z, u.x,u.y,u.z);
            }
        } else if (b.tag == VK_VB) {             /* R_b = other face F2 on b's edge */
            int f2 = sv_other_face(b.a, b.b, FP, facet, nfac);
            if (f2 < 0) goto fallback;
            s_pkey F2 = facet[f2];
            s_point A2=cv->p[F2.a],Q2=cv->p[F2.b],R2=cv->p[F2.c];
            db = sv_det3s(cross_prod(subtract_points(Q2,A2),subtract_points(R2,A2)), nP, nQ);
            if (a.tag == VK_EF) {                /* TSS vs T (gen) */
                int pa[3]; sv_partners(a, cs, pa); int ta = (pa[0]!=t)?pa[0]:pa[1];
                s_point A1=cv->p[a.d],Q1=cv->p[a.e],R1=cv->p[a.f], Ta=sd->p[ta];
                feas = lp3_feasible_TSS_T_gen(A1.x,A1.y,A1.z,Q1.x,Q1.y,Q1.z,R1.x,R1.y,R1.z,
                                              A2.x,A2.y,A2.z,Q2.x,Q2.y,Q2.z,R2.x,R2.y,R2.z,
                                              s.x,s.y,s.z, T.x,T.y,T.z, Ta.x,Ta.y,Ta.z);
            } else if (a.tag == VK_VB) {         /* TTS vs T */
                s_point A=cv->p[a.a],B=cv->p[a.b];
                feas = lp3_feasible_TTS_T(A.x,A.y,A.z,B.x,B.y,B.z, s.x,s.y,s.z, T.x,T.y,T.z,
                                          A2.x,A2.y,A2.z,Q2.x,Q2.y,Q2.z,R2.x,R2.y,R2.z);
            }
        }
    }
    else {                                       /* L3: tet edge (2 faces); verts VB */
        if (b.tag == VK_VB && a.tag == VK_VB) {  /* R_b = bisector B(cs,tb); TTS vs S */
            int pb[3]; sv_partners(b, cs, pb); int tb = pb[0];
            int pa[3]; sv_partners(a, cs, pa); int ta = pa[0];
            db = sv_det3s(subtract_points(sd->p[tb], s), nP, nQ);
            s_point A=cv->p[a.a], B=cv->p[a.b], Ta=sd->p[ta], u=sd->p[tb];
            feas = lp3_feasible_TTS_S(A.x,A.y,A.z,B.x,B.y,B.z, s.x,s.y,s.z,
                                      Ta.x,Ta.y,Ta.z, u.x,u.y,u.z);
        }
    }

    if (feas != 2 && db != 0) return -db * feas;   /* exact position order (may be 0) */

fallback:                                          /* unhandled / degenerate slope */
    {
        s_point d = cross_prod(nP, nQ);
        double v = dot_prod(d, subtract_points(ca, cb));
        return (v>0)-(v<0);
    }
}

/* Total order for sorting: exact position, then the symbolic-perturbation tiebreak
 * (Phase C) for coincident vertices. */
static inline int ord_line(s_vkey a, s_vkey b, s_point ca, s_point cb,
                           s_pkey P, s_pkey Q, int cs,
                           const s_points *sd, const s_points *cv,
                           const s_pkey *facet, int nfac)
{
    int r = ord_line_raw(a, b, ca, cb, P, Q, cs, sd, cv, facet, nfac);
    return r ? r : sv_vkey_cmp(a, b);
}

/* `tris` holds one s_surf_tri per piece face of the cell. Assembles the boundary
 * surface and splits it into connected components (a pinched cell yields >1).
 * On success fills *out_meshes (malloc'd array, caller frees each + the array)
 * and *out_n (>= 1). Returns 1 (built), 0 (empty / nothing valid), -1 (alloc). */
static inline int build_cell_surface(const s_dynarray *tris, int cell_seed,
                                     const s_points *seeds, const s_points *cverts,
                                     double EPS_DEG, s_trimesh **out_meshes, int *out_n)
{
    *out_meshes = NULL; *out_n = 0;
    size_t Ntri = tris->N;
    if (Ntri == 0) return 0;

    /* 1. Weld corners into vertices.  Identity is COMBINATORIAL: two corners
     *    name the same vertex iff their key sets (canonical key + alternates
     *    accumulated by the piece-level merges) INTERSECT, computed as
     *    connected components (union-find) over keys.  The plain same-key weld
     *    is not enough: a piece hull collapses a sub-EPS micro-feature (e.g. a
     *    Voronoi vertex a few ulps off a CDT face together with its incident
     *    edge/face crossings) into one vertex, but different pieces may keep a
     *    DIFFERENT surviving construction as the representative; the shared-key
     *    union makes them one vertex here no matter which one survived. */
    size_t maxk = 3 * Ntri * (size_t)(1 + SV_NALT);
    s_hash_table vh;                                  /* key -> 1-based key id */
    if (!hash_init(&vh, sizeof(s_vkey), sizeof(int), 2 * maxk + 1, maxk,
                   vkey_hash, vkey_equals, NULL)) return -1;
    s_point  *coords = (s_point *)malloc(sizeof(s_point) * 3 * Ntri);
    s_vlabel *vlab   = (s_vlabel *)malloc(sizeof(s_vlabel) * 3 * Ntri);
    int      *tid    = (int *)malloc(sizeof(int) * 3 * Ntri);   /* per-corner vertex id */
    int      *kpar   = (int *)malloc(sizeof(int) * maxk);       /* key-id union-find */
    int      *kvid   = (int *)malloc(sizeof(int) * maxk);       /* key-root -> vertex id */
    if (!coords || !vlab || !tid || !kpar || !kvid) {
        free(coords); free(vlab); free(tid); free(kpar); free(kvid); hash_free(&vh); return -1;
    }
    int Nv = 0, err = 0, nk = 0;
    for (size_t t = 0; t < Ntri && !err; t++) {       /* pass 1: union keys per corner */
        const s_surf_tri *tr = (const s_surf_tri *)dynarray_get_ptr_c(tris, t);
        for (int k = 0; k < 3 && !err; k++) {
            const s_vlabel *L = &tr->v[k].lb;
            int r0 = -1;
            for (int q = -1; q < (int)L->nalt && !err; q++) {
                const s_vkey *kk = (q < 0) ? &L->key : &L->alt[q];
                int *slot = (int *)hash_get_or_create(&vh, kk);
                if (!slot) { err = 1; break; }
                if (*slot == 0) { *slot = ++nk; kpar[nk - 1] = nk - 1; }
                int r = *slot - 1;
                while (kpar[r] != r) r = kpar[r];     /* find */
                if (r0 < 0) r0 = r;
                else if (r != r0) { if (r < r0) { kpar[r0] = r; r0 = r; } else kpar[r] = r0; }
            }
        }
    }
    for (int i = 0; i < nk; i++) kvid[i] = -1;
    for (size_t t = 0; t < Ntri && !err; t++) {       /* pass 2: one vertex per key root */
        const s_surf_tri *tr = (const s_surf_tri *)dynarray_get_ptr_c(tris, t);
        for (int k = 0; k < 3 && !err; k++) {
            int *slot = (int *)hash_get_or_create(&vh, &tr->v[k].lb.key);
            if (!slot || *slot == 0) { err = 1; break; }
            int r = *slot - 1;
            while (kpar[r] != r) r = kpar[r];
            if (kvid[r] < 0) { kvid[r] = Nv; coords[Nv] = tr->v[k].p; vlab[Nv] = tr->v[k].lb; Nv++; }
            else vlabel_merge(&vlab[kvid[r]], &tr->v[k].lb);
            tid[3*t + k] = kvid[r];
        }
    }
    hash_free(&vh);
    free(kpar); free(kvid);
    if (err || Nv == 0) { free(coords); free(vlab); free(tid); return err ? -1 : 0; }

    /* 1.5 TOLERANCE VERTEX MERGE.  The exact re-keying upstream (canonical_vv_label,
     *    canonical_edge_cluster) already collapses every EXACT coincidence (cospherical
     *    Voronoi vertices, cocircular Voronoi-edge crossings) by giving them one key, so
     *    they welded in step 1.  What remains are GENUINE near-coincidences: distinct
     *    constructions of the same cell vertex that no exact predicate can equate -- e.g.
     *    a Voronoi edge exiting a convex domain through two NEAR-COPLANAR hull faces gives
     *    two crossings ~1e-15 apart (different third plane -> exact-distinct), which are
     *    really one vertex of the convex cell.  Left split they form a zero-area spur that
     *    breaks the facet's ear-clip and spawns spurious zero-volume "double-wall" sheets.
     *    Merge vertex groups within tau = tau_rel * bbox_diag.  Because the exact cases are
     *    gone, the only pairs inside tau are these ~1e-15 near-misses (real features are
     *    macro-apart), so any tau in a wide band gives the identical result -- validated by
     *    the tau sweep.  Representative: smallest key (deterministic, order-independent);
     *    labels are unioned so the merged vertex is on_plane for every member's facets.
     *    SV_TAU_REL overrides tau_rel (0 disables the merge; the subdivision stays exact). */
    int *rep = (int *)malloc(sizeof(int) * Nv);
    if (!rep) { free(coords); free(vlab); free(tid); return -1; }
    double tau2;
    {
        double tau_rel = 1e-9;
        if (getenv("SV_TAU_REL")) tau_rel = atof(getenv("SV_TAU_REL"));
        s_point lo = coords[0], hi = coords[0];
        for (int i = 1; i < Nv; i++) {
            if (coords[i].x < lo.x) lo.x = coords[i].x;  if (coords[i].x > hi.x) hi.x = coords[i].x;
            if (coords[i].y < lo.y) lo.y = coords[i].y;  if (coords[i].y > hi.y) hi.y = coords[i].y;
            if (coords[i].z < lo.z) lo.z = coords[i].z;  if (coords[i].z > hi.z) hi.z = coords[i].z;
        }
        tau2 = norm_squared(subtract_points(hi, lo)) * tau_rel * tau_rel;   /* (tau_rel*diag)^2 */
        for (int i = 0; i < Nv; i++) rep[i] = i;
        if (tau_rel <= 0) goto weld_done;                /* SV_TAU_REL=0 disables the merge */
        for (int i = 0; i < Nv; i++) {                   /* tiny union-find */
            for (int j = i + 1; j < Nv; j++) {
                if (norm_squared(subtract_points(coords[i], coords[j])) >= tau2) continue;
                int ri = i, rj = j;
                while (rep[ri] != ri) ri = rep[ri];
                while (rep[rj] != rj) rj = rep[rj];
                if (ri != rj) rep[rj > ri ? rj : ri] = rj > ri ? ri : rj;
            }
        }
        for (int i = 0; i < Nv; i++) {                   /* path-compress */
            int r = i;
            while (rep[r] != r) r = rep[r];
            rep[i] = r;
        }
        for (int r = 0; r < Nv; r++) {                   /* canonical = min key in group */
            if (rep[r] != r) continue;
            int best = r;
            for (int i = r + 1; i < Nv; i++)
                if (rep[i] == r && sv_vkey_cmp(vlab[i].key, vlab[best].key) < 0) best = i;
            for (int i = r; i < Nv; i++) {
                if (rep[i] != r || i == best) continue;
                vlabel_merge(&vlab[best], &vlab[i]);
            }
            if (best != r)
                for (int i = 0; i < Nv; i++) if (rep[i] == r) rep[i] = best;
        }
        for (size_t c = 0; c < 3 * Ntri; c++) tid[c] = rep[tid[c]];
weld_done: ;
    }

    /* 2. Distinct facet planes (piece_to_surface already dropped internal walls,
     *    so every plane here is a real boundary facet: bisector or domain face). */
    s_pkey *facet = (s_pkey *)malloc(sizeof(s_pkey) * Ntri);
    int nfac = 0;
    if (!facet) { free(coords); free(vlab); free(tid); free(rep); return -1; }
    {
        s_hash_table fh;
        if (!hash_init(&fh, sizeof(s_pkey), sizeof(int), 2*Ntri+1, Ntri,
                       pkey_hash, pkey_equals, NULL)) { free(facet); free(coords); free(vlab); free(tid); free(rep); return -1; }
        for (size_t t = 0; t < Ntri && !err; t++) {
            const s_surf_tri *tr = (const s_surf_tri *)dynarray_get_ptr_c(tris, t);
            int *slot = (int *)hash_get_or_create(&fh, &tr->plane);
            if (!slot) { err = 1; break; }
            if (*slot == 0) { *slot = nfac + 1; facet[nfac++] = tr->plane; }
        }
        hash_free(&fh);
        if (err) { free(facet); free(coords); free(vlab); free(tid); free(rep); return -1; }
    }

    /* 3. Per-plane directed-edge XOR, with each triangle edge SUBDIVIDED at the
     *    line-vertices it passes through.  "w lies on segment (u,v)" is decided
     *    on the welded coordinates within the SAME tau as the vertex merge (step
     *    1.5): w is on the line iff dist(w, line uv) < tau, and strictly between
     *    u,v by the dominant-axis coordinate.  Using the merge's tau (not an
     *    exact collinearity test) is essential: near-line vertices left by the
     *    near-hull degeneracy must be treated consistently by both facets sharing
     *    the segment, or a T-junction crack opens.  The XOR supplies coverage:
     *    interior diagonals cancel within a piece and gap edges never appear, so
     *    non-convex facets are correct. */
    s_hash_table eh;
    if (!hash_init(&eh, sizeof(s_edgekey), sizeof(s_edgerec), 2*3*Ntri+1, 3*Ntri,
                   edgekey_hash, edgekey_equals, NULL)) { free(facet); free(coords); free(vlab); free(tid); free(rep); return -1; }
    int *chain = (int *)malloc(sizeof(int) * (Nv + 1));
    if (!chain) { hash_free(&eh); free(facet); free(coords); free(vlab); free(tid); free(rep); return -1; }
    for (size_t t = 0; t < Ntri && !err; t++) {
        const s_surf_tri *tr = (const s_surf_tri *)dynarray_get_ptr_c(tris, t);
        s_pkey P = tr->plane;
        for (int k = 0; k < 3; k++) {
            int u = tid[3*t+k], v = tid[3*t+(k+1)%3];
            if (u == v) continue;
            double dx = fabs(coords[v].x-coords[u].x), dy = fabs(coords[v].y-coords[u].y), dz = fabs(coords[v].z-coords[u].z);
            int ax = (dx>=dy && dx>=dz) ? 0 : (dy>=dz ? 1 : 2);
            double cu = ax==0?coords[u].x:ax==1?coords[u].y:coords[u].z;
            double cv = ax==0?coords[v].x:ax==1?coords[v].y:coords[v].z;
            int nc = 0;
            s_point ev = subtract_points(coords[v], coords[u]);
            double elen2 = norm_squared(ev);
            for (int w = 0; w < Nv; w++) {
                if (w == u || w == v || rep[w] != w) continue;
                s_point cr = cross_prod(ev, subtract_points(coords[w], coords[u]));
                if (norm_squared(cr) >= tau2 * elen2) continue;   /* near-collinear within tau */
                double cw = ax==0?coords[w].x:ax==1?coords[w].y:coords[w].z;
                if ((cu<cw && cw<cv) || (cv<cw && cw<cu)) chain[nc++] = w;   /* strictly between */
            }
            for (int i = 1; i < nc; i++) {               /* sort along u->v (exact for collinear pts) */
                int wi = chain[i], j = i - 1;
                double cwi = ax==0?coords[wi].x:ax==1?coords[wi].y:coords[wi].z;
                while (j >= 0) {
                    double cj = ax==0?coords[chain[j]].x:ax==1?coords[chain[j]].y:coords[chain[j]].z;
                    if ((cu < cv) ? (cj <= cwi) : (cj >= cwi)) break;
                    chain[j+1] = chain[j]; j--;
                }
                chain[j+1] = wi;
            }
            int prev = u;                                /* emit atomic sub-edges into the XOR */
            for (int i = 0; i <= nc; i++) {
                int curv = (i < nc) ? chain[i] : v;
                if (prev != curv) {
                    s_edgekey key = { P, prev<curv?prev:curv, prev<curv?curv:prev };
                    s_edgerec *r = (s_edgerec *)hash_get_or_create(&eh, &key);
                    if (!r) { err = 1; break; }
                    r->plane = P; r->lo = key.lo; r->hi = key.hi;
                    r->net += (prev < curv) ? 1 : -1;
                }
                prev = curv;
            }
        }
    }
    free(chain); free(facet); free(tid); free(rep);
    if (err) { hash_free(&eh); free(vlab); free(coords); return -1; }

    /* 3b. Surviving directed edges (net != 0), grouped by plane. */
    size_t ne = eh.size;
    s_edgerec *edges = (s_edgerec *)malloc(sizeof(s_edgerec) * (ne ? ne : 1));
    if (!edges) { hash_free(&eh); free(coords); return -1; }
    hash_to_array(&eh, edges);
    hash_free(&eh);
    size_t nsurv = 0;
    for (size_t i = 0; i < ne; i++) if (edges[i].net != 0) edges[nsurv++] = edges[i];
    qsort(edges, nsurv, sizeof(s_edgerec), edgerec_cmp_plane);

    free(vlab);

    /* 3c. Trace each plane's surviving directed edges into closed loops.
     * Vertices can have degree > 2 (real pinches; micro-features merged by the
     * weld), where a naive "first edge whose tail matches" walk pairs chains
     * arbitrarily and produces open walks that must be discarded.  The next
     * edge out of a vertex is instead chosen ANGULARLY: rotating CLOCKWISE (in
     * the plane's 2D projection, exact orient2d) from the direction of
     * arrival, take the first available outgoing edge.  In/out degrees are
     * balanced at every vertex (the edges are a XOR of triangle boundaries)
     * and the angular rule cannot cross chains, so every trace closes and
     * every directed edge is consumed exactly once. */
    size_t mtot = 0;
    for (size_t i = 0; i < nsurv; i++) mtot += (size_t)(edges[i].net > 0 ? edges[i].net : -edges[i].net);
    s_dynarray lbuf_da = dynarray_initialize(sizeof(int), 256);
    s_dynarray loff_da = dynarray_initialize(sizeof(int), 64);
    s_dynarray plb_da  = dynarray_initialize(sizeof(int), 256);   /* one plane's raw loops */
    s_dynarray plo_da  = dynarray_initialize(sizeof(int), 16);
    int  *deu   = (int *)malloc(sizeof(int) * (mtot ? mtot : 1));
    int  *dev   = (int *)malloc(sizeof(int) * (mtot ? mtot : 1));
    char *dused = (char *)calloc(mtot ? mtot : 1, 1);
    int  *tmp   = (int *)malloc(sizeof(int) * (mtot + 1));
    if (!lbuf_da.items || !loff_da.items || !plb_da.items || !plo_da.items ||
        !deu || !dev || !dused || !tmp) {
        dynarray_free(&lbuf_da); dynarray_free(&loff_da);
        dynarray_free(&plb_da); dynarray_free(&plo_da);
        free(deu); free(dev); free(dused); free(tmp); free(edges); free(coords); return -1;
    }
    int nloops = 0, zero = 0;
    if (!dynarray_push(&loff_da, &zero)) err = 1;
    size_t md = 0;
    for (size_t b = 0; b < nsurv && !err; ) {
        size_t e = b + 1;
        while (e < nsurv && edgerec_cmp_plane(&edges[b], &edges[e]) == 0) e++;
        size_t m0 = md;                                   /* this plane's directed edges */
        for (size_t i = b; i < e; i++) {
            int n = edges[i].net, cnt = n > 0 ? n : -n;
            for (int k = 0; k < cnt; k++) {
                deu[md] = n > 0 ? edges[i].lo : edges[i].hi;
                dev[md] = n > 0 ? edges[i].hi : edges[i].lo;
                md++;
            }
        }
        s_point nrm = sv_plane_nrm(edges[b].plane, seeds, cverts);
        int drop = coord_with_largest_component_3D(nrm);
        dynarray_clear(&plb_da); dynarray_clear(&plo_da);
        int pl_n = 0;
        if (!dynarray_push(&plo_da, &zero)) err = 1;
        for (size_t st = m0; st < md && !err; st++) {
            if (dused[st]) continue;
            int len = 0, closed = 0;
            size_t cur = st;
            int guard = 0, gmax = (int)(md - m0) + 2;
            while (guard++ < gmax) {
                dused[cur] = 1;
                tmp[len++] = deu[cur];
                int v = dev[cur];
                s_point2d V = surf_proj2d(coords[v], drop);
                s_point2d U = surf_proj2d(coords[deu[cur]], drop);
                size_t best = SIZE_MAX;
                for (size_t q = m0; q < md; q++) {
                    if (deu[q] != v) continue;
                    if (dused[q] && q != st) continue;    /* st re-eligible = closure */
                    if (best == SIZE_MAX) { best = q; continue; }
                    s_point2d A = surf_proj2d(coords[dev[q]], drop);
                    s_point2d B = surf_proj2d(coords[dev[best]], drop);
                    if (sv_cw_before(V, U, A, B)) best = q;
                }
                if (best == SIZE_MAX) {                   /* in/out imbalance: impossible */
                    fprintf(stderr, "build_cell_surface: cell %d open trace (bug)\n", cell_seed);
                    len = 0; closed = 0; break;
                }
                if (best == st) { closed = 1; break; }
                cur = best;
            }
            if (!closed || len < 3) continue;             /* len<3: degenerate sliver */
            for (int j = 0; j < len && !err; j++) if (!dynarray_push(&plb_da, &tmp[j])) err = 1;
            int off = (int)plb_da.N;
            if (!err && !dynarray_push(&plo_da, &off)) err = 1;
            if (!err) pl_n++;
        }
        /* Nested loops (a facet with holes) become keyhole-merged rings here;
         * the single-loop common case passes through untouched. */
        if (!err && pl_n > 0 &&
            sv_emit_plane_loops(coords, drop, (const int *)plb_da.items,
                                (const int *)plo_da.items, pl_n,
                                &lbuf_da, &loff_da, &nloops, cell_seed) < 0) err = 1;
        b = e;
    }
    free(deu); free(dev); free(dused); free(tmp); free(edges);
    dynarray_free(&plb_da); dynarray_free(&plo_da);
    if (err) { dynarray_free(&lbuf_da); dynarray_free(&loff_da); free(coords); return -1; }

    int *lbuf = (int *)lbuf_da.items;
    int *loff = (int *)loff_da.items;

    /* 4. Triangulate each loop by exact ear-clipping. */
    s_dynarray faces = dynarray_initialize(sizeof(int) * 3, 64);
    if (!faces.items) { dynarray_free(&lbuf_da); dynarray_free(&loff_da); free(coords); return -1; }
    for (int L = 0; L < nloops; L++) {
        if (surf_earclip(coords, lbuf + loff[L], loff[L+1] - loff[L], &faces) < 0) {
            dynarray_free(&faces); dynarray_free(&lbuf_da); dynarray_free(&loff_da); free(coords); return -1;
        }
    }
    dynarray_free(&lbuf_da); dynarray_free(&loff_da);

    int Nf = (int)faces.N;
    if (Nf == 0) { dynarray_free(&faces); free(coords); return 0; }
    const int *F = (const int *)faces.items;

    /* 7. Split into connected components. Union triangles across every edge
     *    shared by EXACTLY two triangles (a manifold edge); pinch edges (3+
     *    triangles) are left un-unioned, so the lobes fall into separate
     *    components -- each of which is a clean closed manifold. */
    int *parent = (int *)malloc(sizeof(int) * Nf);
    int *rank   = (int *)malloc(sizeof(int) * Nf);
    int *labels = (int *)malloc(sizeof(int) * Nf);
    if (!parent || !rank || !labels) {
        free(parent); free(rank); free(labels); dynarray_free(&faces); free(coords); return -1;
    }
    UF_initialize(Nf, parent, rank);

    s_hash_table e2h;
    if (!hash_init(&e2h, sizeof(s_e2tkey), sizeof(s_e2tval), 2 * 3 * Nf + 1, 3 * Nf,
                   e2t_hash, e2t_equals, NULL)) {
        free(parent); free(rank); free(labels); dynarray_free(&faces); free(coords); return -1;
    }
    int herr = 0;
    for (int t = 0; t < Nf && !herr; t++) for (int k = 0; k < 3; k++) {
        int a = F[3*t+k], b = F[3*t+(k+1)%3];
        s_e2tkey key = { a < b ? a : b, a < b ? b : a };
        s_e2tval *v = (s_e2tval *)hash_get_or_create(&e2h, &key);
        if (!v) { herr = 1; break; }
        if (v->cnt == 0) { v->lo = key.lo; v->hi = key.hi; }
        if (v->cnt < SV_EMAXT) v->t[v->cnt] = t;
        v->cnt++;
    }
    if (!herr) {
        size_t nk = e2h.size;
        s_e2tval *ev = (s_e2tval *)malloc(sizeof(s_e2tval) * (nk ? nk : 1));
        if (!ev) herr = 1;
        else {
            hash_to_array(&e2h, ev);
            for (size_t i = 0; i < nk; i++)
                if (ev[i].cnt == 2) UF_union(ev[i].t[0], ev[i].t[1], parent, rank);
            free(ev);
        }
    }
    hash_free(&e2h);
    if (herr) { free(parent); free(rank); free(labels); dynarray_free(&faces); free(coords); return -1; }

    int ncomp = UF_count_sets_and_label(Nf, parent, labels);
    free(parent); free(rank);

    /* 8. Build one trimesh per component (compacting its vertices). */
    s_trimesh *meshes = (s_trimesh *)malloc(sizeof(s_trimesh) * ncomp);
    s_point   *cc     = (s_point *)malloc(sizeof(s_point) * Nv);
    int       *cf     = (int *)malloc(sizeof(int) * 3 * Nf);
    int       *l2     = (int *)malloc(sizeof(int) * (Nv ? Nv : 1));
    if (!meshes || !cc || !cf || !l2) {
        free(meshes); free(cc); free(cf); free(l2); free(labels);
        dynarray_free(&faces); free(coords); return -1;
    }
    for (int i = 0; i < Nv; i++) l2[i] = -1;
    int nmesh = 0;
    for (int c = 0; c < ncomp; c++) {
        int nv = 0, nf = 0;
        for (int t = 0; t < Nf; t++) {
            if (labels[t] != c) continue;
            for (int k = 0; k < 3; k++) {
                int g = F[3*t+k];
                if (l2[g] < 0) { l2[g] = nv; cc[nv] = coords[g]; nv++; }
                cf[3*nf+k] = l2[g];
            }
            nf++;
        }
        s_trimesh m = trimesh_from_arrays(cc, nv, cf, nf, EPS_DEG);
        int valid = trimesh_is_valid(&m);
        double cvol = valid ? volume_trimesh(&m) : 0.0;
        int keep = valid, degenerate = 0;
        if (valid) {
            /* Discard degenerate zero-volume "double-wall" sheets: a real lobe
             * has volume ~ size^3, a folded flat sheet ~ machine epsilon.  The
             * near-coplanar hull faces leave such slivers even after the vertex
             * merge; they carry no volume and are not part of the cell surface. */
            s_point clo = cc[0], chi = cc[0];
            for (int i = 1; i < nv; i++) {
                if (cc[i].x<clo.x)clo.x=cc[i].x; if (cc[i].x>chi.x)chi.x=cc[i].x;
                if (cc[i].y<clo.y)clo.y=cc[i].y; if (cc[i].y>chi.y)chi.y=cc[i].y;
                if (cc[i].z<clo.z)clo.z=cc[i].z; if (cc[i].z>chi.z)chi.z=cc[i].z;
            }
            double d2 = norm_squared(subtract_points(chi, clo));
            if (fabs(cvol) < 1e-9 * d2 * sqrt(d2)) { keep = 0; degenerate = 1; }
        }
        if (keep) meshes[nmesh++] = m;                 /* skip bad or degenerate-sheet components */
        else if (degenerate) free_trimesh(&m);
        for (int t = 0; t < Nf; t++) if (labels[t] == c) /* reset only touched l2 entries */
            for (int k = 0; k < 3; k++) l2[F[3*t+k]] = -1;
    }
    free(cc); free(cf); free(l2); free(labels);
    dynarray_free(&faces); free(coords);

    if (nmesh == 0) { free(meshes); return 0; }
    *out_meshes = meshes;
    *out_n = nmesh;
    return 1;
}

#endif
