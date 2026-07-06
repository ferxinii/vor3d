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
#include "gtests.h"     /* exact orient2d (test_orientation_2d) for ear-clipping */

/* ---- key type ---------------------------------------------------------- */

typedef enum { VK_TV = 0, VK_VB, VK_VV, VK_EF } e_vkind;

/* Fixed-size POD, no padding (7 x int32 = 28 bytes): safe to memcmp / hash by
 * bytes.  Unused fields are -1.  Fields are canonicalised (sorted) so equal
 * points always produce byte-identical keys. */
typedef struct vkey {
    int32_t tag;
    int32_t a, b, c, d, e, f;
} s_vkey;

/* Per-hull-vertex label: identity key + plane membership (for face classify).
 * fmask bit i set => vertex lies on local tet face i (0..3).
 * bis[k] = partner seed id of a bisector {s, bis[k]} the vertex lies on (-1 unused). */
typedef struct vlabel {
    s_vkey  key;
    uint8_t fmask;
    int32_t bis[3];
} s_vlabel;

/* One triangle corner emitted to a cell's surface accumulator. */
typedef struct surf_vtx {
    s_vkey  key;
    s_point p;
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

static inline void vlabel_add_bis(s_vlabel *dst, int t)
{
    if (t < 0) return;
    for (int k = 0; k < 3; k++) if (dst->bis[k] == t) return;   /* already present */
    for (int k = 0; k < 3; k++) if (dst->bis[k] < 0) { dst->bis[k] = t; return; }
    /* overflow (>3 distinct partners): ignore; classification only needs a
     * common partner, which is preserved by the first slots. */
}

static inline void vlabel_merge(s_vlabel *dst, const s_vlabel *src)
{
    dst->fmask |= src->fmask;
    for (int k = 0; k < 3; k++) vlabel_add_bis(dst, src->bis[k]);
    /* keys are combinatorial identities: equal points => equal keys, so `dst`'s
     * key is kept as-is. */
}

/* Partner seed id present in all three vertices' bisector lists (-1 if none):
 * the three vertices of a hull face all lie on the bisector {s, that_partner}. */
static inline int vlabel_common_bis(const s_vlabel *a, const s_vlabel *b,
                                    const s_vlabel *c)
{
    for (int i = 0; i < 3; i++) {
        int t = a->bis[i];
        if (t < 0) continue;
        int inb = 0, inc = 0;
        for (int k = 0; k < 3; k++) {
            if (b->bis[k] == t) inb = 1;
            if (c->bis[k] == t) inc = 1;
        }
        if (inb && inc) return t;
    }
    return -1;
}


/* ---- per-cell assembly: boundary-of-union via per-plane edge XOR ---------
 *
 * Every piece face carries its supporting plane (s_pkey). For each plane we XOR
 * the DIRECTED edges of its triangles (an edge and its reverse cancel):
 *   - interior triangulation diagonals cancel within a piece;
 *   - an interior tet face is contributed by both adjacent tets with opposite
 *     orientation -> its whole boundary cancels -> the plane vanishes;
 *   - a boundary tet face or a bisector facet (possibly spanning several tet
 *     pieces, whose seam edges cancel) leaves its outer boundary loop.
 * The surviving directed edges per plane form closed loops (walked by vertex id,
 * purely combinatorial). Each loop is re-triangulated by exact ear-clipping.
 * Connectivity is entirely key/id-based -> no floating point can open a seam;
 * the ear-clip predicates are exact (orient2d on axis-dropped coordinates), so
 * the triangle *geometry* is correct too. */

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

/* Ear-clip a simple loop (vertex ids into `coords`, length n) into `faces_out`
 * (int[3] per triangle). Exact orient2d decisions. Falls back to a fan if no ear
 * is found (still Steiner-free -> still closed). Returns 1 OK, -1 alloc error. */
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

    s_point2d *P = (s_point2d *)malloc(sizeof(s_point2d) * n);
    int       *V = (int *)malloc(sizeof(int) * n);        /* active ring of indices 0..n-1 */
    if (!P || !V) { free(P); free(V); return -1; }
    for (int i = 0; i < n; i++) { P[i] = surf_proj2d(coords[loop[i]], drop); V[i] = i; }
    int m = n, rc = 1;

    /* Polygon orientation from its lexicographically-lowest vertex (always a
     * strictly convex corner): exact via orient2d. */
    int e0 = 0;
    for (int i = 1; i < n; i++)
        if (P[i].x < P[e0].x || (P[i].x == P[e0].x && P[i].y < P[e0].y)) e0 = i;
    int poly_or = test_orientation_2d((s_point2d[]){ P[(e0+n-1)%n], P[e0] }, P[(e0+1)%n]);
    if (poly_or == 0) poly_or = 1;                        /* degenerate: fan still closes */

    /* Clip only STRICTLY-convex ears. A flat (collinear) vertex is never clipped
     * as an ear -- doing so would create a chord skipping it, along the straight
     * cell edge, which collides with the identical chord from the neighbouring
     * facet (both facets share that edge). Leaving it in the ring preserves its
     * two boundary edges (matching the neighbour) and connects it via interior
     * diagonals instead. */
    int guard = 0, guard_max = 3 * n + 4;
    while (m > 3 && guard++ < guard_max) {
        int clipped = 0;
        for (int i = 0; i < m; i++) {
            int a = V[(i+m-1)%m], b = V[i], c = V[(i+1)%m];
            int turn = test_orientation_2d((s_point2d[]){ P[a], P[b] }, P[c]);
            if (turn != poly_or) continue;                /* reflex/collinear: not an ear */
            int ok = 1;                                   /* no other vertex inside a,b,c */
            for (int j = 0; j < m && ok; j++) {
                int vj = V[j];
                if (vj == a || vj == b || vj == c) continue;
                int o1 = test_orientation_2d((s_point2d[]){ P[a], P[b] }, P[vj]);
                int o2 = test_orientation_2d((s_point2d[]){ P[b], P[c] }, P[vj]);
                int o3 = test_orientation_2d((s_point2d[]){ P[c], P[a] }, P[vj]);
                if ((o1 == poly_or || o1 == 0) &&
                    (o2 == poly_or || o2 == 0) &&
                    (o3 == poly_or || o3 == 0)) ok = 0;   /* inside/on -> blocks ear */
            }
            if (!ok) continue;
            int tri[3] = { loop[a], loop[b], loop[c] };
            if (!dynarray_push(faces_out, tri)) { rc = -1; goto done; }
            for (int k = i; k < m - 1; k++) V[k] = V[k+1];
            m--; clipped = 1; break;
        }
        if (!clipped) break;                              /* no ear -> fan the rest */
    }
    if (m == 3) {
        int tri[3] = { loop[V[0]], loop[V[1]], loop[V[2]] };
        if (!dynarray_push(faces_out, tri)) { rc = -1; goto done; }
    } else {
        for (int i = 1; i + 1 < m; i++) {                 /* fan fallback (Steiner-free) */
            int tri[3] = { loop[V[0]], loop[V[i]], loop[V[i+1]] };
            if (!dynarray_push(faces_out, tri)) { rc = -1; goto done; }
        }
    }
done:
    free(P); free(V);
    return rc;
}

/* Undirected edge -> incident-triangle bucket, for splitting the assembled
 * surface into connected components at non-manifold (pinch) edges. */
typedef struct e2tkey { int32_t lo, hi; } s_e2tkey;
typedef struct e2tval { int32_t cnt, t0, t1; } s_e2tval;

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

/* `tris` holds one s_surf_tri per piece face of the cell. Assembles the boundary
 * surface and splits it into connected components (a pinched cell yields >1).
 * On success fills *out_meshes (malloc'd array, caller frees each + the array)
 * and *out_n (>= 1). Returns 1 (built), 0 (empty / nothing valid), -1 (alloc). */
static inline int build_cell_surface(const s_dynarray *tris, double EPS_DEG,
                                     s_trimesh **out_meshes, int *out_n)
{
    *out_meshes = NULL; *out_n = 0;
    size_t Ntri = tris->N;
    if (Ntri == 0) return 0;

    /* 1. Weld corners by combinatorial key -> vertex ids + coords. */
    s_hash_table vh;
    if (!hash_init(&vh, sizeof(s_vkey), sizeof(int), 2 * 3 * Ntri + 1, 3 * Ntri,
                   vkey_hash, vkey_equals, NULL)) return -1;
    s_point *coords = (s_point *)malloc(sizeof(s_point) * 3 * Ntri);
    int     *tid    = (int *)malloc(sizeof(int) * 3 * Ntri);   /* per-corner vertex id */
    if (!coords || !tid) { free(coords); free(tid); hash_free(&vh); return -1; }
    int Nv = 0, err = 0;
    for (size_t t = 0; t < Ntri && !err; t++) {
        const s_surf_tri *tr = (const s_surf_tri *)dynarray_get_ptr_c(tris, t);
        for (int k = 0; k < 3; k++) {
            int *slot = (int *)hash_get_or_create(&vh, &tr->v[k].key);
            if (!slot) { err = 1; break; }
            if (*slot == 0) { *slot = Nv + 1; coords[Nv] = tr->v[k].p; Nv++; }
            tid[3*t + k] = *slot - 1;
        }
    }
    hash_free(&vh);
    if (err) { free(coords); free(tid); return -1; }

    /* 2. Per-plane directed-edge XOR. */
    s_hash_table eh;
    if (!hash_init(&eh, sizeof(s_edgekey), sizeof(s_edgerec), 2 * 3 * Ntri + 1,
                   3 * Ntri, edgekey_hash, edgekey_equals, NULL)) {
        free(coords); free(tid); return -1;
    }
    for (size_t t = 0; t < Ntri && !err; t++) {
        const s_surf_tri *tr = (const s_surf_tri *)dynarray_get_ptr_c(tris, t);
        for (int k = 0; k < 3; k++) {
            int u = tid[3*t + k], v = tid[3*t + (k+1)%3];
            if (u == v) continue;                         /* degenerate edge */
            s_edgekey key = { tr->plane, u < v ? u : v, u < v ? v : u };
            s_edgerec *r = (s_edgerec *)hash_get_or_create(&eh, &key);
            if (!r) { err = 1; break; }
            r->plane = tr->plane; r->lo = key.lo; r->hi = key.hi;
            r->net += (u < v) ? 1 : -1;
        }
    }
    free(tid);
    if (err) { hash_free(&eh); free(coords); return -1; }

    /* 3. Collect surviving directed edges (net != 0), grouped by plane. */
    size_t ne = eh.size;
    s_edgerec *edges = (s_edgerec *)malloc(sizeof(s_edgerec) * (ne ? ne : 1));
    if (!edges) { hash_free(&eh); free(coords); return -1; }
    hash_to_array(&eh, edges);
    hash_free(&eh);
    size_t nsurv = 0;
    for (size_t i = 0; i < ne; i++) if (edges[i].net != 0) edges[nsurv++] = edges[i];
    qsort(edges, nsurv, sizeof(s_edgerec), edgerec_cmp_plane);

    /* 4. Walk each plane's surviving edges into closed loops, stored as
     *    concatenated vertex ids (lbuf) with per-loop offsets (loff). */
    int  *lbuf = (int *)malloc(sizeof(int) * (nsurv ? nsurv : 1));
    int  *loff = (int *)malloc(sizeof(int) * (nsurv + 2));
    int  *lb2  = (int *)malloc(sizeof(int) * (nsurv ? nsurv : 1));   /* scratch for compaction */
    int  *lo2  = (int *)malloc(sizeof(int) * (nsurv + 2));
    char *used = (char *)calloc(nsurv ? nsurv : 1, 1);
    if (!lbuf || !loff || !lb2 || !lo2 || !used) {
        free(lbuf); free(loff); free(lb2); free(lo2); free(used);
        free(edges); free(coords); return -1;
    }
    int nloops = 0, lpos = 0; loff[0] = 0;
    for (size_t b = 0; b < nsurv; ) {
        size_t e = b + 1;
        while (e < nsurv && edgerec_cmp_plane(&edges[b], &edges[e]) == 0) e++;
        /* block [b,e): edges of one plane. Directed u->v = (net>0)?(lo,hi):(hi,lo). */
        for (size_t s = b; s < e; s++) {
            if (used[s]) continue;
            int start = lpos, bad = 0, guard = 0;
            size_t cur = s;
            do {
                used[cur] = 1;
                int u = edges[cur].net > 0 ? edges[cur].lo : edges[cur].hi;
                int v = edges[cur].net > 0 ? edges[cur].hi : edges[cur].lo;
                lbuf[lpos++] = u;
                size_t nxt = e;                           /* next edge whose tail == v */
                for (size_t q = b; q < e; q++) {
                    if (used[q]) continue;
                    int qu = edges[q].net > 0 ? edges[q].lo : edges[q].hi;
                    if (qu == v) { nxt = q; break; }
                }
                if (nxt == e) { if (v != lbuf[start]) bad = 1; break; }
                cur = nxt;
            } while (guard++ < (int)nsurv + 4);
            if (bad || lpos - start < 3) { lpos = start; continue; }  /* drop bad/degenerate */
            loff[++nloops] = lpos;
        }
        b = e;
    }
    free(used); free(edges);

    /* 5. Globally strip pure-artifact collinear vertices: a vertex collinear
     *    exactly TWO distinct neighbours over the whole surface is a false
     *    vertex -- a tet-subdivision point sitting on a straight cell edge
     *    between those two neighbours (real cell corners have degree >= 3, where
     *    >= 3 facets meet). Splicing it out merges its two edges into one; both
     *    facets sharing the edge see the same degree-2 vertex and drop it
     *    identically, so the edge stays matched. This is purely combinatorial --
     *    no coordinates, no tolerance. Iterate to collapse chains. */
    int  *nb0 = (int *)malloc(sizeof(int) * (Nv ? Nv : 1));
    int  *nb1 = (int *)malloc(sizeof(int) * (Nv ? Nv : 1));
    char *many = (char *)malloc(Nv ? Nv : 1);
    char *rem  = (char *)malloc(Nv ? Nv : 1);
    if (!nb0 || !nb1 || !many || !rem) {
        free(nb0); free(nb1); free(many); free(rem);
        free(lbuf); free(loff); free(lb2); free(lo2); free(coords); return -1;
    }
    for (int pass = 0; pass < Nv + 1; pass++) {
        for (int i = 0; i < Nv; i++) { nb0[i] = nb1[i] = -1; many[i] = 0; }
        for (int L = 0; L < nloops; L++) {
            int off = loff[L], len = loff[L+1] - loff[L];
            for (int i = 0; i < len; i++) {
                int v = lbuf[off+i];
                int adj[2] = { lbuf[off+(i-1+len)%len], lbuf[off+(i+1)%len] };
                for (int k = 0; k < 2; k++) {
                    int x = adj[k];
                    if (nb0[v] == -1 || nb0[v] == x) nb0[v] = x;
                    else if (nb1[v] == -1 || nb1[v] == x) nb1[v] = x;
                    else many[v] = 1;
                }
            }
        }
        int nrem = 0;
        for (int i = 0; i < Nv; i++) {
            rem[i] = (!many[i] && nb0[i] != -1 && nb1[i] != -1);   /* exactly 2 neighbours */
            if (rem[i]) nrem++;
        }
        if (nrem == 0) break;
        int wpos = 0, wl = 0; lo2[0] = 0;               /* compact loops, dropping removed */
        for (int L = 0; L < nloops; L++) {
            int off = loff[L], len = loff[L+1] - loff[L], st = wpos;
            for (int i = 0; i < len; i++) { int v = lbuf[off+i]; if (!rem[v]) lb2[wpos++] = v; }
            if (wpos - st >= 3) lo2[++wl] = wpos; else wpos = st;
        }
        { int *tb = lbuf; lbuf = lb2; lb2 = tb; int *to = loff; loff = lo2; lo2 = to; }
        nloops = wl;
    }
    free(nb0); free(nb1); free(many); free(rem);

    /* 6. Triangulate each cleaned loop by exact ear-clipping. */
    s_dynarray faces = dynarray_initialize(sizeof(int) * 3, 64);
    if (!faces.items) { free(lbuf); free(loff); free(lb2); free(lo2); free(coords); return -1; }
    for (int L = 0; L < nloops; L++) {
        if (surf_earclip(coords, lbuf + loff[L], loff[L+1] - loff[L], &faces) < 0) {
            dynarray_free(&faces);
            free(lbuf); free(loff); free(lb2); free(lo2); free(coords); return -1;
        }
    }
    free(lbuf); free(loff); free(lb2); free(lo2);

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
        if (v->cnt == 0) v->t0 = t; else if (v->cnt == 1) v->t1 = t;
        v->cnt++;
    }
    if (!herr) {
        size_t nk = e2h.size;
        s_e2tval *ev = (s_e2tval *)malloc(sizeof(s_e2tval) * (nk ? nk : 1));
        if (!ev) herr = 1;
        else {
            hash_to_array(&e2h, ev);
            for (size_t i = 0; i < nk; i++)
                if (ev[i].cnt == 2) UF_union(ev[i].t0, ev[i].t1, parent, rank);
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
        if (trimesh_is_valid(&m)) meshes[nmesh++] = m;   /* skip any residual bad component */
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
