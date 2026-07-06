/* WARNING: This was built with claude code */
#include "vdiagram.h"
#include "vor3d.h"
#include "delaunay.h"
#include "scplx.h"
#include "dynarray.h"
#include "trimesh.h"
#include "clip_lp.h"     /* Tier-B: 3D LP-feasibility dispatch (lp3_*) */
#include "linalg.h"      /* solve_3x3_ppivot */
#include "surf_vkey.h"   /* combinatorial vertex keys + per-cell surface assembly */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>


/* -----------------------------------------------------------------------
 * Internal validity / free helpers
 * ----------------------------------------------------------------------- */

int ncvx_domain_is_valid(const s_ncvx_domain *d)
{
    return trimesh_is_valid(&d->surface) &&
           d->cdt.head != NULL &&
           convhull_is_valid(&d->bpoly.convh);
}

void free_ncvx_domain(s_ncvx_domain *d)
{
    free_trimesh(&d->surface);
    free_complex(&d->cdt);
    free_bpoly(&d->bpoly);
    memset(d, 0, sizeof(s_ncvx_domain));
}

void free_ncvx_vdiagram(s_ncvx_vdiagram *vd)
{
    int N = vd->seeds.N;   /* read before free_points zeroes the struct */
    free_points(&vd->seeds);
    free_ncvx_domain(&vd->domain);
    if (vd->vcells) {
        for (int i = 0; i < N; i++) free_vcell(&vd->vcells[i]);
        free(vd->vcells);
    }
    memset(vd, 0, sizeof(s_ncvx_vdiagram));
}


/* ---------------------------------------------------------------------- */

s_ncvx_domain ncvx_domain_from_trimesh(const s_trimesh *mesh,
                                        double EPS_DEG, double TOL, int keep_exterior)
{
    if (!trimesh_is_valid(mesh))
        return (s_ncvx_domain){0};

    /* Build the CDT of the interior. keep_exterior=0 returns only interior tets;
     * keep_exterior=1 returns the full convex-hull complex (interior + exterior
     * pockets flagged via nc->interior), which the clip treats identically and
     * which also supports walk-based point location. The scplx is kept as the
     * domain (traversed by the clip); not extracted. */
    s_scplx dt = keep_exterior ? tetrahedralize_domain_flagged(mesh, EPS_DEG, TOL)
                               : tetrahedralize_interior_trimesh(mesh, EPS_DEG, TOL);
    if (!dt.head) {
        fprintf(stderr, "ncvx_domain_from_trimesh: tetrahedralization failed\n");
        return (s_ncvx_domain){0};
    }

    double volume_sum = 0.0;
    for (s_ncell *nc = dt.head; nc; nc = nc->next) {
        if (!nc->interior) continue;              /* skip exterior pockets, if kept */
        s_point v[4];
        extract_vertices_ncell(&dt, nc, v);
        volume_sum += fabs(signed_volume_tetra(v));
    }

    s_trimesh surface = copy_trimesh(mesh);
    if (!trimesh_is_valid(&surface)) { free_complex(&dt); return (s_ncvx_domain){0}; }

    s_bpoly bpoly = bpoly_from_points(&mesh->points, EPS_DEG);
    if (!convhull_is_valid(&bpoly.convh)) {
        free_trimesh(&surface); free_complex(&dt); return (s_ncvx_domain){0};
    }

    s_ncvx_domain domain = {0};
    domain.surface       = surface;
    domain.cdt           = dt;            /* take ownership */
    domain.domain_volume = volume_sum;
    domain.bpoly         = bpoly;
    domain.spatial_index = NULL;
    return domain;
}


/* ---------------------------------------------------------------------- */

/* Cell adjacency (Delaunay dual), reserved for Step 3's canonical bisector cells:
 * cell i and cell j are adjacent iff {i,j} is a Delaunay edge in the seed DT,
 * i.e. they share a Voronoi face F_ij. CSR layout: cell i's neighbors are
 * neighbor[offset[i] .. offset[i+1]-1]. Real seeds occupy DT vertex ids
 * 0..N-1 (see vor3d_in_bp_dt's contract); mirror/auxiliary points (id >= N)
 * are dropped since they are not real Voronoi cells. */
struct ncvx_cell_adjacency {
    int N;
    int *offset;    /* size N+1 */
    int *neighbor;  /* size offset[N] */
    /* Voronoi duals per cell, precomputed together with the adjacency so each
     * seed's incident-cell star is walked once (not again per cell in the clip).
     * CSR: cell i's data is vtx[3*vtx_off[i] ..], edge[4*edge_off[i] ..]. */
    int *vtx_off;   /* size N+1 */
    int *vtx;       /* 3 * vtx_off[N] ints: Voronoi-vertex seed triples */
    int *edge_off;  /* size N+1 */
    int *edge;      /* 4 * edge_off[N] ints: {a, b, apex1, apex2} per Voronoi edge */
};

static int cmp_int_asc(const void *a, const void *b)
{ int x = *(const int *)a, y = *(const int *)b; return (x > y) - (x < y); }

static int cmp_edge_ab(const void *a, const void *b)
{
    const int *x = a, *y = b;
    if (x[0] != y[0]) return (x[0] > y[0]) - (x[0] < y[0]);
    return (x[1] > y[1]) - (x[1] < y[1]);
}

/* Build cell adjacency + Voronoi duals in a SINGLE sweep over the seed-DT tets,
 * scattering each tet's contribution to its four vertices -- replacing the former
 * per-seed incident-star walk (which visited every tet ~4x with mark/scratch
 * churn and an O(deg^2) per-seed edge dedup). Two passes (count, then fill) build
 * duplicate-laden per-seed scatter ranges; neighbours are then sorted+uniqued and
 * edges grouped by (a,b) with their <=2 apex caps merged. Same CSR contents as
 * before, up to ordering (which the clip does not depend on). */
static int build_cell_adjacency(s_scplx *dt, int N, struct ncvx_cell_adjacency *adj)
{
    memset(adj, 0, sizeof(*adj));
    int ok = 0;
    int *nc = NULL, *vc = NULL, *ec = NULL;           /* per-seed start offsets, N+1 */
    int *ncur = NULL, *vcur = NULL, *ecur = NULL;     /* fill cursors, N */
    int *nbr_s = NULL, *edge_s = NULL, *vtx_s = NULL; /* duplicate-laden scatter arrays */
    s_dynarray fnbr = {0}, fedge = {0};

    nc = calloc((size_t)(N + 1), sizeof(int));
    vc = calloc((size_t)(N + 1), sizeof(int));
    ec = calloc((size_t)(N + 1), sizeof(int));
    adj->offset   = malloc(sizeof(int) * (size_t)(N + 1));
    adj->vtx_off  = malloc(sizeof(int) * (size_t)(N + 1));
    adj->edge_off = malloc(sizeof(int) * (size_t)(N + 1));
    if (!nc || !vc || !ec || !adj->offset || !adj->vtx_off || !adj->edge_off) goto done;

    /* Pass 1: per-seed counts (nc/ec are upper bounds incl. duplicates; vc exact). */
    for (s_ncell *t = dt->head; t; t = t->next) {
        for (int vi = 0; vi < 4; vi++) {
            int v = t->vertex_id[vi];
            if (v < 0 || v >= N) continue;
            int o[3], no = 0, nreal = 0;
            for (int k = 0; k < 4; k++)
                if (k != vi) { o[no] = t->vertex_id[k]; if (o[no] < N) nreal++; no++; }
            if (nreal == 3) vc[v]++;
            nc[v] += nreal;
            if (o[0] < N && o[1] < N) ec[v]++;
            if (o[0] < N && o[2] < N) ec[v]++;
            if (o[1] < N && o[2] < N) ec[v]++;
        }
    }
    { int nt = 0, vt = 0, et = 0;
      for (int i = 0; i < N; i++) {
          int a = nc[i]; nc[i] = nt; nt += a;
          int b = vc[i]; vc[i] = vt; vt += b;
          int c = ec[i]; ec[i] = et; et += c;
      }
      nc[N] = nt; vc[N] = vt; ec[N] = et;
      nbr_s  = nt ? malloc(sizeof(int) * (size_t)nt)     : NULL;
      vtx_s  = vt ? malloc(sizeof(int) * (size_t)vt * 3) : NULL;
      edge_s = et ? malloc(sizeof(int) * (size_t)et * 3) : NULL;   /* {a,b,apex} records */
      if ((nt && !nbr_s) || (vt && !vtx_s) || (et && !edge_s)) goto done;
    }

    ncur = malloc(sizeof(int) * (size_t)(N > 0 ? N : 1));
    vcur = malloc(sizeof(int) * (size_t)(N > 0 ? N : 1));
    ecur = malloc(sizeof(int) * (size_t)(N > 0 ? N : 1));
    if (!ncur || !vcur || !ecur) goto done;
    for (int i = 0; i < N; i++) { ncur[i] = nc[i]; vcur[i] = vc[i]; ecur[i] = ec[i]; }

    /* Pass 2: scatter. */
    for (s_ncell *t = dt->head; t; t = t->next) {
        for (int vi = 0; vi < 4; vi++) {
            int v = t->vertex_id[vi];
            if (v < 0 || v >= N) continue;
            int o[3], no = 0;
            for (int k = 0; k < 4; k++) if (k != vi) o[no++] = t->vertex_id[k];
            int r0 = o[0] < N, r1 = o[1] < N, r2 = o[2] < N;
            if (r0 && r1 && r2) { int *d = vtx_s + 3 * vcur[v]++; d[0]=o[0]; d[1]=o[1]; d[2]=o[2]; }
            if (r0) nbr_s[ncur[v]++] = o[0];
            if (r1) nbr_s[ncur[v]++] = o[1];
            if (r2) nbr_s[ncur[v]++] = o[2];
            if (r0 && r1) { int a=o[0], b=o[1]; if (a>b){int t2=a;a=b;b=t2;} int *e=edge_s+3*ecur[v]++; e[0]=a;e[1]=b;e[2]=o[2]; }
            if (r0 && r2) { int a=o[0], b=o[2]; if (a>b){int t2=a;a=b;b=t2;} int *e=edge_s+3*ecur[v]++; e[0]=a;e[1]=b;e[2]=o[1]; }
            if (r1 && r2) { int a=o[1], b=o[2]; if (a>b){int t2=a;a=b;b=t2;} int *e=edge_s+3*ecur[v]++; e[0]=a;e[1]=b;e[2]=o[0]; }
        }
    }

    /* vtx needs no dedup: transfer the scatter array directly. */
    memcpy(adj->vtx_off, vc, sizeof(int) * (size_t)(N + 1));
    adj->vtx = vtx_s; vtx_s = NULL;

    /* neighbours: sort+unique per seed. edges: group by (a,b), merge <=2 apexes. */
    fnbr  = dynarray_initialize(sizeof(int),    (size_t)(nc[N] ? nc[N] : 1));
    fedge = dynarray_initialize(sizeof(int[4]), (size_t)(ec[N] ? ec[N] : 1));
    if (!fnbr.items || !fedge.items) goto done;
    for (int v = 0; v < N; v++) {
        adj->offset[v]   = (int)fnbr.N;
        adj->edge_off[v] = (int)fedge.N;

        int lo = nc[v], hi = nc[v + 1];
        if (hi > lo) {
            qsort(nbr_s + lo, (size_t)(hi - lo), sizeof(int), cmp_int_asc);
            int have = 0, prev = 0;
            for (int k = lo; k < hi; k++)
                if (!have || nbr_s[k] != prev) {
                    if (!dynarray_push(&fnbr, &nbr_s[k])) goto done;
                    prev = nbr_s[k]; have = 1;
                }
        }

        int elo = ec[v], ehi = ec[v + 1];
        if (ehi > elo) {
            qsort(edge_s + 3 * elo, (size_t)(ehi - elo), 3 * sizeof(int), cmp_edge_ab);
            int k = elo;
            while (k < ehi) {
                int a = edge_s[3*k], b = edge_s[3*k+1], ap1 = -1, ap2 = -1, kk = k;
                while (kk < ehi && edge_s[3*kk] == a && edge_s[3*kk+1] == b) {
                    int ap = edge_s[3*kk+2];
                    if (ap1 == -1) ap1 = ap;
                    else if (ap != ap1 && ap2 == -1) ap2 = ap;
                    kk++;
                }
                int rec[4] = { a, b, ap1, ap2 };
                if (!dynarray_push(&fedge, rec)) goto done;
                k = kk;
            }
        }
    }
    adj->offset[N]   = (int)fnbr.N;
    adj->edge_off[N] = (int)fedge.N;
    adj->N = N;

    adj->neighbor = fnbr.N  ? malloc(sizeof(int)    * fnbr.N)  : NULL;
    adj->edge     = fedge.N ? malloc(sizeof(int[4]) * fedge.N) : NULL;
    if ((fnbr.N && !adj->neighbor) || (fedge.N && !adj->edge)) goto done;
    if (fnbr.N)  memcpy(adj->neighbor, fnbr.items,  sizeof(int)    * fnbr.N);
    if (fedge.N) memcpy(adj->edge,     fedge.items, sizeof(int[4]) * fedge.N);
    ok = 1;

done:
    free(nc); free(vc); free(ec);
    free(ncur); free(vcur); free(ecur);
    free(nbr_s); free(edge_s); free(vtx_s);
    dynarray_free(&fnbr); dynarray_free(&fedge);
    if (!ok) {
        free(adj->offset); free(adj->neighbor);
        free(adj->vtx_off); free(adj->vtx);
        free(adj->edge_off); free(adj->edge);
        memset(adj, 0, sizeof(*adj));
    }
    return ok;
}

static void free_cell_adjacency(struct ncvx_cell_adjacency *adj)
{
    free(adj->offset);
    free(adj->neighbor);
    free(adj->vtx_off);
    free(adj->vtx);
    free(adj->edge_off);
    free(adj->edge);
    memset(adj, 0, sizeof(*adj));
}

/* Debug/verification only: a proper Delaunay dual is symmetric -- j in
 * neighbors(i) implies i in neighbors(j). Returns 1 if consistent, 0 if not
 * (indicates a bug in build_cell_adjacency, not a data issue). */
static int cell_adjacency_is_symmetric(const struct ncvx_cell_adjacency *adj)
{
    for (int i = 0; i < adj->N; i++) {
        for (int k = adj->offset[i]; k < adj->offset[i + 1]; k++) {
            int j = adj->neighbor[k];
            int found = 0;
            for (int k2 = adj->offset[j]; k2 < adj->offset[j + 1] && !found; k2++)
                if (adj->neighbor[k2] == i) found = 1;
            if (!found) return 0;
        }
    }
    return 1;
}


/* Nearest real seed (== Voronoi-cell owner) to point p. A seed-DT walk gives a
 * starting candidate (for p outside the seed hull it halts at the nearest hull
 * tet); a greedy descent over the Voronoi adjacency then reaches the true
 * nearest seed. Returns an id in [0,N). The owner of p's cell always contains p,
 * so seeding the flood from it guarantees overlap. */
static double nseed_d2(s_point p, s_point q)
{ double dx=p.x-q.x, dy=p.y-q.y, dz=p.z-q.z; return dx*dx+dy*dy+dz*dz; }

static int nearest_seed(s_scplx *seed_dt, const s_points *rs,
                        const struct ncvx_cell_adjacency *adj, int N, s_point p)
{
    int cur = -1; double cd = 0.0;
    s_ncell *loc = in_ncell_walk(seed_dt, p);
    if (loc)
        for (int k = 0; k < 4; k++) {
            int id = loc->vertex_id[k];
            if (id < 0 || id >= N) continue;
            double d = nseed_d2(p, rs->p[id]);
            if (cur < 0 || d < cd) { cur = id; cd = d; }
        }
    if (cur < 0) {                     /* walk gave nothing -> brute nearest */
        for (int i = 0; i < N; i++) {
            double d = nseed_d2(p, rs->p[i]);
            if (cur < 0 || d < cd) { cur = i; cd = d; }
        }
        return cur;
    }
    for (;;) {                         /* greedy descent to the true nearest */
        int best = cur; double bd = cd;
        for (int k = adj->offset[cur]; k < adj->offset[cur + 1]; k++) {
            int j = adj->neighbor[k];
            if (j < 0 || j >= N) continue;
            double d = nseed_d2(p, rs->p[j]);
            if (d < bd) { best = j; bd = d; }
        }
        if (best == cur) break;
        cur = best; cd = bd;
    }
    return cur;
}

/* 2a: exact owner cell of every CDT vertex, precomputed once (indexed by dense
 * cdt point id). owner[v] = the seed whose Voronoi cell contains cdt point v
 * (== the nearest seed; by the Voronoi-cell theorem, "closer to s than every
 * neighbour of s" == "closer to s than every seed"). tie[v] = 1 iff v is
 * equidistant to owner and some neighbour (lies on a Voronoi face) -- category
 * (b) must then keep v for every incident cell, so it falls back to the exact
 * per-bisector loop. A float greedy descent (nearest_seed) seeds an EXACT
 * monotone descent over the Voronoi adjacency using lp_feasible_T0_T1_S (the
 * very predicate category (b) uses), so owner[]/tie[] reproduce (b) bit-for-bit. */
static int compute_vertex_owners(const s_ncvx_domain *domain, const s_points *rs,
                                 s_scplx *seed_dt, const struct ncvx_cell_adjacency *adj,
                                 int N, int **out_owner, char **out_tie)
{
    int P = domain->cdt.points.N;
    int  *owner = malloc(sizeof(int)  * (size_t)P);
    char *tie   = calloc((size_t)P, 1);
    if (!owner || !tie) { free(owner); free(tie); return 0; }

    for (int v = 0; v < P; v++) {
        s_point V = domain->cdt.points.p[v];
        int m = nearest_seed(seed_dt, rs, adj, N, V);      /* float nearest as start */
        if (m < 0) m = 0;
        for (;;) {                                          /* exact monotone descent */
            int moved = -1;
            s_point sm = rs->p[m];
            for (int k = adj->offset[m]; k < adj->offset[m+1]; k++) {
                int j = adj->neighbor[k];
                if (j < 0 || j >= N) continue;
                s_point sj = rs->p[j];
                if (lp_feasible_T0_T1_S(V.x,V.y,V.z, sm.x,sm.y,sm.z, sj.x,sj.y,sj.z) < 0) {
                    moved = j; break;                       /* V strictly closer to j */
                }
            }
            if (moved < 0) break;
            m = moved;
        }
        char t = 0;                                         /* tie test at final owner */
        s_point sm = rs->p[m];
        for (int k = adj->offset[m]; k < adj->offset[m+1]; k++) {
            int j = adj->neighbor[k];
            if (j < 0 || j >= N) continue;
            s_point sj = rs->p[j];
            if (lp_feasible_T0_T1_S(V.x,V.y,V.z, sm.x,sm.y,sm.z, sj.x,sj.y,sj.z) == 0) { t = 1; break; }
        }
        owner[v] = m; tie[v] = t;
    }
    *out_owner = owner; *out_tie = tie;
    return 1;
}

/* |X-s|^2 - |X-t|^2 (bisector value of s->t at X, in double). */
static double lp_bisval(s_point X, s_point s, s_point t)
{
    s_point xs = subtract_points(X, s), xt = subtract_points(X, t);
    return dot_prod(xs, xs) - dot_prod(xt, xt);
}

/* Exact sign(lam_j - lam_k) for two bisectors s->tj, s->tk that cross tet edge
 * AB, where lam_m = g_mA/(g_mA - g_mB) is the crossing parameter (g_mX the
 * bisector value at X). From lp3_predicates.tex:
 *   sign(lam_j - lam_k) = -sign(det) * slope_j * slope_k
 * with det = lp3_det_TTS, slope = lp3_slope_TTS. Used by the (d) 1D-LP. */
static int lp_cmp_lambda(s_point A, s_point B, s_point s, s_point tj, s_point tk)
{
    int det = lp3_det_TTS(A.x,A.y,A.z, B.x,B.y,B.z, s.x,s.y,s.z,
                          tj.x,tj.y,tj.z, tk.x,tk.y,tk.z);
    int sj  = lp3_slope_TTS(A.x,A.y,A.z, B.x,B.y,B.z, s.x,s.y,s.z, tj.x,tj.y,tj.z);
    int sk  = lp3_slope_TTS(A.x,A.y,A.z, B.x,B.y,B.z, s.x,s.y,s.z, tk.x,tk.y,tk.z);
    return -det * sj * sk;
}

/* ---------------------------------------------------------------------- */
/* Surface extraction helpers (see SURFACE_EXTRACTION_PLAN.md)             */
/* ---------------------------------------------------------------------- */

/* Scatter per-candidate labels (aligned index-for-index with the hull's input
 * point array) onto the hull's output vertices via pmap. Two inputs deduped to
 * one output vertex are merged. Returns malloc'd array of length out_N (caller
 * frees), or NULL on allocation failure. */
static s_vlabel *scatter_labels(const s_dynarray *cand, const int *pmap,
                                int in_N, int out_N)
{
    size_t nn = (size_t)(out_N > 0 ? out_N : 1);
    s_vlabel *labels = calloc(nn, sizeof(s_vlabel));
    char     *set    = calloc(nn, 1);
    if (!labels || !set) { free(labels); free(set); return NULL; }
    for (int ii = 0; ii < in_N; ii++) {
        int v = pmap[ii];
        if (v < 0) continue;
        const s_vlabel *src = (const s_vlabel *)dynarray_get_ptr_c(cand, (size_t)ii);
        if (!set[v]) { labels[v] = *src; set[v] = 1; }
        else vlabel_merge(&labels[v], src);
    }
    free(set);
    return labels;
}

/* Emit every hull face of `piece` (labels aligned with piece->points) into the
 * cell accumulator `acc` (items: s_surf_tri), tagged with its supporting plane:
 *   - a face whose 3 vertices share a common tet face f (exact plane membership)
 *     lies on that tet face -> its GLOBAL CDT face-triangle (from nc->vertex_id);
 *   - otherwise the face lies on the Voronoi bisector {seed_id, common partner}.
 * No surface/internal decision is made here -- interior faces cancel in the
 * per-plane edge XOR of build_cell_surface. `nc` is the tet the piece was carved
 * from (NULL for an interior whole-cell hull; all its faces are bisectors).
 * Returns 1 OK, 0 allocation error. */
static int piece_to_surface(const s_convh *piece, const s_vlabel *labels,
                            const s_ncell *nc, int seed_id, s_dynarray *acc)
{
    for (int fc = 0; fc < piece->Nf; fc++) {
        int v0 = piece->faces[3*fc+0];
        int v1 = piece->faces[3*fc+1];
        int v2 = piece->faces[3*fc+2];
        const s_vlabel *L0 = &labels[v0], *L1 = &labels[v1], *L2 = &labels[v2];

        uint8_t fm = L0->fmask & L1->fmask & L2->fmask;
        s_pkey plane;
        if (fm && nc) {                                   /* common tet face */
            int f = __builtin_ctz(fm);
            int fv[3]; lp3_face_idx(f, fv);
            plane = pkey_face(nc->vertex_id[fv[0]], nc->vertex_id[fv[1]],
                              nc->vertex_id[fv[2]]);
        } else {                                          /* Voronoi bisector wall */
            int t = vlabel_common_bis(L0, L1, L2);
            if (t < 0) continue;                          /* degenerate face: skip */
            plane = pkey_bis(seed_id, t);
        }
        s_surf_tri tri = { plane, {
            { L0->key, piece->points.p[v0] },
            { L1->key, piece->points.p[v1] },
            { L2->key, piece->points.p[v2] },
        } };
        if (!dynarray_push(acc, &tri)) return 0;
    }
    return 1;
}

/* Tier-B exact clip of the Voronoi cell of seed `seed_id` against one tetrahedron.
 * Enumerates the four vertex categories of (cell ^ tet) directly on the constraint
 * set { 4 tet faces T ; neighbour bisectors S }, keeping each candidate iff it is
 * feasible (lp3_feasible) against every other constraint, and hulls the survivors:
 *   (b) tet vertex in cell        (TTT vs every S)
 *   (d) tet edge ^ bisector       (TTS vs every other S; faces automatic)
 *   (a) Voronoi vertex in tet      (SSS vs 4 faces + every other S) -- circumcentre
 *   (c) Voronoi edge ^ tet face    (TSS vs other 3 faces + every other S)
 * (a)/(c) are enumerated brute-force over neighbour triples/pairs; DT-dual
 * (output-sensitive) enumeration is a later optimization. Coordinates are float;
 * every keep/drop decision is exact. Returns 1 (built), 0 (empty), -1 (error). */
static int clip_vcell_clip_tet_lp(int seed_id, const s_points *seeds,
                                  const struct ncvx_cell_adjacency *adj,
                                  const int *vor_vtx, int n_vtx,   /* n_vtx * {t1,t2,t3}     */
                                  const int *vor_edge, int n_edge, /* n_edge * {a,b,ap1,ap2} */
                                  const s_point tet[4], const int tet_vids[4],
                                  const int *owner, const char *tie,
                                  int *resolved_vtx, int tet_id,
                                  double EPS_DEG, s_convh *out,
                                  s_vlabel **out_labels)
{
    *out = convhull_NAN;
    if (out_labels) *out_labels = NULL;
    s_point s = seeds->p[seed_id];
    int NB = adj->offset[seed_id + 1] - adj->offset[seed_id];
    const int *nbr = adj->neighbor + adj->offset[seed_id];   /* NB neighbour seed indices */

    int sigma_f[4];
    lp3_tet_face_orient(tet, sigma_f);

    s_cid3 FC[4];
    for (int f = 0; f < 4; f++) { FC[f].type = CT_TET; FC[f].idx = f; }

    s_dynarray pts  = dynarray_initialize(sizeof(s_point), 64);
    s_dynarray cand = dynarray_initialize(sizeof(s_vlabel), 64);  /* labels aligned with pts */
    if (!pts.items || !cand.items) goto err;

    /* ---- (b) tet vertices inside the cell (2a: owner[] memo) ---- */
    for (int j = 0; j < 4; j++) {
        int V = tet_vids[j];
        int keep;
        if (!tie[V]) {
            keep = (owner[V] == seed_id);   /* exact: V in cell iff its owner is s */
        } else {                            /* boundary vertex -> exact per-bisector fallback */
            keep = 1;
            s_cid3 tri[3]; int n = 0;
            for (int f = 0; f < 4; f++) if (f != j) tri[n++] = FC[f];  /* 3 faces meeting at j */
            for (int k = 0; k < NB && keep; k++) {
                s_cid3 q = { CT_SEED, nbr[k] };
                if (lp3_feasible(tet, sigma_f, s, seeds, tri[0], tri[1], tri[2], q) < 0) keep = 0;
            }
        }
        if (keep) {
            s_vlabel lb = { vkey_TV(tet_vids[j]), (uint8_t)(0xF ^ (1 << j)), {-1, -1, -1} };
            if (!dynarray_push(&pts, &tet[j]) || !dynarray_push(&cand, &lb)) goto err;
        }
    }

    /* ---- (d) tet edge ^ bisector: 1D LP on each of the 6 edges ----
     * On edge [A,B] the cell is an interval; its interior endpoints are the
     * tightest LO (max lambda) and HI (min lambda) bisector crossings. One O(N)
     * pass (exact slope classification + ordering) vs the O(N^2) per-crossing
     * feasibility filter. Ends at lambda=0/1 are tet vertices (category b). */
    for (int f = 0; f < 4; f++) for (int g = f + 1; g < 4; g++) {
        int e[2]; lp3_shared_edge(f, g, e);
        s_point A = tet[e[0]], B = tet[e[1]];
        int lo_k = -1, hi_k = -1, empty = 0;   /* binding LO/HI neighbour index (-1 = segment end) */
        for (int k = 0; k < NB && !empty; k++) {
            s_point t = seeds->p[nbr[k]];
            /* exact signs of gX = |X-s|^2 - |X-t|^2 (>0 == X outside this bisector) */
            int sA = -lp_feasible_T0_T1_S(A.x,A.y,A.z, s.x,s.y,s.z, t.x,t.y,t.z);
            int sB = -lp_feasible_T0_T1_S(B.x,B.y,B.z, s.x,s.y,s.z, t.x,t.y,t.z);
            if (sA < 0 && sB < 0) continue;              /* edge fully inside this bisector */
            if (sA > 0 && sB > 0) { empty = 1; continue; }  /* edge fully outside -> no (d) here */
            if (sA == 0 || sB == 0) continue;            /* endpoint on bisector: no interior crossing */
            /* strict straddle -> crossing in (0,1); slope>0 (g decreasing) => LO bound */
            int slope = lp3_slope_TTS(A.x,A.y,A.z, B.x,B.y,B.z, s.x,s.y,s.z, t.x,t.y,t.z);
            if (slope > 0) {
                if (lo_k < 0 || lp_cmp_lambda(A,B,s, t, seeds->p[nbr[lo_k]]) > 0) lo_k = k;
            } else if (slope < 0) {
                if (hi_k < 0 || lp_cmp_lambda(A,B,s, t, seeds->p[nbr[hi_k]]) < 0) hi_k = k;
            }
        }
        if (empty) continue;
        if (lo_k >= 0 && hi_k >= 0 &&
            lp_cmp_lambda(A,B,s, seeds->p[nbr[lo_k]], seeds->p[nbr[hi_k]]) > 0) continue; /* empty interval */
        for (int side = 0; side < 2; side++) {
            int kk = side ? hi_k : lo_k;
            if (kk < 0) continue;                        /* this end is a tet vertex (category b) */
            s_point t = seeds->p[nbr[kk]];
            double gA = lp_bisval(A, s, t), gB = lp_bisval(B, s, t);
            double lam = gA / (gA - gB);
            s_point Pt = { .x = A.x + lam*(B.x-A.x), .y = A.y + lam*(B.y-A.y), .z = A.z + lam*(B.z-A.z) };
            s_vlabel lb = { vkey_VB(tet_vids[e[0]], tet_vids[e[1]], seed_id, nbr[kk]),
                            (uint8_t)((1 << f) | (1 << g)), { nbr[kk], -1, -1 } };
            if (!dynarray_push(&pts, &Pt) || !dynarray_push(&cand, &lb)) goto err;
        }
    }

    /* ---- (a) Voronoi vertices (circumcentres) inside the tet ----
     * Candidates are the DT-dual Voronoi vertices of this cell (one per all-real
     * incident DT tet), not all C(N,3) neighbour triples. The feasibility filter
     * (inside 4 faces + all other bisectors) is retained. */
    for (int vi = 0; vi < n_vtx; vi++) {
        int j1 = vor_vtx[3*vi+0], j2 = vor_vtx[3*vi+1], j3 = vor_vtx[3*vi+2];
        s_cid3 a = { CT_SEED, j1 }, b = { CT_SEED, j2 }, c = { CT_SEED, j3 };
        /* A DT-dual circumcentre is a Voronoi vertex, so by the empty-sphere
         * property it is inside every bisector -- only the 4 tet faces need
         * testing (no O(N) bisector filter).
         * 2b: a Voronoi vertex is strictly inside exactly one tet. resolved_vtx[vi]
         * caches that tet (0 = unresolved). Resolved elsewhere -> skip, no
         * predicates; resolved here -> keep, no predicates; else test the 4 faces
         * and memoize only on strict containment (a dual ON a face is in 2 tets,
         * left unresolved so both keep it -- exactly as the raw test does). */
        int R = resolved_vtx[vi];
        if (R && R != tet_id) continue;                 /* strictly inside another tet */
        int keep = 1;
        if (R != tet_id) {                              /* R == 0: unresolved -> test */
            int strict = 1;
            for (int f = 0; f < 4 && keep; f++) {
                int fe = lp3_feasible(tet, sigma_f, s, seeds, a, b, c, FC[f]);
                if (fe < 0) keep = 0;
                else if (fe == 0) strict = 0;
            }
            if (keep && strict) resolved_vtx[vi] = tet_id;
        }
        if (keep) {
            s_point tetc[4] = { s, seeds->p[j1], seeds->p[j2], seeds->p[j3] };
            s_point Pt;
            if (circumcentre_tetrahedron(tetc, EPS_DEG, &Pt)) {
                s_vlabel lb = { vkey_VV(seed_id, j1, j2, j3), 0, { j1, j2, j3 } };
                if (!dynarray_push(&pts, &Pt) || !dynarray_push(&cand, &lb)) goto err;
            }
        }
    }

    /* ---- (c) Voronoi edge ^ tet face ----
     * Candidates are the DT-dual Voronoi edges of this cell (one per all-real
     * incident DT triangle), not all C(N,2) neighbour pairs. Feasibility filter
     * (inside other 3 faces + all other bisectors -> confirms the crossing lies on
     * the real edge segment) is retained. */
    for (int ei = 0; ei < n_edge; ei++) {
        int j1 = vor_edge[4*ei+0], j2 = vor_edge[4*ei+1];
        int cap[2] = { vor_edge[4*ei+2], vor_edge[4*ei+3] };
        s_point t1 = seeds->p[j1], t2 = seeds->p[j2];
        s_point n1 = subtract_points(t1, s); double d1 = 0.5*(dot_prod(t1,t1) - dot_prod(s,s));
        s_point n2 = subtract_points(t2, s); double d2 = 0.5*(dot_prod(t2,t2) - dot_prod(s,s));
        for (int f = 0; f < 4; f++) {
            s_cid3 a = FC[f], b = { CT_SEED, j1 }, c = { CT_SEED, j2 };
            int keep = 1;
            for (int g = 0; g < 4 && keep; g++) {
                if (g == f) continue;
                if (lp3_feasible(tet, sigma_f, s, seeds, a, b, c, FC[g]) < 0) keep = 0;
            }
            /* On the real Voronoi edge iff on s's side of its (<=2) capping
             * bisectors (the adjacent DT tets' apexes) -- duality replaces the
             * O(N) all-bisector filter with O(1). */
            for (int ci2 = 0; ci2 < 2 && keep; ci2++) {
                if (cap[ci2] < 0) continue;
                s_cid3 q = { CT_SEED, cap[ci2] };
                if (lp3_feasible(tet, sigma_f, s, seeds, a, b, c, q) < 0) keep = 0;
            }
            if (keep) {
                int fv[3]; lp3_face_idx(f, fv);
                s_point A = tet[fv[0]], Q = tet[fv[1]], R = tet[fv[2]];
                /* Exact degeneracy screen: lp3_D_TSS == sign det[n_T; t1-s; t2-s],
                 * which is 0 iff the Voronoi edge (perpendicular to triangle
                 * {s,t1,t2}) is parallel to the tet face -- then the crossing does
                 * not exist and the float solve would return garbage the feasibility
                 * test still accepts. This robust (indirect) sign is consistent
                 * across the three cells sharing the edge, so it can't leave one
                 * cell's piece overlapping another's. */
                if (lp3_D_TSS(A.x,A.y,A.z, Q.x,Q.y,Q.z, R.x,R.y,R.z,
                              s.x,s.y,s.z, t1.x,t1.y,t1.z, t2.x,t2.y,t2.z) == 0)
                    continue;
                s_point nf = cross_prod(subtract_points(Q,A), subtract_points(R,A));
                double df = dot_prod(nf, A);
                double M[3][3] = {{n1.x,n1.y,n1.z},{n2.x,n2.y,n2.z},{nf.x,nf.y,nf.z}};
                double rhs[3] = { d1, d2, df }, x[3];
                /* EPS_DEG pivot tolerance: numerical backstop. solve_3x3_ppivot
                 * returns the failing column (<3) with x unwritten on a bad pivot,
                 * so success is == 3 (all three pivots usable), not merely non-zero. */
                if (solve_3x3_ppivot(M, rhs, x, EPS_DEG) == 3) {
                    s_point Pt = { .x = x[0], .y = x[1], .z = x[2] };
                    s_vlabel lb = { vkey_EF(seed_id, j1, j2,
                                            tet_vids[fv[0]], tet_vids[fv[1]], tet_vids[fv[2]]),
                                    (uint8_t)(1 << f), { j1, j2, -1 } };
                    if (!dynarray_push(&pts, &Pt) || !dynarray_push(&cand, &lb)) goto err;
                }
            }
        }
    }

    {
        s_points P = { (int)pts.N, (s_point *)pts.items };
        int *pmap = out_labels ? malloc(sizeof(int) * (pts.N ? pts.N : 1)) : NULL;
        if (out_labels && !pmap) goto err;
        int r = convhull_from_points_mapped(&P, EPS_DEG, out, pmap);
        if (r == 1 && out_labels) {
            s_vlabel *labels = scatter_labels(&cand, pmap, (int)pts.N, out->points.N);
            if (!labels) { free(pmap); free_convhull(out); *out = convhull_NAN; goto err; }
            *out_labels = labels;
        }
        free(pmap);
        dynarray_free(&pts);
        dynarray_free(&cand);
        if (r == 0)  return 0;
        if (r == -1) { *out = convhull_NAN; return -1; }
        return 1;
    }
err:
    dynarray_free(&pts);
    dynarray_free(&cand);
    *out = convhull_NAN;
    if (out_labels) *out_labels = NULL;
    return -1;
}

/* Part II: emit an INTERIOR cell as its whole convex Voronoi cell -- the hull of
 * its dual circumcentres (Voronoi vertices). Valid because an interior cell lies
 * entirely inside the domain, so it is not carved by any tet. Returns 1 (built
 * into *out, volume in *vol_out), 0 (degenerate/empty), -1 (allocation error). */
static int emit_interior_cell(int seed_id, const s_points *seeds,
                              const int *vor_vtx, int n_vtx, double EPS_DEG,
                              s_convh *out, double *vol_out, s_vlabel **out_labels)
{
    *out = convhull_NAN; *vol_out = 0.0;
    if (out_labels) *out_labels = NULL;
    s_dynarray pts  = dynarray_initialize(sizeof(s_point), 32);
    s_dynarray cand = dynarray_initialize(sizeof(s_vlabel), 32);
    if (!pts.items || !cand.items) { dynarray_free(&pts); dynarray_free(&cand); return -1; }
    s_point s = seeds->p[seed_id];
    for (int vi = 0; vi < n_vtx; vi++) {
        int j1 = vor_vtx[3*vi+0], j2 = vor_vtx[3*vi+1], j3 = vor_vtx[3*vi+2];
        s_point tetc[4] = { s, seeds->p[j1], seeds->p[j2], seeds->p[j3] };
        s_point Pt;
        if (circumcentre_tetrahedron(tetc, EPS_DEG, &Pt)) {
            /* interior cell: all faces are Voronoi walls, so every vertex is a
             * VV point on the three bisectors {s,j1},{s,j2},{s,j3}. */
            s_vlabel lb = { vkey_VV(seed_id, j1, j2, j3), 0, { j1, j2, j3 } };
            if (!dynarray_push(&pts, &Pt) || !dynarray_push(&cand, &lb)) {
                dynarray_free(&pts); dynarray_free(&cand); return -1;
            }
        }
    }
    s_points P = { (int)pts.N, (s_point *)pts.items };
    int *pmap = out_labels ? malloc(sizeof(int) * (pts.N ? pts.N : 1)) : NULL;
    if (out_labels && !pmap) { dynarray_free(&pts); dynarray_free(&cand); return -1; }
    int r = convhull_from_points_mapped(&P, EPS_DEG, out, pmap);
    if (r == 1 && out_labels) {
        s_vlabel *labels = scatter_labels(&cand, pmap, (int)pts.N, out->points.N);
        if (!labels) { free(pmap); free_convhull(out); *out = convhull_NAN;
                       dynarray_free(&pts); dynarray_free(&cand); return -1; }
        *out_labels = labels;
    }
    free(pmap);
    dynarray_free(&pts);
    dynarray_free(&cand);
    if (r == 1) *vol_out = volume_convhull(out);
    return r;
}

/* A face f of interior tet nc is on the domain boundary iff it faces outside the
 * domain: either no tet lies across it (interior-only CDT) or the tet across it
 * is an exterior pocket (flagged CDT). Makes the clip accept either CDT form. */
static inline int face_is_domain_boundary(const s_ncell *nc, int f)
{
    return !nc->opposite[f] || !nc->opposite[f]->interior;
}

/* Inverse (tet-outer) clip of the whole domain. For each interior CDT tet, find
 * the Voronoi cells that carve it (seed-DT point location of the tet centroid +
 * vertices, then a flood fill over Voronoi adjacency, keeping neighbours of any
 * cell whose clip is non-empty) and clip each, scattering the pieces into
 * per-cell accumulators. Fills vcells[0..N-1]. No AABBs, no candidate lists.
 * The domain CDT may be interior-only or flagged (exterior pockets are skipped).
 * Returns 1 on success, 0 on allocation error. */
static int clip_domain(const s_ncvx_domain *domain, const s_points *real_seeds,
                       s_scplx *seed_dt, const struct ncvx_cell_adjacency *adj,
                       double EPS_DEG, int want_surface, s_vcell *vcells)
{
    int N = real_seeds->N;
    int  *owner    = NULL;                                /* 2a: per-CDT-vertex owner cell */
    char *tie      = NULL;
    int  *resolved = NULL;                                /* 2b: per-dual containing tet id */
    char *cellkind = NULL;                                /* Part II: 0 interior, 1 boundary */
    s_dynarray *pcs   = calloc((size_t)N, sizeof(s_dynarray));
    s_dynarray *surf  = want_surface ? calloc((size_t)N, sizeof(s_dynarray)) : NULL;
    double     *vol   = calloc((size_t)N, sizeof(double));
    int        *stamp = calloc((size_t)N, sizeof(int));   /* per-cell visit stamp */
    s_dynarray  queue = dynarray_initialize(sizeof(int), 64);
    if (!pcs || (want_surface && !surf) || !vol || !stamp || !queue.items) goto err_early;
    for (int i = 0; i < N; i++) {
        pcs[i]  = dynarray_initialize(sizeof(s_convh), 2);
        if (want_surface) surf[i] = dynarray_initialize(sizeof(s_surf_tri), 8);
    }

    if (!compute_vertex_owners(domain, real_seeds, seed_dt, adj, N, &owner, &tie))
        goto err;
    int NVD = adj->vtx_off[N];                            /* total Voronoi-vertex duals */
    if (NVD > 0 && !(resolved = calloc((size_t)NVD, sizeof(int)))) goto err;
    if (!(cellkind = calloc((size_t)N, 1))) goto err;

    /* Part II classification. A cell is BOUNDARY if it is unbounded (an unbounded
     * Voronoi edge -> apex == -1) or it can reach the domain surface; otherwise
     * INTERIOR (its whole convex Voronoi cell lies inside the domain and is emitted
     * directly, skipping all per-tet clips). The reach test uses the one-sided
     * triangle filter against boundary faces, so a cell that truly touches the
     * surface is never mis-classified interior; false positives just take the
     * slow (clipped) path. */
    for (int i = 0; i < N; i++) {
        int ne = adj->edge_off[i+1] - adj->edge_off[i];
        const int *ed = adj->edge + 4 * adj->edge_off[i];
        for (int e = 0; e < ne; e++)
            if (ed[4*e+2] == -1 || ed[4*e+3] == -1) { cellkind[i] = 1; break; }
    }

    /* Debug (NCVX_ALL_BOUNDARY=1): disable the interior-cell shortcut; every
     * cell takes the clipped path, so vol_sum must equal the tet-volume sum. */
    if (getenv("NCVX_ALL_BOUNDARY"))
        for (int i = 0; i < N; i++) cellkind[i] = 1;

    int tet_stamp = 0;

    /* Phase A: flood each boundary tet (one with a surface face) and flag every
     * reachable interior-so-far cell that may touch a boundary face. */
    for (s_ncell *nc = domain->cdt.head; nc; nc = nc->next) {
        if (!nc->interior) continue;              /* skip exterior pockets (flagged CDT) */
        int has_bf = 0;
        for (int f = 0; f < 4; f++) if (face_is_domain_boundary(nc, f)) { has_bf = 1; break; }
        if (!has_bf) continue;
        tet_stamp++;
        s_point tv[4];
        extract_vertices_ncell(&domain->cdt, nc, tv);
        s_point c = { .x = (tv[0].x+tv[1].x+tv[2].x+tv[3].x)*0.25,
                      .y = (tv[0].y+tv[1].y+tv[2].y+tv[3].y)*0.25,
                      .z = (tv[0].z+tv[1].z+tv[2].z+tv[3].z)*0.25 };
        queue.N = 0;
        int o = nearest_seed(seed_dt, real_seeds, adj, N, c);
        if (o >= 0) { stamp[o] = tet_stamp; if (!dynarray_push(&queue, &o)) goto err; }
        for (unsigned qh = 0; qh < queue.N; qh++) {
            int i; dynarray_get_value(&queue, qh, &i);
            const int *nbr = adj->neighbor + adj->offset[i];
            int NB = adj->offset[i+1] - adj->offset[i];
            if (!lp3_cell_maybe_hits_tet(tv, real_seeds->p[i], real_seeds, nbr, NB))
                continue;
            if (!cellkind[i]) {
                for (int f = 0; f < 4; f++) {
                    if (!face_is_domain_boundary(nc, f)) continue;   /* interior face */
                    int fv[3]; lp3_face_idx(f, fv);
                    s_point tri[3] = { tv[fv[0]], tv[fv[1]], tv[fv[2]] };
                    if (lp3_cell_maybe_hits_tri(tri, real_seeds->p[i], real_seeds, nbr, NB)) {
                        cellkind[i] = 1; break;
                    }
                }
            }
            for (int k = adj->offset[i]; k < adj->offset[i+1]; k++) {
                int j = adj->neighbor[k];
                if (j >= 0 && j < N && stamp[j] != tet_stamp) {
                    stamp[j] = tet_stamp;
                    if (!dynarray_push(&queue, &j)) goto err;
                }
            }
        }
    }

    /* Debug (NCVX_TETCHECK=1, combine with NCVX_ALL_BOUNDARY=1): per-tet audit,
     * the pieces carved out of one tet must sum to the tet volume. */
    int tetcheck = getenv("NCVX_TETCHECK") ? 1 : 0;
    double tetcheck_excess = 0.0;
    int    tetcheck_idx = -1;

    /* Phase B: main clip. Interior cells are skipped (emitted whole below) but
     * still expanded through, so the flood stays connected. */
    for (s_ncell *nc = domain->cdt.head; nc; nc = nc->next) {
        if (!nc->interior) continue;              /* skip exterior pockets (flagged CDT) */
        tet_stamp++;
        tetcheck_idx++;
        double tet_acc = 0.0;
        s_point tv[4];
        extract_vertices_ncell(&domain->cdt, nc, tv);

        s_point c = { .x = (tv[0].x+tv[1].x+tv[2].x+tv[3].x)*0.25,
                      .y = (tv[0].y+tv[1].y+tv[2].y+tv[3].y)*0.25,
                      .z = (tv[0].z+tv[1].z+tv[2].z+tv[3].z)*0.25 };
        queue.N = 0;
        int o = nearest_seed(seed_dt, real_seeds, adj, N, c);
        if (o >= 0) { stamp[o] = tet_stamp; if (!dynarray_push(&queue, &o)) goto err; }

        for (unsigned qh = 0; qh < queue.N; qh++) {
            int i; dynarray_get_value(&queue, qh, &i);
            const int *nbr = adj->neighbor + adj->offset[i];
            int NB = adj->offset[i+1] - adj->offset[i];
            if (!lp3_cell_maybe_hits_tet(tv, real_seeds->p[i], real_seeds, nbr, NB))
                continue;                                 /* frontier */

            int do_expand = 0;
            if (!cellkind[i]) {                           /* INTERIOR: emitted whole later */
                do_expand = 1;
            } else {                                      /* BOUNDARY: clip against this tet */
                int n_vtx  = adj->vtx_off[i+1]  - adj->vtx_off[i];
                int n_edge = adj->edge_off[i+1] - adj->edge_off[i];
                const int *vv = adj->vtx  + 3 * adj->vtx_off[i];
                const int *ve = adj->edge + 4 * adj->edge_off[i];
                int *rv = resolved ? resolved + adj->vtx_off[i] : NULL;
                s_convh piece; s_vlabel *labels = NULL;
                int r = clip_vcell_clip_tet_lp(i, real_seeds, adj, vv, n_vtx, ve, n_edge,
                                               tv, nc->vertex_id, owner, tie, rv, tet_stamp,
                                               EPS_DEG, &piece, want_surface ? &labels : NULL);
                if (r == 1) {
                    double pv = volume_convhull(&piece);
                    if (pv < EPS_DEG) { free_convhull(&piece); free(labels); }
                    else {
                        vol[i] += pv;
                        tet_acc += pv;
                        if (want_surface && !piece_to_surface(&piece, labels, nc, i, &surf[i])) {
                            free(labels); free_convhull(&piece); goto err;
                        }
                        free(labels);
                        if (!dynarray_push(&pcs[i], &piece)) { free_convhull(&piece); goto err; }
                        do_expand = 1;
                    }
                } else {
                    free(labels);  /* NULL on r != 1, but be explicit */
                }
            }
            if (do_expand)
                for (int k = adj->offset[i]; k < adj->offset[i+1]; k++) {
                    int j = adj->neighbor[k];
                    if (j >= 0 && j < N && stamp[j] != tet_stamp) {
                        stamp[j] = tet_stamp;
                        if (!dynarray_push(&queue, &j)) goto err;
                    }
                }
        }

        if (tetcheck) {
            double vt = fabs(signed_volume_tetra(tv));
            if (fabs(tet_acc - vt) > 1e-7) {
                int sg[4]; lp3_tet_face_orient(tv, sg);
                fprintf(stderr, "[ncvx] tet %d: pieces=%.9f tet=%.15e diff=%.9f "
                        "verts=(%d,%d,%d,%d) sigma=(%d,%d,%d,%d)\n",
                        tetcheck_idx, tet_acc, vt, tet_acc - vt,
                        nc->vertex_id[0], nc->vertex_id[1],
                        nc->vertex_id[2], nc->vertex_id[3],
                        sg[0], sg[1], sg[2], sg[3]);
                for (int j = 0; j < 4; j++)
                    fprintf(stderr, "[ncvx]   v%d = (%.17g, %.17g, %.17g)\n",
                            j, tv[j].x, tv[j].y, tv[j].z);
                tetcheck_excess += tet_acc - vt;
            }
        }
    }
    if (tetcheck)
        fprintf(stderr, "[ncvx] total per-tet excess: %.9f\n", tetcheck_excess);

    /* Interior cells: build the single whole-cell hull into their accumulator, so
     * the transfer loop below handles interior and boundary cells uniformly. */
    for (int i = 0; i < N; i++) {
        if (cellkind[i]) continue;
        int n_vtx = adj->vtx_off[i+1] - adj->vtx_off[i];
        const int *vv = adj->vtx + 3 * adj->vtx_off[i];
        s_convh hull; double hv; s_vlabel *labels = NULL;
        int r = emit_interior_cell(i, real_seeds, vv, n_vtx, EPS_DEG, &hull, &hv,
                                   want_surface ? &labels : NULL);
        if (r == -1) goto err;
        if (r == 1) {
            vol[i] = hv;
            if (want_surface && !piece_to_surface(&hull, labels, NULL, i, &surf[i])) {
                free(labels); free_convhull(&hull); goto err;
            }
            free(labels);
            if (!dynarray_push(&pcs[i], &hull)) { free_convhull(&hull); goto err; }
        } else {
            free(labels);
        }
    }

    /* Assemble per-cell vcells, transferring each accumulator's buffer directly
     * (no copy / no further allocation that could fail here). */
    for (int i = 0; i < N; i++) {
        vcells[i] = (s_vcell){0};
        vcells[i].seed_id  = i;
        vcells[i].N_pieces = (int)pcs[i].N;
        vcells[i].volume   = vol[i];
        if (pcs[i].N > 0) {
            vcells[i].pieces = (s_convh *)pcs[i].items;     /* transfer ownership */
        } else {
            dynarray_free(&pcs[i]);
            vcells[i].pieces = NULL;
        }

        if (want_surface) {   /* a pinched cell splits into >1 connected component */
            s_trimesh *sm = NULL; int sn = 0;
            if (build_cell_surface(&surf[i], EPS_DEG, &sm, &sn) == 1) {
                vcells[i].surface   = sm;
                vcells[i].N_surface = sn;
            }
            dynarray_free(&surf[i]);
        }
    }
    dynarray_free(&queue);
    free(owner); free(tie); free(resolved); free(cellkind);
    free(pcs); free(surf); free(vol); free(stamp);
    return 1;

err:
    for (int i = 0; i < N; i++) {
        for (unsigned q = 0; q < pcs[i].N; q++) {
            s_convh pp; dynarray_get_value(&pcs[i], q, &pp); free_convhull(&pp);
        }
        dynarray_free(&pcs[i]);
        if (surf) dynarray_free(&surf[i]);
    }
err_early:
    dynarray_free(&queue);
    free(owner); free(tie); free(resolved); free(cellkind);
    free(pcs); free(surf); free(vol); free(stamp);
    return 0;
}


/* ---------------------------------------------------------------------- */

s_ncvx_vdiagram vor3d_in_ncvx_domain(const s_points *seeds,
                                  const s_ncvx_domain *domain,
                                  double vol_max_rel_diff,
                                  double EPS_DEG, double TOL,
                                  int (*randint)(void*, int), void *rctx,
                                  s_dynarray *buff_points, int *out_kept_idx,
                                  int want_surface)
{
    (void)vol_max_rel_diff; (void)randint; (void)rctx; (void)buff_points;
    if (!seeds || seeds->N <= 0 || !ncvx_domain_is_valid(domain))
        return (s_ncvx_vdiagram){0};

    /* Seed Delaunay of the real seeds plus 8 FAR auxiliary points (a box many
     * times the domain size). They guarantee a non-degenerate 3D DT even for 1-3
     * (or coplanar) real seeds -- without which dt_builder_end leaves no real
     * tetrahedra and the diagram is empty. Being far outside every real seed's
     * circumsphere, they never flip a real Delaunay tet, so the all-real dual
     * enumeration the clip relies on is exactly the true real Voronoi structure
     * (unlike mirrors, which sit just outside the hull and corrupt it). Aux ids
     * (>= Nreal) are dropped by build_cell_adjacency; their far bisectors never
     * reach the domain, and the tet clip bounds every cell. */
    if (out_kept_idx)
        for (int i = 0; i < seeds->N; i++) out_kept_idx[i] = i;

    s_point mn = domain->bpoly.min, mx = domain->bpoly.max;
    s_point c  = { .x = 0.5*(mn.x+mx.x), .y = 0.5*(mn.y+mx.y), .z = 0.5*(mn.z+mx.z) };
    double  hx = 0.5*(mx.x-mn.x), hy = 0.5*(mx.y-mn.y), hz = 0.5*(mx.z-mn.z);
    double  h  = hx > hy ? (hx > hz ? hx : hz) : (hy > hz ? hy : hz);
    if (h <= 0) h = 1.0;
    const double FAR = 1000.0;                 /* far enough to clear any real circumsphere */
    s_point aux[8];
    for (int k = 0; k < 8; k++)
        aux[k] = (s_point){ .x = c.x + ((k&1)?FAR:-FAR)*h,
                            .y = c.y + ((k&2)?FAR:-FAR)*h,
                            .z = c.z + ((k&4)?FAR:-FAR)*h };
    s_point big_min = { .x = c.x-(FAR+2)*h, .y = c.y-(FAR+2)*h, .z = c.z-(FAR+2)*h };
    s_point big_max = { .x = c.x+(FAR+2)*h, .y = c.y+(FAR+2)*h, .z = c.z+(FAR+2)*h };

    s_dt_builder builder = dt_builder_begin(seeds, NULL, TOL, &big_min, &big_max);
    if (!builder._stack) {
        if (out_kept_idx) for (int i = 0; i < seeds->N; i++) out_kept_idx[i] = -1;
        return (s_ncvx_vdiagram){0};
    }
    s_points aux_pts = { .N = 8, .p = aux };
    if (!dt_builder_extend(&builder, &aux_pts, TOL)) {
        s_scplx tmp = dt_builder_end(&builder, false, NULL, NULL, 0);
        free_complex(&tmp);
        if (out_kept_idx) for (int i = 0; i < seeds->N; i++) out_kept_idx[i] = -1;
        return (s_ncvx_vdiagram){0};
    }
    int Nreal = seeds->N;                       /* first seeds->N inserted points are real */
    s_scplx seed_dt = dt_builder_end(&builder, false, &Nreal, out_kept_idx,
                                     out_kept_idx ? seeds->N : 0);
    if (!seed_dt.head || Nreal <= 0) { free_complex(&seed_dt); return (s_ncvx_vdiagram){0}; }

    int N = Nreal;
    s_points real_seeds = { .N = N, .p = seed_dt.points.p };  /* compacted: real seeds first */

    struct ncvx_cell_adjacency adj;
    if (!build_cell_adjacency(&seed_dt, N, &adj)) {
        free_complex(&seed_dt); return (s_ncvx_vdiagram){0};
    }
    if (!cell_adjacency_is_symmetric(&adj))
        fprintf(stderr, "vor3d_in_ncvx_domain: cell adjacency is not symmetric "
                        "(bug in build_cell_adjacency)\n");

    s_vcell *vcells = calloc((size_t)N, sizeof(s_vcell));
    if (!vcells) {
        free_cell_adjacency(&adj); free_complex(&seed_dt);
        return (s_ncvx_vdiagram){0};
    }

    /* Inverse (tet-outer) clip fills all vcells. */
    if (!clip_domain(domain, &real_seeds, &seed_dt, &adj, EPS_DEG, want_surface, vcells)) {
        free(vcells); free_cell_adjacency(&adj); free_complex(&seed_dt);
        return (s_ncvx_vdiagram){0};
    }

    /* Assemble result. Output seeds = the compacted real seeds; the output domain
     * carries surface/volume/bpoly only (empty CDT -- nobody reads it). */
    s_ncvx_vdiagram out;
    out.seeds  = copy_points(&real_seeds);
    out.domain = (s_ncvx_domain){
        .surface       = copy_trimesh(&domain->surface),
        .cdt           = (s_scplx){0},
        .domain_volume = domain->domain_volume,
        .bpoly         = bpoly_copy(&domain->bpoly),
        .spatial_index = NULL,
    };

    free_cell_adjacency(&adj);
    free_complex(&seed_dt);

    if (!out.seeds.p || !trimesh_is_valid(&out.domain.surface)) {
        free_points(&out.seeds);
        free_ncvx_domain(&out.domain);
        for (int i = 0; i < N; i++) free_vcell(&vcells[i]);
        free(vcells);
        return (s_ncvx_vdiagram){0};
    }

    out.vcells = vcells;
    return out;
}


s_ncvx_vdiagram vor3d_in_trimesh(const s_points *seeds,
                                   const s_trimesh *mesh,
                                   double vol_max_rel_diff,
                                   double EPS_DEG, double TOL,
                                   int (*randint)(void*, int), void *rctx,
                                   s_dynarray *buff_points, int *out_kept_idx,
                                   int want_surface)
{
    s_ncvx_domain domain = ncvx_domain_from_trimesh(mesh, EPS_DEG, TOL, 0);
    if (!ncvx_domain_is_valid(&domain)) return (s_ncvx_vdiagram){0};
    s_ncvx_vdiagram out = vor3d_in_ncvx_domain(seeds, &domain, vol_max_rel_diff,
                                           EPS_DEG, TOL, randint, rctx,
                                           buff_points, out_kept_idx, want_surface);
    free_ncvx_domain(&domain);  /* result owns its own copy */
    return out;
}
