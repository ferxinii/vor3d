/*
 * Alpha-shape surface extraction for (optionally weighted) point sets.
 *
 * Builds the (weighted) Delaunay triangulation, marks each tetrahedron in/out
 * by the squared power-radius test (r^2 <= alpha, the alpha complex tet test
 * acplx_tet_in), repairs the in-set so its boundary is edge-manifold, and
 * returns that outward-oriented boundary as a closed s_trimesh.  See
 * include/alpha.h and CLAUDE/ALPHA_SHAPE_PLAN.md.
 */

#include "alpha.h"
#include "acplx.h"
#include "scplx.h"
#include "delaunay.h"
#include "gtests.h"
#include "./hlibs/dynarray.h"
#include "./hlibs/hash.h"
#include <float.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


/* ------------------------------------------------------------------ */
/* Tet predicates                                                      */
/* ------------------------------------------------------------------ */

static inline bool tet_touches_sentinel(const s_ncell *nc)
{
    return nc->vertex_id[0] < 4 || nc->vertex_id[1] < 4 ||
           nc->vertex_id[2] < 4 || nc->vertex_id[3] < 4;
}

/* Squared power radius of a real tet: |c - p0|^2 - w0 with c the (power)
 * circumcentre (w0 = 0 and the plain circumcentre when unweighted).  May be
 * negative when weighted.  DBL_MAX if degenerate/flat. */
static double tet_power_radius2(const s_scplx *dt, const s_ncell *nc)
{
    s_point p[4]; extract_vertices_ncell(dt, nc, p);
    double  w[4]; extract_weights_ncell(dt, nc, w);
    s_point cc;
    if (!circumcentre_tetrahedron_w(p, dt->weights ? w : NULL, 1e-20, &cc)) return DBL_MAX;
    return dt->weights ? distance_squared(cc, p[0]) - w[0] : distance_squared(cc, p[0]);
}


/* ------------------------------------------------------------------ */
/* Marking + simple aggregates                                         */
/* ------------------------------------------------------------------ */

/* Set nc->mask_alpha for every tet by the alpha complex tet test; return
 * in-count. */
static int mark_by_alpha(s_scplx *dt, double alpha)
{
    int count = 0;
    for (s_ncell *nc = dt->head; nc; nc = nc->next) {
        nc->mask_alpha = acplx_tet_in(dt, nc, alpha, /*has_big_tetra=*/true);
        if (nc->mask_alpha) count++;
    }
    return count;
}

/* Total volume of all in-tets (mask_alpha). */
static double in_set_volume(const s_scplx *dt)
{
    double vol = 0.0;
    for (s_ncell *nc = dt->head; nc; nc = nc->next) {
        if (!nc->mask_alpha) continue;
        s_point p[4]; extract_vertices_ncell(dt, nc, p);
        vol += fabs(signed_volume_tetra(p));
    }
    return vol;
}


/* ------------------------------------------------------------------ */
/* Hash helpers                                                        */
/* ------------------------------------------------------------------ */

static size_t hash_ptr(const void *key)
{
    uintptr_t p = *(const uintptr_t *)key;
    return (size_t)((p >> 4) * 0x9E3779B97F4A7C15ULL);
}
static bool eq_ptr(const void *a, const void *b)
{
    return *(const uintptr_t *)a == *(const uintptr_t *)b;
}
static size_t hash_edge2(const void *key)
{
    const int *v = key;
    uint64_t h = ((uint64_t)v[0] * 73856093u) ^ ((uint64_t)v[1] * 19349663u);
    h ^= h >> 33; h *= 0xff51afd7ed558ccdULL; h ^= h >> 33;
    return (size_t)h;
}
static bool eq_edge2(const void *a, const void *b)
{
    const int *x = a, *y = b; return x[0] == y[0] && x[1] == y[1];
}


/* ------------------------------------------------------------------ */
/* Manifold repair (Phase 3)                                           */
/* ------------------------------------------------------------------ */

static bool AUX_collect_ring(void *C, const s_ncell *nc)
{
    s_dynarray *ring = C;
    s_ncell *p = (s_ncell *)nc;
    dynarray_push(ring, &p);
    return true;   /* collect the whole ring */
}

/* Make the boundary edge-manifold around edge (la,lb) of tet `start`.
 *
 * The tets around a DT edge form a cyclic ring (a full cycle because the big
 * tetra is kept, so every face has a neighbor).  The boundary of the in-set
 * has 2 faces at this edge per maximal run of in-tets; more than one in-run
 * (equivalently, more than one out-"gap") makes the edge non-manifold.  We
 * dilate: fill every out-gap except one, promoting real out-tets to in.  A gap
 * that contains a sentinel or a locked tet is "blocked" and cannot be filled;
 * if two or more blocked gaps remain, fall back to demoting (and locking) all
 * but the largest in-run so the blocked gaps merge.  Returns #tets changed. */
static int repair_edge(s_ncell *start, int la, int lb, s_hash_table *locked)
{
    int omitted[2];
    for (int k = 0, i = 0; i < 4; i++) if (i != la && i != lb) omitted[k++] = i;

    s_dynarray ring = dynarray_initialize(sizeof(s_ncell *), 16);
    if (!ring.items) return 0;
    walk_ridge_cycle_and_check_ncells(start, omitted, AUX_collect_ring, &ring);
    int R = (int)ring.N;
    s_ncell **rc = ring.items;
    if (R < 3) { dynarray_free(&ring); return 0; }

    bool *isin = malloc(R * sizeof(bool));
    bool *blk  = malloc(R * sizeof(bool));
    if (!isin || !blk) { free(isin); free(blk); dynarray_free(&ring); return 0; }

    int nin = 0;
    for (int i = 0; i < R; i++) {
        isin[i] = rc[i]->mask_alpha;
        if (isin[i]) nin++;
        blk[i] = tet_touches_sentinel(rc[i]) || hash_get(locked, &rc[i]);
    }

    int changed = 0;
    if (nin == 0 || nin == R) goto cleanup;   /* no boundary / all in */

    int s0 = -1; for (int i = 0; i < R; i++) if (isin[i]) { s0 = i; break; }

    /* Enumerate gaps (maximal out-runs) cyclically from an in-tet. */
    int  *gidx   = malloc(R * sizeof(int));
    int  *gstart = malloc(R * sizeof(int));
    int  *glen   = malloc(R * sizeof(int));
    bool *gblk   = malloc(R * sizeof(bool));
    if (!gidx || !gstart || !glen || !gblk) {
        free(gidx); free(gstart); free(glen); free(gblk); goto cleanup;
    }
    int ng = 0, gp = 0; bool prevout = false;
    for (int j = 0; j < R; j++) {
        int idx = (s0 + j) % R;
        if (isin[idx]) { prevout = false; continue; }
        if (!prevout) { gstart[ng] = gp; glen[ng] = 0; gblk[ng] = false; ng++; }
        gidx[gp++] = idx; glen[ng - 1]++;
        if (blk[idx]) gblk[ng - 1] = true;
        prevout = true;
    }

    if (ng > 1) {
        int nblk = 0, largest = 0, largestlen = -1;
        for (int g = 0; g < ng; g++) {
            if (gblk[g]) nblk++;
            if (glen[g] > largestlen) { largestlen = glen[g]; largest = g; }
        }

        if (nblk == 0) {
            /* Keep the largest gap open; fill the others. */
            for (int g = 0; g < ng; g++) {
                if (g == largest) continue;
                for (int t = 0; t < glen[g]; t++) {
                    int idx = gidx[gstart[g] + t];
                    if (!rc[idx]->mask_alpha) { rc[idx]->mask_alpha = true; changed++; }
                }
            }
        } else {
            /* Fill every unblocked gap; keep blocked gaps open. */
            for (int g = 0; g < ng; g++) {
                if (gblk[g]) continue;
                for (int t = 0; t < glen[g]; t++) {
                    int idx = gidx[gstart[g] + t];
                    if (!rc[idx]->mask_alpha) { rc[idx]->mask_alpha = true; changed++; }
                }
            }
            if (nblk >= 2) {
                /* Demotion fallback: after filling, the in-runs are the segments
                 * between the (still-open) blocked gaps.  Keep the largest such
                 * run and demote + lock the rest so the blocked gaps merge. */
                int b0 = -1;
                for (int i = 0; i < R; i++) if (!rc[i]->mask_alpha && blk[i]) { b0 = i; break; }
                if (b0 >= 0) {
                    int *runidx = malloc(R * sizeof(int));
                    int *rstart = malloc(R * sizeof(int));
                    int *rlen   = malloc(R * sizeof(int));
                    if (runidx && rstart && rlen) {
                        int nr = 0, rp = 0; bool previn = false;
                        for (int j = 1; j <= R; j++) {
                            int idx = (b0 + j) % R;
                            if (rc[idx]->mask_alpha) {
                                if (!previn) { rstart[nr] = rp; rlen[nr] = 0; nr++; }
                                runidx[rp++] = idx; rlen[nr - 1]++; previn = true;
                            } else previn = false;
                        }
                        int keep = 0, keeplen = -1;
                        for (int r = 0; r < nr; r++) if (rlen[r] > keeplen) { keeplen = rlen[r]; keep = r; }
                        for (int r = 0; r < nr; r++) {
                            if (r == keep) continue;
                            for (int t = 0; t < rlen[r]; t++) {
                                int idx = runidx[rstart[r] + t];
                                if (rc[idx]->mask_alpha) { rc[idx]->mask_alpha = false; changed++; }
                                bool tv = true; hash_insert(locked, &rc[idx], &tv);
                            }
                        }
                    }
                    free(runidx); free(rstart); free(rlen);
                }
            }
        }
    }

    free(gidx); free(gstart); free(glen); free(gblk);

cleanup:
    free(isin); free(blk);
    dynarray_free(&ring);
    return changed;
}

/* Sweep all in-tet edges repairing non-manifold ones until stable.
 * Returns total number of tet state changes. */
static int repair_manifold(s_scplx *dt)
{
    s_hash_table locked;
    if (!hash_init(&locked, sizeof(s_ncell *), sizeof(bool), 256, 0,
                   hash_ptr, eq_ptr, NULL))
        return 0;

    int total = 0, sweep = 0, cap = dt->N_ncells + 16;
    bool changed = true;
    while (changed && sweep < cap) {
        changed = false; sweep++;
        s_hash_table seen;
        if (!hash_init(&seen, sizeof(int) * 2, sizeof(bool),
                       (size_t)dt->N_ncells * 3 + 16, 0,
                       hash_edge2, eq_edge2, NULL))
            break;
        for (s_ncell *nc = dt->head; nc; nc = nc->next) {
            if (!nc->mask_alpha) continue;
            for (int a = 0; a < 3; a++) for (int b = a + 1; b < 4; b++) {
                int e[2] = { nc->vertex_id[a], nc->vertex_id[b] };
                if (e[0] > e[1]) { int t = e[0]; e[0] = e[1]; e[1] = t; }
                if (hash_get(&seen, e)) continue;
                bool tv = true; hash_insert(&seen, e, &tv);
                int ch = repair_edge(nc, a, b, &locked);
                if (ch) { changed = true; total += ch; }
            }
        }
        hash_free(&seen);
    }
    hash_free(&locked);
    if (sweep >= cap)
        fprintf(stderr, "ashape: manifold repair hit sweep cap (%d)\n", cap);
    return total;
}


/* ------------------------------------------------------------------ */
/* Connected components + coverage diagnostics                         */
/* ------------------------------------------------------------------ */

static int count_components(s_scplx *dt)
{
    s_hash_table seen;   /* ./hlibs/hash.h tables never grow: size for one entry per tet */
    if (!hash_init(&seen, sizeof(s_ncell *), sizeof(bool),
                   (size_t)dt->N_ncells + 16, 0, hash_ptr, eq_ptr, NULL))
        return -1;
    s_dynarray stack = dynarray_initialize(sizeof(s_ncell *), 64);
    int comps = 0;
    for (s_ncell *root = dt->head; root; root = root->next) {
        if (!root->mask_alpha || hash_get(&seen, &root)) continue;
        comps++;
        bool tv = true; hash_insert(&seen, &root, &tv);
        dynarray_push(&stack, &root);
        while (stack.N > 0) {
            s_ncell *nc; dynarray_pop(&stack, &nc);
            for (int i = 0; i < 4; i++) {
                s_ncell *o = nc->opposite[i];
                if (o && o->mask_alpha && !hash_get(&seen, &o)) {
                    hash_insert(&seen, &o, &tv);
                    dynarray_push(&stack, &o);
                }
            }
        }
    }
    dynarray_free(&stack);
    hash_free(&seen);
    return comps;
}

/* Input points that are not a vertex of any in-tet: kept[i] is the DT vid of
 * input i, or -1 (near-duplicate / redundant weighted point -> dropped too). */
static int count_dropped_points(const s_scplx *dt, const int *kept, int N_in)
{
    bool *used = calloc(dt->points.N, sizeof(bool));
    if (!used) return -1;
    for (s_ncell *nc = dt->head; nc; nc = nc->next)
        if (nc->mask_alpha)
            for (int i = 0; i < 4; i++) used[nc->vertex_id[i]] = true;
    int dropped = 0;
    for (int i = 0; i < N_in; i++) if (kept[i] < 0 || !used[kept[i]]) dropped++;
    free(used);
    return dropped;
}


/* ------------------------------------------------------------------ */
/* Auto-alpha (Phase 4)                                                */
/* ------------------------------------------------------------------ */

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/* True iff, at threshold alpha, every point that is a vertex of some real tet
 * is also a vertex of at least one in-tet.  r2cache is aligned with the tet
 * list iteration order (DBL_MAX for sentinel/degenerate tets); comparing the
 * cached circumradii keeps the binary search consistent with the thresholds,
 * which are those same cached values -- and avoids re-running the orthosphere
 * predicate over every tet at every search step (on lattice-like inputs the
 * threshold IS a tet circumradius, so the predicate ties exactly and falls
 * through to exact arithmetic each time). */
static bool covers_all(s_scplx *dt, double alpha, const double *r2cache,
                       const bool *coverable, bool *cov, int Npts)
{
    for (int i = 0; i < Npts; i++) cov[i] = false;
    int idx = 0;
    for (s_ncell *nc = dt->head; nc; nc = nc->next, idx++) {
        if (r2cache[idx] > alpha) continue;
        for (int i = 0; i < 4; i++) cov[nc->vertex_id[i]] = true;
    }
    for (int i = 0; i < Npts; i++) if (coverable[i] && !cov[i]) return false;
    return true;
}

/* Smallest alpha (from the sorted real squared power radii) covering all
 * coverable points.  May be <= 0 when weighted.  Returns NAN if there is no
 * real non-degenerate tet or on allocation failure. */
static double auto_alpha(s_scplx *dt)
{
    int Npts = dt->points.N;
    bool *coverable = calloc(Npts, sizeof(bool));
    bool *cov = malloc(Npts * sizeof(bool));
    double *r2cache = malloc((size_t)dt->N_ncells * sizeof(double));
    s_dynarray vals = dynarray_initialize(sizeof(double), 256);
    if (!coverable || !cov || !r2cache || !vals.items) {
        free(coverable); free(cov); free(r2cache); dynarray_free(&vals);
        return NAN;
    }

    int idx = 0;
    for (s_ncell *nc = dt->head; nc; nc = nc->next, idx++) {
        r2cache[idx] = DBL_MAX;
        if (tet_touches_sentinel(nc)) continue;
        for (int i = 0; i < 4; i++) coverable[nc->vertex_id[i]] = true;
        double r2 = tet_power_radius2(dt, nc);
        if (r2 < DBL_MAX) { dynarray_push(&vals, &r2); r2cache[idx] = r2; }
    }
    if (vals.N == 0) {
        free(coverable); free(cov); free(r2cache); dynarray_free(&vals);
        return NAN;
    }

    qsort(vals.items, vals.N, sizeof(double), cmp_double);
    double *v = vals.items;

    int lo = 0, hi = (int)vals.N - 1, ans = (int)vals.N - 1;
    while (lo <= hi) {
        int mid = (lo + hi) / 2;
        if (covers_all(dt, v[mid], r2cache, coverable, cov, Npts)) { ans = mid; hi = mid - 1; }
        else lo = mid + 1;
    }
    /* Nudge up by a relative 1e-9 so the exact predicate includes the tet
     * whose radius set the threshold (towards +inf also when weighted and
     * negative). */
    double alpha = v[ans] >= 0 ? v[ans] * (1.0 + 1e-9) : v[ans] * (1.0 - 1e-9);

    free(coverable); free(cov); free(r2cache); dynarray_free(&vals);
    return alpha;
}


/* ------------------------------------------------------------------ */
/* Boundary extraction + trimesh assembly                              */
/* ------------------------------------------------------------------ */

/* old2new_out (length dt->points.N, nullable): DT vid -> surface vertex
 * index, or -1 (left untouched on failure). */
static s_trimesh boundary_to_trimesh(s_scplx *dt, int *old2new_out)
{
    s_dynarray faces = dynarray_initialize(sizeof(int) * 3, 64);
    if (!faces.items) return trimesh_NAN;

    for (s_ncell *nc = dt->head; nc; nc = nc->next) {
        if (!nc->mask_alpha) continue;
        for (int i = 0; i < 4; i++) {
            s_ncell *opp = nc->opposite[i];
            if (opp && opp->mask_alpha) continue;   /* interior face, skip */

            int fv[3];
            for (int k = 0, j = 0; j < 4; j++)
                if (j != i) fv[k++] = nc->vertex_id[j];

            /* Wind outward: normal must point away from the apex (vertex i). */
            s_point a = dt->points.p[fv[0]];
            s_point b = dt->points.p[fv[1]];
            s_point c = dt->points.p[fv[2]];
            s_point apex = dt->points.p[nc->vertex_id[i]];
            s_point n = cross_prod(subtract_points(b, a), subtract_points(c, a));
            if (dot_prod(n, subtract_points(apex, a)) > 0.0) {
                int tmp = fv[1]; fv[1] = fv[2]; fv[2] = tmp;
            }
            if (!dynarray_push(&faces, fv)) { dynarray_free(&faces); return trimesh_NAN; }
        }
    }

    if (faces.N == 0) { dynarray_free(&faces); return trimesh_NAN; }

    int Npts = dt->points.N;
    int *old2new = malloc(sizeof(int) * Npts);
    if (!old2new) { dynarray_free(&faces); return trimesh_NAN; }
    for (int i = 0; i < Npts; i++) old2new[i] = -1;

    s_point *verts = malloc(sizeof(s_point) * Npts);   /* upper bound */
    int *tris = malloc(sizeof(int) * 3 * faces.N);
    if (!verts || !tris) { free(old2new); free(verts); free(tris);
                           dynarray_free(&faces); return trimesh_NAN; }

    int Nv = 0;
    for (unsigned f = 0; f < faces.N; f++) {
        int *fv = dynarray_get_ptr(&faces, f);
        for (int k = 0; k < 3; k++) {
            int gid = fv[k];
            if (old2new[gid] < 0) { old2new[gid] = Nv; verts[Nv] = dt->points.p[gid]; Nv++; }
            tris[f * 3 + k] = old2new[gid];
        }
    }

    s_trimesh out = trimesh_from_arrays(verts, Nv, tris, (int)faces.N, 0.0);
    if (old2new_out && out.points.N > 0)   /* points are copied verbatim */
        memcpy(old2new_out, old2new, sizeof(int) * Npts);

    free(old2new); free(verts); free(tris);
    dynarray_free(&faces);
    return out;
}


/* ------------------------------------------------------------------ */
/* Public entry point                                                  */
/* ------------------------------------------------------------------ */

bool alpha_shape_is_valid(const s_alpha_shape *as)
{
    return as != NULL && !isnan(as->alpha);
}

void free_alpha_shape(s_alpha_shape *as)
{
    if (!as) return;
    free_trimesh(&as->surface);
    *as = (s_alpha_shape){ .surface = trimesh_NAN, .alpha = NAN };
}

s_alpha_shape alpha_shape_3d(const s_points *pts, const double *weights,
                             double alpha, double TOL_dup,
                             s_random_context *rng, int *out_kept_idx)
{
    const s_alpha_shape FAIL = { .surface = trimesh_NAN, .alpha = NAN };
    int N_in = pts ? pts->N : 0;
    if (out_kept_idx) for (int i = 0; i < N_in; i++) out_kept_idx[i] = -1;
    if (!pts || pts->N < 4) return FAIL;

    struct timespec _at; clock_gettime(CLOCK_MONOTONIC, &_at);
    double _t0 = _at.tv_sec + 1e-9 * _at.tv_nsec, _t1;
#define ASHAPE_TICK(msg) do { \
        clock_gettime(CLOCK_MONOTONIC, &_at); \
        _t1 = _at.tv_sec + 1e-9 * _at.tv_nsec; \
        if (getenv("ASHAPE_PROFILE")) \
            fprintf(stderr, "[ashape][t] %-10s %.3fs\n", msg, _t1 - _t0); \
        _t0 = _t1; \
    } while (0)

    /* DT with the sentinels kept; kept[i] = DT vid of input i, or -1. */
    int *kept = malloc(pts->N * sizeof(int));
    if (!kept) return FAIL;
    for (int i = 0; i < pts->N; i++) kept[i] = i;
    s_dt_builder b = dt_builder_begin(pts, weights, TOL_dup, NULL, NULL, rng);
    if (!b._stack) { free(kept); return FAIL; }
    s_scplx dt = dt_builder_end(&b, /*keep_big_tetra=*/true, NULL, kept, pts->N);
    ASHAPE_TICK("dt");

    s_alpha_shape as = FAIL;
    int *old2new = NULL;

    double alpha_used = alpha;
    if (isnan(alpha_used)) {
        alpha_used = auto_alpha(&dt);   /* NAN: degenerate (no real tet) / alloc */
        if (isnan(alpha_used)) goto done;
    }
    ASHAPE_TICK("auto_alpha");

    if (mark_by_alpha(&dt, alpha_used) == 0) goto done;   /* empty / degenerate */
    ASHAPE_TICK("mark");

    int promoted = repair_manifold(&dt);
    ASHAPE_TICK("repair");

    old2new = malloc(dt.points.N * sizeof(int));
    if (!old2new) goto done;
    s_trimesh mesh = boundary_to_trimesh(&dt, old2new);
    if (mesh.points.N == 0) goto done;
    ASHAPE_TICK("boundary");

    int comps = count_components(&dt);
    ASHAPE_TICK("components");
    int dropped = count_dropped_points(&dt, kept, pts->N);
    ASHAPE_TICK("dropped");
    if (comps < 0 || dropped < 0) { free_trimesh(&mesh); goto done; }

    as = (s_alpha_shape){
        .surface       = mesh,
        .alpha         = alpha_used,
        .volume        = in_set_volume(&dt),
        .N_components  = comps,
        .N_promoted    = promoted,
        .N_pts_dropped = dropped,
    };
    ASHAPE_TICK("volume");

    if (out_kept_idx)
        for (int i = 0; i < pts->N; i++)
            out_kept_idx[i] = kept[i] >= 0 ? old2new[kept[i]] : -1;
#undef ASHAPE_TICK

done:
    free(old2new);
    free(kept);
    free_complex(&dt);
    return as;
}
