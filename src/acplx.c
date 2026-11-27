#include "acplx.h"
#include "alpha.h"
#include "delaunay.h"
#include "scplx.h"
#include "./hlibs/hash.h"
#include "./hlibs/dynarray.h"
#include "gtests.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


/* Test for every simplex: return true if simplex satisfies radius and conflict conditions */
bool acplx_tet_in(const s_scplx *scplx, const s_ncell *nc, double alpha, bool has_big_tetra)
{   /* Always conflict free, so only check radius condition */
    if (has_big_tetra && (nc->vertex_id[0] < 4 || nc->vertex_id[1] < 4 ||
                          nc->vertex_id[2] < 4 || nc->vertex_id[3] < 4 )) {
        return false;
    }

    s_point p[4];  extract_vertices_ncell(scplx, nc, p);
    double  w[4];  extract_weights_ncell(scplx, nc, w);
    return test_orthosphere_w(4, p, w, alpha) <= 0;
}


static bool face_test(const s_scplx *scplx, const s_ncell *nc, int face_localid,
                      double alpha, bool has_big_tetra)
{
    if (has_big_tetra) {
        int fid[3]; extract_ids_face(nc, 2, &face_localid, fid);
        if (fid[0] < 4 || fid[1] < 4 || fid[2] < 4) return false;
    }

    /* 1) Radius condition. */
    s_point p[3]; double w[3];
    extract_vertices_and_weights_face(scplx, nc, 2, &face_localid, p, w);
    if (test_orthosphere_w(3, p, w, alpha) > 0) return false;

    /* 2) Conflict condition. */
    /* Vertex opposite in nc */
    int vid1 = nc->vertex_id[face_localid];
    if (!has_big_tetra || vid1 >= 4)   
        if (test_orthosphere(3, p, w, scplx->points.p[vid1], 
                             scplx->weights ? scplx->weights[vid1] : 0.0) > 0)
            return false;

    /* Vertex opposite in opp tet */
    s_ncell *opp = nc->opposite[face_localid];
    if (opp) {
        int opp_local;
        face_localid_of_adjacent_ncell(nc, 2, &face_localid, face_localid, &opp_local);
        int vid2 = opp->vertex_id[opp_local];
        if (!has_big_tetra || vid2 >= 4)
            if (test_orthosphere(3, p, w, scplx->points.p[vid2], 
                                 scplx->weights ? scplx->weights[vid2] : 0.0) > 0)
                return false;
    }

    return true;
}


typedef struct {
    s_scplx *scplx;
    int gid0, gid1;
    bool has_big_tetra;
    s_point edge_p[2];
    double edge_w[2];
    bool empty_ball_ok;   /* output */
} edge_test_ctx;

static bool AUX_check_edge_acplx(void *C, const s_ncell *nc)
{   /* Checks conflict condition */
    edge_test_ctx *ctx = C;
    if (ctx->empty_ball_ok) {
        for (int v = 0; v < 4; ++v) {
            int vid = nc->vertex_id[v];
            if (ctx->has_big_tetra && vid < 4) continue;
            if (vid == ctx->gid0 || vid == ctx->gid1) continue;
            
            s_point pv = ctx->scplx->points.p[vid];
            double  wv  = ctx->scplx->weights ? ctx->scplx->weights[vid] : 0.0;
            if (test_orthosphere(2, ctx->edge_p, ctx->edge_w, pv, wv) > 0) {
                ctx->empty_ball_ok = false;
                return false;
            }
        }
    }
    return true;
}

static bool edge_test(s_scplx *scplx, const s_ncell *nc, const int edge_endpoints[2],
                      double alpha, bool has_big_tetra)
{
    int gid0 = nc->vertex_id[edge_endpoints[0]];
    int gid1 = nc->vertex_id[edge_endpoints[1]];
    if (has_big_tetra && (gid0 < 4 || gid1 < 4)) return false;

    int omitted[2];  /* The other two local ids */
    for (int k=0, i=0; i<4; i++) if (i != edge_endpoints[0] && i != edge_endpoints[1])
        omitted[k++] = i;

    s_point p[2]; double w[2];
    extract_vertices_and_weights_face(scplx, nc, 1, omitted, p, w);

    /* 1) Radius condition. */
    if (test_orthosphere_w(2, p, w, alpha) > 0) return false;

    /* 2) Power condition by walking ridge. */
    edge_test_ctx ctx = {
        .scplx = scplx,
        .gid0 = gid0, .gid1 = gid1,
        .has_big_tetra = has_big_tetra,
        .edge_p[0] = p[0], .edge_p[1] = p[1],
        .edge_w[0] = w[0], .edge_w[1] = w[1],
        .empty_ball_ok = true,
    };
    walk_ridge_cycle_and_check_ncells(nc, omitted, AUX_check_edge_acplx, &ctx);
    return ctx.empty_ball_ok;
}


static bool vertex_test(s_scplx *scplx, s_ncell *nc, int v_localid, double alpha,
                        bool has_big_tetra, s_dynarray *buff_ncellPTR)
{
    int vid = nc->vertex_id[v_localid];
    if (has_big_tetra && (vid < 4)) return false;

    /* 1) Radius condition */
    s_point pi_pt = scplx->points.p[vid];
    double  pi_wt = scplx->weights ? scplx->weights[vid] : 0.0;
    if (test_orthosphere_w(1, &pi_pt, &pi_wt, alpha) > 0) return false;

    /* 2) Conflict ("attached") condition (Finding C).  A vertex can be
     * non-redundant in the RT (nonempty power cell) yet be buried inside a
     * bigger neighbouring ball, contributing nothing to the union; such a
     * vertex must be excluded from K.  Vertex a is "attached" to neighbour b iff
     * dist(a,b)^2 + wa - wb < 0, which is exactly test_orthosphere(1, {a}, b) > 0
     * (b in conflict with a's vertex-orthosphere y_a = (a, -wa)).  Walk the star
     * of the vertex and test every neighbouring vertex.
     *
     * ncells_incident_face(dim_face=0) identifies the vertex by the LOCAL ids of
     * the OTHER three tet vertices (it flood-fills across faces opposite those),
     * so pass {0,1,2,3} \ {v_localid} -- NOT extract_ids_face, which returns
     * GLOBAL ids and would break the star walk. */
    int other_localids[3];
    for (int i = 0, k = 0; i < 4; i++)
        if (i != v_localid) other_localids[k++] = i;
    ncells_incident_face(scplx, nc, 0, other_localids, buff_ncellPTR);

    for (unsigned i = 0; i < buff_ncellPTR->N; i++) {
        s_ncell **tc = dynarray_get_ptr(buff_ncellPTR, i);
        for (int j = 0; j < 4; j++) {
            int vjd = (*tc)->vertex_id[j];
            if (vjd == vid) continue;
            if (has_big_tetra && vjd < 4) continue;
            s_point pj = scplx->points.p[vjd];
            double  wj = scplx->weights ? scplx->weights[vjd] : 0.0;
            if (test_orthosphere(1, &pi_pt, &pi_wt, pj, wj) > 0) return false;
        }
    }

    return true;
}



/* Main algorithm to extract alpha complex */

typedef struct vertex_mark {
    bool seen;
    bool in_acplx;
} s_vertex_mark;

static inline void sort4(int *a)
{
    if (a[0] > a[1]) { int t=a[0]; a[0]=a[1]; a[1]=t; }
    if (a[2] > a[3]) { int t=a[2]; a[2]=a[3]; a[3]=t; }
    if (a[0] > a[2]) { int t=a[0]; a[0]=a[2]; a[2]=t; }
    if (a[1] > a[3]) { int t=a[1]; a[1]=a[3]; a[3]=t; }
    if (a[1] > a[2]) { int t=a[1]; a[1]=a[2]; a[2]=t; }
}

static inline void sort3(int *a)
{
    if (a[0] > a[1]) { int t=a[0]; a[0]=a[1]; a[1]=t; }
    if (a[1] > a[2]) { int t=a[1]; a[1]=a[2]; a[2]=t; }
    if (a[0] > a[1]) { int t=a[0]; a[0]=a[1]; a[1]=t; }
}

static inline void sort2(int *a)
{
    if (a[0] > a[1]) { int t=a[0]; a[0]=a[1]; a[1]=t; }
}

static bool hash_eq_faces(const void *k1, const void *k2)
{
    const int *a = k1, *b = k2;
    return (a[0]==b[0] && a[1]==b[1] && a[2]==b[2]);
}

static bool hash_eq_edges(const void *k1, const void *k2)
{
    const int *a = k1, *b = k2;
    return (a[0]==b[0] && a[1]==b[1]);
}

static inline uint64_t mix64(uint64_t x)
{
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
}

static size_t hash_fun_faces(const void *key)
{
    const int *v = key;
    uint64_t h = ((uint64_t)v[0] * 73856093u) ^
                 ((uint64_t)v[1] * 19349663u) ^
                 ((uint64_t)v[2] * 83492791u);
    return (size_t)mix64(h);
}

static size_t hash_fun_edges(const void *key)
{
    const int *v = key;
    uint64_t h = ((uint64_t)v[0] * 73856093u) ^
                 ((uint64_t)v[1] * 19349663u);
    return (size_t)mix64(h);
}


static int mark_recursive_edge(int edge[2], s_hash_table *ht_edges, 
                               s_vertex_mark *vmark)
{
    /* Mark edge */
    sort2(edge);  
    bool t = true;
    if (!hash_insert(ht_edges, edge, &t)) return 0;

    vmark[edge[0]].in_acplx = true;
    vmark[edge[0]].seen = true;
    vmark[edge[1]].in_acplx = true;
    vmark[edge[1]].seen = true;

    return 1;
}


static int mark_recursive_face(int face[3], s_hash_table *ht_faces, 
                               s_hash_table *ht_edges, s_vertex_mark *vmark)
{
    /* Mark face */
    sort3(face);
    bool t = true;
    if (!hash_insert(ht_faces, face, &t)) return 0;

    /* Mark edges */
    if (!hash_insert(ht_edges, (int[2]){face[0], face[1]}, &t)) return 0;
    if (!hash_insert(ht_edges, (int[2]){face[0], face[2]}, &t)) return 0;
    if (!hash_insert(ht_edges, (int[2]){face[1], face[2]}, &t)) return 0;

    /* Mark vertices */
    vmark[face[0]].in_acplx = true;
    vmark[face[0]].seen = true;
    vmark[face[1]].in_acplx = true;
    vmark[face[1]].seen = true;
    vmark[face[2]].in_acplx = true;
    vmark[face[2]].seen = true;

    return 1;
}


static int mark_recursive_ncell(s_ncell *nc, s_hash_table *ht_faces, 
                                s_hash_table *ht_edges, s_vertex_mark *vmark)
{
    nc->mask_alpha = true;

    int tet[4] = {nc->vertex_id[0], nc->vertex_id[1], 
                  nc->vertex_id[2], nc->vertex_id[3]};
    sort4(tet);

    /* Mark faces */
    bool t = true;
    if (!hash_insert(ht_faces, (int[3]){tet[0], tet[1], tet[2]}, &t)) return 0;
    if (!hash_insert(ht_faces, (int[3]){tet[0], tet[1], tet[3]}, &t)) return 0;
    if (!hash_insert(ht_faces, (int[3]){tet[0], tet[2], tet[3]}, &t)) return 0;
    if (!hash_insert(ht_faces, (int[3]){tet[1], tet[2], tet[3]}, &t)) return 0;

    /* Mark edges */
    if (!hash_insert(ht_edges, (int[2]){tet[0], tet[1]}, &t)) return 0;
    if (!hash_insert(ht_edges, (int[2]){tet[0], tet[2]}, &t)) return 0;
    if (!hash_insert(ht_edges, (int[2]){tet[0], tet[3]}, &t)) return 0;
    if (!hash_insert(ht_edges, (int[2]){tet[1], tet[2]}, &t)) return 0;
    if (!hash_insert(ht_edges, (int[2]){tet[1], tet[3]}, &t)) return 0;
    if (!hash_insert(ht_edges, (int[2]){tet[2], tet[3]}, &t)) return 0;

    /* Mark vertices */
    vmark[tet[0]].in_acplx = true;
    vmark[tet[0]].seen = true;
    vmark[tet[1]].in_acplx = true;
    vmark[tet[1]].seen = true;
    vmark[tet[2]].in_acplx = true;
    vmark[tet[2]].seen = true;
    vmark[tet[3]].in_acplx = true;
    vmark[tet[3]].seen = true;

    return 1;
}


int alpha_complex_mark(s_scplx *scplx, bool has_big_tetra, double alpha, s_dynarray *buff_ncellPTR,
                       s_hash_table *out_faces, s_hash_table *out_edges,
                       bool *out_vertices)
{
    int N_tetra = scplx->N_ncells;
    s_hash_table ht_faces, ht_edges;
    s_vertex_mark *vmark = NULL;
    if (!hash_init(&ht_faces, sizeof(int)*3, sizeof(bool), N_tetra*4,
                   N_tetra*4, hash_fun_faces, hash_eq_faces, NULL)) return 0;
    if (!hash_init(&ht_edges, sizeof(int)*2, sizeof(bool), N_tetra*6,
                   N_tetra*6, hash_fun_edges, hash_eq_edges, NULL)) {
        hash_free(&ht_faces);
        return 0;
    }
    vmark = calloc(scplx->points.N, sizeof(s_vertex_mark));
    if (!vmark) goto fail;
    bool f = false;

    /* Detect ncells */
    for (s_ncell *nc = scplx->head; nc; nc = nc->next) {
        if (acplx_tet_in(scplx, nc, alpha, has_big_tetra))
            if (!mark_recursive_ncell(nc, &ht_faces, &ht_edges, vmark)) goto fail;
    }

    /* Detect faces */
    for (s_ncell *nc = scplx->head; nc; nc = nc->next) for (int i=0; i<4; i++) {
        int face_vid[3]; extract_ids_face(nc, 2, &i, face_vid); sort3(face_vid);
        if (!hash_get(&ht_faces, face_vid)) { 
            if (face_test(scplx, nc, i, alpha, has_big_tetra)) {
                if (!mark_recursive_face(face_vid, &ht_faces, &ht_edges, vmark)) goto fail;
            } else {
                if (!hash_insert(&ht_faces, face_vid, &f)) goto fail;
            }
        }
    }

    /* Detect edges */
    for (s_ncell *nc = scplx->head; nc; nc = nc->next) {
        for (int i=0; i<3; i++) for (int j=i+1; j<4; j++) {
            int edge_vid[2] = {nc->vertex_id[i], nc->vertex_id[j]}; sort2(edge_vid);
            if (!hash_get(&ht_edges, edge_vid)) {
                if (edge_test(scplx, nc, (int[2]){i,j}, alpha, has_big_tetra)) {
                    if (!mark_recursive_edge(edge_vid, &ht_edges, vmark)) goto fail;
                } else {
                    if (!hash_insert(&ht_edges, edge_vid, &f)) goto fail;
                }
            }
        }
    }

    /* Detect vertices */
    for (s_ncell *nc = scplx->head; nc; nc = nc->next) {
        for (int i=0; i<4; i++) {
            int vid = nc->vertex_id[i];
            if (!vmark[vid].seen) {
                if (vertex_test(scplx, nc, i, alpha, has_big_tetra, buff_ncellPTR)) {
                    vmark[vid].in_acplx = true;
                    vmark[vid].seen = true;
                } else {
                    vmark[vid].in_acplx = false;
                    vmark[vid].seen = true;
                }
            }
        }
    }

    /* output */
    *out_faces = ht_faces;
    *out_edges = ht_edges;
    for (int i=0; i<scplx->points.N; i++) out_vertices[i] = vmark[i].in_acplx;

    free(vmark);
    return 1;

fail:
    hash_free(&ht_faces);
    hash_free(&ht_edges);
    free(vmark);
    return 0;
}





/* ------------------------------------------------------------------------- */
/* Public API: alpha_complex_3d                                              */
/* ------------------------------------------------------------------------- */

void free_alpha_complex(s_alpha_complex *ac)
{
    if (!ac) return;
    free(ac->points.p);
    free(ac->weights);
    free(ac->tets);
    free(ac->faces);
    free(ac->edges);
    *ac = (s_alpha_complex){ .alpha = NAN };
}

static int cmp_lex(const int *a, const int *b, int k)
{
    for (int i = 0; i < k; i++) if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    return 0;
}
static int cmp_int2(const void *a, const void *b) { return cmp_lex(a, b, 2); }
static int cmp_int3(const void *a, const void *b) { return cmp_lex(a, b, 3); }
static int cmp_int4(const void *a, const void *b) { return cmp_lex(a, b, 4); }

/* Collect the keys (k sorted DT vids) of every entry of ht whose bool value is
 * true, remapped through k_of_vid (monotone, so each record stays ascending),
 * into a fresh array sorted lexicographically.  Returns the count, or -1 on
 * allocation failure.  *out is NULL when the count is 0. */
static int collect_true_keys(const s_hash_table *ht, int k, const int *k_of_vid, int **out)
{
    *out = NULL;
    int n = 0;
    for (size_t b = 0; b < ht->nbuckets; b++)
        for (s_hash_entry *e = ht->buckets[b]; e; e = e->next)
            if (*(bool *)((char *)e + ht->value_offset)) n++;
    if (n == 0) return 0;
    int *arr = malloc((size_t)n * k * sizeof(int));
    if (!arr) return -1;
    int m = 0;
    for (size_t b = 0; b < ht->nbuckets; b++)
        for (s_hash_entry *e = ht->buckets[b]; e; e = e->next) {
            if (!*(bool *)((char *)e + ht->value_offset)) continue;
            const int *key = entry_key(e);
            for (int j = 0; j < k; j++) arr[m*k + j] = k_of_vid[key[j]];
            m++;
        }
    qsort(arr, n, k * sizeof(int), k == 2 ? cmp_int2 : cmp_int3);
    *out = arr;
    return n;
}

s_alpha_complex alpha_complex_3d(const s_points *pts, const double *weights,
                                 double alpha, double TOL_dup,
                                 s_random_context *rng, int *out_kept_idx)
{
    const s_alpha_complex FAIL = { .alpha = NAN };
    int N_in = pts ? pts->N : 0;
    if (out_kept_idx) for (int i = 0; i < N_in; i++) out_kept_idx[i] = -1;
    if (!pts || pts->N < 4 || isnan(alpha)) return FAIL;

    /* 1. DT with the sentinels kept; kept[i] = DT vid of input i, or -1. */
    int *kept = malloc(pts->N * sizeof(int));
    if (!kept) return FAIL;
    for (int i = 0; i < pts->N; i++) kept[i] = i;
    s_dt_builder b = dt_builder_begin(pts, weights, TOL_dup, NULL, NULL, rng);
    if (!b._stack) { free(kept); return FAIL; }
    s_scplx dt = dt_builder_end(&b, /*keep_big_tetra=*/true, NULL, kept, pts->N);

    bool has_real_tet = false;
    for (s_ncell *nc = dt.head; nc && !has_real_tet; nc = nc->next)
        has_real_tet = nc->vertex_id[0] >= 4 && nc->vertex_id[1] >= 4 &&
                       nc->vertex_id[2] >= 4 && nc->vertex_id[3] >= 4;
    if (!has_real_tet) { free(kept); free_complex(&dt); return FAIL; }

    /* 2. Mark K. */
    s_alpha_complex ac = { .alpha = alpha };
    s_hash_table ht_faces, ht_edges;
    bool marked = false;
    int *k_of_vid = NULL;
    bool *vin = malloc(dt.points.N * sizeof(bool));
    s_dynarray buff = dynarray_initialize(sizeof(s_ncell *), 64);
    if (!vin || !buff.items) goto fail;
    if (!alpha_complex_mark(&dt, true, alpha, &buff, &ht_faces, &ht_edges, vin)) goto fail;
    marked = true;

    /* 3. Compact the vertices of K: DT vid -> K index (monotone in vid). */
    k_of_vid = malloc(dt.points.N * sizeof(int));
    if (!k_of_vid) goto fail;
    int nv = 0;
    for (int v = 0; v < dt.points.N; v++) k_of_vid[v] = (v >= 4 && vin[v]) ? nv++ : -1;
    if (nv > 0) {
        ac.points.p = malloc(nv * sizeof(s_point));
        if (!ac.points.p) goto fail;
        if (weights) {
            ac.weights = malloc(nv * sizeof(double));
            if (!ac.weights) goto fail;
        }
        for (int v = 4; v < dt.points.N; v++) {
            if (k_of_vid[v] < 0) continue;
            ac.points.p[k_of_vid[v]] = dt.points.p[v];
            if (weights) ac.weights[k_of_vid[v]] = dt.weights[v];
        }
    }
    ac.points.N = nv;

    /* 4. Simplices.  mark only sets mask_alpha on real tets (sentinel tets fail
     * acplx_tet_in), and every vertex of a K-simplex is in K. */
    int nt = 0;
    for (s_ncell *nc = dt.head; nc; nc = nc->next) nt += nc->mask_alpha;
    if (nt > 0) {
        ac.tets = malloc((size_t)nt * 4 * sizeof(int));
        if (!ac.tets) goto fail;
        int m = 0;
        for (s_ncell *nc = dt.head; nc; nc = nc->next) {
            if (!nc->mask_alpha) continue;
            int t[4];
            for (int j = 0; j < 4; j++) t[j] = k_of_vid[nc->vertex_id[j]];
            sort4(t);
            memcpy(&ac.tets[4*m++], t, sizeof t);
        }
        qsort(ac.tets, nt, 4 * sizeof(int), cmp_int4);
    }
    ac.N_tets = nt;
    ac.N_faces = collect_true_keys(&ht_faces, 3, k_of_vid, &ac.faces);
    ac.N_edges = collect_true_keys(&ht_edges, 2, k_of_vid, &ac.edges);
    if (ac.N_faces < 0 || ac.N_edges < 0) goto fail;

    /* 5. Input -> output map. */
    if (out_kept_idx)
        for (int i = 0; i < pts->N; i++)
            out_kept_idx[i] = kept[i] >= 0 ? k_of_vid[kept[i]] : -1;

    hash_free(&ht_faces);
    hash_free(&ht_edges);
    dynarray_free(&buff);
    free(vin);
    free(k_of_vid);
    free(kept);
    free_complex(&dt);
    return ac;

fail:   /* collect_true_keys leaves its *out NULL on failure, so ac is freeable */
    free_alpha_complex(&ac);
    if (marked) { hash_free(&ht_faces); hash_free(&ht_edges); }
    if (buff.items) dynarray_free(&buff);
    free(vin);
    free(k_of_vid);
    free(kept);
    free_complex(&dt);
    if (out_kept_idx) for (int i = 0; i < pts->N; i++) out_kept_idx[i] = -1;
    return FAIL;
}
