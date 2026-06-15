
#include "scplx.h"
#include "hash.h"
#include "dynarray.h"
#include "gtests.h"
#include <float.h>


/* Test for every simplex: return true if simplex satisfies radius and conflict conditions */
static bool tetra_test(const s_scplx *scplx, const s_ncell *nc, double alpha, bool has_big_tetra)
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

    (void)buff_ncellPTR;
    // /* 2) Conflict condition by checking star of vertex */
    // int v_localid_ARR[3]; extract_ids_face(nc, 2, &v_localid, v_localid_ARR);
    // ncells_incident_face(scplx, nc, 0, v_localid_ARR, buff_ncellPTR);
    //
    // for (unsigned i=0; i<buff_ncellPTR->N; i++) {
    //     s_ncell **tc = dynarray_get_ptr(buff_ncellPTR, i);
    //     for (int j=0; j<4; j++) {
    //         int vjd = (*tc)->vertex_id[j];
    //         if (vjd == vid) continue;
    //         if (has_big_tetra && vjd < 4) continue;
    //         s_point pj = scplx->points.p[vjd];
    //         double  wj  = scplx->weights ? scplx->weights[vjd] : 0.0;
    //         if (test_orthosphere(1, &pi_pt, &pi_wt, pj, wj) > 0) return false;
    //     }
    // }

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


void extract_alpha_complex(s_scplx *scplx, bool has_big_tetra, double alpha, s_dynarray *buff_ncellPTR,
                           s_hash_table *out_faces, s_hash_table *out_edges, 
                           bool *out_vertices)
{
    int N_tetra = scplx->N_ncells;
    s_hash_table ht_faces; hash_init(&ht_faces, sizeof(int)*3, sizeof(bool), N_tetra*4, 
                                     N_tetra*4, hash_fun_faces, hash_eq_faces, NULL);
    s_hash_table ht_edges; hash_init(&ht_edges, sizeof(int)*2, sizeof(bool), N_tetra*6, 
                                     N_tetra*6, hash_fun_edges, hash_eq_edges, NULL);
    s_vertex_mark *vmark = calloc(scplx->points.N, sizeof(s_vertex_mark));
    bool f = false;

    /* Detect ncells */
    for (s_ncell *nc = scplx->head; nc; nc = nc->next) {
        if (tetra_test(scplx, nc, alpha, has_big_tetra))
            mark_recursive_ncell(nc, &ht_faces, &ht_edges, vmark);
    }

    /* Detect faces */
    for (s_ncell *nc = scplx->head; nc; nc = nc->next) for (int i=0; i<4; i++) {
        int face_vid[3]; extract_ids_face(nc, 2, &i, face_vid); sort3(face_vid);
        if (!hash_get(&ht_faces, face_vid)) { 
            if (face_test(scplx, nc, i, alpha, has_big_tetra)) {
                mark_recursive_face(face_vid, &ht_faces, &ht_edges, vmark);
            } else {
                hash_insert(&ht_faces, face_vid, &f);
            }
        }
    }

    /* Detect edges */
    for (s_ncell *nc = scplx->head; nc; nc = nc->next) {
        for (int i=0; i<3; i++) for (int j=i+1; j<4; j++) {
            int edge_vid[2] = {nc->vertex_id[i], nc->vertex_id[j]}; sort2(edge_vid);
            if (!hash_get(&ht_edges, edge_vid)) {
                if (edge_test(scplx, nc, (int[2]){i,j}, alpha, has_big_tetra)) {
                    mark_recursive_edge(edge_vid, &ht_edges, vmark);
                } else {
                    bool f = false;
                    hash_insert(&ht_edges, edge_vid, &f);
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
}


