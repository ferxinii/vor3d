
#include "scplx.h"
#include "linalg.h"
#include "hash.h"
#include "dynarray.h"
#include <float.h>


/* Orthosphere of every simplex. */

static void tetra_orthosphere(const s_point p[4], const double w[4], double EPS_DEGEN, 
                              s_point *out_c, double *out_r2)
{
    double A[3][3], b[3];
    double r0_2 = norm_squared(p[0]);
    for (int i = 1; i <= 3; i++) {
        double ri_2 = norm_squared(p[i]);
        b[i-1] = ri_2 - r0_2 - (w[i] - w[0]);
        for (int j = 0; j < 3; j++)
            A[i-1][j] = 2.0 * (p[i].coords[j] - p[0].coords[j]);
    }

    s_point c;
    if (solve_3x3_ppivot_inplace(A, b, c.coords, EPS_DEGEN) != 3) goto error;

    if (out_c) *out_c = c;
    if (out_r2) *out_r2 = distance_squared(c, p[0]) - w[0];
    return;

    error:
        if (out_c) { (*out_c).x = DBL_MAX; (*out_c).y = DBL_MAX; (*out_c).z = DBL_MAX; }
        *out_r2 = DBL_MAX;
}

static void face_orthosphere(const s_point p[3], const double w[3], double EPS_DEGEN,
                             s_point *out_c, double *out_r2)
{   /* Unique sphere centered on the face with same power to the 3 vertices */
    s_point n, t1, t2;
    if (basis_vectors_plane(p, EPS_DEGEN, &n, &t1, &t2) == 0) goto error;

    /* Project all 3 vertices onto the plane's 2D coordinate system */
    double p2D[2][2];  /* p0 is mapped to (0,0) */
    project_point_to_plane_2D(p[1], p[0], t1, t2, p2D[0]);
    project_point_to_plane_2D(p[2], p[0], t1, t2, p2D[1]);

    /* Build 2x2 system, accounting for p0 = (0,0) */
    double A[2][2], b[2];
    for (int i=0; i<2; i++) {
        A[i][0] = 2.0 * p2D[i][0];
        A[i][1] = 2.0 * p2D[i][1];
        b[i] = p2D[i][0]*p2D[i][0] + p2D[i][1]*p2D[i][1] - (w[i+1] - w[0]);
    }

    double c[2];  /* Circle center, we won't return it */
    if (solve_2x2_cramer(A, b, c, EPS_DEGEN) != 2) goto error;

    if (out_c) *out_c = sum_points(p[0], sum_points(scale_point(t1, c[0]),
                                         scale_point(t2, c[1])));
    if (out_r2) *out_r2 = c[0]*c[0] + c[1]*c[1] - w[0];
    return;

    error:
        if (out_c) { (*out_c).x = DBL_MAX; (*out_c).y = DBL_MAX; (*out_c).z = DBL_MAX; }
        *out_r2 = DBL_MAX;
}

static void edge_orthosphere(const s_point p[2], const double w[2], double EPS_DEGEN,
                             s_point *out_c, double *out_r2)
{   /* Unique orthosphere centered on the edge with same power to both endpoints */
    s_point p01 = subtract_points(p[1], p[0]);
    double d2 = norm_squared(p01);
    if (d2 < EPS_DEGEN) goto error;
    double t = 0.5 + (w[0]-w[1])/(2*d2);
    
    if (out_c) *out_c = sum_points(p[0], scale_point(p01, t));
    if (out_r2) *out_r2 = t*t*d2 - w[0];
    return;

    error:
        if (out_c) { (*out_c).x = DBL_MAX; (*out_c).y = DBL_MAX; (*out_c).z = DBL_MAX; }
        if (out_r2) *out_r2 = DBL_MAX;
}



/* Test for every simplex */

static bool tetra_test(const s_scplx *scplx, const s_ncell *nc, double alpha, bool has_big_tetra,
                       double EPS_DEGEN)
{
    if (has_big_tetra && (nc->vertex_id[0] < 4 || nc->vertex_id[1] < 4 ||
                          nc->vertex_id[2] < 4 || nc->vertex_id[3] < 4 )) {
        return false;
    }

    s_point p[4];  extract_vertices_ncell(scplx, nc, p);
    double  w[4];  extract_weights_ncell(scplx, nc, w);
    double r2; tetra_orthosphere(p, w, EPS_DEGEN, NULL, &r2);
    return (r2 <= alpha);  /* Only radius contiditon necessary */
}


static bool face_test(const s_scplx *scplx, const s_ncell *nc, int face_localid,
                      double alpha, bool has_big_tetra, double EPS_DEGEN)
{
    if (has_big_tetra) {
        int fid[3]; extract_ids_face(nc, 2, &face_localid, fid);
        if (fid[0] < 4 || fid[1] < 4 || fid[2] < 4) return false;
    }

    /* 1) Radius condition. */
    s_point p[3]; double w[3];
    extract_vertices_and_weights_face(scplx, nc, 2, &face_localid, p, w);

    s_point c; double r2;
    face_orthosphere(p, w, EPS_DEGEN, &c, &r2);
    if (r2 > alpha) return false;  

    /* 2) Power condition. */
    /* Vertex opposite in nc */
    int vid1 = nc->vertex_id[face_localid];
    if (!has_big_tetra || vid1 >= 4)
        if (power_distance_point_vertex(scplx, vid1, c) < r2) return false;
    /* Vertex opposite in opp tet */
    s_ncell *opp = nc->opposite[face_localid];
    if (opp) {
        int opp_local;
        face_localid_of_adjacent_ncell(nc, 2, &face_localid, face_localid, &opp_local);
        int vid2 = opp->vertex_id[opp_local];
        if (!has_big_tetra || vid2 >= 4)
            if (power_distance_point_vertex(scplx, vid2, c) < r2) return false;
    }

    return true;
}


typedef struct {
    s_scplx *scplx;
    int gid0, gid1;
    bool has_big_tetra;
    s_point c; double r2;  /* edge's orthosphere */
    bool empty_ball_ok;   /* output */
} edge_test_ctx;

static bool AUX_check_edge_acplx(void *C, const s_ncell *nc)
{   /* Checks power condition */
    edge_test_ctx *ctx = C;
    if (ctx->empty_ball_ok) {
        for (int v = 0; v < 4; ++v) {
            int vid = nc->vertex_id[v];
            if (ctx->has_big_tetra && vid < 4) continue;
            if (vid == ctx->gid0 || vid == ctx->gid1) continue;
            if (power_distance_point_vertex(ctx->scplx, vid, ctx->c) < ctx->r2) {
                ctx->empty_ball_ok = false;
                return false;
            }
        }
    }
    return true;
}

static bool edge_test(s_scplx *scplx, const s_ncell *nc, const int edge_endpoints[2],
                      double alpha, bool has_big_tetra, double EPS_DEGEN)
{
    int gid0 = nc->vertex_id[edge_endpoints[0]];
    int gid1 = nc->vertex_id[edge_endpoints[1]];
    if (has_big_tetra && (gid0 < 4 || gid1 < 4)) return false;

    int omitted[2];  /* The other two local ids */
    for (int k=0, i=0; i<4; i++) if (i != edge_endpoints[0] && i != edge_endpoints[1])
        omitted[k++] = i;

    s_point p[2]; double w[2];
    extract_vertices_and_weights_face(scplx, nc, 1, omitted, p, w);
    s_point c; double r2;
    edge_orthosphere(p, w, EPS_DEGEN, &c, &r2);

    /* 1) Radius condition. */
    if (r2 > alpha) return false;

    /* 2) Power condition by walking ridge. */
    edge_test_ctx ctx = {
        .scplx = scplx,
        .gid0 = gid0, .gid1 = gid1,
        .has_big_tetra = has_big_tetra,
        .c = c, .r2 = r2,
        .empty_ball_ok = true,
    };
    walk_ridge_cycle_and_check_ncells(nc, omitted, AUX_check_edge_acplx, &ctx);
    return ctx.empty_ball_ok;
}


static bool vertex_test(s_scplx *scplx, s_ncell *nc, int v_localid, double alpha,
                        bool has_big_tetra, s_dynarray *buff_ncellPTR)
{
    if (has_big_tetra && (nc->vertex_id[v_localid] < 4)) return false;

    /* 1) Radius condition */
    int vid = nc->vertex_id[v_localid];
    double r2 = scplx->weights ? -scplx->weights[vid] : 0.0;  /* Orthosphere's r2 is -w */
    if (alpha < r2) return false;

    /* 2) Power condition by checking star of vertex */
    int v_localid_ARR[3]; extract_ids_face(nc, 2, &v_localid, v_localid_ARR);
    ncells_incident_face(scplx, nc, 0, v_localid_ARR, buff_ncellPTR);

    s_point pi = scplx->points.p[vid];
    for (unsigned i=0; i<buff_ncellPTR->N; i++) {
        s_ncell **tc = dynarray_get_ptr(buff_ncellPTR, i);
        for (int j=0; j<4; j++) {
            int vjd = (*tc)->vertex_id[j];
            if (vjd == vid) continue;
            if (has_big_tetra && vjd < 4) continue;
            if (power_distance_point_vertex(scplx, vjd, pi) < r2) return false;
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


void extract_alpha_complex(s_scplx *scplx, bool has_big_tetra, double alpha, 
                           double EPS_DEGEN, s_dynarray *buff_ncellPTR,
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
        if (tetra_test(scplx, nc, alpha, has_big_tetra, EPS_DEGEN))
            mark_recursive_ncell(nc, &ht_faces, &ht_edges, vmark);
    }

    /* Detect faces */
    for (s_ncell *nc = scplx->head; nc; nc = nc->next) for (int i=0; i<4; i++) {
        int face_vid[3]; extract_ids_face(nc, 2, &i, face_vid); sort3(face_vid);
        if (!hash_get(&ht_faces, face_vid)) { 
            if (face_test(scplx, nc, i, alpha, has_big_tetra, EPS_DEGEN)) {
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
                if (edge_test(scplx, nc, (int[2]){i,j}, alpha, has_big_tetra, EPS_DEGEN)) {
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


