#include "delaunay.h"
#include "scplx.h"
#include "trimesh.h"
#include "dynarray.h"
#include "hash.h"
#include "gtests.h"
#include "points.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>


static bool edge_in_dt(s_scplx *dt, s_ncell **vertex_start, int *v_local_id,
                        int va, int vb)
{
    if (!vertex_start[va]) return false;

    s_dynarray out_ids = dynarray_initialize(sizeof(int), 16);
    s_dynarray scratch  = dynarray_initialize(sizeof(s_ncell *), 16);
    if (!out_ids.items || !scratch.items) {
        dynarray_free(&out_ids); dynarray_free(&scratch);
        return false;
    }

    vertex_neighbors(dt, va, vertex_start[va], v_local_id[va], 0, &out_ids, &scratch);

    bool found = false;
    for (unsigned i = 0; i < out_ids.N; i++) {
        int nb; dynarray_get_value(&out_ids, i, &nb);
        if (nb == vb) { found = true; break; }
    }
    dynarray_free(&out_ids);
    dynarray_free(&scratch);
    return found;
}

/* Replace every face containing edge (va,vb) with two faces split at new_id. */
static int split_trimesh_edge(s_trimesh *mesh, int va, int vb, int new_id)
{
    int n_split = 0;
    for (int fi = 0; fi < mesh->Nf; fi++) {
        bool has_va = false, has_vb = false;
        for (int j = 0; j < 3; j++) {
            if (mesh->faces[fi*3+j] == va) has_va = true;
            if (mesh->faces[fi*3+j] == vb) has_vb = true;
        }
        if (has_va && has_vb) n_split++;
    }

    int new_Nf = mesh->Nf + n_split;
    int     *new_faces    = malloc((size_t)new_Nf * 3 * sizeof(int));
    s_point *new_fnormals = malloc((size_t)new_Nf     * sizeof(s_point));
    if (!new_faces || !new_fnormals) {
        free(new_faces); free(new_fnormals); return 0;
    }

    int out = 0;
    for (int fi = 0; fi < mesh->Nf; fi++) {
        int v0 = mesh->faces[fi*3], v1 = mesh->faces[fi*3+1], v2 = mesh->faces[fi*3+2];
        bool has_va = (v0==va || v1==va || v2==va);
        bool has_vb = (v0==vb || v1==vb || v2==vb);

        if (!has_va || !has_vb) {
            new_faces[out*3]   = v0;
            new_faces[out*3+1] = v1;
            new_faces[out*3+2] = v2;
            new_fnormals[out]  = mesh->fnormals[fi];
            out++;
        } else {
            int vc = -1;
            int vs[3] = {v0, v1, v2};
            for (int j = 0; j < 3; j++)
                if (vs[j] != va && vs[j] != vb) { vc = vs[j]; break; }
            assert(vc != -1);

            new_faces[out*3]   = va;     new_faces[out*3+1] = new_id; new_faces[out*3+2] = vc;
            new_fnormals[out]  = mesh->fnormals[fi];
            out++;
            new_faces[out*3]   = new_id; new_faces[out*3+1] = vb;     new_faces[out*3+2] = vc;
            new_fnormals[out]  = mesh->fnormals[fi];
            out++;
        }
    }

    free(mesh->faces);    mesh->faces    = new_faces;
    free(mesh->fnormals); mesh->fnormals = new_fnormals;
    free(mesh->adjacency); mesh->adjacency = NULL;
    mesh->Nf = new_Nf;
    return 1;
}

static size_t edge_pair_hash(const void *key)
{
    const int *k = (const int *)key;
    size_t h = 2166136261u;
    h ^= (size_t)(unsigned int)k[0]; h *= 16777619u;
    h ^= (size_t)(unsigned int)k[1]; h *= 16777619u;
    return h;
}

static bool edge_pair_eq(const void *a, const void *b)
{
    return memcmp(a, b, 2 * sizeof(int)) == 0;
}

/* Ensure every trimesh boundary edge is present in the DT.
 * Processes ALL unprotected edges per pass (not one at a time) so each pass
 * is O(N_faces) instead of O(N_faces^2) and the total work is proportional
 * to the number of edges that actually need splitting.
 * mesh->faces holds 0-based canonical indices; the +4 builder offset is
 * applied only at the DT lookup call sites here -- the trimesh is never shifted. */
static int ensure_ridge_protected(s_dt_builder *b, s_trimesh *mesh, double TOL)
{
    int _pass = 0;
    for (;;) {
        _pass++;

        /* Build vertex-to-cell index once per pass. */
        s_ncell **vertex_start = calloc((size_t)b->dt.points.N, sizeof(s_ncell *));
        int      *v_local_id   = calloc((size_t)b->dt.points.N, sizeof(int));
        if (!vertex_start || !v_local_id) {
            free(vertex_start); free(v_local_id); return 0;
        }
        build_vertex_cell_index(&b->dt, vertex_start, v_local_id);

        /* Collect all unique unprotected edges (canonical min<max pairs). */
        s_hash_table edge_set;
        int dummy = 1;
        if (!hash_init(&edge_set, sizeof(int[2]), sizeof(int),
                       (size_t)(mesh->Nf * 4 + 1), (size_t)(mesh->Nf),
                       edge_pair_hash, edge_pair_eq, NULL)) {
            free(vertex_start); free(v_local_id); return 0;
        }

        s_dynarray to_split = dynarray_initialize(sizeof(int[2]), 64);
        if (!to_split.items) {
            hash_free(&edge_set); free(vertex_start); free(v_local_id); return 0;
        }

        for (int fi = 0; fi < mesh->Nf; fi++) {
            for (int ei = 0; ei < 3; ei++) {
                int va = mesh->faces[fi*3 + (ei+1)%3];
                int vb = mesh->faces[fi*3 + (ei+2)%3];
                if (edge_in_dt(&b->dt, vertex_start, v_local_id, va + 4, vb + 4))
                    continue;

                int pair[2] = { va < vb ? va : vb, va < vb ? vb : va };
                if (hash_get(&edge_set, pair)) continue;
                if (!hash_insert(&edge_set, pair, &dummy)) {
                    hash_free(&edge_set); dynarray_free(&to_split);
                    free(vertex_start); free(v_local_id); return 0;
                }
                if (!dynarray_push(&to_split, pair)) {
                    hash_free(&edge_set); dynarray_free(&to_split);
                    free(vertex_start); free(v_local_id); return 0;
                }
            }
        }

        free(vertex_start);
        free(v_local_id);
        hash_free(&edge_set);

        if (to_split.N == 0) {
            dynarray_free(&to_split);
            break;
        }

        fprintf(stderr, "[CDT]   ridge pass %d: Nf=%d Nv=%d  %u edges to split\n",
                _pass, mesh->Nf, b->dt.points.N, to_split.N); fflush(stderr);

        /* Insert midpoints into DT and split trimesh edges in batch. */
        for (unsigned i = 0; i < to_split.N; i++) {
            int *pair = (int *)dynarray_get_ptr(&to_split, i);
            int va = pair[0], vb = pair[1];

            s_point pa = b->dt.points.p[va + 4];
            s_point pb = b->dt.points.p[vb + 4];
            s_point M  = { .x = (pa.x + pb.x) * 0.5,
                           .y = (pa.y + pb.y) * 0.5,
                           .z = (pa.z + pb.z) * 0.5 };

            int new_id_builder = b->dt.points.N;
            s_points single = { .N = 1, .p = &M };
            if (!dt_builder_extend(b, &single, TOL)) {
                dynarray_free(&to_split); return 0;
            }

            if (!split_trimesh_edge(mesh, va, vb, new_id_builder - 4)) {
                dynarray_free(&to_split); return 0;
            }
        }

        dynarray_free(&to_split);
    }
    return 1;
}


/* ----------------------------------------------------------------------- */


static s_ncell *face_in_dt(const s_scplx *dt, int va, int vb, int vc)
{
    for (s_ncell *nc = dt->head; nc; nc = nc->next) {
        bool ha = false, hb = false, hc = false;
        for (int j = 0; j < 4; j++) {
            if (nc->vertex_id[j] == va) ha = true;
            if (nc->vertex_id[j] == vb) hb = true;
            if (nc->vertex_id[j] == vc) hc = true;
        }
        if (ha && hb && hc) return nc;
    }
    return NULL;
}

static size_t face_triple_hash(const void *key)
{
    const int *k = (const int *)key;
    size_t h = 2166136261u;
    for (int i = 0; i < 3; i++) {
        h ^= (size_t)(unsigned int)k[i];
        h *= 16777619u;
    }
    return h;
}

static bool face_triple_eq(const void *a, const void *b)
{
    return memcmp(a, b, 3 * sizeof(int)) == 0;
}

static void sort3(int *a, int *b, int *c)
{
    if (*a > *b) { int t = *a; *a = *b; *b = t; }
    if (*b > *c) { int t = *b; *b = *c; *c = t; }
    if (*a > *b) { int t = *a; *a = *b; *b = t; }
}

static bool is_constrained_face(s_hash_table *ht, int va, int vb, int vc)
{
    int key[3] = {va, vb, vc};
    sort3(&key[0], &key[1], &key[2]);
    return hash_get(ht, key) != NULL;
}

static void mark_constrained_face(s_hash_table *ht, int va, int vb, int vc)
{
    int key[3] = {va, vb, vc};
    sort3(&key[0], &key[1], &key[2]);
    int val = 1;
    hash_insert(ht, key, &val);
}

/* --- Kinetic CDT helpers (Steps 1 & 2) --------------------------------- */

static double det3(double a0, double a1, double a2,
                   double b0, double b1, double b2,
                   double c0, double c1, double c2)
{
    return a0*(b1*c2 - b2*c1) - a1*(b0*c2 - b2*c0) + a2*(b0*c1 - b1*c0);
}

static double det4(const double m[4][4])
{
    return  m[0][0]*det3(m[1][1],m[1][2],m[1][3], m[2][1],m[2][2],m[2][3], m[3][1],m[3][2],m[3][3])
           -m[0][1]*det3(m[1][0],m[1][2],m[1][3], m[2][0],m[2][2],m[2][3], m[3][0],m[3][2],m[3][3])
           +m[0][2]*det3(m[1][0],m[1][1],m[1][3], m[2][0],m[2][1],m[2][3], m[3][0],m[3][1],m[3][3])
           -m[0][3]*det3(m[1][0],m[1][1],m[1][2], m[2][0],m[2][1],m[2][2], m[3][0],m[3][1],m[3][2]);
}

/* Signed 4x4 in-sphere determinant.
 * Rows are (v-q, |v|^2-|q|^2) for v in {a,b,c,p}; q is the test point.
 * Positive result means q is inside the circumsphere of a,b,c,p (not locally Delaunay). */
static double insphere_det_float(s_point a, s_point b, s_point c,
                                  s_point p, s_point q)
{
    double qq = q.x*q.x + q.y*q.y + q.z*q.z;
    s_point verts[4] = {a, b, c, p};
    double m[4][4];
    for (int i = 0; i < 4; i++) {
        m[i][0] = verts[i].x - q.x;
        m[i][1] = verts[i].y - q.y;
        m[i][2] = verts[i].z - q.z;
        m[i][3] = (verts[i].x*verts[i].x + verts[i].y*verts[i].y + verts[i].z*verts[i].z) - qq;
    }
    return det4(m);
}

/* Returns tau* >= 0 at which face g loses local regularity under kinetic weight w(v)=tau*d(v).
 * s is the tet NOT in R (all d(v) <= 0); s->opposite[opp_lid] is the tet IN R.
 * Returns -1.0 if the face never loses regularity. */
static double certify(const s_scplx *dt,
                      const s_ncell *s, int opp_lid,
                      s_point h_normal, s_point h_point)
{
    s_ncell *t = s->opposite[opp_lid];
    if (!t) return -1.0;

    int L_id = s->vertex_id[opp_lid];
    s_point L = dt->points.p[L_id];

    int face_ids[3]; int idx = 0;
    s_point face_pts[3];
    for (int j = 0; j < 4; j++) {
        if (j == opp_lid) continue;
        face_ids[idx] = s->vertex_id[j];
        face_pts[idx] = dt->points.p[s->vertex_id[j]];
        idx++;
    }
    s_point a = face_pts[0], b = face_pts[1], c = face_pts[2];

    int Rv_id = -1;
    for (int j = 0; j < 4; j++) {
        int vid = t->vertex_id[j];
        if (vid != face_ids[0] && vid != face_ids[1] && vid != face_ids[2]) {
            Rv_id = vid; break;
        }
    }
    if (Rv_id < 0) return -1.0;
    s_point Rv = dt->points.p[Rv_id];

    double d_a  = dot_prod(subtract_points(a,  h_point), h_normal);
    double d_b  = dot_prod(subtract_points(b,  h_point), h_normal);
    double d_c  = dot_prod(subtract_points(c,  h_point), h_normal);
    double d_L  = dot_prod(subtract_points(L,  h_point), h_normal);
    double d_Rv = dot_prod(subtract_points(Rv, h_point), h_normal);

    int sign_D0 = test_insphere((s_point[4]){a, b, c, L}, Rv);
    if (sign_D0 == 0) return 0.0;  /* Rv exactly on circumsphere: flip immediately */

    double D0 = insphere_det_float(a, b, c, L, Rv);
    /* Float insphere can have wrong sign when points are far from origin (cancellation
     * in |v|^2 - |q|^2). Correct with the robust sign; tau ratio stays valid. */
    if (D0 != 0.0 && (sign_D0 > 0) != (D0 > 0.0)) D0 = -D0;

    /* D': same first-3-col geometry but 4th col = d(Rv) - d(v) for each row */
    double m[4][4];
    s_point pts[4] = {a, b, c, L};
    double dv[4]   = {d_a, d_b, d_c, d_L};
    for (int i = 0; i < 4; i++) {
        m[i][0] = pts[i].x - Rv.x;
        m[i][1] = pts[i].y - Rv.y;
        m[i][2] = pts[i].z - Rv.z;
        m[i][3] = d_Rv - dv[i];
    }
    double Dprime = det4(m);

    if (Dprime == 0.0) return -1.0;
    double tau = -D0 / Dprime;
    return (tau >= 0.0) ? tau : -1.0;
}

/* --- Min-heap on tau (Step 3) ------------------------------------------- */

typedef struct {
    double   tau;
    int      va, vb, vc;  /* sorted face triple */
    s_ncell *s;           /* tet NOT in R at push time; re-verified on pop */
    int      opp_lid;
} heap_entry;

typedef struct { s_dynarray data; } s_heap;

static s_heap heap_init(void)
{
    return (s_heap){ .data = dynarray_initialize(sizeof(heap_entry), 32) };
}

static void heap_free(s_heap *h) { dynarray_free(&h->data); }
static int  heap_empty(const s_heap *h) { return h->data.N == 0; }

static int heap_push(s_heap *h, double tau, int va, int vb, int vc,
                     s_ncell *s, int opp_lid)
{
    heap_entry e = { tau, va, vb, vc, s, opp_lid };
    if (!dynarray_push(&h->data, &e)) return 0;

    unsigned i = h->data.N - 1;
    while (i > 0) {
        unsigned parent = (i - 1) / 2;
        heap_entry *ei = dynarray_get_ptr(&h->data, i);
        heap_entry *ep = dynarray_get_ptr(&h->data, parent);
        if (ep->tau <= ei->tau) break;
        heap_entry tmp = *ei; *ei = *ep; *ep = tmp;
        i = parent;
    }
    return 1;
}

static heap_entry heap_pop_min(s_heap *h)
{
    heap_entry result;
    dynarray_get_value(&h->data, 0, &result);

    if (h->data.N > 1) {
        /* Read old last element before decrementing — dynarray_get_ptr returns
         * NULL for id >= N, so we must fetch it while it's still in range. */
        heap_entry *last = dynarray_get_ptr(&h->data, h->data.N - 1);
        h->data.N--;
        *(heap_entry *)dynarray_get_ptr(&h->data, 0) = *last;

        unsigned i = 0;
        for (;;) {
            unsigned l = 2*i + 1, r = 2*i + 2, sm = i;
            heap_entry *es = dynarray_get_ptr(&h->data, sm);
            if (l < h->data.N) {
                heap_entry *el = dynarray_get_ptr(&h->data, l);
                if (el->tau < es->tau) { sm = l; es = el; }
            }
            if (r < h->data.N) {
                heap_entry *er = dynarray_get_ptr(&h->data, r);
                if (er->tau < es->tau) sm = r;
            }
            if (sm == i) break;
            heap_entry *ei = dynarray_get_ptr(&h->data, i);
            es = dynarray_get_ptr(&h->data, sm);
            heap_entry tmp = *ei; *ei = *es; *es = tmp;
            i = sm;
        }
    } else {
        h->data.N--;  /* N was 1: just remove the single element */
    }
    return result;
}

/* --- Helpers for the queue-driven flip loop (Step 5) ------------------- */

static bool tet_in_R(const s_scplx *dt, const s_ncell *nc,
                     s_point pa, s_point pb, s_point pc)
{
    bool above = false, below = false;
    for (int j = 0; j < 4; j++) {
        s_point q = dt->points.p[nc->vertex_id[j]];
        int side = test_orientation((s_point[3]){pa, pb, pc}, q);
        if (side > 0) above = true;
        if (side < 0) below = true;
    }
    return above && below;
}

/* --- find_region_R_boundary (Step 4) ------------------------------------ */

/* BFS from a seed R-tet (straddling the constraint plane). Collects:
 *   boundary: sorted int[3] face triples between an R-tet and a non-R-tet
 *   in_R:     s_ncell* pointers for every tet in R
 * Constrained faces are excluded from boundary and not crossed during BFS.
 * Returns 1 on success, 0 on allocation failure. */
static int find_region_R_boundary(s_scplx *dt,
                                   s_point pa, s_point pb, s_point pc,
                                   s_hash_table *constrained,
                                   s_dynarray *boundary,
                                   s_dynarray *in_R)
{
    /* Find any starting tet that straddles the constraint plane */
    s_ncell *seed = NULL;
    for (s_ncell *nc = dt->head; nc; nc = nc->next) {
        if (tet_in_R(dt, nc, pa, pb, pc)) { seed = nc; break; }
    }
    if (!seed) return 1;  /* face already in DT */

    dt->mark_stamp++;
    int stamp = dt->mark_stamp;

    s_dynarray stack = dynarray_initialize(sizeof(s_ncell *), 64);
    if (!stack.items) return 0;

    seed->mark_token = stamp;
    dynarray_push(&stack, &seed);
    dynarray_push(in_R, &seed);

    int ret = 0;
    while (stack.N > 0) {
        s_ncell *cur;
        dynarray_pop(&stack, &cur);

        for (int fi = 0; fi < 4; fi++) {
            s_ncell *nb = cur->opposite[fi];
            if (!nb) continue;

            int fids[3];
            extract_ids_face(cur, 2, &fi, fids);
            sort3(&fids[0], &fids[1], &fids[2]);

            if (is_constrained_face(constrained, fids[0], fids[1], fids[2])) continue;

            if (tet_in_R(dt, nb, pa, pb, pc)) {
                if (nb->mark_token != stamp) {
                    nb->mark_token = stamp;
                    if (!dynarray_push(&stack, &nb)) goto done;
                    if (!dynarray_push(in_R, &nb))   goto done;
                }
            } else {
                if (!dynarray_push(boundary, fids)) goto done;
            }
        }
    }
    ret = 1;
done:
    dynarray_free(&stack);
    return ret;
}

/* Find the s-tet (NOT in R) for a boundary face (fa,fb,fc), and its opp_lid
 * such that s->opposite[opp_lid] is the R-tet.
 * Returns 0 if the face is stale (no longer an R-boundary face). */
static int find_boundary_face_s(const s_scplx *dt,
                                 s_point pa, s_point pb, s_point pc,
                                 int fa, int fb, int fc,
                                 s_ncell **out_s, int *out_opp_lid)
{
    for (s_ncell *nc = dt->head; nc; nc = nc->next) {
        bool ha = false, hb = false, hc = false;
        int opp = -1;
        for (int j = 0; j < 4; j++) {
            int v = nc->vertex_id[j];
            if      (v == fa) ha = true;
            else if (v == fb) hb = true;
            else if (v == fc) hc = true;
            else              opp = j;
        }
        if (!ha || !hb || !hc || opp < 0) continue;

        s_ncell *nb = nc->opposite[opp];
        if (!nb) return 0;

        bool nc_in_R = tet_in_R(dt, nc, pa, pb, pc);
        bool nb_in_R = tet_in_R(dt, nb, pa, pb, pc);
        if (nc_in_R == nb_in_R) return 0;  /* no longer an R-boundary face */

        if (!nc_in_R) {
            *out_s = nc; *out_opp_lid = opp;
        } else {
            int nb_opp = -1;
            for (int j = 0; j < 4; j++) {
                int v = nb->vertex_id[j];
                if (v != fa && v != fb && v != fc) { nb_opp = j; break; }
            }
            if (nb_opp < 0) return 0;
            *out_s = nb; *out_opp_lid = nb_opp;
        }
        return 1;
    }
    return 0;
}

/* ----------------------------------------------------------------------- */

static int flipinsertfacet(s_scplx *dt, int va, int vb, int vc,
                            s_hash_table *constrained,
                            s_point pa, s_point pb, s_point pc)
{
    if (face_in_dt(dt, va, vb, vc)) {
        mark_constrained_face(constrained, va, vb, vc);
        return 1;
    }

    s_point h_normal = normalize_vec(
        cross_prod(subtract_points(pb, pa), subtract_points(pc, pa)), 1e-12);
    s_point h_point  = pa;

    int ret = 0;
    s_dynarray boundary = {0};
    s_dynarray in_R     = {0};
    s_heap Q            = {0};

    boundary = dynarray_initialize(sizeof(int[3]), 32);
    in_R     = dynarray_initialize(sizeof(s_ncell *), 64);
    if (!boundary.items || !in_R.items) goto cleanup;

    if (!find_region_R_boundary(dt, pa, pb, pc, constrained, &boundary, &in_R))
        goto cleanup;

    unsigned diag_n_R = in_R.N, diag_n_bnd = boundary.N;

    Q = heap_init();
    if (!Q.data.items) goto cleanup;

    int diag_pushed = 0;
    for (unsigned i = 0; i < boundary.N; i++) {
        int *fids = dynarray_get_ptr(&boundary, i);
        s_ncell *s; int opp_lid;
        if (!find_boundary_face_s(dt, pa, pb, pc, fids[0], fids[1], fids[2], &s, &opp_lid))
            continue;
        double tau = certify(dt, s, opp_lid, h_normal, h_point);
        if (tau >= 0.0) { heap_push(&Q, tau, fids[0], fids[1], fids[2], s, opp_lid); diag_pushed++; }
    }
    dynarray_free(&boundary);
    dynarray_free(&in_R);

    int diag_flipped = 0, diag_notflipped = 0, diag_stale = 0, diag_constr = 0;
    while (!heap_empty(&Q)) {
        heap_entry e = heap_pop_min(&Q);

        s_ncell *s; int opp_lid;
        if (!find_boundary_face_s(dt, pa, pb, pc, e.va, e.vb, e.vc, &s, &opp_lid))
            { diag_stale++; continue; }
        if (is_constrained_face(constrained, e.va, e.vb, e.vc))
            { diag_constr++; continue; }

        s_ncell *new_tets[4]; int n_new = 0;
        if (dt_flip_face(dt, s, opp_lid, new_tets, &n_new) != 1) { diag_notflipped++; continue; }
        diag_flipped++;

        if (face_in_dt(dt, va, vb, vc)) break;

        /* Enqueue new R-boundary faces from the tets produced by the flip. */
        for (int ni = 0; ni < n_new; ni++) {
            s_ncell *nc = new_tets[ni];
            if (!tet_in_R(dt, nc, pa, pb, pc)) continue;
            for (int fi = 0; fi < 4; fi++) {
                s_ncell *nb = nc->opposite[fi];
                if (!nb || tet_in_R(dt, nb, pa, pb, pc)) continue;

                int fids[3];
                extract_ids_face(nc, 2, &fi, fids);
                sort3(&fids[0], &fids[1], &fids[2]);
                if (is_constrained_face(constrained, fids[0], fids[1], fids[2])) continue;

                int nb_opp = -1;
                for (int j = 0; j < 4; j++) {
                    int v = nb->vertex_id[j];
                    if (v != fids[0] && v != fids[1] && v != fids[2]) { nb_opp = j; break; }
                }
                if (nb_opp < 0) continue;

                double tau = certify(dt, nb, nb_opp, h_normal, h_point);
                if (tau >= 0.0)
                    heap_push(&Q, tau, fids[0], fids[1], fids[2], nb, nb_opp);
            }
        }
    }

    if (face_in_dt(dt, va, vb, vc)) {
        mark_constrained_face(constrained, va, vb, vc);
        ret = 1;
    } else {
        s_point cross = cross_prod(subtract_points(pb, pa), subtract_points(pc, pa));
        double cross_mag = sqrt(cross.x*cross.x + cross.y*cross.y + cross.z*cross.z);
        fprintf(stderr,
            "flipinsertfacet: heap exhausted for face (%d,%d,%d)\n"
            "  pa=(%.3f,%.3f,%.3f) pb=(%.3f,%.3f,%.3f) pc=(%.3f,%.3f,%.3f)\n"
            "  |cross|=%.4e  h_normal=(%.4e,%.4e,%.4e)\n"
            "  in_R=%u  boundary=%u  pushed=%d  flipped=%d  not_flipped=%d  stale=%d  constrained_skip=%d\n",
            va, vb, vc,
            pa.x, pa.y, pa.z, pb.x, pb.y, pb.z, pc.x, pc.y, pc.z,
            cross_mag, h_normal.x, h_normal.y, h_normal.z,
            diag_n_R, diag_n_bnd, diag_pushed, diag_flipped, diag_notflipped, diag_stale, diag_constr);
    }

cleanup:
    heap_free(&Q);
    dynarray_free(&boundary);
    dynarray_free(&in_R);
    return ret;
}



/* -----------------------------------------------------------------------
 * Public entry point
 * ----------------------------------------------------------------------- */

#define POINT_IN_TRIMESH_MAX_RETRIES 10
#define CDT_MAX_ITERATIONS 20

s_scplx tetrahedralize_interior_trimesh(const s_trimesh *mesh,
                                         double EPS_DEG, double TOL)
{
    if (!trimesh_is_valid(mesh)) return (s_scplx){0};

    /* Mutable working copy: ridge protection and face Steiner splits accumulate here. */
    s_trimesh working = copy_trimesh(mesh);
    if (!trimesh_is_valid(&working)) return (s_scplx){0};

    s_scplx dt = {0};

    for (int iter = 0; iter < CDT_MAX_ITERATIONS; iter++) {
        fprintf(stderr, "[CDT] iter %d: Phase A: DT + ridge protection (%d verts, %d faces) ...\n",
                iter, working.points.N, working.Nf); fflush(stderr);

        /* ---- Phase A: build DT and ensure all boundary edges are present ---- */
        s_dt_builder b = dt_builder_begin(&working.points, NULL, TOL, NULL, NULL);
        if (!b._stack) { free_trimesh(&working); return (s_scplx){0}; }

        if (!ensure_ridge_protected(&b, &working, TOL)) {
            s_scplx tmp = dt_builder_end(&b, false, NULL, NULL, 0);
            free_complex(&tmp); free_trimesh(&working); return (s_scplx){0};
        }

        /* Strip sentinels; compaction shifts indices by -4, matching working.faces. */
        dt = dt_builder_end(&b, false, NULL, NULL, 0);
        if (!dt.head) { free_trimesh(&working); return (s_scplx){0}; }

        /* Sync working.points with the full compacted DT point set (which includes
         * any ridge Steiner points added by ensure_ridge_protected to the builder).
         * Without this, the next face-Steiner centroid gets an ID that collides with
         * the ridge Steiner IDs already in working.faces. */
        {
            s_point *p = malloc((size_t)dt.points.N * sizeof(s_point));
            if (!p) { free_complex(&dt); free_trimesh(&working); return (s_scplx){0}; }
            memcpy(p, dt.points.p, (size_t)dt.points.N * sizeof(s_point));
            free(working.points.p);
            working.points.p = p;
            working.points.N = dt.points.N;
        }

        /* ---- Phase B: insert boundary faces; collect any that need Steiner ---- */
        fprintf(stderr, "[CDT] iter %d: Phase B: inserting %d boundary faces ...\n",
                iter, working.Nf); fflush(stderr);

        s_hash_table constrained;
        if (!hash_init(&constrained, sizeof(int[3]), sizeof(int),
                       (size_t)(working.Nf * 2 + 1), (size_t)working.Nf,
                       face_triple_hash, face_triple_eq, NULL)) {
            free_complex(&dt); free_trimesh(&working); return (s_scplx){0};
        }

        /* Process faces; stop at the first failure (running more flips on a partially
         * constrained DT corrupts it into degenerate configurations). */
        int failed_va = -1, failed_vb = -1, failed_vc = -1;

        for (int fi = 0; fi < working.Nf; fi++) {
            if (fi % 500 == 0) { fprintf(stderr, "[CDT]   face %d / %d\n", fi, working.Nf); fflush(stderr); }
            int va = working.faces[fi*3], vb = working.faces[fi*3+1], vc = working.faces[fi*3+2];
            s_point pa = dt.points.p[va], pb = dt.points.p[vb], pc = dt.points.p[vc];

            if (!flipinsertfacet(&dt, va, vb, vc, &constrained, pa, pb, pc)) {
                failed_va = va; failed_vb = vb; failed_vc = vc;
                break;
            }
        }

        hash_free(&constrained);

        if (failed_va < 0) {
            break;  /* All faces inserted — done with Phase B */
        }

        fprintf(stderr, "[CDT] iter %d: face (%d,%d,%d) needs Steiner; restarting ...\n",
                iter, failed_va, failed_vb, failed_vc); fflush(stderr);

        /* Discard current DT; rebuild Phase A with Steiner centroid added. */
        free_complex(&dt);
        dt = (s_scplx){0};

        /* Add centroid to working.points and split face → 3 sub-faces. */
        {
            int va = failed_va, vb = failed_vb, vc = failed_vc;

            /* Find the face in working (by vertex triple — vertices never renumber). */
            int fi_found = -1;
            for (int fi = 0; fi < working.Nf; fi++) {
                if (working.faces[fi*3]==va && working.faces[fi*3+1]==vb && working.faces[fi*3+2]==vc) {
                    fi_found = fi; break;
                }
            }
            if (fi_found < 0) {
                /* Face not found — shouldn't happen */
                free_trimesh(&working); return (s_scplx){0};
            }

            /* Add centroid to working.points. */
            s_point pa = working.points.p[va], pb = working.points.p[vb], pc = working.points.p[vc];
            s_point centroid = {
                .x = (pa.x + pb.x + pc.x) * (1.0/3.0),
                .y = (pa.y + pb.y + pc.y) * (1.0/3.0),
                .z = (pa.z + pb.z + pc.z) * (1.0/3.0),
            };
            int new_id = working.points.N;
            s_point *new_pts = realloc(working.points.p, (size_t)(new_id + 1) * sizeof(s_point));
            if (!new_pts) { free_trimesh(&working); return (s_scplx){0}; }
            working.points.p = new_pts;
            working.points.p[new_id] = centroid;
            working.points.N = new_id + 1;

            /* Split face fi_found → 3 sub-faces.
             * Reuse slot fi_found for sub-face 0; append the other two. */
            int new_Nf = working.Nf + 2;
            int *nf = realloc(working.faces, (size_t)new_Nf * 3 * sizeof(int));
            s_point *nn = realloc(working.fnormals, (size_t)new_Nf * sizeof(s_point));
            if (!nf || !nn) { free_trimesh(&working); return (s_scplx){0}; }
            working.faces = nf; working.fnormals = nn;

            s_point fn = working.fnormals[fi_found];
            working.faces[fi_found*3]     = va;  working.faces[fi_found*3+1]     = vb;
            working.faces[fi_found*3+2]   = new_id; working.fnormals[fi_found]   = fn;

            working.faces[working.Nf*3]   = vb;  working.faces[working.Nf*3+1]   = vc;
            working.faces[working.Nf*3+2] = new_id; working.fnormals[working.Nf] = fn;

            working.faces[(working.Nf+1)*3]   = vc;  working.faces[(working.Nf+1)*3+1] = va;
            working.faces[(working.Nf+1)*3+2] = new_id; working.fnormals[working.Nf+1] = fn;

            free(working.adjacency); working.adjacency = NULL;
            working.Nf = new_Nf;
        }
        /* Continue outer loop → rebuild Phase A with extended working mesh */
    }

    free_trimesh(&working);

    if (!dt.head) {
        fprintf(stderr, "[CDT] Failed: could not insert all boundary faces after %d iterations\n",
                CDT_MAX_ITERATIONS);
        return (s_scplx){0};
    }

    fprintf(stderr, "[CDT] Phase C: filtering %d tets ...\n", dt.N_ncells); fflush(stderr);

    /* Ray-cast each tet centroid against the original mesh to classify interior/exterior.
     * For each exterior tet: null back-references in its neighbours (keeps adjacency valid),
     * then unlink and free it. */
    s_ncell *nc = dt.head;
    while (nc) {
        s_ncell *next = nc->next;

        s_point v[4];
        extract_vertices_ncell(&dt, nc, v);
        s_point centroid = {
            .x = (v[0].x + v[1].x + v[2].x + v[3].x) * 0.25,
            .y = (v[0].y + v[1].y + v[2].y + v[3].y) * 0.25,
            .z = (v[0].z + v[1].z + v[2].z + v[3].z) * 0.25,
        };

        if (point_in_trimesh(mesh, centroid, EPS_DEG, POINT_IN_TRIMESH_MAX_RETRIES) == 1) {
            nc = next;
            continue;
        }

        for (int fi = 0; fi < 4; fi++) {
            s_ncell *nb = nc->opposite[fi];
            if (!nb) continue;
            for (int fj = 0; fj < 4; fj++)
                if (nb->opposite[fj] == nc) { nb->opposite[fj] = NULL; break; }
        }
        if (nc->prev) nc->prev->next = next;
        else          dt.head = next;
        if (next)     next->prev = nc->prev;
        dt.N_ncells--;
        free_ncell(nc);
        nc = next;
    }

    fprintf(stderr, "[CDT] Done. Interior tets kept: %d\n", dt.N_ncells); fflush(stderr);

    return dt;
}
