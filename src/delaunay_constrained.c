#include "delaunay.h"
#include "scplx.h"
#include "trimesh.h"
#include "dynarray.h"
#include "hash.h"
#include "gtests.h"
#include "points.h"
#include "voronoi_predicates.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>


static bool edge_in_dt(s_scplx *dt, int va, int vb)
{
    if (!dt->point2tet[va]) return false;

    s_dynarray out_ids = dynarray_initialize(sizeof(int), 16);
    s_dynarray scratch  = dynarray_initialize(sizeof(s_ncell *), 16);
    if (!out_ids.items || !scratch.items) {
        dynarray_free(&out_ids); dynarray_free(&scratch);
        return false;
    }

    int va_lid = -1;
    for (int k = 0; k < 4; k++)
        if (dt->point2tet[va]->vertex_id[k] == va) { va_lid = k; break; }
    vertex_neighbors(dt, va, dt->point2tet[va], va_lid, 0, &out_ids, &scratch);

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

/* Angle at V in triangle V-A-B, in radians. Returns 0 for degenerate inputs. */
static double __attribute__((unused)) interior_angle(s_point v, s_point a, s_point b, double EPS_DEG)
{
    double ux = a.x - v.x, uy = a.y - v.y, uz = a.z - v.z;
    double wx = b.x - v.x, wy = b.y - v.y, wz = b.z - v.z;
    double lu = sqrt(ux*ux + uy*uy + uz*uz);
    double lw = sqrt(wx*wx + wy*wy + wz*wz);
    if (lu <= EPS_DEG || lw <= EPS_DEG) return 0.0;
    double cosA = (ux*wx + uy*wy + uz*wz) / (lu * lw);
    if (cosA >  1.0) cosA =  1.0;
    if (cosA < -1.0) cosA = -1.0;
    return acos(cosA);
}

/* Parameter t such that A + t*(B-A) is the orthogonal projection of V onto
 * line AB. Returns 0.5 for degenerate (zero-length) segments. */
static double __attribute__((unused)) project_t_onto_segment(s_point v, s_point a, s_point b, double EPS_DEG)
{
    double dx = b.x - a.x, dy = b.y - a.y, dz = b.z - a.z;
    double len2 = dx*dx + dy*dy + dz*dz;
    if (len2 <= EPS_DEG) return 0.5;
    return ((v.x - a.x)*dx + (v.y - a.y)*dy + (v.z - a.z)*dz) / len2;
}

/* Return the DT index of the vertex most deeply inside the diametric circumsphere
 * of segment [va_dt, vb_dt], i.e. the one with the largest interior angle to AB.
 * A vertex V is inside the diametric sphere iff dot(A-V, B-V) < 0.
 * Skips sentinels (indices 0-3) and the endpoints themselves.
 * Returns -1 if no encroaching vertex exists. */
static int scout_refpt(const s_scplx *dt, int va_dt, int vb_dt, double EPS_DEG)
{
    s_point A = dt->points.p[va_dt];
    s_point B = dt->points.p[vb_dt];
    int    refpt  = -1;
    double angmax = 0.0;

    for (int i = 4; i < dt->points.N; i++) {
        if (i == va_dt || i == vb_dt) continue;
        s_point V = dt->points.p[i];
        double ax = A.x - V.x, ay = A.y - V.y, az = A.z - V.z;
        double bx = B.x - V.x, by = B.y - V.y, bz = B.z - V.z;
        if (ax*bx + ay*by + az*bz >= 0.0) continue; /* outside diametric sphere */
        double ang = interior_angle(V, A, B, EPS_DEG);
        if (ang > angmax) { angmax = ang; refpt = i; }
    }
    return refpt;
}

/* Compute the Steiner point to insert on segment [va_dt, vb_dt].
 * refpt_id: result of scout_refpt (-1 → fall back to midpoint).
 * va_mesh/vb_mesh: mesh-space ids of the two endpoints (for adjacent-Steiner check).
 * steiner_origins: dynarray of int[2] indexed by DT point id; {-1,-1} = original vertex.
 *
 * Adjacent-Steiner rule (mirrors TetGen getsteinerptonsegment):
 *   if refpt was itself inserted on an edge sharing endpoint A, place the new
 *   point at the same distance from A as refpt; likewise for B.  This prevents
 *   the symmetric spiral where a Steiner near A keeps encroaching on sub-segments. */
static s_point get_steiner_point(const s_scplx *dt,
                                                          int va_dt, int vb_dt,
                                                          int refpt_id,
                                                          int va_mesh, int vb_mesh,
                                                          const s_dynarray *steiner_origins,
                                                          double EPS_DEG)
{
    s_point A = dt->points.p[va_dt];
    s_point B = dt->points.p[vb_dt];
    double t = 0.5;

    if (refpt_id >= 0) {
        s_point R = dt->points.p[refpt_id];

        /* Check if refpt is a Steiner on an adjacent edge sharing A or B. */
        int adj = 0;
        if (refpt_id < (int)steiner_origins->N) {
            const int *orig = (const int *)dynarray_get_ptr(
                                  (s_dynarray *)steiner_origins, refpt_id);
            if (orig[0] != -1) {
                /* orig[] holds the mesh-space edge (pa, pb) refpt was inserted on. */
                int pa = orig[0], pb = orig[1];
                if (pa == va_mesh || pb == va_mesh) {
                    /* shared endpoint is A: mirror refpt's distance from A */
                    double lab = sqrt((B.x-A.x)*(B.x-A.x) +
                                      (B.y-A.y)*(B.y-A.y) +
                                      (B.z-A.z)*(B.z-A.z));
                    double lar = sqrt((R.x-A.x)*(R.x-A.x) +
                                      (R.y-A.y)*(R.y-A.y) +
                                      (R.z-A.z)*(R.z-A.z));
                    t = (lab > EPS_DEG) ? lar / lab : 0.5;
                    adj = 1;
                } else if (pa == vb_mesh || pb == vb_mesh) {
                    /* shared endpoint is B: mirror from B */
                    double lab = sqrt((B.x-A.x)*(B.x-A.x) +
                                      (B.y-A.y)*(B.y-A.y) +
                                      (B.z-A.z)*(B.z-A.z));
                    double lbr = sqrt((R.x-B.x)*(R.x-B.x) +
                                      (R.y-B.y)*(R.y-B.y) +
                                      (R.z-B.z)*(R.z-B.z));
                    t = (lab > EPS_DEG) ? 1.0 - lbr / lab : 0.5;
                    adj = 1;
                }
            }
        }

        if (!adj)
            t = project_t_onto_segment(R, A, B, EPS_DEG);
    }

    /* Clamp to [0.2, 0.8]; fall back to midpoint if outside. */
    if (t < 0.2 || t > 0.8) t = 0.5;

    return (s_point){ .x = A.x + t*(B.x-A.x),
                      .y = A.y + t*(B.y-A.y),
                      .z = A.z + t*(B.z-A.z) };
}

/* Ensure every trimesh boundary edge is present in the DT.
 * Processes ALL unprotected edges per pass (not one at a time) so each pass
 * is O(N_faces) instead of O(N_faces^2) and the total work is proportional
 * to the number of edges that actually need splitting.
 * mesh->faces holds 0-based canonical indices; the +4 builder offset is
 * applied only at the DT lookup call sites here -- the trimesh is never shifted. */
static int ensure_ridge_protected(s_dt_builder *b, s_trimesh *mesh, double EPS_DEG, double TOL)
{
    /* steiner_origins[dt_point_id] = {va, vb} (mesh-space) if the point is a
     * Steiner inserted on edge (va,vb); {-1,-1} for all original vertices.
     * Kept across passes so scout_refpt can identify adjacent Steiner points. */
    s_dynarray steiner_origins = dynarray_initialize(sizeof(int[2]), 64);
    if (!steiner_origins.items) return 0;
    int neg[2] = {-1, -1};
    for (int i = 0; i < b->dt.points.N; i++) {
        if (!dynarray_push(&steiner_origins, neg)) {
            dynarray_free(&steiner_origins); return 0;
        }
    }

    int _pass = 0;
    for (;;) {
        _pass++;

        /* Collect all unique unprotected edges (canonical min<max pairs). */
        s_hash_table edge_set;
        int dummy = 1;
        if (!hash_init(&edge_set, sizeof(int[2]), sizeof(int),
                       (size_t)(mesh->Nf * 4 + 1), (size_t)(mesh->Nf),
                       edge_pair_hash, edge_pair_eq, NULL)) {
            dynarray_free(&steiner_origins); return 0;
        }

        s_dynarray to_split = dynarray_initialize(sizeof(int[2]), 64);
        if (!to_split.items) {
            hash_free(&edge_set);
            dynarray_free(&steiner_origins); return 0;
        }

        for (int fi = 0; fi < mesh->Nf; fi++) {
            for (int ei = 0; ei < 3; ei++) {
                int va = mesh->faces[fi*3 + (ei+1)%3];
                int vb = mesh->faces[fi*3 + (ei+2)%3];
                if (edge_in_dt(&b->dt, va + 4, vb + 4))
                    continue;

                int pair[2] = { va < vb ? va : vb, va < vb ? vb : va };
                if (hash_get(&edge_set, pair)) continue;
                if (!hash_insert(&edge_set, pair, &dummy)) {
                    hash_free(&edge_set); dynarray_free(&to_split);
                    dynarray_free(&steiner_origins); return 0;
                }
                if (!dynarray_push(&to_split, pair)) {
                    hash_free(&edge_set); dynarray_free(&to_split);
                    dynarray_free(&steiner_origins); return 0;
                }
            }
        }

        hash_free(&edge_set);

        if (to_split.N == 0) {
            dynarray_free(&to_split);
            break;
        }

        fprintf(stderr, "[CDT]   ridge pass %d: Nf=%d Nv=%d  %u edges to split\n",
                _pass, mesh->Nf, b->dt.points.N, to_split.N); fflush(stderr);

        /* Insert scout+project Steiner points and split trimesh edges in batch. */
        for (unsigned i = 0; i < to_split.N; i++) {
            int *pair = (int *)dynarray_get_ptr(&to_split, i);
            int va = pair[0], vb = pair[1];

            int refpt_id = scout_refpt(&b->dt, va + 4, vb + 4, EPS_DEG);
            s_point M = get_steiner_point(&b->dt, va + 4, vb + 4,
                                          refpt_id, va, vb,
                                          &steiner_origins, EPS_DEG);

            int new_id_builder = b->dt.points.N;
            s_points single = { .N = 1, .p = &M };
            if (!dt_builder_extend(b, &single, TOL)) {
                dynarray_free(&to_split);
                dynarray_free(&steiner_origins); return 0;
            }

            /* Record origin edge for the newly inserted point (if not deduped). */
            if (b->dt.points.N > (int)steiner_origins.N) {
                int origin[2] = {va, vb};
                if (!dynarray_push(&steiner_origins, origin)) {
                    dynarray_free(&to_split);
                    dynarray_free(&steiner_origins); return 0;
                }
            }

            if (!split_trimesh_edge(mesh, va, vb, new_id_builder - 4)) {
                dynarray_free(&to_split);
                dynarray_free(&steiner_origins); return 0;
            }
        }

        dynarray_free(&to_split);
    }
    dynarray_free(&steiner_origins);
    return 1;
}


/* ----------------------------------------------------------------------- */


static s_ncell *face_in_dt(const s_scplx *dt, int va, int vb, int vc)
{
    s_ncell *start = dt->point2tet[va];
    if (!start) return NULL;
    int va_lid = -1;
    for (int k = 0; k < 4; k++)
        if (start->vertex_id[k] == va) { va_lid = k; break; }
    int comp[3]; for (int i = 0, k = 0; i < 4; i++) if (i != va_lid) comp[k++] = i;

    s_dynarray ring = dynarray_initialize(sizeof(s_ncell *), 16);
    if (!ring.items) return NULL;
    ncells_incident_face((s_scplx *)dt, start, 0, comp, &ring);

    s_ncell *found = NULL;
    for (unsigned i = 0; i < ring.N; i++) {
        s_ncell *nc; dynarray_get_value(&ring, i, &nc);
        bool has_vb = false, has_vc = false;
        for (int j = 0; j < 4; j++) {
            if (nc->vertex_id[j] == vb) has_vb = true;
            if (nc->vertex_id[j] == vc) has_vc = true;
        }
        if (has_vb && has_vc) { found = nc; break; }
    }
    dynarray_free(&ring);
    return found;
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



/* Signed 4x4 in-sphere determinant.
 * Rows are (v-q, |v|^2-|q|^2) for v in {a,b,c,p}; q is the test point.
 * Positive result means q is inside the circumsphere of a,b,c,p (not locally Delaunay). */


/* --- SoS predicates ------------------------------------------------------- */

/* Wraps test_insphere with SoS perturbation.
 * When test_insphere returns 0, iterates through the 5 vertices in ascending
 * global index order.  For the vertex at row r in {a(0),b(1),c(2),L(3),Rv(4)},
 * the SoS sign is (-1)^r * orient3d(4 remaining points in row order). */
static int insphere_sos(s_point a, s_point b, s_point c, s_point L, s_point Rv,
                        int ia, int ib, int ic, int iL, int iRv)
{
    int s = test_insphere((s_point[4]){a, b, c, L}, Rv);
    if (s != 0) return s;

    typedef struct { int idx; int row; } ie;
    ie order[5] = {{ia,0},{ib,1},{ic,2},{iL,3},{iRv,4}};
    for (int i = 1; i < 5; i++) {
        ie tmp = order[i]; int j = i - 1;
        while (j >= 0 && order[j].idx > tmp.idx) { order[j+1] = order[j]; j--; }
        order[j+1] = tmp;
    }

    s_point pts[5] = {a, b, c, L, Rv};
    for (int k = 0; k < 5; k++) {
        int row = order[k].row;
        s_point rem[4]; int ri = 0;
        for (int r = 0; r < 5; r++) if (r != row) rem[ri++] = pts[r];
        int o = test_orientation((s_point[3]){rem[0], rem[1], rem[2]}, rem[3]);
        if (o != 0) return (row % 2 == 0 ? 1 : -1) * o;
    }
    return 0; /* all five points coincident — degenerate input */
}

/* Wraps cdt_dprime with SoS perturbation.
 * D' has rows for {a(0),b(1),c(2),L(3)} with Rv as the reference; Rv has no
 * row of its own.  When cdt_dprime returns 0, iterates {a,b,c,L} in ascending
 * global index order and evaluates the cofactor of the dot-product column:
 *   (-1)^(row+3) * orient3d(remaining 3 rows, Rv). */
static int dprime_sos(s_point a, s_point b, s_point c, s_point L, s_point Rv,
                      int ia, int ib, int ic, int iL,
                      s_point pa, s_point pb, s_point pc)
{
    int s = cdt_dprime(a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z,
                       L.x, L.y, L.z, Rv.x, Rv.y, Rv.z,
                       pa.x, pa.y, pa.z, pb.x, pb.y, pb.z, pc.x, pc.y, pc.z);
    if (s != 0) return s;

    typedef struct { int idx; int row; } ie;
    ie order[4] = {{ia,0},{ib,1},{ic,2},{iL,3}};
    for (int i = 1; i < 4; i++) {
        ie tmp = order[i]; int j = i - 1;
        while (j >= 0 && order[j].idx > tmp.idx) { order[j+1] = order[j]; j--; }
        order[j+1] = tmp;
    }

    s_point pts[4] = {a, b, c, L};
    for (int k = 0; k < 4; k++) {
        int row = order[k].row;
        s_point rem[3]; int ri = 0;
        for (int r = 0; r < 4; r++) if (r != row) rem[ri++] = pts[r];
        int o = test_orientation((s_point[3]){rem[0], rem[1], rem[2]}, Rv);
        /* (-1)^(row+3) = (-1)^(row+1): alternates -,+,-,+ starting at row 0 */
        if (o != 0) return (row % 2 == 0 ? -1 : 1) * o;
    }
    return 0; /* all five points coplanar — degenerate input */
}

#define CERTIFY_NO_EVENT (-1)
#define CERTIFY_EVENT      0

/* Returns CERTIFY_EVENT if face loses regularity at τ >= 0 and fills *out_L_id,
 * *out_Rv_id.  Returns CERTIFY_NO_EVENT otherwise. */
static int certify(const s_scplx *dt,
                   const s_ncell *s, int opp_lid,
                   s_point pa, s_point pb, s_point pc,
                   int *out_L_id, int *out_Rv_id)
{
    s_ncell *t = s->opposite[opp_lid];
    if (!t) return CERTIFY_NO_EVENT;

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
    if (Rv_id < 0) return CERTIFY_NO_EVENT;

    s_point Rv = dt->points.p[Rv_id];

    int s_D0     = insphere_sos(a, b, c, L, Rv,
                                face_ids[0], face_ids[1], face_ids[2], L_id, Rv_id);
    int s_Dprime = dprime_sos(a, b, c, L, Rv,
                              face_ids[0], face_ids[1], face_ids[2], L_id,
                              pa, pb, pc);

    if (s_D0 == 0 || s_Dprime == 0) return CERTIFY_NO_EVENT;

    if (-s_D0 * s_Dprime < 0) return CERTIFY_NO_EVENT;  /* τ < 0 */

    *out_L_id  = L_id;
    *out_Rv_id = Rv_id;
    return CERTIFY_EVENT;
}

/* --- Min-heap on tau (Step 3) ------------------------------------------- */

typedef struct {
    int      va, vb, vc;  /* sorted face triple */
    int      L_id;        /* opposite vertex in non-R tet s */
    int      Rv_id;       /* opposite vertex in R tet t */
    s_ncell *s;           /* tet NOT in R at push time; re-verified on pop */
    int      opp_lid;
} heap_entry;

typedef struct {
    s_dynarray     data;
    const s_scplx *dt;
    s_point        pa, pb, pc;
} s_heap;

static s_heap heap_init(const s_scplx *dt, s_point pa, s_point pb, s_point pc)
{
    return (s_heap){ .data = dynarray_initialize(sizeof(heap_entry), 32),
                     .dt = dt, .pa = pa, .pb = pb, .pc = pc };
}

/* Returns true if event A has a smaller flip time than event B.
 * τ_A < τ_B  ⟺  sign(D′_A)·sign(D′_B)·sign(D0_B·D′_A − D0_A·D′_B) < 0 */
/* SoS wrapper for sign(D0_B·D'_A − D0_A·D'_B).
 * When the unweighted predicate returns 0, iterates all vertex IDs from both
 * events in ascending order, perturbing one weight at a time until nonzero. */
static int sign_cross_dprime_sos(
    s_point aA, s_point bA, s_point cA, s_point LA, s_point RvA,
    int iaA, int ibA, int icA, int iLA, int iRvA,
    s_point aB, s_point bB, s_point cB, s_point LB, s_point RvB,
    int iaB, int ibB, int icB, int iLB, int iRvB,
    s_point pa, s_point pb, s_point pc)
{
    int s = cdt_sign_cross_dprime(
        aA.x, aA.y, aA.z, bA.x, bA.y, bA.z, cA.x, cA.y, cA.z,
        LA.x, LA.y, LA.z, RvA.x, RvA.y, RvA.z,
        aB.x, aB.y, aB.z, bB.x, bB.y, bB.z, cB.x, cB.y, cB.z,
        LB.x, LB.y, LB.z, RvB.x, RvB.y, RvB.z,
        pa.x, pa.y, pa.z, pb.x, pb.y, pb.z, pc.x, pc.y, pc.z);
    if (s != 0) return s;

    /* collect, sort and deduplicate all vertex IDs from both events */
    int ids[10] = { iaA, ibA, icA, iLA, iRvA, iaB, ibB, icB, iLB, iRvB };
    for (int i = 1; i < 10; i++) {
        int tmp = ids[i], j = i - 1;
        while (j >= 0 && ids[j] > tmp) { ids[j+1] = ids[j]; j--; }
        ids[j+1] = tmp;
    }
    int unique[10]; int n = 0;
    for (int i = 0; i < 10; i++)
        if (n == 0 || unique[n-1] != ids[i]) unique[n++] = ids[i];

    const int idsA[5] = { iaA, ibA, icA, iLA, iRvA };
    const int idsB[5] = { iaB, ibB, icB, iLB, iRvB };

    for (int k = 0; k < n; k++) {
        int v = unique[k];
        double kA[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        double kB[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
        for (int i = 0; i < 5; i++) {
            if (idsA[i] == v) kA[i] = 1.0;
            if (idsB[i] == v) kB[i] = 1.0;
        }
        int sv = cdt_sign_cross_dprime_weighted(
            aA.x, aA.y, aA.z, bA.x, bA.y, bA.z, cA.x, cA.y, cA.z,
            LA.x, LA.y, LA.z, RvA.x, RvA.y, RvA.z,
            aB.x, aB.y, aB.z, bB.x, bB.y, bB.z, cB.x, cB.y, cB.z,
            LB.x, LB.y, LB.z, RvB.x, RvB.y, RvB.z,
            pa.x, pa.y, pa.z, pb.x, pb.y, pb.z, pc.x, pc.y, pc.z,
            kA, kB);
        if (sv != 0) return sv;
    }
    return 0; /* fully degenerate — shouldn't happen in practice */
}

static bool tau_lt(const heap_entry *A, const heap_entry *B,
                   const s_scplx *dt, s_point pa, s_point pb, s_point pc)
{
    s_point aA  = dt->points.p[A->va],  bA  = dt->points.p[A->vb],
            cA  = dt->points.p[A->vc],  LA  = dt->points.p[A->L_id],
            RvA = dt->points.p[A->Rv_id];
    s_point aB  = dt->points.p[B->va],  bB  = dt->points.p[B->vb],
            cB  = dt->points.p[B->vc],  LB  = dt->points.p[B->L_id],
            RvB = dt->points.p[B->Rv_id];

    int raw_dpA = cdt_dprime(aA.x, aA.y, aA.z, bA.x, bA.y, bA.z, cA.x, cA.y, cA.z,
                             LA.x, LA.y, LA.z, RvA.x, RvA.y, RvA.z,
                             pa.x, pa.y, pa.z, pb.x, pb.y, pb.z, pc.x, pc.y, pc.z);
    int raw_dpB = cdt_dprime(aB.x, aB.y, aB.z, bB.x, bB.y, bB.z, cB.x, cB.y, cB.z,
                             LB.x, LB.y, LB.z, RvB.x, RvB.y, RvB.z,
                             pa.x, pa.y, pa.z, pb.x, pb.y, pb.z, pc.x, pc.y, pc.z);

    /* D'=0 means τ=∞: those events go last */
    if (raw_dpA == 0 && raw_dpB != 0) return false;
    if (raw_dpA != 0 && raw_dpB == 0) return true;
    if (raw_dpA == 0 && raw_dpB == 0) {
        if (A->Rv_id != B->Rv_id) return A->Rv_id < B->Rv_id;
        if (A->L_id  != B->L_id)  return A->L_id  < B->L_id;
        return false;
    }

    int sDpA = dprime_sos(aA, bA, cA, LA, RvA, A->va, A->vb, A->vc, A->L_id, pa, pb, pc);
    int sDpB = dprime_sos(aB, bB, cB, LB, RvB, B->va, B->vb, B->vc, B->L_id, pa, pb, pc);
    int sP   = sign_cross_dprime_sos(
        aA, bA, cA, LA, RvA, A->va, A->vb, A->vc, A->L_id, A->Rv_id,
        aB, bB, cB, LB, RvB, B->va, B->vb, B->vc, B->L_id, B->Rv_id,
        pa, pb, pc);

    return sDpA * sDpB * sP < 0;
}

static void heap_free(s_heap *h) { dynarray_free(&h->data); }
static int  heap_empty(const s_heap *h) { return h->data.N == 0; }

static int heap_push(s_heap *h, int va, int vb, int vc,
                     int L_id, int Rv_id, s_ncell *s, int opp_lid)
{
    heap_entry e = { va, vb, vc, L_id, Rv_id, s, opp_lid };
    if (!dynarray_push(&h->data, &e)) return 0;

    unsigned i = h->data.N - 1;
    while (i > 0) {
        unsigned parent = (i - 1) / 2;
        heap_entry *ei = dynarray_get_ptr(&h->data, i);
        heap_entry *ep = dynarray_get_ptr(&h->data, parent);
        if (!tau_lt(ei, ep, h->dt, h->pa, h->pb, h->pc)) break;
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
            if (l < h->data.N) {
                heap_entry *el = dynarray_get_ptr(&h->data, l);
                heap_entry *es = dynarray_get_ptr(&h->data, sm);
                if (tau_lt(el, es, h->dt, h->pa, h->pb, h->pc)) sm = l;
            }
            if (r < h->data.N) {
                heap_entry *er = dynarray_get_ptr(&h->data, r);
                heap_entry *es = dynarray_get_ptr(&h->data, sm);
                if (tau_lt(er, es, h->dt, h->pa, h->pb, h->pc)) sm = r;
            }
            if (sm == i) break;
            heap_entry *ei = dynarray_get_ptr(&h->data, i);
            heap_entry *es = dynarray_get_ptr(&h->data, sm);
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
    s_point tri[3] = {pa, pb, pc};
    s_point tet[4];
    for (int j = 0; j < 4; j++) tet[j] = dt->points.p[nc->vertex_id[j]];
    return test_triangle_tetrahedron_intersect(tri, tet, 0.0)
           == INTERSECT_NONDEGENERATE;
}

/* --- find_region_R_boundary (Step 4) ------------------------------------ */

/* BFS from a seed R-tet (straddling the constraint plane). Collects:
 *   boundary: sorted int[3] face triples between an R-tet and a non-R-tet
 *   in_R:     s_ncell* pointers for every tet in R
 * Constrained faces are excluded from boundary and not crossed during BFS.
 * Returns 1 on success, 0 on allocation failure. */
static int find_region_R_boundary(s_scplx *dt,
                                   int va,
                                   s_point pa, s_point pb, s_point pc,
                                   s_hash_table *constrained,
                                   s_dynarray *boundary,
                                   s_dynarray *in_R)
{
    /* Find any starting tet that straddles the constraint plane.
     * Ring-walk from point2tet[va] first; fall back to head scan if needed. */
    s_ncell *seed = NULL;
    if (dt->point2tet[va]) {
        s_ncell *start = dt->point2tet[va];
        int va_lid = -1;
        for (int k = 0; k < 4; k++)
            if (start->vertex_id[k] == va) { va_lid = k; break; }
        int comp[3]; for (int i = 0, k = 0; i < 4; i++) if (i != va_lid) comp[k++] = i;

        s_dynarray ring = dynarray_initialize(sizeof(s_ncell *), 16);
        if (ring.items) {
            ncells_incident_face(dt, start, 0, comp, &ring);
            for (unsigned i = 0; i < ring.N && !seed; i++) {
                s_ncell *nc; dynarray_get_value(&ring, i, &nc);
                if (tet_in_R(dt, nc, pa, pb, pc)) seed = nc;
            }
            dynarray_free(&ring);
        }
    }
    if (!seed) {
        for (s_ncell *nc = dt->head; nc; nc = nc->next) {
            if (tet_in_R(dt, nc, pa, pb, pc)) { seed = nc; break; }
        }
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
    if (!dt->point2tet[fa]) return 0;
    s_ncell *start = dt->point2tet[fa];
    int fa_lid = -1;
    for (int k = 0; k < 4; k++)
        if (start->vertex_id[k] == fa) { fa_lid = k; break; }
    int comp[3]; for (int i = 0, k = 0; i < 4; i++) if (i != fa_lid) comp[k++] = i;

    s_dynarray ring = dynarray_initialize(sizeof(s_ncell *), 16);
    if (!ring.items) return 0;
    ncells_incident_face((s_scplx *)dt, start, 0, comp, &ring);

    int result = 0;
    for (unsigned i = 0; i < ring.N; i++) {
        s_ncell *nc; dynarray_get_value(&ring, i, &nc);
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
        if (!nb) break;

        bool nc_in_R = tet_in_R(dt, nc, pa, pb, pc);
        bool nb_in_R = tet_in_R(dt, nb, pa, pb, pc);
        if (nc_in_R == nb_in_R) break;  /* no longer an R-boundary face */

        if (!nc_in_R) {
            *out_s = nc; *out_opp_lid = opp;
        } else {
            int nb_opp = -1;
            for (int j = 0; j < 4; j++) {
                int v = nb->vertex_id[j];
                if (v != fa && v != fb && v != fc) { nb_opp = j; break; }
            }
            if (nb_opp < 0) break;
            *out_s = nb; *out_opp_lid = nb_opp;
        }
        result = 1;
        break;
    }
    dynarray_free(&ring);
    return result;
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

    int ret = 0;
    s_dynarray boundary = {0};
    s_dynarray in_R     = {0};
    s_heap Q            = {0};

    boundary = dynarray_initialize(sizeof(int[3]), 32);
    in_R     = dynarray_initialize(sizeof(s_ncell *), 64);
    if (!boundary.items || !in_R.items) goto cleanup;

    if (!find_region_R_boundary(dt, va, pa, pb, pc, constrained, &boundary, &in_R))
        goto cleanup;

    unsigned diag_n_R = in_R.N, diag_n_bnd = boundary.N;

    Q = heap_init(dt, pa, pb, pc);
    if (!Q.data.items) goto cleanup;

    int diag_pushed = 0;
    for (unsigned i = 0; i < boundary.N; i++) {
        int *fids = dynarray_get_ptr(&boundary, i);
        s_ncell *s; int opp_lid;
        if (!find_boundary_face_s(dt, pa, pb, pc, fids[0], fids[1], fids[2], &s, &opp_lid))
            continue;
        int L_id = -1, Rv_id = -1;
        if (certify(dt, s, opp_lid, pa, pb, pc, &L_id, &Rv_id) == CERTIFY_EVENT)
            { heap_push(&Q, fids[0], fids[1], fids[2], L_id, Rv_id, s, opp_lid); diag_pushed++; }
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

                int L_id = -1, Rv_id = -1;
                if (certify(dt, nb, nb_opp, pa, pb, pc, &L_id, &Rv_id) == CERTIFY_EVENT)
                    heap_push(&Q, fids[0], fids[1], fids[2], L_id, Rv_id, nb, nb_opp);
            }
        }
    }

    if (face_in_dt(dt, va, vb, vc)) {
        mark_constrained_face(constrained, va, vb, vc);
        ret = 1;
    } else {
        (void)diag_n_R; (void)diag_n_bnd; (void)diag_pushed; (void)diag_flipped; (void)diag_notflipped; 
        (void)diag_stale; (void)diag_constr;
        // s_point cross = cross_prod(subtract_points(pb, pa), subtract_points(pc, pa));
        // double cross_mag = sqrt(cross.x*cross.x + cross.y*cross.y + cross.z*cross.z);
        // fprintf(stderr,
            // "flipinsertfacet: heap exhausted for face (%d,%d,%d)\n"
            // "  pa=(%.3f,%.3f,%.3f) pb=(%.3f,%.3f,%.3f) pc=(%.3f,%.3f,%.3f)\n"
            // "  |cross|=%.4e\n"
            // "  in_R=%u  boundary=%u  pushed=%d  flipped=%d  not_flipped=%d  stale=%d  constrained_skip=%d\n",
            // va, vb, vc,
            // pa.x, pa.y, pa.z, pb.x, pb.y, pb.z, pc.x, pc.y, pc.z,
            // cross_mag,
            // diag_n_R, diag_n_bnd, diag_pushed, diag_flipped, diag_notflipped, diag_stale, diag_constr);
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
s_scplx tetrahedralize_interior_trimesh(const s_trimesh *mesh,
                                         double EPS_DEG, double TOL)
{
    if (!trimesh_is_valid(mesh)) return (s_scplx){0};

    s_trimesh working = copy_trimesh(mesh);
    if (!trimesh_is_valid(&working)) return (s_scplx){0};

    s_scplx dt = {0};

#define CDT_MAX_ITERATIONS 64
    for (int iter = 0; iter < CDT_MAX_ITERATIONS; iter++) {
        /* ---- Phase A: build DT and ensure all boundary edges are present ---- */
        fprintf(stderr, "[CDT] iter %d: Phase A: DT + ridge protection (%d verts, %d faces) ...\n",
                iter, working.points.N, working.Nf); fflush(stderr);

        s_dt_builder b = dt_builder_begin(&working.points, NULL, TOL, NULL, NULL);
        if (!b._stack) { free_trimesh(&working); return (s_scplx){0}; }

        if (!ensure_ridge_protected(&b, &working, EPS_DEG, TOL)) {
            s_scplx tmp = dt_builder_end(&b, false, NULL, NULL, 0);
            free_complex(&tmp); free_trimesh(&working); return (s_scplx){0};
        }

        dt = dt_builder_end(&b, false, NULL, NULL, 0);
        if (!dt.head) { free_trimesh(&working); return (s_scplx){0}; }

        /* Sync working.points with the compacted DT point set (includes ridge Steiners). */
        {
            s_point *p = malloc((size_t)dt.points.N * sizeof(s_point));
            if (!p) { free_complex(&dt); free_trimesh(&working); return (s_scplx){0}; }
            memcpy(p, dt.points.p, (size_t)dt.points.N * sizeof(s_point));
            free(working.points.p);
            working.points.p = p;
            working.points.N = dt.points.N;
        }

        /* ---- Phase B: FIFO queue of faces ------------------------------------ */
        fprintf(stderr, "[CDT] iter %d: Phase B: inserting %d boundary faces ...\n",
                iter, working.Nf); fflush(stderr);

        s_hash_table constrained;
        if (!hash_init(&constrained, sizeof(int[3]), sizeof(int),
                       (size_t)(working.Nf * 4 + 1), (size_t)(working.Nf * 2),
                       face_triple_hash, face_triple_eq, NULL)) {
            free_complex(&dt); free_trimesh(&working); return (s_scplx){0};
        }

        s_dynarray face_queue = dynarray_initialize(sizeof(int[3]), (size_t)working.Nf + 64);
        if (!face_queue.items) {
            hash_free(&constrained); free_complex(&dt); free_trimesh(&working); return (s_scplx){0};
        }

        /* Pre-mark faces already present in the DT (from previous iterations). */
        for (int fi = 0; fi < working.Nf; fi++) {
            int va = working.faces[fi*3], vb = working.faces[fi*3+1], vc = working.faces[fi*3+2];
            if (face_in_dt(&dt, va, vb, vc)) {
                mark_constrained_face(&constrained, va, vb, vc);
            } else {
                int triple[3] = {va, vb, vc};
                if (!dynarray_push(&face_queue, triple)) {
                    dynarray_free(&face_queue); hash_free(&constrained);
                    free_complex(&dt); free_trimesh(&working); return (s_scplx){0};
                }
            }
        }

        int failed_va = -1, failed_vb = -1, failed_vc = -1;
        size_t front = 0;

        while (front < face_queue.N) {
            int triple[3];
            memcpy(triple, dynarray_get_ptr(&face_queue, front++), sizeof(triple));
            int va = triple[0], vb = triple[1], vc = triple[2];

            if (is_constrained_face(&constrained, va, vb, vc)) continue;

            if ((front - 1) % 500 == 0)
                fprintf(stderr, "[CDT]   face %u / %u\n",
                        (unsigned)(front - 1), (unsigned)face_queue.N);

            s_point pa = dt.points.p[va], pb = dt.points.p[vb], pc = dt.points.p[vc];

            if (flipinsertfacet(&dt, va, vb, vc, &constrained, pa, pb, pc)) continue;

            failed_va = va; failed_vb = vb; failed_vc = vc;
            break;
        }

        dynarray_free(&face_queue);
        hash_free(&constrained);

        if (failed_va < 0) break; /* All faces inserted — done */

        fprintf(stderr, "[CDT] iter %d: face (%d,%d,%d) needs Steiner; restarting ...\n",
                iter, failed_va, failed_vb, failed_vc); fflush(stderr);

        /* Discard DT; add centroid Steiner to working and split the failed face → 3. */
        free_complex(&dt);
        dt = (s_scplx){0};

        {
            int va = failed_va, vb = failed_vb, vc = failed_vc;

            int fi_found = -1;
            for (int fi = 0; fi < working.Nf; fi++) {
                if (working.faces[fi*3]==va && working.faces[fi*3+1]==vb && working.faces[fi*3+2]==vc) {
                    fi_found = fi; break;
                }
            }
            if (fi_found < 0) { free_trimesh(&working); return (s_scplx){0}; }

            s_point pa = working.points.p[va], pb = working.points.p[vb], pc = working.points.p[vc];
            s_point centroid = { .x = (pa.x+pb.x+pc.x)*(1.0/3.0),
                                 .y = (pa.y+pb.y+pc.y)*(1.0/3.0),
                                 .z = (pa.z+pb.z+pc.z)*(1.0/3.0) };
            int new_id = working.points.N;
            s_point *new_pts = realloc(working.points.p, (size_t)(new_id+1) * sizeof(s_point));
            if (!new_pts) { free_trimesh(&working); return (s_scplx){0}; }
            working.points.p = new_pts;
            working.points.p[new_id] = centroid;
            working.points.N = new_id + 1;

            int new_Nf = working.Nf + 2;
            int     *nf = realloc(working.faces,    (size_t)new_Nf * 3 * sizeof(int));
            s_point *nn = realloc(working.fnormals, (size_t)new_Nf     * sizeof(s_point));
            if (!nf || !nn) { free_trimesh(&working); return (s_scplx){0}; }
            working.faces = nf; working.fnormals = nn;

            s_point fn = working.fnormals[fi_found];
            working.faces[fi_found*3]       = va;     working.faces[fi_found*3+1]       = vb;
            working.faces[fi_found*3+2]     = new_id; working.fnormals[fi_found]        = fn;
            working.faces[working.Nf*3]     = vb;     working.faces[working.Nf*3+1]     = vc;
            working.faces[working.Nf*3+2]   = new_id; working.fnormals[working.Nf]      = fn;
            working.faces[(working.Nf+1)*3]   = vc;   working.faces[(working.Nf+1)*3+1] = va;
            working.faces[(working.Nf+1)*3+2] = new_id; working.fnormals[working.Nf+1]  = fn;
            free(working.adjacency); working.adjacency = NULL;
            working.Nf = new_Nf;
        }
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
