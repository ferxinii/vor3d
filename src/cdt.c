// TODO: the reference uses isTetLocallyDelaunay for condition iii (circumsphere/
// locally-Delaunay check). It differs from our condition iii in that it determines
// vertex visibility against the original initial boundary (original_C_bnd), not the
// current active faces — this is a more accurate occlusion test. But the current
// failures are all FAIL(ii), not FAIL(iii).

#include "delaunay.h"
#include "scplx.h"
#include "trimesh.h"
#include "dynarray.h"
#include "hash.h"
#include "gtests.h"
#include "points.h"
#include "voronoi_predicates.h"
#include "cdt_predicates.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>

/* Per-point LNC info. v1==-1 means explicit (original) vertex.
 * For Steiners: represents t*V[v1] + (1-t)*V[v2], v1/v2 are DT indices. */
typedef struct {
    int    v1, v2;
    double t;
} s_lnc_info;

/* Half-cavity from splitting the region Tf by the constraint plane.
 * All indices are DT indices. */
typedef struct {
    s_dynarray tets;      /* s_ncell* — tets with ≥1 vertex on this side */
    s_dynarray verts;     /* int     — DT indices of vertices on/above (C1) or on/below (C2) */
    s_dynarray boundary;  /* int[4]  — {sorted face triple, interior_v} of ∂Ci */
} s_half_cavity;

static void half_cavity_free(s_half_cavity *c)
{
    dynarray_free(&c->tets);
    dynarray_free(&c->verts);
    dynarray_free(&c->boundary);
}

static size_t ncell_ptr_hash(const void *key)
{
    uintptr_t p; memcpy(&p, key, sizeof(p));
    p ^= p >> 16; p *= (uintptr_t)0x45d9f3b3u; p ^= p >> 16;
    return (size_t)p;
}

static bool ncell_ptr_eq(const void *a, const void *b)
{
    return memcmp(a, b, sizeof(s_ncell *)) == 0;
}

static size_t int_hash(const void *key)
{
    unsigned v; memcpy(&v, key, sizeof(v));
    return (size_t)(v * 2654435761u);
}

static bool int_eq(const void *a, const void *b)
{
    return memcmp(a, b, sizeof(int)) == 0;
}


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
static double interior_angle(s_point v, s_point a, s_point b, double EPS_DEG)
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
static double project_t_onto_segment(s_point v, s_point a, s_point b, double EPS_DEG)
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

/* Compute the parameter t in (0,1) for a Steiner on segment [va_dt, vb_dt].
 * Paper convention: position = t*V[va_dt] + (1-t)*V[vb_dt].
 * refpt_id: result of scout_refpt (-1 → midpoint).
 * lnc_array/lnc_n: per-point LNC info indexed by DT id; v1==-1 means explicit.
 *
 * Adjacent-Steiner rule: if refpt was itself a Steiner on an edge sharing
 * endpoint A or B, mirror its distance from that endpoint to avoid spirals. */
static double compute_steiner_t(const s_scplx *dt,
                                 int va_dt, int vb_dt,
                                 int refpt_id,
                                 const s_lnc_info *lnc_array, int lnc_n,
                                 double EPS_DEG)
{
    s_point A = dt->points.p[va_dt];
    s_point B = dt->points.p[vb_dt];
    double t = 0.5;

    if (refpt_id >= 0) {
        s_point R = dt->points.p[refpt_id];
        int adj = 0;

        if (refpt_id < lnc_n && lnc_array[refpt_id].v1 != -1) {
            int pa = lnc_array[refpt_id].v1;  /* DT indices of refpt's segment */
            int pb = lnc_array[refpt_id].v2;
            if (pa == va_dt || pb == va_dt) {
                double lab = sqrt((B.x-A.x)*(B.x-A.x) +
                                  (B.y-A.y)*(B.y-A.y) +
                                  (B.z-A.z)*(B.z-A.z));
                double lar = sqrt((R.x-A.x)*(R.x-A.x) +
                                  (R.y-A.y)*(R.y-A.y) +
                                  (R.z-A.z)*(R.z-A.z));
                t = (lab > EPS_DEG) ? lar / lab : 0.5;
                adj = 1;
            } else if (pa == vb_dt || pb == vb_dt) {
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

        if (!adj)
            t = project_t_onto_segment(R, A, B, EPS_DEG);
    }

    if (t < 0.2 || t > 0.8) t = 0.5;
    return t;
}

/* Ensure every trimesh boundary edge is present in the DT.
 * Processes ALL unprotected edges per pass (not one at a time) so each pass
 * is O(N_faces) instead of O(N_faces^2) and the total work is proportional
 * to the number of edges that actually need splitting.
 * mesh->faces holds 0-based canonical indices; the +4 builder offset is
 * applied only at the DT lookup call sites here -- the trimesh is never shifted. */
static int ensure_ridge_protected(s_dt_builder *b, s_trimesh *mesh, double EPS_DEG, double TOL)
{
    /* Register all initial DT points (including sentinels 0-3) as explicit. */
    cdt_predicates_init();
    for (int i = 0; i < b->dt.points.N; i++) {
        s_point *p = &b->dt.points.p[i];
        cdt_point_set_explicit(i, p->x, p->y, p->z);
    }

    /* lnc_info[dt_point_id]: v1==-1 for explicit vertices, DT indices for Steiners. */
    s_dynarray lnc_info = dynarray_initialize(sizeof(s_lnc_info), 64);
    if (!lnc_info.items) return 0;
    s_lnc_info explicit_entry = { .v1 = -1, .v2 = 0, .t = 0.0 };
    for (int i = 0; i < b->dt.points.N; i++) {
        if (!dynarray_push(&lnc_info, &explicit_entry)) {
            dynarray_free(&lnc_info); return 0;
        }
    }

    for (;;) {
        /* Collect all unique unprotected edges (canonical min<max pairs). */
        s_hash_table edge_set;
        int dummy = 1;
        if (!hash_init(&edge_set, sizeof(int[2]), sizeof(int),
                       (size_t)(mesh->Nf * 4 + 1), (size_t)(mesh->Nf),
                       edge_pair_hash, edge_pair_eq, NULL)) {
            dynarray_free(&lnc_info); return 0;
        }

        s_dynarray to_split = dynarray_initialize(sizeof(int[2]), 64);
        if (!to_split.items) {
            hash_free(&edge_set);
            dynarray_free(&lnc_info); return 0;
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
                    dynarray_free(&lnc_info); return 0;
                }
                if (!dynarray_push(&to_split, pair)) {
                    hash_free(&edge_set); dynarray_free(&to_split);
                    dynarray_free(&lnc_info); return 0;
                }
            }
        }

        hash_free(&edge_set);

        if (to_split.N == 0) {
            dynarray_free(&to_split);
            break;
        }

        for (unsigned i = 0; i < to_split.N; i++) {
            int *pair = (int *)dynarray_get_ptr(&to_split, i);
            int va = pair[0], vb = pair[1];
            int va_dt = va + 4, vb_dt = vb + 4;

            int refpt_id = scout_refpt(&b->dt, va_dt, vb_dt, EPS_DEG);
            double t = compute_steiner_t(&b->dt, va_dt, vb_dt, refpt_id,
                                          (const s_lnc_info *)lnc_info.items,
                                          (int)lnc_info.N, EPS_DEG);

            s_point A = b->dt.points.p[va_dt];
            s_point B = b->dt.points.p[vb_dt];
            s_point M = { .x = A.x + t*(B.x - A.x),
                          .y = A.y + t*(B.y - A.y),
                          .z = A.z + t*(B.z - A.z) };

            int new_id = b->dt.points.N;
            s_points single = { .N = 1, .p = &M };
            if (!dt_builder_extend(b, &single, TOL)) {
                dynarray_free(&to_split);
                dynarray_free(&lnc_info); return 0;
            }

            /* If the point was not deduped, register it as LNC. */
            if (b->dt.points.N > (int)lnc_info.N) {
                cdt_point_set_lnc(new_id, va_dt, vb_dt, t);
                s_lnc_info info = { .v1 = va_dt, .v2 = vb_dt, .t = t };
                if (!dynarray_push(&lnc_info, &info)) {
                    dynarray_free(&to_split);
                    dynarray_free(&lnc_info); return 0;
                }
            }

            if (!split_trimesh_edge(mesh, va, vb, new_id - 4)) {
                dynarray_free(&to_split);
                dynarray_free(&lnc_info); return 0;
            }
        }

        dynarray_free(&to_split);
    }
    dynarray_free(&lnc_info);
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
    if (va_lid < 0) return NULL;
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


/* ----------------------------------------------------------------------- */

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
static int __attribute__((unused)) find_boundary_face_s(const s_scplx *dt,
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

/* Split Tf by the constraint plane, producing two half-cavities C1 (above)
 * and C2 (below).  va_dt/vb_dt/vc_dt are DT indices of the constraint face.
 * C1->tets: Tf tets with ≥1 vertex strictly above.
 * C2->tets: Tf tets with ≥1 vertex strictly below.
 * C1->verts: DT indices of vertices on or above the plane.
 * C2->verts: DT indices of vertices on or below the plane.
 * C1->boundary / C2->boundary: sorted int[3] face triples of ∂Ci, including
 * the constraint face itself. */
static int build_half_cavities(s_scplx *dt,
                                int va_dt, int vb_dt, int vc_dt,
                                s_hash_table *constrained,
                                s_half_cavity *C1, s_half_cavity *C2)
{
    s_point pa = dt->points.p[va_dt];
    s_point pb = dt->points.p[vb_dt];
    s_point pc = dt->points.p[vc_dt];

    s_dynarray in_R     = dynarray_initialize(sizeof(s_ncell *), 64);
    s_dynarray boundary = dynarray_initialize(sizeof(int[3]), 32);
    if (!in_R.items || !boundary.items) {
        dynarray_free(&in_R); dynarray_free(&boundary); return 0;
    }

    if (!find_region_R_boundary(dt, va_dt, pa, pb, pc, constrained, &boundary, &in_R)) {
        dynarray_free(&in_R); dynarray_free(&boundary); return 0;
    }
    dynarray_free(&boundary);

    C1->tets     = dynarray_initialize(sizeof(s_ncell *), 32);
    C1->verts    = dynarray_initialize(sizeof(int), 32);
    C1->boundary = dynarray_initialize(sizeof(int[4]), 32);
    C2->tets     = dynarray_initialize(sizeof(s_ncell *), 32);
    C2->verts    = dynarray_initialize(sizeof(int), 32);
    C2->boundary = dynarray_initialize(sizeof(int[4]), 32);
    if (!C1->tets.items || !C1->verts.items || !C1->boundary.items ||
        !C2->tets.items || !C2->verts.items || !C2->boundary.items) {
        dynarray_free(&in_R); return 0;
    }

    /* Build a membership set for ALL of R so we can identify the true R-vs-nonR
     * boundary faces.  The old code used per-half sets (c1_set/c2_set), which
     * caused faces between a straddling tet and a purely-opposite-side R-tet to
     * be treated as boundary faces — they are interior to R and can never be
     * filled by gift-wrapping, producing the spurious FAIL(ii)/NO APEX. */
    size_t cap = in_R.N * 2 + 1;
    int dummy = 1;
    s_hash_table in_R_set;
    if (!hash_init(&in_R_set, sizeof(s_ncell *), sizeof(int), cap, in_R.N,
                   ncell_ptr_hash, ncell_ptr_eq, NULL)) {
        dynarray_free(&in_R); return 0;
    }
    for (unsigned i = 0; i < in_R.N; i++) {
        s_ncell *nc; dynarray_get_value(&in_R, i, &nc);
        hash_insert(&in_R_set, &nc, &dummy);
    }

    s_hash_table v1_set, v2_set;
    if (!hash_init(&v1_set, sizeof(int), sizeof(int),
                   in_R.N * 8 + 1, in_R.N * 4, int_hash, int_eq, NULL)) {
        hash_free(&in_R_set); dynarray_free(&in_R); return 0;
    }
    if (!hash_init(&v2_set, sizeof(int), sizeof(int),
                   in_R.N * 8 + 1, in_R.N * 4, int_hash, int_eq, NULL)) {
        hash_free(&v1_set); hash_free(&in_R_set); dynarray_free(&in_R); return 0;
    }

    /* Split R into half-cavities and collect per-side vertices in one pass. */
    for (unsigned i = 0; i < in_R.N; i++) {
        s_ncell *nc; dynarray_get_value(&in_R, i, &nc);
        bool above = false, below = false;
        for (int j = 0; j < 4; j++) {
            int o = cdt_orient3d(va_dt, vb_dt, vc_dt, nc->vertex_id[j]);
            if (o > 0) above = true;
            if (o < 0) below = true;
        }
        if (above) dynarray_push(&C1->tets, &nc);
        if (below) dynarray_push(&C2->tets, &nc);

        for (int j = 0; j < 4; j++) {
            int v = nc->vertex_id[j];
            int o = cdt_orient3d(va_dt, vb_dt, vc_dt, v);
            if (o >= 0 && !hash_get(&v1_set, &v)) {
                hash_insert(&v1_set, &v, &dummy);
                dynarray_push(&C1->verts, &v);
            }
            if (o <= 0 && !hash_get(&v2_set, &v)) {
                hash_insert(&v2_set, &v, &dummy);
                dynarray_push(&C2->verts, &v);
            }
        }
    }
    hash_free(&v1_set); hash_free(&v2_set);

    /* Build half-cavity boundaries.
     * A face is on the true R boundary iff its neighbor is outside in_R_set.
     * Non-plane faces: add to C1 if all verts orient >= 0, C2 if all <= 0.
     * Mixed-orientation faces (straddle the plane) are skipped.
     * Plane faces (all 3 verts orient == 0): adding to both halves would give each
     * toggle the same void_sign, so the half whose side doesn't match would search
     * for an apex in the wrong direction and fail. Instead route to exactly one half
     * by looking at the interior vertex (fids4[3]): above → C1, below → C2.
     * Already-constrained faces are skipped; the constraint face cap is added below. */
    for (unsigned i = 0; i < in_R.N; i++) {
        s_ncell *nc; dynarray_get_value(&in_R, i, &nc);
        for (int fi = 0; fi < 4; fi++) {
            s_ncell *nb = nc->opposite[fi];
            if (nb && hash_get(&in_R_set, &nb)) continue; /* internal to R */

            int fids[3]; int k = 0;
            for (int j = 0; j < 4; j++) if (j != fi) fids[k++] = nc->vertex_id[j];
            sort3(&fids[0], &fids[1], &fids[2]);
            if (is_constrained_face(constrained, fids[0], fids[1], fids[2])) continue;

            int o0 = cdt_orient3d(va_dt, vb_dt, vc_dt, fids[0]);
            int o1 = cdt_orient3d(va_dt, vb_dt, vc_dt, fids[1]);
            int o2 = cdt_orient3d(va_dt, vb_dt, vc_dt, fids[2]);
            int fids4[4] = {fids[0], fids[1], fids[2], nc->vertex_id[fi]};

            bool all_ge0 = (o0 >= 0 && o1 >= 0 && o2 >= 0);
            bool all_le0 = (o0 <= 0 && o1 <= 0 && o2 <= 0);

            if (all_ge0 && all_le0) {
                /* Plane face: route to one half only based on interior vertex side. */
                int oi = cdt_orient3d(va_dt, vb_dt, vc_dt, nc->vertex_id[fi]);
                if      (oi > 0) { if (!dynarray_push(&C1->boundary, fids4)) goto err; }
                else if (oi < 0) { if (!dynarray_push(&C2->boundary, fids4)) goto err; }
                /* oi == 0: all 4 verts coplanar — degenerate, skip. */
            } else {
                if (all_ge0) if (!dynarray_push(&C1->boundary, fids4)) goto err;
                if (all_le0) if (!dynarray_push(&C2->boundary, fids4)) goto err;
            }
        }
    }

    {
        int key[4] = {va_dt, vb_dt, vc_dt, -1};
        sort3(&key[0], &key[1], &key[2]);
        if (!dynarray_push(&C1->boundary, key)) goto err;
        if (!dynarray_push(&C2->boundary, key)) goto err;
    }

    hash_free(&in_R_set);
    dynarray_free(&in_R);
    return 1;

err:
    hash_free(&in_R_set);
    dynarray_free(&in_R);
    return 0;
}

/* -----------------------------------------------------------------------
 * Step 5 — Cavity Expansion
 * ----------------------------------------------------------------------- */

/* Adjacency info for a pre-removal boundary face → external global neighbor. */
typedef struct { s_ncell *nb; int opp_lid; } s_face_adj_entry;
/* Used during inner adjacency pairing of new tets. */
typedef struct { s_ncell *nc; int fi;     } s_face_nc_entry;
/* Toggle hash value: count (odd=active) + void_sign (precomputed orient3d side of the cavity). */
typedef struct { int count; int void_sign; } s_bnd_face_val;

/* Rebuild C->boundary by face-count of all C->tets (count==1 → boundary).
 * The constraint face is always appended at the end. */
static int rebuild_half_boundary(s_half_cavity *C,
                                  int va_dt, int vb_dt, int vc_dt)
{
    typedef struct { int count; int interior_v; } s_cnt_val;
    s_hash_table cnt;
    if (!hash_init(&cnt, sizeof(int[3]), sizeof(s_cnt_val),
                   C->tets.N * 8 + 1, C->tets.N * 4,
                   face_triple_hash, face_triple_eq, NULL)) return 0;

    for (unsigned i = 0; i < C->tets.N; i++) {
        s_ncell *nc; dynarray_get_value(&C->tets, i, &nc);
        for (int fi = 0; fi < 4; fi++) {
            int fids[3]; int k = 0;
            for (int j = 0; j < 4; j++) if (j != fi) fids[k++] = nc->vertex_id[j];
            sort3(&fids[0], &fids[1], &fids[2]);
            s_cnt_val *p = (s_cnt_val *)hash_get_or_create(&cnt, fids);
            if (!p) { hash_free(&cnt); return 0; }
            if (p->count == 0) p->interior_v = nc->vertex_id[fi];
            p->count++;
        }
    }

    dynarray_free(&C->boundary);
    C->boundary = dynarray_initialize(sizeof(int[4]), 32);
    if (!C->boundary.items) { hash_free(&cnt); return 0; }

    for (size_t b = 0; b < cnt.nbuckets; b++) {
        for (s_hash_entry *e = cnt.buckets[b]; e; e = e->next) {
            s_cnt_val *val = (s_cnt_val *)entry_value(&cnt, e);
            if (val->count == 1) {
                const int *key3 = (const int *)entry_key(e);
                int fids4[4] = {key3[0], key3[1], key3[2], val->interior_v};
                if (!dynarray_push(&C->boundary, fids4)) {
                    hash_free(&cnt); return 0;
                }
            }
        }
    }
    hash_free(&cnt);

    int key[4] = {va_dt, vb_dt, vc_dt, -1};
    sort3(&key[0], &key[1], &key[2]);
    return dynarray_push(&C->boundary, key);
}

/* Try to expand half-cavity C until every face in ∂C appears in the local DT
 * of Vi, or fail if a wrong-side vertex is encountered.
 * side: +1 = C1 (above), -1 = C2 (below).
 * On success writes *ldt_out (caller must free_complex), *l2g_out (caller must
 * free), *n_out = seed count (local DT real vertices = [4 .. 4+n-1]).
 * On failure returns 0; output pointers are untouched. */
static int try_cavity_expansion(s_scplx *dt,
                                 s_half_cavity *C,
                                 int va_dt, int vb_dt, int vc_dt,
                                 int side,
                                 double TOL,
                                 s_scplx *ldt_out,
                                 int **l2g_out,
                                 int *n_out)
{
    int n = (int)C->verts.N;
    int *l2g = malloc((size_t)n * sizeof(int));
    if (!l2g) return 0;

    s_hash_table g2l;
    if (!hash_init(&g2l, sizeof(int), sizeof(int), (size_t)(n*4+1), 0,
                   int_hash, int_eq, NULL)) {
        free(l2g); return 0;
    }
    for (int i = 0; i < n; i++) {
        int g; dynarray_get_value(&C->verts, i, &g);
        l2g[i] = g;
        hash_insert(&g2l, &g, &i);
    }

    /* Track the growing set of C-tets to identify external neighbors. */
    s_hash_table removal_set;
    if (!hash_init(&removal_set, sizeof(s_ncell *), sizeof(int),
                   C->tets.N * 4 + 1, 0, ncell_ptr_hash, ncell_ptr_eq, NULL)) {
        hash_free(&g2l); free(l2g); return 0;
    }
    int dummy = 1;
    for (unsigned i = 0; i < C->tets.N; i++) {
        s_ncell *nc; dynarray_get_value(&C->tets, i, &nc);
        hash_insert(&removal_set, &nc, &dummy);
    }

    int constraint_key[3] = {va_dt, vb_dt, vc_dt};
    sort3(&constraint_key[0], &constraint_key[1], &constraint_key[2]);

    s_scplx ldt = {0};
    int ret = 0;

    for (;;) {
        /* Build local DT of the current Vi. */
        if (ldt.head) { free_complex(&ldt); ldt = (s_scplx){0}; }

        s_point *lpts_arr = malloc((size_t)n * sizeof(s_point));
        if (!lpts_arr) goto done;
        for (int i = 0; i < n; i++) lpts_arr[i] = dt->points.p[l2g[i]];

        s_points lpts = { .N = n, .p = lpts_arr };
        s_dt_builder lb = dt_builder_begin(&lpts, NULL, TOL, NULL, NULL);
        free(lpts_arr);
        if (!lb._stack) goto done;
        ldt = dt_builder_end(&lb, true, NULL, NULL, 0);
        if (!ldt.head) goto done;

        /* Check every boundary face (except the constraint face). */
        int expand_face[3] = {-1, -1, -1};
        for (unsigned bi = 0; bi < C->boundary.N; bi++) {
            int *f = (int *)dynarray_get_ptr(&C->boundary, bi);
            if (f[0]==constraint_key[0] && f[1]==constraint_key[1] && f[2]==constraint_key[2])
                continue;
            int *pi0 = (int *)hash_get(&g2l, &f[0]);
            int *pi1 = (int *)hash_get(&g2l, &f[1]);
            int *pi2 = (int *)hash_get(&g2l, &f[2]);
            if (!pi0 || !pi1 || !pi2) goto done;
            if (face_in_dt(&ldt, *pi0+4, *pi1+4, *pi2+4)) continue;
            expand_face[0] = f[0]; expand_face[1] = f[1]; expand_face[2] = f[2];
            break;
        }

        if (expand_face[0] < 0) {
            *ldt_out = ldt; ldt = (s_scplx){0};
            *l2g_out = l2g; l2g = NULL;
            *n_out   = n;
            ret = 1;
            goto done;
        }

        /* Find which C-tet has expand_face with an external opposite. */
        s_ncell *expand_nb = NULL;
        for (unsigned ti = 0; ti < C->tets.N && !expand_nb; ti++) {
            s_ncell *nc; dynarray_get_value(&C->tets, ti, &nc);
            for (int fi = 0; fi < 4; fi++) {
                int fids[3]; int k = 0;
                for (int j = 0; j < 4; j++) if (j != fi) fids[k++] = nc->vertex_id[j];
                sort3(&fids[0], &fids[1], &fids[2]);
                if (fids[0]!=expand_face[0] || fids[1]!=expand_face[1] || fids[2]!=expand_face[2])
                    continue;
                s_ncell *nb = nc->opposite[fi];
                if (nb && !hash_get(&removal_set, &nb)) { expand_nb = nb; break; }
            }
        }
        if (!expand_nb) goto done;

        /* New vertex: 4th vertex of expand_nb not in expand_face. */
        int new_v = -1;
        for (int j = 0; j < 4; j++) {
            int v = expand_nb->vertex_id[j];
            if (v!=expand_face[0] && v!=expand_face[1] && v!=expand_face[2]) { new_v=v; break; }
        }
        if (new_v < 0) goto done;

        /* Side check: new vertex must not be on the wrong side of the plane. */
        int o = cdt_orient3d(va_dt, vb_dt, vc_dt, new_v);
        if (o != 0 && o != side) goto done;

        /* Grow Vi and the C-tet set. */
        if (!hash_get(&g2l, &new_v)) {
            int *nl2g = realloc(l2g, (size_t)(n+1) * sizeof(int));
            if (!nl2g) goto done;
            l2g = nl2g;
            l2g[n] = new_v;
            hash_insert(&g2l, &new_v, &n);
            n++;
            if (!dynarray_push(&C->verts, &new_v)) goto done;
        }
        if (!dynarray_push(&C->tets, &expand_nb)) goto done;
        hash_insert(&removal_set, &expand_nb, &dummy);

        if (!rebuild_half_boundary(C, va_dt, vb_dt, vc_dt)) goto done;
    }

done:
    if (ldt.head) free_complex(&ldt);
    hash_free(&g2l);
    hash_free(&removal_set);
    free(l2g);
    return ret;
}

/* Commit both half-cavities: remove old Tf tets, insert local-DT tets, fix adjacency.
 * On failure, the DT is inconsistent; the caller must abort. */
static int commit_cavity_expansion(s_scplx *dt,
                                    s_half_cavity *C1, s_scplx *ldt1, int *l2g1,
                                    s_half_cavity *C2, s_scplx *ldt2, int *l2g2)
{
    size_t total_tets = C1->tets.N + C2->tets.N;

    /* Build unified removal set. */
    s_hash_table removal_set;
    int dummy = 1;
    if (!hash_init(&removal_set, sizeof(s_ncell *), sizeof(int),
                   total_tets*2+1, total_tets, ncell_ptr_hash, ncell_ptr_eq, NULL)) return 0;
    for (int half = 0; half < 2; half++) {
        s_half_cavity *C = half ? C2 : C1;
        for (unsigned i = 0; i < C->tets.N; i++) {
            s_ncell *nc; dynarray_get_value(&C->tets, i, &nc);
            hash_insert(&removal_set, &nc, &dummy);
        }
    }

    /* Pre-collect external adjacency: for each C-tet face with non-C opposite,
     * store the neighbor and its opp_lid so we can reconnect after removal. */
    s_hash_table adj_map;
    if (!hash_init(&adj_map, sizeof(int[3]), sizeof(s_face_adj_entry),
                   (C1->boundary.N + C2->boundary.N)*2+1, 0,
                   face_triple_hash, face_triple_eq, NULL)) {
        hash_free(&removal_set); return 0;
    }
    for (int half = 0; half < 2; half++) {
        s_half_cavity *C = half ? C2 : C1;
        for (unsigned ti = 0; ti < C->tets.N; ti++) {
            s_ncell *nc; dynarray_get_value(&C->tets, ti, &nc);
            for (int fi = 0; fi < 4; fi++) {
                s_ncell *nb = nc->opposite[fi];
                if (!nb || hash_get(&removal_set, &nb)) continue;
                int fids[3]; int k = 0;
                for (int j = 0; j < 4; j++) if (j != fi) fids[k++] = nc->vertex_id[j];
                sort3(&fids[0], &fids[1], &fids[2]);
                if (hash_get(&adj_map, fids)) continue;
                int opp_lid = -1;
                for (int k2 = 0; k2 < 4; k2++)
                    if (nb->opposite[k2] == nc) { opp_lid = k2; break; }
                if (opp_lid < 0) continue;
                s_face_adj_entry fadj = {nb, opp_lid};
                hash_insert(&adj_map, fids, &fadj);
            }
        }
    }

    /* Fix external neighbors' opposite pointers and point2tet, then unlink/free. */
    for (int half = 0; half < 2; half++) {
        s_half_cavity *C = half ? C2 : C1;
        for (unsigned i = 0; i < C->tets.N; i++) {
            s_ncell *nc; dynarray_get_value(&C->tets, i, &nc);
            for (int fi = 0; fi < 4; fi++) {
                s_ncell *nb = nc->opposite[fi];
                if (!nb || hash_get(&removal_set, &nb)) continue;
                for (int k = 0; k < 4; k++)
                    if (nb->opposite[k] == nc) { nb->opposite[k] = NULL; break; }
            }
            for (int j = 0; j < 4; j++) {
                int v = nc->vertex_id[j];
                if (dt->point2tet[v] != nc) continue;
                s_ncell *rep = NULL;
                for (int fi = 0; fi < 4 && !rep; fi++) {
                    if (fi == j) continue;
                    s_ncell *nb = nc->opposite[fi];
                    if (nb && !hash_get(&removal_set, &nb)) rep = nb;
                }
                dt->point2tet[v] = rep;
            }
        }
    }
    /* Unlink and free using removal_set to avoid double-free of straddling tets
     * (tets with vertices on both sides appear in both C1->tets and C2->tets). */
    for (size_t rb = 0; rb < removal_set.nbuckets; rb++) {
        for (s_hash_entry *e = removal_set.buckets[rb]; e; e = e->next) {
            s_ncell *nc; memcpy(&nc, entry_key(e), sizeof(s_ncell *));
            if (nc->prev) nc->prev->next = nc->next; else dt->head = nc->next;
            if (nc->next) nc->next->prev = nc->prev;
            dt->N_ncells--;
            free_ncell(nc);
        }
    }
    hash_free(&removal_set);

    /* Insert real tets from local DTs into global DT. */
    s_dynarray new_tets = dynarray_initialize(sizeof(s_ncell *), total_tets + 8);
    if (!new_tets.items) { hash_free(&adj_map); return 0; }

    for (int half = 0; half < 2; half++) {
        s_scplx *ldt = half ? ldt2 : ldt1;
        int     *l2g = half ? l2g2 : l2g1;
        for (s_ncell *lnc = ldt->head; lnc; lnc = lnc->next) {
            bool sent = false;
            for (int j = 0; j < 4; j++) if (lnc->vertex_id[j] < 4) { sent = true; break; }
            if (sent) continue;

            s_ncell *nc = malloc_ncell();
            if (!nc) { dynarray_free(&new_tets); hash_free(&adj_map); return 0; }
            memset(nc, 0, sizeof(s_ncell));
            for (int j = 0; j < 4; j++) nc->vertex_id[j] = l2g[lnc->vertex_id[j] - 4];

            nc->next = dt->head; nc->prev = NULL;
            if (dt->head) dt->head->prev = nc;
            dt->head = nc;
            dt->N_ncells++;

            for (int j = 0; j < 4; j++)
                if (!dt->point2tet[nc->vertex_id[j]]) dt->point2tet[nc->vertex_id[j]] = nc;

            if (!dynarray_push(&new_tets, &nc)) {
                dynarray_free(&new_tets); hash_free(&adj_map); return 0;
            }
        }
    }

    /* Fix opposite pointers: inner pairing first, then external boundary. */
    s_hash_table face_nc_map;
    if (!hash_init(&face_nc_map, sizeof(int[3]), sizeof(s_face_nc_entry),
                   new_tets.N * 8 + 1, 0, face_triple_hash, face_triple_eq, NULL)) {
        dynarray_free(&new_tets); hash_free(&adj_map); return 0;
    }
    for (unsigned ti = 0; ti < new_tets.N; ti++) {
        s_ncell *nc; dynarray_get_value(&new_tets, ti, &nc);
        for (int fi = 0; fi < 4; fi++) {
            int fids[3]; int k = 0;
            for (int j = 0; j < 4; j++) if (j != fi) fids[k++] = nc->vertex_id[j];
            sort3(&fids[0], &fids[1], &fids[2]);
            s_face_nc_entry *ex = (s_face_nc_entry *)hash_get(&face_nc_map, fids);
            if (ex) {
                ex->nc->opposite[ex->fi] = nc;
                nc->opposite[fi]         = ex->nc;
            } else {
                s_face_nc_entry entry = {nc, fi};
                hash_insert(&face_nc_map, fids, &entry);
            }
        }
    }
    /* Connect remaining NULL opposite[fi] to pre-stored external neighbors. */
    for (unsigned ti = 0; ti < new_tets.N; ti++) {
        s_ncell *nc; dynarray_get_value(&new_tets, ti, &nc);
        for (int fi = 0; fi < 4; fi++) {
            if (nc->opposite[fi]) continue;
            int fids[3]; int k = 0;
            for (int j = 0; j < 4; j++) if (j != fi) fids[k++] = nc->vertex_id[j];
            sort3(&fids[0], &fids[1], &fids[2]);
            s_face_adj_entry *fadj = (s_face_adj_entry *)hash_get(&adj_map, fids);
            if (fadj) {
                nc->opposite[fi]                  = fadj->nb;
                fadj->nb->opposite[fadj->opp_lid] = nc;
            }
        }
    }

    hash_free(&face_nc_map);
    hash_free(&adj_map);
    dynarray_free(&new_tets);

    /* Final safety pass: any vertex whose representative was removed but that
     * still lives in a surviving tet not reached above (e.g. a sentinel ghost
     * tet pulled into the cavity by expansion) must get a valid point2tet. */
    for (s_ncell *nc = dt->head; nc; nc = nc->next)
        for (int j = 0; j < 4; j++)
            if (!dt->point2tet[nc->vertex_id[j]]) dt->point2tet[nc->vertex_id[j]] = nc;

    return 1;
}

/* -----------------------------------------------------------------------
 * Step 6 — Gift-Wrapping Fallback (§4.4)
 * ----------------------------------------------------------------------- */

/* Returns 1 if tet (t0,t1,t2,t3) and triangle (u0,u1,u2) intersect only at
 * a shared sub-simplex (no spurious crossing). */
/* Returns 1 = no spurious intersection (OK to use this apex), 0 = intersection found.
 * t0,t1,t2 = sigma face (already removed from boundary); t3 = candidate apex.
 * u0,u1,u2 = another active boundary face to test against.
 *
 * Mirrors the reference isTetIntersecting logic:
 * 1. Sigma half-space early exit: if all tri vertices are on the non-apex side of
 *    sigma (or on sigma), the triangle cannot intersect the apex-side geometry → OK.
 * 2. For n_sh==2: same early exit suffices (tri's unshared vertex on non-apex side).
 * 3. For n_sh==1: check unshared tri vertices inside tet; check only apex edges vs tri.
 * 4. For n_sh==0: check tri vertices in tet; check apex edges vs tri; check tri edges
 *    vs the three non-sigma tet faces only (skip sigma face). */
/* Exact segment-triangle intersection test covering all cases:
 * both endpoints coplanar, one endpoint coplanar, and non-coplanar.
 * Returns 1 if the interior of s1-s2 intersects the interior of <v0,v1,v2>
 * (including coplanar interior crossings). */
static int seg_tri_exact(int s1, int s2, int v0, int v1, int v2)
{
    int o1 = cdt_orient3d(v0, v1, v2, s1);
    int o2 = cdt_orient3d(v0, v1, v2, s2);

    if (o1 == 0 && o2 == 0) {
        /* Both coplanar: check if either endpoint is inside the triangle,
         * or if the segment crosses any triangle edge. */
        if (cdt_point_in_inner_tri(s1, v0, v1, v2)) return 1;
        if (cdt_point_in_inner_tri(s2, v0, v1, v2)) return 1;
        if (cdt_inner_segs_cross(s1, s2, v0, v1))   return 1;
        if (cdt_inner_segs_cross(s1, s2, v1, v2))   return 1;
        if (cdt_inner_segs_cross(s1, s2, v2, v0))   return 1;
        return 0;
    }
    if (o1 == 0) return cdt_point_in_inner_tri(s1, v0, v1, v2);
    if (o2 == 0) return cdt_point_in_inner_tri(s2, v0, v1, v2);
    if (o1 == o2) return 0;
    return cdt_inner_seg_crosses_inner_tri(s1, s2, v0, v1, v2);
}

/* Returns true if vertex v is occluded from tet <s0,s1,s2,w> by any triangle
 * in original_bnd (array of sorted int[3] face triples).  Mirrors the reference
 * isTetLocallyDelaunay Test 2: for each boundary triangle, check whether it
 * blocks the line of sight between v and the tet by testing:
 *   (a) 4 segments v→{s0,s1,s2,w} vs the triangle interior
 *   (b) 18 boundary-triangle edge vs cone-face {sigma-edge,v}/{apex-edge,v} checks
 *   (c) 9 boundary-triangle vertex inside non-sigma sub-tet {face,v} checks
 * Returns false (not occluded / visible) if no boundary triangle blocks the view. */
static bool is_occluded(int v, int s0, int s1, int s2, int w,
                         const s_dynarray *original_bnd)
{
    for (unsigned i = 0; i < original_bnd->N; i++) {
        const int *u = (const int *)dynarray_get_ptr((s_dynarray *)original_bnd, i);
        int u0 = u[0], u1 = u[1], u2 = u[2];

        /* Skip if tri is a face of the tet (all 3 verts are tet vertices). */
        if ((u0==s0||u0==s1||u0==s2||u0==w) &&
            (u1==s0||u1==s1||u1==s2||u1==w) &&
            (u2==s0||u2==s1||u2==s2||u2==w)) continue;

        /* Skip if v is coplanar with this boundary triangle (mirrors reference
         * bnd_o3dc==0 check: coplanar vertex is not occluded by the triangle). */
        int hv = cdt_orient3d(u0, u1, u2, v);
        if (hv == 0) continue;

        /* Half-space early exit: if {s0,s1,s2,w,v} are all on the same strict side
         * of the boundary triangle's plane, no intersection is possible. */
        int hs0 = cdt_orient3d(u0, u1, u2, s0);
        int hs1 = cdt_orient3d(u0, u1, u2, s1);
        int hs2 = cdt_orient3d(u0, u1, u2, s2);
        int hw  = cdt_orient3d(u0, u1, u2, w);
        if (hv==hs0 && hv==hs1 && hv==hs2 && hv==hw) continue;

        /* (a) Segments v→{s0,s1,s2,w} vs inner triangle (u0,u1,u2). */
        if (seg_tri_exact(v, s0, u0, u1, u2)) return true;
        if (seg_tri_exact(v, s1, u0, u1, u2)) return true;
        if (seg_tri_exact(v, s2, u0, u1, u2)) return true;
        if (seg_tri_exact(v, w,  u0, u1, u2)) return true;

        /* (b) Tri edges vs "cone" faces: sigma edges paired with v, and apex
         *     edges paired with v.  These are the faces of the cone from v to
         *     the tet, excluding sigma itself. */
        if (seg_tri_exact(u0,u1, s0,s1,v)) return true;
        if (seg_tri_exact(u0,u1, s1,s2,v)) return true;
        if (seg_tri_exact(u0,u1, s2,s0,v)) return true;
        if (seg_tri_exact(u0,u1, s0,w, v)) return true;
        if (seg_tri_exact(u0,u1, s1,w, v)) return true;
        if (seg_tri_exact(u0,u1, s2,w, v)) return true;
        if (seg_tri_exact(u1,u2, s0,s1,v)) return true;
        if (seg_tri_exact(u1,u2, s1,s2,v)) return true;
        if (seg_tri_exact(u1,u2, s2,s0,v)) return true;
        if (seg_tri_exact(u1,u2, s0,w, v)) return true;
        if (seg_tri_exact(u1,u2, s1,w, v)) return true;
        if (seg_tri_exact(u1,u2, s2,w, v)) return true;
        if (seg_tri_exact(u2,u0, s0,s1,v)) return true;
        if (seg_tri_exact(u2,u0, s1,s2,v)) return true;
        if (seg_tri_exact(u2,u0, s2,s0,v)) return true;
        if (seg_tri_exact(u2,u0, s0,w, v)) return true;
        if (seg_tri_exact(u2,u0, s1,w, v)) return true;
        if (seg_tri_exact(u2,u0, s2,w, v)) return true;

        /* (c) Tri vertices inside non-sigma sub-tets {face,v} (mirrors
         *     intersectionTEST_3 for the three non-sigma tet faces). */
        if (cdt_point_in_tet(u0, s0,s1,w,v)) return true;
        if (cdt_point_in_tet(u1, s0,s1,w,v)) return true;
        if (cdt_point_in_tet(u2, s0,s1,w,v)) return true;
        if (cdt_point_in_tet(u0, s1,s2,w,v)) return true;
        if (cdt_point_in_tet(u1, s1,s2,w,v)) return true;
        if (cdt_point_in_tet(u2, s1,s2,w,v)) return true;
        if (cdt_point_in_tet(u0, s2,s0,w,v)) return true;
        if (cdt_point_in_tet(u1, s2,s0,w,v)) return true;
        if (cdt_point_in_tet(u2, s2,s0,w,v)) return true;
    }
    return false;
}

/* Returns 1 if the tetrahedron <t0,t1,t2,t3> does NOT spuriously intersect
 * the boundary triangle <u0,u1,u2>. The sigma face is <t0,t1,t2>.
 * Mirrors the reference isTetIntersecting: 5 half-space early exits first,
 * then vertex-in-tet, apex edges vs triangle, triangle edges vs non-sigma faces. */
static int no_spurious_intersect(int t0, int t1, int t2, int t3,
                                  int u0, int u1, int u2)
{
    int uv[3] = {u0, u1, u2};

    /* Apex sign: orient3D(sigma, apex). Positive = apex on cavity side. */
    int apex_sign = cdt_orient3d(t0, t1, t2, t3);

    /* Half-space early exit 1 (sigma face, reference check 1, uses <=):
     * All u vertices on the non-apex side of sigma → no intersection possible.
     * Shared u vertices have orient=0 which satisfies the condition. */
    {
        int s0 = cdt_orient3d(t0, t1, t2, u0);
        int s1 = cdt_orient3d(t0, t1, t2, u1);
        int s2 = cdt_orient3d(t0, t1, t2, u2);
        if (s0 * apex_sign <= 0 && s1 * apex_sign <= 0 && s2 * apex_sign <= 0)
            return 1;
    }

    /* Half-space early exits 2-4 (non-sigma faces, reference checks 2-4, use <):
     * For each face the opposite tet vertex lies on the apex_sign side;
     * all u on -apex_sign side means they're strictly outside T through that face.
     * orient3D(t2,t1,t3,t0) = apex_sign  (face opp t0)
     * orient3D(t0,t2,t3,t1) = apex_sign  (face opp t1)
     * orient3D(t1,t0,t3,t2) = apex_sign  (face opp t2) */
    {
        int e0 = cdt_orient3d(t2,t1,t3,u0), e1 = cdt_orient3d(t2,t1,t3,u1), e2 = cdt_orient3d(t2,t1,t3,u2);
        if (e0 * apex_sign < 0 && e1 * apex_sign < 0 && e2 * apex_sign < 0) return 1;
    }
    {
        int f0 = cdt_orient3d(t0,t2,t3,u0), f1 = cdt_orient3d(t0,t2,t3,u1), f2 = cdt_orient3d(t0,t2,t3,u2);
        if (f0 * apex_sign < 0 && f1 * apex_sign < 0 && f2 * apex_sign < 0) return 1;
    }
    {
        int g0 = cdt_orient3d(t1,t0,t3,u0), g1 = cdt_orient3d(t1,t0,t3,u1), g2 = cdt_orient3d(t1,t0,t3,u2);
        if (g0 * apex_sign < 0 && g1 * apex_sign < 0 && g2 * apex_sign < 0) return 1;
    }

    /* Half-space early exit 5 (triangle plane, reference check at lines 415-416):
     * All T vertices on the same side (or plane) of u → no intersection. */
    {
        int h0 = cdt_orient3d(u0,u1,u2,t0), h1 = cdt_orient3d(u0,u1,u2,t1);
        int h2 = cdt_orient3d(u0,u1,u2,t2), h3 = cdt_orient3d(u0,u1,u2,t3);
        if ((h0 <= 0 && h1 <= 0 && h2 <= 0 && h3 <= 0) ||
            (h0 >= 0 && h1 >= 0 && h2 >= 0 && h3 >= 0)) return 1;
    }

    /* Triangle vertices strictly inside tet (skip vertices shared with tet). */
    for (int i = 0; i < 3; i++) {
        if (uv[i] == t0 || uv[i] == t1 || uv[i] == t2 || uv[i] == t3) continue;
        if (cdt_point_in_tet(uv[i], t0, t1, t2, t3)) return 0;
    }

    /* Apex edges (t0-t3, t1-t3, t2-t3) vs triangle. */
    if (seg_tri_exact(t0, t3, u0, u1, u2)) return 0;
    if (seg_tri_exact(t1, t3, u0, u1, u2)) return 0;
    if (seg_tri_exact(t2, t3, u0, u1, u2)) return 0;

    /* Triangle edges vs 3 non-sigma tet faces {t1,t2,t3}, {t0,t2,t3}, {t0,t1,t3}. */
    for (int ei = 0; ei < 3; ei++) {
        int a = uv[ei], b = uv[(ei+1)%3];
        if (seg_tri_exact(a, b, t1, t2, t3)) return 0;
        if (seg_tri_exact(a, b, t0, t2, t3)) return 0;
        if (seg_tri_exact(a, b, t0, t1, t3)) return 0;
    }

    return 1;
}

/* a,b,c must already be sorted. void_sign = precomputed orient3d side of the cavity.
 * Only stored on first insertion (count==0 → 1). */
static int toggle_boundary_face(s_hash_table *toggle, int a, int b, int c, int void_sign)
{
    int f[3] = {a, b, c}; /* caller must pre-sort */
    s_bnd_face_val *p = (s_bnd_face_val *)hash_get_or_create(toggle, f);
    if (!p) return 0;
    if (p->count == 0) p->void_sign = void_sign;
    p->count++;
    return 1;
}

/* Gift-wrap one half-cavity.
 * toggle: face → int count (odd = active).
 * usable_v: dynarray of int DT-indices; negative means removed.
 * original_bnd: snapshot of the initial boundary (int[3] sorted triples) used
 *   for the occlusion test in condition (iii). */
static int gw_half(s_scplx *dt,
                    s_hash_table *toggle,
                    s_dynarray *usable_v,
                    int va_dt, int vb_dt, int vc_dt,
                    s_dynarray *new_tets_out,
                    const s_dynarray *original_bnd)
{
    int ck[3] = {va_dt, vb_dt, vc_dt};
    sort3(&ck[0], &ck[1], &ck[2]);

    for (;;) {
        /* Collect active non-constraint boundary faces (sorted int[3] keys). */
        s_dynarray active = dynarray_initialize(sizeof(int[3]), 16);
        if (!active.items) return 0;

        for (size_t b = 0; b < toggle->nbuckets; b++) {
            for (s_hash_entry *e = toggle->buckets[b]; e; e = e->next) {
                const s_bnd_face_val *val = (const s_bnd_face_val *)entry_value(toggle, e);
                if (val->count % 2 == 0) continue;
                const int *k = (const int *)entry_key(e);
                if (k[0]==ck[0] && k[1]==ck[1] && k[2]==ck[2]) continue;
                if (!dynarray_push(&active, k)) { dynarray_free(&active); return 0; }
            }
        }
        if (active.N == 0) { dynarray_free(&active); return 1; }

        bool made_progress = false;
        for (unsigned si = 0; si < active.N && !made_progress; si++) {
            const int *sigma = (const int *)dynarray_get_ptr(&active, si);

            /* int_sign: the cavity (void) side sign for this face.
             * Stored as void_sign in the toggle at seeding/toggle time. */
            const s_bnd_face_val *sv = (const s_bnd_face_val *)hash_get(toggle, sigma);
            if (!sv || sv->void_sign == 0) continue;
            int int_sign = sv->void_sign;

            int apex = -1;
            for (unsigned vi = 0; vi < usable_v->N; vi++) {
                int w; dynarray_get_value(usable_v, vi, &w);
                if (w < 0) continue;
                /* (i) Correct side. */
                int o_i = cdt_orient3d(sigma[0], sigma[1], sigma[2], w);
                if (o_i != int_sign) continue;

                /* (ii) No spurious intersections with active faces or the
                 * constraint face.  The constraint face is excluded from
                 * `active` (we never try to connect it), but a candidate tet
                 * that pokes through the constraint plane must still be
                 * rejected, so we test against ck explicitly. */
                bool ok = true;
                for (unsigned ai = 0; ai < active.N && ok; ai++) {
                    const int *u = (const int *)dynarray_get_ptr(&active, ai);
                    if (!no_spurious_intersect(sigma[0], sigma[1], sigma[2], w,
                                               u[0], u[1], u[2]))
                        ok = false;
                }
                if (ok && !no_spurious_intersect(sigma[0], sigma[1], sigma[2], w,
                                                  ck[0], ck[1], ck[2]))
                    ok = false;
                if (!ok) continue;

                /* (iii) No VISIBLE vertex inside circumsphere (SoS).
                 * v is visible if (a) it is on the same side of sigma as w,
                 * and (b) no original boundary triangle occludes the line of
                 * sight from v to the tet <sigma, w>. */
                for (unsigned ui = 0; ui < usable_v->N && ok; ui++) {
                    int v; dynarray_get_value(usable_v, ui, &v);
                    if (v < 0 || v==sigma[0] || v==sigma[1] || v==sigma[2] || v==w)
                        continue;
                    int ov = cdt_orient3d(sigma[0], sigma[1], sigma[2], v);
                    if (ov != 0 && ov != int_sign) continue; /* wrong side → not visible */
                    if (is_occluded(v, sigma[0], sigma[1], sigma[2], w, original_bnd))
                        continue; /* blocked by cavity wall → not visible */
                    if (cdt_perturbed_insphere(sigma[0], sigma[1], sigma[2], w, v) == 1)
                        ok = false;
                }
                if (!ok) continue;

                apex = w;
                break;
            }
            if (apex < 0) continue;

            /* Create tet and link into dt. */
            s_ncell *nc = malloc_ncell();
            if (!nc) { dynarray_free(&active); return 0; }
            memset(nc, 0, sizeof(s_ncell));
            nc->vertex_id[0] = sigma[0]; nc->vertex_id[1] = sigma[1];
            nc->vertex_id[2] = sigma[2]; nc->vertex_id[3] = apex;
            nc->next = dt->head; nc->prev = NULL;
            if (dt->head) dt->head->prev = nc;
            dt->head = nc;
            dt->N_ncells++;
            for (int j = 0; j < 4; j++)
                if (!dt->point2tet[nc->vertex_id[j]])
                    dt->point2tet[nc->vertex_id[j]] = nc;
            if (!dynarray_push(new_tets_out, &nc)) { dynarray_free(&active); return 0; }

            /* Toggle boundary: remove sigma (void_sign irrelevant for inactive face).
             * For each new face, void_sign = -orient3d(sorted_new_face, sigma_opposite_vertex).
             * The sigma vertex NOT in the new face is inside the new tet (filled side),
             * so negating its orient3d gives the cavity (void) side. */
            int f01[3] = {sigma[0], sigma[1], apex}; sort3(&f01[0],&f01[1],&f01[2]);
            int f12[3] = {sigma[1], sigma[2], apex}; sort3(&f12[0],&f12[1],&f12[2]);
            int f20[3] = {sigma[2], sigma[0], apex}; sort3(&f20[0],&f20[1],&f20[2]);
            int vs01 = -cdt_orient3d(f01[0],f01[1],f01[2], sigma[2]);
            int vs12 = -cdt_orient3d(f12[0],f12[1],f12[2], sigma[0]);
            int vs20 = -cdt_orient3d(f20[0],f20[1],f20[2], sigma[1]);
            if (!toggle_boundary_face(toggle, sigma[0], sigma[1], sigma[2], 0) ||
                !toggle_boundary_face(toggle, f01[0],   f01[1],   f01[2],   vs01) ||
                !toggle_boundary_face(toggle, f12[0],   f12[1],   f12[2],   vs12) ||
                !toggle_boundary_face(toggle, f20[0],   f20[1],   f20[2],   vs20)) {
                dynarray_free(&active); return 0;
            }
            made_progress = true;
        }
        dynarray_free(&active);
        if (!made_progress) return 0;
    }
}

/* Gift-wrapping fallback for when cavity expansion fails.
 * C1/C2 still hold original Tf tets; they are removed here.
 * On success, Tf tets are replaced and the constraint face is filled. */
static int gift_wrap_face(s_scplx *dt,
                           s_half_cavity *C1, s_half_cavity *C2,
                           int va_dt, int vb_dt, int vc_dt)
{
    size_t total = C1->tets.N + C2->tets.N;
    int dummy = 1;

    /* Build unified removal set. */
    s_hash_table removal_set;
    if (!hash_init(&removal_set, sizeof(s_ncell *), sizeof(int),
                   total*2+1, total, ncell_ptr_hash, ncell_ptr_eq, NULL)) return 0;
    for (int half = 0; half < 2; half++) {
        s_half_cavity *C = half ? C2 : C1;
        for (unsigned i = 0; i < C->tets.N; i++) {
            s_ncell *nc; dynarray_get_value(&C->tets, i, &nc);
            hash_insert(&removal_set, &nc, &dummy);
        }
    }

    /* Pre-collect external adjacency. */
    s_hash_table adj_map;
    if (!hash_init(&adj_map, sizeof(int[3]), sizeof(s_face_adj_entry),
                   (C1->boundary.N + C2->boundary.N)*2+1, 0,
                   face_triple_hash, face_triple_eq, NULL)) {
        hash_free(&removal_set); return 0;
    }
    for (int half = 0; half < 2; half++) {
        s_half_cavity *C = half ? C2 : C1;
        for (unsigned ti = 0; ti < C->tets.N; ti++) {
            s_ncell *nc; dynarray_get_value(&C->tets, ti, &nc);
            for (int fi = 0; fi < 4; fi++) {
                s_ncell *nb = nc->opposite[fi];
                if (!nb || hash_get(&removal_set, &nb)) continue;
                int fids[3]; int k = 0;
                for (int j = 0; j < 4; j++) if (j != fi) fids[k++] = nc->vertex_id[j];
                sort3(&fids[0], &fids[1], &fids[2]);
                if (hash_get(&adj_map, fids)) continue;
                int opp_lid = -1;
                for (int k2 = 0; k2 < 4; k2++)
                    if (nb->opposite[k2] == nc) { opp_lid = k2; break; }
                if (opp_lid < 0) continue;
                s_face_adj_entry fadj = {nb, opp_lid};
                hash_insert(&adj_map, fids, &fadj);
            }
        }
    }

    /* Fix external neighbors, point2tet, then unlink+free Tf tets. */
    for (int half = 0; half < 2; half++) {
        s_half_cavity *C = half ? C2 : C1;
        for (unsigned i = 0; i < C->tets.N; i++) {
            s_ncell *nc; dynarray_get_value(&C->tets, i, &nc);
            for (int fi = 0; fi < 4; fi++) {
                s_ncell *nb = nc->opposite[fi];
                if (!nb || hash_get(&removal_set, &nb)) continue;
                for (int k = 0; k < 4; k++)
                    if (nb->opposite[k] == nc) { nb->opposite[k] = NULL; break; }
            }
            for (int j = 0; j < 4; j++) {
                int v = nc->vertex_id[j];
                if (dt->point2tet[v] != nc) continue;
                s_ncell *rep = NULL;
                for (int fi = 0; fi < 4 && !rep; fi++) {
                    if (fi == j) continue;
                    s_ncell *nb = nc->opposite[fi];
                    if (nb && !hash_get(&removal_set, &nb)) rep = nb;
                }
                dt->point2tet[v] = rep;
            }
        }
    }
    /* Unlink and free via removal_set to avoid double-free of straddling tets. */
    for (size_t rb = 0; rb < removal_set.nbuckets; rb++) {
        for (s_hash_entry *e = removal_set.buckets[rb]; e; e = e->next) {
            s_ncell *nc; memcpy(&nc, entry_key(e), sizeof(s_ncell *));
            if (nc->prev) nc->prev->next = nc->next; else dt->head = nc->next;
            if (nc->next) nc->next->prev = nc->prev;
            dt->N_ncells--;
            free_ncell(nc);
        }
    }
    hash_free(&removal_set);

    s_dynarray new_tets = dynarray_initialize(sizeof(s_ncell *), (int)(total + 8));
    if (!new_tets.items) { hash_free(&adj_map); return 0; }

    /* Gift-wrap C1 (above). */
    s_hash_table toggle1;
    if (!hash_init(&toggle1, sizeof(int[3]), sizeof(s_bnd_face_val),
                   C1->boundary.N * 4 + 1, 0,
                   face_triple_hash, face_triple_eq, NULL)) {
        dynarray_free(&new_tets); hash_free(&adj_map); return 0;
    }
    for (unsigned i = 0; i < C1->boundary.N; i++) {
        const int *f = (const int *)dynarray_get_ptr(&C1->boundary, i);
        /* f[0..2] already sorted; f[3] = interior vertex on the void (cavity) side.
         * void_sign = orient3d(face, f[3]) — f[3] is always on the cavity side. */
        int vs = 0;
        if (f[3] >= 0)
            vs = cdt_orient3d(f[0],f[1],f[2],f[3]);
        if (!toggle_boundary_face(&toggle1, f[0], f[1], f[2], vs)) {
            hash_free(&toggle1); dynarray_free(&new_tets); hash_free(&adj_map); return 0;
        }
    }
    /* Snapshot initial C1 boundary as int[3] triples for the occlusion test. */
    s_dynarray orig_bnd1 = dynarray_initialize(sizeof(int[3]), (int)C1->boundary.N + 1);
    if (!orig_bnd1.items) {
        hash_free(&toggle1); dynarray_free(&new_tets); hash_free(&adj_map); return 0;
    }
    for (unsigned i = 0; i < C1->boundary.N; i++) {
        const int *f = (const int *)dynarray_get_ptr(&C1->boundary, i);
        int key[3] = {f[0], f[1], f[2]};
        dynarray_push(&orig_bnd1, key);
    }

    s_dynarray usable1 = dynarray_initialize(sizeof(int), (int)C1->verts.N + 1);
    if (!usable1.items) {
        dynarray_free(&orig_bnd1);
        hash_free(&toggle1); dynarray_free(&new_tets); hash_free(&adj_map); return 0;
    }
    for (unsigned i = 0; i < C1->verts.N; i++) {
        int v; dynarray_get_value(&C1->verts, i, &v);
        dynarray_push(&usable1, &v);
    }
    int ok = gw_half(dt, &toggle1, &usable1, va_dt, vb_dt, vc_dt, &new_tets, &orig_bnd1);
    hash_free(&toggle1); dynarray_free(&usable1); dynarray_free(&orig_bnd1);
    if (!ok) { dynarray_free(&new_tets); hash_free(&adj_map); return 0; }

    /* Gift-wrap C2 (below). */
    s_hash_table toggle2;
    if (!hash_init(&toggle2, sizeof(int[3]), sizeof(s_bnd_face_val),
                   C2->boundary.N * 4 + 1, 0,
                   face_triple_hash, face_triple_eq, NULL)) {
        dynarray_free(&new_tets); hash_free(&adj_map); return 0;
    }
    for (unsigned i = 0; i < C2->boundary.N; i++) {
        const int *f = (const int *)dynarray_get_ptr(&C2->boundary, i);
        /* void_sign = orient3d(face, f[3]) — f[3] is always on the cavity side. */
        int vs = 0;
        if (f[3] >= 0)
            vs = cdt_orient3d(f[0],f[1],f[2],f[3]);
        if (!toggle_boundary_face(&toggle2, f[0], f[1], f[2], vs)) {
            hash_free(&toggle2); dynarray_free(&new_tets); hash_free(&adj_map); return 0;
        }
    }
    /* Snapshot initial C2 boundary as int[3] triples for the occlusion test. */
    s_dynarray orig_bnd2 = dynarray_initialize(sizeof(int[3]), (int)C2->boundary.N + 1);
    if (!orig_bnd2.items) {
        hash_free(&toggle2); dynarray_free(&new_tets); hash_free(&adj_map); return 0;
    }
    for (unsigned i = 0; i < C2->boundary.N; i++) {
        const int *f = (const int *)dynarray_get_ptr(&C2->boundary, i);
        int key[3] = {f[0], f[1], f[2]};
        dynarray_push(&orig_bnd2, key);
    }

    s_dynarray usable2 = dynarray_initialize(sizeof(int), (int)C2->verts.N + 1);
    if (!usable2.items) {
        dynarray_free(&orig_bnd2);
        hash_free(&toggle2); dynarray_free(&new_tets); hash_free(&adj_map); return 0;
    }
    for (unsigned i = 0; i < C2->verts.N; i++) {
        int v; dynarray_get_value(&C2->verts, i, &v);
        dynarray_push(&usable2, &v);
    }
    ok = gw_half(dt, &toggle2, &usable2, va_dt, vb_dt, vc_dt, &new_tets, &orig_bnd2);
    hash_free(&toggle2); dynarray_free(&usable2); dynarray_free(&orig_bnd2);
    if (!ok) { dynarray_free(&new_tets); hash_free(&adj_map); return 0; }

    /* Fix opposite pointers: inner pairing then external. */
    s_hash_table face_nc_map;
    if (!hash_init(&face_nc_map, sizeof(int[3]), sizeof(s_face_nc_entry),
                   new_tets.N * 8 + 1, 0, face_triple_hash, face_triple_eq, NULL)) {
        dynarray_free(&new_tets); hash_free(&adj_map); return 0;
    }
    for (unsigned ti = 0; ti < new_tets.N; ti++) {
        s_ncell *nc; dynarray_get_value(&new_tets, ti, &nc);
        for (int fi = 0; fi < 4; fi++) {
            int fids[3]; int k = 0;
            for (int j = 0; j < 4; j++) if (j != fi) fids[k++] = nc->vertex_id[j];
            sort3(&fids[0], &fids[1], &fids[2]);
            s_face_nc_entry *ex = (s_face_nc_entry *)hash_get(&face_nc_map, fids);
            if (ex) {
                ex->nc->opposite[ex->fi] = nc;
                nc->opposite[fi]         = ex->nc;
            } else {
                s_face_nc_entry entry = {nc, fi};
                hash_insert(&face_nc_map, fids, &entry);
            }
        }
    }
    for (unsigned ti = 0; ti < new_tets.N; ti++) {
        s_ncell *nc; dynarray_get_value(&new_tets, ti, &nc);
        for (int fi = 0; fi < 4; fi++) {
            if (nc->opposite[fi]) continue;
            int fids[3]; int k = 0;
            for (int j = 0; j < 4; j++) if (j != fi) fids[k++] = nc->vertex_id[j];
            sort3(&fids[0], &fids[1], &fids[2]);
            s_face_adj_entry *fadj = (s_face_adj_entry *)hash_get(&adj_map, fids);
            if (fadj) {
                nc->opposite[fi]                  = fadj->nb;
                fadj->nb->opposite[fadj->opp_lid] = nc;
            }
        }
    }
    hash_free(&face_nc_map);
    hash_free(&adj_map);
    dynarray_free(&new_tets);

    /* Final safety pass: restore point2tet for any vertex orphaned during
     * removal but still living in a surviving tet (see commit_cavity_expansion). */
    for (s_ncell *nc = dt->head; nc; nc = nc->next)
        for (int j = 0; j < 4; j++)
            if (!dt->point2tet[nc->vertex_id[j]]) dt->point2tet[nc->vertex_id[j]] = nc;
    return 1;
}

/* -----------------------------------------------------------------------
 * Step 7 — Ghost-Vertex Flood-Fill for Interior/Exterior (§4.5)
 * ----------------------------------------------------------------------- */

static void classify_interior_exterior(s_scplx *dt, s_hash_table *constrained)
{
    int ext_stamp = ++dt->mark_stamp;
    int int_stamp = ++dt->mark_stamp;

    s_dynarray queue = dynarray_initialize(sizeof(s_ncell *), 64);
    if (!queue.items) return;

    /* Phase 1: label all sentinel-incident tets as external. */
    for (s_ncell *nc = dt->head; nc; nc = nc->next) {
        for (int j = 0; j < 4; j++) {
            if (nc->vertex_id[j] < 4) { nc->mark_token = ext_stamp; break; }
        }
    }

    /* Phase 2: seed queue from real neighbors of ghost tets. */
    for (s_ncell *nc = dt->head; nc; nc = nc->next) {
        if (nc->mark_token != ext_stamp) continue;
        for (int fi = 0; fi < 4; fi++) {
            s_ncell *nb = nc->opposite[fi];
            if (!nb || nb->mark_token == ext_stamp || nb->mark_token == int_stamp) continue;
            int fids[3]; int k = 0;
            for (int j = 0; j < 4; j++) if (j != fi) fids[k++] = nc->vertex_id[j];
            sort3(&fids[0], &fids[1], &fids[2]);
            nb->mark_token = is_constrained_face(constrained, fids[0], fids[1], fids[2])
                             ? int_stamp : ext_stamp;
            dynarray_push(&queue, &nb);
        }
    }

    /* Phase 3: BFS — propagate label, flip at constrained faces. */
    while (queue.N > 0) {
        s_ncell *cur; dynarray_pop(&queue, &cur);
        int cur_label = cur->mark_token;
        for (int fi = 0; fi < 4; fi++) {
            s_ncell *nb = cur->opposite[fi];
            if (!nb || nb->mark_token == ext_stamp || nb->mark_token == int_stamp) continue;
            int fids[3]; int k = 0;
            for (int j = 0; j < 4; j++) if (j != fi) fids[k++] = cur->vertex_id[j];
            sort3(&fids[0], &fids[1], &fids[2]);
            nb->mark_token = is_constrained_face(constrained, fids[0], fids[1], fids[2])
                             ? (cur_label == ext_stamp ? int_stamp : ext_stamp)
                             : cur_label;
            dynarray_push(&queue, &nb);
        }
    }
    dynarray_free(&queue);

    /* Phase 4: delete all external tets. */
    s_ncell *nc = dt->head;
    while (nc) {
        s_ncell *next = nc->next;
        if (nc->mark_token == ext_stamp) {
            /* Null out opposite pointers in internal neighbors. */
            for (int fi = 0; fi < 4; fi++) {
                s_ncell *nb = nc->opposite[fi];
                if (!nb || nb->mark_token == ext_stamp) continue;
                for (int k = 0; k < 4; k++)
                    if (nb->opposite[k] == nc) { nb->opposite[k] = NULL; break; }
            }
            /* Fix point2tet. */
            for (int j = 0; j < 4; j++) {
                int v = nc->vertex_id[j];
                if (v < 4) { dt->point2tet[v] = NULL; continue; }
                if (dt->point2tet[v] != nc) continue;
                s_ncell *rep = NULL;
                for (int fi = 0; fi < 4 && !rep; fi++) {
                    s_ncell *nb = nc->opposite[fi];
                    if (nb && nb->mark_token == int_stamp) rep = nb;
                }
                dt->point2tet[v] = rep;
            }
            if (nc->prev) nc->prev->next = nc->next; else dt->head = nc->next;
            if (nc->next) nc->next->prev = nc->prev;
            dt->N_ncells--;
            free_ncell(nc);
        }
        nc = next;
    }
}

/* -----------------------------------------------------------------------
 * Public entry point
 * ----------------------------------------------------------------------- */

s_scplx tetrahedralize_interior_trimesh(const s_trimesh *mesh,
                                         double EPS_DEG, double TOL)
{
    if (!trimesh_is_valid(mesh)) return (s_scplx){0};

    s_trimesh working = copy_trimesh(mesh);
    if (!trimesh_is_valid(&working)) return (s_scplx){0};

    /* Phase A: DT + segment recovery */
    s_dt_builder b = dt_builder_begin(&working.points, NULL, TOL, NULL, NULL);
    if (!b._stack) { free_trimesh(&working); return (s_scplx){0}; }

    if (!ensure_ridge_protected(&b, &working, EPS_DEG, TOL)) {
        s_scplx tmp = dt_builder_end(&b, false, NULL, NULL, 0);
        free_complex(&tmp); free_trimesh(&working); return (s_scplx){0};
    }

    /* Keep the big tetra: the CDT relies on sentinels staying at indices 0-3
     * (real vertices at working_vertex+4) for the +4 indexing convention and
     * for Phase C's ghost-tet detection (vertex_id < 4).  Building with
     * keep_big_tetra=false would compact point indices to 0-based and shrink
     * point2tet, making every +4 access out of bounds. */
    s_scplx dt = dt_builder_end(&b, true, NULL, NULL, 0);
    if (!dt.head) { free_trimesh(&working); return (s_scplx){0}; }

    /* Sync working.points with the compacted DT point set. */
    {
        s_point *p = malloc((size_t)dt.points.N * sizeof(s_point));
        if (!p) { free_complex(&dt); free_trimesh(&working); return (s_scplx){0}; }
        memcpy(p, dt.points.p, (size_t)dt.points.N * sizeof(s_point));
        free(working.points.p);
        working.points.p = p;
        working.points.N = dt.points.N;
    }

    /* Phase B: face recovery (Steps 4-6) — Algorithm 2 outer loop.
     * Gift-wrapping guarantees every face is recovered in one pass; the
     * while terminates immediately on the second iteration. */
    s_hash_table constrained;
    if (!hash_init(&constrained, sizeof(int[3]), sizeof(int),
                   (size_t)(working.Nf * 4 + 1), (size_t)(working.Nf * 2),
                   face_triple_hash, face_triple_eq, NULL)) {
        free_complex(&dt); free_trimesh(&working); return (s_scplx){0};
    }

#define PHASE_B_ABORT() do { \
        half_cavity_free(&C1); half_cavity_free(&C2); \
        hash_free(&constrained); free_trimesh(&working); \
        free_complex(&dt); cdt_predicates_clear(); return (s_scplx){0}; \
    } while (0)

    bool any_missing = true;
    while (any_missing) {
        any_missing = false;

        for (int fi = 0; fi < working.Nf; fi++) {
            /* working.faces is 0-based; DT vertex indices are +4 (sentinels at 0-3). */
            int va = working.faces[fi*3]   + 4;
            int vb = working.faces[fi*3+1] + 4;
            int vc = working.faces[fi*3+2] + 4;

            if (face_in_dt(&dt, va, vb, vc)) {
                mark_constrained_face(&constrained, va, vb, vc);
                continue;
            }
            any_missing = true;

            s_half_cavity C1 = {0}, C2 = {0};
            if (!build_half_cavities(&dt, va, vb, vc, &constrained, &C1, &C2))
                PHASE_B_ABORT();

            s_scplx ldt1 = {0}, ldt2 = {0};
            int *l2g1 = NULL, *l2g2 = NULL, n1 = 0, n2 = 0;

            int ok1 = try_cavity_expansion(&dt, &C1, va, vb, vc, +1, TOL, &ldt1, &l2g1, &n1);
            int ok2 = ok1
                      ? try_cavity_expansion(&dt, &C2, va, vb, vc, -1, TOL, &ldt2, &l2g2, &n2)
                      : 0;

            if (ok1 && ok2) {
                int commit_ok = commit_cavity_expansion(&dt,
                                    &C1, &ldt1, l2g1,
                                    &C2, &ldt2, l2g2);
                if (ldt1.head) free_complex(&ldt1); free(l2g1);
                if (ldt2.head) free_complex(&ldt2); free(l2g2);
                if (commit_ok)
                    mark_constrained_face(&constrained, va, vb, vc);
                else
                    PHASE_B_ABORT();
            } else {
                if (ldt1.head) free_complex(&ldt1); free(l2g1);
                if (ldt2.head) free_complex(&ldt2); free(l2g2);
                if (gift_wrap_face(&dt, &C1, &C2, va, vb, vc))
                    mark_constrained_face(&constrained, va, vb, vc);
                else
                    PHASE_B_ABORT();
            }

            half_cavity_free(&C1);
            half_cavity_free(&C2);
        }
    }

#undef PHASE_B_ABORT

    /* Phase C: interior/exterior classification (Step 7) */
    classify_interior_exterior(&dt, &constrained);

    hash_free(&constrained);
    free_trimesh(&working);
    cdt_predicates_clear();

    return dt;
}
