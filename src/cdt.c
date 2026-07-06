/* WARNING: This file was written by claude code */

#include "delaunay.h"
#include "scplx.h"
#include "trimesh.h"
#include "dynarray.h"
#include "hash.h"
#include "gtests.h"
#include "points.h"
#include "cdt_predicates.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>

/* Internal Phase B recovery diagnostics.  Off by default (these fire on benign
 * events like cavity-expansion giving up on a face before gift-wrap recovers
 * it); set CDT_VERBOSE=1 to see them.  The user-facing centroid-split fallback
 * warning is a separate, always-on fprintf. */
static int g_cdt_verbose = -1;
#define CDT_DBG(...) do { \
        if (g_cdt_verbose < 0) g_cdt_verbose = getenv("CDT_VERBOSE") ? 1 : 0; \
        if (g_cdt_verbose) fprintf(stderr, "[cdt] " __VA_ARGS__); \
    } while (0)
static int g_dbg_this_face = 0;

/* Structural validation (CDT_VALIDATE=1): adjacency reciprocity and shared
 * vertices.  Returns the number of violations found (0 = consistent). */
static int validate_complex(const s_scplx *dt)
{
    int bad = 0, nulls = 0;
    for (s_ncell *nc = dt->head; nc; nc = nc->next) {
        for (int fi = 0; fi < 4; fi++) {
            s_ncell *nb = nc->opposite[fi];
            if (!nb) { nulls++; continue; }
            int back = -1;
            for (int k = 0; k < 4; k++) if (nb->opposite[k] == nc) { back = k; break; }
            if (back < 0) { bad++; continue; }
            for (int j = 0; j < 4; j++) {
                if (j == fi) continue;
                bool found = false;
                for (int k = 0; k < 4; k++)
                    if (nb->vertex_id[k] == nc->vertex_id[j]) { found = true; break; }
                if (!found) { bad++; break; }
            }
        }
    }
    /* With the big tetra kept, exactly 4 NULL opposites (its hull) are
     * expected; any extra NULL is a hole left by an incomplete refill. */
    if (nulls > 4) bad += nulls - 4;
    return bad;
}


/* LNC: LiNear Combination, a point implicitly defined as a lin. comb. of 2 explicit points */
typedef struct {
    int    v1, v2;  /* DT indices of explicit points, v1==-1 means explicit vertex. */
    double t;
} s_lnc_info;

/* Half-cavity from splitting the region Tf by the constraint plane.
 * All indices are DT indices. */
typedef struct {
    s_dynarray tets;      /* s_ncell*, tets with >=1 vertex on this side */
    s_dynarray verts;     /* int, DT indices of vertices on/above (C1) or on/below (C2) */
    s_dynarray boundary;  /* int[4], {sorted face triple, interior_v} of bound(Ci) */
} s_half_cavity;

static void half_cavity_free(s_half_cavity *c)
{
    dynarray_free(&c->tets);
    dynarray_free(&c->verts);
    dynarray_free(&c->boundary);
}

/* Memory buffers used throughout */
typedef struct {
    s_dynarray ring;      /* s_ncell* -- face_in_dt, find_region_R_boundary seed, edge_in_dt */
    s_dynarray out_ids;   /* int      -- edge_in_dt neighbor ids */
    s_dynarray stack;     /* s_ncell* -- find_region_R_boundary BFS */
    s_dynarray in_R;      /* s_ncell* -- build_half_cavities */
    s_dynarray queue;     /* s_ncell* -- classify_interior_exterior */
    s_dynarray new_tets;  /* s_ncell* -- commit_cavity_expansion, gift_wrap_face */
    s_dynarray orig_bnd;  /* int[3]   -- snapshot of initial cavity boundary for is_occluded */
    s_dynarray usable;    /* int      -- usable1/2 in gift_wrap_face */
    s_dynarray active;    /* int[3]   -- active faces per gw_half iteration */
    s_dynarray lnc_info;  /* s_lnc_info -- ensure_ridge_protected */
    s_dynarray to_split;  /* int[2]   -- ensure_ridge_protected inner loop */
    s_dynarray comm_cur;  /* int[3]   -- gw_half current deferred on-plane faces */
} s_cdt_scratch;

static void cdt_scratch_free(s_cdt_scratch *sc)
{
    dynarray_free(&sc->ring);
    dynarray_free(&sc->out_ids);
    dynarray_free(&sc->stack);
    dynarray_free(&sc->in_R);
    dynarray_free(&sc->queue);
    dynarray_free(&sc->new_tets);
    dynarray_free(&sc->orig_bnd);
    dynarray_free(&sc->usable);
    dynarray_free(&sc->active);
    dynarray_free(&sc->lnc_info);
    dynarray_free(&sc->to_split);
    dynarray_free(&sc->comm_cur);
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

static inline int ncell_local_id(const s_ncell *nc, int v)
{
    for (int k = 0; k < 4; k++) if (nc->vertex_id[k] == v) return k;
    return -1;
}


static bool edge_in_dt(s_scplx *dt, int va, int vb, s_cdt_scratch *sc)
{
    if (!dt->point2tet[va]) return false;

    dynarray_clear(&sc->out_ids);
    dynarray_clear(&sc->ring);

    int va_lid = ncell_local_id(dt->point2tet[va], va);
    vertex_neighbors(dt, va, dt->point2tet[va], va_lid, 0, &sc->out_ids, &sc->ring);

    for (unsigned i = 0; i < sc->out_ids.N; i++) {
        int nb; dynarray_get_value(&sc->out_ids, i, &nb);
        if (nb == vb) return true;
    }
    return false;
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
 * refpt_id: result of scout_refpt (-1 -> midpoint).
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
            double lab = sqrt((B.x-A.x)*(B.x-A.x) +
                              (B.y-A.y)*(B.y-A.y) +
                              (B.z-A.z)*(B.z-A.z));
            if (pa == va_dt || pb == va_dt) {
                double lar = sqrt((R.x-A.x)*(R.x-A.x) +
                                  (R.y-A.y)*(R.y-A.y) +
                                  (R.z-A.z)*(R.z-A.z));
                t = (lab > EPS_DEG) ? lar / lab : 0.5;
                adj = 1;
            } else if (pa == vb_dt || pb == vb_dt) {
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
static int ensure_ridge_protected(s_dt_builder *b, s_trimesh *mesh, double EPS_DEG, double TOL,
                                   s_cdt_scratch *sc)
{
    /* Registry: in exact mode dt_builder_begin_exact already cleared and
     * registered every initial point; in coordinate mode we must populate it
     * here (Phase B consumes it regardless of how the DT was built). */
    if (!b->dt.exact_ids) {
        cdt_predicates_clear();
        for (int i = 0; i < b->dt.points.N; i++) {
            s_point *p = &b->dt.points.p[i];
            cdt_point_set_explicit(i, p->x, p->y, p->z);
        }
    }

    /* lnc_info[dt_point_id]: v1==-1 for explicit vertices, DT indices for Steiners. */
    dynarray_clear(&sc->lnc_info);
    s_lnc_info explicit_entry = { .v1 = -1, .v2 = 0, .t = 0.0 };
    for (int i = 0; i < b->dt.points.N; i++) {
        if (!dynarray_push(&sc->lnc_info, &explicit_entry)) return 0;
    }

    for (;;) {
        s_hash_table edge_set;
        int dummy = 1;
        if (!hash_init(&edge_set, sizeof(int[2]), sizeof(int),
                       (size_t)(mesh->Nf * 4 + 1), (size_t)(mesh->Nf),
                       edge_pair_hash, edge_pair_eq, NULL)) return 0;

        dynarray_clear(&sc->to_split);

        for (int fi = 0; fi < mesh->Nf; fi++) {
            for (int ei = 0; ei < 3; ei++) {
                int va = mesh->faces[fi*3 + (ei+1)%3];
                int vb = mesh->faces[fi*3 + (ei+2)%3];
                if (edge_in_dt(&b->dt, va + 4, vb + 4, sc)) continue;

                int pair[2] = { va < vb ? va : vb, va < vb ? vb : va };
                if (hash_get(&edge_set, pair)) continue;
                if (!hash_insert(&edge_set, pair, &dummy)) {
                    hash_free(&edge_set); return 0;
                }
                if (!dynarray_push(&sc->to_split, pair)) {
                    hash_free(&edge_set); return 0;
                }
            }
        }

        hash_free(&edge_set);

        if (sc->to_split.N == 0) break;

        for (unsigned i = 0; i < sc->to_split.N; i++) {
            int *pair = (int *)dynarray_get_ptr(&sc->to_split, i);
            int va = pair[0], vb = pair[1];
            int va_dt = va + 4, vb_dt = vb + 4;

            int refpt_id = scout_refpt(&b->dt, va_dt, vb_dt, EPS_DEG);
            double t = compute_steiner_t(&b->dt, va_dt, vb_dt, refpt_id,
                                          (const s_lnc_info *)sc->lnc_info.items,
                                          (int)sc->lnc_info.N, EPS_DEG);

            int new_id = b->dt.points.N;
            /* Steiner at fraction t from A: M = (1-t)*A + t*B.  In the
             * cdt_point_set_lnc convention s*V[va] + (1-s)*V[vb] this is s = 1-t,
             * so both paths store M in points.p AND register the matching exact
             * LNC (consistent -- no reflection across the midpoint). */
            if (b->dt.exact_ids) {
                if (!dt_builder_extend_lnc(b, va_dt, vb_dt, 1.0 - t, TOL)) return 0;
            } else {
                s_point A = b->dt.points.p[va_dt];
                s_point B = b->dt.points.p[vb_dt];
                s_point M = { .x = (1.0-t)*A.x + t*B.x,
                              .y = (1.0-t)*A.y + t*B.y,
                              .z = (1.0-t)*A.z + t*B.z };
                s_points single = { .N = 1, .p = &M };
                if (!dt_builder_extend(b, &single, TOL)) return 0;
                if (b->dt.points.N > new_id)
                    cdt_point_set_lnc(new_id, va_dt, vb_dt, 1.0 - t);
            }

            if (b->dt.points.N > (int)sc->lnc_info.N) {
                s_lnc_info info = { .v1 = va_dt, .v2 = vb_dt, .t = t };
                if (!dynarray_push(&sc->lnc_info, &info)) return 0;
            }

            if (!split_trimesh_edge(mesh, va, vb, new_id - 4)) return 0;
        }
    }
    return 1;
}


/* ----------------------------------------------------------------------- */


static s_ncell *face_in_dt(const s_scplx *dt, int va, int vb, int vc, s_cdt_scratch *sc)
{
    s_ncell *start = dt->point2tet[va];
    if (!start) return NULL;
    int va_lid = ncell_local_id(start, va);
    if (va_lid < 0) return NULL;
    int comp[3]; for (int i = 0, k = 0; i < 4; i++) if (i != va_lid) comp[k++] = i;

    dynarray_clear(&sc->ring);
    ncells_incident_face((s_scplx *)dt, start, 0, comp, &sc->ring);

    for (unsigned i = 0; i < sc->ring.N; i++) {
        s_ncell *nc; dynarray_get_value(&sc->ring, i, &nc);
        bool has_vb = false, has_vc = false;
        for (int j = 0; j < 4; j++) {
            if (nc->vertex_id[j] == vb) has_vb = true;
            if (nc->vertex_id[j] == vc) has_vc = true;
        }
        if (has_vb && has_vc) return nc;
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

/* Collect the 3 vertex ids of face fi of nc into fids[3]. */
#define FACE_IDS(nc, fi, fids) do { \
    int _k = 0; \
    for (int _j = 0; _j < 4; _j++) if (_j != (fi)) (fids)[_k++] = (nc)->vertex_id[_j]; \
} while (0)

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



/* ----------------------------------------------------------------------- */

/* --- find_region_R (Step 4) --------------------------------------------- */

/* Exact test: does tet nc "block" constraint triangle (va,vb,vc), i.e. does
 * it improperly intersect the closed triangle so that f cannot be a face of
 * the DT while nc exists?  Three exact sub-tests (paper's Tf, extended with
 * the coplanar contacts our per-triangle recovery requires):
 *   1. a tet edge strictly straddling the plane pierces f's interior;
 *   2. a coplanar tet edge crosses f's interior (2D);
 *   3. a coplanar tet vertex lies strictly inside f (2D).
 * Any other improper contact would imply a DT edge/face crossing one of f's
 * protected edges, which cannot happen in a valid complex. */
static bool tet_blocks_face(int va, int vb, int vc, const s_ncell *nc)
{
    int ids[4], o[4];
    for (int j = 0; j < 4; j++) {
        ids[j] = nc->vertex_id[j];
        o[j]   = cdt_orient3d(va, vb, vc, ids[j]);
    }

    for (int j = 0; j < 4; j++) {
        for (int k = j + 1; k < 4; k++) {
            int u = ids[j], v = ids[k];
            bool u_corner = (u==va || u==vb || u==vc);
            bool v_corner = (v==va || v==vb || v==vc);
            if (u_corner && v_corner) continue; /* an edge/chord of f itself */

            if (o[j] * o[k] < 0) {
                /* Straddling edge piercing f's interior. */
                if (cdt_inner_seg_crosses_inner_tri(u, v, va, vb, vc))
                    return true;
            } else if (o[j] == 0 && o[k] == 0) {
                /* Coplanar edge vs f's interior (2D). */
                if (!u_corner && cdt_point_in_inner_tri(u, va, vb, vc)) return true;
                if (!v_corner && cdt_point_in_inner_tri(v, va, vb, vc)) return true;
                if (cdt_inner_segs_cross(u, v, va, vb)) return true;
                if (cdt_inner_segs_cross(u, v, vb, vc)) return true;
                if (cdt_inner_segs_cross(u, v, vc, va)) return true;
            }
        }
    }

    /* Coplanar vertex strictly inside f (vertex-only contact). */
    for (int j = 0; j < 4; j++) {
        if (o[j] != 0) continue;
        if (ids[j]==va || ids[j]==vb || ids[j]==vc) continue;
        if (cdt_point_in_inner_tri(ids[j], va, vb, vc)) return true;
    }
    return false;
}

/* Collect the tets improperly intersecting constraint triangle (va,vb,vc).
 * Seeds: blockers among the tets incident to f's corners.  Growth: neighbors
 * that block.  The blocked part of f's region is tiled by blockers, so the
 * set is face-connected across f's region and reachable from the corners.
 * Returns 1 on success, 0 on allocation failure. */
static int find_region_R(s_scplx *dt,
                          int va, int vb, int vc,
                          s_dynarray *in_R,
                          s_cdt_scratch *sc)
{
    dt->mark_stamp++;
    int stamp = dt->mark_stamp;

    /* sc->stack doubles as the result vector (iterated by index). */
    dynarray_clear(&sc->stack);

    int corners[3] = {va, vb, vc};
    for (int e = 0; e < 3; e++) {
        int p = corners[e];
        if (!dt->point2tet[p]) continue;
        s_ncell *start = dt->point2tet[p];
        int p_lid = ncell_local_id(start, p);
        if (p_lid < 0) continue;
        int comp[3]; for (int i = 0, k = 0; i < 4; i++) if (i != p_lid) comp[k++] = i;

        dynarray_clear(&sc->ring);
        ncells_incident_face(dt, start, 0, comp, &sc->ring);
        for (unsigned i = 0; i < sc->ring.N; i++) {
            s_ncell *nc; dynarray_get_value(&sc->ring, i, &nc);
            if (nc->mark_token == stamp) continue;
            if (!tet_blocks_face(va, vb, vc, nc)) continue;
            nc->mark_token = stamp;
            if (!dynarray_push(&sc->stack, &nc)) return 0;
        }
    }

    /* Fallback for a blocked region not touching any corner ring. */
    if (sc->stack.N == 0) {
        for (s_ncell *nc = dt->head; nc; nc = nc->next) {
            if (tet_blocks_face(va, vb, vc, nc)) {
                nc->mark_token = stamp;
                if (!dynarray_push(&sc->stack, &nc)) return 0;
                break;
            }
        }
    }

    for (unsigned k = 0; k < sc->stack.N; k++) {
        s_ncell *cur; dynarray_get_value(&sc->stack, k, &cur);
        for (int fi = 0; fi < 4; fi++) {
            s_ncell *nb = cur->opposite[fi];
            if (!nb || nb->mark_token == stamp) continue;
            if (!tet_blocks_face(va, vb, vc, nb)) continue;
            nb->mark_token = stamp;
            if (!dynarray_push(&sc->stack, &nb)) return 0;
        }
    }

    for (unsigned k = 0; k < sc->stack.N; k++) {
        s_ncell *nc; dynarray_get_value(&sc->stack, k, &nc);
        if (!dynarray_push(in_R, &nc)) return 0;
    }
    return 1;
}


/* Split Tf by the constraint plane, producing two half-cavities C1 (above)
 * and C2 (below).  va_dt/vb_dt/vc_dt are DT indices of the constraint face.
 * C1->tets: Tf tets with >=1 vertex strictly above.
 * C2->tets: Tf tets with >=1 vertex strictly below.
 * C1->verts: DT indices of vertices on or above the plane.
 * C2->verts: DT indices of vertices on or below the plane.
 * C1->boundary / C2->boundary: sorted int[3] face triples of dCi, including
 * the constraint face itself. */
static int build_half_cavities(s_scplx *dt,
                                int va_dt, int vb_dt, int vc_dt,
                                s_hash_table *constrained,
                                s_half_cavity *C1, s_half_cavity *C2,
                                s_cdt_scratch *sc)
{
    (void)constrained; /* region growth is bounded by the blocker test alone */

    dynarray_clear(&sc->in_R);

    if (!find_region_R(dt, va_dt, vb_dt, vc_dt, &sc->in_R, sc))
        return 0;

    C1->tets     = dynarray_initialize(sizeof(s_ncell *), 32);
    C1->verts    = dynarray_initialize(sizeof(int), 32);
    C1->boundary = dynarray_initialize(sizeof(int[4]), 32);
    C2->tets     = dynarray_initialize(sizeof(s_ncell *), 32);
    C2->verts    = dynarray_initialize(sizeof(int), 32);
    C2->boundary = dynarray_initialize(sizeof(int[4]), 32);
    if (!C1->tets.items || !C1->verts.items || !C1->boundary.items ||
        !C2->tets.items || !C2->verts.items || !C2->boundary.items) {
        return 0;
    }

    size_t cap = sc->in_R.N * 2 + 1;
    int dummy = 1;
    s_hash_table in_R_set;
    if (!hash_init(&in_R_set, sizeof(s_ncell *), sizeof(int), cap, sc->in_R.N,
                   ncell_ptr_hash, ncell_ptr_eq, NULL)) return 0;
    for (unsigned i = 0; i < sc->in_R.N; i++) {
        s_ncell *nc; dynarray_get_value(&sc->in_R, i, &nc);
        hash_insert(&in_R_set, &nc, &dummy);
    }

    s_hash_table v1_set, v2_set;
    if (!hash_init(&v1_set, sizeof(int), sizeof(int),
                   sc->in_R.N * 8 + 1, sc->in_R.N * 4, int_hash, int_eq, NULL)) {
        hash_free(&in_R_set); return 0;
    }
    if (!hash_init(&v2_set, sizeof(int), sizeof(int),
                   sc->in_R.N * 8 + 1, sc->in_R.N * 4, int_hash, int_eq, NULL)) {
        hash_free(&v1_set); hash_free(&in_R_set); return 0;
    }

    /* Split R into half-cavities and collect per-side vertices in one pass. */
    for (unsigned i = 0; i < sc->in_R.N; i++) {
        s_ncell *nc; dynarray_get_value(&sc->in_R, i, &nc);
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

    for (unsigned i = 0; i < sc->in_R.N; i++) {
        s_ncell *nc; dynarray_get_value(&sc->in_R, i, &nc);
        for (int fi = 0; fi < 4; fi++) {
            s_ncell *nb = nc->opposite[fi];
            if (nb && hash_get(&in_R_set, &nb)) continue; /* internal to R */

            int fids[3]; FACE_IDS(nc, fi, fids);
            sort3(&fids[0], &fids[1], &fids[2]);
            /* NOTE: constrained faces are NOT skipped here.  A previously
             * recovered face adjacent to the cavity is still a wall of the
             * cavity; dropping it leaves the boundary open and gift-wrapping
             * can never close it (this caused the lung-lobe Phase B aborts). */

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
                /* oi == 0: all 4 verts coplanar -- degenerate, skip. */
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
    return 1;

err:
    hash_free(&in_R_set);
    return 0;
}

/* -----------------------------------------------------------------------
 * Step 5 -- Cavity Expansion
 * ----------------------------------------------------------------------- */

/* Adjacency info for a pre-removal boundary face -> external global neighbor. */
typedef struct { s_ncell *nb; int opp_lid; } s_face_adj_entry;
/* Used during inner adjacency pairing of new tets. */
typedef struct { s_ncell *nc; int fi;     } s_face_nc_entry;
/* Toggle hash value: count (odd=active) + void_sign (precomputed orient3d side of the cavity). */
typedef struct { int count; int void_sign; } s_bnd_face_val;

/* Rebuild C->boundary by face-count of all C->tets (count==1 -> boundary).
 * Only faces on this half's side of the constraint plane are kept (side=+1
 * keeps all-o>=0 faces, side=-1 all-o<=0; plane faces are routed by the
 * interior vertex), mirroring build_half_cavities: the wrong-side faces of
 * straddling tets belong to the OTHER half-cavity.
 * The constraint face is always appended at the end. */
static int rebuild_half_boundary(s_half_cavity *C,
                                  int va_dt, int vb_dt, int vc_dt, int side)
{
    typedef struct { int count; int interior_v; } s_cnt_val;
    s_hash_table cnt;
    if (!hash_init(&cnt, sizeof(int[3]), sizeof(s_cnt_val),
                   C->tets.N * 8 + 1, C->tets.N * 4,
                   face_triple_hash, face_triple_eq, NULL)) return 0;

    for (unsigned i = 0; i < C->tets.N; i++) {
        s_ncell *nc; dynarray_get_value(&C->tets, i, &nc);
        for (int fi = 0; fi < 4; fi++) {
            int fids[3]; FACE_IDS(nc, fi, fids);
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
                int o0 = cdt_orient3d(va_dt, vb_dt, vc_dt, key3[0]);
                int o1 = cdt_orient3d(va_dt, vb_dt, vc_dt, key3[1]);
                int o2 = cdt_orient3d(va_dt, vb_dt, vc_dt, key3[2]);
                bool keep;
                if (o0 == 0 && o1 == 0 && o2 == 0) {
                    /* Plane face: belongs to the half its interior vertex is on. */
                    keep = (cdt_orient3d(va_dt, vb_dt, vc_dt, val->interior_v) == side);
                } else if (side > 0) {
                    keep = (o0 >= 0 && o1 >= 0 && o2 >= 0);
                } else {
                    keep = (o0 <= 0 && o1 <= 0 && o2 <= 0);
                }
                if (!keep) continue;
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

/* Try to expand half-cavity C until every face in dC appears in the local DT
 * of Vi, or fail if a wrong-side vertex is encountered.
 * side: +1 = C1 (above), -1 = C2 (below).
 * On success writes *ldt_out (caller must free_complex), *l2g_out (caller must
 * free), *n_out = vertex count of the expanded local DT (real vertices = [4 .. 4+n-1]).
 * On failure returns 0; output pointers are untouched. */
static int try_cavity_expansion(s_scplx *dt,
                                 s_half_cavity *C,
                                 int va_dt, int vb_dt, int vc_dt,
                                 int side,
                                 double TOL,
                                 s_hash_table *constrained,
                                 s_scplx *ldt_out,
                                 int **l2g_out,
                                 int *n_out,
                                 s_cdt_scratch *sc)
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
    int dbg_iter = 0;

    for (;;) {
        /* Build local DT of the current Vi. */
        if (ldt.head) { free_complex(&ldt); ldt = (s_scplx){0}; }

        s_point *lpts_arr = malloc((size_t)n * sizeof(s_point));
        if (!lpts_arr) goto done;
        for (int i = 0; i < n; i++) lpts_arr[i] = dt->points.p[l2g[i]];

        s_points lpts = { .N = n, .p = lpts_arr };
        /* Match the global DT's geometry: an exact global DT gets an exact local
         * cavity DT sharing the same registry (l2g maps local real ids -> global;
         * sentinels use scratch slots past the global point count). */
        s_dt_builder lb = dt->exact_ids
            ? dt_builder_begin_exact_local(&lpts, TOL, l2g, n, dt->points.N)
            : dt_builder_begin(&lpts, NULL, TOL, NULL, NULL);
        free(lpts_arr);
        if (!lb._stack) { CDT_DBG("expand side=%d iter=%d: lb begin fail\n", side, dbg_iter); goto done; }
        ldt = dt_builder_end(&lb, true, NULL, NULL, 0);
        if (!ldt.head) { CDT_DBG("expand side=%d iter=%d: ldt build fail (n=%d)\n", side, dbg_iter, n); goto done; }

        /* Check every boundary face.  The constraint face itself MUST also be
         * present in the local DT: committing a local DT that lacks it would
         * rebuild the cavity without recovering f, and the Phase B outer loop
         * would retry the exact same cavity forever. */
        int expand_face[3] = {-1, -1, -1};
        bool missing_constrained = false;
        for (unsigned bi = 0; bi < C->boundary.N; bi++) {
            int *f = (int *)dynarray_get_ptr(&C->boundary, bi);
            bool is_ck = (f[0]==constraint_key[0] && f[1]==constraint_key[1] &&
                          f[2]==constraint_key[2]);
            int *pi0 = (int *)hash_get(&g2l, &f[0]);
            int *pi1 = (int *)hash_get(&g2l, &f[1]);
            int *pi2 = (int *)hash_get(&g2l, &f[2]);
            if (!pi0 || !pi1 || !pi2) {
                CDT_DBG("expand side=%d iter=%d: bnd face (%d,%d,%d) has vertex not in Vi\n",
                        side, dbg_iter, f[0], f[1], f[2]);
                goto done;
            }
            if (face_in_dt(&ldt, *pi0+4, *pi1+4, *pi2+4, sc)) continue;
            /* Neither the constraint face nor a previously recovered
             * constrained wall face may be expanded across; try the other
             * missing faces first and fail (-> gift-wrapping) if they never
             * conform. */
            if (is_ck || is_constrained_face(constrained, f[0], f[1], f[2])) {
                missing_constrained = true;
                continue;
            }
            expand_face[0] = f[0]; expand_face[1] = f[1]; expand_face[2] = f[2];
            break;
        }

        if (expand_face[0] < 0) {
            if (missing_constrained) {
                CDT_DBG("expand side=%d iter=%d: constrained wall face not in local DT\n",
                        side, dbg_iter);
                goto done;
            }
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
                int fids[3]; FACE_IDS(nc, fi, fids);
                sort3(&fids[0], &fids[1], &fids[2]);
                if (fids[0]!=expand_face[0] || fids[1]!=expand_face[1] || fids[2]!=expand_face[2])
                    continue;
                s_ncell *nb = nc->opposite[fi];
                if (nb && !hash_get(&removal_set, &nb)) { expand_nb = nb; break; }
            }
        }
        if (!expand_nb) {
            CDT_DBG("expand side=%d iter=%d: no external neighbor for face (%d,%d,%d)\n",
                    side, dbg_iter, expand_face[0], expand_face[1], expand_face[2]);
            goto done;
        }

        /* New vertex: 4th vertex of expand_nb not in expand_face. */
        int new_v = -1;
        for (int j = 0; j < 4; j++) {
            int v = expand_nb->vertex_id[j];
            if (v!=expand_face[0] && v!=expand_face[1] && v!=expand_face[2]) { new_v=v; break; }
        }
        if (new_v < 0) goto done;

        /* Side check: new vertex must not be on the wrong side of the plane. */
        int o = cdt_orient3d(va_dt, vb_dt, vc_dt, new_v);
        if (o != 0 && o != side) {
            CDT_DBG("expand side=%d iter=%d: wrong-side vertex %d (o=%d)\n",
                    side, dbg_iter, new_v, o);
            goto done;
        }
        if (++dbg_iter > 2000) {
            CDT_DBG("expand side=%d: iteration cap reached, giving up\n", side);
            goto done;
        }

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

        if (!rebuild_half_boundary(C, va_dt, vb_dt, vc_dt, side)) goto done;
    }

done:
    if (ldt.head) free_complex(&ldt);
    hash_free(&g2l);
    hash_free(&removal_set);
    free(l2g);
    return ret;
}

/* Remove C1+C2 tets from dt, recording external adjacency in *adj_map (caller
 * must hash_free it on success).  Uses removal_set to avoid double-free of tets
 * that straddle both halves.  Returns 1 on success, 0 on allocation failure. */
static int remove_cavity_tets(s_scplx *dt,
                               s_half_cavity *C1, s_half_cavity *C2,
                               s_hash_table *adj_map)
{
    size_t total = C1->tets.N + C2->tets.N;
    int dummy = 1;

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

    if (!hash_init(adj_map, sizeof(int[3]), sizeof(s_face_adj_entry),
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
                int fids[3]; FACE_IDS(nc, fi, fids);
                sort3(&fids[0], &fids[1], &fids[2]);
                if (hash_get(adj_map, fids)) continue;
                int opp_lid = -1;
                for (int k2 = 0; k2 < 4; k2++)
                    if (nb->opposite[k2] == nc) { opp_lid = k2; break; }
                if (opp_lid < 0) continue;
                s_face_adj_entry fadj = {nb, opp_lid};
                hash_insert(adj_map, fids, &fadj);
            }
        }
    }

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
    /* Use removal_set to avoid double-free of tets straddling both halves. */
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
    return 1;
}

/* Wire opposite pointers for new_tets: inner pairing first, then external via
 * adj_map.  Finishes with a safety pass that restores any still-null point2tet
 * entries (covers ghost tets pulled in by expansion).
 * Returns 1 on success, 0 on allocation failure. */
static int link_new_tets(s_scplx *dt, s_dynarray *new_tets, s_hash_table *adj_map)
{
    s_hash_table face_nc_map;
    if (!hash_init(&face_nc_map, sizeof(int[3]), sizeof(s_face_nc_entry),
                   new_tets->N * 8 + 1, 0, face_triple_hash, face_triple_eq, NULL))
        return 0;

    for (unsigned ti = 0; ti < new_tets->N; ti++) {
        s_ncell *nc; dynarray_get_value(new_tets, ti, &nc);
        for (int fi = 0; fi < 4; fi++) {
            int fids[3]; FACE_IDS(nc, fi, fids);
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
    for (unsigned ti = 0; ti < new_tets->N; ti++) {
        s_ncell *nc; dynarray_get_value(new_tets, ti, &nc);
        for (int fi = 0; fi < 4; fi++) {
            if (nc->opposite[fi]) continue;
            int fids[3]; FACE_IDS(nc, fi, fids);
            sort3(&fids[0], &fids[1], &fids[2]);
            s_face_adj_entry *fadj = (s_face_adj_entry *)hash_get(adj_map, fids);
            if (fadj) {
                nc->opposite[fi]                  = fadj->nb;
                fadj->nb->opposite[fadj->opp_lid] = nc;
            }
        }
    }
    hash_free(&face_nc_map);

    for (s_ncell *nc = dt->head; nc; nc = nc->next)
        for (int j = 0; j < 4; j++)
            if (!dt->point2tet[nc->vertex_id[j]]) dt->point2tet[nc->vertex_id[j]] = nc;

    return 1;
}

/* Commit both half-cavities: remove old Tf tets, insert local-DT tets, fix adjacency.
 * On failure, the DT is inconsistent; the caller must abort. */
static int commit_cavity_expansion(s_scplx *dt,
                                    s_half_cavity *C1, s_scplx *ldt1, int *l2g1,
                                    s_half_cavity *C2, s_scplx *ldt2, int *l2g2,
                                    s_cdt_scratch *sc)
{
    s_hash_table adj_map;
    if (!remove_cavity_tets(dt, C1, C2, &adj_map)) return 0;

    /* Insert real tets from local DTs into global DT. */
    dynarray_clear(&sc->new_tets);

    for (int half = 0; half < 2; half++) {
        s_scplx *ldt = half ? ldt2 : ldt1;
        int     *l2g = half ? l2g2 : l2g1;
        for (s_ncell *lnc = ldt->head; lnc; lnc = lnc->next) {
            bool sent = false;
            for (int j = 0; j < 4; j++) if (lnc->vertex_id[j] < 4) { sent = true; break; }
            if (sent) continue;

            s_ncell *nc = malloc_ncell();
            if (!nc) { hash_free(&adj_map); return 0; }
            memset(nc, 0, sizeof(s_ncell));
            for (int j = 0; j < 4; j++) nc->vertex_id[j] = l2g[lnc->vertex_id[j] - 4];

            nc->next = dt->head; nc->prev = NULL;
            if (dt->head) dt->head->prev = nc;
            dt->head = nc;
            dt->N_ncells++;

            for (int j = 0; j < 4; j++)
                if (!dt->point2tet[nc->vertex_id[j]]) dt->point2tet[nc->vertex_id[j]] = nc;

            if (!dynarray_push(&sc->new_tets, &nc)) {
                hash_free(&adj_map); return 0;
            }
        }
    }

    if (!link_new_tets(dt, &sc->new_tets, &adj_map)) {
        hash_free(&adj_map); return 0;
    }
    hash_free(&adj_map);
    return 1;
}

/* -----------------------------------------------------------------------
 * Step 6 -- Gift-Wrapping Fallback (S.4.4)
 * ----------------------------------------------------------------------- */

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
 *   (a) 4 segments v->{s0,s1,s2,w} vs the triangle interior
 *   (b) 18 boundary-triangle edge vs cone-face {sigma-edge,v}/{apex-edge,v} checks
 *   (c) 9 boundary-triangle vertex inside non-sigma sub-tet {face,v} checks
 * Returns false (not occluded / visible) if no boundary triangle blocks the view. */
static bool is_occluded(int v, int s0, int s1, int s2, int w,
                         const s_dynarray *original_bnd)
{
    for (unsigned i = 0; i < original_bnd->N; i++) {
        const int *u = (const int *)dynarray_get_ptr_c(original_bnd, i);
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

        /* (a) Segments v->{s0,s1,s2,w} vs inner triangle (u0,u1,u2). */
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
     * All u vertices on the non-apex side of sigma -> no intersection possible.
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
     * All T vertices on the same side (or plane) of u -> no intersection. */
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

/* 2D overlap test for two triangles lying in the same plane, by vertex ids.
 * Returns 1 if their interiors overlap improperly: identical triples and
 * triangles sharing only an edge/vertex return 0.  Needed because
 * no_spurious_intersect's half-space early exit cannot see coplanar overlap. */
static int coplanar_tris_improper_overlap(int a0, int a1, int a2,
                                           int b0, int b1, int b2)
{
    int a[3] = {a0, a1, a2}, b[3] = {b0, b1, b2};
    int shared = 0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            if (a[i] == b[j]) { shared++; break; }
    if (shared == 3) return 0;

    for (int i = 0; i < 3; i++) {
        if (cdt_point_in_inner_tri(a[i], b0, b1, b2)) return 1;
        if (cdt_point_in_inner_tri(b[i], a0, a1, a2)) return 1;
    }
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            if (cdt_inner_segs_cross(a[i], a[(i+1)%3], b[j], b[(j+1)%3]))
                return 1;
    return 0;
}

/* a,b,c must already be sorted. void_sign = precomputed orient3d side of the cavity.
 * Only stored on first insertion (count==0 -> 1). */
static int toggle_boundary_face(s_hash_table *toggle, int a, int b, int c, int void_sign)
{
    int f[3] = {a, b, c}; /* caller must pre-sort */
    s_bnd_face_val *p = (s_bnd_face_val *)hash_get_or_create(toggle, f);
    if (!p) return 0;
    if (p->count == 0) p->void_sign = void_sign;
    p->count++;
    return 1;
}

/* Relative orientation of a triangle (f0,f1,f2) lying in the constraint plane
 * wrt the constraint triangle (va,vb,vc): returns k in {-1,+1} such that
 * orient3d(f0,f1,f2,x) == k * orient3d(va,vb,vc,x) for any point x. */
static int coplanar_rel_orientation(int f0, int f1, int f2,
                                     int va_dt, int vb_dt, int vc_dt)
{
    for (int s = 0; s < 4; s++) {  /* sentinels: at least one is off-plane */
        int os = cdt_orient3d(va_dt, vb_dt, vc_dt, s);
        if (os == 0) continue;
        int cs = cdt_orient3d(f0, f1, f2, s);
        return (cs == os) ? 1 : -1;
    }
    return 0; /* unreachable: the 4 sentinels cannot all lie in one plane */
}

/* Gift-wrap one half-cavity.
 * toggle: face -> s_bnd_face_val {count, void_sign}; odd count = active.
 * usable_v: dynarray of int DT-indices of candidate apex vertices.
 * original_bnd: snapshot of the initial boundary (int[3] sorted triples) used
 *   for the occlusion test in condition (iii).
 * comm_out: if non-NULL (upper-cavity fill), active faces lying entirely in
 *   the constraint plane with their void BELOW the plane (this includes the
 *   constraint face itself) are not connected here; they are collected into
 *   comm_out and form the shared middle wall that the lower fill must consume
 *   (paper S.4.4: the triangulation induced on f by C1 is inherited by C2). */
static int gw_half(s_scplx *dt,
                    s_hash_table *toggle,
                    s_dynarray *usable_v,
                    int va_dt, int vb_dt, int vc_dt,
                    s_dynarray *new_tets_out,
                    const s_dynarray *original_bnd,
                    s_dynarray *comm_out,
                    s_cdt_scratch *sc)
{
    int dbg = 0;
    int guard = 0;

    for (;;) {
        if (++guard > 20000) {
            CDT_DBG("gw: pass cap reached, giving up\n");
            return 0;
        }
        dynarray_clear(&sc->active);
        dynarray_clear(&sc->comm_cur);

        for (size_t b = 0; b < toggle->nbuckets; b++) {
            for (s_hash_entry *e = toggle->buckets[b]; e; e = e->next) {
                const s_bnd_face_val *val = (const s_bnd_face_val *)entry_value(toggle, e);
                if (val->count % 2 == 0) continue;
                const int *k = (const int *)entry_key(e);
                if (comm_out &&
                    cdt_orient3d(va_dt, vb_dt, vc_dt, k[0]) == 0 &&
                    cdt_orient3d(va_dt, vb_dt, vc_dt, k[1]) == 0 &&
                    cdt_orient3d(va_dt, vb_dt, vc_dt, k[2]) == 0) {
                    /* On-plane face: defer to the lower fill unless its void
                     * is above the plane (then this fill must consume it). */
                    int krel = coplanar_rel_orientation(k[0], k[1], k[2],
                                                        va_dt, vb_dt, vc_dt);
                    if (val->void_sign != krel) {
                        if (!dynarray_push(&sc->comm_cur, k)) return 0;
                        continue;
                    }
                }
                if (!dynarray_push(&sc->active, k)) return 0;
            }
        }
        if (sc->active.N == 0) {
            if (comm_out) {
                for (unsigned i = 0; i < sc->comm_cur.N; i++) {
                    const int *k = (const int *)dynarray_get_ptr(&sc->comm_cur, i);
                    if (g_dbg_this_face)
                        CDT_DBG("gw: comm face (%d,%d,%d)\n", k[0], k[1], k[2]);
                    if (!dynarray_push(comm_out, k)) return 0;
                }
            }
            return 1;
        }

        bool made_progress = false;
        for (unsigned si = 0; si < sc->active.N && !made_progress; si++) {
            const int *sigma = (const int *)dynarray_get_ptr(&sc->active, si);

            /* int_sign: the cavity (void) side sign for this face.
             * Stored as void_sign in the toggle at seeding/toggle time. */
            const s_bnd_face_val *sv = (const s_bnd_face_val *)hash_get(toggle, sigma);
            if (!sv || sv->void_sign == 0) {
                if (dbg) CDT_DBG("gw: sigma (%d,%d,%d) skipped, void_sign=%d\n",
                                 sigma[0], sigma[1], sigma[2], sv ? sv->void_sign : -99);
                continue;
            }
            int int_sign = sv->void_sign;

            int n_i = 0, n_ii = 0, n_iii = 0;
            int apex = -1;
            for (unsigned vi = 0; vi < usable_v->N; vi++) {
                int w; dynarray_get_value(usable_v, vi, &w);
                if (w < 0) continue;
                /* (i) Correct side. */
                int o_i = cdt_orient3d(sigma[0], sigma[1], sigma[2], w);
                if (dbg) CDT_DBG("gw:   candidate w=%d o_sigma=%d (want %d) o_constraint=%d\n",
                                 w, o_i, int_sign,
                                 cdt_orient3d(va_dt, vb_dt, vc_dt, w));
                if (o_i != int_sign) { n_i++; continue; }

                /* (ii) No spurious intersections with active faces or the
                 * deferred on-plane (comm) faces -- both are cavity walls a
                 * candidate tet must not poke through. */
                bool ok = true;
                for (unsigned ai = 0; ai < sc->active.N && ok; ai++) {
                    const int *u = (const int *)dynarray_get_ptr(&sc->active, ai);
                    if (!no_spurious_intersect(sigma[0], sigma[1], sigma[2], w,
                                               u[0], u[1], u[2]))
                        ok = false;
                }
                for (unsigned ci = 0; ci < sc->comm_cur.N && ok; ci++) {
                    const int *u = (const int *)dynarray_get_ptr(&sc->comm_cur, ci);
                    if (!no_spurious_intersect(sigma[0], sigma[1], sigma[2], w,
                                               u[0], u[1], u[2]))
                        ok = false;
                }

                /* (ii-b) Any new face of the candidate that lies entirely in
                 * the constraint plane must not improperly overlap the middle
                 * wall: neither the constraint face itself nor another
                 * on-plane wall face.  This steers the wall triangulation to
                 * contain the constraint face exactly (our pipeline recovers
                 * individual triangles, not whole coplanar PLC patches). */
                if (ok && cdt_orient3d(va_dt, vb_dt, vc_dt, w) == 0) {
                    for (int nf = 0; nf < 3 && ok; nf++) {
                        int p0 = sigma[nf], p1 = sigma[(nf+1)%3], p2 = w;
                        if (cdt_orient3d(va_dt, vb_dt, vc_dt, p0) != 0) continue;
                        if (cdt_orient3d(va_dt, vb_dt, vc_dt, p1) != 0) continue;
                        if (coplanar_tris_improper_overlap(p0, p1, p2,
                                                           va_dt, vb_dt, vc_dt)) {
                            ok = false; break;
                        }
                        for (unsigned ai = 0; ai < sc->active.N && ok; ai++) {
                            const int *u = (const int *)dynarray_get_ptr(&sc->active, ai);
                            if (cdt_orient3d(va_dt, vb_dt, vc_dt, u[0]) != 0 ||
                                cdt_orient3d(va_dt, vb_dt, vc_dt, u[1]) != 0 ||
                                cdt_orient3d(va_dt, vb_dt, vc_dt, u[2]) != 0)
                                continue;
                            if (coplanar_tris_improper_overlap(p0, p1, p2,
                                                               u[0], u[1], u[2]))
                                ok = false;
                        }
                        for (unsigned ci = 0; ci < sc->comm_cur.N && ok; ci++) {
                            const int *u = (const int *)dynarray_get_ptr(&sc->comm_cur, ci);
                            if (coplanar_tris_improper_overlap(p0, p1, p2,
                                                               u[0], u[1], u[2]))
                                ok = false;
                        }
                    }
                }
                if (!ok) { n_ii++; continue; }

                /* (iii) No VISIBLE vertex inside circumsphere (SoS).
                 * v is visible if (a) it is on the same side of sigma as w,
                 * and (b) no original boundary triangle occludes the line of
                 * sight from v to the tet <sigma, w>. */
                for (unsigned ui = 0; ui < usable_v->N && ok; ui++) {
                    int v; dynarray_get_value(usable_v, ui, &v);
                    if (v < 0 || v==sigma[0] || v==sigma[1] || v==sigma[2] || v==w)
                        continue;
                    int ov = cdt_orient3d(sigma[0], sigma[1], sigma[2], v);
                    if (ov != 0 && ov != int_sign) continue; /* wrong side -> not visible */
                    if (is_occluded(v, sigma[0], sigma[1], sigma[2], w, original_bnd))
                        continue; /* blocked by cavity wall -> not visible */
                    if (cdt_perturbed_insphere(sigma[0], sigma[1], sigma[2], w, v) == 1)
                        ok = false;
                }
                if (!ok) { n_iii++; continue; }

                apex = w;
                break;
            }
            if (apex < 0) {
                if (dbg) CDT_DBG("gw: sigma (%d,%d,%d) void_sign=%d NO APEX "
                                 "(usable=%u rej_side=%d rej_intersect=%d rej_insphere=%d)\n",
                                 sigma[0], sigma[1], sigma[2], int_sign,
                                 usable_v->N, n_i, n_ii, n_iii);
                continue;
            }

            if (g_dbg_this_face)
                CDT_DBG("gw: create tet (%d,%d,%d)+%d (comm_mode=%d void=%d)\n",
                        sigma[0], sigma[1], sigma[2], apex,
                        comm_out != NULL, int_sign);

            /* Create tet and link into dt. */
            s_ncell *nc = malloc_ncell();
            if (!nc) return 0;
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
            if (!dynarray_push(new_tets_out, &nc)) return 0;

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
                !toggle_boundary_face(toggle, f20[0],   f20[1],   f20[2],   vs20))
                return 0;
            made_progress = true;
        }
        if (!made_progress) {
            if (dbg) return 0;
            CDT_DBG("gw: STUCK with %u active faces, rerunning pass with diagnostics\n",
                    sc->active.N);
            dbg = 1;
        }
    }
}

/* Seed toggle and vertex lists for one half-cavity, then run gw_half.
 * comm_out: non-NULL for the upper fill -- receives the middle-wall faces.
 * comm_in:  non-NULL for the lower fill -- middle-wall faces produced by the
 *           upper fill, seeded here with their void on the lower side.
 * Returns 1 on success, 0 on allocation failure. */
static int gw_half_cavity(s_scplx *dt,
                           s_half_cavity *C,
                           int va_dt, int vb_dt, int vc_dt,
                           s_dynarray *new_tets_out,
                           s_dynarray *comm_out,
                           const s_dynarray *comm_in,
                           s_cdt_scratch *sc)
{
    size_t n_seed = C->boundary.N + (comm_in ? comm_in->N : 0);
    s_hash_table toggle;
    if (!hash_init(&toggle, sizeof(int[3]), sizeof(s_bnd_face_val),
                   n_seed * 4 + 1, 0,
                   face_triple_hash, face_triple_eq, NULL)) return 0;
    for (unsigned i = 0; i < C->boundary.N; i++) {
        const int *f = (const int *)dynarray_get_ptr(&C->boundary, i);
        /* The constraint-face entry (f[3]==-1) is never seeded.  The middle
         * wall (f included) is not known a priori: it EMERGES as the on-plane
         * faces created by the upper fill, which are deferred to comm_out and
         * seeded into the lower fill via comm_in. */
        if (f[3] < 0) continue;
        int vs = cdt_orient3d(f[0],f[1],f[2],f[3]);
        if (!toggle_boundary_face(&toggle, f[0], f[1], f[2], vs)) {
            hash_free(&toggle); return 0;
        }
    }
    if (comm_in) {
        for (unsigned i = 0; i < comm_in->N; i++) {
            const int *f = (const int *)dynarray_get_ptr_c(comm_in, i);
            /* Void side of a middle-wall face for the lower fill: below the
             * constraint plane, i.e. -krel (orient3d(f,w)=krel*orient3d(ck,w)). */
            int vs = -coplanar_rel_orientation(f[0], f[1], f[2],
                                               va_dt, vb_dt, vc_dt);
            if (!toggle_boundary_face(&toggle, f[0], f[1], f[2], vs)) {
                hash_free(&toggle); return 0;
            }
        }
    }
    dynarray_clear(&sc->orig_bnd);
    for (unsigned i = 0; i < C->boundary.N; i++) {
        const int *f = (const int *)dynarray_get_ptr(&C->boundary, i);
        if (f[3] < 0) continue; /* occlusion walls: original boundary only */
        int key[3] = {f[0], f[1], f[2]};
        dynarray_push(&sc->orig_bnd, key);
    }
    dynarray_clear(&sc->usable);
    for (unsigned i = 0; i < C->verts.N; i++) {
        int v; dynarray_get_value(&C->verts, i, &v);
        dynarray_push(&sc->usable, &v);
    }
    int ok = gw_half(dt, &toggle, &sc->usable, va_dt, vb_dt, vc_dt,
                     new_tets_out, &sc->orig_bnd, comm_out, sc);
    hash_free(&toggle);
    return ok;
}

/* Gift-wrapping fallback for when cavity expansion fails.
 * C1/C2 still hold original Tf tets; they are removed here.
 * On success, Tf tets are replaced and the constraint face is filled. */
static int gift_wrap_face(s_scplx *dt,
                           s_half_cavity *C1, s_half_cavity *C2,
                           int va_dt, int vb_dt, int vc_dt,
                           s_cdt_scratch *sc)
{
    if (g_dbg_this_face) {
        CDT_DBG("gift_wrap ENTRY constraint=(%d,%d,%d)\n", va_dt, vb_dt, vc_dt);
        for (int half = 0; half < 2; half++) {
            s_half_cavity *C = half ? C2 : C1;
            CDT_DBG("C%d: %u tets, %u verts, %u bnd\n", half+1, C->tets.N, C->verts.N, C->boundary.N);
            for (unsigned i = 0; i < C->tets.N; i++) {
                s_ncell *nc; dynarray_get_value(&C->tets, i, &nc);
                CDT_DBG("  tet (%d[%d],%d[%d],%d[%d],%d[%d])\n",
                        nc->vertex_id[0], cdt_orient3d(va_dt,vb_dt,vc_dt,nc->vertex_id[0]),
                        nc->vertex_id[1], cdt_orient3d(va_dt,vb_dt,vc_dt,nc->vertex_id[1]),
                        nc->vertex_id[2], cdt_orient3d(va_dt,vb_dt,vc_dt,nc->vertex_id[2]),
                        nc->vertex_id[3], cdt_orient3d(va_dt,vb_dt,vc_dt,nc->vertex_id[3]));
            }
            for (unsigned i = 0; i < C->boundary.N; i++) {
                int *f = (int *)dynarray_get_ptr(&C->boundary, i);
                CDT_DBG("  bnd (%d[%d],%d[%d],%d[%d]) int_v=%d\n",
                        f[0], cdt_orient3d(va_dt,vb_dt,vc_dt,f[0]),
                        f[1], cdt_orient3d(va_dt,vb_dt,vc_dt,f[1]),
                        f[2], cdt_orient3d(va_dt,vb_dt,vc_dt,f[2]), f[3]);
            }
        }
    }

    s_hash_table adj_map;
    if (!remove_cavity_tets(dt, C1, C2, &adj_map)) return 0;

    dynarray_clear(&sc->new_tets);

    /* Upper fill first: it produces the middle-wall triangulation (comm),
     * which the lower fill inherits (paper S.4.4). */
    s_dynarray comm = dynarray_initialize(sizeof(int[3]), 16);
    if (!comm.items) { hash_free(&adj_map); return 0; }

    if (!gw_half_cavity(dt, C1, va_dt, vb_dt, vc_dt, &sc->new_tets, &comm, NULL, sc)) {
        CDT_DBG("gift_wrap: C1 failed (tets=%u verts=%u bnd=%u)\n",
                C1->tets.N, C1->verts.N, C1->boundary.N);
        dynarray_free(&comm); hash_free(&adj_map); return 0;
    }
    if (!gw_half_cavity(dt, C2, va_dt, vb_dt, vc_dt, &sc->new_tets, NULL, &comm, sc)) {
        CDT_DBG("gift_wrap: C2 failed (tets=%u verts=%u bnd=%u comm=%u)\n",
                C2->tets.N, C2->verts.N, C2->boundary.N, comm.N);
        dynarray_free(&comm); hash_free(&adj_map); return 0;
    }
    dynarray_free(&comm);

    if (!link_new_tets(dt, &sc->new_tets, &adj_map)) {
        hash_free(&adj_map); return 0;
    }
    hash_free(&adj_map);
    return 1;
}

/* -----------------------------------------------------------------------
 * Step 7 -- Ghost-Vertex Flood-Fill for Interior/Exterior (S.4.5)
 * ----------------------------------------------------------------------- */

static inline int ncell_has_sentinel(const s_ncell *nc)
{
    for (int j = 0; j < 4; j++) if (nc->vertex_id[j] < 4) return 1;
    return 0;
}

/* Tets removed by classify's Phase 4:
 *   keep_exterior=false: every externally-labelled tet (ghosts + pockets).
 *   keep_exterior=true : only ghost/sentinel tets, so the surviving complex is
 *                        the full convex-hull tetrahedralization (pockets kept,
 *                        flagged interior=false). */
static inline int ncell_is_deleted(const s_ncell *nc, bool keep_exterior, int ext_stamp)
{
    return keep_exterior ? ncell_has_sentinel(nc) : (nc->mark_token == ext_stamp);
}

static void classify_interior_exterior(s_scplx *dt, s_hash_table *constrained,
                                       s_cdt_scratch *sc, bool keep_exterior)
{
    int ext_stamp = ++dt->mark_stamp;
    int int_stamp = ++dt->mark_stamp;

    dynarray_clear(&sc->queue);

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
            int fids[3]; FACE_IDS(nc, fi, fids);
            sort3(&fids[0], &fids[1], &fids[2]);
            nb->mark_token = is_constrained_face(constrained, fids[0], fids[1], fids[2])
                             ? int_stamp : ext_stamp;
            dynarray_push(&sc->queue, &nb);
        }
    }

    /* Phase 3: BFS -- propagate label, flip at constrained faces. */
    while (sc->queue.N > 0) {
        s_ncell *cur; dynarray_pop(&sc->queue, &cur);
        int cur_label = cur->mark_token;
        for (int fi = 0; fi < 4; fi++) {
            s_ncell *nb = cur->opposite[fi];
            if (!nb || nb->mark_token == ext_stamp || nb->mark_token == int_stamp) continue;
            int fids[3]; FACE_IDS(cur, fi, fids);
            sort3(&fids[0], &fids[1], &fids[2]);
            nb->mark_token = is_constrained_face(constrained, fids[0], fids[1], fids[2])
                             ? (cur_label == ext_stamp ? int_stamp : ext_stamp)
                             : cur_label;
            dynarray_push(&sc->queue, &nb);
        }
    }

    /* Phase 4: record the interior flag on every tet, then delete tets per the
     * keep_exterior policy (ghosts only, or ghosts + pockets). */
    s_ncell *nc = dt->head;
    while (nc) {
        s_ncell *next = nc->next;
        nc->interior = (nc->mark_token == int_stamp);
        if (ncell_is_deleted(nc, keep_exterior, ext_stamp)) {
            /* Null out opposite pointers in surviving neighbors. */
            for (int fi = 0; fi < 4; fi++) {
                s_ncell *nb = nc->opposite[fi];
                if (!nb || ncell_is_deleted(nb, keep_exterior, ext_stamp)) continue;
                for (int k = 0; k < 4; k++)
                    if (nb->opposite[k] == nc) { nb->opposite[k] = NULL; break; }
            }
            /* Fix point2tet: repoint to any surviving incident neighbor. */
            for (int j = 0; j < 4; j++) {
                int v = nc->vertex_id[j];
                if (v < 4) { dt->point2tet[v] = NULL; continue; }
                if (dt->point2tet[v] != nc) continue;
                s_ncell *rep = NULL;
                for (int fi = 0; fi < 4 && !rep; fi++) {
                    s_ncell *nb = nc->opposite[fi];
                    if (nb && !ncell_is_deleted(nb, keep_exterior, ext_stamp)) rep = nb;
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

/* Split face fi of mesh into 3 subfaces at its centroid.  The centroid lies
 * on the input surface, so the described domain is unchanged.
 * Returns 1 on success, 0 on allocation failure. */
static int split_trimesh_face_centroid(s_trimesh *mesh, int fi)
{
    int v0 = mesh->faces[fi*3], v1 = mesh->faces[fi*3+1], v2 = mesh->faces[fi*3+2];
    s_point p0 = mesh->points.p[v0], p1 = mesh->points.p[v1], p2 = mesh->points.p[v2];
    s_point c = { .x = (p0.x + p1.x + p2.x) / 3.0,
                  .y = (p0.y + p1.y + p2.y) / 3.0,
                  .z = (p0.z + p1.z + p2.z) / 3.0 };

    s_point *np = realloc(mesh->points.p, (size_t)(mesh->points.N + 1) * sizeof(s_point));
    if (!np) return 0;
    np[mesh->points.N] = c;
    mesh->points.p = np;
    int cid = mesh->points.N++;

    int new_Nf = mesh->Nf + 2;
    int     *nf  = realloc(mesh->faces,    (size_t)new_Nf * 3 * sizeof(int));
    if (!nf) return 0;
    mesh->faces = nf;
    s_point *nn  = realloc(mesh->fnormals, (size_t)new_Nf     * sizeof(s_point));
    if (!nn) return 0;
    mesh->fnormals = nn;

    s_point fn = mesh->fnormals[fi];
    /* Replace fi by (v0,v1,c); append (v1,v2,c) and (v2,v0,c). */
    mesh->faces[fi*3]   = v0; mesh->faces[fi*3+1] = v1; mesh->faces[fi*3+2] = cid;
    mesh->faces[mesh->Nf*3]   = v1; mesh->faces[mesh->Nf*3+1] = v2; mesh->faces[mesh->Nf*3+2] = cid;
    mesh->fnormals[mesh->Nf]  = fn;
    mesh->Nf++;
    mesh->faces[mesh->Nf*3]   = v2; mesh->faces[mesh->Nf*3+1] = v0; mesh->faces[mesh->Nf*3+2] = cid;
    mesh->fnormals[mesh->Nf]  = fn;
    mesh->Nf++;

    free(mesh->adjacency);
    mesh->adjacency = NULL;
    return 1;
}

static s_scplx cdt_attempt(s_trimesh *working, double EPS_DEG, double TOL,
                            s_cdt_scratch *sc, int *fail_fi, bool keep_exterior);

/* Shared body for both public entries. keep_exterior selects the classify
 * policy (see classify_interior_exterior / ncell_is_deleted). */
static s_scplx cdt_build(const s_trimesh *mesh, double EPS_DEG, double TOL,
                         bool keep_exterior)
{
    if (!trimesh_is_valid(mesh)) return (s_scplx){0};

    s_trimesh working = copy_trimesh(mesh);
    if (!trimesh_is_valid(&working)) return (s_scplx){0};

    s_cdt_scratch sc = {0};
    sc.ring     = dynarray_initialize(sizeof(s_ncell *),   16);
    sc.out_ids  = dynarray_initialize(sizeof(int),         16);
    sc.stack    = dynarray_initialize(sizeof(s_ncell *),   64);
    sc.in_R     = dynarray_initialize(sizeof(s_ncell *),   64);
    sc.queue    = dynarray_initialize(sizeof(s_ncell *),   64);
    sc.new_tets = dynarray_initialize(sizeof(s_ncell *),   64);
    sc.orig_bnd = dynarray_initialize(sizeof(int[3]),      32);
    sc.usable   = dynarray_initialize(sizeof(int),         32);
    sc.active   = dynarray_initialize(sizeof(int[3]),      16);
    sc.lnc_info = dynarray_initialize(sizeof(s_lnc_info),  64);
    sc.to_split = dynarray_initialize(sizeof(int[2]),      64);
    sc.comm_cur = dynarray_initialize(sizeof(int[3]),      16);
    if (!sc.ring.items || !sc.out_ids.items || !sc.stack.items || !sc.in_R.items ||
        !sc.queue.items || !sc.new_tets.items || !sc.orig_bnd.items || !sc.usable.items ||
        !sc.active.items || !sc.lnc_info.items || !sc.to_split.items || !sc.comm_cur.items) {
        cdt_scratch_free(&sc); free_trimesh(&working); return (s_scplx){0};
    }

    /* Steiner-split retry loop: if one face resists recovery (near-degenerate
     * local configuration), split it at its centroid -- a new point ON the
     * input surface, so the domain geometry is unchanged -- and rebuild.
     * (TetGen-style fallback; each round only perturbs the failing spot.)
     * BACKSTOP ONLY: with the exact one-geometry CDT this should never fire.
     * A split is a LOUD warning -- it flags a genuine recovery bug to fix at
     * source, not a routine degeneracy to paper over. */
    s_scplx result = {0};
    const int MAX_ROUNDS = 8;
    for (int round = 0; round < MAX_ROUNDS; round++) {
        int fail_fi = -1;
        result = cdt_attempt(&working, EPS_DEG, TOL, &sc, &fail_fi, keep_exterior);
        if (result.head) break;
        if (fail_fi < 0) break;   /* failure not tied to a specific face */
        fprintf(stderr,
            "\n**** UNINTENDED CDT FALLBACK (band-aid) ****\n"
            "  The exact one-geometry CDT FAILED to recover constraint face "
            "fi=%d and is\n  splitting it at its centroid (round %d/8) to force "
            "recovery.\n  This is NOT expected and indicates a real recovery "
            "bug that should be\n  fixed at the source, not papered over.  The "
            "output remains geometrically\n  correct but the triangulation now "
            "contains extra Steiner points on the\n  surface (a refinement of "
            "the input mesh).  PLEASE REPORT with the input.\n"
            "********************************************\n\n",
            fail_fi, round);
        if (!split_trimesh_face_centroid(&working, fail_fi)) break;
    }

    cdt_scratch_free(&sc);
    free_trimesh(&working);

    /* The build ran in exact_ids mode against the cdt_predicates registry, which
     * is cleared once construction finishes. Reset the returned complex to
     * coordinate predicates so it is safe to walk / point-locate afterwards
     * (dtp_orient in exact mode would dereference the freed registry). */
    result.exact_ids = 0;
    return result;
}

s_scplx tetrahedralize_interior_trimesh(const s_trimesh *mesh,
                                         double EPS_DEG, double TOL)
{
    return cdt_build(mesh, EPS_DEG, TOL, /*keep_exterior=*/false);
}

s_scplx tetrahedralize_domain_flagged(const s_trimesh *mesh,
                                       double EPS_DEG, double TOL)
{
    return cdt_build(mesh, EPS_DEG, TOL, /*keep_exterior=*/true);
}

/* One full CDT attempt (Phases A-C) on the given working mesh.  On failure
 * returns an empty complex; if the failure is pinned to one constraint face,
 * *fail_fi is set to its index in working->faces (else left at -1).
 * The working mesh is mutated (Phase A edge splits, point sync) but remains a
 * valid refinement of the same surface, so the caller can retry with it. */
static s_scplx cdt_attempt(s_trimesh *working, double EPS_DEG, double TOL,
                            s_cdt_scratch *sc, int *fail_fi, bool keep_exterior)
{
    /* Phase A: DT + segment recovery, in exact-id mode -- ONE exact geometry
     * end-to-end.  The builder and Phase B (blocker test, cavity expansion,
     * gift-wrap) all reason on the cdt_predicates registry by vertex id, over
     * the true implicit LNC Steiner positions.  (The float DT builder remains
     * for the Voronoi seed / weighted paths; only the CDT is exact.) */
    cdt_predicates_init();  /* FPU rounding mode; idempotent */
    s_dt_builder b = dt_builder_begin_exact(&working->points, TOL, NULL, NULL);
    if (!b._stack) return (s_scplx){0};

    if (!ensure_ridge_protected(&b, working, EPS_DEG, TOL, sc)) {
        s_scplx tmp = dt_builder_end(&b, false, NULL, NULL, 0);
        free_complex(&tmp); return (s_scplx){0};
    }

    /* Keep the big tetra: the CDT relies on sentinels staying at indices 0-3
     * (real vertices at working_vertex+4) for the +4 indexing convention and
     * for Phase C's ghost-tet detection (vertex_id < 4).  Building with
     * keep_big_tetra=false would compact point indices to 0-based and shrink
     * point2tet, making every +4 access out of bounds. */
    s_scplx dt = dt_builder_end(&b, true, NULL, NULL, 0);
    if (!dt.head) return (s_scplx){0};

    /* Sync working->points with the DT point set (0-based: skip sentinels),
     * so faces indices stay aligned and Phase A Steiner coords are present. */
    {
        int n_real = dt.points.N - 4;
        s_point *p = malloc((size_t)n_real * sizeof(s_point));
        if (!p) { free_complex(&dt); return (s_scplx){0}; }
        memcpy(p, dt.points.p + 4, (size_t)n_real * sizeof(s_point));
        free(working->points.p);
        working->points.p = p;
        working->points.N = n_real;
    }

    /* Phase B: face recovery (Steps 4-6) -- Algorithm 2 outer loop. */
    s_hash_table constrained;
    if (!hash_init(&constrained, sizeof(int[3]), sizeof(int),
                   (size_t)(working->Nf * 4 + 1), (size_t)(working->Nf * 2),
                   face_triple_hash, face_triple_eq, NULL)) {
        free_complex(&dt); return (s_scplx){0};
    }

#define PHASE_B_ABORT() do { \
        half_cavity_free(&C1); half_cavity_free(&C2); \
        *fail_fi = fi; \
        hash_free(&constrained); \
        free_complex(&dt); cdt_predicates_clear(); \
        return (s_scplx){0}; \
    } while (0)

    int dbg_n_recovered = 0;
    int pass = 0;
    bool any_missing = true;
    while (any_missing) {
        any_missing = false;

        /* Guard against surgeries ping-ponging (face A's recovery breaking
         * face B and vice versa): give up on the first missing face so the
         * caller splits it. */
        if (++pass > 30) {
            for (int fi = 0; fi < working->Nf; fi++) {
                int va = working->faces[fi*3]   + 4;
                int vb = working->faces[fi*3+1] + 4;
                int vc = working->faces[fi*3+2] + 4;
                if (!face_in_dt(&dt, va, vb, vc, sc)) {
                    CDT_DBG("PHASE B: no convergence after %d passes; "
                            "giving up on fi=%d\n", pass, fi);
                    *fail_fi = fi;
                    break;
                }
            }
            hash_free(&constrained);
            free_complex(&dt); cdt_predicates_clear();
            return (s_scplx){0};
        }

        for (int fi = 0; fi < working->Nf; fi++) {
            /* working->faces is 0-based; DT vertex indices are +4 (sentinels at 0-3). */
            int va = working->faces[fi*3]   + 4;
            int vb = working->faces[fi*3+1] + 4;
            int vc = working->faces[fi*3+2] + 4;

            if (face_in_dt(&dt, va, vb, vc, sc)) {
                mark_constrained_face(&constrained, va, vb, vc);
                continue;
            }
            any_missing = true;

            {
                static int dbg_face[3] = {-2, -2, -2};
                if (dbg_face[0] == -2) {
                    dbg_face[0] = -1;
                    const char *s = getenv("CDT_DBG_FACE");
                    if (s) sscanf(s, "%d,%d,%d", &dbg_face[0], &dbg_face[1], &dbg_face[2]);
                }
                g_dbg_this_face = (va==dbg_face[0] && vb==dbg_face[1] && vc==dbg_face[2]);
            }

            if (g_dbg_this_face) {
                CDT_DBG("face (%d,%d,%d) missing; edges in DT: ab=%d bc=%d ca=%d\n",
                        va, vb, vc,
                        edge_in_dt(&dt, va, vb, sc),
                        edge_in_dt(&dt, vb, vc, sc),
                        edge_in_dt(&dt, vc, va, sc));
            }

            s_half_cavity C1 = {0}, C2 = {0};
            if (!build_half_cavities(&dt, va, vb, vc, &constrained, &C1, &C2, sc))
                PHASE_B_ABORT();

            s_scplx ldt1 = {0}, ldt2 = {0};
            int *l2g1 = NULL, *l2g2 = NULL, n1 = 0, n2 = 0;

            int ok1 = try_cavity_expansion(&dt, &C1, va, vb, vc, +1, TOL, &constrained,
                                           &ldt1, &l2g1, &n1, sc);
            int ok2 = ok1
                      ? try_cavity_expansion(&dt, &C2, va, vb, vc, -1, TOL, &constrained,
                                             &ldt2, &l2g2, &n2, sc)
                      : 0;

            bool face_ok;
            if (ok1 && ok2) {
                int commit_ok = commit_cavity_expansion(&dt,
                                    &C1, &ldt1, l2g1,
                                    &C2, &ldt2, l2g2, sc);
                face_ok = commit_ok;
            } else {
                /* try_cavity_expansion GROWS C1/C2 in place (appends verts/tets,
                 * rebuilds boundary) -- so after a partial success (e.g. ok1=1,
                 * ok2=0) they no longer hold the original Tf region gift_wrap_face
                 * expects; an expansion that crossed a mesh-surface face even
                 * pulls a ghost/sentinel vertex into the usable set and dead-ends
                 * the gift-wrap.  Rebuild the pristine half-cavities first. */
                half_cavity_free(&C1); half_cavity_free(&C2);
                C1 = (s_half_cavity){0}; C2 = (s_half_cavity){0};
                if (!build_half_cavities(&dt, va, vb, vc, &constrained, &C1, &C2, sc))
                    PHASE_B_ABORT();
                face_ok = gift_wrap_face(&dt, &C1, &C2, va, vb, vc, sc);
            }
            if (ldt1.head) free_complex(&ldt1); free(l2g1);
            if (ldt2.head) free_complex(&ldt2); free(l2g2);
            /* Safety net: a "successful" recovery that did not actually
             * produce the face would make the outer loop retry the same
             * cavity forever.  Abort loudly instead. */
            if (face_ok && !face_in_dt(&dt, va, vb, vc, sc)) {
                CDT_DBG("PHASE B: face fi=%d (%d,%d,%d) still missing after "
                        "recovery reported success (ok1=%d ok2=%d)\n",
                        fi, va, vb, vc, ok1, ok2);
                face_ok = false;
            }
            {
                static int validate = -1;
                if (validate < 0) validate = getenv("CDT_VALIDATE") ? 1 : 0;
                if (validate && face_ok) {
                    int bad = validate_complex(&dt);
                    if (bad) {
                        CDT_DBG("VALIDATE: %d adjacency violations after recovering "
                                "fi=%d (%d,%d,%d) via %s\n",
                                bad, fi, va, vb, vc,
                                (ok1 && ok2) ? "commit" : "gift-wrap");
                        face_ok = false;
                    }
                }
            }
            if (face_ok) {
                dbg_n_recovered++;
                mark_constrained_face(&constrained, va, vb, vc);
            } else {
                CDT_DBG("PHASE B ABORT: fi=%d (va,vb,vc)=(%d,%d,%d) ok1=%d ok2=%d "
                        "(%d faces surgically recovered before this)\n",
                        fi, va, vb, vc, ok1, ok2, dbg_n_recovered);
                PHASE_B_ABORT();
            }

            half_cavity_free(&C1);
            half_cavity_free(&C2);
        }
    }

#undef PHASE_B_ABORT

    /* Phase C: interior/exterior classification (Step 7) */
    classify_interior_exterior(&dt, &constrained, sc, keep_exterior);

    hash_free(&constrained);
    cdt_predicates_clear();

    return dt;
}
