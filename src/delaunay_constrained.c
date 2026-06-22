#include "delaunay.h"
#include "scplx.h"
#include "trimesh.h"
#include "dynarray.h"
#include "hash.h"
#include "gtests.h"
#include "points.h"
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

/* Ensure every trimesh boundary edge is present in the DT.
 * mesh->faces holds 0-based canonical indices; the +4 builder offset is
 * applied only at the DT lookup call sites here — the trimesh is never shifted. */
static int ensure_ridge_protected(s_dt_builder *b, s_trimesh *mesh, double TOL)
{
    bool changed = true;
    while (changed) {
        changed = false;

        s_ncell **vertex_start = calloc((size_t)b->dt.points.N, sizeof(s_ncell *));
        int      *v_local_id   = calloc((size_t)b->dt.points.N, sizeof(int));
        if (!vertex_start || !v_local_id) {
            free(vertex_start); free(v_local_id); return 0;
        }
        build_vertex_cell_index(&b->dt, vertex_start, v_local_id);

        for (int fi = 0; fi < mesh->Nf && !changed; fi++) {
            for (int ei = 0; ei < 3 && !changed; ei++) {
                int va = mesh->faces[fi*3 + (ei+1)%3];  /* canonical (0-based) */
                int vb = mesh->faces[fi*3 + (ei+2)%3];

                /* Builder indices are canonical + 4 (sentinels occupy 0..3). */
                if (edge_in_dt(&b->dt, vertex_start, v_local_id, va + 4, vb + 4))
                    continue;

                s_point pa = b->dt.points.p[va + 4];
                s_point pb = b->dt.points.p[vb + 4];
                s_point M  = { .x = (pa.x + pb.x) * 0.5,
                               .y = (pa.y + pb.y) * 0.5,
                               .z = (pa.z + pb.z) * 0.5 };

                /* new_id is in builder space; canonical = new_id - 4. */
                int new_id_builder = b->dt.points.N;
                s_points single = { .N = 1, .p = &M };
                if (!dt_builder_extend(b, &single, TOL)) {
                    free(vertex_start); free(v_local_id); return 0;
                }

                if (!split_trimesh_edge(mesh, va, vb, new_id_builder - 4)) {
                    free(vertex_start); free(v_local_id); return 0;
                }

                changed = true;
                free(vertex_start);
                free(v_local_id);
            }
        }

        if (!changed) {
            free(vertex_start);
            free(v_local_id);
        }
    }
    return 1;
}


/* ----------------------------------------------------------------------- */

static s_ncell *find_tet_with_edge(const s_scplx *dt, int va, int vb,
                                    int *out_va_local, int *out_vb_local)
{
    for (s_ncell *nc = dt->head; nc; nc = nc->next) {
        int la = -1, lb = -1;
        for (int j = 0; j < 4; j++) {
            if (nc->vertex_id[j] == va) la = j;
            if (nc->vertex_id[j] == vb) lb = j;
        }
        if (la >= 0 && lb >= 0) {
            if (out_va_local) *out_va_local = la;
            if (out_vb_local) *out_vb_local = lb;
            return nc;
        }
    }
    return NULL;
}

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

static int flipinsertfacet(s_scplx *dt, int va, int vb, int vc,
                            s_hash_table *constrained)
{
    if (face_in_dt(dt, va, vb, vc)) {
        mark_constrained_face(constrained, va, vb, vc);
        return 1;
    }

    s_point tri[3] = {dt->points.p[va], dt->points.p[vb], dt->points.p[vc]};

    for (int iter = 0; iter < 2000; iter++) {
        if (face_in_dt(dt, va, vb, vc)) {
            mark_constrained_face(constrained, va, vb, vc);
            return 1;
        }

        int eu = -1, ev = -1;
        for (s_ncell *nc = dt->head; nc && eu < 0; nc = nc->next) {
            for (int i = 0; i < 3 && eu < 0; i++) {
                for (int j = i+1; j < 4 && eu < 0; j++) {
                    int u = nc->vertex_id[i], v = nc->vertex_id[j];
                    if (u==va || u==vb || u==vc || v==va || v==vb || v==vc) continue;
                    s_point seg[2] = {dt->points.p[u], dt->points.p[v]};
                    e_intersect_type r =
                        test_segment_triangle_intersect_3D(seg, tri, 0.0, 0.0);
                    if (r == INTERSECT_NONDEGENERATE) {
                        eu = u; ev = v;
                    } else if (r == INTERSECT_DEGENERATE) {
                        fprintf(stderr, "flipinsertfacet: degenerate crossing edge "
                                "(%d,%d) vs face (%d,%d,%d) — coplanar input?\n",
                                u, v, va, vb, vc);
                    } else {
                        fprintf(stderr, "flipinsertfacet: whoah! Unconsidered intersection type.\n");
                    }
                }
            }
        }

        if (eu < 0) {
            fprintf(stderr, "flipinsertfacet: no crossing edges but face (%d,%d,%d) "
                    "not in DT\n", va, vb, vc);
            return 0;
        }

        int eu_l, ev_l;
        s_ncell *start = find_tet_with_edge(dt, eu, ev, &eu_l, &ev_l);
        if (!start) continue;

        bool flipped = false;
        s_ncell *cur  = start;
        int cur_eu = eu_l, cur_ev = ev_l;
        int fan_steps = 0;

        do {
            for (int opp = 0; opp < 4 && !flipped; opp++) {
                if (opp == cur_eu || opp == cur_ev) continue;
                if (!cur->opposite[opp]) continue;

                int fids[3]; extract_ids_face(cur, 2, &opp, fids);
                if (is_constrained_face(constrained, fids[0], fids[1], fids[2])) continue;

                if (dt_flip_face(dt, cur, opp) == 1) flipped = true;
            }

            if (flipped) break;

            int new_eu, new_ev;
            s_ncell *nxt = next_ncell_ridge_cycle(cur, cur_eu, cur_ev, &new_eu, &new_ev);
            if (!nxt || nxt == start) break;
            cur = nxt; cur_eu = new_eu; cur_ev = new_ev;
        } while (++fan_steps < 64);

        if (!flipped) {
            fprintf(stderr, "flipinsertfacet: no flippable face around crossing edge "
                    "(%d,%d) for constraining face (%d,%d,%d)\n", eu, ev, va, vb, vc);
            return 0;
        }
    }

    if (face_in_dt(dt, va, vb, vc)) {
        mark_constrained_face(constrained, va, vb, vc);
        return 1;
    }
    fprintf(stderr, "flipinsertfacet: iteration limit reached for face (%d,%d,%d)\n",
            va, vb, vc);
    return 0;
}


/* ---------------------------------------------------------------------- */

static int restore_delaunay(s_scplx *dt, s_hash_table *constrained)
{
    bool progress = true;
    int rounds_left = dt->N_ncells * 20 + 100;

    while (progress && rounds_left-- > 0) {
        progress = false;
        for (s_ncell *nc = dt->head; nc; nc = nc->next) {
            for (int fi = 0; fi < 4; fi++) {
                if (!nc->opposite[fi]) continue;

                int fids[3]; extract_ids_face(nc, 2, &fi, fids);
                if (is_constrained_face(constrained, fids[0], fids[1], fids[2])) continue;

                if (!are_locally_delaunay(dt, nc, fi, DELAUNAY_TEST_NONSTRICT)) {
                    if (dt_flip_face(dt, nc, fi) == 1) {
                        progress = true;
                        goto restart;
                    }
                }
            }
        }
        restart:;
    }
    return 1;
}


/* -----------------------------------------------------------------------
 * Public entry point
 * ----------------------------------------------------------------------- */

#define POINT_IN_TRIMESH_MAX_RETRIES 10

s_scplx tetrahedralize_interior_trimesh(const s_trimesh *mesh,
                                         double EPS_DEG, double TOL)
{
    if (!trimesh_is_valid(mesh)) return (s_scplx){0};

    /* Build initial DT with big tetra kept for incremental Steiner insertion. */
    s_dt_builder b = dt_builder_begin(&mesh->points, NULL, TOL, NULL, NULL);
    if (!b._stack) return (s_scplx){0};

    /* Mutable working copy: ridge protection may split boundary edges. */
    s_trimesh working = copy_trimesh(mesh);
    if (!trimesh_is_valid(&working)) {
        s_scplx tmp = dt_builder_end(&b, false, NULL, NULL, 0);
        free_complex(&tmp);
        return (s_scplx){0};
    }

    /* Sub-phase A: ridge protection using incremental insertion.
     * ensure_ridge_protected applies the +4 builder offset internally. */
    if (!ensure_ridge_protected(&b, &working, TOL)) {
        s_scplx tmp = dt_builder_end(&b, false, NULL, NULL, 0);
        free_complex(&tmp);
        free_trimesh(&working);
        return (s_scplx){0};
    }

    /* Strip sentinels and compact. Assuming no seed duplicates, compaction
     * shifts all non-sentinel indices by exactly -4, matching the canonical
     * 0-based indices already stored in working.faces. */
    s_scplx dt = dt_builder_end(&b, false, NULL, NULL, 0);
    if (!dt.head) { free_trimesh(&working); return (s_scplx){0}; }

    /* Sub-phase B: force each boundary face into the DT. */
    s_hash_table constrained;
    if (!hash_init(&constrained, sizeof(int[3]), sizeof(int),
                   (size_t)(working.Nf * 2 + 1), (size_t)working.Nf,
                   face_triple_hash, face_triple_eq, NULL)) {
        free_complex(&dt); free_trimesh(&working); return (s_scplx){0};
    }

    for (int fi = 0; fi < working.Nf; fi++) {
        int va = working.faces[fi*3], vb = working.faces[fi*3+1], vc = working.faces[fi*3+2];
        if (!flipinsertfacet(&dt, va, vb, vc, &constrained)) {
            hash_free(&constrained); free_complex(&dt); free_trimesh(&working);
            return (s_scplx){0};
        }
    }

    /* Sub-phase C: restore Delaunay on non-constrained faces. */
    restore_delaunay(&dt, &constrained);
    hash_free(&constrained);
    free_trimesh(&working);

    /* Filter to interior tets by centroid test.
     * For each exterior tet: null back-references in its interior neighbours
     * (so adjacency among kept tets stays valid), then unlink and free it. */
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

        /* Null back-references from interior neighbours. */
        for (int fi = 0; fi < 4; fi++) {
            s_ncell *nb = nc->opposite[fi];
            if (!nb) continue;
            for (int fj = 0; fj < 4; fj++)
                if (nb->opposite[fj] == nc) { nb->opposite[fj] = NULL; break; }
        }

        /* Unlink from doubly-linked list. */
        if (nc->prev) nc->prev->next = next;
        else          dt.head = next;
        if (next)     next->prev = nc->prev;
        dt.N_ncells--;

        free_ncell(nc);
        nc = next;
    }

    return dt;
}
