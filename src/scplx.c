
#include "scplx.h"
#include "points.h"
#include "gtests.h"
#include "dt_predseam.h"  /* Phase 1: id-based predicate seam (walk-by-id) */
#include "gnuplotc.h"
#include "dynarray.h"
#include "random.h"       /* caller-owned PRNG for deterministic walk tie-breaking */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <float.h>

static int id_where_equal_int(const int *arr, int N, int entry) 
{
    for (int ii=0; ii<N; ii++) if (arr[ii] == entry) return ii;
    fprintf(stderr, "id_where_equal_int: Could not find id.\n"); 
    assert(1==0);
    return -10;
}

static int inarray(const int *arr1, int N, int a)
{
    for (int ii=0; ii<N; ii++) if (arr1[ii] == a) return 1;
    return 0;
}



/* ---- ncell slab pool (see the pool field in scplx.h) ---------------------
 * Bump allocation from fixed-size blocks + a freelist threaded through the
 * `next` field of dead cells.  Blocks never move, so every s_ncell* stays
 * stable for the complex's lifetime; free_complex releases whole blocks.
 * Cells allocated consecutively are memory-adjacent, so insertion-order
 * traversals (the head->next walks all over dt/cdt/medax) get cache locality
 * that per-cell malloc cannot provide. */
#ifndef NCELL_NO_POOL

#define NCELL_BLOCK 1024

typedef struct ncell_block {
    struct ncell_block *next;
    int used;                      /* bump cursor into cells[] */
    s_ncell cells[NCELL_BLOCK];
} s_ncell_block;

typedef struct ncell_pool {
    s_ncell_block *blocks;         /* newest block first */
    s_ncell       *free_head;      /* freelist of returned cells */
} s_ncell_pool;

s_ncell *malloc_ncell(s_scplx *scplx)
{
    s_ncell_pool *P = scplx->pool;
    if (!P) {
        P = calloc(1, sizeof *P);
        assert(P && "Could not malloc ncell pool");
        scplx->pool = P;
    }
    s_ncell *out;
    if (P->free_head) {
        out = P->free_head;
        P->free_head = out->next;
    } else {
        if (!P->blocks || P->blocks->used == NCELL_BLOCK) {
            s_ncell_block *b = malloc(sizeof *b);
            assert(b && "Could not malloc ncell block");
            b->used = 0;
            b->next = P->blocks;
            P->blocks = b;
        }
        out = &P->blocks->cells[P->blocks->used++];
    }
    memset(out, 0, sizeof(s_ncell));
    return out;
}

void free_ncell(s_scplx *scplx, s_ncell *ncell)
{
    s_ncell_pool *P = scplx->pool;
    ncell->next = P->free_head;
    P->free_head = ncell;
}

static void free_ncell_storage(s_scplx *scplx)
{
    if (!scplx->pool) return;
    s_ncell_block *b = scplx->pool->blocks;
    while (b) {
        s_ncell_block *next = b->next;
        free(b);
        b = next;
    }
    free(scplx->pool);
}

#else  /* NCELL_NO_POOL: plain per-cell malloc, ASAN/leaks-friendly */

s_ncell *malloc_ncell(s_scplx *scplx)
{
    (void)scplx;
    s_ncell *out = malloc(sizeof(s_ncell));
    assert(out && "Could not malloc ncell");
    memset(out, 0, sizeof(s_ncell));
    return out;
}

void free_ncell(s_scplx *scplx, s_ncell *ncell)
{
    (void)scplx;
    free(ncell);
}

static void free_ncell_storage(s_scplx *scplx)
{   /* only the still-live cells are in the list; freed ones are already gone */
    s_ncell *current = scplx->head;
    while (current) {
        s_ncell *next = current->next;
        free(current);
        current = next;
    }
}

#endif


void free_complex(s_scplx *scplx)
{
    free_points(&scplx->points);
    if (scplx->weights) free(scplx->weights);
    if (scplx->point2tet) free(scplx->point2tet);

    free_ncell_storage(scplx);

    memset(scplx, 0, sizeof(s_scplx));
}


void assert_point2tet(const s_scplx *scplx)
{
    if (!scplx->point2tet) return;

    /* Every tet's vertices must have a non-NULL entry that contains that vertex. */
    for (const s_ncell *c = scplx->head; c; c = c->next) {
        for (int k = 0; k < 4; k++) {
            int v = c->vertex_id[k];
            if (!scplx->point2tet[v]) {
                fprintf(stderr, "assert_point2tet: vertex %d in tet %p has NULL point2tet\n",
                        v, (const void *)c);
                assert(scplx->point2tet[v]);
            }
        }
    }

    /* Every non-NULL entry must point to a tet that contains that vertex. */
    for (int v = 0; v < scplx->points.N; v++) {
        s_ncell *t = scplx->point2tet[v];
        if (!t) continue;
        int found = 0;
        for (int k = 0; k < 4; k++) found |= (t->vertex_id[k] == v);
        if (!found) {
            fprintf(stderr, "assert_point2tet: point2tet[%d] = %p does not contain vertex %d"
                    " (has %d %d %d %d)\n",
                    v, (void *)t,
                    v, t->vertex_id[0], t->vertex_id[1], t->vertex_id[2], t->vertex_id[3]);
            assert(found);
        }
    }
}


void print_ncell(const s_ncell *ncell)
{
    printf("%p, ( ", (void*)ncell);
    for (int ii=0; ii<4; ii++) {
        printf("%d ", ncell->vertex_id[ii]);
    }
    printf(") \n");
}


void print_scomplex(const s_scplx *scplx)
{
    puts("NCELLS");
    s_ncell *current = scplx->head;
    int ii = 0;
    while (current) {
        printf("%d  |  p : %p  |  prev : %p  |  marked : %d  |  vertex_ids :", ii, (void*)current, (void*)current->prev, current->mark_token);
        for (int jj=0; jj<4; jj++) {
            printf(" %d", current->vertex_id[jj]);
        }
        printf("  |  opposite :");
        for (int jj=0; jj<4; jj++) {
            printf(" %p", (void*)current->opposite[jj]);
        }
        printf("\n");
        ii++;
        current = current->next;
    }
    puts("");
}


// void write_ncell3d_file(const s_scplx *scplx, const s_ncell *ncell, FILE *file)
// {
//     fprintf(file, "%f %f %f\n", scplx->points.p[ncell->vertex_id[0]].x,
//                                 scplx->points.p[ncell->vertex_id[0]].y,
//                                 scplx->points.p[ncell->vertex_id[0]].z);
//     fprintf(file, "%f %f %f\n", scplx->points.p[ncell->vertex_id[1]].x,
//                                 scplx->points.p[ncell->vertex_id[1]].y,
//                                 scplx->points.p[ncell->vertex_id[1]].z);
//     fprintf(file, "%f %f %f\n\n", scplx->points.p[ncell->vertex_id[2]].x,
//                                   scplx->points.p[ncell->vertex_id[2]].y,
//                                   scplx->points.p[ncell->vertex_id[2]].z);
//
//     fprintf(file, "%f %f %f\n", scplx->points.p[ncell->vertex_id[0]].x,
//                                 scplx->points.p[ncell->vertex_id[0]].y,
//                                 scplx->points.p[ncell->vertex_id[0]].z);
//     fprintf(file, "%f %f %f\n", scplx->points.p[ncell->vertex_id[1]].x,
//                                 scplx->points.p[ncell->vertex_id[1]].y,
//                                 scplx->points.p[ncell->vertex_id[1]].z);
//     fprintf(file, "%f %f %f\n\n", scplx->points.p[ncell->vertex_id[3]].x,
//                                   scplx->points.p[ncell->vertex_id[3]].y,
//                                   scplx->points.p[ncell->vertex_id[3]].z);
//
//     fprintf(file, "%f %f %f\n", scplx->points.p[ncell->vertex_id[3]].x,
//                                 scplx->points.p[ncell->vertex_id[3]].y,
//                                 scplx->points.p[ncell->vertex_id[3]].z);
//     fprintf(file, "%f %f %f\n", scplx->points.p[ncell->vertex_id[1]].x,
//                                 scplx->points.p[ncell->vertex_id[1]].y,
//                                 scplx->points.p[ncell->vertex_id[1]].z);
//     fprintf(file, "%f %f %f\n\n", scplx->points.p[ncell->vertex_id[2]].x,
//                                   scplx->points.p[ncell->vertex_id[2]].y,
//                                   scplx->points.p[ncell->vertex_id[2]].z);
//
//     fprintf(file, "%f %f %f\n", scplx->points.p[ncell->vertex_id[0]].x,
//                                 scplx->points.p[ncell->vertex_id[0]].y,
//                                 scplx->points.p[ncell->vertex_id[0]].z);
//     fprintf(file, "%f %f %f\n", scplx->points.p[ncell->vertex_id[3]].x,
//                                 scplx->points.p[ncell->vertex_id[3]].y,
//                                 scplx->points.p[ncell->vertex_id[3]].z);
//     fprintf(file, "%f %f %f\n\n\n", scplx->points.p[ncell->vertex_id[2]].x,
//                                     scplx->points.p[ncell->vertex_id[2]].y,
//                                     scplx->points.p[ncell->vertex_id[2]].z);
// }
// void write_scomplex_file(const s_scplx *scplx, FILE *file)
// {
//     s_ncell *current = scplx->head;
//     while (current) {
//         write_ncell3d_file(scplx, current, file);
//         fprintf(file, "\n\n");
//         current = current->next;
//     }
// }



void extract_vertices_ncell(const s_scplx *scplx, const s_ncell *ncell, s_point out[4])
{
    for (int ii=0; ii<4; ii++)
        out[ii] = scplx->points.p[ncell->vertex_id[ii]];
}


void extract_weights_ncell(const s_scplx *scplx, const s_ncell *ncell, double out[4])
{
    if (!scplx->weights) { out[0]=0; out[1]=0; out[2]=0; out[3]=0; return; }
    for (int ii=0; ii<4; ii++)
        out[ii] = scplx->weights[ncell->vertex_id[ii]];
}



void extract_vertices_and_weights_ncell(const s_scplx *scplx, const s_ncell *ncell, s_point p_out[4], double w_out[4])
{
    for (int ii=0; ii<4; ii++) {
        int vid = ncell->vertex_id[ii];

        p_out[ii] = scplx->points.p[vid];
        w_out[ii] = scplx->weights ? scplx->weights[vid] : 0.0;
    }
}


void extract_ids_face(const s_ncell *ncell, int dim_face, const int *v_localid, int *out)
{   // v_localid[dim-dim_face], out[dim_face+1]
    if (dim_face < 0 || dim_face >= 3) {
        fprintf(stderr, "scplx.c: dim_face must be >= 0 && < 3\n");
        exit(1);
    }

    int kk = 0;
    for (int ii=0; ii<4; ii++) {
        if (!inarray(v_localid, 3 - dim_face, ii)) {
            out[kk] = ncell->vertex_id[ii];
            kk++;
        }
    }
}


void extract_vertices_face(const s_scplx *scplx, const s_ncell *ncell, int dim_face, 
                           const int v_localid[3-dim_face], s_point out[dim_face+1])
{   // v_localid[dim-dim_face], out[dim_face+1]
    if (dim_face < 0 || dim_face >= 3) {
        fprintf(stderr, "scplx.c: dim_face must be >= 0 && < 3\n");
        exit(1);
    }

    int face_vertex_id[dim_face+1];
    extract_ids_face(ncell, dim_face, v_localid, face_vertex_id);
    for (int ii=0; ii<dim_face+1; ii++) {
        out[ii] = scplx->points.p[face_vertex_id[ii]];
    }
}


void extract_weights_face(const s_scplx *scplx, const s_ncell *ncell, int dim_face,
                          const int v_localid[3-dim_face], double out[dim_face+1])
{
    if (dim_face < 0 || dim_face >= 3) {
        fprintf(stderr, "scplx.c: dim_face must be >= 0 && < 3\n");
        exit(1);
    }

    if (!scplx->weights) {
        for (int i=0; i<dim_face+1; i++) out[i] = 0;
        return;
    }

    int ids[dim_face+1];
    extract_ids_face(ncell, dim_face, v_localid, ids);
    for (int i=0; i<dim_face+1; i++)
        out[i] = scplx->weights ? scplx->weights[ids[i]] : 0.0;
}


void extract_vertices_and_weights_face(const s_scplx *scplx, const s_ncell *ncell, int dim_face,
                                       const int v_localid[3-dim_face],
                                       s_point pts_out[dim_face+1], double wts_out[dim_face+1])
{
    if (dim_face < 0 || dim_face >= 3) {
        fprintf(stderr, "scplx.c: dim_face must be >= 0 && < 3\n");
        exit(1);
    }

    int ids[dim_face+1];
    extract_ids_face(ncell, dim_face, v_localid, ids);
    for (int i = 0; i < dim_face+1; i++) {
        pts_out[i] = scplx->points.p[ids[i]];
        wts_out[i] = scplx->weights ? scplx->weights[ids[i]] : 0.0;
    }
}


void face_localid_of_adjacent_ncell(const s_ncell *ncell, int dim_face, const int v_localid[3-dim_face], int id_adjacent, int out_v_localid[3-dim_face])
{
    if (dim_face < 0 || dim_face >= 3) {
        fprintf(stderr, "dim_face must be >= 0 && < 3\n");
        exit(1);
    }

    s_ncell *adjacent = ncell->opposite[id_adjacent];

    int vertex_id_face[dim_face+1];
    extract_ids_face(ncell, dim_face, v_localid, vertex_id_face);

    int kk = 0;
    for (int ii=0; ii<4; ii++) {
        if (!inarray(vertex_id_face, dim_face+1, adjacent->vertex_id[ii])) {
            out_v_localid[kk] = ii;
            kk++;
        }
    }
}


s_ncell *next_ncell_ridge_cycle(const s_ncell *ncell, int v_localid_main, int v_localid_2,
                                int *new_v_localid_main, int *new_v_localid_2)
{
    s_ncell *next = ncell->opposite[v_localid_main];
    if (!next) return NULL;
    *new_v_localid_main = id_where_equal_int(next->vertex_id, 4, ncell->vertex_id[v_localid_2]);
    for (int ii=0; ii<4; ii++) {
        if (next->opposite[ii] == ncell) {
            *new_v_localid_2 = ii;
            return next;
        }
    }
    return NULL;
}


void walk_ridge_cycle_and_check_ncells(const s_ncell *start, const int omitted[2],
                                       bool (*check)(void *ctx, const s_ncell *), void *ctx)
{
    const s_ncell *cur = start;
    int cur_main = omitted[0], cur_2 = omitted[1];
    bool hit_boundary = false;
    do {
        if (!check(ctx, cur)) return;
        if (!cur->opposite[cur_main]) { hit_boundary = true; break; }
        int next_main, next_2;
        const s_ncell *next = next_ncell_ridge_cycle(cur, cur_main, cur_2,
                                                     &next_main, &next_2);
        if (next == start || !next) break;
        cur = next; cur_main = next_main; cur_2 = next_2;
    } while (cur != start);

    if (hit_boundary) {
        cur = start; cur_main = omitted[1]; cur_2 = omitted[0];
        if (!cur->opposite[cur_main]) return;  // boundary on both sides
        int next_main, next_2;
        cur = next_ncell_ridge_cycle(cur, cur_main, cur_2, &next_main, &next_2);
        cur_main = next_main; cur_2 = next_2;
        while (cur && cur != start) {
            if (!check(ctx, cur)) return;
            if (!cur->opposite[cur_main]) break;
            const s_ncell *next = next_ncell_ridge_cycle(cur, cur_main, cur_2,
                                                          &next_main, &next_2);
            if (next == start || !next) break;
            cur = next; cur_main = next_main; cur_2 = next_2;
        }
    }
}



static int new_mark_stamp(s_scplx *scplx)
{
    scplx->mark_stamp++;
    return scplx->mark_stamp;
}

static int mark_ncells_incident_face_STEP(int mark_stamp, s_ncell *ncell, int dim_face, const int v_localid[3-dim_face], s_dynarray *out_ncell_ptrs)
{
    for (int ii=0; ii<4; ii++) if (inarray(v_localid, 3-dim_face, ii)) {
        s_ncell *adjacent_ncell = ncell->opposite[ii];
        if (adjacent_ncell && adjacent_ncell->mark_token2 != mark_stamp) {
            adjacent_ncell->mark_token2 = mark_stamp;
            if (!dynarray_push(out_ncell_ptrs, &adjacent_ncell)) return 0;
            
            int new_v_localid[3-dim_face];
            face_localid_of_adjacent_ncell(ncell, dim_face, v_localid, ii, new_v_localid);

            /* Recursion: */
            if (!mark_ncells_incident_face_STEP(mark_stamp, adjacent_ncell, dim_face, new_v_localid, out_ncell_ptrs)) return 0;
        }
    }
    return 1;
}

int ncells_incident_face(s_scplx *scplx, s_ncell *ncell, int dim_face, const int v_localid[3-dim_face], s_dynarray *out_ncell_ptrs)
{   /* 0 ERROR, 1 OK */
    if (dim_face < 0 || dim_face >= 3) {
        fprintf(stderr, "dim_face must be >= 0 && < 3\n");
        exit(1);
    }

    int mark_stamp = new_mark_stamp(scplx);
    ncell->mark_token2 = mark_stamp;

    out_ncell_ptrs->N = 0;
    if (!dynarray_push(out_ncell_ptrs, &ncell)) return 0;
    
    /* Recursion, Flood-Fill */
    if (!mark_ncells_incident_face_STEP(mark_stamp, ncell, dim_face, v_localid, out_ncell_ptrs)) return 0;
    return 1;
}


int vertex_neighbors(s_scplx *scplx, int v_global,
                     s_ncell *start_cell, int v_local,
                     int skip_below, s_dynarray *out_ids,
                     s_dynarray *scratch_cells)
{
    /* The 3 local vertex slots NOT equal to v_local are the "opposite" face of v_global.
     * ncells_incident_face with dim_face=0 traverses all cells sharing those 3 slots,
     * which is exactly the star (all incident cells) of v_global. */
    int opp_localids[3], k = 0;
    for (int j = 0; j < 4; j++) if (j != v_local) opp_localids[k++] = j;

    scratch_cells->N = 0;
    if (!ncells_incident_face(scplx, start_cell, 0, opp_localids, scratch_cells))
        return -1;

    /* Dedup by scanning the output list itself: a vertex star is a few tens
     * of entries, so O(k) per candidate beats clearing any O(points.N) seen
     * structure per call (callers loop this over many vertices).  Anything
     * already in out_ids (callers may pre-seed it) is also never repeated. */
    int count = 0;
    s_ncell **cells = (s_ncell **)scratch_cells->items;
    for (unsigned ci = 0; ci < scratch_cells->N; ci++)
        for (int vi = 0; vi < 4; vi++) {
            int gid = cells[ci]->vertex_id[vi];
            if (gid == v_global || gid < skip_below) continue;
            const int *ids = (const int *)out_ids->items;
            bool dup = false;
            for (unsigned q = 0; q < out_ids->N && !dup; q++)
                dup = (ids[q] == gid);
            if (dup) continue;
            if (!dynarray_push(out_ids, &gid)) return -1;
            count++;
        }
    return count;
}



int test_point_in_ncell(const s_scplx *scplx, const s_ncell *ncell, s_point query)
{   // TODO Assumes consistent ordering of vertices?
    s_point v0 = scplx->points.p[ncell->vertex_id[0]];
    s_point v1 = scplx->points.p[ncell->vertex_id[1]];
    s_point v2 = scplx->points.p[ncell->vertex_id[2]];
    s_point v3 = scplx->points.p[ncell->vertex_id[3]];
    s_point tetra[4] = {v0, v1, v2, v3};
    return test_point_in_tetrahedron(tetra, query, 0, 0);
 }


s_ncell *bruteforce_find_ncell_containing(const s_scplx *scplx, s_point p)
{
    s_ncell *current = scplx->head;
    while (current) {
        e_geom_test test = test_point_in_ncell(scplx, current, p);
        if (test == TEST_IN || test == TEST_BOUNDARY) return current;
        if (test == TEST_ERROR) fprintf(stderr, "Warning: test error\n");
        // if (test == TEST_DEGENERATE);  /* These are expected! Just skip them */
        current = current->next;
    }
    /* No tet contains p: for a convex complex this means p is outside it. Return
     * NULL (a legitimate "outside" answer) rather than aborting the process. */
    return NULL;
}


static s_point ncell_centroid(const s_scplx *scplx, const s_ncell *ncell) {
    s_point vertices[4];
    extract_vertices_ncell(scplx, ncell, vertices);
    s_point out = { .x = (vertices[0].x + vertices[1].x + vertices[2].x + vertices[3].x) / 4.0,
                    .y = (vertices[0].y + vertices[1].y + vertices[2].y + vertices[3].y) / 4.0, 
                    .z = (vertices[0].z + vertices[1].z + vertices[2].z + vertices[3].z) / 4.0 };
    return out;
}


/* Fisher-Yates over {0,1,2,3} using the complex's caller-owned PRNG. When rng is
 * NULL the identity order {0,1,2,3} is used: still deterministic (the walk's
 * loop-detection + brute-force backstop handle any degenerate cycle). */
static void random_order_04(s_random_context *rng, int out[4])
{
    for (int ii=0; ii<4; ii++) out[ii] = ii;

    if (!rng) return;
    for (int ii=3; ii>0; ii--) {
        int jj = (int)random_uniform_range_u64(rng, (uint64_t)(ii + 1));
        int tmp = out[ii];
        out[ii] = out[jj];
        out[jj] = tmp;
    }
}

/* Shared point-location walk.  If query_id >= 0 the query is the registered
 * point points.p[query_id] and its orientation against each face is computed
 * through the predicate seam BY ID (exact-mode capable; used by
 * insert_one_point).  If query_id < 0 the query is the arbitrary coordinate p
 * and uses coordinate orientation (Voronoi callers, always exact_ids==0).  The
 * two are bit-identical when scplx->exact_ids == 0.  The o1 (opposite-vertex)
 * test is always a vertex, so it always goes through the seam by id. */
static s_ncell *in_ncell_walk_core(const s_scplx *scplx, s_point p, int query_id)
{
    // return bruteforce_find_ncell_containing(scplx, p);

    assert(scplx->N_ncells >= 1 && "N_ncells < 1");
    s_ncell *current = scplx->head;
    int hint = scplx->walk_hint_vid;
    /* rng is caller-owned; mutating *rng through a const complex is legal (the
     * pointee is non-const). NULL => start at index 0 / head (deterministic). */
    s_random_context *rng = scplx->rng;
    if (scplx->point2tet && hint >= 0 && hint < scplx->points.N &&
        scplx->point2tet[hint]) {
        current = scplx->point2tet[hint];
    } else if (scplx->point2tet) {
        int randi = rng ? (int)random_uniform_range_u64(rng, (uint64_t)scplx->points.N) : 0;
        for (int ii = 0; ii < scplx->points.N; ii++) {
            int idx = (randi + ii) % scplx->points.N;
            if (scplx->point2tet[idx]) { current = scplx->point2tet[idx]; break; }
        }
    } else {
        int randi = rng ? (int)random_uniform_range_u64(rng, (uint64_t)scplx->N_ncells) : 0;
        for (int ii = 0; ii < randi; ii++) current = current->next;
    }

    s_point facet_vertices[3];
    s_ncell *prev = current;

    /* Loop detection: stamp each cell as we enter it; if we re-enter a cell
     * already visited THIS walk, the walk is oscillating (degenerate/coplanar
     * geometry, or p on a boundary face) -- stop and resolve by brute force.
     * walk_stamp/walk_token are private to the walk (see scplx.h), so this never
     * disturbs mark_token/mark_stamp used by CDT and other traversals. The stamp
     * is mutable scratch on a logically-const complex, hence the cast. */
    const int wstamp = ++((s_scplx *)scplx)->walk_stamp;

    STEP:
    if (current->walk_token == wstamp)
        return bruteforce_find_ncell_containing(scplx, p);
    current->walk_token = wstamp;

    int order[4]; random_order_04(rng, order);
    for (int kk=0; kk<4; kk++) {  /* Visit faces in random order to prevent loops ? */
        int ii = order[kk];

        s_ncell *next = current->opposite[ii];  /* NULL == boundary (hull) face */

        int fids[3]; extract_ids_face(current, 2, &ii, fids);
        extract_vertices_face(scplx, current, 2, &ii, facet_vertices);

        int o1 = dtp_orient(scplx, fids[0], fids[1], fids[2], current->vertex_id[ii]);
        int o2 = (query_id >= 0)
                     ? dtp_orient(scplx, fids[0], fids[1], fids[2], query_id)
                     : test_orientation(facet_vertices, p);

        if (o1 == 0) {  /* Tetrahedron is degenerate */
            if (next && next != prev) {  /* come from a different adjacent one -> walk towards */
                prev = current;
                current = next;
                goto STEP;
            }
            else continue;
        }

        if (o2 == 0) {  /* Query is coplanar with this face's plane */
            if (scplx->exact_ids && query_id >= 0) {
                /* Exact walk (reference-style, CDT insertExistingVertex): if the
                 * query lies IN this face it is on the tet boundary -> found;
                 * otherwise it is coplanar-but-outside, so it is strictly beyond
                 * some OTHER face -- do NOT step here (no coordinate centroid
                 * tie-break, which drifts/loops), let another face step. */
                if (dtp_point_in_triangle(scplx, query_id, fids[0], fids[1], fids[2]) >= 0)
                    return current;
                continue;
            }
            e_geom_test test = test_point_in_triangle_3D(facet_vertices, p, 0, 0);
            if (test == TEST_IN || test == TEST_BOUNDARY) { return current; }
            if (next && next != prev) {
                // TIE-BREAKING: move to the neighbor whose centroid is closer to p.
                double d_curr = distance_squared(ncell_centroid(scplx, current), p);
                double d_next = distance_squared(ncell_centroid(scplx, next), p);
                if (d_next < d_curr) {
                    prev = current;
                    current = next;
                    goto STEP;
                }
            }
            continue;
        } else if (o1 != o2) {  // Query strictly beyond this face -> step across it
            if (!next) return NULL;   /* beyond a boundary face -> outside the domain
                                       * (exact only for a convex complex; walking is
                                       * only sound there anyway) */
            prev = current;
            current = next;
            goto STEP;
        }
    }

    return current;
}

/* Wrapper: run the walk, then remember where it ended (by VERTEX id, never by
 * tet pointer -- flips free tets, but point2tet is maintained through every
 * flip) so the next walk starts there.  Consecutive spatially coherent queries
 * (BRIO-style insertion, per-vertex owner scans) then walk O(1) tets instead
 * of O(n^1/3) from a random start.  The stamp is mutable scratch on a
 * logically-const complex, same pattern as walk_stamp above. */
static s_ncell *in_ncell_walk_impl(const s_scplx *scplx, s_point p, int query_id)
{
    s_ncell *res = in_ncell_walk_core(scplx, p, query_id);
    if (res)
        ((s_scplx *)scplx)->walk_hint_vid = res->vertex_id[0];
    return res;
}

s_ncell *in_ncell_walk(const s_scplx *scplx, s_point p)
{
    return in_ncell_walk_impl(scplx, p, -1);
}

s_ncell *in_ncell_walk_id(const s_scplx *scplx, int point_id)
{
    return in_ncell_walk_impl(scplx, scplx->points.p[point_id], point_id);
}


double power_distance_point_vertex(const s_scplx *scplx, int vid, s_point p)
{
    if (!scplx->weights) return distance_squared(scplx->points.p[vid], p);
    else return distance_squared(scplx->points.p[vid], p) - scplx->weights[vid];
}




/* PLOTS */
void plot_add_ncell(s_gnuplot *interface, const s_scplx *scplx, const s_ncell *ncell, char *config)
{
    s_point face_vertices[3];
    for (int ii=0; ii<4; ii++) {
        extract_vertices_face(scplx, ncell, 2, &ii, face_vertices);
        draw_solid_triangle_3d(interface, face_vertices[0].coords, face_vertices[1].coords,
                               face_vertices[2].coords, config);
    }
}


void plot_ncell_3d(const s_scplx *scplx, const s_ncell *ncell, char *f_name, s_point ranges[2])
{
    s_gnuplot *interface = gnuplot_start(PNG_3D, f_name, (int[2]){1080, 1080}, 18);
    gnuplot_config(interface, "set pm3d depthorder",
                              "set pm3d border lc 'black' lw 0.5",
                              "set view 100, 60",
                              "set xyplane at 0");
    if (ranges) {
        char buff[1024];
        snprintf(buff, 1024, "set xrange [%f:%f]\n set yrange [%f:%f]\n set zrange [%f:%f]", 
                 ranges[0].x, ranges[1].x, ranges[0].y, ranges[1].y, ranges[0].z, ranges[1].z); 
        gnuplot_config(interface, buff);
    }

    plot_add_ncell(interface, scplx, ncell, "fs transparent solid 0.2 fc rgb '#000090'");
    gnuplot_end(interface);
}


void plot_all_ncells_3d(const s_scplx *scplx, char *f_name, s_point ranges[2], char *view_command)
{
     char colors[][20] = { "#000090", "#000fff", "#0090ff", "#0fffee", 
        "#90ff70", "#ffee00", "#ff7000", "#ee0000", "#7f0000" };
    char buff[1024];

    s_gnuplot *interface = gnuplot_start(PNG_3D, f_name, (int[2]){1080, 1080}, 18);
    gnuplot_config(interface, "set pm3d depthorder",
                              "set pm3d border lc 'black' lw 0.5",
                              "set xyplane at 0",
                              view_command);
    if (ranges) {
        snprintf(buff, 1024, "set xrange [%f:%f]\n set yrange [%f:%f]\n set zrange [%f:%f]", 
                 ranges[0].x, ranges[1].x, ranges[0].y, ranges[1].y, ranges[0].z, ranges[1].z); 
    }
    if (ranges) gnuplot_config(interface, buff);

    int it = 0;
    s_ncell *current = scplx->head;
    while (current) {
        snprintf(buff, 1024, "fs transparent solid 0.2 fc rgb '%s'", colors[it%8]);
        plot_add_ncell(interface, scplx, current, buff);
        current = current->next;
        it++;
    }
    gnuplot_end(interface);
}


void plot_dt_3d_differentviews(const s_scplx *scplx, char *f_name, s_point ranges[2])
{
    char final_name[512];
    snprintf(final_name, 512, "%s_v1.png", f_name);
    plot_all_ncells_3d(scplx, final_name, ranges, "set view 100, 60, 1.5");

    snprintf(final_name, 512, "%s_v2.png", f_name);
    plot_all_ncells_3d(scplx, final_name, ranges, "set view 100, 90, 1.5");

    snprintf(final_name, 512, "%s_v3.png", f_name);
    plot_all_ncells_3d(scplx, final_name, ranges, "set view 100, 180, 1.5");

    snprintf(final_name, 512, "%s_v4.png", f_name);
    plot_all_ncells_3d(scplx, final_name, ranges, "set view 100, 270, 1.5");

}

