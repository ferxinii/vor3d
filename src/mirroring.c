
#include "mirroring.h"
#include "bpoly.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>


static s_point mirror_plane(s_point normal, double d_plane, s_point p)
{
    double factor = 2 * (d_plane - dot_prod(normal, p));
    s_point out;
    out.x = p.x + factor * normal.x;
    out.y = p.y + factor * normal.y;
    out.z = p.z + factor * normal.z;
    return out;
}


/* Memory handling for lists */
s_vertex_list initialize_vertex_list(int Nmax) {
    if (Nmax <= 0) Nmax = 50;
    s_vertex_list out = { .Nmax = Nmax, 
                          .list = malloc(Nmax * sizeof(s_extruding_vertex)),
                          .N = 0 };
    return out;
}

int increase_memory_vertex_list_if_needed(s_vertex_list *v_list, int N_needed) {
    while (N_needed >= v_list->Nmax) {
        s_extruding_vertex *tmp = realloc(v_list->list, 2 * v_list->Nmax * sizeof(s_extruding_vertex));
        if (!tmp) {
            return 0;
        }
        v_list->list = tmp;
        v_list->Nmax *= 2;
    }
    return 1;
}

void free_vertex_list(s_vertex_list *v_list) {
    free(v_list->list);
    memset(v_list, 0, sizeof(s_vertex_list));
}


typedef struct point_list {
    s_point *list;
    int N;
    int Nmax;
} s_point_list;

static s_point_list initialize_point_list(int Nmax) {
    if (Nmax <= 0) Nmax = 50;
    s_point_list out = { .Nmax = Nmax, 
                         .list = malloc(Nmax * sizeof(s_point)),
                         .N = 0 };
    return out;
}

static int increase_memory_point_list_if_needed(s_point_list *point_list, int N_needed) {
    while (N_needed >= point_list->Nmax) {
        s_point *tmp = realloc(point_list->list, 2 * point_list->Nmax * sizeof(s_point));
        if (!tmp) return 0;
        point_list->list = tmp;
        point_list->Nmax *= 2;
    }
    return 1;
}

static void free_point_list(s_point_list *point_list) {
    free(point_list->list);
    memset(point_list, 0, sizeof(s_point_list));
}


typedef struct int_list {
    int *list;
    int N;
    int Nmax;
} s_int_list;

static s_int_list initialize_int_list(int Nmax) {
    if (Nmax <= 0) Nmax = 50;
    s_int_list out = { .Nmax = Nmax, 
                       .list = malloc(Nmax * sizeof(int)),
                       .N = 0 };
    return out;
}

static int increase_memory_int_list_if_needed(s_int_list *int_list, int N_needed) {
    while (N_needed >= int_list->Nmax) {
        int *tmp = realloc(int_list->list, 2 * int_list->Nmax * sizeof(int));
        if (!tmp) {
            return 0;
        }
        int_list->list = tmp;
        int_list->Nmax *= 2;
    }
    return 1;
}

static void free_int_list(s_int_list *int_list) {
    free(int_list->list);
    memset(int_list, 0, sizeof(s_int_list));
}



/* Initial mirroring */
static int witness_points_from_face(const s_points *seeds, const s_point face[3], double EPS_degenerate, s_point out_witness[seeds->N+3])
{   /* 0 ERROR, 1 OK */
    int N = 3;
    out_witness[0] = face[0];  out_witness[1] = face[1];  out_witness[2] = face[2];

    for (int ii=0; ii<seeds->N; ii++) {
        s_point closest = closest_point_on_triangle(face, EPS_degenerate, seeds->p[ii]);
        if (!point_is_valid(closest)) return 0;
        out_witness[N++] = closest;
    }
    return 1;
}

static int find_closest_seeds_to_witness(s_point witness, const s_points *seeds, double TOL2, s_int_list *buff_closest)
{   /* Find those seeds which are closest to the witness, might be more than 1 */
    /* 0 ERROR, 1 OK */
    double dmin2 = distance_squared(witness, seeds->p[0]);
    buff_closest->N = 1;
    buff_closest->list[0] = 0;
    for (int jj=1; jj<seeds->N; jj++) {
        double d2 = distance_squared(witness, seeds->p[jj]);
        if (fabs(d2-dmin2) < TOL2) {  /* Same distance than previous min */
            if (!increase_memory_int_list_if_needed(buff_closest, buff_closest->N+1)) return 0;
            buff_closest->list[buff_closest->N++] = jj;
            // printf("DEBUG: TIE!\n");
        } else if (d2 < dmin2) {  /* Found new minimum */
            dmin2 = d2;
            buff_closest->N = 1;
            buff_closest->list[0] = jj;
        }
    }
    return 1;
}

static int mirror_seeds_face(const s_points *seeds, const s_point face[3], s_point normal, double d_plane, double EPS_degenerate, double TOL2, bool buff_mask[seeds->N], s_point buff_witness[seeds->N+3], s_int_list *buff_closest, s_point_list *out_mirrored)
{   /* 0 ERROR, 1 OK */
    memset(buff_mask, 0, sizeof(bool) * seeds->N);
    
    if (!witness_points_from_face(seeds, face, EPS_degenerate, buff_witness)) return 0;
    /* Loop over witness points (N+3) */
    for (int ii=0; ii<seeds->N+3; ii++) {  
        if (!find_closest_seeds_to_witness(buff_witness[ii], seeds, TOL2, buff_closest)) return 0;
        /* mirror if not already mirrored */
        for (int jj=0; jj<buff_closest->N; jj++) {
            if (!buff_mask[buff_closest->list[jj]]) {
                if (!increase_memory_point_list_if_needed(out_mirrored, out_mirrored->N+1)) return 0;
                out_mirrored->list[out_mirrored->N++] = mirror_plane(normal, d_plane, seeds->p[buff_closest->list[jj]]);
                buff_mask[buff_closest->list[jj]] = true;
            }
        }
    }
    return 1;
}


int extend_sites_mirroring_initial(const s_bpoly *bp, double EPS_degenerate, double TOL, s_points *inout_seeds)
{   /* 0 ERROR, 1 OK */
    s_point_list mirrored = {0};
    s_int_list buff_closest = {0};
    s_point *buff_witness = NULL;
    bool *buff_mask_mirrored = {0};
    
    mirrored = initialize_point_list(0);
    buff_closest = initialize_int_list(10);
    buff_witness = malloc(sizeof(s_point) * (inout_seeds->N + 3)); 
    buff_mask_mirrored = calloc(inout_seeds->N, sizeof(bool));
    if (!mirrored.list || !buff_closest.list || !buff_witness || !buff_mask_mirrored) goto error;
    
    /* Compute all mirrored seeds */
    const double TOL2 = TOL * TOL;
    for (int ii=0; ii<bp->convh.Nf; ii++) {
        s_point face[3];
        convh_get_face(&bp->convh, ii, face);
        s_point normal = normalize_vec(bp->convh.fnormals[ii], EPS_degenerate);
        if (!point_is_valid(normal)) goto error;
        double d_plane =  dot_prod(normal, face[0]);
        if (!mirror_seeds_face(inout_seeds, face, normal, d_plane, EPS_degenerate, TOL2, buff_mask_mirrored, buff_witness, &buff_closest, &mirrored)) goto error;
    }
    printf("DEBUG: mirrored: %d\n", mirrored.N);

    /* Add them to inout_seeds */
    s_point *tmp = realloc(inout_seeds->p, sizeof(s_point) * (inout_seeds->N + mirrored.N));
    if (!tmp) goto error;
    for (int ii=0; ii< mirrored.N; ii++) 
        tmp[inout_seeds->N + ii] = mirrored.list[ii];

    inout_seeds->p = tmp;
    inout_seeds->N += mirrored.N;
    printf("DEBUG: mirrored->N = %d, total = %d\n", mirrored.N, inout_seeds->N);

    free_point_list(&mirrored);
    free_int_list(&buff_closest);
    free(buff_witness);
    free(buff_mask_mirrored);
    return 1;
    
    error:
        free_point_list(&mirrored);
        free_int_list(&buff_closest);
        free(buff_witness);
        free(buff_mask_mirrored);
        return 0;
}


static int mirror_from_extruding_vertex(const s_points *seeds, const s_bpoly *bp, double EPS_degenerate, double TOL2, s_point vertex, s_int_list *buff_closest_faces, s_int_list *buff_closest_seeds, s_point_list *mirrored)
{
    s_int_list *closest_faces = buff_closest_faces;

    /* Find faces which are closest to the extruding vertex */
    double dmin2 = distance_squared(vertex, bp->convh.points.p[bp->convh.faces[0]]);
    closest_faces->N = 1;
    closest_faces->list[0] = 0;
    for (int ii=0; ii<bp->convh.Nf; ii++) {
        s_point face[3];
        convh_get_face(&bp->convh, ii, face);

        s_point closest = closest_point_on_triangle(face, EPS_degenerate, vertex);
        double d2 = distance_squared(vertex, closest);

        if (fabs(d2-dmin2) < TOL2) {  /* Same distance than previous min */
            if (!increase_memory_int_list_if_needed(closest_faces, closest_faces->N+1)) return 0;
            closest_faces->list[closest_faces->N++] = ii;
        } else if (d2 < dmin2) {  /* Found new minimum */
            dmin2 = d2;
            closest_faces->N = 1;
            closest_faces->list[0] = ii;
        }
    }

    // printf("DEBUG: closest_faces->N = %d\n", closest_faces->N);
    for (int ii=0; ii<closest_faces->N; ii++) {
        int face_id = closest_faces->list[ii];
        s_point face[3];
        convh_get_face(&bp->convh, face_id, face);
        s_point normal = normalize_vec(bp->convh.fnormals[face_id], EPS_degenerate);
        if (!point_is_valid(normal)) return 0;
        double d_plane =  dot_prod(normal, face[0]);

        s_point witness = closest_point_on_triangle(face, EPS_degenerate, vertex);
        s_int_list *closest_seeds = buff_closest_seeds;
        find_closest_seeds_to_witness(witness, seeds, TOL2, closest_seeds);
        if (closest_seeds->N > 1) printf("DEBUG: MULTIPLE SEEDS! %d\n", closest_seeds->N);
        for (int ii=0; ii<closest_seeds->N; ii++) {
            if (!increase_memory_point_list_if_needed(mirrored, mirrored->N+1)) return 0;
            mirrored->list[mirrored->N++] = mirror_plane(normal, d_plane, seeds->p[closest_seeds->list[ii]]);
        }
    }
    return 1;
}


int extend_sites_mirroring_extruding(const s_bpoly *bp, double EPS_degenerate, double TOL, const s_vertex_list *extruding, int Nreal, s_points *inout_seeds)
{   /* 0 ERROR, 1 OK */
    s_int_list buff_closest_faces = {0}, buff_closest_seeds = {0};
    s_point_list mirrored = {0};

    buff_closest_faces = initialize_int_list(3);
    buff_closest_seeds = initialize_int_list(3);
    mirrored = initialize_point_list(extruding->N);
    if (!buff_closest_faces.list || !buff_closest_seeds.list || !mirrored.list) goto error;

    /* For each extruding vertex, mirror necessary seeds accross necessary faces */
    const int Naux = inout_seeds->N;
    inout_seeds->N = Nreal;
    const double TOL2 = TOL*TOL;
    for (int ii=0; ii<extruding->N; ii++) {
        s_point vertex = extruding->list[ii].vertex;
        if (!mirror_from_extruding_vertex(inout_seeds, bp, EPS_degenerate, TOL2, vertex, &buff_closest_faces, &buff_closest_seeds, &mirrored)) goto error;
    }

    /* Add mirrored to inout_seeds */
    inout_seeds->N = Naux;
    s_point *tmp = realloc(inout_seeds->p, sizeof(s_point) * (inout_seeds->N + mirrored.N));
    if (!tmp) goto error;
    for (int ii=0; ii< mirrored.N; ii++) 
        tmp[inout_seeds->N + ii] = mirrored.list[ii];

    inout_seeds->p = tmp;
    inout_seeds->N += mirrored.N;
    printf("DEBUG: EXTRUDING: mirrored->N = %d, total = %d\n", mirrored.N, inout_seeds->N);

    free_int_list(&buff_closest_faces);
    free_int_list(&buff_closest_seeds);
    free_point_list(&mirrored);
    return 1;

    error:
        free_int_list(&buff_closest_faces);
        free_int_list(&buff_closest_seeds);
        free_point_list(&mirrored);
        return 0;
}


