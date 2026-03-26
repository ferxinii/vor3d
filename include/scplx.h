#ifndef VOR3D_SIMPLICAL_COMPLEX_H
#define VOR3D_SIMPLICAL_COMPLEX_H

#include <stdio.h>
#include "points.h"

typedef struct simplical_complex {  // May in stack
    s_points points;  
    double *weights;  // size point.N, May be NULL if unweighted
    int N_ncells;
    struct ncell *head;  // Linked list of ncells
    int mark_stamp;
} s_scplx;


typedef struct ncell {  // Must live in heap
    int vertex_id[4];
    struct ncell *opposite[4];
    struct ncell *next;  // Linked list of cells
    struct ncell *prev;
    int mark_token;  // Used to mark particular ncells: mark_token == mark_stamp
    bool mask_alpha;  // if it belongs to the alpha complex (for a given alpha)
} s_ncell;

typedef enum delaunay_test_type {
    DELAUNAY_TEST_STRICT,
    DELAUNAY_TEST_NONSTRICT
} e_delaunay_test_type;

typedef struct dynarray s_dynarray;
typedef struct hash_table s_hash_table;

void free_complex(s_scplx *scplx);
void print_ncell(const s_ncell *ncell);
void print_scomplex(const s_scplx *scplx);


int are_locally_delaunay(const s_scplx *scplx, const s_ncell *ncell, int id_opposite, 
                         e_delaunay_test_type type);
int is_delaunay_3d(const s_scplx *scplx, e_delaunay_test_type type);
int test_point_in_ncell(const s_scplx *scplx, const s_ncell *ncell, s_point query);
s_ncell *bruteforce_find_ncell_containing(const s_scplx *scplx, s_point p);
s_ncell *in_ncell_walk(const s_scplx *scplx, s_point p);
void plot_dt_3d_differentviews(const s_scplx *scplx, char *f_name, s_point ranges[2]);


/* HELPERS */
s_ncell *malloc_ncell();
void free_ncell(s_ncell *ncell);

int ncells_incident_face(s_scplx *scplx, s_ncell *ncell, int dim_face,
                         const int *v_localid, s_dynarray *out);
void face_localid_of_adjacent_ncell(const s_ncell *ncell, int dim_face,
                                    const int v_localid[3-dim_face], int id_adjacent,
                                    int out_v_localid[3-dim_face]);
s_ncell *next_ncell_ridge_cycle(const s_ncell *ncell, int v_localid_main, int v_localid_2, 
                                int *new_v_localid_main, int *new_v_localid_2);
void walk_ridge_cycle_and_check_ncells(const s_ncell *start, const int omitted[2],
                                       bool (*check)(void *ctx, const s_ncell *), void *ctx);

double power_distance_point_vertex(const s_scplx *scplx, int vid, s_point p);

void extract_vertices_ncell(const s_scplx *scplx, const s_ncell *ncell, s_point out[4]);
void extract_weights_ncell(const s_scplx *scplx, const s_ncell *ncell, double out[4]);
void extract_vertices_and_weights_ncell(const s_scplx *scplx, const s_ncell *ncell,
                                        s_point p_out[4], double w_out[4]);
void extract_ids_face(const s_ncell *ncell, int dim_face, const int *v_localid, int *out);
void extract_vertices_face(const s_scplx *scplx, const s_ncell *ncell, int dim_face,
                           const int v_localid[3-dim_face], s_point out[dim_face+1]);
void extract_weights_face(const s_scplx *scplx, const s_ncell *ncell, int dim_face,
                          const int v_localid[3-dim_face], double out[dim_face+1]);
void extract_vertices_and_weights_face(const s_scplx *scplx, const s_ncell *ncell, 
                                       int dim_face, const int v_localid[3-dim_face],
                                       s_point p_out[dim_face+1], double w_out[dim_face+1]);


#endif
