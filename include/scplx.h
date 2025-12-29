#ifndef VOR3D_SIMPLICAL_COMPLEX_H
#define VOR3D_SIMPLICAL_COMPLEX_H

#include <stdio.h>
#include "points.h"
#include "lists.h"

typedef struct simplical_complex {  // May live in stack
    s_points points;  // N = (4 + n_points) 
                      // The first 4 points correspond to the big n_cell
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
    // int count;   // ID
} s_ncell;


// ncells are dynamically allocated, complex is static
s_ncell *malloc_ncell();
void free_ncell(s_ncell *ncell);
void free_complex(s_scplx *setup);

void print_ncell(const s_ncell *ncell);
void print_scomplex(const s_scplx *setup);

int ncells_incident_face(s_scplx *setup, s_ncell *ncell, int dim_face, const int *v_localid, s_list *out);

void extract_vertices_ncell(const s_scplx *setup, const s_ncell *ncell, s_point out[4]);
void extract_ids_face(const s_ncell *ncell, int dim_face, const int *v_localid, int *out);
void extract_vertices_face(const s_scplx *setup, const s_ncell *ncell, int dim_face, const int v_localid[3-dim_face], s_point out[dim_face+1]);
void face_localid_of_adjacent_ncell(const s_ncell *ncell, int dim_face, const int v_localid[3-dim_face], int id_adjacent, int out_v_localid[3-dim_face]);
s_ncell *next_ncell_ridge_cycle(const s_ncell *ncell, int v_localid_main, int v_localid_2, int *new_v_localid_main, int *new_v_localid_2);

int are_locally_delaunay_strict(const s_scplx *setup, const s_ncell *ncell, int id_opposite);
int are_locally_delaunay_nonstrict(const s_scplx *setup, const s_ncell *ncell, int id_opposite);
int test_point_in_ncell(const s_scplx *setup, const s_ncell *ncell, s_point query);
s_ncell *bruteforce_find_ncell_containing(const s_scplx *setup, s_point p);
s_ncell *in_ncell_walk(const s_scplx *setup, s_point p);
int is_delaunay_3d(const s_scplx *setup);

void plot_dt_3d_differentviews(const s_scplx *setup, char *f_name, s_point ranges[2]);

#endif
