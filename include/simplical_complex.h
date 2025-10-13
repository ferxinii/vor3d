#ifndef VOR3D_SIMPLICAL_COMPLEX_H
#define VOR3D_SIMPLICAL_COMPLEX_H

#include <stdio.h>
#include "external/geometry/geometry.h"

typedef struct simplical_complex {
    int dim;
    int N_points;
    s_point *points;  // (n_points + dim + 1) x dim
                      // The first? (dim+1) points correspond to the big n_cell
    int N_ncells;
    struct ncell *head;  // Linked list of ncells
} s_scplx;


typedef struct ncell {
    int *vertex_id;
    struct ncell **opposite;
    struct ncell *next;  // Linked list of cells
    struct ncell *prev;
    int mark;  // Used to mark particular ncells
    int count;  // ID
    double volume;  // DEBUGGING
} s_ncell;


s_ncell *malloc_ncell(const s_scplx *setup);
void free_ncell(s_ncell *ncell);
void free_complex(s_scplx *setup);

void print_ncell(const s_scplx *setup, const s_ncell *ncell);
void print_ncells(const s_scplx *setup);
void write_ncell3d_file(s_scplx *setup, s_ncell *ncell, FILE *file);
void write_dt3d_file(s_scplx *setup, FILE *file);

void initialize_ncells_counter(const s_scplx *setup);
void initialize_ncells_mark(const s_scplx *setup);
void print_marked(const s_scplx *setup);
int count_marked(const s_scplx *setup);
void mark_ncells_incident_face(const s_scplx *setup, s_ncell *ncell, const int *v_localid, int dim_face);

void extract_vertices_ncell(const s_scplx *setup, const s_ncell *ncell, s_point *out);
void extract_ids_face(const s_scplx *setup, const s_ncell *ncell, const int *v_localid, int dim_face, int *out);
void extract_vertices_face(const s_scplx *setup, const s_ncell *ncell, const int *v_localid, int dim_face, s_point *out);
void extract_ace_center_and_normal(const s_scplx *setup, const s_ncell *ncell, int face_localid, s_point *fc, s_point *n);
void face_localid_of_adjacent_ncell(const s_scplx *setup, const s_ncell *ncell, const int *v_localid,
                                    int dim_face, int id_adjacent, int *out_v_localid);
s_ncell *next_ncell_ridge_cycle(const s_scplx *setup, const s_ncell *ncell, int v_localid_main, int v_localid_2, 
                                int *new_v_localid_main, int *new_v_localid_2);
int count_cycle_ridge(const s_scplx *setup, const s_ncell *ncell, int v_localid_main, int v_localid_2);

int are_locally_delaunay_strict(const s_scplx *setup, const s_ncell *ncell, int id_opposite);
int are_locally_delaunay_nonstrict(const s_scplx *setup, const s_ncell *ncell, int id_opposite);
int point_in_tetra(const s_scplx *setup, s_point x, const s_ncell *nc);  // Assumes consisten ordering of vertices?
s_ncell *in_ncell_walk(const s_scplx *setup, s_point p);
int is_delaunay_3d(const s_scplx *setup);

void add_ncell_volume_3d(const s_scplx *setup, s_ncell *ncell);
double compute_volume_complex(s_scplx *setup);

void plot_dt_3d_differentviews(const s_scplx *setup, char *f_name, s_point ranges[2]);

#endif
