#ifndef SIMPLICAL_COMPLEX_H
#define SIMPLICAL_COMPLEX_H

#include <stdio.h>

typedef struct setup {
    int dim;
    int N_points;
    double **points;  // (n_points + dim + 1) x dim
                      // The last (dim+1) points correspond to the big n_cell
    int N_ncells;
    struct ncell *head;  // Linked list of ncells
} s_setup;


typedef struct ncell {
    int *vertex_id;
    struct ncell **opposite;
    struct ncell *next;  // Linked list of cells
    struct ncell *prev;
    int mark;  // Used to mark particular ncells
    int count;
    double volume;  // DEBUGGING
} s_ncell;


s_ncell *malloc_ncell(const s_setup *setup);
void free_ncell(s_ncell *ncell);
void free_complex(s_setup *setup);
void print_ncell(const s_setup *setup, const s_ncell *ncell);
void print_ncells(const s_setup *setup);
void write_ncell3d_file(s_setup *setup, s_ncell *ncell, FILE *file);
void write_dt3d_file(s_setup *setup, FILE *file);
void initialize_ncells_counter(const s_setup *setup);
void initialize_ncells_mark(const s_setup *setup);
void print_marked(const s_setup *setup);
int count_marked(const s_setup *setup);
void extract_vertices_ncell(const s_setup *setup, const s_ncell *ncell, double **out);
void extract_ids_face(const s_setup *setup, const s_ncell *ncell, const int *v_localid, int dim_face, int *out);
void extract_vertices_face(const s_setup *setup, const s_ncell *ncell, const int *v_localid, int dim_face, double **out);
void extract_face_center_and_normal(const s_setup *setup, const s_ncell *ncell, int face_localid, double *fc, double *n);
void face_localid_of_adjacent_ncell(const s_setup *setup, const s_ncell *ncell, const int *v_localid,
                                         int dim_face, int id_adjacent, int *out_v_localid);
s_ncell *next_ncell_ridge_cycle(const s_setup *setup, const s_ncell *ncell, int v_localid_main, int v_localid_2, 
                                int *new_v_localid_main, int *new_v_localid_2);
int count_cycle_ridge(const s_setup *setup, const s_ncell *ncell, int v_localid_main, int v_localid_2);
void mark_ncells_incident_face(const s_setup *setup, s_ncell *ncell, const int *v_localid, int dim_face);
int are_locally_delaunay_strict(const s_setup *setup, const s_ncell *ncell, int id_opposite);
int are_locally_delaunay_nonstrict(const s_setup *setup, const s_ncell *ncell, int id_opposite);
int point_in_tetra(s_setup *setup, double *x, s_ncell *nc);
s_ncell *in_ncell_walk(s_setup *setup, double *p);
int is_delaunay_3d(const s_setup *setup);
void add_ncell_volume_3d(s_setup *setup, s_ncell *ncell);
double compute_volume_complex(s_setup *setup);
void plot_dt_3d(s_setup *setup, char *f_name, double *ranges, int max_files);

#endif
