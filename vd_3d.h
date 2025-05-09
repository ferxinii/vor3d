#ifndef VD_3D_H
#define VD_3D_H

#include <float.h>
#include "simplical_complex.h"
// #define CONVHULL_3D_ENABLE
// #include "convhull_3d.h"
#include "bpoly.h"

#define VCELL_BLOCK_VERTICES 1000
#define MAX_PLANES 1000

typedef struct plane {
    double A[3];
    double b;
} s_plane;

typedef struct planes_poly {
    s_plane planes[MAX_PLANES];
    int id[MAX_PLANES];
    int Nplanes;
} s_planes_poly;


typedef struct vdiagram {
    int N_vcells;
    struct vcell **vcells;  // Array of pointers to the cells, N_vcells x 1
    const struct bound_poly *bpoly;
    double **seeds;
} s_vdiagram;

typedef struct vcell {
    int seed_id;
    // struct vcell *next;     // linked list of vcells
    int Nv;
    int Nv_capacity;
    double **vertices;  // Nv x 3
    int **origin_vertices;  // Nv x 4, LAST column indicates the dual ncell if POSITIVE,
                            // if -1: comes from circumcenter delaunay. The rest of 
                            //        the columns indicate the delaunay indices of the 
                            //        face whose normal was extended, 
                            // if -2: ARTIFICIAL extension, first column is id of face of
                            //        convhull of setup points
                            // if -3: it is from the bounding polyhedron, and the first 
                            //        index is the point id 
                            // if -4: it comes from some coords, supposedly from the ray 
                            //        intersection with bounding convex hull, and the first
                            //        columns correspond to face's vertex_id
    int Nf;
    int *faces;
    double **fnormals;
    double volume;
    s_planes_poly *planes_poly;
} s_vcell;


void free_vdiagram(s_vdiagram *vdiagram);

void write_vd_file(s_vdiagram *vd, FILE *file);

s_vdiagram *malloc_vdiagram(const s_setup *setup, int Nreal);

void print_vdiagram(const s_vdiagram *vdiagram);

s_vcell *malloc_vcell(int seed_id);

void compute_vcell_volume(s_vcell *vcell);

s_vdiagram *voronoi_from_delaunay_3d(const s_setup *setup, s_bound_poly *bpoly, int Nreal);

int find_inside_which_vcell(s_vdiagram *vd, double *x);

void plot_vcell(s_vdiagram *vdiag, s_vcell *vcell, char *f_name, double *ranges);

void plot_vdiagram(s_vdiagram *vdiagram, char *f_name, double *ranges, int max_files, double **aux_points, int *N_aux);

void plot_vdiagram_auto(s_vdiagram *vdiagram, char *f_name, int max_files);

void clear_volumes_file(char *fname);

void append_volumes_to_file(s_vdiagram *vdiagram, char *fname, int id);

#endif
