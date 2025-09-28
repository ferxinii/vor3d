#ifndef VOR3D_VDIAGRAM_H
#define VOR3D_VDIAGRAM_H

#include <float.h>
#include "simplical_complex.h"
#include "bpoly.h"


#define VCELL_BLOCK_VERTICES 1000


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
void plot_vdiagram_differentviews(s_vdiagram *vdiagram, char *f_name, double *ranges);

void clear_volumes_file(char *fname);
void append_volumes_to_file(s_vdiagram *vdiagram, char *fname, int id);

#endif
