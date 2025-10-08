#ifndef VOR3D_VDIAGRAM_H
#define VOR3D_VDIAGRAM_H

#include <float.h>
#include "simplical_complex.h"
#include "bpoly.h"


#define VCELL_BLOCK_VERTICES 1000  // Used for efficient mallocing, incrementing by blocks if necessary


typedef struct vdiagram {
    int N;  // Number of vcells
    struct vcell **vcells;  // Array of pointers to the cells, N_vcells x 1
    const struct bounding_polyhedron *bpoly;
    s_point *seeds;  // TODO Ns? Or there should be as many seeds as vcells?
} s_vdiagram;


typedef struct vcell {
    int seed_id;
    // struct vcell *next;     // linked list of vcells
    int Nv;
    int Nv_capacity;  // malloc limit
    s_point *vertices;  // Nv x 3
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
    int *faces;  // Convex hull, 1 x 3*Nf
    s_point *fnormals;
    double volume;
} s_vcell;


void free_vdiagram(s_vdiagram *vdiagram);
void write_vd_file(const s_vdiagram *vd, FILE *file);
s_vdiagram *malloc_vdiagram(const s_scplx *setup, int Nreal);  // TODO what is Nreal?
void print_vdiagram(const s_vdiagram *vdiagram);

s_vcell *malloc_vcell(int seed_id);
void compute_vcell_volume(s_vcell *vcell);
s_vdiagram *voronoi_from_delaunay_3d(const s_scplx *setup, const s_bpoly *bpoly, int Nreal);
int find_inside_which_vcell(const s_vdiagram *vd, s_point x);
void plot_vcell(const s_vdiagram *vdiag, const s_vcell *vcell, char *f_name, const s_point ranges[2]);
void plot_vdiagram_differentviews(const s_vdiagram *vdiagram, char *f_name, const s_point ranges[2]);

void clear_volumes_file(char *fname);
void append_volumes_to_file(s_vdiagram *vdiagram, char *fname, int id);

#endif
