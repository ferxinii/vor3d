#ifndef VOR3D_VDIAGRAM_H
#define VOR3D_VDIAGRAM_H

#include <float.h>
#include "simplical_complex.h"
#include "bpoly.h"

#define VCELL_BLOCK_VERTICES 1000  // Used for efficient mallocing, incrementing by blocks if necessary

typedef struct vdiagram {
    int N;  // Number of vcells
    s_point *seeds;  
    struct vcell **vcells;  // Array of pointers to the cells
    const struct bounding_polyhedron *bpoly;
} s_vdiagram;


typedef struct vcell {
    int seed_id;
    int Nv;
    int Nv_capacity;  // malloc limit
    s_point *vertices;  // Nv x 3
    int Nf;
    int *faces;  // Convex hull, 1 x 3*Nf
    s_point *fnormals;
    double volume;
} s_vcell;


void free_vdiagram(s_vdiagram *vdiagram);
void write_vd_file(const s_vdiagram *vd, FILE *file);
s_vdiagram *malloc_vdiagram(const s_scplx *setup, int Nreal);  // Nreal is the N seeds without reflection!
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
