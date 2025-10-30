#ifndef VOR3D_VDIAGRAM_H
#define VOR3D_VDIAGRAM_H

#include <float.h>
#include "scplx.h"
#include "bpoly.h"
#include "convh.h"

#define VCELL_BLOCK_VERTICES 1000  // Used for efficient mallocing, incrementing by blocks if necessary

typedef struct vdiagram {
    s_points seeds;  
    struct vcell *vcells;  // Array of vcells (size: seeds.N)
    s_bpoly bpoly;
} s_vdiagram;


typedef struct vcell {
    int seed_id;
    s_convhull convh;
    double volume;
} s_vcell;


void free_vdiagram(s_vdiagram *vdiagram);
void write_vcell_file(const s_vcell *vcell, FILE *file);
void write_vd_file(const s_vdiagram *vd, FILE *file);
void print_vcell(const s_vcell *vcell);
void print_vdiagram(const s_vdiagram *vdiagram);

s_vdiagram voronoi_from_delaunay_3d(const s_scplx *setup, const s_bpoly *bpoly, int Nreal);  // copy of bpoly inside
int find_inside_which_vcell(const s_vdiagram *vd, s_point x);

void plot_vcell(const s_vdiagram *vdiag, const s_vcell *vcell, char *f_name, const s_point ranges[2]);
void plot_all_vcells(const s_vdiagram *vdiagram, char *f_name, const s_point ranges[2], char *view_command);
void plot_vdiagram_differentviews(const s_vdiagram *vdiagram, char *f_name, const s_point ranges[2]);

void clear_volumes_file(char *fname);
void append_volumes_to_file(s_vdiagram *vdiagram, char *fname, int id);

#endif
