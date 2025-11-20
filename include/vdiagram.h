#ifndef VOR3D_VDIAGRAM_H
#define VOR3D_VDIAGRAM_H

#include <float.h>
#include "scplx.h"
#include "bpoly.h"
#include "convh.h"

#define VBUFF_N_INIT 1000  // Used for efficient mallocing, incrementing by blocks if necessary

typedef struct vdiagram {
    s_points seeds;  
    struct vcell *vcells;  // Array of vcells (size: seeds.N)
    s_bpoly bpoly;
} s_vdiagram;


typedef struct vcell {
    int seed_id;
    s_convh convh;
    double volume;
} s_vcell;


void free_vdiagram(s_vdiagram *vdiagram);
void print_vcell(const s_vcell *vcell);
void print_vdiagram(const s_vdiagram *vdiagram);

s_vdiagram voronoi_from_delaunay_3d(const s_scplx *setup, const s_bpoly *bpoly, int Nreal, double EPS_degenerate, double TOL);  // copy of bpoly inside
int find_inside_which_vcell(const s_vdiagram *vd, s_point x, double EPS_degenerate, double TOL);

void plot_vcell(const s_vdiagram *vdiag, const s_vcell *vcell, char *f_name, const s_point ranges[2]);
void plot_all_vcells(const s_vdiagram *vdiagram, char *f_name, const s_point ranges[2], char *view_command);
void plot_vdiagram_differentviews(const s_vdiagram *vdiagram, char *f_name, const s_point ranges[2]);

#endif
