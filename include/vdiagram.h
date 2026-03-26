#ifndef VOR3D_VDIAGRAM_H
#define VOR3D_VDIAGRAM_H

#include <float.h>
#include <stdbool.h>
#include "scplx.h"
#include "convh.h"

// #define VBUFF_N_INIT 1000  // Used for efficient mallocing, incrementing by blocks if necessary

typedef struct bounding_polyhedron {
    s_convh convh;
    // double dmax;  // Max distance between two pairs of points
    s_point CM;
    s_point min;
    s_point max;
    double volume;
} s_bpoly;

typedef struct vcell {
    int seed_id;
    s_convh convh;
    double volume;
} s_vcell;

typedef struct vdiagram {
    s_points seeds;  
    s_bpoly bpoly;
    struct vcell *vcells;  // Array of vcells (size: seeds.N)
} s_vdiagram;



/* BOUNDING POLYHEDRON */
/* (New copy of points inside) */

s_bpoly bpoly_from_points(const s_points *points, double EPS_degenerate, double TOL);  
s_bpoly bpoly_from_csv(const char *fname, double EPS_degenerate, double TOL);
s_bpoly bpoly_from_convh(const s_convh *convh);
s_bpoly bpoly_from_convh_scaled(const s_convh *convh, double s, s_point pivot,
                                double EPS_degenerate, double TOL);
s_bpoly bpoly_copy(const s_bpoly *in);
s_bpoly bpoly_copy_scaled(const s_bpoly *bp, double scale, s_point pivot,
                          double EPS_degenerate, double TOL);
s_bpoly bpoly_copy_scaled_volume(const s_bpoly *bp, double objective_volume,
                                 double EPS_degenerate, double TOL);
void free_bpoly(s_bpoly *bpoly);

s_points generate_poisson_dist_inside(const s_bpoly *bpoly, 
                                      double (*rmax)(double*, void*), void *rmax_params,
                                      double (*randd01)(void*), int (*randint)(void*, int),
                                      void *rctx, double EPS_DEG);
double find_closest_point_on_bp(const s_bpoly *bp, s_point p, double EPS_degenerate, s_point *out);

void plot_bpoly_differentviews(s_bpoly *bpoly, char *f_name, s_point ranges[2], char *color);
void plot_bpoly(s_bpoly *bpoly, char *f_name, s_point ranges[2], char *color, char *view_command);

int extend_sites_mirroring(const s_bpoly *bp, double EPS_degenerate, double TOL, 
                           s_points *inout_seeds, int (*randint)(void* rctx, int), void* rctx,
                           s_dynarray *buff_points, s_dynarray *buff_LPconstraints2D);


/* VORONOI DIAGRAM */
int serialize_vdiagram(const s_vdiagram *vd, uint8_t *buff_write, size_t *size,
                       uint8_t **out);
int deserialize_vdiagram(const uint8_t *data, s_vdiagram *out, size_t *bytes_read);
int write_serialized_vdiagram(const char *file, const uint8_t *data, size_t size);
int read_serialized_vdiagram(const char *file, uint8_t **outbuf, size_t *outsize);

int vdiagram_is_valid(const s_vdiagram *vd);
void free_vdiagram(s_vdiagram *vdiagram);
void print_vcell(const s_vcell *vcell);
void print_vdiagram(const s_vdiagram *vdiagram);
int find_inside_which_vcell(const s_vdiagram *vd, s_point x, double EPS_degenerate, double TOL);

s_vdiagram voronoi_from_delaunay_3d(const s_scplx *setup, const s_bpoly *bpoly, int Nreal, 
                                    double EPS_degenerate, double TOL);

void plot_vcell(const s_vdiagram *vdiag, const s_vcell *vcell, char *f_name, const s_point ranges[2]);
void plot_all_vcells(const s_vdiagram *vdiagram, char *f_name, const s_point ranges[2], char *view_command);
void plot_vdiagram_differentviews(const s_vdiagram *vdiagram, char *f_name, const s_point ranges[2]);

#endif
