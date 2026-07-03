#ifndef VOR3D_VDIAGRAM_H
#define VOR3D_VDIAGRAM_H

#include <float.h>
#include <stdbool.h>
#include "scplx.h"
#include "convh.h"
#include "trimesh.h"

// #define VBUFF_N_INIT 1000  // Used for efficient mallocing, incrementing by blocks if necessary

typedef struct bounding_polyhedron {
    s_convh convh;
    s_point CM;
    s_point min;
    s_point max;
    double volume;
} s_bpoly;

/* A Voronoi cell is represented as an array of convex pieces.
 * For convex domains N_pieces is always 1; for non-convex domains N_pieces >= 1. */
typedef struct vcell {
    int      seed_id;
    s_convh *pieces;
    int      N_pieces;
    double   volume;
} s_vcell;

void free_vcell(s_vcell *c);

typedef struct vdiagram {
    s_points seeds;
    s_bpoly bpoly;
    struct vcell *vcells;  // Array of vcells (size: seeds.N)
} s_vdiagram;


/* ---- Non-convex domain and diagram ---- */

typedef struct ncvx_domain {
    s_trimesh   surface;        /* boundary mesh (owns memory) */
    s_scplx     cdt;            /* CDT of the interior (interior tets; owns memory) */
    double      domain_volume;  /* sum of interior tet volumes */
    s_bpoly     bpoly;          /* convex hull of surface.points; used for mirroring */
    void       *spatial_index;  /* NULL in MVP */
} s_ncvx_domain;

typedef struct ncvx_vdiagram {
    s_points      seeds;    /* copy of seeds */
    s_ncvx_domain domain;   /* copy of domain */
    s_vcell      *vcells;   /* array of length seeds.N */
} s_ncvx_vdiagram;

/*
 * Build a domain from a trimesh.
 * Tetrahedralizes the interior and computes the convex hull for mirroring.
 * Returns a zero-initialised s_ncvx_domain on error.
 */
s_ncvx_domain ncvx_domain_from_trimesh(const s_trimesh *mesh,
                                        double EPS_DEG, double TOL);
int  ncvx_domain_is_valid(const s_ncvx_domain *d);
void free_ncvx_domain(s_ncvx_domain *d);
void free_ncvx_vdiagram(s_ncvx_vdiagram *vd);


/* BOUNDING POLYHEDRON */
/* (New copy of points inside) */
s_bpoly bpoly_from_points(const s_points *points, double EPS_DEG);
s_bpoly bpoly_from_csv(const char *fname, double EPS_DEG);
s_bpoly bpoly_from_convh(const s_convh *convh);
s_bpoly bpoly_from_convh_scaled(const s_convh *convh, double s, s_point pivot,
                                double EPS_DEG);
s_bpoly bpoly_copy(const s_bpoly *in);
s_bpoly bpoly_copy_scaled(const s_bpoly *bp, double scale, s_point pivot,
                          double EPS_DEG);
s_bpoly bpoly_copy_scaled_volume(const s_bpoly *bp, double objective_volume,
                                 double EPS_DEG);
void free_bpoly(s_bpoly *bpoly);

s_points generate_poisson_dist_inside(const s_bpoly *bpoly,
                                      double (*rmax)(double*, void*), void *rmax_params,
                                      double (*randd01)(void*), int (*randint)(void*, int),
                                      void *rctx, double EPS_DEG);
double find_closest_point_on_bp(const s_bpoly *bp, s_point p, double EPS_DEG, s_point *out);

void plot_bpoly_differentviews(s_bpoly *bpoly, char *f_name, s_point ranges[2], char *color);
void plot_bpoly(s_bpoly *bpoly, char *f_name, s_point ranges[2], char *color, char *view_command);

int extract_mirrored_points(const s_bpoly *bp, double EPS_degenerate,
                            s_scplx *dt, int N_seeds,
                            int (*randint)(void *rctx, int), void *rctx,
                            s_dynarray *out_points);

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
int find_inside_which_vcell(const s_vdiagram *vd, s_point x, double EPS_DEG, double TOL);

s_vdiagram voronoi_from_delaunay_3d(const s_scplx *setup, const s_bpoly *bpoly, int Nreal,
                                    double EPS_DEG, double TOL);

void plot_vcell(const s_vdiagram *vdiag, const s_vcell *vcell, char *f_name, const s_point ranges[2]);
void plot_all_vcells(const s_vdiagram *vdiagram, char *f_name, const s_point ranges[2], char *view_command);
void plot_vdiagram_differentviews(const s_vdiagram *vdiagram, char *f_name, const s_point ranges[2]);

#endif
