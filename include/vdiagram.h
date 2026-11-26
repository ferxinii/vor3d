#ifndef VOR3D_VDIAGRAM_H
#define VOR3D_VDIAGRAM_H

#include <float.h>
#include <stdbool.h>
#include <stdint.h>
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

/* A Voronoi cell is represented as an array of convex pieces: 1 piece in a
 * convex domain, >= 1 in a non-convex one. An EMPTY cell (seed outside the
 * domain whose cell misses it, or a weighted seed dominated everywhere in the
 * domain) keeps its slot with N_pieces = 0, pieces = NULL, volume = 0 and no
 * surface. */
typedef struct vcell {
    int        seed_id;
    s_convh   *pieces;
    int        N_pieces;
    double     volume;
    s_trimesh *surface;     /* array of the cell's boundary meshes, one per connected
                               component (a non-convex cell straddling a domain fold can
                               pinch into >1 lobe). NULL/0 if extraction failed. */
    int        N_surface;   /* number of surface components */
} s_vcell;

void free_vcell(s_vcell *c);

/* ---- Domain ---- */

/* The region a diagram is clipped to. Built by the vor3d_domain_from_* functions
 * (vor3d.h), which also decide is_convex; the convex kernel is used iff
 * is_convex, the non-convex (CDT) kernel otherwise. */
typedef struct vor3d_domain {
    bool       is_convex;
    s_trimesh  surface;   /* boundary mesh (owns memory), always present */
    s_bpoly    bpoly;     /* convex: the domain itself; non-convex: its convex hull */
    s_scplx    cdt;       /* non-convex only: CDT of the interior (interior tets,
                             plus flagged exterior pockets if built with
                             keep_exterior); empty (head == NULL) when is_convex */
    bool       cdt_has_exterior;  /* cdt was built with keep_exterior = 1 */
    double     volume;    /* convex: bpoly.volume; non-convex: sum of interior tets */
} s_vor3d_domain;


/* ---- Diagram ---- */

typedef struct vdiagram {
    s_points  seeds;      /* kept seeds, compacted (see out_kept_idx in vor3d.h) */
    s_vcell  *vcells;     /* length seeds.N, index-aligned with seeds */
    /* Combinatorial cell adjacency (Delaunay dual of the seed DT): cells i and j
     * are adjacent iff they share a Voronoi/power facet. CSR layout: cell i's
     * neighbours are adj_neighbor[adj_offset[i] .. adj_offset[i+1]-1]. */
    int      *adj_offset;     /* size seeds.N+1 */
    int      *adj_neighbor;   /* size adj_offset[seeds.N] */
    double    domain_volume;
    double    vol_rel_err;    /* (domain_volume - sum of cell volumes) / domain_volume */
    bool      is_convex;      /* kernel that built it */
} s_vdiagram;

void free_vdiagram(s_vdiagram *vd);


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

double find_closest_point_on_bp(const s_bpoly *bp, s_point p, double EPS_DEG, s_point *out);

void plot_bpoly_differentviews(s_bpoly *bpoly, char *f_name, s_point ranges[2], char *color);
void plot_bpoly(s_bpoly *bpoly, char *f_name, s_point ranges[2], char *color, char *view_command);


/* VORONOI DIAGRAM */
int vdiagram_is_valid(const s_vdiagram *vd);
void print_vcell(const s_vcell *vcell);
void print_vdiagram(const s_vdiagram *vdiagram);
/* Index of the cell containing x (any piece, boundary included), or -1. */
int find_inside_which_vcell(const s_vdiagram *vd, s_point x, double EPS_DEG, double TOL);

void plot_vcell(const s_vor3d_domain *domain, const s_vcell *vcell, char *f_name, const s_point ranges[2]);
void plot_all_vcells(const s_vdiagram *vdiagram, char *f_name, const s_point ranges[2], char *view_command);
void plot_vdiagram_differentviews(const s_vdiagram *vdiagram, char *f_name, const s_point ranges[2]);

#endif
