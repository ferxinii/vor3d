#ifndef VOR3D_BPOLY_H
#define VOR3D_BPOLY_H

#include "geometry.h"  
#include "convh.h"

// Used for Poisson disc sampling inside bpoly:
#define MAX_TRIAL_POINTS 10000
#define MAX_TRIAL_TESTS 50


typedef struct bounding_polyhedron {
    s_convhull convh;
    double dmax;  // Max distance between two pairs of points
    s_point CM;
    s_point min;
    s_point max;
    double volume;
} s_bpoly;


s_bpoly bpoly_from_points(const s_points *points);  // (Different copy of points inside)
s_bpoly bpoly_from_csv(const char *fname);
s_bpoly bpoly_copy(const s_bpoly *in);
s_bpoly copy_bpoly_scaled(const s_bpoly *bp, double factor);
s_bpoly copy_bpoly_scaled_volume(const s_bpoly *bp, double objective_volume);
void free_bpoly(s_bpoly *bpoly);

s_points generate_poisson_dist_inside(const s_bpoly *bpoly, double (*rmax)(double *));
void extend_sites_mirroring(const s_bpoly *bp, s_points *seeds_inout);
double find_closest_point_on_bp(const s_bpoly *bp, s_point p, s_point *out);

void generate_file_cube_bp(const char *filename, double length);
void generate_file_tetrahedron_bp(const char *filename, double length);
void generate_file_sphere_bp(const char *filename, double radius, int nTheta, int nPhi);

void plot_bpoly_differentviews(s_bpoly *bpoly, char *f_name, s_point ranges[2], char *color);
void plot_bpoly(s_bpoly *bpoly, char *f_name, s_point ranges[2], char *color, char *view_command);

#endif
