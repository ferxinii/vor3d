#ifndef BPOLY_H
#define BPOLY_H

#include "float.h"
// #define TPH_POISSON_IMPLEMENTATION
#include "poisson_disk.h"


#define MAX_TRIAL_POINTS 10000
#define MAX_TRIAL_TESTS 50
typedef struct point {
    double coords[3];
} s_point;


typedef struct bound_poly {
    int Np;
    double **points;
    int Nf;
    int *faces;  // Its flat! Nf x 3
    double **fnormals;
    double dmax;  // Max distance between two pairs of points
    double CM[3];
    double min[3];
    double max[3];
    double volume;
} s_bound_poly;


void free_bpoly(s_bound_poly *bpoly);
void add_noise_to_bp(s_bound_poly *bpoly);
void extract_dmax_bp(s_bound_poly *bpoly);
void extract_CM_bp(s_bound_poly *bpoly);
void extract_min_max_coord(s_bound_poly *bpoly, double *min, double *max);
void extract_convhull_bp(s_bound_poly *bpoly);
s_bound_poly *new_bpoly_from_points(double **points, double Np, int add_noise);
void new_bpoly_from_txt(const char *fname, double ***OUT_points, int *OUT_Np, s_bound_poly **OUT_bpoly, int add_noise);
s_bound_poly *new_bpoly_copy(s_bound_poly *in);
void extract_vertices_face_bpoly(const s_bound_poly *bpoly, int *face, double **out);
void scale_bpoly_vertices(double **points, int Np, double s);
void scale_bpoly(s_bound_poly **bp, double objective_volume);



void find_closest_point_on_bp(s_bound_poly *bp, double *p, double *OUT);



int should_mirror(double *n, double *s, double d, double *f1, double *f2, double *f3, double **all_seeds, int Ns, int seed_id);
int extend_sites_mirroring(s_bound_poly *bp, double ***s, int Ns);


void random_point_uniform(double *min, double *max, s_point *out);
void random_point_around(double *x, double r, double *out);
int is_valid(s_bound_poly *bpoly, double *q, s_point *samples, int Nsamples, double (*rmax)(double *));
double **generate_nonuniform_poisson_dist_inside(s_bound_poly *bpoly, double (*rmax)(double *), int *Np_generated);


void generate_file_cube_bp(const char *filename, double length);
void generate_file_tetrahedron_bp(const char *filename, double length);
void generate_file_sphere_bp(const char *filename, double radius, int nTheta, int nPhi);
void plot_bpoly_with_points(s_bound_poly *bpoly, double **points, int Np, char *f_name, double *ranges, char *color);
#endif
