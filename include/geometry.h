#ifndef VOR3D_GEOMETRY_H
#define VOR3D_GEOMETRY_H

#include "convhull_3d.h"

typedef struct point {
    union {
        double coords[3];
        struct {
            double x, y, z;
        };
    };
} s_point;


int orientation(const s_point *p3, s_point q);
int in_sphere(const s_point *p4, s_point q);

int segment_crosses_triangle_3d(const s_point *triangle, s_point a, s_point b);
int between_1d(double x, double a, double b, double eps);
int segments_intersect_2d(const s_point*AB, const s_point *pd);

int are_in_general_position_3d(double **points, int N);  // TODO remove or upgrade this?
                                                         //
s_point find_center_mass(const s_point *in, int N_points);
int coord_with_largest_component_3d(double *n);

s_point cross_prod(s_point u, s_point v);
double dot_prod(s_point u, s_point v);
s_point subtract_points(s_point u, s_point v);
double distance_squared(s_point a, s_point b);

s_point closest_point_on_triangle(const s_point *triangle, s_point p);
s_point closest_point_on_segment(const s_point *segment, s_point p);
int point_in_triangle_2d(const s_point *triangle, s_point p);

ch_vertex *malloc_points_to_chvertex(const s_point *points, int Np);
s_point *extract_normals_from_ch(const ch_vertex *vertices, int *faces, int Nf, s_point ch_CM, int NORMALIZE);

int is_inside_convhull(s_point query, const s_point *pch, const int *faces, int Nf);
int is_in_boundary_convhull(const int *faces, int Nf, int vid);
s_point random_point_uniform_3d(s_point min, s_point max);
s_point random_point_inside_convhull(const s_point *pch, const int *faces, int Nf, s_point min, s_point max);
int mark_inside_convhull(const s_point *query, int Np, const s_point *pch, const int *faces, int Nf, int *mark);

double volume_tetrahedron_approx(s_point p1, s_point p2, s_point p3, s_point p4);
double compute_volume_convhull(const s_point *points, const int *faces, const s_point *fnormals_UNNORMALIZED, int Nf);
double compute_volume_convhull_from_points(const s_point *points, int Np);

void extract_faces_convhull_from_points(const s_point *points, int Np, int **faces, s_point **fnormals, int *Nf);

#endif
