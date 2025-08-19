#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "convhull_3d.h"

int orientation(double **p, double *q, int dim);
int in_sphere(double **p, double *q, int dim);
int segment_crosses_triangle_3d(double **triangle, double *a, double *b);
int between_1d(double x, double a, double b, double eps);
int segments_intersect_2d(double *A, double *B, double *p, double *d);
int are_in_general_position_3d(double **points, int N);
void find_center_mass(double **in, int N_points, int dim, double *out);
int coord_with_largest_component_3d(double *n);
void cross_3d(const double *u, const double *v, double *out);
double dot_3d(const double *u, const double *v);
void subtract_3d(const double *u, const double *v, double *result);
double distance_squared(const double *a, const double *b);
void closest_point_on_triangle(const double *A, const double *B, const double *C, const double *p, double *c_out);
void closest_point_on_segment(double *p, double *A, double *B, double *OUT);
int point_in_triangle_2d_NEW(double *v1, double *v2, double *v3, double *p);
int point_in_triangle_2d(double *v1, double *v2, double *v3, double *p);
int point_in_triangle_2d_candegenerate(double *v1, double *v2, double *v3, double *p);
int ray_triangle_intersection_3d(double **triangle, const double *origin, const double *dir, double *intersection);
ch_vertex *convert_points_to_chvertex(double **points, int Np);
double **extract_normals_from_ch(ch_vertex *vertices, int *faces, int Nf, double *ch_CM);
double **extract_normals_from_ch_UNNORMALIZED(ch_vertex *vertices, int *faces, int Nf, double *ch_CM);
int is_inside_convhull(double *query, double **pch, int *faces, double **fnormals, int Nf);
void inside_ray_convhull_intersection(double **pch, int *faces, double **fnormals, int Nf, double *origin, double *dir, double *OUT);
int is_in_boundary_convhull(int *faces, int Nf, int vid);
void random_point_uniform_3d(double *min, double *max, double *out);
void random_point_inside_convhull(double **pch, int *faces, double **fnormals, int Nf, double *min, double *max, double *out);
int mark_inside_convhull(double **points, int Np, double **pch, int *faces, double **fnormals, int Nf, int *mark);
double volume_tetrahedron_approx(double *p1, double *p2, double *p3, double *p4);
double compute_volume_convhull(double **points, int *faces, double **fnormals, int Nf);
double compute_volume_convhull_from_points(double **points, int Np);
void extract_faces_convhull_from_points(double **points, int Np, int **faces, double ***fnormals, int *Nf);
void segment_convex_hull_intersection(const double *p0, const double *p1, double **pch, const int *faces, double **fnormals, int Nf, int *ind1, double *i1, int *ind2, double *i2);

#endif
