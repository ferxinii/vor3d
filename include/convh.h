#ifndef VOR3D_CONVH_H
#define VOR3D_CONVH_H

#include "external/geometry/geometry.h"

void convhull_from_points(const s_point *points, int Np, int **faces, s_point **fnormals, int *Nf);  // fnormals out is optional
// s_point *extract_normals_from_ch(const ch_vertex *vertices, int *faces, int Nf, s_point ch_CM, int NORMALIZE);  // TODO remove visibility of this! only used in externally in vdiagram.
int is_inside_convhull(s_point query, const s_point *pch, const int *faces, int Nf);  // 1: inside, 0: outise, -1: in boundary
int is_in_boundary_convhull(const int *faces, int Nf, int vid);
int mark_inside_convhull(const s_point *query, int Np, const s_point *pch, const int *faces, int Nf, int *mark);

s_point random_point_inside_convhull(const s_point *pch, const int *faces, int Nf, s_point min, s_point max);

double volume_tetrahedron_approx(s_point p1, s_point p2, s_point p3, s_point p4);
double compute_volume_convhull(const s_point *points, const int *faces, const s_point *fnormals_UNNORMALIZED, int Nf);
double compute_volume_convhull_from_points(const s_point *points, int Np);

#endif
