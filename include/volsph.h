#ifndef VOR3D_VOLSPH_H
#define VOR3D_VOLSPH_H

#include "points.h"
typedef struct dynarray s_dynarray;

double volume_intersection_2_spheres(double rA, double rB, double dAB,
                                     bool check_intersection);
double volume_intersection_3_spheres(double rA, double rB, double rC, 
                                     double dAB, double dAC, double dBC,
                                     bool check_intersection, double EPS_DEGEN);


double volume_sphere(double r);
double volume_union_2_spheres(s_point pA, double rA,
                              s_point pB, double rB,
                              bool check_intersection);

double volume_union_3_spheres(s_point pA, double rA,
                              s_point pB, double rB,
                              s_point pC, double rC,
                              bool check_intersection, double EPS_DEGEN);

double volume_union_spheres(const s_points *centers, const double radii[centers->N],
                            double EPS_DEGEN, double TOL_dup, struct dynarray *buff_ncellPTR);

#endif
