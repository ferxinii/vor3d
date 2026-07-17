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

/* Per-ball contributions to Vol(union): out_contrib has length centers->N,
 * indexed by ORIGINAL ball id, with sum_i out_contrib[i] == volume_union_spheres(
 * centers, radii, EPS_DEGEN, TOL_dup, buff_ncellPTR) (up to rounding).  The
 * per-ball unit is the Laguerre/power partition Vol(union INTERSECT
 * power-cell(i)).  Redundant/fully-contained balls (dropped by the RT) get 0.
 * An individual contribution is always >= 0; the K-tet power split may produce
 * negative intermediate pieces that cancel, but the total per ball is a volume.
 * (The split uses the trig-free "quad form" of Mach A.16 and is applied to
 * EVERY K-tet, including near-flat ones, whose signed pieces are load-bearing;
 * see VOLSPH_SPLIT_PLAN.md.) */
void volume_contribution_spheres(const s_points *centers, const double radii[centers->N],
                                 double EPS_DEGEN, double TOL_dup,
                                 struct dynarray *buff_ncellPTR,
                                 double out_contrib[centers->N]);

#endif
