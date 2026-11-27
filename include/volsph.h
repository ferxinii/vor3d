#ifndef VOR3D_VOLSPH_H
#define VOR3D_VOLSPH_H

#include "points.h"
#include "vor3d_opts.h"

double volume_intersection_2_spheres(double rA, double rB, double dAB,
                                     bool check_intersection);
double volume_intersection_3_spheres(double rA, double rB, double rC, 
                                     double dAB, double dAC, double dBC,
                                     bool check_intersection);


double volume_sphere(double r);

/* Analytic unions of 2 / 3 balls (inclusion-exclusion).
 * check_intersection: handle containment / disjoint / no-triple-point cases;
 *   pass false only when the balls are known to overlap properly.
 * out_contrib (nullable, length 2 / 3, in argument order): per-ball power
 *   (Laguerre) contributions, same unit as volume_union_spheres; a contained
 *   ball gets 0.  Passing it never changes the returned volume. */
double volume_union_2_spheres(s_point pA, double rA,
                              s_point pB, double rB,
                              bool check_intersection, double *out_contrib /* nullable */);

double volume_union_3_spheres(s_point pA, double rA,
                              s_point pB, double rB,
                              s_point pC, double rC,
                              bool check_intersection,
                              double *out_contrib /* nullable */);

/* Volume of the union of N balls (centers, radii).  N = 1, 2, 3: analytic;
 * N >= 4: weighted Delaunay (regular) triangulation + alpha complex at alpha = 0
 * (Edelsbrunner's inclusion-exclusion).  Returns NAN on allocation failure.
 * opts (nullable): tol_dup, forwarded to the RT.
 * out_contrib (nullable, length centers->N, indexed by ORIGINAL ball id):
 *   per-ball contributions, sum_i out_contrib[i] == the returned volume (up to
 *   rounding; passing it never changes the returned value).  The per-ball unit
 *   is the Laguerre/power partition Vol(union INTERSECT power-cell(i)).
 *   Redundant/fully-contained balls (dropped by the RT) get 0.  An individual
 *   contribution is always >= 0; the K-tet power split may produce negative
 *   intermediate pieces that cancel, but the total per ball is a volume.  (The
 *   split uses the trig-free "quad form" of Mach A.16 and is applied to EVERY
 *   K-tet, including near-flat ones, whose signed pieces are load-bearing; see
 *   VOLSPH_SPLIT_PLAN.md.) */
double volume_union_spheres(const s_points *centers, const double radii[centers->N],
                            const s_vor3d_opts *opts, double *out_contrib /* nullable */);

#endif
