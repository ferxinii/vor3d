#ifndef CDT_PREDICATES_H
#define CDT_PREDICATES_H

#ifdef __cplusplus
extern "C" {
#endif

/* Must be called once before any predicate (sets FPU rounding mode). */
void cdt_predicates_init(void);

/* Reset the point registry (call before each DT build). */
void cdt_predicates_clear(void);

/* Register a new explicit point at DT index id. */
void cdt_point_set_explicit(int id, double x, double y, double z);

/* Register a Steiner point at DT index id.
 * Represents t*V[v1] + (1-t)*V[v2]  (paper convention).
 * v1 and v2 must already be registered as explicit points. */
void cdt_point_set_lnc(int id, int v1, int v2, double t);

/* orient3D -- returns -1, 0, or +1.
 * Arguments are DT point indices; any may be explicit or LNC. */
int cdt_orient3d(int a, int b, int c, int d);

/* inSphere -- returns -1 (inside), 0 (on), +1 (outside).
 * Arguments are DT point indices; any may be explicit or LNC. */
int cdt_insphere(int a, int b, int c, int d, int e);

/* SoS-wrapped inSphere for gift-wrapping (Algorithm 1 from paper).
 * i1..i4 are the tet vertices, i5 is the query point.
 * Returns -1 if i5 is inside or on circumsphere, +1 otherwise.
 * Never returns 0 (ties broken by global index order). */
int cdt_perturbed_insphere(int i1, int i2, int i3, int i4, int i5);

/* TRUE if point p is inside or on the boundary of tet <t0,t1,t2,t3>. */
int cdt_point_in_tet(int p, int t0, int t1, int t2, int t3);

/* TRUE if segment <s1,s2> intersects triangle <v0,v1,v2>
 * (not coplanar; interior of segment only). */
int cdt_segment_crosses_triangle(int s1, int s2, int v0, int v1, int v2);

/* TRUE if interior of s1-s2 crosses interior of <v0,v1,v2> at a single point.
 * Returns false if either endpoint is coplanar with the triangle. */
int cdt_inner_seg_crosses_inner_tri(int s1, int s2, int v0, int v1, int v2);

/* TRUE if interior of coplanar segments a-b and p-q cross at a single point.
 * Points are assumed coplanar; tries all three axis projections. */
int cdt_inner_segs_cross(int a, int b, int p, int q);

/* TRUE if point p is strictly inside (not on boundary of) triangle <a,b,c>.
 * Points are assumed coplanar; tries all three axis projections. */
int cdt_point_in_inner_tri(int p, int a, int b, int c);

/* Point p vs triangle <a,b,c>: -1 OUTSIDE, 0 on BOUNDARY, +1 strictly INSIDE.
 * Mirrors test_point_in_triangle_3D(.,0,0) (OUT if not coplanar). */
int cdt_point_in_triangle(int p, int a, int b, int c);

/* TRUE if the CLOSED segments <s1,s2> and <p,q> cross (coplanar inputs). */
int cdt_segments_cross(int s1, int s2, int p, int q);

/* Get approximate float coordinates of any registered point (explicit or LNC).
 * For DT insertion of Steiner points. Returns 0 if id is out of range. */
int cdt_get_approx_coords(int id, double *x, double *y, double *z);

#ifdef __cplusplus
}
#endif

#endif /* CDT_PREDICATES_H */
