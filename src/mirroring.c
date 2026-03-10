
#include "mirroring.h"
#include "dynarray.h"
#include "bpoly.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>


static s_point mirror_plane(s_point normal, double d_plane, s_point p)
{
    double factor = 2 * (d_plane - dot_prod(normal, p));
    s_point out;
    out.x = p.x + factor * normal.x;
    out.y = p.y + factor * normal.y;
    out.z = p.z + factor * normal.z;
    return out;
}


/* Linear programming in 2D. */
typedef struct {  /* Constraint: a*x + b*y <= c */
    double a, b, c;
} s_constraint;

typedef struct {
    double x, y;
} s_point_2D;

static bool constraint_valid(s_constraint con, double TOL2)
{
    double n2 = con.a*con.a + con.b*con.b;
    if (n2 <= TOL2) return false;
    else return true;
}

static bool point_satisfies_constraint(s_constraint con, s_point_2D p, double TOL) 
{   /* a*x + b*y <= c + TOL */
    double val = con.a * p.x + con.b * p.y - con.c;
    return (val <= TOL);
}

static void line_from_constraint(s_constraint con, s_point_2D *out_p, s_point_2D *out_v) 
{   /* Returns point on line p and direction v */
    *out_v = (s_point_2D){-con.b, con.a}; 
    /* To find p, set either x=0 or y=0 */
    if (fabs(con.b) > fabs(con.a)) *out_p = (s_point_2D){0, con.c/con.b};
    else *out_p = (s_point_2D){con.c/con.a, 0}; 
}

static void constraints_from_triangle(const double t1[2], const double t2[2], const double t3[2], s_constraint out[3])
{
    double v1[2] = { t1[0], t1[1] };
    double v2[2] = { t2[0], t2[1] };
    double v3[2] = { t3[0], t3[1] };

    if (orient2d(t1, t2, t3) < 0) {  /* Check proper orientation */
        double tmp[2] = {v1[0], v1[1]};
        v1[0] = v2[0];  v1[1] = v2[1];
        v2[0] = tmp[0]; v2[1] = tmp[1];
    }

    out[0] = (s_constraint){ v2[1]-v1[1], -(v2[0]-v1[0]), (v2[1]-v1[1])*v1[0] - (v2[0]-v1[0])*v1[1] };
    out[1] = (s_constraint){ v3[1]-v2[1], -(v3[0]-v2[0]), (v3[1]-v2[1])*v2[0] - (v3[0]-v2[0])*v2[1] };
    out[2] = (s_constraint){ v1[1]-v3[1], -(v1[0]-v3[0]), (v1[1]-v3[1])*v3[0] - (v1[0]-v3[0])*v3[1] };
}

static void shuffle_constraints(int N, s_constraint constraints[N])
{
	for (int i = N - 1; i > 0; --i) {
		int j = rand() % (i + 1);
		s_constraint tmp = constraints[i];
		constraints[i] = constraints[j];
		constraints[j] = tmp;
	}
}

static bool LP_is_feasible(int N, s_constraint constraints[N], double TOL) 
{   /* Seidel algorithm */
    const double TOL2 = TOL * TOL;
    shuffle_constraints(N, constraints);
    s_point_2D x = {0.0, 0.0};  /* Initial guess */
    for (int ii=0; ii<N; ii++) {  
        if (!constraint_valid(constraints[ii], TOL2)) {
            if (constraints[ii].c < -TOL) return false;
            else continue;
        }
        if (point_satisfies_constraint(constraints[ii], x, TOL)) continue;

        /* Find new x on boundary of ii satisfiying all previous constraints */
        s_point_2D p, v;
        line_from_constraint(constraints[ii], &p, &v);
        double norm_v = sqrt(v.x*v.x + v.y*v.y);

        double t_min = -INFINITY;
        double t_max = +INFINITY;
        double TOL_t = (norm_v > 0.0) ? (TOL / norm_v) : TOL;
        for (int jj=0; jj<ii; jj++) {
            double alpha = constraints[jj].a * v.x + constraints[jj].b * v.y;
            double beta  = constraints[jj].c - (constraints[jj].a * p.x + constraints[jj].b * p.y);
            beta += TOL;  /* To be consistent with inclusion test */

            if (fabs(alpha) <= TOL * norm_v) {  /* alpha ~ 0 */
                if (beta < -TOL) return false;  /* Impossible */
                else continue;  /* jj places no restriction on this line */
            }

            double t_bound = beta / alpha;
            if (alpha > 0.0) {  /* t <= t_bound */
                if (t_bound < t_max) t_max = t_bound;
            } else {  /* alpha < 0 -> t >= t_bound */
                if (t_bound > t_min) t_min = t_bound;
            }

            /* check emptiness with t-tolerance mapped from coordinate TOL */
            if (t_min > t_max + TOL_t) return false;
        } 

        double t;
        if (isinf(t_min) && isinf(t_max)) t = 0.0;
        else if (isinf(t_min)) t = t_max - 1.0;
        else if (isinf(t_max)) t = t_min + 1.0;
        else t = 0.5 * (t_min + t_max);

        x.x = p.x + t * v.x;
        x.y = p.y + t * v.y;
    } 
    // *out_point = x;
    return true;
}


/* Mirroring criteria based on weighted 2D Voronoi diagram */
static s_point project_to_plane_3D(s_point p, s_point plane_n, double plane_d)
{
    double dist = dot_prod(p, plane_n) - plane_d;
    return subtract_points(p, scale_point(plane_n, dist));
}

static inline void plane_3D_to_2D(s_point X, s_point O, s_point u, s_point v, double out[2])
{
    s_point d = subtract_points(X, O);
    out[0] = dot_prod(d, u);
    out[1] = dot_prod(d, v);
}

static bool should_mirror(const s_point face[3], const s_points *seeds, int best_n, int *best_idx, int id, double TOL, s_dynarray *buff_constraints)
{   /* Checks if the weighted voronoi diagram in the plane of face intersects the face */
    /* Project onto 2D */
    s_point n, u, v;
    if (!basis_vectors_plane(face, TOL, &n, &u, &v)) {
        printf("DEBUG: Could not construct basis for plane of face.\n");
        return false;
    }
    double plane_d = dot_prod(n, face[0]);
    s_point O = face[0];  /* Plane origin */
    double t1[2];  plane_3D_to_2D(face[0], O, u, v, t1);
    double t2[2];  plane_3D_to_2D(face[1], O, u, v, t2);
    double t3[2];  plane_3D_to_2D(face[2], O, u, v, t3);
    double s[2];   plane_3D_to_2D(project_to_plane_3D(seeds->p[best_idx[id]], n, plane_d), O, u, v, s);

    /* Find closest point to plane to determine appropriate plane offset */
    // s_point closest = project_to_plane_3D(seeds->p[0], n, plane_d);
    // double d2min = distance_squared(seeds->p[0], closest);
    // for (int ii=1; ii<seeds->N; ii++) {
    //     closest = project_to_plane_3D(seeds->p[ii], n, plane_d);
    //     double d2 = distance_squared(seeds->p[ii], closest);
    //     if (d2 < d2min) d2min = d2;
    // }
    // double plane_offset = sqrt(d2min) / 3;
    double plane_offset = 0;

    double ds = signed_distance_point_to_plane(seeds->p[best_idx[id]], face, TOL);
    double ws = ds*ds - 2*fabs(ds)*plane_offset;

    /* Add constraints */
    s_dynarray *C = buff_constraints;
    C->N = 0;
    for (int ii=0; ii<best_n; ii++) {
        if (ii == id) continue;
        double t[2];   plane_3D_to_2D(project_to_plane_3D(seeds->p[best_idx[ii]], n, plane_d), O, u, v, t);
        double dt = signed_distance_point_to_plane(seeds->p[best_idx[ii]], face, TOL);
        double wt = dt*dt - 2*fabs(dt)*plane_offset;
        s_constraint con = { .a = 2*(t[0]-s[0]),
                             .b = 2*(t[1]-s[1]),
                             .c = (t[0]*t[0]+t[1]*t[1] + wt)-(s[0]*s[0]+s[1]*s[1] + ws) };
        if (!dynarray_push(C, &con)) goto error;
    }

    /* Add triangle as constraints */
    s_constraint triangle_con[3];  constraints_from_triangle(t1, t2, t3, triangle_con);
    if (!dynarray_push(C, &triangle_con[0])) goto error;
    if (!dynarray_push(C, &triangle_con[1])) goto error;
    if (!dynarray_push(C, &triangle_con[2])) goto error;

    /* Check feasability */
    return LP_is_feasible(C->N, C->items, TOL);

    error:
        return 0;
}


static int should_mirror_OLD(s_point normal, double d_plane, s_point face[3], const s_points *all_seeds, int seed_id, double EPS_degenerate)
{
    s_point s = all_seeds->p[seed_id];
    double dist = d_plane - dot_prod(normal, s);
    s_point proj_fplane = {{{s.x + dist*normal.x, s.y + dist*normal.y, s.z + dist*normal.z}}};
    s_point c = closest_point_on_triangle(face, EPS_degenerate, proj_fplane);
    if (!point_is_valid(c)) return 0;
    
    // Check nearest‐neighbor at c
    double d_s = distance_squared(c, s);
    for (int jj=0; jj<all_seeds->N; jj++) {
        if (jj != seed_id && distance(all_seeds->p[jj], c) + 1e-8 < d_s)
            return 0;   // someone else is nearer
    }
    return 1;
}


static double distance_plane(s_point n, double d, s_point p)
{
    return fabs(dot_prod(n, p) - d);
}

static void insert_topk(int idx[], double val[], int *n, int k, int new_idx, double new_val)
{   /* Keep k smallest distances (descending order: worst at 0). Only efficient for small k. */
	if (*n < k) {
		int pos = (*n)++;
		idx[pos] = new_idx;
		val[pos] = new_val;
		/* bubble up */
		while (pos > 0 && val[pos] > val[pos-1]) {
			double tv = val[pos]; val[pos] = val[pos-1]; val[pos-1] = tv;
			int ti = idx[pos]; idx[pos] = idx[pos-1]; idx[pos-1] = ti;
			pos--;
		}
	} else if (new_val < val[0]) {
		/* replace worst */
		idx[0] = new_idx;
		val[0] = new_val;
		/* bubble down */
		int pos = 0;
		while (pos+1 < *n && val[pos] < val[pos+1]) {
			double tv = val[pos]; val[pos] = val[pos+1]; val[pos+1] = tv;
			int ti = idx[pos]; idx[pos] = idx[pos+1]; idx[pos+1] = ti;
			pos++;
		}
	}
}



int extend_sites_mirroring_initial(const s_bpoly *bp, double EPS_degenerate, double TOL, s_points *inout_seeds)
{   /* 0 ERROR, 1 OK */
    (void)should_mirror, (void)should_mirror_OLD, (void)TOL;  // Here for when debugging avoiding warnings
    s_dynarray buff_constraints = dynarray_initialize(sizeof(s_constraint), 20);
    s_dynarray mirrored = dynarray_initialize(sizeof(s_point), bp->convh.Nf);

    for (int ff=0; ff<bp->convh.Nf; ff++) {
        // printf("ff: %d / %d. inout_seeds.N : %d\n", ff, bp->convh.Nf, inout_seeds->N);
        s_point face[3]; 
        convh_get_face(&bp->convh, ff, face);
        s_point normal = normalize_vec(bp->convh.fnormals[ff], EPS_degenerate);
        if (!point_is_valid(normal)) continue;
        double d_plane =  dot_prod(normal, face[0]);

        /* Find k points closest to face, and only these are candidates to mirror */
        int k = 10;
        int best_idx[k];
        double best_dist[k];
        int best_n = 0;
        for (int pp=0; pp < inout_seeds->N; pp++) {
            double d = distance_plane(normal, d_plane, inout_seeds->p[pp]);
            insert_topk(best_idx, best_dist, &best_n, k, pp, d);
        }
        
        for (int ii=0; ii<best_n; ii++) {
            if (should_mirror(face, inout_seeds, best_n, best_idx, ii, TOL, &buff_constraints)) {
                s_point p_mirror = mirror_plane(normal, d_plane, inout_seeds->p[best_idx[ii]]);
                if (!dynarray_push(&mirrored, &p_mirror)) goto error;
            }
        }
    }
    // printf("Nadded %d\n", mirrored.N);


    /* Add mirrored to inout_seeds */
    // printf("DEBUG: Mirrored.N = %d, total = %d\n", mirrored.N, inout_seeds->N + mirrored.N);
    s_point *tmp = realloc(inout_seeds->p, sizeof(s_point) * (inout_seeds->N + mirrored.N));
    if (!tmp) goto error;
    memcpy(&tmp[inout_seeds->N], mirrored.items, mirrored.N * sizeof(s_point));
    inout_seeds->p = tmp;
    inout_seeds->N += mirrored.N;

    dynarray_free(&buff_constraints);
    dynarray_free(&mirrored);
    return 1;
    
    error:
        fprintf(stderr, "Error in exten_sites_mirroring!\n");
        dynarray_free(&buff_constraints);
        dynarray_free(&mirrored);
        return 0;
}


