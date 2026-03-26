
#include "dynarray.h"
#include "LP.h"
#include "points.h"
#include "gtests.h"
#include "vdiagram.h"
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


/* Linear Programming. */
static void constraints_from_triangle(const double t1[2], const double t2[2], const double t3[2], 
                                      s_LPconstraint2D out[3])
{
    double v1[2] = { t1[0], t1[1] };
    double v2[2] = { t2[0], t2[1] };
    double v3[2] = { t3[0], t3[1] };

    if (orient2d(t1, t2, t3) < 0) {  /* Check proper orientation */
        double tmp[2] = {v1[0], v1[1]};
        v1[0] = v2[0];  v1[1] = v2[1];
        v2[0] = tmp[0]; v2[1] = tmp[1];
    }

    out[0] = (s_LPconstraint2D){ v2[1]-v1[1], -(v2[0]-v1[0]), (v2[1]-v1[1])*v1[0] - (v2[0]-v1[0])*v1[1] };
    out[1] = (s_LPconstraint2D){ v3[1]-v2[1], -(v3[0]-v2[0]), (v3[1]-v2[1])*v2[0] - (v3[0]-v2[0])*v2[1] };
    out[2] = (s_LPconstraint2D){ v1[1]-v3[1], -(v1[0]-v3[0]), (v1[1]-v3[1])*v3[0] - (v1[0]-v3[0])*v3[1] };
}



/* Mirroring criteria based on weighted 2D Voronoi diagram */
static bool should_mirror(const s_point face[3], const s_points *seeds, int best_n, int *best_idx,
                          int id, double TOL, int (*randint)(void* rctx, int), void* rctx,
                          s_dynarray *buff_3dbls)
{   /* Checks if the weighted voronoi diagram in the plane of face intersects the face */
    /* Project onto 2D */
    s_point n, u, v;
    if (!basis_vectors_plane(face, TOL, &n, &u, &v)) {
        fprintf(stderr, "should_mirror: Could not construct basis for plane of face. Ignoring face.\n");
        return false;
    }
    double t1[2] = {0, 0};  /* Plane origin is face[0] */
    double t2[2];  project_point_to_plane_2D(face[1], face[0], u, v, t2);
    double t3[2];  project_point_to_plane_2D(face[2], face[0], u, v, t3);
    double s[2];   project_point_to_plane_2D(seeds->p[best_idx[id]], face[0], u, v, s);

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
    s_dynarray *C = buff_3dbls;
    C->N = 0;
    for (int ii=0; ii<best_n; ii++) {
        if (ii == id) continue;
        double t[2];  project_point_to_plane_2D(seeds->p[best_idx[ii]], face[0], u, v, t);
        double dt = signed_distance_point_to_plane(seeds->p[best_idx[ii]], face, TOL);
        double wt = dt*dt - 2*fabs(dt)*plane_offset;
        s_LPconstraint2D con = { .a = 2*(t[0]-s[0]),
                                .b = 2*(t[1]-s[1]),
                                .c = (t[0]*t[0]+t[1]*t[1] + wt)-(s[0]*s[0]+s[1]*s[1] + ws) };
        if (!dynarray_push(C, &con)) goto error;
    }

    /* Add triangle as constraints */
    s_LPconstraint2D triangle_con[3];  constraints_from_triangle(t1, t2, t3, triangle_con);
    if (!dynarray_push(C, &triangle_con[0])) goto error;
    if (!dynarray_push(C, &triangle_con[1])) goto error;
    if (!dynarray_push(C, &triangle_con[2])) goto error;

    /* Check feasability */
    return LP_is_feasible_2D(randint, rctx, C->N, C->items, TOL, NULL, NULL);

    error:
        fprintf(stderr, "should_mirror: Fatal error. (maybe no memory?)\n");
        return 0;
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


int extend_sites_mirroring(const s_bpoly *bp, double EPS_degenerate, double TOL, 
                           s_points *inout_seeds, int (*randint)(void* rctx, int), void* rctx,
                           s_dynarray *buff_points, s_dynarray *buff_3dbls)
{   /* 0 ERROR, 1 OK */
    s_dynarray *mirrored = buff_points;  mirrored->N = 0;

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
            double d = fabs(signed_distance_point_to_plane_v2(inout_seeds->p[pp], normal, d_plane));
            insert_topk(best_idx, best_dist, &best_n, k, pp, d);
        }
        
        for (int ii=0; ii<best_n; ii++) if (should_mirror(face, inout_seeds, best_n, best_idx, ii,
                                                          TOL, randint, rctx, buff_3dbls)) {
            s_point p_mirror = mirror_plane(normal, d_plane, inout_seeds->p[best_idx[ii]]);
            if (!dynarray_push(mirrored, &p_mirror)) goto error;
        }
    }
    // printf("DEBUG: Mirrored.N = %d, total = %d\n", mirrored.N, inout_seeds->N + mirrored.N);


    /* Add mirrored to inout_seeds */
    s_point *tmp = realloc(inout_seeds->p, sizeof(s_point) * (inout_seeds->N + mirrored->N));
    if (!tmp) goto error;
    memcpy(&tmp[inout_seeds->N], mirrored->items, mirrored->N * sizeof(s_point));
    inout_seeds->p = tmp;
    inout_seeds->N += mirrored->N;

    return 1;
    
    error:
        fprintf(stderr, "Error in exten_sites_mirroring!\n");
        return 0;
}


// static int should_mirror_OLD(s_point normal, double d_plane, s_point face[3], const s_points *all_seeds, int seed_id, double EPS_degenerate)
// {
//     s_point s = all_seeds->p[seed_id];
//     double dist = d_plane - dot_prod(normal, s);
//     s_point proj_fplane = {{{s.x + dist*normal.x, s.y + dist*normal.y, s.z + dist*normal.z}}};
//     s_point c = closest_point_on_triangle(face, EPS_degenerate, proj_fplane);
//     if (!point_is_valid(c)) return 0;
//     
//     // Check nearest‐neighbor at c
//     double d_s = distance_squared(c, s);
//     for (int jj=0; jj<all_seeds->N; jj++) {
//         if (jj != seed_id && distance(all_seeds->p[jj], c) + 1e-8 < d_s)
//             return 0;   // someone else is nearer
//     }
//     return 1;
// }

