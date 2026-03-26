
#include "points.h"
#include "convh.h"
#include "gnuplotc.h"
#include "vdiagram.h"
#include <math.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>


void free_bpoly(s_bpoly *bpoly)
{
    free_convhull(&bpoly->convh);
    memset(bpoly, 0, sizeof(s_bpoly));
}


double find_closest_point_on_bp(const s_bpoly *bp, s_point p, double EPS_DEG, s_point *out)
{
    *out = closest_point_on_convhull_boundary(&bp->convh, p, EPS_DEG);
    return distance(p, *out);
}


// static void extract_dmax_bp(s_bpoly *bpoly)
// {
//     double dmax = 0; 
//     for (int ii=0; ii<bpoly->convh.points.N-1; ii++) {
//         for (int jj=ii+1; jj<bpoly->convh.points.N; jj++) {
//             double d = distance_squared(bpoly->convh.points.p[ii], 
//                                         bpoly->convh.points.p[jj]);
//             if (d > dmax) 
//                 dmax = d;
//         }
//     }
//     bpoly->dmax = sqrt(dmax);
// }


static void extract_min_max_coord(s_bpoly *bpoly)
{
    s_point min = bpoly->convh.points.p[0], max = bpoly->convh.points.p[0];

    for (int ii=1; ii<bpoly->convh.points.N; ii++) {
        if (bpoly->convh.points.p[ii].x < min.x) 
            min.x = bpoly->convh.points.p[ii].x;
        if (bpoly->convh.points.p[ii].y < min.y) 
            min.y = bpoly->convh.points.p[ii].y;
        if (bpoly->convh.points.p[ii].z < min.z) 
            min.z = bpoly->convh.points.p[ii].z;

        if (bpoly->convh.points.p[ii].x > max.x)
            max.x = bpoly->convh.points.p[ii].x;
        if (bpoly->convh.points.p[ii].y > max.y) 
            max.y = bpoly->convh.points.p[ii].y;
        if (bpoly->convh.points.p[ii].z > max.z) 
            max.z = bpoly->convh.points.p[ii].z;
    }

    bpoly->min = min;
    bpoly->max = max;
}


s_bpoly bpoly_from_points(const s_points *points, double EPS_DEG, double TOL)
{   /* New copy of points inside! */
    s_bpoly bpoly;
    if (convhull_from_points(points, EPS_DEG, TOL, &bpoly.convh) != 1) goto error;
    // extract_dmax_bp(&bpoly);
    bpoly.CM = point_average(&bpoly.convh.points);
    extract_min_max_coord(&bpoly);
    bpoly.volume = volume_convhull(&bpoly.convh);
    return bpoly;

    error:
        memset(&bpoly, 0, sizeof(s_bpoly));
        return bpoly;
}


s_bpoly bpoly_from_csv(const char *fname, double EPS_DEG, double TOL)
{
    s_points points = read_points_from_csv(fname);
    s_bpoly bpoly = bpoly_from_points(&points, EPS_DEG, TOL);

    free_points(&points);
    return bpoly;
}


s_bpoly bpoly_from_convh(const s_convh *convh)
{
    s_bpoly bpoly;
    bpoly.convh = copy_convhull(convh);
    if (!convhull_is_valid(&bpoly.convh)) goto error;
    // extract_dmax_bp(&bpoly);
    bpoly.CM = point_average(&bpoly.convh.points);
    extract_min_max_coord(&bpoly);
    bpoly.volume = volume_convhull(&bpoly.convh);
    return bpoly;

    error:
        fprintf(stderr, "Error bpoly_from_convh\n");
        memset(&bpoly, 0, sizeof(s_bpoly));
        return bpoly;
}


s_bpoly bpoly_from_convh_scaled(const s_convh *convh, double s, s_point pivot, double EPS_DEG, double TOL)
{
    s_points new_p = copy_points(&convh->points);
    homotethy_points(&new_p, s, pivot);

    s_bpoly out = bpoly_from_points(&new_p, EPS_DEG, TOL);
    free_points(&new_p);
    return out;
}



s_bpoly bpoly_copy(const s_bpoly *in)
{
    s_bpoly out;
    out.convh = copy_convhull(&in->convh); 
    // out.dmax = in->dmax;
    out.volume = in->volume;
    out.CM = in->CM;
    out.min = in->min;
    out.max = in->max;
    return out;
}

s_bpoly bpoly_copy_scaled(const s_bpoly *bp, double factor, s_point pivot, double EPS_DEG, double TOL)
{
    s_points new_p = copy_points(&bp->convh.points);
    homotethy_points(&new_p, factor, pivot);
    s_bpoly out = bpoly_from_points(&new_p, EPS_DEG, TOL);
    free_points(&new_p);
    return out;
}


s_bpoly bpoly_copy_scaled_volume(const s_bpoly *bp, double objective_volume, double EPS_DEG, double TOL)
{
    double F = objective_volume / (bp)->volume;
    double s = cbrt(F);
    return bpoly_copy_scaled(bp, s, bp->CM, EPS_DEG, TOL);
}



/* IMPLEMENTATION OF POISSON DISK SAMPLING WITH WEIGHT FUNCTION
 * Generate a random candidate point around a given point p.
 * The candidate is generated uniformly in the spherical shell [r, 2r],
 * where r = r_of_x(p). */
static s_point random_point_around(double (*randd01)(void *rctx), void *rctx,
                                   s_point x, double r)
{
    double radius = r + r * randd01(rctx);

    // Generate a random direction uniformly over the sphere:
    double aux = 1 - 2.0 * randd01(rctx);
    if (aux >= 1) aux = 1; if (aux <= -1) aux = -1;
    double theta = acos(aux);  // polar angle, 0 <= theta <= pi.
    double phi = 2.0 * M_PI * randd01(rctx);  // azimuthal, 0 <= phi < 2pi.
    
    s_point out;
    out.x = x.x + radius * sin(theta) * cos(phi);
    out.y = x.y + radius * sin(theta) * sin(phi);
    out.z = x.z + radius * cos(theta);
    return out;
}


static int poisson_is_valid(const s_bpoly *bpoly, s_point query, const s_points *samples,
                            double (*rmax)(double*, void*), void *rmax_params, double EPS_DEG)
{
    if (test_point_in_convhull(&bpoly->convh, query, EPS_DEG, 0) != TEST_IN) return 0;

    double rq = rmax(query.coords, rmax_params);
    for (int ii = 0; ii<samples->N; ii++) {
        double rx = rmax(samples->p[ii].coords, rmax_params);
        double minDist = fmin(rq, rx);
        if (distance_squared(query, samples->p[ii]) < (minDist * minDist))
            return 0;  // candidate too close to an existing sample
    }
    return 1;  // valid candidate
}


#define MAX_TRIAL_POINTS 10000
#define MAX_TRIAL_TESTS 50

s_points generate_poisson_dist_inside(const s_bpoly *bpoly, 
                                      double (*rmax)(double*, void*), void *rmax_params,
                                      double (*randd01)(void*), int (*randint)(void*, int),
                                      void *rctx, double EPS_DEG)
{
    s_point _samples[MAX_TRIAL_POINTS], _active[MAX_TRIAL_POINTS];
    s_points samples = {.N = 0, .p = _samples};
    s_points active= {.N = 0, .p = _active};

    s_point x = random_point_inside_convhull(randd01, rctx, &bpoly->convh, EPS_DEG, bpoly->min, bpoly->max); 
    samples.p[samples.N++] = x;
    active.p[active.N++] = x;

    while (active.N > 0 && samples.N < MAX_TRIAL_POINTS) {
        int random_id = randint(rctx, active.N);
        s_point p = active.p[random_id];
        int found = 0;

        double rp = rmax(p.coords, rmax_params);
        for (int ii=0; ii<MAX_TRIAL_TESTS; ii++) {
            s_point q = random_point_around(randd01, rctx, p, rp);
            if (poisson_is_valid(bpoly, q, &samples, rmax, rmax_params, EPS_DEG)) {
                samples.p[samples.N++] = q;
                active.p[active.N++] = q;
                found = 1;
                break;
            }
        }

        if (samples.N >= MAX_TRIAL_POINTS-1) {
            fprintf(stderr, "ERROR! Max_trial_points in poisson sampling.\n");
            return points_NAN;
        }

        // Replace activeList[idx] with the last active point and decrease count.
        if (found == 0) active.p[random_id] = active.p[--active.N];
    }

    s_points out = copy_points(&samples);
    return out;
}






void plot_bpoly(s_bpoly *bpoly, char *f_name, s_point ranges[2], char *color, char *view_command)
{
    s_gnuplot *interface = gnuplot_start(PNG_3D, f_name, (int[2]){1080, 1080}, 18);
    gnuplot_config(interface, "set pm3d depthorder", 
                              "set pm3d border lc 'black' lw 0.1",
                              "set xyplane at 0",
                              "unset border",
                              "unset xtics",
                              "unset ytics",
                              "unset ztics",
                              view_command);
    if (ranges) {
        char buff[1024];
        snprintf(buff, 1024, "set xrange [%f:%f]\n set yrange [%f:%f]\n set zrange [%f:%f]", 
                 ranges[0].x, ranges[1].x, ranges[0].x, ranges[1].y, ranges[0].z, ranges[1].z); 
        gnuplot_config(interface, buff);
    }

    for (int ii=0; ii<bpoly->convh.Nf; ii++) {
        char buff[256];
        if (color) snprintf(buff, 256, "fs transparent solid 0.2 fc '%s'", color);
        else snprintf(buff, 256, "fs transparent solid 0.2 fc 'blue'");

        s_point face[3];
        convh_get_face(&bpoly->convh, ii, face);
        draw_solid_triangle_3d(interface, face[0].coords, face[1].coords, face[2].coords, buff);
    }
    gnuplot_end(interface);
}


void plot_bpoly_differentviews(s_bpoly *bpoly, char *f_name, s_point ranges[2], char *color)
{   
    char real_name[256];
    snprintf(real_name, 256, "%s_v1.png", f_name);
    plot_bpoly(bpoly, real_name, ranges, color, "set view 100, 60, 1.5");

    snprintf(real_name, 256, "%s_v2.png", f_name);
    plot_bpoly(bpoly, real_name, ranges, color, "set view 100, 90, 1.5");

    snprintf(real_name, 256, "%s_v3.png", f_name);
    plot_bpoly(bpoly, real_name, ranges, color, "set view 100, 180, 1.5");

    snprintf(real_name, 256, "%s_v4.png", f_name);
    plot_bpoly(bpoly, real_name, ranges, color, "set view 100, 270, 1.5");
}


