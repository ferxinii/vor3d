
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


s_bpoly bpoly_from_points(const s_points *points, double EPS_DEG)
{   /* New copy of points inside! */
    s_bpoly bpoly;
    if (convhull_from_points(points, EPS_DEG, &bpoly.convh) != 1) goto error;
    // extract_dmax_bp(&bpoly);
    bpoly.CM = point_average(&bpoly.convh.points);
    extract_min_max_coord(&bpoly);
    bpoly.volume = volume_convhull(&bpoly.convh);
    return bpoly;

    error:
        memset(&bpoly, 0, sizeof(s_bpoly));
        return bpoly;
}


s_bpoly bpoly_from_csv(const char *fname, double EPS_DEG)
{
    s_points points = read_points_from_csv(fname);
    s_bpoly bpoly = bpoly_from_points(&points, EPS_DEG);

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


s_bpoly bpoly_from_convh_scaled(const s_convh *convh, double s, s_point pivot, double EPS_DEG)
{
    s_points new_p = copy_points(&convh->points);
    homotethy_points(&new_p, s, pivot);

    s_bpoly out = bpoly_from_points(&new_p, EPS_DEG);
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

s_bpoly bpoly_copy_scaled(const s_bpoly *bp, double factor, s_point pivot, double EPS_DEG)
{
    s_points new_p = copy_points(&bp->convh.points);
    homotethy_points(&new_p, factor, pivot);
    s_bpoly out = bpoly_from_points(&new_p, EPS_DEG);
    free_points(&new_p);
    return out;
}


s_bpoly bpoly_copy_scaled_volume(const s_bpoly *bp, double objective_volume, double EPS_DEG)
{
    double F = objective_volume / (bp)->volume;
    double s = cbrt(F);
    return bpoly_copy_scaled(bp, s, bp->CM, EPS_DEG);
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


