
#include "bpoly.h"
#include "points.h"
#include "convh.h"
#include "gnuplotc.h"
#include <math.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>


static void extract_dmax_bp(s_bpoly *bpoly)
{
    double dmax = 0; 
    for (int ii=0; ii<bpoly->convh.points.N-1; ii++) {
        for (int jj=ii+1; jj<bpoly->convh.points.N; jj++) {
            double d = distance_squared(bpoly->convh.points.p[ii], 
                                        bpoly->convh.points.p[jj]);
            if (d > dmax) 
                dmax = d;
        }
    }
    bpoly->dmax = sqrt(dmax);
}


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


s_bpoly bpoly_from_points(const s_points *points, double EPS_degenerate, double TOL)
{   // New copy of points inside!
    s_bpoly bpoly;
    if (convhull_from_points(points, EPS_degenerate, TOL, &bpoly.convh) != 1) goto error;
    extract_dmax_bp(&bpoly);
    bpoly.CM = point_average(&bpoly.convh.points);
    extract_min_max_coord(&bpoly);
    bpoly.volume = volume_convhull(&bpoly.convh);
    return bpoly;

    error:
        bpoly.dmax = 0;
        bpoly.CM = (s_point){0};
        bpoly.min = (s_point){0};
        bpoly.max = (s_point){0};
        bpoly.volume = 0;
        return bpoly;
}


s_bpoly bpoly_from_csv(const char *fname, double EPS_degenerate, double TOL)
{
    FILE *f = fopen(fname, "r");
    if (!f) {
        puts("new_bpoly_from_txt: Could not open file");
        exit(1);
    }
    
    s_points points = read_points_from_csv(fname);
    s_bpoly bpoly = bpoly_from_points(&points, EPS_degenerate, TOL);

    free_points(&points);
    return bpoly;
}


s_bpoly bpoly_from_convh(const s_convh *convh)
{
    s_bpoly bpoly;
    bpoly.convh = copy_convhull(convh);
    if (!convhull_is_valid(&bpoly.convh)) goto error;
    extract_dmax_bp(&bpoly);
    bpoly.CM = point_average(&bpoly.convh.points);
    extract_min_max_coord(&bpoly);
    bpoly.volume = volume_convhull(&bpoly.convh);
    return bpoly;

    error:
        bpoly.dmax = 0;
        bpoly.CM = (s_point){0};
        bpoly.min = (s_point){0};
        bpoly.max = (s_point){0};
        bpoly.volume = 0;
        return bpoly;
}



s_bpoly bpoly_copy(const s_bpoly *in)
{
    s_bpoly out;
    out.convh = copy_convhull(&in->convh); 
    out.dmax = in->dmax;
    out.volume = in->volume;
    out.CM = in->CM;
    out.min = in->min;
    out.max = in->max;
    return out;
}


void free_bpoly(s_bpoly *bpoly)
{
    free_convhull(&bpoly->convh);
    memset(bpoly, 0, sizeof(s_bpoly));
}


static void scale_points(s_points *points, double s)
{
    s_point CM = point_average(points);
    s_point b = {{{(1-s)*CM.x, (1-s)*CM.y, (1-s)*CM.z}}};
    for (int ii=0; ii<points->N; ii++) {
        points->p[ii].x = s*points->p[ii].x + b.x;
        points->p[ii].y = s*points->p[ii].y + b.y;
        points->p[ii].z = s*points->p[ii].z + b.z;
    }
}


s_bpoly copy_bpoly_scaled(const s_bpoly *bp, double factor, double EPS_degenerate, double TOL)
{
    s_points new_p = copy_points(&bp->convh.points);
    scale_points(&new_p, factor);
    s_bpoly out = bpoly_from_points(&new_p, EPS_degenerate, TOL);
    free_points(&new_p);
    return out;
}


s_bpoly copy_bpoly_scaled_volume(const s_bpoly *bp, double objective_volume, double EPS_degenerate, double TOL)
{
    double F = objective_volume / (bp)->volume;
    double s = cbrt(F);
    return copy_bpoly_scaled(bp, s, EPS_degenerate, TOL);
}


static int should_mirror(s_point normal, double d_plane, s_point face[3], const s_points *all_seeds, int seed_id, double EPS_degenerate)
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


static int sites_mirroring_CORE(const s_bpoly *bp, const s_points *sites, double EPS_degenerate, s_points *out)
{
    // out should be malloced beforehand with space for Ns+N_mirror. Two-passes!
    int N_mirror = 0;
    for (int ff=0; ff<bp->convh.Nf; ff++) {
        s_point face[3];
        convh_get_face(&bp->convh, ff, face);
        s_point normal = normalize_vec(bp->convh.fnormals[ff], EPS_degenerate);
        if (!point_is_valid(normal)) continue;
        double plane_d = dot_prod(normal, face[0]);
        for (int jj=0; jj<sites->N; jj++) {
            if (!should_mirror(normal, plane_d, face, sites, jj, EPS_degenerate)) continue;

            if (out) {
                s_point sj = sites->p[jj];
                double factor = 2 * (plane_d - dot_prod(normal, sj));
                // if (fabs(factor / 2.0) < 1e-9) continue;
                out->p[sites->N+N_mirror].x = sj.x + factor * normal.x;
                out->p[sites->N+N_mirror].y = sj.y + factor * normal.y;
                out->p[sites->N+N_mirror].z = sj.z + factor * normal.z;
            }
            N_mirror++;
        }
    }
    return N_mirror;
}


void extend_sites_mirroring(const s_bpoly *bp, double EPS_degenerate, s_points *inout)
{
    int N_mirror = sites_mirroring_CORE(bp, inout, EPS_degenerate, NULL);
    inout->p = realloc(inout->p, sizeof(s_point) * (inout->N + N_mirror));
    sites_mirroring_CORE(bp, inout, EPS_degenerate, inout);
    inout->N += N_mirror;  // Need to update AFTER second pass!
    // printf("DEBUG: Reflected %d sites, total sites now: %d\n", N_mirror, Ns+N_mirror);
}


// MY IMPLEMENTATION FOR POISSON DISK SAMPLING WITH WEIGHT FUNCTION

// Generate a random candidate point around a given point p.
// The candidate is generated uniformly in the spherical shell [r, 2r],
// where r = r_of_x(p). We also perturb in all 3 dimensions.
static s_point random_point_around(s_point x, double r)
{
    double radius = r + r * rand()/((double) RAND_MAX + 1);

    // Generate a random direction uniformly over the sphere:
    double aux = 1 - 2.0 * rand()/((double) RAND_MAX + 1);
    if (aux >= 1) aux = 1;
    if (aux <= -1) aux = -1;
    double theta = acos(aux);  // polar angle, 0 <= theta <= pi.
    double phi = 2.0 * M_PI * rand()/((double) RAND_MAX + 1);  // azimuthal, 0 <= phi < 2pi.
    
    s_point out;
    out.x = x.x + radius * sin(theta) * cos(phi);
    out.y = x.y + radius * sin(theta) * sin(phi);
    out.z = x.z + radius * cos(theta);
    return out;
}


static int poisson_is_valid(const s_bpoly *bpoly, s_point query, const s_points *samples, double (*rmax)(double*, void*), void *rmax_params, double EPS_degenerate)
{
    if (test_point_in_convhull(&bpoly->convh, query, EPS_degenerate, 0) != TEST_IN) return 0;

    double rq = rmax(query.coords, rmax_params);
    for (int ii = 0; ii<samples->N; ii++) {
        double rx = rmax(samples->p[ii].coords, rmax_params);
        double minDist = fmin(rq, rx);
        if (distance_squared(query, samples->p[ii]) < (minDist * minDist))
            return 0; // candidate too close to an existing sample
    }
    return 1; // valid candidate
}


s_points generate_poisson_dist_inside(const s_bpoly *bpoly, double (*rmax)(double*, void*), void *rmax_params, double EPS_degenerate)
{
    s_point _samples[MAX_TRIAL_POINTS], _active[MAX_TRIAL_POINTS];
    s_points samples = {.N = 0, .p = _samples};
    s_points active= {.N = 0, .p = _active};

    s_point x = random_point_inside_convhull(&bpoly->convh, EPS_degenerate, bpoly->min, bpoly->max); 
    samples.p[samples.N++] = x;
    active.p[active.N++] = x;

    while (active.N > 0 && samples.N < MAX_TRIAL_POINTS) {
        int random_id = rand() / (RAND_MAX / active.N + 1);
        s_point p = active.p[random_id];
        int found = 0;

        double rp = rmax(p.coords, rmax_params);
        for (int ii=0; ii<MAX_TRIAL_TESTS; ii++) {
            s_point q = random_point_around(p, rp);
            if (poisson_is_valid(bpoly, q, &samples, rmax, rmax_params, EPS_degenerate)) {
                samples.p[samples.N++] = q;
                active.p[active.N++] = q;
                found = 1;
                break;
            }
        }

        if (samples.N >= MAX_TRIAL_POINTS-1) {
            puts("ERROR! Max_trial_points in poisson sampling.\n");
            exit(1);
        }
        
        if (found == 0) {
            // Replace activeList[idx] with the last active point and decrease count.
            active.p[random_id] = active.p[active.N--];
        }
    }

    while (samples.N < 4) {
        x = random_point_inside_convhull(&bpoly->convh, EPS_degenerate, bpoly->min, bpoly->max);
        samples.p[samples.N++] = x;
    }

    s_points out = copy_points(&samples);
    return out;
}


double find_closest_point_on_bp(const s_bpoly *bp, s_point p, double EPS_degenerate, s_point *out)
{
    *out = closest_point_on_convhull_boundary(&bp->convh, p, EPS_degenerate);
    return distance(p, *out);
}


void generate_file_cube_bp(const char *filename, double length)
{
    double s = length / 2;
    
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "%f, %f, %f\n", -s, -s, -s);
    fprintf(fp, "%f, %f, %f\n", -s, -s, s);
    fprintf(fp, "%f, %f, %f\n", -s, s, -s);
    fprintf(fp, "%f, %f, %f\n", s, -s, -s);
    fprintf(fp, "%f, %f, %f\n", -s, s, s);
    fprintf(fp, "%f, %f, %f\n", s, -s, s);
    fprintf(fp, "%f, %f, %f\n", s, s, -s);
    fprintf(fp, "%f, %f, %f\n", s, s, s);
    fclose(fp);
}


void generate_file_tetrahedron_bp(const char *filename, double length)
{
    double s = length / 2;
    
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "%f, %f, %f\n", -s, -s, -s);
    fprintf(fp, "%f, %f, %f\n", -s, -s, s);
    fprintf(fp, "%f, %f, %f\n", -s, s, -s);
    fprintf(fp, "%f, %f, %f\n", s, s, s);
    fclose(fp);
}


void generate_file_sphere_bp(const char *filename, double radius, int nTheta, int nPhi)
{
    // ntheta: example 18; // Number of steps in the polar angle
    // nphi: example 36;   // Number of steps in the azimuthal angle
    FILE *fp = fopen(filename, "w");
    // fprintf(fp, "%d\n\n", 2 + nPhi * (nTheta-1));

    for (int i = 0; i <= nTheta; i++) {
        double theta = M_PI * i / nTheta;
        
        if (i == 0 || i == nTheta) { // If at a pole, compute and write the coordinate once.
            double x = 0.0;
            double y = 0.0;
            double z = radius * cos(theta);  // will be +radius or -radius
            fprintf(fp, "%f, %f, %f\n", x, y, z);
        } else {
            for (int j = 0; j < nPhi; j++) {
                double phi = 2 * M_PI * j / nPhi;
                double x = radius * sin(theta) * cos(phi);
                double y = radius * sin(theta) * sin(phi);
                double z = radius * cos(theta);
                fprintf(fp, "%f, %f, %f\n", x, y, z);
            }
        }
    }
    fclose(fp);
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


