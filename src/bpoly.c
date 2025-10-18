
#include "bpoly.h"
#include "geometry.h"
#include "convh.h"
#include "float.h"
#include "math.h"
#include "gnuplotc.h"
#include <assert.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>


void extract_dmax_bp(s_bpoly *bpoly)
{
    double dmax = 0; 
    for (int ii=0; ii<bpoly->Np-1; ii++) {
        for (int jj=ii+1; jj<bpoly->Np; jj++) {
            double d = distance_squared(bpoly->points[ii], bpoly->points[jj]);
            if (d > dmax) 
                dmax = d;
        }
    }
    bpoly->dmax = sqrt(dmax);
}


void extract_CM_bp(s_bpoly *bpoly)
{
    bpoly->CM = find_center_mass(bpoly->points, bpoly->Np);
}


void extract_min_max_coord(const s_bpoly *bpoly, s_point *min, s_point *max)
{
    min->x = DBL_MAX;   min->y = DBL_MAX;   min->z = DBL_MAX;
    max->x = -DBL_MAX;  max->y = -DBL_MAX;  max->z = -DBL_MAX;
    for (int ii=0; ii<bpoly->Np; ii++) {
        if (bpoly->points[ii].x < min->x) min->x = bpoly->points[ii].x;
        if (bpoly->points[ii].y < min->y) min->y = bpoly->points[ii].y;
        if (bpoly->points[ii].z < min->z) min->z = bpoly->points[ii].z;

        if (bpoly->points[ii].x > max->x) max->x = bpoly->points[ii].x;
        if (bpoly->points[ii].y > max->y) max->y = bpoly->points[ii].y;
        if (bpoly->points[ii].z > max->z) max->z = bpoly->points[ii].z;
    }
}


void extract_convhull_bp(s_bpoly *bp)
{   
    convhull_from_points(bp->points, bp->Np, &bp->faces, &bp->fnormals, &bp->Nf);
}


double compute_volume_bpoly(s_bpoly *bpoly)
{
    double vol = 0;
    for (int ii=0; ii<bpoly->Nf; ii++) {
        int i0 = 0;
        int i1 = 1;
        int i2 = 2;

        // Normals should be unnormalized!
        double Nx = (bpoly->points[bpoly->faces[ii*3 + 1]].coords[i1] -
                     bpoly->points[bpoly->faces[ii*3 + 0]].coords[i1]) *
                    (bpoly->points[bpoly->faces[ii*3 + 2]].coords[i2] -
                     bpoly->points[bpoly->faces[ii*3 + 0]].coords[i2]) 
                    -
                    (bpoly->points[bpoly->faces[ii*3 + 1]].coords[i2] -
                     bpoly->points[bpoly->faces[ii*3 + 0]].coords[i2]) *
                    (bpoly->points[bpoly->faces[ii*3 + 2]].coords[i1] -
                     bpoly->points[bpoly->faces[ii*3 + 0]].coords[i1]);

        vol += Nx * (bpoly->points[bpoly->faces[ii*3 + 0]].coords[i0] +
                     bpoly->points[bpoly->faces[ii*3 + 1]].coords[i0] +
                     bpoly->points[bpoly->faces[ii*3 + 2]].coords[i0]);
    }
    return vol/6;
}


s_bpoly *new_bpoly_from_points(const s_point *points, double Np)
{   // New copy of points inside!
    s_bpoly *bpoly = malloc(sizeof(s_bpoly));
    bpoly->points = malloc(sizeof(s_point) * Np);
    memcpy(bpoly->points, points, Np * sizeof(s_point));
    bpoly->Np = Np;
    
    extract_dmax_bp(bpoly);
    extract_convhull_bp(bpoly);
    extract_CM_bp(bpoly);
    extract_min_max_coord(bpoly, &bpoly->min, &bpoly->max);

    bpoly->volume = compute_volume_bpoly(bpoly);
    
    return bpoly;
}


s_bpoly *new_bpoly_from_txt(const char *fname)
{
    FILE *f = fopen(fname, "r");
    if (!f) {
        puts("new_bpoly_from_txt: Could not open file");
        exit(1);
    }
    
    int Np;
    fscanf(f, "%d\n\n", &Np);

    s_point *points = malloc(Np * sizeof(s_point));
    for (int ii=0; ii<Np; ii++) {
        fscanf(f, "%lf, %lf, %lf\n", &points[ii].x, &points[ii].y, &points[ii].z); 
    }

     s_bpoly *out = new_bpoly_from_points(points, Np);
     free(points);
     return out;
}


s_bpoly *new_bpoly_copy(s_bpoly *in)
{
    s_bpoly *out = malloc(sizeof(s_bpoly));

    out->Np = in->Np;
    out->points = malloc(in->Np * sizeof(s_point));
    memcpy(out->points, in->points, sizeof(s_point) * in->Np);

    out->Nf = in->Nf;
    out->faces = malloc(sizeof(int) * in->Nf * 3);
    memcpy(out->faces, in->faces, sizeof(int) * 3 * in->Nf);

    out->fnormals = malloc(sizeof(s_point) * in->Nf);
    memcpy(out->fnormals, in->fnormals, sizeof(s_point) * in->Nf);

    out->dmax = in->dmax;
    out->volume = in->volume;
    out->CM = in->CM;
    out->min = in->min;
    out->max = in->max;

    return out;
}


void free_bpoly(s_bpoly *bpoly)
{
    free(bpoly->points);
    free(bpoly->faces);
    free(bpoly->fnormals);
    free(bpoly);
}


void scale_points(s_point *points, int Np, double s)
{
    s_point CM = find_center_mass(points, Np);
    s_point b = {{{(1-s)*CM.x, (1-s)*CM.y, (1-s)*CM.z}}};
    for (int ii=0; ii<Np; ii++) {
        points[ii].x = s*points[ii].x + b.x;
        points[ii].y = s*points[ii].y + b.y;
        points[ii].z = s*points[ii].z + b.z;
    }
}


s_bpoly *copy_bpoly_scaled(const s_bpoly *bp, double factor)
{
    // double OLD_VOL = bp->volume;
    s_point *new_p = malloc(sizeof(s_point) * bp->Np);
    memcpy(new_p, bp->points, sizeof(s_point) * bp->Np);
    scale_points(new_p, bp->Np, factor);
    s_bpoly *out = new_bpoly_from_points(new_p, bp->Np);
    free(new_p);
    return out;
}


s_bpoly *copy_bpoly_scaled_volume(const s_bpoly *bp, double objective_volume)
{
    assert(fabs(bp->volume) > 1e-6);
    double F = objective_volume / (bp)->volume;
    double s = cbrt(F);
    return copy_bpoly_scaled(bp, s);
}


int should_mirror(s_point normal, double d_plane, s_point face[3], const s_point *all_seeds, int Ns, int seed_id)
{
    s_point s = all_seeds[seed_id];
    double dist = d_plane - dot_prod(normal, s);
    s_point proj_fplane = {{{s.x + dist*normal.x, s.y + dist*normal.y, s.z + dist*normal.z}}};
    s_point c = closest_point_on_triangle(face, proj_fplane);
    
    // Check nearest‐neighbor at c
    double d_s = distance_squared(c, s);
    for (int j = 0; j < Ns; j++) {
        if (j != seed_id && distance(all_seeds[j], c) + 1e-6 < d_s)
            return 0;   // someone else is nearer
    }
    return 1;
}


int sites_mirroring_CORE(const s_bpoly *bp, const s_point *s, int Ns, s_point *out)
{
    // out should be malloced beforehand with space for Ns+N_mirror. Two-passes!
    int N_mirror = 0;
    for (int ff=0; ff<bp->Nf; ff++) {
        s_point face[3] = {bp->points[bp->faces[ff*3+0]],
                           bp->points[bp->faces[ff*3+1]],
                           bp->points[bp->faces[ff*3+2]]};
        s_point normal = bp->fnormals[ff];
        double plane_d = dot_prod(normal, face[0]);
        for (int jj=0; jj<Ns; jj++) {
            if (!should_mirror(normal, plane_d, face, s, Ns, jj)) continue;

            if (out) {
                s_point sj = s[jj];
                double factor = 2 * (plane_d - dot_prod(normal, sj));
                // if (fabs(factor / 2.0) < 1e-9) continue;
                out[Ns+N_mirror].x = sj.x + factor * normal.x;
                out[Ns+N_mirror].y = sj.y + factor * normal.y;
                out[Ns+N_mirror].z = sj.z + factor * normal.z;
            }
            N_mirror++;
        }
    }
    return N_mirror;
}


int extend_sites_mirroring(const s_bpoly *bp, s_point **s, int Ns)
{
    int N_mirror = sites_mirroring_CORE(bp, *s, Ns, NULL);
    *s = realloc(*s, sizeof(s_point) * (Ns+N_mirror));
    sites_mirroring_CORE(bp, *s, Ns, *s);
    // printf("DEBUG: Reflected %d sites, total sites now: %d\n", N_mirror, Ns+N_mirror);
    return Ns+N_mirror;
}



// MY IMPLEMENTATION FOR POISSON DISK SAMPLING WITH WEIGHT FUNCTION

// Generate a random candidate point around a given point p.
// The candidate is generated uniformly in the spherical shell [r, 2r],
// where r = r_of_x(p). We also perturb in all 3 dimensions.
s_point random_point_around(s_point x, double r)
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


int poisson_is_valid(const s_bpoly *bpoly, s_point query, s_point *samples, int Ns, double (*rmax)(double *))
{
    if (is_inside_convhull(query, bpoly->points, bpoly->faces,bpoly->Nf) != 1) 
        return 0;

    double rq = rmax(query.coords);
    for (int ii = 0; ii<Ns; ii++) {
        double rx = rmax(samples[ii].coords);
        double minDist = fmin(rq, rx);
        if (distance_squared(query, samples[ii]) < (minDist * minDist))
            return 0; // candidate too close to an existing sample
    }
    return 1; // valid candidate
}


s_point *generate_poisson_dist_inside(const s_bpoly *bpoly, double (*rmax)(double *), int *Np_generated)
{
    s_point samples[MAX_TRIAL_POINTS], active_list[MAX_TRIAL_POINTS];
    int Nsamples = 0, Nactive = 0;

    s_point x = random_point_uniform_3d(bpoly->min, bpoly->max);
    while (is_inside_convhull(x, bpoly->points, bpoly->faces, bpoly->Nf) != 1) {
        x = random_point_uniform_3d(bpoly->min, bpoly->max);
    }

    samples[Nsamples++] = x;
    active_list[Nactive++] = x;

    while (Nactive > 0 && Nsamples < MAX_TRIAL_POINTS) {
        int random_id = rand() % Nactive;
        s_point p = active_list[random_id];
        int found = 0;

        double rp = rmax(p.coords);
        for (int ii=0; ii<MAX_TRIAL_TESTS; ii++) {
            s_point q = random_point_around(p, rp);
            if (poisson_is_valid(bpoly, q, samples, Nsamples, rmax)) {
                samples[Nsamples++] = q;
                active_list[Nactive++] = q;
                found = 1;
                break;
            }
        }

        if (Nsamples >= MAX_TRIAL_POINTS-1) {
            puts("ERROR! Max_trial_points in poisson sampling.\n");
            exit(1);
        }
        
        if (found == 0) {
            // Replace activeList[idx] with the last active point and decrease count.
            active_list[random_id] = active_list[Nactive - 1];
            Nactive--;
        }
    }

    while (Nsamples < 4) {
        x = random_point_uniform_3d(bpoly->min, bpoly->max);
        while (is_inside_convhull(x, bpoly->points, bpoly->faces, bpoly->Nf) != 1) {
            x = random_point_uniform_3d(bpoly->min, bpoly->max);
        }
        samples[Nsamples++] = x;
    }

    s_point *out = malloc(sizeof(s_point) * Nsamples);
    for (int ii=0; ii<Nsamples; ii++) {
        out[ii] = samples[ii];
    }
    
    *Np_generated = Nsamples;
    return out;
}


s_point find_closest_point_on_bp(const s_bpoly *bp, s_point p)
{
    assert(bp->Nf != 0 && bp->Np != 0);
    s_point out;
    double best_d2 = DBL_MAX;

    for (int f=0; f<bp->Nf; f++) {
        s_point triangle[3] = {bp->points[bp->faces[3*f+0]],
                               bp->points[bp->faces[3*f+1]],
                               bp->points[bp->faces[3*f+2]]};
        s_point tmp = closest_point_on_triangle(triangle, p);
        double d2 = distance_squared(p, tmp);
        if ( d2 < best_d2) {
            out = tmp;
            best_d2 = d2; 
        }
    }
    return out;
}


void generate_file_cube_bp(const char *filename, double length)
{
    double s = length / 2;
    
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "%d\n\n", 8);
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
    fprintf(fp, "%d\n\n", 4);
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
    fprintf(fp, "%d\n\n", 2 + nPhi * (nTheta-1));

    for (int i = 0; i <= nTheta; i++) {
        double theta = M_PI * i / nTheta;
        
        if (i == 0 || i == nTheta) { // If at a pole, compute and write the coordinate once.
            double x = 0.0;
            double y = 0.0;
            double z = radius * cos(theta);  // will be +radius or -radius
            fprintf(fp, "%f %f %f\n", x, y, z);
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

    for (int ii=0; ii<bpoly->Nf; ii++) {
        char buff[256];
        if (color) snprintf(buff, 256, "fs transparent solid 0.2 fc '%s'", color);
        else snprintf(buff, 256, "fs transparent solid 0.2 fc 'blue'");

        draw_solid_triangle_3d(interface, bpoly->points[bpoly->faces[ii*3]].coords, 
                bpoly->points[bpoly->faces[ii*3+1]].coords, 
                bpoly->points[bpoly->faces[ii*3+2]].coords, buff);
    }
    gnuplot_end(interface);
}


void plot_bpoly_differentviews(s_bpoly *bpoly, char *f_name, s_point ranges[2], char *color)
{   
    char real_name[256];
    snprintf(real_name, 256, "%s_v1.png", f_name);
    plot_bpoly(bpoly, f_name, ranges, color, "set view 100, 60, 1.5");

    snprintf(real_name, 256, "%s_v2.png", f_name);
    plot_bpoly(bpoly, f_name, ranges, color, "set view 100, 90, 1.5");

    snprintf(real_name, 256, "%s_v3.png", f_name);
    plot_bpoly(bpoly, f_name, ranges, color, "set view 100, 180, 1.5");

    snprintf(real_name, 256, "%s_v4.png", f_name);
    plot_bpoly(bpoly, f_name, ranges, color, "set view 100, 270, 1.5");
}


