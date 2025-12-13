
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
#include <stdbool.h>


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


/* Mirroring seeds */
typedef struct int_list {
    int *list;
    int N;
    int Nmax;
} s_int_list;

static s_int_list initialize_int_list(int Nmax) {
    if (Nmax <= 0) Nmax = 50;
    s_int_list out = { .Nmax = Nmax, 
                       .list = malloc(Nmax * sizeof(int)),
                       .N = 0 };
    return out;
}

static int increase_memory_int_list_if_needed(s_int_list *int_list, int N_needed) {
    while (N_needed >= int_list->Nmax) {
        int *tmp = realloc(int_list->list, 2 * int_list->Nmax * sizeof(int));
        if (!tmp) {
            return 0;
        }
        int_list->list = tmp;
        int_list->Nmax *= 2;
    }
    return 1;
}

static void free_int_list(s_int_list *int_list) {
    free(int_list->list);
    memset(int_list, 0, sizeof(s_int_list));
}

typedef struct point_list {
    s_point *list;
    int N;
    int Nmax;
} s_point_list;

static s_point_list initialize_point_list(int Nmax) {
    if (Nmax <= 0) Nmax = 50;
    s_point_list out = { .Nmax = Nmax, 
                         .list = malloc(Nmax * sizeof(s_point)),
                         .N = 0 };
    return out;
}

static int increase_memory_point_list_if_needed(s_point_list *point_list, int N_needed) {
    while (N_needed >= point_list->Nmax) {
        s_point *tmp = realloc(point_list->list, 2 * point_list->Nmax * sizeof(s_point));
        if (!tmp) return 0;
        point_list->list = tmp;
        point_list->Nmax *= 2;
    }
    return 1;
}

static void free_point_list(s_point_list *point_list) {
    free(point_list->list);
    memset(point_list, 0, sizeof(s_point_list));
}



static s_point mirror_plane(s_point normal, double d_plane, s_point p)
{
    double factor = 2 * (d_plane - dot_prod(normal, p));
    s_point out;
    out.x = p.x + factor * normal.x;
    out.y = p.y + factor * normal.y;
    out.z = p.z + factor * normal.z;
    return out;
}


// static int witness_points_from_face(const s_points *seeds, const s_point face[3], double EPS_degenerate, s_point out_witness[seeds->N+3])
// {   /* 0 ERROR, 1 OK */
//     int N = 3;
//     out_witness[0] = face[0];  out_witness[1] = face[1];  out_witness[2] = face[2];
//
//     for (int ii=0; ii<seeds->N; ii++) {
//         s_point closest = closest_point_on_triangle(face, EPS_degenerate, seeds->p[ii]);
//         if (!point_is_valid(closest)) return 0;
//         out_witness[N++] = closest;
//     }
//     return 1;
// }
//
// static int mirror_seeds_face(const s_points *seeds, const s_point face[3], s_point normal, double d_plane, double EPS_degenerate, double TOL2, bool buff_mask[seeds->N], s_point buff_witness[seeds->N+3], s_int_list *buff_closest, s_point_list *out_mirrored)
// {   /* 0 ERROR, 1 OK */
//     memset(buff_mask, 0, sizeof(bool) * seeds->N);
//     
//     if (!witness_points_from_face(seeds, face, EPS_degenerate, buff_witness)) return 0;
//     for (int ii=0; ii<seeds->N+3; ii++) {  /* Loop over witness: N+3 */
//         /* Find those seeds which are closest to the witness, might be more than 1 */
//         double dmin2 = distance_squared(buff_witness[ii], seeds->p[0]);
//         buff_closest->N = 1;
//         buff_closest->list[0] = 0;
//         for (int jj=1; jj<seeds->N; jj++) {
//             double d2 = distance_squared(buff_witness[ii], seeds->p[jj]);
//             if (fabs(d2-dmin2) < TOL2) {  /* Same distance than previous min */
//                 if (!increase_memory_int_list_if_needed(buff_closest, buff_closest->N+1)) return 0;
//                 buff_closest->list[buff_closest->N++] = jj;
//                 printf("DEBUG: TIE!\n");
//             } else if (d2 < dmin2) {  /* Found new minimum */
//                 dmin2 = d2;
//                 buff_closest->N = 1;
//                 buff_closest->list[0] = jj;
//             }
//             
//         }
//         /* mirror if not already mirrored */
//         for (int jj=0; jj<buff_closest->N; jj++) {
//             if (!buff_mask[buff_closest->list[jj]]) {
//                 if (!increase_memory_point_list_if_needed(out_mirrored, out_mirrored->N+1)) return 0;
//                 out_mirrored->list[out_mirrored->N++] = mirror_plane(normal, d_plane, seeds->p[buff_closest->list[jj]]);
//                 buff_mask[buff_closest->list[jj]] = true;
//             }
//         }
//     }
//     return 1;
// }

int extend_sites_mirroring(const s_bpoly *bp, double EPS_degenerate, double TOL, s_points *inout_seeds)
{   /* 0 ERROR, 1 OK */
    s_point_list mirrored = initialize_point_list(0);
    s_int_list buff_closest = initialize_int_list(0);
    s_point *buff_witness = malloc(sizeof(s_point) * (inout_seeds->N + 3)); 
    bool *buff_mask_mirrored = calloc(inout_seeds->N, sizeof(bool));
    if (!mirrored.list || !buff_witness || !buff_mask_mirrored) return 0;
    
    /* Compute all mirrored seeds */
    const double TOL2 = TOL * TOL;
    for (int ii=0; ii<bp->convh.Nf; ii++) {
        s_point face[3];
        convh_get_face(&bp->convh, ii, face);
        s_point normal = normalize_vec(bp->convh.fnormals[ii], EPS_degenerate);
        if (!point_is_valid(normal)) goto error;
        double d_plane =  dot_prod(normal, face[0]);
        if (!mirror_seeds_face(inout_seeds, face, normal, d_plane, EPS_degenerate, TOL2, buff_mask_mirrored, buff_witness, &buff_closest, &mirrored)) goto error;
    }
    printf("DEBUG: mirrored: %d\n", mirrored.N);

    /* Add them to inout_seeds */
    s_point *tmp = realloc(inout_seeds->p, sizeof(s_point) * (inout_seeds->N + mirrored.N));
    if (!tmp) goto error;
    for (int ii=0; ii< mirrored.N; ii++) 
        tmp[inout_seeds->N + ii] = mirrored.list[ii];

    inout_seeds->p = tmp;
    inout_seeds->N += mirrored.N;
    printf("DEBUG: mirrored->N = %d, total = %d\n", mirrored.N, inout_seeds->N);

    free_point_list(&mirrored);
    free_int_list(&buff_closest);
    free(buff_witness);
    free(buff_mask_mirrored);
    return 1;

    error:
        free_point_list(&mirrored);
        free_int_list(&buff_closest);
        free(buff_witness);
        free(buff_mask_mirrored);
        return 0;
}


// static int should_mirror(s_point normal, double d_plane, const s_point face[3], const s_points *all_seeds, int seed_id, double EPS_degenerate, double TOL)
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
//         if (jj != seed_id && distance(all_seeds->p[jj], c) + TOL < d_s)
//             return 0;   // someone else is nearer
//     }
//     return 1;
// }
//
//
// static int sites_mirroring_CORE(const s_bpoly *bp, const s_points *sites, double EPS_degenerate, s_points *out)
// {
//     // out should be malloced beforehand with space for Ns+N_mirror. Two-passes!
//     int N_mirror = 0;
//     for (int ff=0; ff<bp->convh.Nf; ff++) {
//         s_point face[3];
//         convh_get_face(&bp->convh, ff, face);
//         s_point normal = normalize_vec(bp->convh.fnormals[ff], EPS_degenerate);
//         if (!point_is_valid(normal)) continue;
//         double plane_d = dot_prod(normal, face[0]);
//         for (int jj=0; jj<sites->N; jj++) {
//             if (!should_mirror(normal, plane_d, face, sites, jj, EPS_degenerate)) continue;
//
//             if (out) {
//                 s_point sj = sites->p[jj];
//                 double factor = 2 * (plane_d - dot_prod(normal, sj));
//                 // if (fabs(factor / 2.0) < 1e-9) continue;
//                 out->p[sites->N+N_mirror].x = sj.x + factor * normal.x;
//                 out->p[sites->N+N_mirror].y = sj.y + factor * normal.y;
//                 out->p[sites->N+N_mirror].z = sj.z + factor * normal.z;
//             }
//             N_mirror++;
//         }
//     }
//     return N_mirror;
// }





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


