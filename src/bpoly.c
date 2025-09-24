
#include "bpoly.h"
#include "algebra.h"
#include "geometry.h"
#include "float.h"
#include "math.h"
#include <assert.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>


void free_bpoly(s_bound_poly *bpoly)
{
    free_matrix(bpoly->points, bpoly->Np);
    free(bpoly->faces);
    free_matrix(bpoly->fnormals, bpoly->Nf);
    free(bpoly);
}


void add_noise_to_bp(s_bound_poly *bpoly)
{   // ADD SOME NOISE TO AVOID COLINEARITIES... Necessary??
    const double s = 0.01;
    for (int ii=0; ii<bpoly->Np; ii++) {
        for (int jj=0; jj<3; jj++) {
            double aux = 2.0 * rand() / RAND_MAX - 1;  
            bpoly->points[ii][jj] += s * bpoly->dmax * aux;
        }
    }
}


void extract_dmax_bp(s_bound_poly *bpoly)
{
    double dmax = 0; 
    for (int ii=0; ii<bpoly->Np-1; ii++) {
        for (int jj=ii+1; jj<bpoly->Np; jj++) {
            double d = norm_difference_squared(bpoly->points[ii], bpoly->points[jj], 3);
            if (d > dmax) 
                dmax = d;
        }
    }
    bpoly->dmax = sqrt(dmax);
}


void extract_CM_bp(s_bound_poly *bpoly)
{
    find_center_mass(bpoly->points, bpoly->Np, 3, bpoly->CM);
}


void extract_min_max_coord(s_bound_poly *bpoly, double *min, double *max)
{
    min[0] = DBL_MAX;   min[1] = DBL_MAX;   min[2] = DBL_MAX;
    max[0] = -DBL_MAX;  max[1] = -DBL_MAX;  max[2] = -DBL_MAX;
    for (int ii=0; ii<bpoly->Np; ii++) {
        if (bpoly->points[ii][0] < min[0]) min[0] = bpoly->points[ii][0];
        if (bpoly->points[ii][1] < min[1]) min[1] = bpoly->points[ii][1];
        if (bpoly->points[ii][2] < min[2]) min[2] = bpoly->points[ii][2];

        if (bpoly->points[ii][0] > max[0]) max[0] = bpoly->points[ii][0];
        if (bpoly->points[ii][1] > max[1]) max[1] = bpoly->points[ii][1];
        if (bpoly->points[ii][2] > max[2]) max[2] = bpoly->points[ii][2];
    }
}


void extract_convhull_bp(s_bound_poly *bpoly)
{   
    ch_vertex *pch = convert_points_to_chvertex(bpoly->points, bpoly->Np);

    convhull_3d_build(pch, bpoly->Np, &bpoly->faces, &bpoly->Nf);
    // printf("DEBUG: Nf = %d, Np = %d\n", bpoly->Nf, bpoly->Np);
    
    double CM[3];
    find_center_mass(bpoly->points, bpoly->Np, 3, CM);
    double **normals = extract_normals_from_ch(pch, bpoly->faces, bpoly->Nf, CM);
    bpoly->fnormals = normals;
    free(pch);
}


double compute_volume_bpoly(s_bound_poly *bpoly)
{
    double vol = 0;
    for (int ii=0; ii<bpoly->Nf; ii++) {
        int i0 = 0;
        int i1 = 1;
        int i2 = 2;

        double Nx = (bpoly->points[bpoly->faces[ii*3 + 1]][i1] -
                     bpoly->points[bpoly->faces[ii*3 + 0]][i1]) *
                    (bpoly->points[bpoly->faces[ii*3 + 2]][i2] -
                     bpoly->points[bpoly->faces[ii*3 + 0]][i2]) 
                    -
                    (bpoly->points[bpoly->faces[ii*3 + 1]][i2] -
                     bpoly->points[bpoly->faces[ii*3 + 0]][i2]) *
                    (bpoly->points[bpoly->faces[ii*3 + 2]][i1] -
                     bpoly->points[bpoly->faces[ii*3 + 0]][i1]);

        // Nx = bpoly->fnormals[ii][0];
        vol += Nx * (bpoly->points[bpoly->faces[ii*3 + 0]][i0] +
                     bpoly->points[bpoly->faces[ii*3 + 1]][i0] +
                     bpoly->points[bpoly->faces[ii*3 + 2]][i0]);
    }

    return vol/6;
}


s_bound_poly *new_bpoly_from_points(double **points, double Np, int add_noise)
{
    s_bound_poly *bpoly = malloc(sizeof(s_bound_poly));
    bpoly->points = malloc_matrix(Np, 3);
    bpoly->Np = Np;
    copy_matrix(points, bpoly->points, Np, 3);
    
    extract_dmax_bp(bpoly);
    if (add_noise != 0) add_noise_to_bp(bpoly);
    extract_convhull_bp(bpoly);
    extract_CM_bp(bpoly);
    extract_min_max_coord(bpoly, bpoly->min, bpoly->max);

    bpoly->volume = compute_volume_bpoly(bpoly);
    
    return bpoly;
}


s_bound_poly *new_bpoly_from_txt(const char *fname)
{
    FILE *f = fopen(fname, "r");
    if (!f) {
        puts("new_bpoly_from_txt: Could not open file");
        exit(1);
    }
    
    int Np;
    fscanf(f, "%d\n\n", &Np);

    double **points = malloc_matrix(Np, 3);
    for (int ii=0; ii<Np; ii++) {
        fscanf(f, "%lf, %lf, %lf\n", &points[ii][0], &points[ii][1], &points[ii][2]); 
    }

     return new_bpoly_from_points(points, Np, 0);
}


s_bound_poly *new_bpoly_copy(s_bound_poly *in)
{
    s_bound_poly *out = malloc(sizeof(s_bound_poly));
    out->Np = in->Np;

    out->points = malloc_matrix(in->Np, 3);
    copy_matrix(in->points, out->points, in->Np, 3);

    out->Nf = in->Nf;
    out->faces = malloc(sizeof(int) * in->Nf * 3);
    copy_matrix_int(&in->faces, &out->faces, 1, in->Nf * 3);

    out->fnormals = malloc_matrix(in->Nf, 3);
    copy_matrix(in->fnormals, out->fnormals, in->Nf, 3);

    out->dmax = in->dmax;

    out->volume = in->volume;
    
    out->CM[0] = in->CM[0];
    out->CM[1] = in->CM[1];
    out->CM[2] = in->CM[2];

    out->min[0] = in->min[0];
    out->min[1] = in->min[1];
    out->min[2] = in->min[2];

    out->max[0] = in->max[0];
    out->max[1] = in->max[1];
    out->max[2] = in->max[2];

    return out;
}


void extract_vertices_face_bpoly(const s_bound_poly *bpoly, int *face, double **out)
{
    out[0][0] = bpoly->points[face[0]][0];
    out[0][1] = bpoly->points[face[0]][1];
    out[0][2] = bpoly->points[face[0]][2];

    out[1][0] = bpoly->points[face[1]][0];
    out[1][1] = bpoly->points[face[1]][1];
    out[1][2] = bpoly->points[face[1]][2];

    out[2][0] = bpoly->points[face[2]][0];
    out[2][1] = bpoly->points[face[2]][1];
    out[2][2] = bpoly->points[face[2]][2];
}


void scale_bpoly_vertices(double **points, int Np, double s)
{
    double CM[3];
    find_center_mass(points, Np, 3, CM);
    double b[3] = {(1-s)*CM[0],
                   (1-s)*CM[1],
                   (1-s)*CM[2]};
    
    for (int ii=0; ii<Np; ii++) {
        points[ii][0] = s*points[ii][0] + b[0];
        points[ii][1] = s*points[ii][1] + b[1];
        points[ii][2] = s*points[ii][2] + b[2];
    }
}


s_bound_poly *scale_bpoly(s_bound_poly *bp, double factor)
{
    double OLD_VOL = (bp)->volume;
    double **new_p = malloc_matrix((bp)->Np, 3);
    copy_matrix((bp)->points, new_p, (bp)->Np, 3);
    scale_bpoly_vertices(new_p, (bp)->Np, factor);

    int Np = (bp)->Np;
    return(new_bpoly_from_points(new_p, Np, 0));
    printf("DEBUG SCALE_BPOLY: new volume = %f, old volume = %f\n", (bp)->volume, OLD_VOL);
}


s_bound_poly *scale_bpoly_objective_volume(s_bound_poly *bp, double objective_volume)
{
    assert(fabs((bp)->volume) > 1e-6);
    double F = objective_volume / (bp)->volume;
    double s = cbrt(F);

    double **new_p = malloc_matrix((bp)->Np, 3);
    copy_matrix((bp)->points, new_p, (bp)->Np, 3);
    scale_bpoly_vertices(new_p, (bp)->Np, s);

    int Np = (bp)->Np;
    return(new_bpoly_from_points(new_p, Np, 0));
}


int should_mirror(double *n, double *s, double d, double *f1, double *f2, double *f3, 
                  double **all_seeds, int Ns, int seed_id)
{
    double dist = d - dot_3d(n, s);
    // Projection onto the face's plane:
    double p[3] = { s[0] + dist*n[0],
                    s[1] + dist*n[1],
                    s[2] + dist*n[2] };

    // 2) find closest point c on triangle to p
    double c[3];
    closest_point_on_triangle(f1, f2, f3, p, c);
    
    // 3) check nearest‐neighbor at c
    double d_s = norm_difference_squared(c, s, 3);
    for (int j = 0; j < Ns; j++) {
        if (j != seed_id && norm_difference_squared(all_seeds[j], c, 3) + 1e-6 < d_s)
            return 0;   // someone else is nearer
    }

    return 1;
}


int extend_sites_mirroring(s_bound_poly *bp, double ***s, int Ns)
{
    // worst‐case 
    double **out = malloc_matrix(Ns * (1 + bp->Nf), 3);

    for (int ii=0; ii<Ns; ii++) {
        out[ii][0] = (*s)[ii][0];
        out[ii][1] = (*s)[ii][1];
        out[ii][2] = (*s)[ii][2];
    }
    int kk = Ns;

    // 1st order reflections
    for (int ff=0; ff<bp->Nf; ff++) {
        int i0 = bp->faces[ff*3 + 0],
            i1 = bp->faces[ff*3 + 1],
            i2 = bp->faces[ff*3 + 2];
        double *f1 = bp->points[i0],
               *f2 = bp->points[i1],
               *f3 = bp->points[i2];

        double *n = bp->fnormals[ff];
        // assert(fabs(norm_squared(n, 3)-1) < 1e-9);
        double d = dot_3d(n, f1);

        for (int jj=0; jj<Ns; jj++) {
            double *sj = (*s)[jj];
            if (!should_mirror(n, sj, d, f1, f2, f3, *s, Ns, jj))
                continue;
            
            double factor = 2 * (d - dot_3d(n, sj));
            // if (fabs(factor / 2.0) < 1e-9) continue;

            out[kk][0] = sj[0] + factor * n[0];
            out[kk][1] = sj[1] + factor * n[1];
            out[kk][2] = sj[2] + factor * n[2];
            kk++;
        }
    }
    
    // 2nd order reflections:
    // precompute face adjacency: list of (fi,fj, shared_edge_vertices)
    // FaceNeighbor **f_adjacency = compute_face_adjacency(bp);
    // kk = extend_sites_double_reflection(bp, f_adjacency, out, Ns, kk);

    // 3rd order reflections:
    // kk = extend_sites_triple_reflection(bp, out, Ns, kk);

    free_matrix(*s, Ns);
    *s = realloc_matrix(out, Ns * (1 + bp->Nf), kk, 3);
    // printf("DEBUG: Reflected %d sites\n", kk-Ns);
    return kk;
}



// MY IMPLEMENTATION FOR POISSON DISK SAMPLING WITH WEIGHT FUNCTION

void random_point_uniform(double *min, double *max, s_point *out)
{
    double ux = rand() / ((double) RAND_MAX + 1.0);
    double uy = rand() / ((double) RAND_MAX + 1.0);
    double uz = rand() / ((double) RAND_MAX + 1.0);

    out->coords[0] = min[0] + (max[0] - min[0]) * ux;
    out->coords[1] = min[1] + (max[1] - min[1]) * uy;
    out->coords[2] = min[2] + (max[2] - min[2]) * uz;
}


// Generate a random candidate point around a given point p.
// The candidate is generated uniformly in the spherical shell [r, 2r],
// where r = r_of_x(p). We also perturb in all 3 dimensions.
void random_point_around(double *x, double r, double *out)
{
    double radius = r + r * rand()/((double) RAND_MAX + 1);

    // Generate a random direction uniformly over the sphere:
    double aux = 1 - 2.0 * rand()/((double) RAND_MAX + 1);
    if (aux >= 1) aux = 1;
    if (aux <= -1) aux = -1;
    double theta = acos(aux);  // polar angle, 0 <= theta <= pi.
    double phi = 2.0 * M_PI * rand()/((double) RAND_MAX + 1);         // azimuthal, 0 <= phi < 2pi.
    
    out[0] = x[0] + radius * sin(theta) * cos(phi);
    out[1] = x[1] + radius * sin(theta) * sin(phi);
    out[2] = x[2] + radius * cos(theta);
}


int is_valid(s_bound_poly *bpoly, double *q, s_point *samples, int Nsamples, double (*rmax)(double *)) {
    if (!is_inside_convhull(q, bpoly->points, bpoly->faces, bpoly->fnormals, bpoly->Nf)) 
        return 0;

    double rq = rmax(q);
    for (int ii = 0; ii<Nsamples; ii++) {
        double rx = rmax(samples[ii].coords);
        double minDist = fmin(rq, rx);
        if (distance_squared(q, samples[ii].coords) < (minDist * minDist))
            return 0; // candidate too close to an existing sample
    }
    return 1; // valid candidate
}


double **generate_nonuniform_poisson_dist_inside(s_bound_poly *bpoly, double (*rmax)(double *), int *Np_generated)
{
    s_point samples[MAX_TRIAL_POINTS];
    int Nsamples = 0;

    s_point active_list[MAX_TRIAL_POINTS];
    int Nactive = 0;

    s_point x; // RANDOM!!
    random_point_uniform(bpoly->min, bpoly->max, &x);
    while (!is_inside_convhull(x.coords, bpoly->points, bpoly->faces, bpoly->fnormals, bpoly->Nf)) {
        random_point_uniform(bpoly->min, bpoly->max, &x);
    }

    samples[Nsamples++] = x;
    active_list[Nactive++] = x;

    while (Nactive > 0 && Nsamples < MAX_TRIAL_POINTS) {
        int random_id = rand() % Nactive;
        s_point p = active_list[random_id];
        int found = 0;

        double rp = rmax(p.coords);
        for (int ii=0; ii<MAX_TRIAL_TESTS; ii++) {
            s_point q;
            random_point_around(p.coords, rp, q.coords);
            if (is_valid(bpoly, q.coords, samples, Nsamples, rmax)) {
                samples[Nsamples++] = q;
                active_list[Nactive++] = q;
                found = 1;
                break;
            }
        }

        if (Nsamples >= MAX_TRIAL_POINTS-1) {
            fprintf(stderr, "WARNING! Max_trial_points in poisson sampling.\n");
        }
        
        if (found == 0) {
            // Replace activeList[idx] with the last active point and decrease count.
            active_list[random_id] = active_list[Nactive - 1];
            Nactive--;
        }
    }

    while (Nsamples < 4) {
        random_point_uniform(bpoly->min, bpoly->max, &x);
        while (!is_inside_convhull(x.coords, bpoly->points, bpoly->faces, bpoly->fnormals, bpoly->Nf)) {
            random_point_uniform(bpoly->min, bpoly->max, &x);
        }
        samples[Nsamples++] = x;
    }

    double **out_points = malloc_matrix(Nsamples, 3);
    for (int ii=0; ii<Nsamples; ii++) {
        out_points[ii][0] = samples[ii].coords[0];
        out_points[ii][1] = samples[ii].coords[1];
        out_points[ii][2] = samples[ii].coords[2];
    }
    
    // printf("DEBUG POISSON: Nsamples = %d\n", Nsamples);
    // if (Nsamples < 3) {
    //     printf("WARNING! TOO FEW SAMPLES..., N = %d\n", Nsamples);
    //     exit(1);
    // }
    *Np_generated = Nsamples;
    return out_points;
}


void find_closest_point_on_bp(s_bound_poly *bp, double *p, double *OUT)
{
    OUT[0] = DBL_MAX; OUT[1] = DBL_MAX; OUT[2] = DBL_MAX;
    for (int f=0; f<bp->Nf; f++) {
        double tmp[3];

        double *A = bp->points[bp->faces[3*f + 0]];
        double *B = bp->points[bp->faces[3*f + 1]];
        double *C = bp->points[bp->faces[3*f + 2]];
        closest_point_on_triangle(A, B, C, p, tmp);

        if (distance_squared(p, tmp) < distance_squared(p, OUT)) {
            OUT[0] = tmp[0]; OUT[1] = tmp[1]; OUT[2] = tmp[2];
        }
    }
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


void plot_bpoly_with_points(s_bound_poly *bpoly, double **points, int Np, char *f_name, double *ranges, char *color)
{
    FILE *pipe = popen("gnuplot -persistent 2>&1", "w");
    fprintf(pipe, "set terminal pdfcairo enhanced font 'Arial,18' size 4,4 enhanced \n");
    fprintf(pipe, "set pm3d depthorder\n");
    fprintf(pipe, "set pm3d border lc 'black' lw 0.1\n");
    fprintf(pipe, "set xyplane at 0\n");
    fprintf(pipe, "unset border\n");
    fprintf(pipe, "unset xtics\n");
    fprintf(pipe, "unset ytics\n");
    fprintf(pipe, "unset ztics\n");
    if (ranges) {
        fprintf(pipe, "set xrange [%f:%f]\n", ranges[0], ranges[1]);
        fprintf(pipe, "set yrange [%f:%f]\n", ranges[2], ranges[3]);
        fprintf(pipe, "set zrange [%f:%f]\n", ranges[4], ranges[5]);
    }
    fflush(pipe);

    fprintf(pipe, "set view 100, 10, 1.5\n");
    fprintf(pipe, "set output '%s_v0.pdf'\n", f_name);
    fprintf(pipe, "splot ");

    for (int ii=0; ii<bpoly->Nf; ii++) {
            fprintf(pipe, "\"<echo \'");
            fprintf(pipe, "%f %f %f\\n", bpoly->points[bpoly->faces[ii*3]][0],
                                         bpoly->points[bpoly->faces[ii*3]][1], 
                                         bpoly->points[bpoly->faces[ii*3]][2]);
            fprintf(pipe, "%f %f %f\\n", bpoly->points[bpoly->faces[ii*3+1]][0],
                                         bpoly->points[bpoly->faces[ii*3+1]][1], 
                                         bpoly->points[bpoly->faces[ii*3+1]][2]);
            fprintf(pipe, "%f %f %f'\"", bpoly->points[bpoly->faces[ii*3+2]][0],
                                         bpoly->points[bpoly->faces[ii*3+2]][1], 
                                         bpoly->points[bpoly->faces[ii*3+2]][2]);
            if (color) {
                char buff[256];
                snprintf(buff, 256, "w polygons fs transparent solid 0.2 fc '%s' notitle, ", color);
                fprintf(pipe, "%s", buff);
            } else {
                fprintf(pipe, "w polygons fs transparent solid 0.2 fc 'blue' notitle, ");
            }
    }
    fprintf(pipe, "\n");

    fprintf(pipe, "set output '%s_v1.pdf'\n", f_name);
    fprintf(pipe, "set view 100, 90, 1.5\n");
    fprintf(pipe, "replot\n");

    fprintf(pipe, "set output '%s_v2.pdf'\n", f_name);
    fprintf(pipe, "set view 100, 180, 1.5\n");
    fprintf(pipe, "replot\n");

    fprintf(pipe, "set output '%s_v3.pdf'\n", f_name);
    fprintf(pipe, "set view 100, 270, 1.5\n");
    fprintf(pipe, "replot\n");

    // fprintf(pipe, "\"<echo \'");
    // for (int ii=0; ii<bpoly->Np; ii++) {
    //     fprintf(pipe, "%f %f %f\\n", bpoly->points[ii][0], bpoly->points[ii][1], bpoly->points[ii][2]);
    // }
    // fprintf(pipe, "'\" pt 7 lc rgb 'black' notitle, ");

    
    if (points) {
        fprintf(pipe, "\"<echo \'");
        for (int ii=0; ii<Np; ii++) {
            fprintf(pipe, "%f %f %f\\n", points[ii][0], points[ii][1], points[ii][2]);
        }
        fprintf(pipe, "'\" pt 3 lc rgb 'red' notitle, ");
    }

    fprintf(pipe, "\n");
    pclose(pipe);
}


