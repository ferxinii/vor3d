#ifndef BPOLY_C
#define BPOLY_C

#include "algebra.c"
#include "geometry.c"
#include "float.h"
#include <time.h>
// REQUIRES CONVHULL_3D.H
#define TPH_POISSON_IMPLEMENTATION
#include "poisson_disk.h"



typedef struct bound_poly {
    int Np;
    double **points;
    int Nf;
    int *faces;  // Its flat! Nf x 3
    double **fnormals;
    double dmax;  // Max distance between two pairs of points
    double CM[3];
    double min[3];
    double max[3];
} s_bound_poly;


void free_bpoly(s_bound_poly *bpoly)
{
    printf("Nf = %d\n", bpoly->Nf);
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
            double d = norm_difference(bpoly->points[ii], bpoly->points[jj], 3);
            if (d > dmax) 
                dmax = d;
        }
    }
    bpoly->dmax = dmax;
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
    
    double CM[3];
    find_center_mass(bpoly->points, bpoly->Np, 3, CM);
    double **normals = extract_normals_from_ch(pch, bpoly->faces, bpoly->Nf, CM);
    bpoly->fnormals = normals;
    free(pch);
}


s_bound_poly *new_bpoly_from_points(double **points, double Np)
{
    s_bound_poly *bpoly = malloc(sizeof(s_bound_poly));
    bpoly->points = malloc_matrix(Np, 3);
    bpoly->Np = Np;
    copy_matrix(points, bpoly->points, Np, 3);
    
    extract_dmax_bp(bpoly);
    add_noise_to_bp(bpoly);
    extract_convhull_bp(bpoly);
    extract_CM_bp(bpoly);
    extract_min_max_coord(bpoly, bpoly->min, bpoly->max);
    return bpoly;
}


void new_bpoly_from_txt(const char *fname, double ***OUT_points, int *OUT_Np, s_bound_poly **OUT_bpoly)
{
    FILE *f = fopen(fname, "r");
    if (!f) {
        puts("new_bpoly_from_txt: Could not open file");
        exit(1);
    }
    
    fscanf(f, "%d\n\n", OUT_Np);

    *OUT_points = malloc_matrix(*OUT_Np, 3);
    for (int ii=0; ii<*OUT_Np; ii++) {
        fscanf(f, "%lf, %lf, %lf\n", &(*OUT_points)[ii][0], &(*OUT_points)[ii][1], &(*OUT_points)[ii][2]); 
    }
    
    *OUT_bpoly = new_bpoly_from_points(*OUT_points, *OUT_Np);
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


void find_intersection_with_bounding_poly(const s_bound_poly *bpoly, const double *origin, double *dir, double *intersection, double *extension)
{   // in extension, we extend the intersection point in the ray direction a distance of 2*bpoly->dmax
    normalize_inplace(dir, 3);

    static double **face_v = NULL;
    if (!face_v) face_v = malloc_matrix(3, 3);

    double dmin = DBL_MAX;
    int intersected = 0;

    for (int ii=0; ii<bpoly->Nf; ii++) {
        extract_vertices_face_bpoly(bpoly, &bpoly->faces[ii*3], face_v);
        double intersection_aux[3];
        int indicator = ray_triangle_intersection_3d(face_v, origin, dir, intersection_aux);
        if (indicator == 1) {
            double d = norm_difference(intersection_aux, origin, 3);
            if (d < dmin) {
                intersected = 1;
                dmin = d;
                intersection[0] = intersection_aux[0];
                intersection[1] = intersection_aux[1];
                intersection[2] = intersection_aux[2];
            }
        }
    }
    
    assert(bpoly->dmax > 0.1 && "Is bpoly->dmax initialized?");
    double s = 2;
    extension[0] = intersection[0] + s * bpoly->dmax * dir[0];
    extension[1] = intersection[1] + s * bpoly->dmax * dir[1];
    extension[2] = intersection[2] + s * bpoly->dmax * dir[2];

    assert(intersected != 0 && "Did not find any intersection with bounding poly.");
}


double **generate_poisson_dist_inside(s_bound_poly *bpoly, double rmax, int *Np_generated)
{
    double s = 0.1;

    const tph_poisson_real bounds_min[3] = { (tph_poisson_real)(bpoly->min[0] - s * bpoly->dmax),
                                             (tph_poisson_real)(bpoly->min[1] - s * bpoly->dmax), 
                                             (tph_poisson_real)(bpoly->min[2] - s * bpoly->dmax)};
    const tph_poisson_real bounds_max[3] = { (tph_poisson_real)(bpoly->max[0] + s * bpoly->dmax), 
                                             (tph_poisson_real)(bpoly->max[1] + s * bpoly->dmax), 
                                             (tph_poisson_real)(bpoly->max[2] + s * bpoly->dmax)};
    const tph_poisson_args args = { .bounds_min = bounds_min,
                                    .bounds_max = bounds_max,
                                    .radius = (tph_poisson_real)rmax,
                                    .ndims = INT32_C(3),
                                    .max_sample_attempts = UINT32_C(30),
                                    .seed = (uint64_t)time(NULL) };

    tph_poisson_allocator *alloc = NULL;

    tph_poisson_sampling sampling;
    memset(&sampling, 0, sizeof(tph_poisson_sampling));

    puts("DEBUG: NOW CREATING DISTRIBUTION");
    int ret = tph_poisson_create(&args, alloc, &sampling);
    if (ret != TPH_POISSON_SUCCESS) {
      // No need to destroy sampling here!
      printf("Failed creating Poisson sampling! Error code: %d", ret);
      exit(1);
    }

    puts("DEBUG: NOW GETTING SAMPLES");
    const tph_poisson_real *samples = tph_poisson_get_samples(&sampling);

    double **samples_arr = malloc_matrix(sampling.nsamples, 3);
    for (int ii = 0; ii < sampling.nsamples; ii++) {
        samples_arr[ii][0] = samples[ii * sampling.ndims];
        samples_arr[ii][1] = samples[ii * sampling.ndims + 1];
        samples_arr[ii][2] = samples[ii * sampling.ndims + 2];
    }


    // Filter all samples and get only those inside the bpoly
    int *mark_inside = malloc(sizeof(int) * sampling.nsamples);
    mark_inside_convhull(samples_arr, sampling.nsamples, bpoly->points, bpoly->faces, 
                         bpoly->fnormals, bpoly->Nf, mark_inside);
    int N_in = 0;
    for (int ii=0; ii<sampling.nsamples; ii++) {
        if (mark_inside[ii] == 1) N_in++;
    }
    
    double **out = malloc_matrix(N_in, 3);
    int kk = 0;
    for (int ii=0; ii<sampling.nsamples; ii++) {
        if (mark_inside[ii] == 1) {
            copy_matrix(&samples_arr[ii], &out[kk], 1, 3);
            kk++;
        }
    }

    free(mark_inside);
    free_matrix(samples_arr, sampling.nsamples);
    tph_poisson_destroy(&sampling);
    *Np_generated = N_in;
    return out;
}


void plot_bpoly_with_points(s_bound_poly *bpoly, double **points, int Np, char *f_name, double *ranges)
{
    FILE *pipe = popen("gnuplot -persistent 2>&1", "w");
    fprintf(pipe, "set terminal pngcairo enhanced font 'Arial,18' size 1080,1080 enhanced \n");
    fprintf(pipe, "set output '%s.png'\n", f_name);
    fprintf(pipe, "set pm3d depthorder\n");
    fprintf(pipe, "set view 100, 60, \n");
    fprintf(pipe, "set xyplane at 0\n");
    fprintf(pipe, "set xrange [%f:%f]\n", ranges[0], ranges[1]);
    fprintf(pipe, "set yrange [%f:%f]\n", ranges[2], ranges[3]);
    fprintf(pipe, "set zrange [%f:%f]\n", ranges[4], ranges[5]);
    fprintf(pipe, "set xlabel 'x'\n");
    fprintf(pipe, "set ylabel 'y'\n");
    fprintf(pipe, "set zlabel 'z'\n");
    fflush(pipe);
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
            fprintf(pipe, "w polygons fs transparent solid 0.05 fc 'black' notitle, ");
    }

    fprintf(pipe, "\"<echo \'");
    for (int ii=0; ii<Np; ii++) {
        fprintf(pipe, "%f %f %f\\n", points[ii][0], points[ii][1], points[ii][2]);
    }
    fprintf(pipe, "'\" pt 7 lc rgb 'black' notitle, ");

    fprintf(pipe, "\n");

}

#endif
