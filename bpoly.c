#ifndef BPOLY_C
#define BPOLY_C

#include "algebra.c"
#include "geometry.c"
#include "float.h"
#include <time.h>
// REQUIRES CONVHULL_3D.H
#define TPH_POISSON_IMPLEMENTATION
#include "poisson_disk.h"


#define MAX_TRIAL_POINTS 10000
#define MAX_TRIAL_TESTS 100
typedef struct point {
    double coords[3];
} s_point;


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
    double volume;
} s_bound_poly;


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
    printf("DEBUG: Nf = %d, Np = %d\n", bpoly->Nf, bpoly->Np);
    
    double CM[3];
    find_center_mass(bpoly->points, bpoly->Np, 3, CM);
    double **normals = extract_normals_from_ch(pch, bpoly->faces, bpoly->Nf, CM);
    bpoly->fnormals = normals;
    free(pch);
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

    // COMPUTE VOLUME
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

        vol += Nx * (bpoly->points[bpoly->faces[ii*3 + 0]][i0] +
                     bpoly->points[bpoly->faces[ii*3 + 1]][i0] +
                     bpoly->points[bpoly->faces[ii*3 + 2]][i0]);
    }
    bpoly->volume = vol / 6;

    return bpoly;
}


void new_bpoly_from_txt(const char *fname, double ***OUT_points, int *OUT_Np, s_bound_poly **OUT_bpoly, int add_noise)
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


    *OUT_bpoly = new_bpoly_from_points(*OUT_points, *OUT_Np, add_noise);
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


double **generate_uniform_poisson_dist_inside(s_bound_poly *bpoly, double rmax, int *Np_generated)
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

    // puts("DEBUG: NOW CREATING DISTRIBUTION");
    int ret = tph_poisson_create(&args, alloc, &sampling);
    if (ret != TPH_POISSON_SUCCESS) {
      // No need to destroy sampling here!
      printf("Failed creating Poisson sampling! Error code: %d", ret);
      exit(1);
    }

    // puts("DEBUG: NOW GETTING SAMPLES");
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


int is_plane_saved(double **planes, int N, double *n, double d)
{
    double EPS = 1e-9;
    for (int ii=0; ii<N; ii++) {
        if (norm_difference(planes[ii], n, 3) < EPS && 
            fabs(planes[ii][3] - d) < EPS)
            return 1;
    }
    return 0;
}


void closest_point_on_triangle(const double *A, const double *B, const double *C, const double *p,
                               double *c_out) 
{
    double AB[3] = {B[0]-A[0], B[1]-A[1], B[2]-A[2]},
           AC[3] = {C[0]-A[0], C[1]-A[1], C[2]-A[2]},
           AP[3] = {p [0]-A[0], p [1]-A[1], p [2]-A[2]};

    double d1 = dot_3d(AB, AP),
           d2 = dot_3d(AC, AP);
    if (d1 <= 0 && d2 <= 0) {          // Vertex region A
        memcpy(c_out, A, 3*sizeof(double));
        return;
    }

    // Check vertex B region
    double BP[3] = { p[0]-B[0], p[1]-B[1], p[2]-B[2] };
    double d3 = dot_3d(AB, BP),
           d4 = dot_3d(AC, BP);
    if (d3 >= 0 && d4 <= d3) {         // Vertex region B
        memcpy(c_out, B, 3*sizeof(double));
        return;
    }

    // Check edge AB region
    double vc = d1*d4 - d3*d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) {
        double v = d1 / (d1 - d3);
        c_out[0] = A[0] + v * AB[0];
        c_out[1] = A[1] + v * AB[1];
        c_out[2] = A[2] + v * AB[2];
        return;
    }

    // Check vertex C region
    double CP[3] = { p[0]-C[0], p[1]-C[1], p[2]-C[2] };
    double d5 = dot_3d(AB, CP),
           d6 = dot_3d(AC, CP);
    if (d6 >= 0 && d5 <= d6) {         // Vertex region C
        memcpy(c_out, C, 3*sizeof(double));
        return;
    }

    // Check edge AC region
    double vb = d5*d2 - d1*d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {
        double w = d2 / (d2 - d6);
        c_out[0] = A[0] + w * AC[0];
        c_out[1] = A[1] + w * AC[1];
        c_out[2] = A[2] + w * AC[2];
        return;
    }

    // Check edge BC region
    double va = d3*d6 - d5*d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        // projection onto BC
        double BC[3] = {C[0]-B[0], C[1]-B[1], C[2]-B[2]};
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        c_out[0] = B[0] + w * BC[0];
        c_out[1] = B[1] + w * BC[1];
        c_out[2] = B[2] + w * BC[2];
        return;
    }

    // Inside face region
    // Barycentric coordinates (u,v,w)
    double denom = va + vb + vc;
    assert(fabs(denom) > 1e-9);
    double v = vb / denom;
    double w = vc / denom;
    c_out[0] = A[0] + AB[0]*v + AC[0]*w;
    c_out[1] = A[1] + AB[1]*v + AC[1]*w;
    c_out[2] = A[2] + AB[2]*v + AC[2]*w;
}



int should_mirror(double *n, double *s, double d, double *f1, double *f2, double *f3, 
                  double **all_seeds, int Ns, int seed_id)
{
    double dist = dot_3d(n, s) - d;
    // Projection onto the face's plane:
    double p[3] = { s[0] - dist*n[0],
                    s[1] - dist*n[1],
                    s[2] - dist*n[2] };

    // 2) find closest point c on triangle to p
    double c[3];
    closest_point_on_triangle(f1, f2, f3, p, c);
    
    // 3) check nearest‐neighbor at c
    double d_s = norm_difference(c, s, 3);
    for (int j = 0; j < Ns; j++) {
        if (j == seed_id) continue;
        if (norm_difference(all_seeds[j], c, 3) < d_s)
            return 0;   // someone else is nearer
    }
    return 1;
}


int extend_sites_mirroring_NEW(s_bound_poly *bp, double ***s, int Ns)
{
    // worst‐case 
    double **out = malloc_matrix(Ns * (1 + bp->Nf), 3);

    for (int ii=0; ii<Ns; ii++) {
        out[ii][0] = (*s)[ii][0];
        out[ii][1] = (*s)[ii][1];
        out[ii][2] = (*s)[ii][2];
    }
    int kk = Ns;

    for (int ff=0; ff<bp->Nf; ff++) {
        int   i0 = bp->faces[ff*3 + 0],
              i1 = bp->faces[ff*3 + 1],
              i2 = bp->faces[ff*3 + 2];
        double *f1 = bp->points[i0],
               *f2 = bp->points[i1],
               *f3 = bp->points[i2];

        double *n = bp->fnormals[ff];
        double d = dot_3d(n, f1);

        for (int jj=0; jj<Ns; jj++) {
            double *sj = (*s)[jj];
            if (!should_mirror(n, sj, d, f1, f2, f3, *s, Ns, jj))
                continue;

            double dist   = dot_3d(n, sj) - d;
            double factor = -2.0 * dist;
            out[kk][0] = sj[0] + factor * n[0];
            out[kk][1] = sj[1] + factor * n[1];
            out[kk][2] = sj[2] + factor * n[2];
            kk++;
        }
    }

    free_matrix(*s, Ns);
    *s = realloc_matrix(out, Ns * (1 + bp->Nf), kk, 3);
    return kk;
}


int extend_sites_mirroring(s_bound_poly *bp, double ***s, int Ns)
{
    double **out = malloc_matrix(Ns*(1+bp->Nf), 3);
    for (int ii=0; ii<Ns; ii++) {
        out[ii][0] = (*s)[ii][0];
        out[ii][1] = (*s)[ii][1];
        out[ii][2] = (*s)[ii][2];
    }

    double **stored_planes = malloc_matrix(bp->Nf, 4);
    int kk = Ns, Np_unique = 0;
    for (int ii=0; ii<bp->Nf; ii++) {
        double ni[3] = {bp->fnormals[ii][0],
                        bp->fnormals[ii][1],
                        bp->fnormals[ii][2]};
        double di = dot_3d(bp->fnormals[ii], 
                           bp->points[bp->faces[ii*3]]);

        assert(fabs(norm_squared(ni, 3) - 1) < 1e-9);
        

        if (!is_plane_saved(stored_planes, Np_unique, ni, di)) {
            stored_planes[Np_unique][0] = ni[0];
            stored_planes[Np_unique][1] = ni[1];
            stored_planes[Np_unique][2] = ni[2];
            stored_planes[Np_unique][3] = di;
            Np_unique++;
            for (int jj=0; jj<Ns; jj++) {
                double factor = 2 * (di - dot_3d(ni, (*s)[jj]));
                assert(fabs(factor) > 1e-7); 
                out[kk][0] = (*s)[jj][0] +  factor * ni[0];
                out[kk][1] = (*s)[jj][1] +  factor * ni[1];
                out[kk][2] = (*s)[jj][2] +  factor * ni[2];
                kk++;
            }
        }
    }

    free_matrix(*s, Ns);
    free_matrix(stored_planes, bp->Nf);
    *s = realloc_matrix(out, Ns*(1+bp->Nf), kk, 3);
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
    puts("DEBUG POISSON: Finding initial point...");
    random_point_uniform(bpoly->min, bpoly->max, &x);
    while (!is_inside_convhull(x.coords, bpoly->points, bpoly->faces, bpoly->fnormals, bpoly->Nf)) {
        printf("DEBUG POISSON: (%f, %f, %f) : %d\n", x.coords[0], x.coords[1], x.coords[2], 
                is_inside_convhull(x.coords, bpoly->points, bpoly->faces, bpoly->fnormals, bpoly->Nf));
        random_point_uniform(bpoly->min, bpoly->max, &x);
    }
    puts("DEBUG POISSON: Found!");

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

    if (Nsamples < 2) {
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
    
    printf("DEBUG POISSON: Nsamples = %d\n", Nsamples);
    if (Nsamples < 4) {
        printf("WARNING! TOO FEW SAMPLES..., N = %d\n", Nsamples);
        exit(1);
    }
    *Np_generated = Nsamples;
    return out_points;
}







void plot_bpoly_with_points(s_bound_poly *bpoly, double **points, int Np, char *f_name, double *ranges)
{
    FILE *pipe = popen("gnuplot -persistent 2>&1", "w");
    fprintf(pipe, "set terminal pngcairo enhanced font 'Arial,18' size 1080,1080 enhanced \n");
    fprintf(pipe, "set output '%s.png'\n", f_name);
    fprintf(pipe, "set pm3d depthorder\n");
    fprintf(pipe, "set pm3d border lc 'black' lw 0.5\n");
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
    for (int ii=0; ii<bpoly->Np; ii++) {
        fprintf(pipe, "%f %f %f\\n", bpoly->points[ii][0], bpoly->points[ii][1], bpoly->points[ii][2]);
    }
    fprintf(pipe, "'\" pt 7 lc rgb 'black' notitle, ");


    fprintf(pipe, "\"<echo \'");
    for (int ii=0; ii<Np; ii++) {
        fprintf(pipe, "%f %f %f\\n", points[ii][0], points[ii][1], points[ii][2]);
    }
    fprintf(pipe, "'\" pt 3 lc rgb 'red' notitle, ");

    fprintf(pipe, "\n");
    pclose(pipe);
}

#endif
