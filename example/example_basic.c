#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <omp.h>
#include <math.h>
#include "vor3d.h"
#include "random.h"
#include "dynarray.h"

#define FILE_BP "bp_points.txt.aux"
#define PLOT_VOLUMES(name) system("./plot_together.plt " name)

double EPS_degenerate = 1e-12, TOL = 1e-12;
double vol_max_rel_diff = 1e-3;
s_random_context *rctx_omp;

double r_const(double *x, void *params)
{   
    (void)params;
    (void)x;
    // return NAN;
    return 1;
}

int randint(void *ctx, int N)
{
    s_random_context *C = ctx;
    return random_uniform_range_u64(C, N);
}

double randd01(void *ctx)
{
    s_random_context *C = ctx;
    return random_uniform_double(C);
}

typedef struct {
    double z0, zf, r_mean;
} s_params;

s_params r_fun_params;

double r_fun(double *x, void *params)
{   
    s_params *p = params;
    double K = + 2 * p->r_mean * 0.1 / (p->zf - p->z0);
    return K * (x[2] - (p->z0 + p->zf)/2) + p->r_mean;
}


static void check_volume(const s_vdiagram *vd)
{
    double sum_vol = 0;
    for (int ii=0; ii<vd->seeds.N; ii++) {
        sum_vol += vd->vcells[ii].volume;
        if (vd->vcells[ii].volume <= 0) {
            printf("VOL %d : %f\n", ii, vd->vcells[ii].volume);
        }
    }
    printf("Vol = %f, Diff = %.16f, rel_diff = %.16f\n", vd->bpoly.volume, vd->bpoly.volume - sum_vol, (vd->bpoly.volume - sum_vol) / vd->bpoly.volume);
}

static void iterate(int iterations, bool print_progress, const s_bpoly *bp,
                    s_dynarray *buff_points_omp)
{
    int fail = 0;
    #pragma omp parallel for
    for (int ii=0; ii<iterations; ii++) {
        // printf("%d\n", ii);
        if (print_progress && ii%100 == 0) printf("%d / %d\n", ii, iterations);
        s_vdiagram vd = vor3d_inside_bp_PDS(&r_const, NULL, bp, vol_max_rel_diff,
                                          EPS_degenerate, TOL,
                                          randint, randd01, &rctx_omp[omp_get_thread_num()],
                                          &buff_points_omp[omp_get_thread_num()]);
        if (vd.seeds.N == 0) {fail++; continue;}
        // check_volume(&vd);
        free_vdiagram(&vd);
    }

    printf("FAILED %d / %d\n", fail, iterations);
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


int main(void)
{
    rctx_omp = malloc(sizeof(s_random_context) * omp_get_max_threads());
    random_initialize_threads(time(NULL), 1, omp_get_max_threads(), rctx_omp);

    s_dynarray *buff_points_omp = malloc(sizeof(s_dynarray) * omp_get_max_threads());
    for (int i=0; i<omp_get_max_threads(); i++) {
        buff_points_omp[i] = dynarray_initialize(sizeof(s_point), 0);
    }
    
    bool PLOT = false;
    (void)PLOT;

    puts("\nTETRAHEDON:");
    generate_file_tetrahedron_bp(FILE_BP, 3);
    s_bpoly bp_tet = bpoly_from_csv(FILE_BP, EPS_degenerate);
    s_vdiagram vd_tet = vor3d_inside_bp_PDS(&r_const, NULL, &bp_tet, vol_max_rel_diff,
                                          EPS_degenerate, TOL, randint, randd01, &rctx_omp[0],
                                          buff_points_omp);
    check_volume(&vd_tet);
    free_vdiagram(&vd_tet);
    iterate(10000, false, &bp_tet, buff_points_omp);
    if (PLOT) plot_vdiagram_differentviews(&vd_tet, "plots/tet", NULL);
    free_bpoly(&bp_tet);
    // exit(1);

    // /* Testing serialization */
    // size_t size;
    // uint8_t *serialized;
    // serialize_vdiagram(&vd_tet, NULL, &size, &serialized);
    // write_serialized_vdiagram("tet_serialized.bin", serialized, size);
    //
    // size_t new_size;
    // uint8_t *new_serialized;
    // read_serialized_vdiagram("tet_serialized.bin", &new_serialized, &new_size);
    // s_vdiagram new_vd_tet;
    // deserialize_vdiagram(new_serialized, &new_vd_tet, NULL);
    // plot_vdiagram_differentviews(&vd_tet, "plots/new_tet", NULL);
    // printf("%d\n%d\n", (int)size, (int)new_size);
    // free_vdiagram(&vd_tet);

    puts("\nCUBE:");
    generate_file_cube_bp(FILE_BP, 2);
    s_bpoly bp_cube = bpoly_from_csv(FILE_BP, EPS_degenerate);
    s_vdiagram vd_cube = vor3d_inside_bp_PDS(&r_const, NULL, &bp_cube, vol_max_rel_diff,
                                           EPS_degenerate, TOL, randint, randd01, &rctx_omp[0],
                                           buff_points_omp);
    if (!vdiagram_is_valid(&vd_cube)) { puts("Could not construct vd."); exit(1); }
    check_volume(&vd_cube);
    if (PLOT) plot_vdiagram_differentviews(&vd_cube, "plots/cube", NULL);
    free_vdiagram(&vd_cube);
    iterate(10000, false, &bp_cube, buff_points_omp);
    free_bpoly(&bp_cube);
    // exit(1);


    puts("\nSPHERE:");
    generate_file_sphere_bp(FILE_BP, 1.5, 15, 20);
    s_bpoly bp_sph = bpoly_from_csv(FILE_BP, EPS_degenerate);
    s_vdiagram vd_sph = vor3d_inside_bp_PDS(&r_const, NULL, &bp_sph, vol_max_rel_diff,
                                          EPS_degenerate, TOL, randint, randd01, &rctx_omp[0],
                                          buff_points_omp);
    write_convhull_to_m(&vd_sph.bpoly.convh, "test_bp.m");
    if (!vdiagram_is_valid(&vd_sph)) { puts("Could not construct vd."); exit(1); }
    check_volume(&vd_sph);
    if (PLOT) plot_vdiagram_differentviews(&vd_sph, "plots/sph", NULL);
    free_vdiagram(&vd_sph);
    iterate(100, true, &bp_sph, buff_points_omp);
    free_bpoly(&bp_sph);


    free(rctx_omp);
}

