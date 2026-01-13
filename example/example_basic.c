#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <omp.h>
#include "vor3d.h"

#define FILE_BP "bp_points.txt.aux"
#define PLOT_VOLUMES(name) system("./plot_together.plt " name)

double EPS_degenerate = 1e-12, TOL = 1e-12;
double vol_max_rel_diff = 1e-3;

double r_const(double *x, void *params)
{   
    (void)params;
    (void)x;
    // return NAN;
    return 1;
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

static void iterate(int iterations, int MAX_TRIES, bool print_progress)
{
    int fail = 0;
    #pragma omp parallel for
    for (int ii=0; ii<iterations; ii++) {
        if (print_progress && ii%100 == 0) printf("%d / %d\n", ii, iterations);
        s_vdiagram vd = vor3d_from_txt_PDS(&r_const, NULL, FILE_BP, vol_max_rel_diff, MAX_TRIES, EPS_degenerate, TOL);
        if (vd.seeds.N == 0) {fail++; continue;}
        // check_volume(&vd);
        free_vdiagram(&vd);
    }

    printf("FAILED %d / %d\n", fail, iterations);
}

int main(void)
{
    srand(time(NULL));
    // system("rm -f plots/*");
    
    int MAX_TRIES = 1;
    bool PLOT = false;
    (void)PLOT, (void)MAX_TRIES;

    // generate_file_tetrahedron_bp(FILE_BP, 3);
    // // generate_file_sphere_bp(FILE_BP, 1.5, 15, 20);
    // s_points ptest = read_points_from_csv("../seeds.csv");
    // s_vdiagram vd_test = vor3d_from_txt(&ptest, FILE_BP, vol_max_rel_diff, MAX_TRIES, 1e-14, TOL);
    // check_volume(&vd_test);
    // // plot_vdiagram_differentviews(&vd_tet, "plots/tet", NULL);
    // exit(1);

    puts("\nTETRAHEDON:");
    generate_file_tetrahedron_bp(FILE_BP, 3);
    s_vdiagram vd_tet = vor3d_from_txt_PDS(&r_const, NULL, FILE_BP, vol_max_rel_diff, MAX_TRIES, EPS_degenerate, TOL);
    check_volume(&vd_tet);
    iterate(10000, MAX_TRIES, false);
    if (PLOT) plot_vdiagram_differentviews(&vd_tet, "plots/tet", NULL);
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
    s_vdiagram vd_cube = vor3d_from_txt_PDS(&r_const, NULL, FILE_BP, vol_max_rel_diff, 5, EPS_degenerate, TOL);
    if (vd_cube.seeds.N == 0) { puts("Could not construct vd in max_tries."); exit(1); }
    check_volume(&vd_cube);
    if (PLOT) plot_vdiagram_differentviews(&vd_cube, "plots/cube", NULL);
    free_vdiagram(&vd_cube);
    iterate(10000, MAX_TRIES, false);
    // exit(1);
    

    puts("\nSPHERE:");
    generate_file_sphere_bp(FILE_BP, 1.5, 15, 20);
    s_vdiagram vd_sph = vor3d_from_txt_PDS(&r_const, NULL, FILE_BP, vol_max_rel_diff, MAX_TRIES, EPS_degenerate, TOL);
    write_convhull_to_m(&vd_sph.bpoly.convh, "test_bp.m");
    if (vd_sph.seeds.N == 0) { puts("Could not construct vd in max_tries."); exit(1); }
    check_volume(&vd_sph);
    if (PLOT) plot_vdiagram_differentviews(&vd_sph, "plots/sph", NULL);
    free_vdiagram(&vd_sph);
    iterate(100, MAX_TRIES, true);
}

