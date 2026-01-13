#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <omp.h>
#include <math.h>
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

void clear_volumes_file(char *fname)
{
    FILE *f = fopen(fname, "w");
    fclose(f);
}

void append_volumes_to_file(s_vdiagram *vdiagram, char *fname, int id)
{
    FILE *f = fopen(fname, "a");
    assert(f && "Could not open file to write volumes to");

    for (int ii=0; ii<vdiagram->seeds.N; ii++) {
        fprintf(f, "%d, %f, %f, %f, %f\n", id, vdiagram->seeds.p[vdiagram->vcells[ii].seed_id].x,
                                               vdiagram->seeds.p[vdiagram->vcells[ii].seed_id].y, 
                                               vdiagram->seeds.p[vdiagram->vcells[ii].seed_id].z, 
                                               vdiagram->vcells[ii].volume);   
    }
    fclose(f);
}

static void generate_statistics(const s_bpoly *bp, int N_simu, char *FILE_VOLS)
{
    clear_volumes_file(FILE_VOLS);
    for (int ii=0; ii<N_simu; ii++) {
        s_vdiagram vd = vor3d_from_bp_PDS(&r_fun, &r_fun_params, bp, vol_max_rel_diff, 5, EPS_degenerate, TOL);
        append_volumes_to_file(&vd, FILE_VOLS, ii);
        check_volume(&vd);
        free_vdiagram(&vd);
    }
}

int main(void)
{
    srand(time(NULL));
    // system("rm -f plots/*");
    
    int Nsimu = 10;  /* Number of sample lungs for each configuration */

    // --------------------------------------------------
    // ------------ LUNG VIRTUAL POPULATION -------------
    // --------------------------------------------------
    system("rm -f virtual_population/*png");
    system("rm -f virtual_population/*txt");

    puts("\nLEFT LUNG:");
    s_bpoly bp_L = bpoly_from_csv("lobes/L.txt", EPS_degenerate, TOL);
    printf("volume: %f\n", bp_L.volume);
    printf("min: (%f, %f, %f)\n max: (%f, %f, %f)\n", bp_L.min.x, bp_L.min.y, bp_L.min.z, bp_L.max.x, bp_L.max.y, bp_L.max.z);
    // plot_bpoly_differentviews(&bp_L, "plots/bp_L.png", NULL, "blue");
    s_point CM_L = point_average(&bp_L.convh.points);

    puts("\nRIGHT LUNG");
    s_bpoly bp_R = bpoly_from_csv("lobes/R.txt", EPS_degenerate, TOL);
    printf("volume: %f\n", bp_R.volume);
    printf("min: (%f, %f, %f)\n max: (%f, %f, %f)\n", bp_R.min.x, bp_R.min.y, bp_R.min.z, bp_R.max.x, bp_R.max.y, bp_R.max.z);
    // plot_bpoly_differentviews(&bp_R, "plots/bp_R.png", NULL, "orange");
    s_point CM_R = point_average(&bp_R.convh.points);

    double TOT_VOL = bp_L.volume + bp_R.volume;


    // ------------------ LEFT LUNG ---------------------
    puts("\n\n Now generating stats for right lung.\n");

    // ADULT:
    double s = cbrt(5000.0 / TOT_VOL);
    s_bpoly bp_L_adult = copy_bpoly_scaled(&bp_L, s, CM_L, EPS_degenerate, TOL);
    r_fun_params.z0 = bp_L_adult.min.z;
    r_fun_params.zf = bp_L_adult.max.z;
    r_fun_params.r_mean = 1.1;
    printf("ADULT:    s: %g,    Vol %g,    params: (%g, %g, %g)\n", s, bp_L_adult.volume, r_fun_params.z0, r_fun_params.zf, r_fun_params.r_mean);
    generate_statistics(&bp_L_adult, Nsimu, "virtual_population/L_adult.txt");
        // PLOT
        // s_vdiagram *vd_L = vor3d_from_txt_PDS(&r_fun, "lobes/L.txt", 5);
        // check_volume(vd_L);
        // plot_vdiagram_differentviews(vd_L, "plots/L", NULL);
        // free_vdiagram(vd_L);

    // 13 Y.O.:
    s = cbrt(3000.0 / TOT_VOL);
    s_bpoly bp_L_13 = copy_bpoly_scaled(&bp_L, s, CM_L, EPS_degenerate, TOL);
    r_fun_params.z0 = bp_L_13.min.z;
    r_fun_params.zf = bp_L_13.max.z;
    r_fun_params.r_mean = 1;
    printf("13YO:    s: %g,    Vol %g,    params: (%g, %g, %g)\n", s, bp_L_13.volume, r_fun_params.z0, r_fun_params.zf, r_fun_params.r_mean);
    generate_statistics(&bp_L_13, Nsimu, "virtual_population/L_13yo.txt");

    // 8 Y.O.:
    s = cbrt(1500.0 / TOT_VOL);
    s_bpoly bp_L_8 = copy_bpoly_scaled(&bp_L, s, CM_L, EPS_degenerate, TOL);
    r_fun_params.z0 = bp_L_8.min.z;
    r_fun_params.zf = bp_L_8.max.z;
    r_fun_params.r_mean = 0.9;
    printf("8YO:    s: %g,    Vol %g,    params: (%g, %g, %g)\n", s, bp_L_8.volume, r_fun_params.z0, r_fun_params.zf, r_fun_params.r_mean);
    generate_statistics(&bp_L_8, Nsimu, "virtual_population/L_8yo.txt");

    // 3 Y.O.:
    s = cbrt(600.0 / TOT_VOL);
    s_bpoly bp_L_3 = copy_bpoly_scaled(&bp_L, s, CM_L, EPS_degenerate, TOL);
    r_fun_params.z0 = bp_L_3.min.z;
    r_fun_params.zf = bp_L_3.max.z;
    r_fun_params.r_mean = 0.7;
    printf("3YO:    s: %g,    Vol %g,    params: (%g, %g, %g)\n", s, bp_L_3.volume, r_fun_params.z0, r_fun_params.zf, r_fun_params.r_mean);
    generate_statistics(&bp_L_3, Nsimu, "virtual_population/L_3yo.txt");



    // ------------------ RIGHT LUNG ---------------------
    puts("\n\n Now generating stats for right lung.\n");

    // ADULT:
    s = cbrt(5000.0 / TOT_VOL);
    s_bpoly bp_R_adult = copy_bpoly_scaled(&bp_R, 1.06, CM_R, EPS_degenerate, TOL);
    r_fun_params.z0 = bp_R_adult.min.z;
    r_fun_params.zf = bp_R_adult.max.z;
    r_fun_params.r_mean = 1.1;
    printf("ADULT:    s: %g,    Vol %g,    params: (%g, %g, %g)\n", s, bp_R_adult.volume, r_fun_params.z0, r_fun_params.zf, r_fun_params.r_mean);
    generate_statistics(&bp_R_adult, Nsimu, "virtual_population/R_adult.txt");
        // PLOT
        // s_vdiagram *vd_R = vor3d_from_txt_PDS(&r_fun, "lobes/R.txt", 5);
        // check_volume(vd_R);
        // plot_vdiagram_differentviews(vd_R, "plots/R", NULL);
        // free_vdiagram(vd_R);

    // 13 Y.O.:
    s = cbrt(3000.0 / TOT_VOL);
    s_bpoly bp_R_13 = copy_bpoly_scaled(&bp_R, s, CM_R, EPS_degenerate, TOL);
    r_fun_params.z0 = bp_R_13.min.z;
    r_fun_params.zf = bp_R_13.max.z;
    r_fun_params.r_mean = 1;
    printf("13YO:    s: %g,    Vol %g,    params: (%g, %g, %g)\n", s, bp_R_13.volume, r_fun_params.z0, r_fun_params.zf, r_fun_params.r_mean);
    generate_statistics(&bp_R_13, Nsimu, "virtual_population/R_13yo.txt");

    // 8 Y.O.:
    s = cbrt(1500.0 / TOT_VOL);
    s_bpoly bp_R_8 = copy_bpoly_scaled(&bp_R, s, CM_R, EPS_degenerate, TOL);
    r_fun_params.z0 = bp_R_8.min.z;
    r_fun_params.zf = bp_R_8.max.z;
    r_fun_params.r_mean = 0.9;
    printf("8YO:    s: %g,    Vol %g,    params: (%g, %g, %g)\n", s, bp_R_8.volume, r_fun_params.z0, r_fun_params.zf, r_fun_params.r_mean);
    generate_statistics(&bp_R_8, Nsimu, "virtual_population/R_8yo.txt");

    // 3 Y.O.:
    s = cbrt(600.0 / TOT_VOL);
    s_bpoly bp_R_3 = copy_bpoly_scaled(&bp_R, s, CM_R, EPS_degenerate, TOL);
    r_fun_params.z0 = bp_R_3.min.z;
    r_fun_params.zf = bp_R_3.max.z;
    r_fun_params.r_mean = 0.7;
    printf("3YO:    s: %g,    Vol %g,    params: (%g, %g, %g)\n", s, bp_R_3.volume, r_fun_params.z0, r_fun_params.zf, r_fun_params.r_mean);
    generate_statistics(&bp_R_3, Nsimu, "virtual_population/R_3yo.txt");
}

