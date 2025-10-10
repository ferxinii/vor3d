#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "vor3d.h"

#define FILE_BP "bp_points.txt.aux"
#define PLOT_VOLUMES(name) system("./plot_together.plt " name)

double z0, zf, r_mean;

double r_const(double *x)
{   
    (void)x;
    return 1;
}

double r_fun(double *x)
{   
    double K = + 2 * r_mean * 0.1 / (zf - z0);
    return K * (x[2] - (z0 + zf)/2) + r_mean;
}


void check_volume(s_vdiagram *vd)
{
    double sum_vol = 0;
    for (int ii=0; ii<vd->N; ii++) {
        sum_vol += vd->vcells[ii]->volume;
        if (vd->vcells[ii]->volume <= 0) {
            printf("VOL %d : %f\n", ii, vd->vcells[ii]->volume);
        }
    }
    printf("Vol = %f, Diff = %.16f, rel_diff = %.16f\n", vd->bpoly->volume, vd->bpoly->volume - sum_vol, (vd->bpoly->volume - sum_vol) / vd->bpoly->volume);
}


void generate_statistics(s_bpoly *bp, int N_simu, char *FILE_VOLS)
{
    clear_volumes_file(FILE_VOLS);
    for (int ii=0; ii<N_simu; ii++) {
        s_bpoly *bp_tmp = new_bpoly_copy(bp);
        printf("%d, Nf = %d\n", ii, bp_tmp->Nf);
        s_vdiagram *vd = vor3d_from_bp_PDS(&r_fun, bp_tmp, 5);
        append_volumes_to_file(vd, FILE_VOLS, ii);
        check_volume(vd);
        free_vdiagram(vd);
    }
}


int main(void)
{
    srand(time(NULL));
    // system("rm -f plots/*");
    
    int MAX_TRIES = 5;

    puts("\nTETRAHEDON:");
    generate_file_tetrahedron_bp(FILE_BP, 3);
    s_vdiagram *vd_tet = vor3d_from_txt_PDS(&r_const, FILE_BP, MAX_TRIES);
    // s_vdiagram *vd_tet = vor3d_from_txt(test_s, Ns, FILE_BP, MAX_TRIES);
    if (!vd_tet) { puts("Could not construct vd in max_tries."); exit(1); }
    check_volume(vd_tet);
    plot_vdiagram_differentviews(vd_tet, "plots/tet", NULL);
    free_vdiagram(vd_tet);

    puts("\nCUBE:");
    generate_file_cube_bp(FILE_BP, 2);
    s_vdiagram *vd_cube = vor3d_from_txt_PDS(&r_const, FILE_BP, 5);
    if (!vd_cube) { puts("Could not construct vd in max_tries."); exit(1); }
    check_volume(vd_cube);
    plot_vdiagram_differentviews(vd_cube, "plots/cube", NULL);
    free_vdiagram(vd_cube);

    puts("\nSPHERE:");
    generate_file_sphere_bp(FILE_BP, 1.5, 15, 20);
    s_vdiagram *vd_sph = vor3d_from_txt_PDS(&r_const, FILE_BP, 5);
    if (!vd_sph) { puts("Could not construct vd in max_tries."); exit(1); }
    check_volume(vd_sph);
    plot_vdiagram_differentviews(vd_sph, "plots/sph", NULL);
    free_vdiagram(vd_sph);
    
    // Remove tmp file
    char cmd[256];
    snprintf(cmd, 256, "rm -f %s", FILE_BP);
    system(cmd);
    

    // --------------------------------------------------
    // ------------ LUNG VIRTUAL POPULATION -------------
    // --------------------------------------------------
    system("rm -f virtual_population/*png");
    system("rm -f virtual_population/*txt");

    // ------------------ LEFT LUNG ---------------------
    puts("\nLEFT LUNG:");
    s_bpoly *bp_L = new_bpoly_from_txt("lobes/L.txt");
    printf("volume: %f\n", bp_L->volume);
    printf("min: (%f, %f, %f)\n max: (%f, %f, %f)\n", bp_L->min.x, bp_L->min.y, bp_L->min.z, bp_L->max.x, bp_L->max.y, bp_L->max.z);
    plot_bpoly_differentviews(bp_L, "plots/bp_L.png", NULL, "blue");

    
    int Nsimu = 10;
    // ADULT:
    s_bpoly *bp_L_adult = copy_bpoly_scaled(bp_L, 1.06);
    z0 = bp_L_adult->min.z;
    zf = bp_L_adult->max.z;
    r_mean = 1.1;
    generate_statistics(bp_L_adult, Nsimu, "virtual_population/L_adult.txt");
        // PLOT
        // s_vdiagram *vd_L = vor3d_from_txt_PDS(&r_fun, "lobes/L.txt", 5);
        // check_volume(vd_L);
        // plot_vdiagram_differentviews(vd_L, "plots/L", NULL);
        // free_vdiagram(vd_L);

    // 13 Y.O.:
    s_bpoly *bp_L_13 = copy_bpoly_scaled(bp_L, 0.94);
    z0 = bp_L_13->min.z;
    zf = bp_L_13->max.z;
    r_mean = 1;
    generate_statistics(bp_L_13, Nsimu, "virtual_population/L_13yo.txt");

    // 8 Y.O.:
    s_bpoly *bp_L_8 = copy_bpoly_scaled(bp_L, 0.71);
    z0 = bp_L_8->min.z;
    zf = bp_L_8->max.z;
    r_mean = 0.9;
    generate_statistics(bp_L_8, Nsimu, "virtual_population/L_8yo.txt");

    // 3 Y.O.:
    s_bpoly *bp_L_3 = copy_bpoly_scaled(bp_L, 0.52);
    z0 = bp_L_3->min.z;
    zf = bp_L_3->max.z;
    r_mean = 0.7;
    generate_statistics(bp_L_3, Nsimu, "virtual_population/L_3yo.txt");



    // ------------------ LEFT LUNG ---------------------
    puts("\nRIGHT LUNG");
    s_bpoly *bp_R = new_bpoly_from_txt("lobes/R.txt");
    printf("volume: %f\n", bp_R->volume);
    printf("min: (%f, %f, %f)\n max: (%f, %f, %f)\n", bp_R->min.x, bp_R->min.y, bp_R->min.z, bp_R->max.x, bp_R->max.y, bp_R->max.z);
    plot_bpoly_differentviews(bp_R, "plots/bp_R.png", NULL, "orange");

    // ADULT:
    s_bpoly *bp_R_adult = copy_bpoly_scaled(bp_R, 1.06);
    z0 = bp_R_adult->min.z;
    zf = bp_R_adult->max.z;
    r_mean = 1.1;
    generate_statistics(bp_L_adult, Nsimu, "virtual_population/R_adult.txt");
        // PLOT
        // s_vdiagram *vd_R = vor3d_from_txt_PDS(&r_fun, "lobes/R.txt", 5);
        // check_volume(vd_R);
        // plot_vdiagram_differentviews(vd_R, "plots/R", NULL);
        // free_vdiagram(vd_R);

    // 13 Y.O.:
    s_bpoly *bp_R_13 = copy_bpoly_scaled(bp_R, 0.94);
    z0 = bp_R_13->min.z;
    zf = bp_R_13->max.z;
    r_mean = 1;
    generate_statistics(bp_R_13, Nsimu, "virtual_population/R_13yo.txt");

    // 8 Y.O.:
    s_bpoly *bp_R_8 = copy_bpoly_scaled(bp_R, 0.72);
    z0 = bp_R_8->min.z;
    zf = bp_R_8->max.z;
    r_mean = 0.9;
    generate_statistics(bp_R_8, Nsimu, "virtual_population/R_8yo.txt");

    // 3 Y.O.:
    s_bpoly *bp_R_3 = copy_bpoly_scaled(bp_R, 0.52);
    z0 = bp_R_3->min.z;
    zf = bp_R_3->max.z;
    r_mean = 0.7;
    generate_statistics(bp_R_3, Nsimu, "virtual_population/R_3yo.txt");

}
