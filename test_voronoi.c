#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "voronoi.h"


#define FILE_BP "bp_points.txt"
#define PLOT_VOLUMES(name) system("./plot_volumes.plt " name)

double z0, zf, r_mean;

double r_fun(double *x)
{   
    double K = + 2 * r_mean * 0.1 / (zf - z0);
    // printf("%f\n", K * (x[2] - (z0 + zf)/2) + r_mean);
    return K * (x[2] - (z0 + zf)/2) + r_mean;
}


void check_volume(s_vdiagram *vd)
{
    double sum_vol = 0;
    for (int ii=0; ii<vd->N_vcells; ii++) {
        sum_vol += vd->vcells[ii]->volume;
        if (vd->vcells[ii]->volume <= 0) {
            printf("VOL %d : %f\n", ii, vd->vcells[ii]->volume);
        }
    }
    printf("V lung = %f, Diff = %.16f, rel_diff = %.16f\n", vd->bpoly->volume, vd->bpoly->volume - sum_vol, (vd->bpoly->volume - sum_vol) / vd->bpoly->volume);
}


void generate_statistics(s_bound_poly *bp, int N_simu, char *FILE_VOLS)
{
    clear_volumes_file(FILE_VOLS);
    for (int ii=0; ii<N_simu; ii++) {
        s_bound_poly *bp_tmp = new_bpoly_copy(bp);
        printf("%d, Nf = %d\n", ii, bp_tmp->Nf);
        s_vdiagram *vd = construct_vd_from_bp(&r_fun, bp_tmp, 5);
        append_volumes_to_file(vd, FILE_VOLS, ii);
        check_volume(vd);
        free_vdiagram(vd);
    }
}


int main(void)
{
    srand(time(NULL));
    system("rm -f plot_vd/*");
    

    puts("\nTETRAHEDON:");
    // generate_file_tetrahedron_bp(FILE_BP, 3);
    // s_vdiagram *vd_tet = construct_vd_from_txt(&r_fun, FILE_BP, 5);
    // if (!vd_tet) { puts("Could not construct vd in max_tries."); exit(1); }
    // check_volume(vd_tet);
    // plot_vdiagram_auto(vd_tet, "plot_vd/tet", 0);
    // free_vdiagram(vd_tet);

    puts("\nCUBE:");
    // generate_file_cube_bp(FILE_BP, 2);
    // s_vdiagram *vd_cube = construct_vd_from_txt(&r_fun, FILE_BP, 5);
    // if (!vd_cube) { puts("Could not construct vd in max_tries."); exit(1); }
    // check_volume(vd_cube);
    // FILE *f_vcells = fopen("test_vcells.txt", "w");
    // write_vd_file(vd_cube, f_vcells);
    // fclose(f_vcells);
    // plot_vdiagram_auto(vd_cube, "plot_vd/cube", 5);
    // free_vdiagram(vd_cube);

    puts("\nSPHERE:");
    // generate_file_sphere_bp(FILE_BP, 1.5, 15, 20);
    // s_vdiagram *vd_sph = construct_vd_from_txt(&r_fun, FILE_BP, 5);
    // if (!vd_sph) { puts("Could not construct vd in max_tries."); exit(1); }
    // check_volume(vd_sph);
    // fclose(f_vcells);
    // plot_vdiagram_auto(vd_sph, "plot_vd/sph", 5);
    // free_vdiagram(vd_sph);
    


    // ------------------ LEFT LUNG ---------------------
    puts("\nLEFT LUNG:");
    double **points_bp_L;
    int Np_L;
    s_bound_poly *bp_L;
    new_bpoly_from_txt("lobes/L.txt", &points_bp_L, &Np_L, &bp_L, 0);
    printf("volume: %f\n", bp_L->volume);
    printf("min: (%f, %f, %f)\n max: (%f, %f, %f)\n", bp_L->min[0], bp_L->min[1], bp_L->min[2], bp_L->max[0], bp_L->max[1], bp_L->max[2]);
    // plot_bpoly_with_points(bp_L, NULL, 0, "plot_vd/bp_L", NULL, "blue");

    
    int Nsimu = 10;
    // ADULT:
    s_bound_poly *bp_L_adult = scale_bpoly(bp_L, 1.06);
    z0 = bp_L_adult->min[2];
    zf = bp_L_adult->max[2];
    r_mean = 1.1;
    // PLOT
    s_vdiagram *vd_L = construct_vd_from_txt(&r_fun, "lobes/L.txt", 5);
    check_volume(vd_L);
    // plot_vdiagram_auto(vd_L, "plot_vd/L", 0);
    free_vdiagram(vd_L);
    // STATS
    generate_statistics(bp_L_adult, Nsimu, "volumes/L_adult.txt");
    PLOT_VOLUMES("volumes/L_adult");

    // 13 Y.O.:
    s_bound_poly *bp_L_13 = scale_bpoly(bp_L, 0.94);
    z0 = bp_L_13->min[2];
    zf = bp_L_13->max[2];
    r_mean = 1;
    generate_statistics(bp_L_13, Nsimu, "volumes/L_13yo.txt");
    PLOT_VOLUMES("volumes/L_13yo");

    // // 8 Y.O.:
    s_bound_poly *bp_L_8 = scale_bpoly(bp_L, 0.71);
    z0 = bp_L_8->min[2];
    zf = bp_L_8->max[2];
    r_mean = 0.9;
    generate_statistics(bp_L_8, Nsimu, "volumes/L_8yo.txt");
    PLOT_VOLUMES("volumes/L_8yo");

    // // 3 Y.O.:
    s_bound_poly *bp_L_3 = scale_bpoly(bp_L, 0.52);
    z0 = bp_L_3->min[2];
    zf = bp_L_3->max[2];
    r_mean = 0.7;
    generate_statistics(bp_L_3, Nsimu, "volumes/L_3yo.txt");
    PLOT_VOLUMES("volumes/L_3yo");





    puts("\nRIGHT LUNG");
    double **points_bp_R;
    int Np_R;
    s_bound_poly *bp_R;
    new_bpoly_from_txt("lobes/R.txt", &points_bp_R, &Np_R, &bp_R, 0);
    printf("volume: %f\n", bp_R->volume);
    printf("min: (%f, %f, %f)\n max: (%f, %f, %f)\n", bp_R->min[0], bp_R->min[1], bp_R->min[2], bp_R->max[0], bp_R->max[1], bp_R->max[2]);
    // plot_bpoly_with_points(bp_R, NULL, 0, "plot_vd/bp_R", NULL, "orange");

    // ADULT:
    s_bound_poly *bp_R_adult = scale_bpoly(bp_R, 1.06);
    z0 = bp_R_adult->min[2];
    zf = bp_R_adult->max[2];
    r_mean = 1.1;
    // PLOT
    s_vdiagram *vd_R = construct_vd_from_txt(&r_fun, "lobes/R.txt", 5);
    check_volume(vd_R);
    // plot_vdiagram_auto(vd_R, "plot_vd/R", 0);
    free_vdiagram(vd_R);
    // STATS
    generate_statistics(bp_L_adult, Nsimu, "volumes/R_adult.txt");
    PLOT_VOLUMES("volumes/R_adult");

    // 13 Y.O.:
    s_bound_poly *bp_R_13 = scale_bpoly(bp_R, 0.94);
    z0 = bp_R_13->min[2];
    zf = bp_R_13->max[2];
    r_mean = 1;
    generate_statistics(bp_R_13, Nsimu, "volumes/R_13yo.txt");
    PLOT_VOLUMES("volumes/R_13yo");

    // 8 Y.O.:
    s_bound_poly *bp_R_8 = scale_bpoly(bp_R, 0.72);
    z0 = bp_R_8->min[2];
    zf = bp_R_8->max[2];
    r_mean = 0.9;
    generate_statistics(bp_R_8, Nsimu, "volumes/R_8yo.txt");
    PLOT_VOLUMES("volumes/R_8yo");

    // // 3 Y.O.:
    s_bound_poly *bp_R_3 = scale_bpoly(bp_R, 0.52);
    z0 = bp_R_3->min[2];
    zf = bp_R_3->max[2];
    r_mean = 0.7;
    generate_statistics(bp_R_3, Nsimu, "volumes/R_3yo.txt");
    PLOT_VOLUMES("volumes/R_3yo");



    exit(1);


    
    // // DO A LOOP TO GENERATE A LOT OF LUNGS
    // int N_simu = 10;
    // for (int ii=0; ii<N_simu; ii++) {
    //     printf("%d\n", ii);
    //     bp = new_bpoly_copy(bp_L);
    //     vd = construct_cells_nonuniform(bp);
    //     append_volumes_to_file(vd, FILE_VOLS_SPHERE);
    //     check_volume(bp, vd);
    //     free_vdiagram(vd);
    // }

    // PLOTTING HISTOGRAM
    system("gnuplot sphere.plt");
}
