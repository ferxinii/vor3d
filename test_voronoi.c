#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "voronoi.h"
#include "geometry.h"


#define FILE_BP "bp_points.txt"
#define FILE_VOLS "volumes.txt"


double r_fun(double *x)
{   
    double r0 = 1.5;
    double z0 = -2;
    double K = 0.04;

    return r0 + K * fabs(x[2] - z0);
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


int main(void)
{
    srand(time(NULL));
    system("rm -f plot_vd/*");
    

    puts("\nTETRAHEDON:");
    generate_file_tetrahedron_bp(FILE_BP, 3);
    s_vdiagram *vd_tet = construct_vd(&r_fun, FILE_BP, 5);
    if (!vd_tet) { puts("Could not construct vd in max_tries."); exit(1); }
    check_volume(vd_tet);
    // plot_vdiagram_auto(vd_tet, "plot_vd/tet", 0);
    // free_vdiagram(vd_tet);

    puts("\nCUBE:");
    generate_file_cube_bp(FILE_BP, 2);
    s_vdiagram *vd_cube = construct_vd(&r_fun, FILE_BP, 5);
    if (!vd_cube) { puts("Could not construct vd in max_tries."); exit(1); }
    check_volume(vd_cube);
    // FILE *f_vcells = fopen("test_vcells.txt", "w");
    // write_vd_file(vd_cube, f_vcells);
    // fclose(f_vcells);
    // plot_vdiagram_auto(vd_cube, "plot_vd/cube", 5);
    // free_vdiagram(vd_cube);

    puts("\nSPHERE:");
    generate_file_sphere_bp(FILE_BP, 1.5, 15, 20);
    s_vdiagram *vd_sph = construct_vd(&r_fun, FILE_BP, 5);
    if (!vd_sph) { puts("Could not construct vd in max_tries."); exit(1); }
    check_volume(vd_sph);
    // FILE *f_vcells = fopen("test_vcells.txt", "w");
    // write_vd_file(vd_sph, f_vcells);
    // fclose(f_vcells);
    // plot_vdiagram_auto(vd_sph, "plot_vd/sph", 5);
    // free_vdiagram(vd_sph);

    puts("\nLEFT LUNG:");
    double **points_bp_L;
    int Np_L;
    s_bound_poly *bp_L;
    new_bpoly_from_txt("lobes/L.txt", &points_bp_L, &Np_L, &bp_L, 0);
    printf("volume: %f\n", compute_volume_convhull(bp_L->points, bp_L->faces, bp_L->fnormals, bp_L->Nf));
    // plot_bpoly_with_points(bp_L, NULL, 0, "plot_vd/bp_L", NULL, "blue");
    s_vdiagram *vd_L = construct_vd(&r_fun, "lobes/L.txt", 5);
    if (!vd_L) { puts("Could not construct vd in max_tries."); exit(1); }
    check_volume(vd_L);
    // FILE *f_vcells = fopen("test_vcells.txt", "w");
    // write_vd_file(vd_L, f_vcells);
    // fclose(f_vcells);
    // plot_vdiagram_auto(vd_L, "plot_vd/L", 0);
    // free_vdiagram(vd_L);


    puts("\nRIGHT LUNG");
    double **points_bp_R;
    int Np_R;
    s_bound_poly *bp_R;
    new_bpoly_from_txt("lobes/R.txt", &points_bp_R, &Np_R, &bp_R, 0);
    printf("volume: %f\n", compute_volume_convhull(bp_R->points, bp_R->faces, bp_R->fnormals, bp_R->Nf));
    // plot_bpoly_with_points(bp_R, NULL, 0, "plot_vd/bp_R", NULL, "orange");
    s_vdiagram *vd_R = construct_vd(&r_fun, "lobes/L.txt", 5);
    if (!vd_R) { puts("Could not construct vd in max_tries."); exit(1); }
    check_volume(vd_R);
    // FILE *f_vcells = fopen("test_vcells.txt", "w");
    // write_vd_file(vd_L, f_vcells);
    // fclose(f_vcells);
    plot_vdiagram_auto(vd_L, "plot_vd/L", 0);
    // free_vdiagram(vd_R);
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
