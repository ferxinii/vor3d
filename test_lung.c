

#include <stdio.h>
#include <time.h>
#include "simplical_complex.c"
#include "dt_3d_incremental.c"
#include "vd_3d.c"
#include "bpoly.c"


#define FILE_VOLS_1 "test_volumes.txt"
#define FILE_VOLS_2 "test_volumes_2.txt"


double r_fun(double *x)
{   
    double r0 = 0.4;
    double z0 = -30;
    double K = 0.04;
    // printf("DEBUG R_FUN: (%f,%f,%f)  r: %f\n", x[0], x[1], x[2], r0 + K * fabs(x[2] - z0));

    return r0 + K * fabs(x[2] - z0);
}


void check_volume(s_bound_poly *bp, s_vdiagram *vd)
{
    double sum_vol = 0;
    for (int ii=0; ii<vd->N_vcells; ii++) {
        sum_vol += vd->vcells[ii]->volume;
    }
    printf("V lung = %f, Diff = %f\n", bp->volume, bp->volume - sum_vol);

}


s_vdiagram *construct_lung_cells_nonuniform(s_bound_poly *bp, double rmax)
{
    int N_points_poiss;
    (void)rmax;
    // double **points_poiss = generate_uniform_poisson_dist_inside(bp, rmax, &N_points_poiss);
    double **points_poiss = generate_nonuniform_poisson_dist_inside(bp, &r_fun, &N_points_poiss);
    printf("NPOINTS: %d\n", N_points_poiss);

    // puts("Constructing dt...");
    s_setup *dt = construct_dt_3d(points_poiss, N_points_poiss);

    // puts("Constructing vd...");
    s_vdiagram *vd = voronoi_from_delaunay_3d(dt, bp);

    free_complex(dt);
    return vd;
}


s_vdiagram *construct_lung_cells_uniform(s_bound_poly *bp, double rmax)
{
    int N_points_poiss;
    double **points_poiss = generate_uniform_poisson_dist_inside(bp, rmax, &N_points_poiss);
    // double **points_poiss = generate_nonuniform_poisson_dist_inside(bp, &r_fun, &N_points_poiss);
    printf("NPOINTS: %d\n", N_points_poiss);

    // puts("Constructing dt...");
    s_setup *dt = construct_dt_3d(points_poiss, N_points_poiss);

    // puts("Constructing vd...");
    s_vdiagram *vd = voronoi_from_delaunay_3d(dt, bp);

    free_complex(dt);
    return vd;
}


int main(void) {
    srand(time(NULL));
    fopen(FILE_VOLS_1, "w");
    fopen(FILE_VOLS_2, "w");

    double rmax = 1;

    // READ BP AND STORE THEM IN A POINTER, THAT I WILL COPY EACH TIME I CREATE A NEW VD
    double **points_bp_L;
    int N_points_bp_L;
    s_bound_poly *bp_L;
    new_bpoly_from_txt("lobes/L.txt", &points_bp_L, &N_points_bp_L, &bp_L, 0); 

    double **points_bp_R;
    int N_points_bp_R;
    s_bound_poly *bp_R;
    new_bpoly_from_txt("lobes/R.txt", &points_bp_R, &N_points_bp_R, &bp_R, 0);
    

    s_bound_poly *bp;
    bp = new_bpoly_copy(bp_L);
    printf("BP_L: MIN: (%f, %f, %f)\n", bp_L->min[0], bp_L->min[1], bp_L->min[2]);
    printf("      MAX: (%f, %f, %f)\n", bp_L->max[0], bp_L->max[1], bp_L->max[2]);
    s_vdiagram *vd = construct_lung_cells_nonuniform(bp, rmax);
    // append_volumes_to_file(vd, FILE_VOLS_1);
    // double ranges_plot[6];
    // ranges_plot[0] = bp->min[0];     ranges_plot[1] = bp->max[0];
    // ranges_plot[2] = bp->min[1];     ranges_plot[3] = bp->max[1];
    // ranges_plot[4] = bp->min[2];     ranges_plot[5] = bp->max[2];
    // plot_vdiagram(vd, "plot_lobes/L", ranges_plot, 0);
    free_vdiagram(vd);

    bp = new_bpoly_copy(bp_R);
    printf("BP_R: MIN: (%f, %f, %f)\n", bp_R->min[0], bp_R->min[1], bp_R->min[2]);
    printf("      MAX: (%f, %f, %f)\n", bp_R->max[0], bp_R->max[1], bp_R->max[2]);
    vd = construct_lung_cells_nonuniform(bp, rmax);
    // append_volumes_to_file(vd, FILE_VOLS_2);
    // double ranges_plot[6];
    // ranges_plot[0] = bp->min[0];     ranges_plot[1] = bp->max[0];
    // ranges_plot[2] = bp->min[1];     ranges_plot[3] = bp->max[1];
    // ranges_plot[4] = bp->min[2];     ranges_plot[5] = bp->max[2];
    // plot_vdiagram(vd, "plot_lobes/L", ranges_plot, 0);
    free_vdiagram(vd);

    
    // DO A LOOP TO GENERATE A LOT OF LUNGS
    int N_simu = 10;
    for (int ii=0; ii<N_simu; ii++) {
        printf("%d\n", ii);
        bp = new_bpoly_copy(bp_L);
        vd = construct_lung_cells_nonuniform(bp, rmax);
        append_volumes_to_file(vd, FILE_VOLS_1);
        check_volume(bp, vd);
        free_vdiagram(vd);

        bp = new_bpoly_copy(bp_R);
        vd = construct_lung_cells_nonuniform(bp, rmax);
        append_volumes_to_file(vd, FILE_VOLS_1);
        check_volume(bp, vd);
        free_vdiagram(vd);
    }

    for (int ii=0; ii<N_simu; ii++) {
        printf("%d\n", ii);
        bp = new_bpoly_copy(bp_L);
        vd = construct_lung_cells_uniform(bp, rmax);
        append_volumes_to_file(vd, FILE_VOLS_2);
        check_volume(bp, vd);
        free_vdiagram(vd);

        bp = new_bpoly_copy(bp_R);
        vd = construct_lung_cells_uniform(bp, rmax);
        append_volumes_to_file(vd, FILE_VOLS_2);
        check_volume(bp, vd);
        free_vdiagram(vd);
    }

    // PLOTTING HISTOGRAM
    system("gnuplot hist_vol.plt");
}
