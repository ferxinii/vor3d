#include <stdio.h>
#include <time.h>
#include "simplical_complex.c"
#include "dt_3d_incremental.c"
#include "vd_3d.c"
#include "bpoly.c"


#define FILE_COORDS_SPHERE "sphere_coords.txt"
#define FILE_VOLS_SPHERE "sphere_volumes.txt"


double r_fun(double *x)
{   
    double r0 = 1;
    double z0 = -2;
    double K = 0.04;
    // printf("DEBUG R_FUN: (%f,%f,%f)  r: %f\n", x[0], x[1], x[2], r0 + K * fabs(x[2] - z0));

    return r0 + K * fabs(x[2] - z0);
}


void write_sphere_txt(void)
{
    // Parameters for the sphere
    double radius = 2;
    int nTheta = 18; // Number of steps in the polar angle (θ)
    int nPhi = 36;   // Number of steps in the azimuthal angle (φ)

    // Open the file for writing coordinates
    FILE *fp = fopen(FILE_COORDS_SPHERE, "w");
    fprintf(fp, "%d\n\n", 2 + nPhi * (nTheta-1));

    // Generate coordinates on the surface of the sphere
    // θ: angle from the positive z-axis (0 to PI)
    // φ: angle in the x-y plane (0 to 2*PI)
        // Loop over theta (0 to PI)
    for (int i = 0; i <= nTheta; i++) {
        double theta = M_PI * i / nTheta;
        
        // Check if we are at a pole (theta = 0 or PI)
        // If at a pole, compute and write the coordinate once.
        if (i == 0 || i == nTheta) {
            double x = 0.0;
            double y = 0.0;
            double z = radius * cos(theta);  // will be +radius or -radius
            fprintf(fp, "%f %f %f\n", x, y, z);
        } else {
            // If not a pole, loop over φ (0 to 2*PI) normally
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




void check_volume(s_bound_poly *bp, s_vdiagram *vd)
{
    double sum_vol = 0;
    for (int ii=0; ii<vd->N_vcells; ii++) {
        sum_vol += vd->vcells[ii]->volume;
        if (vd->vcells[ii]->volume <= 0) {
            printf("VOL %d : %f\n", ii, vd->vcells[ii]->volume);
        }
    }
    printf("V lung = %f, Diff = %f\n", bp->volume, bp->volume - sum_vol);

}


s_vdiagram *construct_cells_nonuniform(s_bound_poly *bp)
{
    int N_points_poiss;
    puts("Now generating poisson inside, nonuniform");
    double **points_poiss = generate_nonuniform_poisson_dist_inside(bp, &r_fun, &N_points_poiss);
    printf("NPOINTS: %d\n", N_points_poiss);

    // puts("Constructing dt...");
    s_setup *dt = construct_dt_3d(points_poiss, N_points_poiss);

    double sum_vol = 0;
    s_ncell *current = dt->head;
    int it = 0;
    while (current) {
        sum_vol += current->volume;
        if (current->volume < 0) {
            printf("vol dt %d: %f\n", it++, current->volume);
        }
        current = current->next;
    }


    printf("DT, SUM VOLUMES = %f, CONVHULL VOLUME = %f\n", sum_vol, compute_volume_complex(dt));
    // double ranges_plot[6];
    // ranges_plot[0] = bp->min[0];     ranges_plot[1] = bp->max[0];
    // ranges_plot[2] = bp->min[1];     ranges_plot[3] = bp->max[1];
    // ranges_plot[4] = bp->min[2];     ranges_plot[5] = bp->max[2];
    // plot_dt_3d(dt, "test_dt.png", ranges_plot);

    // puts("Constructing vd...");
    s_vdiagram *vd = voronoi_from_delaunay_3d(dt, bp);

    free_complex(dt);
    return vd;
}


int main(void) {
    srand(time(NULL));
    fopen(FILE_VOLS_SPHERE, "w");

    // READ BP AND STORE THEM IN A POINTER, THAT I WILL COPY EACH TIME I CREATE A NEW VD
    double **points_bp_L;
    int N_points_bp_L;
    s_bound_poly *bp_L;
    write_sphere_txt();
    new_bpoly_from_txt(FILE_COORDS_SPHERE, &points_bp_L, &N_points_bp_L, &bp_L, 0); 

    

    s_bound_poly *bp;
    bp = new_bpoly_copy(bp_L);
    printf("BP_L: MIN: (%f, %f, %f)\n", bp_L->min[0], bp_L->min[1], bp_L->min[2]);
    printf("      MAX: (%f, %f, %f)\n", bp_L->max[0], bp_L->max[1], bp_L->max[2]);
    printf("      Np = %d\n", bp->Np);
    // for (int ii=0; ii<bp->Np; ii++) {
    //     printf("(%f, %f, %f)\n", bp->points[ii][0], bp->points[ii][1], bp->points[ii][2]);
    // }
    s_vdiagram *vd = construct_cells_nonuniform(bp);

    append_volumes_to_file(vd, FILE_VOLS_SPHERE);
    printf("RESULT: N_vcells = %d\n", vd->N_vcells);

    check_volume(bp, vd);


    // TESTING FOR LOST VOLUME???
    FILE *f_vcells = fopen("test_vcells.txt", "w");
    write_vd_file(vd, f_vcells);
    fclose(f_vcells);



    FILE *ftest = fopen("lost_volume.txt", "w");
    int Ntest = 1000;
    double **ptest = malloc_matrix(Ntest, 3);
    for (int ii=0; ii<Ntest; ii++) {
        double x[3]; 
        random_point_inside_convhull(bp->points, bp->faces, bp->fnormals, bp->Nf, bp->min, bp->max, x);
        while (find_inside_which_vcell(vd, x) != -1) {
            random_point_inside_convhull(bp->points, bp->faces, bp->fnormals, bp->Nf, bp->min, bp->max, x);
        }
        ptest[ii][0] = x[0];
        ptest[ii][1] = x[1];
        ptest[ii][2] = x[2];
        fprintf(ftest, "%f, %f, %f\n", x[0], x[1], x[2]);
    }
    fclose(ftest);


    // double ranges_plot[6];
    // ranges_plot[0] = bp->min[0];     ranges_plot[1] = bp->max[0];
    // ranges_plot[2] = bp->min[1];     ranges_plot[3] = bp->max[1];
    // ranges_plot[4] = bp->min[2];     ranges_plot[5] = bp->max[2];
    // puts("PLOTTING...");
    // plot_vdiagram(vd, "plot_sphere/sph", ranges_plot, 0, ptest, &Ntest);



    exit(1);


    free_vdiagram(vd);


    
    // DO A LOOP TO GENERATE A LOT OF LUNGS
    int N_simu = 10;
    for (int ii=0; ii<N_simu; ii++) {
        printf("%d\n", ii);
        bp = new_bpoly_copy(bp_L);
        vd = construct_cells_nonuniform(bp);
        append_volumes_to_file(vd, FILE_VOLS_SPHERE);
        check_volume(bp, vd);
        free_vdiagram(vd);
    }

    // PLOTTING HISTOGRAM
    system("gnuplot sphere.plt");
}
