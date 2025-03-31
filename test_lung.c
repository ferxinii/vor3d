

#include <stdio.h>
#include <time.h>
#include "geometry.c"
#include "simplical_complex.c"
#include "dt_3d_incremental.c"
#include "vd_3d.c"
#include "bpoly.c"


int main(void) {
    srand(time(NULL));

    double rmax = 2;

    int Np_LS;
    double **LS_p;
    s_bound_poly *bp_LS;
    new_bpoly_from_txt("lobes/R.txt", &LS_p, &Np_LS, &bp_LS);
    printf("\n\n\n --- READING FROM TXT, NP R: %d\n", Np_LS);

    int Np_in_LS;
    double **p_in_LS = generate_poisson_dist_inside(bp_LS, rmax, &Np_in_LS);

    printf("\n\n\n --- POISSON INSIDE LS, NP IN R: %d\n", Np_in_LS);
    puts("Now constructing dt");
    s_setup *dt_LS = construct_dt_3d(p_in_LS, Np_in_LS);

    puts("--- NOW CONSTRUCTING VD");
    s_vdiagram *vd_LS = voronoi_from_delaunay_3d(dt_LS, bp_LS);

    puts("--- NOW PLOTTING VD");
    double ranges_plot[6];
    ranges_plot[0] = bp_LS->min[0];     ranges_plot[1] = bp_LS->max[0];
    ranges_plot[2] = bp_LS->min[1];     ranges_plot[3] = bp_LS->max[1];
    ranges_plot[4] = bp_LS->min[2];     ranges_plot[5] = bp_LS->max[2];

    plot_vdiagram(vd_LS, "plot_lobes/LS", ranges_plot);
}
