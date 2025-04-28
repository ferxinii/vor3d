#include "voronoi.h"
#include "simplical_complex.h"
#include "dt_3d_incremental.h"
#include "vd_3d.h"
#include "bpoly.h"
#include "geometry.h"
#include <stdlib.h>
#include <assert.h>


// BP.C
extern void generate_file_cube_bp(const char *filename, double length);
extern void generate_file_tetrahedron_bp(const char *filename, double length);
extern void generate_file_sphere_bp(const char *filename, double radius, int nTheta, int nPhi);


// VD_3D.C
extern void clear_volumes_file(char *fname);
extern void append_volumes_to_file(s_vdiagram *vdiagram, char *fname);
extern void plot_vdiagram_auto(s_vdiagram *vdiagram, char *f_name, int max_files);


int valid_volumes(s_bound_poly *bp, s_vdiagram *vd)
{
    double sum_vol = 0;
    for (int ii=0; ii<vd->N_vcells; ii++) {
        if (vd->vcells[ii]->volume <= 0) {
            printf("Volume is not positive? ii = %d, Vol = %.16f\n", ii, vd->vcells[ii]->volume);
            printf("orient(vertices): %d\n", orientation(&vd->vcells[ii]->vertices[0], 
                                                         vd->vcells[ii]->vertices[3], 3));
            return 0;
        }
        sum_vol += vd->vcells[ii]->volume;
    }

    double relative_diff =  (bp->volume - sum_vol) / bp->volume;
    if (fabs(relative_diff) < 1e-3) return 1;
    else return 0;
}


s_vdiagram *construct_vd(double (*f_radius_poiss)(double *), char *file_bounding_polyhedron, int max_tries)
{
    for (int ii=0; ii<max_tries; ii++) {
        double **points_bp;
        int N_points_bp;
        s_bound_poly *bp;
        new_bpoly_from_txt(file_bounding_polyhedron, &points_bp, &N_points_bp, &bp, 0);

        int Ns;
        double **seeds = generate_nonuniform_poisson_dist_inside(bp, f_radius_poiss, &Ns);
        
        int Ns_extended = extend_sites_mirroring(bp, &seeds, Ns);
        
        s_setup *dt = construct_dt_3d(seeds, Ns_extended);
        
        s_vdiagram *vd = voronoi_from_delaunay_3d(dt, bp, Ns);
        if (!vd) continue;
        
        if (valid_volumes(bp, vd)) return vd;
        else free_vdiagram(vd);
    }
    return NULL;
}



