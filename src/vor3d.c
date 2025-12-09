#include "vor3d.h"
#include "scplx.h"
#include "delaunay.h"
#include "vdiagram.h"
#include "bpoly.h"
#include "points.h"
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>


static int valid_volumes(const s_bpoly *bp, const s_vdiagram *vd, double max_rel_diff)
{
    double sum_vol = 0;
    for (int ii=0; ii<vd->seeds.N; ii++) {
        assert(vd->vcells[ii].volume > 0);
        sum_vol += vd->vcells[ii].volume;
    }

    double relative_diff =  (bp->volume - sum_vol) / bp->volume;
    if (fabs(relative_diff) < max_rel_diff) return 1;
    else {
        fprintf(stderr, "vor3d: vdiagram invalid volume (max_rel_diff = %g):\n", max_rel_diff);
        fprintf(stderr, "       bp->volume = %g, sum_vol = %g, rel_diff = %g\n", bp->volume, sum_vol, relative_diff);
        return 0;
    }
}



// GENERALISED SEED HANDLING
typedef struct {
    double (*f_radius_poiss)(double*, void*);
    void *f_params;
}   PDS_userdata;

typedef union {
    s_points seeds;
    PDS_userdata pds;
} u_userdata;

typedef struct {
    u_userdata generator;
    double EPS_degenerate;
    double TOL;
} seed_userdata;

typedef s_points (*f_seed_generator)(const s_bpoly *bp, seed_userdata ud);

static s_points fixed_generator(const s_bpoly *bp, seed_userdata ud) 
{
    s_point *seed_copy = malloc(ud.generator.seeds.N * sizeof(s_point));

    int count = 0;
    for (int ii = 0; ii < ud.generator.seeds.N; ii++) {
        if (test_point_in_convhull(&bp->convh, ud.generator.seeds.p[ii], ud.EPS_degenerate, ud.TOL) == TEST_IN) {
            seed_copy[count++] = ud.generator.seeds.p[ii];
        } else {
            // printf("Point %d not strictly inside\n", ii);
        }
    }
    seed_copy = realloc(seed_copy, count * sizeof(s_point));
    return (s_points){count, seed_copy};
}

static s_points PDS_generator(const s_bpoly *bp, seed_userdata ud) {
    double *bp_centroid = (double*)bp->CM.coords;
    if (isnan(ud.generator.pds.f_radius_poiss(bp_centroid, ud.generator.pds.f_params))) {  /* Check if single cell */
        s_point *p = malloc(sizeof(s_point)); 
        (*p).x = bp_centroid[0]; (*p).y = bp_centroid[1]; (*p).z = bp_centroid[2];
        return (s_points){ .N = 1, .p = p };
    }
    return generate_poisson_dist_inside(bp, ud.generator.pds.f_radius_poiss, ud.generator.pds.f_params, ud.EPS_degenerate);
}

static s_vdiagram vor3d_core(const s_bpoly *bp, double vol_max_rel_diff, int max_tries, f_seed_generator f_seeds, seed_userdata ud)
{
    s_vdiagram vd = {0};
    for (int ii = 0; ii < max_tries; ii++) {
        if (ii > 0) fprintf(stderr, "Retrying to build voronoi diagram. (%d / %d)\n", ii+1, max_tries);
        s_points seeds = f_seeds(bp, ud);
        int Nreal = seeds.N;

        extend_sites_mirroring(bp, ud.EPS_degenerate, &seeds);

        s_scplx dt = construct_dt_3d(&seeds, ud.TOL);
        vd = voronoi_from_delaunay_3d(&dt, bp, Nreal, ud.EPS_degenerate, ud.TOL);
        free_complex(&dt);

        if (vd.seeds.N == 0) {  /* Retry */
            fprintf(stderr, "vor3d: vdiagram is invalid.");
            continue;  
        }

        if (valid_volumes(bp, &vd, vol_max_rel_diff)) return vd;
        else free_vdiagram(&vd);
    }
    return vd;
}




s_vdiagram vor3d_from_bp(const s_points *seeds, const s_bpoly *bp, double vol_max_rel_diff, int max_tries, double EPS_degenerate, double TOL)
{
    seed_userdata ud = {.generator.seeds = copy_points(seeds), .EPS_degenerate = EPS_degenerate, .TOL = TOL};
    return vor3d_core(bp, vol_max_rel_diff, max_tries, &fixed_generator, ud);
}


s_vdiagram vor3d_from_bp_PDS(double (*f_radius_poiss)(double*, void*), void *f_params, const s_bpoly *bp, double vol_max_rel_diff, int max_tries,  double EPS_degenerate, double TOL)
{   
    seed_userdata ud = {.generator.pds = { .f_radius_poiss = f_radius_poiss, .f_params = f_params },
                        .EPS_degenerate = EPS_degenerate, .TOL = TOL};
    return vor3d_core(bp, vol_max_rel_diff, max_tries, &PDS_generator, ud);
}


s_vdiagram vor3d_from_txt(const s_points *seeds, char *file_bounding_polyhedron, double vol_max_rel_diff, int max_tries, double EPS_degenerate, double TOL)
{   
    s_bpoly bp = bpoly_from_csv(file_bounding_polyhedron, EPS_degenerate, TOL);
    if (bp.convh.Nf == 0) { puts("Error: bp->Nf == 0..."); exit(1); }
    
    s_vdiagram out = vor3d_from_bp(seeds, &bp, vol_max_rel_diff, max_tries, EPS_degenerate, TOL);

    free_bpoly(&bp);
    return out;
}


s_vdiagram vor3d_from_txt_PDS(double (*f_radius_poiss)(double*, void*), void *f_params, char *file_bounding_polyhedron, double vol_max_rel_diff, int max_tries, double EPS_degenerate, double TOL)
{   
    s_bpoly bp = bpoly_from_csv(file_bounding_polyhedron, EPS_degenerate, TOL);
    if (bp.convh.Nf == 0) { puts("Error: bp->Nf == 0..."); exit(1); }
    
    s_vdiagram out = vor3d_from_bp_PDS(f_radius_poiss, f_params, &bp, vol_max_rel_diff, max_tries, EPS_degenerate, TOL);

    free_bpoly(&bp);

    return out;
}


