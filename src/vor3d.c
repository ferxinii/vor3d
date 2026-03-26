#include "vor3d.h"
#include "scplx.h"
#include "delaunay.h"
#include "vdiagram.h"
#include "points.h"
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>


static int valid_volumes(const s_bpoly *bp, const s_vdiagram *vd, double max_rel_diff)
{
    double sum_vol = 0;
    for (int ii=0; ii<vd->seeds.N; ii++) {
        if (vd->vcells[ii].volume <= 0) {
            fprintf(stderr, "vcell %d has zero or negative volume: %g\n", ii, vd->vcells[ii].volume);
            return 0;
        }
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
    int (*randint)(void*, int);
    double (*randd01)(void*);
    void *rctx;
}   PDS_userdata;

typedef union {
    s_points seeds;
    PDS_userdata pds;
} u_userdata;

typedef struct {
    u_userdata generator;
    double EPS_DEG;
    double TOL;
} seed_userdata;

typedef s_points (*f_seed_generator)(const s_bpoly *bp, seed_userdata ud);

static s_points fixed_generator(const s_bpoly *bp, seed_userdata ud) 
{
    s_point *seed_copy = malloc(ud.generator.seeds.N * sizeof(s_point));

    int count = 0;
    for (int ii = 0; ii < ud.generator.seeds.N; ii++) {
        if (test_point_in_convhull(&bp->convh, ud.generator.seeds.p[ii], ud.EPS_DEG, ud.TOL) == TEST_IN)
            seed_copy[count++] = ud.generator.seeds.p[ii];
    }
    seed_copy = realloc(seed_copy, count * sizeof(s_point));
    return (s_points){count, seed_copy};
}

static s_points PDS_generator(const s_bpoly *bp, seed_userdata ud) {
    double *bp_centroid = (double*)bp->CM.coords;
    /* Check if function returns NAN. If so, generate a single seed / cell */
    if (isnan(ud.generator.pds.f_radius_poiss(bp_centroid, ud.generator.pds.f_params))) {  
        s_point *p = malloc(sizeof(s_point)); 
        (*p).x = bp_centroid[0]; (*p).y = bp_centroid[1]; (*p).z = bp_centroid[2];
        return (s_points){ .N = 1, .p = p };
    }
    return generate_poisson_dist_inside(bp, ud.generator.pds.f_radius_poiss, ud.generator.pds.f_params,
                                        ud.generator.pds.randd01, ud.generator.pds.randint, 
                                        ud.generator.pds.rctx, ud.EPS_DEG);
}



static s_vdiagram vor3d_core(const s_bpoly *bp, double vol_max_rel_diff, int max_tries, 
                             f_seed_generator f_seeds, seed_userdata ud, int (*randint)(void*, int),
                             void *rctx, s_dynarray *buff_points, s_dynarray *buff_3dbls)
{
    s_vdiagram vd = {0};
    for (int ii = 0; ii < max_tries; ii++) {
        if (ii > 0) { 
            fprintf(stderr, "Retrying to build voronoi diagram. (%d / %d).\n", ii+1, max_tries);
            fprintf(stderr, "Incrementing TOL and EPS by a factor of 10.\n");
            ud.TOL *= 10;  ud.EPS_DEG *= 10; 
        }
        s_points seeds = f_seeds(bp, ud);
        int Nreal = seeds.N;

        if (!extend_sites_mirroring(bp, ud.EPS_DEG, ud.TOL, &seeds, randint, rctx, 
                                    buff_points, buff_3dbls)) {
            fprintf(stderr, "vor3d: error extending sites (initial).");
            free_points(&seeds);
        }

        s_scplx dt = construct_dt_3d(&seeds, NULL, false, ud.TOL);
        vd = voronoi_from_delaunay_3d(&dt, bp, Nreal, ud.EPS_DEG, ud.TOL);
        free_complex(&dt);
        if (vd.seeds.N == 0) {  /* Retry */
            fprintf(stderr, "vor3d: vdiagram is invalid.");
            free_points(&seeds);
            continue;  
        }
        
        if (valid_volumes(bp, &vd, vol_max_rel_diff)) {
            free_points(&seeds);
            return vd;;
        }
        // else {
        //     write_points_to_csv("seeds.csv", "w", &seeds);
        //     for (int ii=0; ii<vd.seeds.N; ii++) {
        //         char buff[256];
        //         snprintf(buff, 256, "v%d.m", ii);
        //         write_convhull_to_m(&vd.vcells[ii].convh, buff);
        //     }
        //     write_convhull_to_m(&vd.bpoly.convh, "bp.m");
        //     // return (s_vdiagram){0};
        //     exit(1);
        // }
        free_points(&seeds);
        free_vdiagram(&vd);
    }

    return (s_vdiagram){0};
}




s_vdiagram vor3d_from_bp(const s_points *seeds, const s_bpoly *bp, double vol_max_rel_diff,
                         int max_tries, double EPS_DEG, double TOL, int (*randint)(void*, int),
                         void *rctx, s_dynarray *buff_points, s_dynarray *buff_3dbls)
{
    seed_userdata ud = {.generator.seeds = copy_points(seeds),
                        .EPS_DEG = EPS_DEG, .TOL = TOL};
    return vor3d_core(bp, vol_max_rel_diff, max_tries, &fixed_generator, ud, 
                      randint, rctx, buff_points, buff_3dbls);
}


s_vdiagram vor3d_from_bp_PDS(double (*f_radius_poiss)(double*, void*), void *f_params,
                             const s_bpoly *bp, double vol_max_rel_diff, int max_tries,
                             double EPS_DEG, double TOL,
                             int (*randint)(void*, int), double (*randd01)(void*),
                             void *rctx, s_dynarray *buff_points, s_dynarray *buff_3dbls)
{   
    seed_userdata ud = {.generator.pds = { .f_radius_poiss = f_radius_poiss, .f_params = f_params,
                                           .randint = randint, .randd01 = randd01, .rctx = rctx },
                        .EPS_DEG = EPS_DEG, .TOL = TOL};
    return vor3d_core(bp, vol_max_rel_diff, max_tries, &PDS_generator, ud,
                      randint, rctx, buff_points, buff_3dbls);
}


s_vdiagram vor3d_from_txt(const s_points *seeds, char *file_bounding_polyhedron,
                          double vol_max_rel_diff, int max_tries, double EPS_DEG, double TOL,
                          int (*randint)(void*, int), void *rctx,
                          s_dynarray *buff_points, s_dynarray *buff_3dbls)
{   
    s_bpoly bp = bpoly_from_csv(file_bounding_polyhedron, EPS_DEG, TOL);
    if (bp.convh.Nf == 0) { puts("Error: bp->Nf == 0..."); exit(1); }
    
    s_vdiagram out = vor3d_from_bp(seeds, &bp, vol_max_rel_diff, max_tries, EPS_DEG, TOL,
                                   randint, rctx, buff_points, buff_3dbls);

    free_bpoly(&bp);
    return out;
}


s_vdiagram vor3d_from_txt_PDS(double (*f_radius_poiss)(double*, void*), void *f_params,
                              char *file_bounding_polyhedron, double vol_max_rel_diff,
                              int max_tries, double EPS_DEG, double TOL,
                              int (*randint)(void*, int), double (*randd01)(void*), void *rctx,
                              s_dynarray *buff_points, s_dynarray *buff_3dbls)
{   
    s_bpoly bp = bpoly_from_csv(file_bounding_polyhedron, EPS_DEG, TOL);
    if (bp.convh.Nf == 0) { puts("Error: bp->Nf == 0..."); exit(1); }
    
    s_vdiagram out = vor3d_from_bp_PDS(f_radius_poiss, f_params, &bp, vol_max_rel_diff,
                                       max_tries, EPS_DEG, TOL,
                                       randint, randd01, rctx, buff_points, buff_3dbls);

    free_bpoly(&bp);

    return out;
}


