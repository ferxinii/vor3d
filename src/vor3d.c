#include "vor3d.h"
#include "simplical_complex.h"
#include "delaunay.h"
#include "vdiagram.h"
#include "bpoly.h"
#include "geometry.h"
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>


int valid_volumes(const s_bpoly *bp, const s_vdiagram *vd)
{
    double sum_vol = 0;
    for (int ii=0; ii<vd->seeds.N; ii++) {
        if (vd->vcells[ii].volume <= 0) {
            printf("Volume is not positive? ii = %d, Vol = %.16f\n", ii, vd->vcells[ii].volume);
            printf("orient(vertices): %d\n", orientation(&vd->vcells[ii].convh.points.p[0], 
                                                         vd->vcells[ii].convh.points.p[3]));
            return 0;
        }
        sum_vol += vd->vcells[ii].volume;
    }

    double relative_diff =  (bp->volume - sum_vol) / bp->volume;
    printf("bp->volume = %f, sum_vol = %f, rel_diff = %.16f\n", bp->volume, sum_vol, relative_diff);
    if (fabs(relative_diff) < 1e-3) return 1;
    else return 0;
}



// GENERALISED SEED HANDLING
typedef union {
    s_points seeds;
    double (*f_radius_poiss)(double *);
} seed_userdata;

typedef s_points (*f_seed_generator)(const s_bpoly *bp, seed_userdata ud);

static s_points fixed_generator(const s_bpoly *bp, seed_userdata ud) 
{
    s_point *seed_copy = malloc(ud.seeds.N * sizeof(s_point));

    int count = 0;
    for (int ii = 0; ii < ud.seeds.N; ii++) {
        if (is_inside_convhull(&bp->convh, ud.seeds.p[ii]) == 1) {
            seed_copy[count++] = ud.seeds.p[ii];
        } else {
            // printf("Point %d not strictly inside\n", ii);
        }
    }
    seed_copy = realloc(seed_copy, count * sizeof(s_point));
    return (s_points){count, seed_copy};
}

static s_points PDS_generator(const s_bpoly *bp, seed_userdata ud) {
    return generate_poisson_dist_inside(bp, ud.f_radius_poiss);
}

static s_vdiagram vor3d_core(const s_bpoly *bp, int max_tries, f_seed_generator f_seeds, seed_userdata ud)
{
    s_vdiagram vd = {0};
    for (int ii = 0; ii < max_tries; ii++) {
        if (ii > 0) puts("Retrying to build voronoi diagram.");
        s_points seeds = f_seeds(bp, ud);
        int Nreal = seeds.N;

        extend_sites_mirroring(bp, &seeds);

        s_scplx dt = construct_dt_3d(&seeds);

        vd = voronoi_from_delaunay_3d(&dt, bp, Nreal);
        if (vd.seeds.N == 0) continue;  // Retry

        if (valid_volumes(bp, &vd)) return vd;
        else free_vdiagram(&vd);
    }
    return vd;
}




s_vdiagram vor3d_from_bp(const s_points *seeds, const s_bpoly *bp, int max_tries)
{
    seed_userdata ud = {.seeds = copy_points(seeds)};
    return vor3d_core(bp, max_tries, &fixed_generator, ud);
}


s_vdiagram vor3d_from_bp_PDS(double (*f_radius_poiss)(double *), const s_bpoly *bp, int max_tries)
{   
    seed_userdata ud = {.f_radius_poiss = f_radius_poiss};
    return vor3d_core(bp, max_tries, &PDS_generator, ud);
}


s_vdiagram vor3d_from_txt(const s_points *seeds, char *file_bounding_polyhedron, int max_tries)
{   
    s_bpoly bp = bpoly_from_csv(file_bounding_polyhedron);
    if (bp.convh.Nf == 0) { puts("Error: bp->Nf == 0..."); exit(1); }
    
    s_vdiagram out = vor3d_from_bp(seeds, &bp, max_tries);

    free_bpoly(&bp);
    return out;
}


s_vdiagram vor3d_from_txt_PDS(double (*f_radius_poiss)(double *), char *file_bounding_polyhedron, int max_tries)
{   
    s_bpoly bp = bpoly_from_csv(file_bounding_polyhedron);
    if (bp.convh.Nf == 0) { puts("Error: bp->Nf == 0..."); exit(1); }
    
    s_vdiagram out = vor3d_from_bp_PDS(f_radius_poiss, &bp, max_tries);

    free_bpoly(&bp);

    return out;
}


