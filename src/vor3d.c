#include "vor3d.h"
#include "scplx.h"
#include "vdiagram.h"
#include "points.h"
#include "./hlibs/dynarray.h"

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

typedef s_points (*f_seed_generator)(const s_bpoly *bp, seed_userdata ud, int *out_kept_idx);

/*
 * out_kept_idx, if non-NULL, must have ud.generator.seeds.N entries: filled
 * with the kept seed's index in the returned s_points (parallel to the
 * eventual vd.seeds/vd.vcells), or -1 for a seed dropped here (not strictly
 * interior to bp).
 */
static s_points fixed_generator(const s_bpoly *bp, seed_userdata ud, int *out_kept_idx)
{
    s_point *seed_copy = malloc(ud.generator.seeds.N * sizeof(s_point));

    int count = 0;
    for (int ii = 0; ii < ud.generator.seeds.N; ii++) {
        if (test_point_in_convhull(&bp->convh, ud.generator.seeds.p[ii], ud.EPS_DEG, ud.TOL) == TEST_IN) {
            seed_copy[count] = ud.generator.seeds.p[ii];
            if (out_kept_idx) out_kept_idx[ii] = count;
            count++;
        } else if (out_kept_idx) {
            out_kept_idx[ii] = -1;
        }
    }
    seed_copy = realloc(seed_copy, count * sizeof(s_point));
    return (s_points){count, seed_copy};
}

/* PDS-generated seeds have no caller-supplied identity to track, so
 * out_kept_idx is unused here -- callers of vor3d_inside_bp_PDS pass NULL. */
static s_points PDS_generator(const s_bpoly *bp, seed_userdata ud, int *out_kept_idx) {
    (void)out_kept_idx;
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



/*
 * Single attempt: builds the diagram once (vor3d_convex_clip) and checks volume conservation -- no retries, no
 * TOL/EPS_DEG escalation. out_kept_idx: see vor3d.h.
 */
static s_vdiagram vor3d_core(const s_bpoly *bp, double vol_max_rel_diff,
                             f_seed_generator f_seeds, seed_userdata ud,
                             int *out_kept_idx)
{
    int n_kept_idx = out_kept_idx ? ud.generator.seeds.N : 0;
    s_points seeds = f_seeds(bp, ud, out_kept_idx);
    if (seeds.N <= 0) {
        free_points(&seeds);
        return (s_vdiagram){0};
    }

    /* weights always NULL here: vor3d_in_bp / vor3d_in_convh reject non-NULL
     * weights up front (see vor3d_in_bp). */
    s_vdiagram vd = vor3d_convex_clip(&seeds, NULL, bp, ud.EPS_DEG, ud.TOL,
                                      out_kept_idx, n_kept_idx);
    free_points(&seeds);
    if (vd.seeds.N == 0) {
        fprintf(stderr, "vor3d: convex build failed.\n");
        if (out_kept_idx)
            for (int i = 0; i < n_kept_idx; i++) out_kept_idx[i] = -1;
        return (s_vdiagram){0};
    }
    if (!valid_volumes(bp, &vd, vol_max_rel_diff)) {
        free_vdiagram(&vd);
        if (out_kept_idx)
            for (int i = 0; i < n_kept_idx; i++) out_kept_idx[i] = -1;
        return (s_vdiagram){0};
    }
    return vd;
}



s_vdiagram vor3d_in_bp(const s_points *seeds, const double *weights,
                      const s_bpoly *bp, double vol_max_rel_diff,
                      double EPS_DEG, double TOL, int *out_kept_idx)
{
    /* Not wired through yet: fixed_generator drops non-interior seeds, so the
     * weights would need the same filtering. Reject loudly rather than build an
     * unweighted diagram -- weighted convex builds use vor3d_convex_clip. */
    if (weights) {
        fprintf(stderr, "vor3d_in_bp: non-NULL weights are not supported here; "
                        "use vor3d_convex_clip directly.\n");
        return (s_vdiagram){0};
    }
    /* Shallow view of the caller's seeds: fixed_generator only reads them and
     * returns its own filtered copy. (The former deep copy here was leaked.) */
    seed_userdata ud = {.generator.seeds = *seeds, .EPS_DEG = EPS_DEG, .TOL = TOL};
    return vor3d_core(bp, vol_max_rel_diff, &fixed_generator, ud, out_kept_idx);
}


s_vdiagram vor3d_in_convh(const s_points *seeds, const double *weights,
                          const s_convh *convh, double vol_max_rel_diff,
                          double EPS_DEG, double TOL, int *out_kept_idx)
{
    s_bpoly bp = bpoly_from_convh(convh);
    s_vdiagram out = vor3d_in_bp(seeds, weights, &bp, vol_max_rel_diff, EPS_DEG, TOL,
                                  out_kept_idx);
    free_bpoly(&bp);
    return out;
}


s_vdiagram vor3d_in_bp_pds(double (*f_radius_poiss)(double*, void*), void *f_params,
                           const s_bpoly *bp, double vol_max_rel_diff,
                           double EPS_DEG, double TOL,
                           int (*randint)(void*, int), double (*randd01)(void*),
                           void *rctx)
{
    seed_userdata ud = {.generator.pds = { .f_radius_poiss = f_radius_poiss, .f_params = f_params,
                                           .randint = randint, .randd01 = randd01, .rctx = rctx },
                        .EPS_DEG = EPS_DEG, .TOL = TOL};
    return vor3d_core(bp, vol_max_rel_diff, &PDS_generator, ud, NULL);
}




