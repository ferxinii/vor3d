#include "vor3d.h"
#include "simplical_complex.h"
#include "delaunay.h"
#include "vdiagram.h"
#include "bpoly.h"
#include "geometry.h"
#include <stdlib.h>
#include <assert.h>
#include <math.h>


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
    // printf("rel_diff = %.16f\n", relative_diff);
    if (fabs(relative_diff) < 1e-3) return 1;
    else return 0;
}



// GENERALISED SEED HANDLING

typedef struct {
    int Ns;
    double **seeds;
} s_seed_set;

typedef union {
    s_seed_set seed_set;
    double (*f_radius_poiss)(double *);
} seed_userdata;

typedef s_seed_set (*f_seed_generator)(s_bound_poly *bp, seed_userdata ud);

static s_seed_set fixed_generator(s_bound_poly *bp, seed_userdata ud) {
    (void)bp;
    return ud.seed_set;
}

static s_seed_set PDS_generator(s_bound_poly *bp, seed_userdata ud) {
    int Ns;
    double **seeds = generate_nonuniform_poisson_dist_inside(bp, ud.f_radius_poiss, &Ns);
    return (s_seed_set){Ns, seeds};
}

static s_vdiagram *vor3d_core(s_bound_poly *bp, int max_tries, f_seed_generator f_seeds, seed_userdata ud)
{
    s_bound_poly *bp_tmp = new_bpoly_copy(bp);
    for (int ii = 0; ii < max_tries; ii++) {
        s_seed_set seed_set = f_seeds(bp_tmp, ud);
        double **seeds = seed_set.seeds;
        int Ns = seed_set.Ns;
        int Ns_extended = extend_sites_mirroring(bp_tmp, &seeds, Ns);

        s_setup *dt = construct_dt_3d(seeds, Ns_extended);

        s_vdiagram *vd = voronoi_from_delaunay_3d(dt, bp_tmp, Ns);
        if (!vd) {  // Retry
            free_bpoly(bp_tmp);
            bp_tmp = new_bpoly_copy(bp);
            continue;
        }

        if (valid_volumes(bp, vd)) {  // Return, with original bp
            free_bpoly(bp_tmp);
            vd->bpoly = bp;
            return vd;  
        } else {
            free_vdiagram(vd);
            bp_tmp = new_bpoly_copy(bp);
        }
    }
    free_bpoly(bp_tmp);
    return NULL;
}


s_vdiagram *vor3d_from_txt(double **seeds, int Ns, char *file_bounding_polyhedron, int max_tries)
{   
    s_bound_poly *bp = new_bpoly_from_txt(file_bounding_polyhedron);
    if (bp->Nf == 0) { free_bpoly(bp); puts("Error: bp->Nf == 0..."); exit(1); }
    
    seed_userdata ud = {.seed_set = {Ns, seeds}};
    return vor3d_core(bp, max_tries, &fixed_generator, ud);
}


s_vdiagram *vor3d_from_txt_PDS(double (*f_radius_poiss)(double *), char *file_bounding_polyhedron, int max_tries)
{   
    s_bound_poly *bp = new_bpoly_from_txt(file_bounding_polyhedron);
    if (bp->Nf == 0) { free_bpoly(bp); puts("Error: bp->Nf == 0..."); exit(1); }
    
    seed_userdata ud = {.f_radius_poiss = f_radius_poiss};
    return vor3d_core(bp, max_tries, &PDS_generator, ud);
}


s_vdiagram *vor3d_from_bp(double **seeds, int Ns, s_bound_poly *bp, int max_tries)
{
    seed_userdata ud = {.seed_set = {Ns, seeds}};
    return vor3d_core(bp, max_tries, &fixed_generator, ud);
}


s_vdiagram *vor3d_from_bp_PDS(double (*f_radius_poiss)(double *), s_bound_poly *bp, int max_tries)
{   
    seed_userdata ud = {.f_radius_poiss = f_radius_poiss};
    return vor3d_core(bp, max_tries, &PDS_generator, ud);
}

