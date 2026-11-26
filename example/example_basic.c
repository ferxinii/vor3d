#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <omp.h>
#include <math.h>
#include "vor3d.h"
#include "random.h"
#include "dynarray.h"

#define FILE_BP "bp_points.txt.aux"
#define PLOT_VOLUMES(name) system("./plot_together.plt " name)

/* Seeds are uniform random points inside the domain. The count roughly matches
 * the density of the Poisson-disk sampler (radius 1) this example used before
 * PDS was removed from the library: ~1.25 seeds per unit volume. */
#define SEEDS_PER_UNIT_VOLUME 1.25

double EPS_degenerate = 1e-12, TOL = 1e-12;
s_random_context *rctx_omp;


/* N uniform random points inside dom (rejection sampling on its bounding box). */
static s_points random_seeds_inside(const s_vor3d_domain *dom, s_random_context *rng)
{
    int N = (int)lround(SEEDS_PER_UNIT_VOLUME * dom->volume);
    if (N < 1) N = 1;
    s_points s = { .N = 0, .p = malloc(sizeof(s_point) * (size_t)N) };
    s_point mn = dom->bpoly.min, mx = dom->bpoly.max;
    while (s.N < N) {
        s_point p = {{{ mn.x + (mx.x - mn.x) * random_uniform_double(rng),
                        mn.y + (mx.y - mn.y) * random_uniform_double(rng),
                        mn.z + (mx.z - mn.z) * random_uniform_double(rng) }}};
        if (vor3d_domain_contains(dom, p, EPS_degenerate, TOL)) s.p[s.N++] = p;
    }
    return s;
}

static s_vdiagram build_random(const s_vor3d_domain *dom, s_random_context *rng)
{
    s_points seeds = random_seeds_inside(dom, rng);
    s_vdiagram vd = vor3d_diagram(&seeds, NULL, dom, EPS_degenerate, TOL, rng,
                                  NULL, false, false);
    free_points(&seeds);
    return vd;
}

static void check_volume(const s_vdiagram *vd)
{
    for (int ii=0; ii<vd->seeds.N; ii++)
        if (vd->vcells[ii].volume <= 0)
            printf("VOL %d : %f\n", ii, vd->vcells[ii].volume);
    printf("%d cells, Vol = %f, rel_diff = %.16f\n",
           vd->seeds.N, vd->domain_volume, vd->vol_rel_err);
}

/* Domain from a csv of points: the convex hull of the points. */
static s_vor3d_domain domain_from_csv(const char *fname)
{
    s_bpoly bp = bpoly_from_csv(fname, EPS_degenerate);
    s_vor3d_domain dom = vor3d_domain_from_convh(&bp.convh);
    free_bpoly(&bp);
    return dom;
}

static void iterate(int iterations, bool print_progress, const s_vor3d_domain *dom)
{
    int fail = 0;
    #pragma omp parallel for reduction(+:fail)
    for (int ii=0; ii<iterations; ii++) {
        if (print_progress && ii%100 == 0) printf("%d / %d\n", ii, iterations);
        s_vdiagram vd = build_random(dom, &rctx_omp[omp_get_thread_num()]);
        if (!vdiagram_is_valid(&vd) || fabs(vd.vol_rel_err) > 1e-6) fail++;
        free_vdiagram(&vd);
    }

    printf("FAILED %d / %d\n", fail, iterations);
}

void generate_file_cube_bp(const char *filename, double length)
{
    double s = length / 2;
    
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "%f, %f, %f\n", -s, -s, -s);
    fprintf(fp, "%f, %f, %f\n", -s, -s, s);
    fprintf(fp, "%f, %f, %f\n", -s, s, -s);
    fprintf(fp, "%f, %f, %f\n", s, -s, -s);
    fprintf(fp, "%f, %f, %f\n", -s, s, s);
    fprintf(fp, "%f, %f, %f\n", s, -s, s);
    fprintf(fp, "%f, %f, %f\n", s, s, -s);
    fprintf(fp, "%f, %f, %f\n", s, s, s);
    fclose(fp);
}


void generate_file_tetrahedron_bp(const char *filename, double length)
{
    double s = length / 2;
    
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "%f, %f, %f\n", -s, -s, -s);
    fprintf(fp, "%f, %f, %f\n", -s, -s, s);
    fprintf(fp, "%f, %f, %f\n", -s, s, -s);
    fprintf(fp, "%f, %f, %f\n", s, s, s);
    fclose(fp);
}


void generate_file_sphere_bp(const char *filename, double radius, int nTheta, int nPhi)
{
    // ntheta: example 18; // Number of steps in the polar angle
    // nphi: example 36;   // Number of steps in the azimuthal angle
    FILE *fp = fopen(filename, "w");
    // fprintf(fp, "%d\n\n", 2 + nPhi * (nTheta-1));

    for (int i = 0; i <= nTheta; i++) {
        double theta = M_PI * i / nTheta;
        
        if (i == 0 || i == nTheta) { // If at a pole, compute and write the coordinate once.
            double x = 0.0;
            double y = 0.0;
            double z = radius * cos(theta);  // will be +radius or -radius
            fprintf(fp, "%f, %f, %f\n", x, y, z);
        } else {
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


int main(void)
{
    rctx_omp = malloc(sizeof(s_random_context) * omp_get_max_threads());
    random_initialize_threads(time(NULL), 1, omp_get_max_threads(), rctx_omp);

    
    bool PLOT = false;
    (void)PLOT;

    puts("\nTETRAHEDON:");
    generate_file_tetrahedron_bp(FILE_BP, 3);
    s_vor3d_domain dom_tet = domain_from_csv(FILE_BP);
    if (!vor3d_domain_is_valid(&dom_tet)) { puts("Could not build domain."); exit(1); }
    s_vdiagram vd_tet = build_random(&dom_tet, &rctx_omp[0]);
    if (!vdiagram_is_valid(&vd_tet)) { puts("Could not construct vd."); exit(1); }
    check_volume(&vd_tet);
    if (PLOT) plot_vdiagram_differentviews(&vd_tet, "plots/tet", NULL);
    free_vdiagram(&vd_tet);
    iterate(10000, false, &dom_tet);
    free_vor3d_domain(&dom_tet);

    puts("\nCUBE:");
    generate_file_cube_bp(FILE_BP, 2);
    s_vor3d_domain dom_cube = domain_from_csv(FILE_BP);
    if (!vor3d_domain_is_valid(&dom_cube)) { puts("Could not build domain."); exit(1); }
    s_vdiagram vd_cube = build_random(&dom_cube, &rctx_omp[0]);
    if (!vdiagram_is_valid(&vd_cube)) { puts("Could not construct vd."); exit(1); }
    check_volume(&vd_cube);
    if (PLOT) plot_vdiagram_differentviews(&vd_cube, "plots/cube", NULL);
    free_vdiagram(&vd_cube);
    iterate(10000, false, &dom_cube);
    free_vor3d_domain(&dom_cube);

    puts("\nSPHERE:");
    generate_file_sphere_bp(FILE_BP, 1.5, 15, 20);
    s_vor3d_domain dom_sph = domain_from_csv(FILE_BP);
    if (!vor3d_domain_is_valid(&dom_sph)) { puts("Could not build domain."); exit(1); }
    write_convhull_to_m(&dom_sph.bpoly.convh, "test_bp.m");
    s_vdiagram vd_sph = build_random(&dom_sph, &rctx_omp[0]);
    if (!vdiagram_is_valid(&vd_sph)) { puts("Could not construct vd."); exit(1); }
    check_volume(&vd_sph);
    if (PLOT) plot_vdiagram_differentviews(&vd_sph, "plots/sph", NULL);
    free_vdiagram(&vd_sph);
    iterate(100, true, &dom_sph);
    free_vor3d_domain(&dom_sph);

    free(rctx_omp);
}

