#include "volsph.h"
#include "dynarray.h"
#include "points.h"
#include "scplx.h"
#include "delaunay.h"
#include <stdio.h>
#include <math.h>
#include <time.h>

#define TOL     1e-6    /* tolerance for exact (analytic) comparisons */
#define EPS     1e-10
#define TOL_DUP 1e-12

static int n_passed = 0;
static int n_failed = 0;

static void check(const char *name, double got, double expected, double tol)
{
    double err = fabs(got - expected);
    double rel = fabs(expected) > 1e-14 ? err / fabs(expected) : err;
    if (rel < tol) {
        printf("  PASS  %-34s got %.8f  expected %.8f\n", name, got, expected);
        n_passed++;
    } else {
        printf("  FAIL  %-34s got %.8f  expected %.8f  rel %.2e\n",
               name, got, expected, rel);
        n_failed++;
    }
}

/* Monte-Carlo estimate of the union volume, for validation only. */
static double mc_volume_union_spheres(const s_points *centers, const double *radii,
                                      int N_samples, unsigned int seed)
{
    double lo[3], hi[3];
    for (int k = 0; k < 3; k++) { lo[k] = 1e300; hi[k] = -1e300; }
    for (int i = 0; i < centers->N; i++)
        for (int k = 0; k < 3; k++) {
            lo[k] = fmin(lo[k], centers->p[i].coords[k] - radii[i]);
            hi[k] = fmax(hi[k], centers->p[i].coords[k] + radii[i]);
        }
    double L[3], box_vol = 1.0;
    for (int k = 0; k < 3; k++) { L[k] = hi[k] - lo[k]; box_vol *= L[k]; }

    srand(seed);
    long hits = 0;
    for (int s = 0; s < N_samples; s++) {
        double x[3];
        for (int k = 0; k < 3; k++) x[k] = lo[k] + L[k] * ((double)rand() / RAND_MAX);
        for (int i = 0; i < centers->N; i++) {
            double d2 = 0.0;
            for (int k = 0; k < 3; k++) {
                double dk = x[k] - centers->p[i].coords[k];
                d2 += dk * dk;
            }
            if (d2 <= radii[i] * radii[i]) { hits++; break; }
        }
    }
    return box_vol * ((double)hits / N_samples);
}

static void check_mc(const char *name, double exact, double mc, double mc_tol)
{
    double rel = fabs(exact - mc) / fmax(fabs(exact), 1e-10);
    if (rel < mc_tol) {
        printf("  PASS  %-34s exact %.6f  mc %.6f  rel %.2e\n", name, exact, mc, rel);
        n_passed++;
    } else {
        printf("  FAIL  %-34s exact %.6f  mc %.6f  rel %.2e\n", name, exact, mc, rel);
        n_failed++;
    }
}

int main(void)
{
    s_dynarray buff = dynarray_initialize(sizeof(s_ncell *), 64);
    const int    N_MC   = 4000000;   /* ~0.1-0.2% MC accuracy */
    const double MC_TOL = 0.015;     /* 1.5% relative (MC noise + boundary conventions) */
    const unsigned int SEED = 20260713u;

    printf("=== Analytic cases (N=1,2,3: geometric primitives, no weighted RT) ===\n");
    {   /* N=1: single sphere. */
        s_point c[1] = { {{{0, 0, 0}}} };
        double  r[1] = { 1.3 };
        s_points pts = { 1, c };
        check("N=1 single sphere",
              volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff), volume_sphere(1.3), TOL);
    }
    {   /* N=2: disjoint -> sum of volumes. */
        s_point c[2] = { {{{0,0,0}}}, {{{10,0,0}}} };
        double  r[2] = { 1.0, 1.2 };
        s_points pts = { 2, c };
        check("N=2 disjoint",
              volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff),
              volume_sphere(1.0) + volume_sphere(1.2), TOL);
    }
    {   /* N=2: containment -> volume of the big sphere. */
        s_point c[2] = { {{{0,0,0}}}, {{{0.1,0,0}}} };
        double  r[2] = { 2.0, 0.5 };
        s_points pts = { 2, c };
        check("N=2 containment",
              volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff), volume_sphere(2.0), TOL);
    }
    {   /* N=2: overlap -> matches the closed-form 2-ball union. */
        s_point c[2] = { {{{0,0,0}}}, {{{2,0,0}}} };
        double  r[2] = { 1.5, 1.1 };
        s_points pts = { 2, c };
        check("N=2 overlap vs closed form",
              volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff),
              volume_union_2_spheres(c[0], r[0], c[1], r[1], true), TOL);
    }

    printf("=== Monte-Carlo validation (dihedral / full pipeline) ===\n");
    {   /* N=3 mutually intersecting incl. triple lens -> exercises the
         * Cayley-Menger dihedral (the resolved duplicated block). */
        s_point c[3] = { {{{0,0,0}}}, {{{2,0,0}}}, {{{1.0, 1.7320508075688772, 0}}} };
        double  r[3] = { 1.5, 1.6, 1.9 };
        s_points pts = { 3, c };
        check_mc("N=3 overlapping (triple lens)",
                 volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff),
                 mc_volume_union_spheres(&pts, r, N_MC, SEED), MC_TOL);
    }

    printf("=== N>=4 (weighted RT + alpha complex) ===\n");
    {   /* N=4 disjoint -> sum; validates isolated balls contribute (Finding B). */
        s_point c[4] = { {{{0,0,0}}}, {{{10,0,0}}}, {{{0,10,0}}}, {{{0,0,10}}} };
        double  r[4] = { 1.0, 1.2, 0.8, 1.1 };
        s_points pts = { 4, c };
        double sum = 0.0; for (int i = 0; i < 4; i++) sum += volume_sphere(r[i]);
        check("N=4 disjoint (isolated balls)",
              volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff), sum, TOL);
    }
    {   /* N=4 overlapping on a tetrahedron. */
        s_point c[4] = { {{{1,1,1}}}, {{{1,-1,-0.5}}}, {{{-1,1,-1}}}, {{{-1,-1,1}}} };
        double  r[4] = { 2.4, 1.8, 1.9, 1.5 };
        s_points pts = { 4, c };
        check_mc("N=4 overlapping (tetrahedral)",
                 volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff),
                 mc_volume_union_spheres(&pts, r, N_MC, SEED), MC_TOL);
    }
    {   /* N=5 partial overlap. */
        s_point c[5] = { {{{0,0,0}}}, {{{1.2,0,0}}}, {{{0.6,1.0,0}}},
                         {{{0,0,1.2}}}, {{{1.2,0,1.2}}} };
        double  r[5] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
        s_points pts = { 5, c };
        check_mc("N=5 partial overlap",
                 volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff),
                 mc_volume_union_spheres(&pts, r, N_MC, SEED), MC_TOL);
    }

    dynarray_free(&buff);
    printf("\n%d passed, %d failed\n", n_passed, n_failed);
    return n_failed ? 1 : 0;
}
