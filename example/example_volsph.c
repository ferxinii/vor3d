#include "volsph.h"
#include "dynarray.h"
#include "points.h"
#include "scplx.h"
#include "delaunay.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#define PI M_PI
#define TOL 1e-6
#define EPS 1e-10
#define TOL_DUP 1e-12

static int n_passed = 0;
static int n_failed = 0;

static void check(const char *name, double got, double expected, double tol)
{
    double err = fabs(got - expected);
    double rel = fabs(expected) > 1e-14 ? err / fabs(expected) : err;
    if (rel < tol) {
        printf("  PASS  %s  (got %.10f, expected %.10f)\n", name, got, expected);
        n_passed++;
    } else {
        printf("  FAIL  %s  (got %.10f, expected %.10f, rel_err %.2e)\n",
               name, got, expected, rel);
        n_failed++;
    }
}




/* Simple MC estimate of union volume, for validation only */
static double mc_volume_union_spheres(const s_points *centers,
                                      const double *radii,
                                      int N_samples, unsigned int seed)
{
    /* Compute bounding box */
    double xmin = centers->p[0].x - radii[0];
    double xmax = centers->p[0].x + radii[0];
    double ymin = centers->p[0].y - radii[0];
    double ymax = centers->p[0].y + radii[0];
    double zmin = centers->p[0].z - radii[0];
    double zmax = centers->p[0].z + radii[0];
    for (int i = 1; i < centers->N; i++) {
        xmin = fmin(xmin, centers->p[i].x - radii[i]);
        xmax = fmax(xmax, centers->p[i].x + radii[i]);
        ymin = fmin(ymin, centers->p[i].y - radii[i]);
        ymax = fmax(ymax, centers->p[i].y + radii[i]);
        zmin = fmin(zmin, centers->p[i].z - radii[i]);
        zmax = fmax(zmax, centers->p[i].z + radii[i]);
    }
    double Lx = xmax - xmin;
    double Ly = ymax - ymin;
    double Lz = zmax - zmin;
    double box_vol = Lx * Ly * Lz;

    srand(seed);
    int hits = 0;
    for (int s = 0; s < N_samples; s++) {
        double x = xmin + Lx * ((double)rand() / RAND_MAX);
        double y = ymin + Ly * ((double)rand() / RAND_MAX);
        double z = zmin + Lz * ((double)rand() / RAND_MAX);
        for (int i = 0; i < centers->N; i++) {
            double dx = x - centers->p[i].x;
            double dy = y - centers->p[i].y;
            double dz = z - centers->p[i].z;
            if (dx*dx + dy*dy + dz*dz <= radii[i]*radii[i]) {
                hits++;
                break;  /* union: count once */
            }
        }
    }
    return box_vol * ((double)hits / N_samples);
}

static void check_mc(const char *name, double exact, double mc,
                     double mc_tol, int *n_passed, int *n_failed)
{
    double err = fabs(exact - mc) / fmax(fabs(exact), 1e-10);
    if (err < mc_tol) {
        printf("  PASS  %s  exact=%.6f  mc=%.6f  rel_err=%.2e\n",
               name, exact, mc, err);
        (*n_passed)++;
    } else {
        printf("  FAIL  %s  exact=%.6f  mc=%.6f  rel_err=%.2e\n",
               name, exact, mc, err);
        (*n_failed)++;
    }
}



void TEST_DIHEDRAL(const s_point p[4], double EPS_DEGEN);
void TEST_CAYLEY_MENGER(const s_point p[4]);

int main(void)
{
    srand(time(NULL));
    s_dynarray buff = dynarray_initialize(sizeof(s_ncell*), 64);
    (void)buff;
    
    // printf("=== DIHEDRAL ANGLES ===\n");
    // {
    //     s_point c[4] = { {{{0, 0, 0}}},
    //                      {{{2, 0, 0}}},
    //                      {{{2, 2, 0}}},
    //                      {{{1, 1, 1}}} };
    //     TEST_DIHEDRAL(c, EPS);
    //
    //     puts("");
    //     printf("%f\n", M_PI/4);
    //     printf("%f\n", M_PI/2);
    //     printf("%f\n", M_PI/3);
    //     printf("%f\n", M_PI/4);
    //     printf("%f\n", 2*M_PI/3);
    //     printf("%f\n", M_PI/3);
    // }
    //
    // printf("=== DIHEDRAL ANGLES ===\n");
    // {
    //     s_point c[4] = { {{{ 1, 1, 1}}},
    //                      {{{ 1,-1,-1}}},
    //                      {{{-1, 1,-1}}},
    //                      {{{-1,-1, 1}}} };
    //     TEST_DIHEDRAL(c, EPS);
    // }
    //
    // printf("=== DIHEDRAL ANGLES ===\n");
    // {
    //     s_point c[4] = { {{{ 1, 1, 0.5}}},
    //                      {{{ 1,-1,-1}}},
    //                      {{{-1, 1,-1}}},
    //                      {{{-1,-1, 2}}} };
    //     TEST_DIHEDRAL(c, EPS);
    // }
    //
    //
    //
    // printf("=== DT sanity check ===\n");
    // {
    //     s_point c[4] = { {{{ 1, 1, 1}}},
    //                      {{{ 1,-1,-1}}},
    //                      {{{-1, 1,-1}}},
    //                      {{{-1,-1, 1}}} };
    //     double  r[4]  = {0.8, 0.8, 0.8, 0.8};
    //     double  w[4];
    //     for (int i=0; i<4; i++) w[i] = r[i]*r[i];
    //     s_points pts = {4, c};
    //     s_scplx dt = construct_dt_3d(&pts, w, false, TOL_DUP);
    //     printf("N_ncells = %d (expected 1)\n", dt.N_ncells);
    //     // int ii = 0;
    //     // for (s_ncell *nc = dt.head; nc; nc = nc->next)
    //     //     printf("ncell[%d]: vids = %d %d %d %d\n", ii++,
    //     //            nc->vertex_id[0], nc->vertex_id[1],
    //     //            nc->vertex_id[2], nc->vertex_id[3]);
    //     free_complex(&dt);
    // }
    //
    // printf("=== Test 2: two non-overlapping spheres ===\n");
    // {
    //     /* Centres 10 apart, radii 1 — no overlap */
    //     s_point c[2] = { {{{0,0,0}}}, {{{10,0,0}}} };
    //     double  r[2] = {1.0, 1.0};
    //     s_points pts = {2, c};
    //     double got      = volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff);
    //     double expected = 2.0 * volume_sphere(1.0);
    //     check("N=2 non-overlapping", got, expected, TOL);
    // }
    //
    // printf("=== Test 3: two identical coincident spheres ===\n");
    // {
    //     /* Same centre — union is just one sphere. Use TOL_dup > 0 to avoid
    //      * degenerate DT, so we place them epsilon apart and use large radii. */
    //     s_point c[2] = {{ {{0,0,0}}}, {{{1e-4,0,0}}} };
    //     double  r[2] = {1.0, 1.0};
    //     s_points pts = {2, c};
    //     double got      = volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff);
    //     double expected = volume_union_2_spheres(c[0], r[0], c[1], r[1], true);
    //     check("N=2 nearly coincident vs volume_union_2_spheres", got, expected, TOL);
    // }
    //
    // printf("=== Test 4: three mutually intersecting spheres ===\n");
    // {
    //     /* Three equal spheres whose centres form an equilateral triangle in z = 0.
    //      * With r = 1.5 and side length d = 2.0:
    //      *   - pairwise overlaps exist (d < 2r)
    //      *   - triple intersection is non-empty (r > d / sqrt(3))
    //      *
    //      * Union = 3*V1 - 3*V2 + V3.
    //      */
    //     s_point c[3] = {
    //         {{{ 0.0, 0.0, 0.0}}},
    //         {{{ 2.0, 0.0, 0.0}}},
    //         {{{ 1.0, 1.7320508075688772, 0.0}}}  /* sqrt(3) */
    //     };
    //     double r[3] = {1.5, 1.5, 1.5};
    //     s_points pts = {3, c};
    //
    //     double got = volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff);
    //
    //     double d01 = distance(c[0], c[1]);
    //     double d02 = distance(c[0], c[2]);
    //     double d12 = distance(c[1], c[2]);
    //
    //     double V1 = volume_sphere(r[0]);
    //     double V2 = volume_intersection_2_spheres(r[0], r[1], d01, true);
    //     double V3 = volume_intersection_3_spheres(
    //         r[0], r[1], r[2],
    //         d01, d02, d12,
    //         true, EPS
    //     );
    //
    //     double expected = 3.0 * V1 - 3.0 * V2 + V3;
    //     check("N=3 mutually intersecting spheres", got, expected, TOL);
    // }
    //
    // printf("=== Test 6: four non-overlapping spheres ===\n");
    // {
    //     s_point c[4] = { {{{0,0,0}}}, {{{10,0,0}}}, {{{0,10,0}}}, {{{0,0,10}}} };
    //     double  r[4] = {1.0, 1.2, 0.8, 1.1};
    //     s_points pts = {4, c};
    //     double got      = volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff);
    //     double expected = volume_sphere(r[0]) + volume_sphere(r[1])
    //                       + volume_sphere(r[2]) + volume_sphere(r[3]);
    //     check("N=4 non-overlapping", got, expected, TOL);
    // }
    //
    // printf("=== Test 7: four spheres on regular tetrahedron vertices, all same radius ===\n");
    // {
    //     /* Regular tetrahedron with edge length L, all radii r.
    //      * Choose r small enough so only adjacent pairs overlap. 
    //      * Expected = inclusion-exclusion with 4 pairs, 4 triples, 1 quadruple. */
    //     double L = 2.0, r0 = 0.8;
    //     s_point c[4] = {
    //         {{{ 1,  1,  1}}},
    //         {{{ 1, -1, -1}}},
    //         {{{-1,  1, -1}}},
    //         {{{-1, -1,  1}}}
    //     };
    //     double  r[4] = {r0, r0, r0, r0};
    //     s_points pts = {4, c};
    //     double got = volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff);
    //
    //     /* Build expected via inclusion-exclusion manually */
    //     double dAB = distance(c[0], c[1]);
    //     double V1  = volume_sphere(r0);
    //     double Vii = volume_intersection_2_spheres(r0, r0, dAB, true);  /* same for all 6 edges */
    //     /* All triples have same pairwise distances */
    //     double Viii = volume_intersection_3_spheres(r0, r0, r0, dAB, dAB, dAB, true, EPS);
    //     double expected = 4*V1 - 6*Vii + 4*Viii;
    //     /* Note: 4-sphere intersection negligible if r0 < L/2 */
    //     check("N=4 regular tetrahedron, symmetric, expected has some error", got, expected, TOL);
    //     (void)L;
    // }
    //
    // printf("=== Test 8: one sphere contained inside another ===\n");
    // {
    //     s_point c[2] = { {{{0,0,0}}}, {{{0.1,0,0}}} };
    //     double  r[2] = {2.0, 0.5};  /* small sphere fully inside large */
    //     s_points pts = {2, c};
    //     double got      = volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff);
    //     double expected = volume_union_2_spheres(c[0], r[0], c[1], r[1], true);
    //     check("N=2 containment vs volume_union_2_spheres", got, expected, TOL);
    // }
    //
    // printf("=== Test 9: N=5 non-overlapping, sum of volumes ===\n");
    // {
    //     s_point c[5] = {
    //                     {{{0,0,0}}},
    //                     {{{10,0,0}}},
    //                     {{{0,10,0}}},
    //                     {{{0,0,10}}},
    //                     {{{10,10,0}}}
    //     };
    //     double  r[5] = {1.0, 0.9, 1.1, 0.8, 1.2};
    //     s_points pts = {5, c};
    //     double got      = volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff);
    //     double expected = 0;
    //     for (int i = 0; i < 5; i++) expected += volume_sphere(r[i]);
    //     check("N=5 non-overlapping", got, expected, TOL);
    // }
    //
    // printf("=== Test 10: all spheres mutually overlapping, large radii ===\n");
    // {
    //     /* All centres close together, large radii — union ≈ largest enclosing sphere.
    //      * Just check result is between max single sphere and sum of spheres. */
    //     s_point c[4] = { {{{0,0,0}}}, 
    //                      {{{0.1,0,0}}},
    //                      {{{0,0.1,0}}},
    //                      {{{0,0,0.1}}}
    //     };
    //     double  r[4] = {2.0, 2.0, 2.0, 2.0};
    //     s_points pts = {4, c};
    //     double got  = volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff);
    //     double vmin = volume_sphere(2.0);           /* at least one sphere */
    //     double vmax = 4.0 * volume_sphere(2.0);     /* at most sum */
    //     if (got >= vmin*(1-TOL) && got <= vmax*(1+TOL)) {
    //         printf("  PASS  N=4 large overlap (%.6f in [%.6f, %.6f])\n", got, vmin, vmax);
    //         n_passed++;
    //     } else {
    //         printf("  FAIL  N=4 large overlap (%.6f not in [%.6f, %.6f])\n", got, vmin, vmax);
    //         n_failed++;
    //     }
    // }
    //
    // dynarray_free(&buff);
    // printf("\n%d passed, %d failed\n", n_passed, n_failed);
    //

    (void)check;
    (void)mc_volume_union_spheres;
    (void)check_mc;
    printf("=== Monte Carlo validation tests ===\n");
    {
        /* MC tolerance: 1M samples gives ~0.1% accuracy for simple cases */
        const int    N_MC  = 100000000;
        const double MC_TOL = 0.01;  /* 1% relative tolerance */
        const unsigned int SEED = time(NULL);

        {
            s_point c[2] = {
                {{{ 0.0, 0.0, 0.0}}},
                {{{ 2.0, 0.0, 0.0}}},
            };
            double r[2] = {1.5, 1.1};
            s_points pts = {2, c};

            double exact = volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff);
            double mc    = mc_volume_union_spheres(&pts, r, N_MC, SEED);
            check_mc("2 overlapping", exact, mc, MC_TOL, &n_passed, &n_failed);
        }

        {
            s_point c[3] = {
                {{{ 0.0, 0.0, 0.0}}},
                {{{ 2.0, 0.0, 0.0}}},
                {{{ 1.0, 1.7320508075688772, 0.0}}}  /* sqrt(3) */
            };
            double r[3] = {1.5, 1.6, 1.9};
            s_points pts = {3, c};

            double exact = volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff);
            double mc    = mc_volume_union_spheres(&pts, r, N_MC, SEED);
            check_mc("3 overlapping", exact, mc, MC_TOL, &n_passed, &n_failed);
        }

        {
            s_point c[4] = {
                {{{ 1,  1,  1}}},
                {{{ 1, -1, -0.5}}},
                {{{-1,  1, -1}}},
                {{{-1, -1,  1}}}
            };

            double r[4] = {2.4, 1.8, 1.9, 1.5};
            s_points pts = {4, c};

            double exact = volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff);
            double mc    = mc_volume_union_spheres(&pts, r, N_MC, SEED);

            check_mc("4 overlapping (tetrahedral)", exact, mc, MC_TOL, &n_passed, &n_failed);
        }

        /* Test A: 5 non-overlapping spheres */
        {
            s_point c[5] = { {{{0,0,0}}}, {{{10,0,0}}}, {{{0,10,0}}},
                             {{{0,0,10}}}, {{{10,10,0}}} };
            double  r[5] = {1.0, 0.9, 1.1, 0.8, 1.2};
            s_points pts = {5, c};
            double exact = volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff);
            double mc    = mc_volume_union_spheres(&pts, r, N_MC, SEED);
            check_mc("5 non-overlapping", exact, mc, MC_TOL, &n_passed, &n_failed);
        }

        /* Test B: 5 partially overlapping spheres */
        {
            s_point c[5] = { {{{0,0,0}}}, {{{1.2,0,0}}}, {{{0.6,1.0,0}}},
                             {{{0,0,1.2}}}, {{{1.2,0,1.2}}} };
            double  r[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
            s_points pts = {5, c};
            double exact = volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff);
            double mc    = mc_volume_union_spheres(&pts, r, N_MC, SEED);
            check_mc("5 partial overlap", exact, mc, MC_TOL, &n_passed, &n_failed);
        }


        /* Test C: 6 spheres, mixed radii, some containment */
        {
            s_point c[6] = { {{{0,0,0}}}, {{{0.1,0,0}}}, {{{3,0,0}}},
                             {{{3,3,0}}}, {{{0,3,0}}}, {{{1.5,1.5,1.5}}} };
            double  r[6] = {2.0, 0.5, 1.0, 1.0, 1.0, 1.5};
            s_points pts = {6, c};
            double exact = volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff);
            double mc    = mc_volume_union_spheres(&pts, r, N_MC, SEED);
            check_mc("6 mixed radii+containment", exact, mc, MC_TOL, &n_passed, &n_failed);
        }


        /* Test D: 8 spheres on cube vertices, same radius */
        {
            s_point c[8] = { {{{0,0,0}}}, {{{2,0,0}}}, {{{0,2,0}}}, {{{2,2,0}}},
                             {{{0,0,2}}}, {{{2,0,2}}}, {{{0,2,2}}}, {{{2,2,2}}} };
            double r[8]; for (int i=0; i<8; i++) r[i] = 1.2;
            s_points pts = {8, c};
            double exact = volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff);
            double mc    = mc_volume_union_spheres(&pts, r, N_MC, SEED);
            check_mc("8 cube vertices r=1.2", exact, mc, MC_TOL, &n_passed, &n_failed);
        }

        /* Test E: 10 random spheres — stress test */
        {
            s_point c[10] = {
                {{{0.0, 0.0, 0.0}}}, {{{1.5, 0.3, 0.1}}}, {{{0.2, 1.6, 0.4}}},
                {{{1.4, 1.5, 0.2}}}, {{{0.7, 0.8, 1.5}}}, {{{2.1, 0.1, 1.3}}},
                {{{0.1, 2.2, 1.4}}}, {{{2.0, 2.1, 1.5}}}, {{{1.0, 1.0, 0.5}}},
                {{{1.5, 0.5, 2.0}}}
            };
            double r[10] = {0.8, 0.7, 0.9, 0.6, 0.8, 0.7, 0.6, 0.8, 0.5, 0.7};
            s_points pts = {10, c};
            double exact = volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff);
            double mc    = mc_volume_union_spheres(&pts, r, N_MC, SEED);
            check_mc("10 random spheres", exact, mc, MC_TOL, &n_passed, &n_failed);
        }

        {
            s_point c[4] = {
                {{{ 0,  0,  0}}},
                {{{ 1, 0, 0}}},
                {{{0,  1, 0}}},
                {{{0, 0,  1}}}
            };

            double r[4] = {10, 1, 1.1, 0.8};
            s_points pts = {4, c};

            double exact = volume_union_spheres(&pts, r, EPS, TOL_DUP, &buff);
            double mc    = mc_volume_union_spheres(&pts, r, N_MC, SEED);

            check_mc("3 contained inside big", exact, mc, MC_TOL, &n_passed, &n_failed);
        }
    

      }

    return 1;
}


