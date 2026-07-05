/*
 * Test harness for alpha_shape_3d (src/ashape.c).
 *
 * Builds several point-cloud fixtures, extracts the alpha-shape surface, and
 * checks: mesh validity, positive volume, the volume invariant (enclosed
 * volume == sum of in-tet volumes), the Euler characteristic / genus of closed
 * fixtures, connected-component counts, and auto-alpha coverage.  Writes each
 * surface to a .obj in the working directory.
 */

#include "ashape.h"
#include "trimesh.h"
#include "points.h"
#include "random.h"
#include "dynarray.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

static int g_fail = 0;

static void check(const char *name, int ok)
{
    printf("    [%s] %s\n", ok ? "PASS" : "FAIL", name);
    if (!ok) g_fail = 1;
}

/* Legacy-VTK PolyData writer for the input cloud (POINTS + VERTICES).
 * Open in ParaView; use the "Points"/"Point Gaussian" representation. */
static void write_vtk_points(const char *path, const s_points *pts)
{
    FILE *f = fopen(path, "w");
    if (!f) return;
    fprintf(f, "# vtk DataFile Version 3.0\ninput points\nASCII\n");
    fprintf(f, "DATASET POLYDATA\nPOINTS %d double\n", pts->N);
    for (int i = 0; i < pts->N; i++)
        fprintf(f, "%.17g %.17g %.17g\n", pts->p[i].x, pts->p[i].y, pts->p[i].z);
    fprintf(f, "VERTICES %d %d\n", pts->N, 2 * pts->N);
    for (int i = 0; i < pts->N; i++) fprintf(f, "1 %d\n", i);
    fclose(f);
}

/* Legacy-VTK PolyData writer for the alpha-shape surface (POINTS + POLYGONS). */
static void write_vtk_surface(const char *path, const s_trimesh *m)
{
    FILE *f = fopen(path, "w");
    if (!f) return;
    fprintf(f, "# vtk DataFile Version 3.0\nalpha shape surface\nASCII\n");
    fprintf(f, "DATASET POLYDATA\nPOINTS %d double\n", m->points.N);
    for (int i = 0; i < m->points.N; i++)
        fprintf(f, "%.17g %.17g %.17g\n",
                m->points.p[i].x, m->points.p[i].y, m->points.p[i].z);
    fprintf(f, "POLYGONS %d %d\n", m->Nf, 4 * m->Nf);
    for (int i = 0; i < m->Nf; i++)
        fprintf(f, "3 %d %d %d\n", m->faces[i*3], m->faces[i*3+1], m->faces[i*3+2]);
    fclose(f);
}

/* Write both the input cloud and the resulting surface for a case. */
static void dump_case(const char *name, const s_points *pts, const s_trimesh *m)
{
    char path[256];
    snprintf(path, sizeof path, "%s_points.vtk", name);
    write_vtk_points(path, pts);
    if (trimesh_is_valid(m)) {
        snprintf(path, sizeof path, "%s_surface.vtk", name);
        write_vtk_surface(path, m);
    }
}

/* Euler characteristic of a closed triangle mesh: V - E + F, E = 3F/2. */
static int euler_char(const s_trimesh *m)
{
    return m->points.N - (m->Nf * 3) / 2 + m->Nf;
}

static void print_info(const s_ashape_info *info)
{
    printf("        alpha=%.5g  tets_in=%d  promoted=%d  comps=%d  dropped=%d  vol=%.6g\n",
           info->alpha_used, info->N_tets_in, info->N_promoted,
           info->N_components, info->N_pts_dropped, info->volume);
}

/* Common checks for every mesh expected to be a valid closed solid. */
static void basic_checks(const s_trimesh *m, const s_ashape_info *info)
{
    print_info(info);
    check("mesh is valid", trimesh_is_valid(m));
    if (!trimesh_is_valid(m)) return;
    double vt = volume_trimesh(m);
    check("positive volume", vt > 0.0);
    double rel = fabs(vt - info->volume) / (info->volume > 0 ? info->volume : 1.0);
    check("volume invariant (mesh == tets)", rel < 1e-9);
}


/* ---- Fixtures ---------------------------------------------------------- */

static s_points ball_cloud(s_random_context *rng, int N, double radius, s_point center)
{
    s_points pts = { .N = N, .p = malloc(sizeof(s_point) * N) };
    for (int i = 0; i < N; ) {
        double x = 2*random_uniform_double(rng) - 1;
        double y = 2*random_uniform_double(rng) - 1;
        double z = 2*random_uniform_double(rng) - 1;
        if (x*x + y*y + z*z > 1.0) continue;
        pts.p[i].x = center.x + radius*x;
        pts.p[i].y = center.y + radius*y;
        pts.p[i].z = center.z + radius*z;
        i++;
    }
    return pts;
}

static s_points torus_cloud(s_random_context *rng, int N, double R, double r)
{
    s_points pts = { .N = N, .p = malloc(sizeof(s_point) * N) };
    for (int i = 0; i < N; ) {
        double u = 2*random_uniform_double(rng) - 1;
        double v = 2*random_uniform_double(rng) - 1;
        if (u*u + v*v > 1.0) continue;
        double theta = 2*M_PI*random_uniform_double(rng);
        double rr = R + r*u;
        pts.p[i].x = rr*cos(theta);
        pts.p[i].y = rr*sin(theta);
        pts.p[i].z = r*v;
        i++;
    }
    return pts;
}

/* Jittered grid filling the L-shape [0,1]^3 minus the octant x>0.5 && y>0.5. */
static s_points lshape_cloud(s_random_context *rng, double h)
{
    int n = (int)(1.0 / h) + 1;
    s_dynarray buf = dynarray_initialize(sizeof(s_point), 1024);
    for (int i = 0; i <= n; i++) for (int j = 0; j <= n; j++) for (int k = 0; k <= n; k++) {
        double x = i * h, y = j * h, z = k * h;
        if (x > 1.0 || y > 1.0 || z > 1.0) continue;
        if (x > 0.5 && y > 0.5) continue;            /* carve the notch */
        double jit = 0.15 * h;
        x += jit*(2*random_uniform_double(rng)-1);
        y += jit*(2*random_uniform_double(rng)-1);
        z += jit*(2*random_uniform_double(rng)-1);
        s_point p = {{{x, y, z}}};
        dynarray_push(&buf, &p);
    }
    s_points pts = { .N = (int)buf.N, .p = malloc(sizeof(s_point) * buf.N) };
    for (unsigned i = 0; i < buf.N; i++) pts.p[i] = *(s_point*)dynarray_get_ptr(&buf, i);
    dynarray_free(&buf);
    return pts;
}

/* Two balls joined by a thin bridge of points -- stresses manifold repair. */
static s_points dumbbell_cloud(s_random_context *rng, int seed_scale)
{
    s_points a = ball_cloud(rng, 600 + seed_scale, 1.0, (s_point){{{0,0,0}}});
    s_points b = ball_cloud(rng, 600 + seed_scale, 1.0, (s_point){{{4,0,0}}});
    int Nbridge = 40;
    int N = a.N + b.N + Nbridge;
    s_points pts = { .N = N, .p = malloc(sizeof(s_point) * N) };
    int idx = 0;
    for (int i = 0; i < a.N; i++) pts.p[idx++] = a.p[i];
    for (int i = 0; i < b.N; i++) pts.p[idx++] = b.p[i];
    for (int i = 0; i < Nbridge; i++) {
        double t = (double)i / (Nbridge - 1);
        pts.p[idx].x = t * 4.0;
        pts.p[idx].y = 0.2 * (2*random_uniform_double(rng)-1);
        pts.p[idx].z = 0.2 * (2*random_uniform_double(rng)-1);
        idx++;
    }
    free_points(&a); free_points(&b);
    return pts;
}


/* Time just the extraction call (excludes cloud generation and file I/O). */
static s_trimesh timed_ashape(const s_points *pts, double alpha, double tol,
                              s_ashape_info *info)
{
    clock_t t0 = clock();
    s_trimesh m = alpha_shape_3d(pts, alpha, tol, info);
    double ms = 1000.0 * (double)(clock() - t0) / CLOCKS_PER_SEC;
    printf("        [extract] %.1f ms  (%d input pts)\n", ms, pts->N);
    return m;
}


/* ---- Cases ------------------------------------------------------------- */

static void case_sphere(s_random_context *rng)
{
    printf("  Case: sphere cloud (2000 pts, r=1)\n");
    s_points pts = ball_cloud(rng, 2000, 1.0, (s_point){{{0,0,0}}});
    s_ashape_info info;
    s_trimesh m = timed_ashape(&pts, 0.04, 1e-9, &info);  /* alpha = 0.2^2 */
    basic_checks(&m, &info);
    if (trimesh_is_valid(&m)) {
        check("sphere genus 0 (Euler 2)", euler_char(&m) == 2);
        check("one component", info.N_components == 1);
        check("volume below hull", volume_trimesh(&m) < 4.19);   /* 4/3 pi */
    }
    dump_case("ashape_sphere", &pts, &m);
    free_trimesh(&m);
    free_points(&pts);
}

static void case_torus(s_random_context *rng)
{
    printf("  Case: torus cloud (6000 pts, R=1, r=0.3)\n");
    s_points pts = torus_cloud(rng, 6000, 1.0, 0.3);
    s_ashape_info info;
    s_trimesh m = timed_ashape(&pts, 0.02, 1e-9, &info);  /* alpha ~ 0.14^2 */
    basic_checks(&m, &info);
    if (trimesh_is_valid(&m)) {
        check("torus genus 1 (Euler 0)", euler_char(&m) == 0);
    }
    dump_case("ashape_torus", &pts, &m);
    free_trimesh(&m);
    free_points(&pts);
}

/* Auto-alpha picks the SMALLEST alpha covering every point: the tightest
 * wrap.  It guarantees coverage, validity, and connectivity, but the surface
 * is deliberately jagged (under-fills volume, may carry small handles), so we
 * do NOT assert a particular genus here. */
static void case_torus_auto(s_random_context *rng)
{
    printf("  Case: torus cloud, AUTO alpha\n");
    s_points pts = torus_cloud(rng, 6000, 1.0, 0.3);
    s_ashape_info info;
    s_trimesh m = timed_ashape(&pts, -1.0, 1e-9, &info);
    basic_checks(&m, &info);
    if (trimesh_is_valid(&m)) {
        check("auto covers all points (dropped 0)", info.N_pts_dropped == 0);
        check("single component", info.N_components == 1);
    }
    dump_case("ashape_torus_auto", &pts, &m);
    free_trimesh(&m);
    free_points(&pts);
}

static void case_sphere_auto(s_random_context *rng)
{
    printf("  Case: sphere cloud, AUTO alpha\n");
    s_points pts = ball_cloud(rng, 2000, 1.0, (s_point){{{0,0,0}}});
    s_ashape_info info;
    s_trimesh m = timed_ashape(&pts, -1.0, 1e-9, &info);
    basic_checks(&m, &info);
    if (trimesh_is_valid(&m)) {
        check("auto covers all points (dropped 0)", info.N_pts_dropped == 0);
        check("single component", info.N_components == 1);
    }
    dump_case("ashape_sphere_auto", &pts, &m);
    free_trimesh(&m);
    free_points(&pts);
}

static void case_two_clusters(s_random_context *rng)
{
    printf("  Case: two separated clusters\n");
    s_points a = ball_cloud(rng, 1200, 1.0, (s_point){{{0,0,0}}});
    s_points b = ball_cloud(rng, 1200, 1.0, (s_point){{{10,0,0}}});
    s_points pts = { .N = a.N + b.N, .p = malloc(sizeof(s_point)*(a.N+b.N)) };
    for (int i = 0; i < a.N; i++) pts.p[i] = a.p[i];
    for (int i = 0; i < b.N; i++) pts.p[a.N+i] = b.p[i];
    s_ashape_info info;
    s_trimesh m = timed_ashape(&pts, 0.04, 1e-9, &info);
    basic_checks(&m, &info);
    if (trimesh_is_valid(&m)) {
        check("two components", info.N_components == 2);
    }
    dump_case("ashape_two_clusters", &pts, &m);
    free_trimesh(&m);
    free_points(&pts); free_points(&a); free_points(&b);
}

static void case_lshape(s_random_context *rng)
{
    printf("  Case: L-shape (jittered grid, notched cube)\n");
    double h = 0.04;
    s_points pts = lshape_cloud(rng, h);
    s_ashape_info info;
    /* Fixed alpha ~ (2.5h)^2 fills the interior solidly while leaving the
     * empty notch carved out. */
    s_trimesh m = timed_ashape(&pts, (2.5*h)*(2.5*h), 1e-9, &info);
    basic_checks(&m, &info);
    if (trimesh_is_valid(&m)) {
        double vt = volume_trimesh(&m);
        printf("        L-volume=%.4g (analytic 0.75, cube hull ~1)\n", vt);
        /* Crisp non-convexity test: the notch centre is outside the solid,
         * the bulk centre is inside. */
        int notch_out = point_in_trimesh(&m, (s_point){{{0.75,0.75,0.5}}}, 1e-9, 8);
        int bulk_in   = point_in_trimesh(&m, (s_point){{{0.25,0.25,0.5}}}, 1e-9, 8);
        check("notch centre is OUTSIDE (concavity carved)", notch_out == 0);
        check("bulk centre is INSIDE (solid)", bulk_in == 1);
        check("volume in solid range [0.55,0.85]", vt > 0.55 && vt < 0.85);
    }
    dump_case("ashape_lshape", &pts, &m);
    free_trimesh(&m);
    free_points(&pts);
}

static void case_repair_stress(s_random_context *rng)
{
    printf("  Case: dumbbell repair stress (12 seeds)\n");
    int all_ok = 1, comp_min = 1e9, comp_max = 0;
    for (int s = 0; s < 12; s++) {
        s_points pts = dumbbell_cloud(rng, s * 13);
        s_ashape_info info;
        s_trimesh m = timed_ashape(&pts, 0.09, 1e-9, &info);   /* alpha = 0.3^2 */
        int ok = trimesh_is_valid(&m);
        if (ok) {
            double vt = volume_trimesh(&m);
            double rel = fabs(vt - info.volume) / (info.volume > 0 ? info.volume : 1.0);
            if (!(vt > 0.0 && rel < 1e-9)) ok = 0;
        }
        if (info.N_components < comp_min) comp_min = info.N_components;
        if (info.N_components > comp_max) comp_max = info.N_components;
        printf("        seed %2d: comps=%d promoted=%d vol=%.4g%s\n",
               s, info.N_components, info.N_promoted, info.volume,
               ok ? "" : "  <-- FAILED");
        if (!ok) all_ok = 0;
        if (s == 0) dump_case("ashape_dumbbell", &pts, &m);
        free_trimesh(&m);
        free_points(&pts);
    }
    printf("        components across seeds: min=%d max=%d\n", comp_min, comp_max);
    check("all dumbbell seeds valid + volume invariant", all_ok);
}

static void case_degenerate(void)
{
    printf("  Case: degenerate inputs -> trimesh_NAN\n");
    s_points p3 = { .N = 3, .p = malloc(sizeof(s_point)*3) };
    p3.p[0]=(s_point){{{0,0,0}}}; p3.p[1]=(s_point){{{1,0,0}}}; p3.p[2]=(s_point){{{0,1,0}}};
    s_trimesh m3 = alpha_shape_3d(&p3, 1.0, 1e-9, NULL);
    check("3 points -> invalid", !trimesh_is_valid(&m3));
    free_trimesh(&m3); free_points(&p3);

    s_points pc = { .N = 25, .p = malloc(sizeof(s_point)*25) };
    for (int i = 0; i < 25; i++) pc.p[i] = (s_point){{{ (double)(i%5), (double)(i/5), 0.0 }}};
    s_trimesh mc = alpha_shape_3d(&pc, 100.0, 1e-9, NULL);
    check("coplanar -> invalid", !trimesh_is_valid(&mc));
    free_trimesh(&mc); free_points(&pc);
}


int main(void)
{
    s_random_context rng = random_initialize(12345);

    printf("=== alpha_shape_3d tests ===\n");
    case_sphere(&rng);
    case_sphere_auto(&rng);
    case_torus(&rng);
    case_torus_auto(&rng);
    case_two_clusters(&rng);
    case_lshape(&rng);
    case_repair_stress(&rng);
    case_degenerate();

    printf("\n%s\n", g_fail ? "SOME TESTS FAILED" : "ALL TESTS PASSED");
    return g_fail;
}
