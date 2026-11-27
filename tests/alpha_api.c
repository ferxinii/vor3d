/* -----------------------------------------------------------------------------
 * alpha_api.c -- tests for the public alpha API (include/alpha.h):
 * alpha_complex_3d / alpha_shape_3d.  See CLAUDE/PLAN_ALPHA_API.md step 5.
 *
 *   (1) unweighted: complex tets == shape in-tets before repair (volume, when
 *       repair changed nothing); empty K at tiny alpha is not a failure.
 *   (2) closure: every face of a listed tet and every edge of a listed face is
 *       listed; records ascending, lists sorted, indices in range.
 *   (3) out_kept_idx round trip (coordinates match, each output point hit
 *       once) for both calls; weights copied.
 *   (4) weighted, alpha = 0:
 *       - every ball whose centre lies in its own power cell (brute force) is
 *         a vertex of K; every vertex of K is a non-redundant ball;
 *       - per-ball union contributions are ~0 for balls not in K;
 *       - shape volume <= volume_union_spheres.
 *   (5) NAN auto: dropped == duplicates + redundant points (weighted) / 0
 *       (unweighted), i.e. every DT vertex touches an in-tet.
 *   (6) degenerate inputs (N < 4, coplanar, all duplicates, NaN alpha for the
 *       complex) fail with alpha = NAN and out_kept_idx all -1.
 *
 * Build (from repo root, after `cmake --build build`):
 *   CLAUDE/alpha_api_baseline/cc_test.sh tests/alpha_api.c /tmp/alpha_api && /tmp/alpha_api
 * --------------------------------------------------------------------------- */

#include "alpha.h"
#include "volsph.h"
#include "scplx.h"
#include "dynarray.h"
#include "random.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int g_fail = 0;
#define CHECK(c, ...) do { if (!(c)) { g_fail = 1; printf("  FAIL: " __VA_ARGS__); printf("\n"); } } while (0)

static unsigned long g_seed = 0x1234567UL;
static double fr(void)
{
    g_seed = g_seed * 6364136223846793005UL + 1442695040888963407UL;
    return (double)(g_seed >> 11) / (double)(1UL << 53);
}

static int cmp_lex(const int *a, const int *b, int k)
{
    for (int i = 0; i < k; i++) if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    return 0;
}

/* Binary search in a lexicographically sorted list of k-records. */
static int has_rec(const int *arr, int n, int k, const int *key)
{
    int lo = 0, hi = n - 1;
    while (lo <= hi) {
        int mid = (lo + hi) / 2, c = cmp_lex(&arr[mid*k], key, k);
        if (c == 0) return 1;
        if (c < 0) lo = mid + 1; else hi = mid - 1;
    }
    return 0;
}

static void check_list(const int *arr, int n, int k, int Np, const char *what)
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) CHECK(arr[i*k+j] >= 0 && arr[i*k+j] < Np, "%s index range", what);
        for (int j = 0; j + 1 < k; j++) CHECK(arr[i*k+j] < arr[i*k+j+1], "%s record not ascending", what);
        if (i) CHECK(cmp_lex(&arr[(i-1)*k], &arr[i*k], k) < 0, "%s list not sorted/unique", what);
    }
}

static void check_complex(const s_alpha_complex *ac, const s_points *pts, const double *w,
                          const int *kept, const char *tag)
{
    check_list(ac->tets, ac->N_tets, 4, ac->points.N, "tets");
    check_list(ac->faces, ac->N_faces, 3, ac->points.N, "faces");
    check_list(ac->edges, ac->N_edges, 2, ac->points.N, "edges");
    for (int t = 0; t < ac->N_tets; t++)
        for (int o = 0; o < 4; o++) {
            int f[3], m = 0;
            for (int j = 0; j < 4; j++) if (j != o) f[m++] = ac->tets[4*t+j];
            CHECK(has_rec(ac->faces, ac->N_faces, 3, f), "%s: face of a tet missing", tag);
        }
    for (int t = 0; t < ac->N_faces; t++) {
        const int *F = &ac->faces[3*t];
        int e[3][2] = {{F[0],F[1]},{F[0],F[2]},{F[1],F[2]}};
        for (int j = 0; j < 3; j++) CHECK(has_rec(ac->edges, ac->N_edges, 2, e[j]), "%s: edge of a face missing", tag);
    }
    CHECK((ac->weights != NULL) == (w != NULL), "%s: weights presence", tag);
    int *hit = calloc(ac->points.N + 1, sizeof(int));
    for (int i = 0; i < pts->N; i++) {
        if (kept[i] < 0) continue;
        CHECK(kept[i] < ac->points.N, "%s: kept range", tag);
        s_point a = pts->p[i], b = ac->points.p[kept[i]];
        CHECK(a.x == b.x && a.y == b.y && a.z == b.z, "%s: kept coords", tag);
        if (w) CHECK(ac->weights[kept[i]] == w[i], "%s: weight copy", tag);
        hit[kept[i]]++;
    }
    for (int k = 0; k < ac->points.N; k++) CHECK(hit[k] == 1, "%s: K point %d hit %d times", tag, k, hit[k]);
    free(hit);
}

static double complex_tet_volume(const s_alpha_complex *ac)
{
    double v = 0;
    for (int t = 0; t < ac->N_tets; t++) {
        s_point T[4];
        for (int j = 0; j < 4; j++) T[j] = ac->points.p[ac->tets[4*t+j]];
        v += fabs(signed_volume_tetra(T));
    }
    return v;
}

static void check_shape_kept(const s_alpha_shape *as, const s_points *pts, const int *kept, const char *tag)
{
    for (int i = 0; i < pts->N; i++) {
        if (kept[i] < 0) continue;
        CHECK(kept[i] < as->surface.points.N, "%s: shape kept range", tag);
        s_point a = pts->p[i], b = as->surface.points.p[kept[i]];
        CHECK(a.x == b.x && a.y == b.y && a.z == b.z, "%s: shape kept coords", tag);
    }
}

int main(void)
{
    const int N = 1500;
    s_point *p = malloc(N * sizeof *p);
    double *w = malloc(N * sizeof *w);
    int *kept = malloc(N * sizeof(int));
    for (int i = 0; i < N; i++) { p[i] = (s_point){{{ fr(), fr(), fr() }}}; w[i] = 0.004 * fr(); }
    s_points pts = { N, p };

    /* (1)-(3) unweighted */
    printf("[unweighted] complex vs shape, closure, kept\n");
    double alphas[] = { 1e-6, 0.002, 0.005, 0.01, 0.05, 1e3 };
    for (int a = 0; a < 6; a++) {
        s_random_context r1 = random_initialize(3), r2 = random_initialize(3);
        s_alpha_complex ac = alpha_complex_3d(&pts, NULL, alphas[a], 1e-12, &r1, kept);
        CHECK(!isnan(ac.alpha) && ac.alpha == alphas[a], "complex alpha=%g failed", alphas[a]);
        check_complex(&ac, &pts, NULL, kept, "unweighted");
        s_alpha_shape as = alpha_shape_3d(&pts, NULL, alphas[a], 1e-12, &r2, kept);
        double vk = complex_tet_volume(&ac);
        if (ac.N_tets == 0) CHECK(isnan(as.alpha), "empty in-set must fail the shape");
        else {
            CHECK(!isnan(as.alpha), "shape alpha=%g failed", alphas[a]);
            check_shape_kept(&as, &pts, kept, "unweighted");
            if (as.N_promoted == 0)
                CHECK(fabs(vk - as.volume) <= 1e-12 * (1 + vk), "alpha=%g complex vol %.17g vs shape %.17g",
                      alphas[a], vk, as.volume);
        }
        printf("  alpha=%-6g V=%d E=%d F=%d T=%d vol=%.6g | shape vol=%.6g promoted=%d\n", alphas[a],
               ac.points.N, ac.N_edges, ac.N_faces, ac.N_tets, vk, as.volume, as.N_promoted);
        free_alpha_shape(&as);
        free_alpha_complex(&ac);
    }

    /* (4) weighted */
    printf("[weighted] closure, alpha=0 vertex set, union volume\n");
    double walphas[] = { -0.001, 0.0, 0.002, 1e3 };
    int N_nonredundant = -1;
    for (int a = 0; a < 4; a++) {
        s_alpha_complex ac = alpha_complex_3d(&pts, w, walphas[a], 1e-12, NULL, kept);
        CHECK(!isnan(ac.alpha), "weighted complex alpha=%g failed", walphas[a]);
        check_complex(&ac, &pts, w, kept, "weighted");
        printf("  alpha=%-6g V=%d E=%d F=%d T=%d\n", walphas[a], ac.points.N, ac.N_edges, ac.N_faces, ac.N_tets);
        if (walphas[a] == 1e3) N_nonredundant = ac.points.N;   /* every DT vertex is in K */
        free_alpha_complex(&ac);
    }
    {
        int *kept_big = malloc(N * sizeof(int));
        s_alpha_complex big = alpha_complex_3d(&pts, w, 1e3, 1e-12, NULL, kept_big);
        s_alpha_complex ac0 = alpha_complex_3d(&pts, w, 0.0, 1e-12, NULL, kept);
        double *rad = malloc(N * sizeof(double)), *contrib = malloc(N * sizeof(double));
        for (int i = 0; i < N; i++) rad[i] = sqrt(w[i]);
        s_dynarray buff = dynarray_initialize(sizeof(s_ncell *), 64);
        double Vu = volume_union_spheres(&pts, rad, 1e-10, 1e-12, &buff, contrib);
        int n_own = 0, n_bad_own = 0, n_bad_red = 0, n_bad_contrib = 0;
        for (int i = 0; i < N; i++) {
            /* centre in its own power cell: pi_j(p_i) >= pi_i(p_i) for all j */
            bool own = true;
            for (int j = 0; j < N && own; j++)
                if (j != i && distance_squared(p[i], p[j]) - w[j] < -w[i]) own = false;
            n_own += own;
            if (own && kept[i] < 0) n_bad_own++;
            if (kept[i] >= 0 && kept_big[i] < 0) n_bad_red++;
            if (kept[i] < 0 && fabs(contrib[i]) > 1e-12 * Vu) n_bad_contrib++;
        }
        CHECK(n_bad_own == 0, "%d balls with centre in own power cell missing from K(0)", n_bad_own);
        CHECK(n_bad_red == 0, "%d redundant balls in K(0)", n_bad_red);
        CHECK(n_bad_contrib == 0, "%d balls outside K(0) with nonzero union contribution", n_bad_contrib);
        s_alpha_shape as = alpha_shape_3d(&pts, w, 0.0, 1e-12, NULL, kept);
        CHECK(!isnan(as.alpha), "weighted shape alpha=0 failed");
        CHECK(as.volume <= Vu, "shape vol %g > union %g", as.volume, Vu);
        check_shape_kept(&as, &pts, kept, "weighted");
        printf("  K(0): V=%d (centre-in-own-cell %d, non-redundant %d); shape vol=%.6g <= union %.6g\n",
               ac0.points.N, n_own, big.points.N, as.volume, Vu);
        free_alpha_shape(&as);
        dynarray_free(&buff);
        free(rad); free(contrib); free(kept_big);
        free_alpha_complex(&ac0);
        free_alpha_complex(&big);
    }

    /* (5) NAN auto */
    printf("[auto] dropped points\n");
    {
        s_alpha_shape as = alpha_shape_3d(&pts, w, NAN, 1e-12, NULL, kept);
        CHECK(!isnan(as.alpha), "weighted auto failed");
        CHECK(as.N_pts_dropped == N - N_nonredundant, "weighted auto dropped %d, redundant %d",
              as.N_pts_dropped, N - N_nonredundant);
        check_shape_kept(&as, &pts, kept, "weighted auto");
        printf("  weighted:   alpha=%g dropped=%d (redundant %d)\n", as.alpha, as.N_pts_dropped, N - N_nonredundant);
        free_alpha_shape(&as);

        /* unweighted, with 5 exact duplicates appended: dropped == 5 */
        s_point *pd = malloc((N + 5) * sizeof *pd);
        memcpy(pd, p, N * sizeof *pd);
        for (int i = 0; i < 5; i++) pd[N+i] = p[7*i];
        int *kd = malloc((N + 5) * sizeof(int));
        as = alpha_shape_3d(&(s_points){ N + 5, pd }, NULL, NAN, 1e-12, NULL, kd);
        CHECK(!isnan(as.alpha) && as.N_pts_dropped == 5, "unweighted auto + 5 dups: dropped %d", as.N_pts_dropped);
        printf("  unweighted: alpha=%g dropped=%d (5 duplicates)\n", as.alpha, as.N_pts_dropped);
        free_alpha_shape(&as);
        free(pd); free(kd);
    }

    /* (6) failures */
    printf("[degenerate] failures\n");
    {
        s_point cop[20], dup[10];
        for (int i = 0; i < 20; i++) cop[i] = (s_point){{{ fr(), fr(), 0.0 }}};
        for (int i = 0; i < 10; i++) dup[i] = p[0];
        s_points cases[3] = { { 3, p }, { 20, cop }, { 10, dup } };
        const char *names[3] = { "N=3", "coplanar", "all duplicates" };
        for (int c = 0; c < 3; c++) {
            for (int i = 0; i < cases[c].N; i++) kept[i] = 123;
            s_alpha_complex f = alpha_complex_3d(&cases[c], NULL, 1.0, 1e-12, NULL, kept);
            CHECK(isnan(f.alpha) && f.points.N == 0 && !f.tets, "complex %s should fail", names[c]);
            for (int i = 0; i < cases[c].N; i++) CHECK(kept[i] == -1, "complex %s kept not -1", names[c]);
            free_alpha_complex(&f);
            for (int i = 0; i < cases[c].N; i++) kept[i] = 123;
            s_alpha_shape fs = alpha_shape_3d(&cases[c], NULL, NAN, 1e-12, NULL, kept);
            CHECK(isnan(fs.alpha) && fs.surface.points.N == 0, "shape %s should fail", names[c]);
            for (int i = 0; i < cases[c].N; i++) CHECK(kept[i] == -1, "shape %s kept not -1", names[c]);
            free_alpha_shape(&fs);
        }
        s_alpha_complex f = alpha_complex_3d(&pts, NULL, NAN, 1e-12, NULL, NULL);
        CHECK(isnan(f.alpha), "complex NaN alpha should fail");
        free_alpha_complex(&f);
    }

    printf("%s\n", g_fail ? "SOME TESTS FAILED" : "ALL PASS");
    free(p); free(w); free(kept);
    return g_fail;
}
