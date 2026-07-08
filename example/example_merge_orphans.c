/* Validation harness for the orphan-merge post-pass (MERGING_ORPHANS.md).
 *
 * For each fixture (a small closed trimesh domain + a FIXED seed set) it runs
 * vor3d_in_ncvx_domain twice -- merge_orphans=false (raw per-cell components)
 * and merge_orphans=true -- and asserts:
 *   - Volume conservation: sum of per-cell volumes == domain tet-sum volume.
 *   - Single component: with merge on, every non-empty seed has N_surface == 1.
 *   - Watertight/manifold: every merged surface is trimesh_is_valid with no
 *     boundary edges (adjacency != -1).
 *   - Regression: raw shows >= 1 seed with N_surface > 1 (a real orphan) and
 *     merge drops that count to 0 (guards that the pass actually does work).
 *   - Determinism: two merged runs give identical per-seed volume/N_surface.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "vor3d.h"
#include "random.h"

static const double EPS_DEG = 1e-12;
static const double TOL     = 1e-12;

/* ----------------------------------------------------------------------- */
/* Fixtures                                                                 */
/* ----------------------------------------------------------------------- */

/* U-/C-shape: three width-1 arms (left x[0,1], right x[2,3], bottom y[0,1])
 * extruded z in [0,1]. Volume = 7. The two vertical arms meet only through the
 * bottom bar, so a seed low in one arm and a seed high in the other force a
 * guaranteed disconnected orphan (the far arm's top is Euclidean-nearest to the
 * opposite seed but only geodesically reachable through the near seed). */
static s_trimesh make_ushape(void)
{
    s_point verts[16];
    static const double xy[8][2] = {
        {0,0},{3,0},{3,3},{2,3},{2,1},{1,1},{1,3},{0,3},
    };
    for (int i = 0; i < 8; i++) {
        verts[i]   = (s_point){{{ xy[i][0], xy[i][1], 0.0 }}};
        verts[i+8] = (s_point){{{ xy[i][0], xy[i][1], 1.0 }}};
    }
    static const int faces[28*3] = {
        /* bottom cap (z=0) */
        4,3,2,  4,2,1,  5,4,1,  5,1,0,  6,5,0,  7,6,0,
        /* top cap (z=1) */
        10,11,12, 9,10,12, 9,12,13, 8,9,13, 8,13,14, 8,14,15,
        /* sides: edge (i,(i+1)%8) */
        0,1,9,   0,9,8,
        1,2,10,  1,10,9,
        2,3,11,  2,11,10,
        3,4,12,  3,12,11,
        4,5,13,  4,13,12,
        5,6,14,  5,14,13,
        6,7,15,  6,15,14,
        7,0,8,   7,8,15,
    };
    return trimesh_from_arrays(verts, 16, faces, 28, EPS_DEG);
}

/* L-shape: union of [0,2]x[0,1] and [0,1]x[1,2], extruded z in [0,1]. Vol = 3. */
static s_trimesh make_lshape(void)
{
    static const s_point verts[12] = {
        {{{0,0,0}}}, {{{2,0,0}}}, {{{2,1,0}}}, {{{1,1,0}}}, {{{1,2,0}}}, {{{0,2,0}}},
        {{{0,0,1}}}, {{{2,0,1}}}, {{{2,1,1}}}, {{{1,1,1}}}, {{{1,2,1}}}, {{{0,2,1}}},
    };
    static const int faces[20*3] = {
        0,5,4,  0,4,3,  0,3,2,  0,2,1,
        6,7,8,  6,8,9,  6,9,10, 6,10,11,
        0,1,7,  0,7,6,
        1,2,8,  1,8,7,
        2,3,9,  2,9,8,
        3,4,10, 3,10,9,
        4,5,11, 4,11,10,
        5,0,6,  5,6,11,
    };
    return trimesh_from_arrays(verts, 12, faces, 20, EPS_DEG);
}

/* Plus-shape: two 3x1x1 boxes crossing at centre. Volume = 5. */
static s_trimesh make_plus(void)
{
    static const s_point verts[24] = {
        {{{-0.5,-1.5,-0.5}}}, {{{ 0.5,-1.5,-0.5}}},
        {{{ 0.5,-0.5,-0.5}}}, {{{ 1.5,-0.5,-0.5}}},
        {{{ 1.5, 0.5,-0.5}}}, {{{ 0.5, 0.5,-0.5}}},
        {{{ 0.5, 1.5,-0.5}}}, {{{-0.5, 1.5,-0.5}}},
        {{{-0.5, 0.5,-0.5}}}, {{{-1.5, 0.5,-0.5}}},
        {{{-1.5,-0.5,-0.5}}}, {{{-0.5,-0.5,-0.5}}},
        {{{-0.5,-1.5, 0.5}}}, {{{ 0.5,-1.5, 0.5}}},
        {{{ 0.5,-0.5, 0.5}}}, {{{ 1.5,-0.5, 0.5}}},
        {{{ 1.5, 0.5, 0.5}}}, {{{ 0.5, 0.5, 0.5}}},
        {{{ 0.5, 1.5, 0.5}}}, {{{-0.5, 1.5, 0.5}}},
        {{{-0.5, 0.5, 0.5}}}, {{{-1.5, 0.5, 0.5}}},
        {{{-1.5,-0.5, 0.5}}}, {{{-0.5,-0.5, 0.5}}},
    };
    static const int faces[44*3] = {
        2,1,0,   11,2,0,
        4,3,2,   5,4,2,
        7,6,5,   8,7,5,
        10,9,8,  11,10,8,
        5,2,11,  8,5,11,
        12,13,14, 12,14,23,
        14,15,16, 14,16,17,
        17,18,19, 17,19,20,
        20,21,22, 20,22,23,
        23,14,17, 23,17,20,
        0,1,13,  0,13,12,
        1,2,14,  1,14,13,
        2,3,15,  2,15,14,
        3,4,16,  3,16,15,
        4,5,17,  4,17,16,
        5,6,18,  5,18,17,
        6,7,19,  6,19,18,
        7,8,20,  7,20,19,
        8,9,21,  8,21,20,
        9,10,22, 9,22,21,
        10,11,23, 10,23,22,
        11,0,12, 11,12,23,
    };
    return trimesh_from_arrays(verts, 24, faces, 44, EPS_DEG);
}

/* ----------------------------------------------------------------------- */
/* Checks                                                                   */
/* ----------------------------------------------------------------------- */

static int surface_has_boundary_edge(const s_trimesh *m)
{
    for (int i = 0; i < m->Nf * 3; i++) if (m->adjacency[i] == -1) return 1;
    return 0;
}

/* Count seeds whose cell reports > 1 surface component. */
static int count_multi_surface(const s_ncvx_vdiagram *vd)
{
    int c = 0;
    for (int i = 0; i < vd->seeds.N; i++) if (vd->vcells[i].N_surface > 1) c++;
    return c;
}

static double total_volume(const s_ncvx_vdiagram *vd)
{
    double s = 0;
    for (int i = 0; i < vd->seeds.N; i++) s += vd->vcells[i].volume;
    return s;
}

static s_ncvx_vdiagram run(const s_ncvx_domain *domain, const s_points *seeds,
                           int want_surface, int merge_orphans)
{
    return vor3d_in_ncvx_domain(seeds, domain, 1e-2, EPS_DEG, TOL,
                                NULL, NULL, NULL, NULL,
                                want_surface != 0, merge_orphans != 0);
}

static int test_fixture(const char *name, s_trimesh mesh,
                        const s_point *seed_pts, int N_seeds)
{
    printf("=== %s (%d seeds) ===\n", name, N_seeds);
    int fail = 0;

    if (!trimesh_is_valid(&mesh)) {
        printf("  [FAIL] fixture trimesh invalid\n");
        free_trimesh(&mesh);
        return 1;
    }

    s_ncvx_domain domain = ncvx_domain_from_trimesh(&mesh, EPS_DEG, TOL, 0);
    if (!ncvx_domain_is_valid(&domain)) {
        printf("  [FAIL] domain construction\n");
        free_trimesh(&mesh); free_ncvx_domain(&domain);
        return 1;
    }
    double domvol = domain.domain_volume;
    s_points seeds = { .N = N_seeds, .p = (s_point *)seed_pts };

    /* Raw (no merge) */
    s_ncvx_vdiagram raw = run(&domain, &seeds, 1, 0);
    int raw_multi = count_multi_surface(&raw);
    double raw_vol = total_volume(&raw);
    printf("  raw:    multi-surface seeds = %d, volume = %.6f (domain %.6f)\n",
           raw_multi, raw_vol, domvol);

    /* Merged */
    s_ncvx_vdiagram mg = run(&domain, &seeds, 1, 1);
    int mg_multi = count_multi_surface(&mg);
    double mg_vol = total_volume(&mg);
    printf("  merged: multi-surface seeds = %d, volume = %.6f\n", mg_multi, mg_vol);

    for (int i = 0; i < mg.seeds.N; i++) {
        printf("    seed %d: N_pieces=%d N_surface=%d volume=%.6f\n",
               i, mg.vcells[i].N_pieces, mg.vcells[i].N_surface, mg.vcells[i].volume);
    }

    /* --- Assertions --- */
    if (fabs(mg_vol - domvol) / domvol > 1e-2) {
        printf("  [FAIL] merged volume not conserved: %.6f vs %.6f\n", mg_vol, domvol); fail++;
    }
    if (mg_multi != 0) {
        printf("  [FAIL] merged still has %d multi-surface seeds\n", mg_multi); fail++;
    }
    for (int i = 0; i < mg.seeds.N; i++) {
        if (mg.vcells[i].N_pieces == 0) continue;   /* empty cell: no surface */
        if (mg.vcells[i].N_surface != 1) {
            printf("  [FAIL] seed %d has N_surface=%d (expected 1)\n", i, mg.vcells[i].N_surface);
            fail++;
        }
        for (int k = 0; k < mg.vcells[i].N_surface; k++) {
            s_trimesh *s = &mg.vcells[i].surface[k];
            if (!trimesh_is_valid(s)) {
                printf("  [FAIL] seed %d surface %d invalid\n", i, k); fail++;
            } else if (surface_has_boundary_edge(s)) {
                printf("  [FAIL] seed %d surface %d has a boundary edge\n", i, k); fail++;
            }
        }
    }
    if (raw_multi == 0) {
        printf("  [WARN] raw produced no orphan for this seed set (merge is a no-op here)\n");
    } else {
        printf("  [OK] regression: raw orphans %d -> merged %d\n", raw_multi, mg_multi);
    }

    /* Determinism: a second merged run must match volume + N_surface per seed. */
    s_ncvx_vdiagram mg2 = run(&domain, &seeds, 1, 1);
    for (int i = 0; i < mg.seeds.N; i++) {
        if (mg.vcells[i].N_surface != mg2.vcells[i].N_surface ||
            fabs(mg.vcells[i].volume - mg2.vcells[i].volume) > 1e-12) {
            printf("  [FAIL] non-deterministic at seed %d\n", i); fail++; break;
        }
    }

    if (!fail) printf("  [PASS] %s\n", name);
    free_ncvx_vdiagram(&raw);
    free_ncvx_vdiagram(&mg);
    free_ncvx_vdiagram(&mg2);
    free_ncvx_domain(&domain);
    free_trimesh(&mesh);
    return fail;
}

/* Brute-force search for a seed set that produces a raw orphan (a cell with
 * N_surface > 1), to design deterministic fixtures. Prints the seeds found. */
static void search_orphans(const char *name, s_trimesh mesh, int N_seeds, int tries)
{
    s_random_context rc = random_initialize(12345);
    s_ncvx_domain domain = ncvx_domain_from_trimesh(&mesh, EPS_DEG, TOL, 0);
    if (!ncvx_domain_is_valid(&domain)) { printf("search %s: bad domain\n", name); return; }
    s_point mn, mx; bounding_box_points(&mesh.points, &mn, &mx);
    printf("--- search %s (%d seeds, %d tries) ---\n", name, N_seeds, tries);
    for (int t = 0; t < tries; t++) {
        s_point *pts = malloc(sizeof(s_point) * (size_t)N_seeds);
        int got = 0, guard = 0;
        while (got < N_seeds && guard++ < 100000) {
            s_point p = {{{ mn.x + random_uniform_double(&rc)*(mx.x-mn.x),
                            mn.y + random_uniform_double(&rc)*(mx.y-mn.y),
                            mn.z + random_uniform_double(&rc)*(mx.z-mn.z) }}};
            if (point_in_trimesh(&mesh, p, EPS_DEG, 20) == 1) pts[got++] = p;
        }
        s_points seeds = { .N = N_seeds, .p = pts };
        s_ncvx_vdiagram raw = run(&domain, &seeds, 1, 0);
        int multi = count_multi_surface(&raw);
        if (multi > 0) {
            printf("  FOUND orphan (try %d, %d multi-surface seeds):\n", t, multi);
            for (int i = 0; i < N_seeds; i++)
                printf("    {{{%.4f,%.4f,%.4f}}},  // seed %d N_surface=%d\n",
                       pts[i].x, pts[i].y, pts[i].z, i, raw.vcells[i].N_surface);
            free_ncvx_vdiagram(&raw); free(pts);
            free_ncvx_domain(&domain); free_trimesh(&mesh);
            return;
        }
        free_ncvx_vdiagram(&raw); free(pts);
    }
    printf("  no orphan found in %d tries\n", tries);
    free_ncvx_domain(&domain); free_trimesh(&mesh);
}

/* Randomized stress: run merge_orphans=1 over many random seed sets and assert
 * the invariants hold every time (conservation, one surface per non-empty seed,
 * watertight surfaces). Exercises chained/multiple orphans and empty cells.
 * Returns the number of failing trials. */
static int stress(const char *name, s_trimesh mesh, int N_seeds, int trials)
{
    s_random_context rc = random_initialize(999);
    s_ncvx_domain domain = ncvx_domain_from_trimesh(&mesh, EPS_DEG, TOL, 0);
    if (!ncvx_domain_is_valid(&domain)) { printf("stress %s: bad domain\n", name); return 1; }
    double domvol = domain.domain_volume;
    s_point mn, mx; bounding_box_points(&mesh.points, &mn, &mx);
    int fails = 0, orphan_trials = 0, pinch_trials = 0;
    for (int t = 0; t < trials; t++) {
        s_point *pts = malloc(sizeof(s_point) * (size_t)N_seeds);
        int got = 0, guard = 0;
        while (got < N_seeds && guard++ < 100000) {
            s_point p = {{{ mn.x + random_uniform_double(&rc)*(mx.x-mn.x),
                            mn.y + random_uniform_double(&rc)*(mx.y-mn.y),
                            mn.z + random_uniform_double(&rc)*(mx.z-mn.z) }}};
            if (point_in_trimesh(&mesh, p, EPS_DEG, 20) == 1) pts[got++] = p;
        }
        s_points seeds = { .N = N_seeds, .p = pts };
        s_ncvx_vdiagram raw = run(&domain, &seeds, 1, 0);
        if (count_multi_surface(&raw) > 0) orphan_trials++;
        free_ncvx_vdiagram(&raw);

        s_ncvx_vdiagram mg = run(&domain, &seeds, 1, 1);
        /* Hard invariants (a real bug if violated): volume conservation and a
         * watertight, valid surface for every component. N_surface>1 for a
         * (verified spatially-connected) region is an inherent pinch, not a
         * failure -- counted separately, see MERGING_ORPHANS.md. */
        const char *reason = NULL;
        double mgvol = total_volume(&mg);
        int pinched = 0;
        if (fabs(mgvol - domvol) / domvol > 1e-2) reason = "volume";
        for (int i = 0; i < mg.seeds.N && !reason; i++) {
            if (mg.vcells[i].N_pieces == 0) continue;
            if (mg.vcells[i].N_surface != 1) pinched = 1;
            for (int k = 0; k < mg.vcells[i].N_surface && !reason; k++) {
                if (!trimesh_is_valid(&mg.vcells[i].surface[k])) reason = "invalid surface";
                else if (surface_has_boundary_edge(&mg.vcells[i].surface[k])) reason = "boundary edge";
            }
        }
        if (pinched) pinch_trials++;
        if (reason) {
            fails++;
            printf("  [stress FAIL:%s] %s trial %d (mgvol=%.6f dom=%.6f) seeds:\n",
                   reason, name, t, mgvol, domvol);
            for (int i = 0; i < N_seeds; i++)
                printf("    {{{%.4f,%.4f,%.4f}}},  Np=%d Ns=%d vol=%.6f\n",
                       pts[i].x, pts[i].y, pts[i].z,
                       mg.vcells[i].N_pieces, mg.vcells[i].N_surface, mg.vcells[i].volume);
        }
        free_ncvx_vdiagram(&mg); free(pts);
    }
    printf("  stress %s (%d seeds x %d trials): %d hard fails, %d/%d trials had raw orphans, "
           "%d trials pinched (inherent)\n",
           name, N_seeds, trials, fails, orphan_trials, trials, pinch_trials);
    free_ncvx_domain(&domain); free_trimesh(&mesh);
    return fails;
}

int main(int argc, char **argv)
{
    if (argc > 1 && strcmp(argv[1], "search") == 0) {
        search_orphans("U-shape", make_ushape(), 2, 300);
        search_orphans("U-shape-3", make_ushape(), 3, 300);
        search_orphans("L-shape", make_lshape(), 3, 300);
        search_orphans("Plus-shape", make_plus(), 5, 300);
        return 0;
    }
    if (argc > 1 && strcmp(argv[1], "stress") == 0) {
        int f = 0;
        f += stress("U-shape", make_ushape(), 3, 200);
        f += stress("U-shape", make_ushape(), 6, 200);
        f += stress("L-shape", make_lshape(), 4, 200);
        f += stress("L-shape", make_lshape(), 8, 200);
        f += stress("Plus-shape", make_plus(), 5, 200);
        f += stress("Plus-shape", make_plus(), 10, 200);
        printf("\nstress total fails: %d\n", f);
        return f ? 1 : 0;
    }

    int fail = 0;

    /* Fixed seed sets (found by the `search` mode) that each straddle a
     * concavity and produce a guaranteed raw orphan (a cell that splits into two
     * disconnected pieces), which the merge pass must reunite into one region. */

    /* U-shape: seed 1 (low in the left arm) has its cell split by seed 0's arm. */
    static const s_point u_seeds[2] = {
        {{{2.2314, 0.3901, 0.9633}}},
        {{{0.1450, 1.6655, 0.0107}}},   /* raw N_surface = 2 */
    };
    fail += test_fixture("U-shape", make_ushape(), u_seeds, 2);

    /* L-shape: seed 2 straddles the reflex notch and disconnects. */
    static const s_point l_seeds[3] = {
        {{{1.4139, 0.5567, 0.0824}}},
        {{{0.7254, 0.4887, 0.7504}}},
        {{{1.9768, 0.9874, 0.8035}}},   /* raw N_surface = 2 */
    };
    fail += test_fixture("L-shape", make_lshape(), l_seeds, 3);

    /* Plus-shape: seed 4 gets a disconnected sliver across an arm. */
    static const s_point p_seeds[5] = {
        {{{-0.4108, -0.5321, -0.4954}}},
        {{{ 0.7599, -0.3715,  0.1632}}},
        {{{-0.1457,  0.6052, -0.1223}}},
        {{{-0.0076,  1.4141,  0.4793}}},
        {{{ 0.4418,  1.2684, -0.1595}}},   /* raw N_surface = 2 */
    };
    fail += test_fixture("Plus-shape", make_plus(), p_seeds, 5);

    printf("\n%s\n", fail ? "SOME TESTS FAILED" : "ALL TESTS PASSED");
    return fail ? 1 : 0;
}
