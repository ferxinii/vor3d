#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "vor3d.h"
#include "trimesh.h"
#include "random.h"
#include "dynarray.h"

#define N_LOBES     5
#define N_SEEDS       200   /* test 1: baseline seeds per lobe */
#define N_SEEDS_DENSE 2000  /* test 2: 10x seed density (Voronoi stress test) */

static const char *LOBE_OBJS[N_LOBES] = {
    "lobes/lung_upper_lobe_left.obj",
    "lobes/lung_lower_lobe_left.obj",
    "lobes/lung_upper_lobe_right.obj",
    "lobes/lung_middle_lobe_right.obj",
    "lobes/lung_lower_lobe_right.obj",
};

static double EPS_DEG    = 1e-9;
static double TOL        = 1e-9;
static double VOL_REL_DIFF = 1e-2;

static s_random_context rctx;

static int randint_fn(void *ctx, int N)
{
    return (int)random_uniform_range_u64(ctx, (uint64_t)N);
}

static double randd01(void *ctx)
{
    return random_uniform_double(ctx);
}

static s_points random_seeds_inside(const s_trimesh *m, int N)
{
    s_point bb_min, bb_max;
    bounding_box_points(&m->points, &bb_min, &bb_max);

    s_points seeds = { .N = N, .p = malloc(sizeof(s_point) * (size_t)N) };
    if (!seeds.p) return (s_points){0};

    int got = 0;
    while (got < N) {
        s_point p;
        p.x = bb_min.x + randd01(&rctx) * (bb_max.x - bb_min.x);
        p.y = bb_min.y + randd01(&rctx) * (bb_max.y - bb_min.y);
        p.z = bb_min.z + randd01(&rctx) * (bb_max.z - bb_min.z);
        if (point_in_trimesh(m, p, EPS_DEG, 20) == 1)
            seeds.p[got++] = p;
    }
    return seeds;
}

/* Write all lobes' diagrams into ONE vtk. Cells carry two scalars:
 * lobe_id (0..N_LOBES-1) and a globally-offset cell_id (unique across lobes). */
static void write_vtk_all(const s_ncvx_vdiagram *vds, int n_lobes, const char *path)
{
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "Cannot open %s for writing\n", path); return; }

    int total_tris = 0;
    for (int l = 0; l < n_lobes; l++)
        for (int i = 0; i < vds[l].seeds.N; i++)
            for (int p = 0; p < vds[l].vcells[i].N_pieces; p++)
                total_tris += vds[l].vcells[i].pieces[p].Nf;

    fprintf(f, "# vtk DataFile Version 2.0\n");
    fprintf(f, "Lung lobes Voronoi (%d lobes)\n", n_lobes);
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(f, "POINTS %d float\n", total_tris * 3);
    for (int l = 0; l < n_lobes; l++)
        for (int i = 0; i < vds[l].seeds.N; i++)
            for (int p = 0; p < vds[l].vcells[i].N_pieces; p++) {
                const s_convh *ch = &vds[l].vcells[i].pieces[p];
                for (int fi = 0; fi < ch->Nf; fi++)
                    for (int v = 0; v < 3; v++)
                        fprintf(f, "%f %f %f\n",
                                ch->points.p[ch->faces[fi*3+v]].x,
                                ch->points.p[ch->faces[fi*3+v]].y,
                                ch->points.p[ch->faces[fi*3+v]].z);
            }

    fprintf(f, "CELLS %d %d\n", total_tris, total_tris * 4);
    for (int t = 0; t < total_tris; t++)
        fprintf(f, "3 %d %d %d\n", t*3, t*3+1, t*3+2);

    fprintf(f, "CELL_TYPES %d\n", total_tris);
    for (int t = 0; t < total_tris; t++)
        fprintf(f, "5\n");

    fprintf(f, "CELL_DATA %d\n", total_tris);
    fprintf(f, "SCALARS cell_id int 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    int base = 0;
    for (int l = 0; l < n_lobes; l++) {
        for (int i = 0; i < vds[l].seeds.N; i++)
            for (int p = 0; p < vds[l].vcells[i].N_pieces; p++)
                for (int fi = 0; fi < vds[l].vcells[i].pieces[p].Nf; fi++)
                    fprintf(f, "%d\n", base + i);
        base += vds[l].seeds.N;
    }
    fprintf(f, "SCALARS lobe_id int 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int l = 0; l < n_lobes; l++)
        for (int i = 0; i < vds[l].seeds.N; i++)
            for (int p = 0; p < vds[l].vcells[i].N_pieces; p++)
                for (int fi = 0; fi < vds[l].vcells[i].pieces[p].Nf; fi++)
                    fprintf(f, "%d\n", l);

    fclose(f);
    printf("  Wrote %d triangles (%d cells, %d lobes) to %s\n",
           total_tris, base, n_lobes, path);
}

/* Build the CDT domain for one lobe (heavy part, done once).  Keeps the surface
 * for seed generation; *mesh_volume_out returns the reference volume. */
static s_ncvx_domain build_domain(const char *obj_path, double *mesh_volume_out)
{
    printf("\n=== %s ===\n", obj_path);
    s_trimesh mesh = trimesh_from_obj(obj_path, EPS_DEG);
    if (!trimesh_is_valid(&mesh)) {
        fprintf(stderr, "Failed to load mesh %s\n", obj_path);
        return (s_ncvx_domain){0};
    }
    double mesh_volume = volume_trimesh(&mesh);
    printf("  Vertices: %d   Faces: %d   Volume: %.4f\n",
           mesh.points.N, mesh.Nf, mesh_volume);

    clock_t t0 = clock();
    s_ncvx_domain domain = ncvx_domain_from_trimesh(&mesh, EPS_DEG, TOL, 0);
    clock_t t1 = clock();
    free_trimesh(&mesh);
    if (!ncvx_domain_is_valid(&domain)) {
        fprintf(stderr, "Failed to build domain for %s\n", obj_path);
        return (s_ncvx_domain){0};
    }
    printf("  Domain (CDT) volume: %.4f   built in %.2f s\n",
           domain.domain_volume, (double)(t1 - t0) / CLOCKS_PER_SEC);
    *mesh_volume_out = mesh_volume;
    return domain;
}

/* Run the Voronoi diagram of n_seeds random seeds inside an already-built
 * domain.  Returns a zeroed struct on failure. */
static s_ncvx_vdiagram run_voronoi(const s_ncvx_domain *domain, double mesh_volume,
                                   int n_seeds, const char *obj_path)
{
    clock_t tsg0 = clock();
    s_points seeds = random_seeds_inside(&domain->surface, n_seeds);
    clock_t tsg1 = clock();
    if (!seeds.p) {
        fprintf(stderr, "Failed to generate seeds for %s\n", obj_path);
        return (s_ncvx_vdiagram){0};
    }

    s_dynarray buff = dynarray_initialize(sizeof(s_point), 0);
    clock_t t0 = clock();
    s_ncvx_vdiagram vd = vor3d_in_ncvx_domain(&seeds, domain, VOL_REL_DIFF,
                                              EPS_DEG, TOL, &rctx,
                                              &buff, NULL, false, false);
    clock_t t1 = clock();
    dynarray_free(&buff);
    free_points(&seeds);

    if (vd.seeds.N == 0) {
        fprintf(stderr, "Voronoi computation failed for %s (n_seeds=%d)\n",
                obj_path, n_seeds);
        return (s_ncvx_vdiagram){0};
    }

    double vol_sum = 0.0, vol_min = 1e300, vol_max = 0.0;
    for (int i = 0; i < vd.seeds.N; i++) {
        double v = vd.vcells[i].volume;
        vol_sum += v;
        if (v < vol_min) vol_min = v;
        if (v > vol_max) vol_max = v;
    }
    double vol_rel_err = (vol_sum - mesh_volume) / mesh_volume;
    printf("  [%4d seeds] %d cells  vol_sum=%.4f  vol_min=%.4f  vol_max=%.4f  "
           "rel_err=%.4e   (seedgen %.2fs, voronoi %.2fs)\n",
           n_seeds, vd.seeds.N, vol_sum, vol_min, vol_max, vol_rel_err,
           (double)(tsg1 - tsg0) / CLOCKS_PER_SEC,
           (double)(t1 - t0) / CLOCKS_PER_SEC);
    return vd;
}

int main(void)
{
    rctx = random_initialize((uint64_t)time(NULL));

    /* Build each lobe's CDT domain ONCE; the Voronoi tests below reuse it. */
    s_ncvx_domain domains[N_LOBES];
    double mesh_vols[N_LOBES];
    int n_domains_ok = 0;
    for (int l = 0; l < N_LOBES; l++) {
        domains[l] = build_domain(LOBE_OBJS[l], &mesh_vols[l]);
        if (ncvx_domain_is_valid(&domains[l])) n_domains_ok++;
    }
    if (n_domains_ok != N_LOBES) {
        fprintf(stderr, "\n%d / %d domains FAILED to build\n",
                N_LOBES - n_domains_ok, N_LOBES);
        for (int l = 0; l < N_LOBES; l++)
            if (ncvx_domain_is_valid(&domains[l])) free_ncvx_domain(&domains[l]);
        return 1;
    }

    struct { const char *label; int n_seeds; const char *vtk; } tests[2] = {
        { "TEST 1 -- baseline",         N_SEEDS,       "lobes_voronoi.vtk" },
        { "TEST 2 -- 10x seed density", N_SEEDS_DENSE, "lobes_voronoi_dense.vtk" },
    };

    int total_failed = 0;
    for (int t = 0; t < 2; t++) {
        printf("\n########################################################\n");
        printf("# %s: %d seeds/lobe (%d total)\n",
               tests[t].label, tests[t].n_seeds, tests[t].n_seeds * N_LOBES);
        printf("########################################################\n");

        s_ncvx_vdiagram vds[N_LOBES];
        int failed = 0;
        for (int l = 0; l < N_LOBES; l++) {
            printf("--- %s ---\n", LOBE_OBJS[l]);
            vds[l] = run_voronoi(&domains[l], mesh_vols[l], tests[t].n_seeds,
                                 LOBE_OBJS[l]);
            if (vds[l].seeds.N == 0) failed++;
        }

        int any_ok = 0;
        for (int l = 0; l < N_LOBES; l++) if (vds[l].seeds.N > 0) any_ok = 1;
        if (any_ok) {
            clock_t tv0 = clock();
            write_vtk_all(vds, N_LOBES, tests[t].vtk);
            printf("  VTK write: %.2f s\n", (double)(clock() - tv0) / CLOCKS_PER_SEC);
        }

        for (int l = 0; l < N_LOBES; l++)
            if (vds[l].seeds.N > 0) free_ncvx_vdiagram(&vds[l]);

        if (failed) {
            fprintf(stderr, "%s: %d / %d lobes FAILED\n",
                    tests[t].label, failed, N_LOBES);
            total_failed += failed;
        } else {
            printf("%s: all %d lobes OK.\n", tests[t].label, N_LOBES);
        }
    }

    for (int l = 0; l < N_LOBES; l++) free_ncvx_domain(&domains[l]);

    if (total_failed) {
        fprintf(stderr, "\nSOME TESTS FAILED\n");
        return 1;
    }
    printf("\nAll tests OK (baseline + 10x density).\n");
    return 0;
}
