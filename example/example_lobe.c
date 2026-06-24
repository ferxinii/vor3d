#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "vor3d.h"
#include "trimesh.h"
#include "random.h"
#include "dynarray.h"

#define LOBE_OBJ   "lobes/lung_middle_lobe_right.obj"
#define OUTPUT_VTK "lobe_voronoi.vtk"
#define N_SEEDS    200

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

static void write_vtk(const s_ncvx_vdiagram *vd, const char *path)
{
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "Cannot open %s for writing\n", path); return; }

    int total_tris = 0;
    for (int i = 0; i < vd->seeds.N; i++)
        for (int p = 0; p < vd->vcells[i].N_pieces; p++)
            total_tris += vd->vcells[i].pieces[p].Nf;

    fprintf(f, "# vtk DataFile Version 2.0\n");
    fprintf(f, "Lung middle lobe right Voronoi\n");
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(f, "POINTS %d float\n", total_tris * 3);
    for (int i = 0; i < vd->seeds.N; i++)
        for (int p = 0; p < vd->vcells[i].N_pieces; p++) {
            const s_convh *ch = &vd->vcells[i].pieces[p];
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
    for (int i = 0; i < vd->seeds.N; i++)
        for (int p = 0; p < vd->vcells[i].N_pieces; p++)
            for (int fi = 0; fi < vd->vcells[i].pieces[p].Nf; fi++)
                fprintf(f, "%d\n", i);

    fclose(f);
    printf("Wrote %d triangles (%d cells) to %s\n", total_tris, vd->seeds.N, path);
}

int main(void)
{
    rctx = random_initialize((uint64_t)time(NULL));

    /* Load mesh */
    printf("Loading %s ...\n", LOBE_OBJ);
    s_trimesh mesh = trimesh_from_obj(LOBE_OBJ, EPS_DEG);
    if (!trimesh_is_valid(&mesh)) {
        fprintf(stderr, "Failed to load mesh\n");
        return 1;
    }
    double mesh_volume = volume_trimesh(&mesh);
    printf("  Vertices: %d   Faces: %d   Volume: %.4f\n",
           mesh.points.N, mesh.Nf, mesh_volume);

    /* Build domain (tetrahedralises interior) */
    printf("Building domain ...\n");
    s_ncvx_domain domain = ncvx_domain_from_trimesh(&mesh, EPS_DEG, TOL);
    free_trimesh(&mesh);
    if (!ncvx_domain_is_valid(&domain)) {
        fprintf(stderr, "Failed to build domain\n");
        return 1;
    }
    printf("  Domain volume (tet sum): %.4f   N_tets: %d\n",
           domain.domain_volume, domain.N_tets);

    /* Random seeds inside the lobe */
    printf("Generating %d seeds ...\n", N_SEEDS);
    s_points seeds = random_seeds_inside(&domain.surface, N_SEEDS);
    if (!seeds.p) { fprintf(stderr, "Failed to generate seeds\n"); free_ncvx_domain(&domain); return 1; }

    /* Voronoi diagram */
    printf("Computing Voronoi diagram ...\n");
    s_dynarray buff = dynarray_initialize(sizeof(s_point), 0);
    s_ncvx_vdiagram vd = vor3d_in_ncvx_domain(&seeds, &domain, VOL_REL_DIFF,
                                               EPS_DEG, TOL, randint_fn, &rctx,
                                               &buff, NULL);
    dynarray_free(&buff);
    free_points(&seeds);
    free_ncvx_domain(&domain);

    if (vd.seeds.N == 0) {
        fprintf(stderr, "Voronoi computation failed\n");
        return 1;
    }

    /* Print stats */
    double vol_sum = 0.0, vol_min = 1e300, vol_max = 0.0;
    for (int i = 0; i < vd.seeds.N; i++) {
        double v = vd.vcells[i].volume;
        vol_sum += v;
        if (v < vol_min) vol_min = v;
        if (v > vol_max) vol_max = v;
    }
    double vol_rel_err = (vol_sum - mesh_volume) / mesh_volume;
    printf("Done: %d cells   vol_sum=%.4f   vol_min=%.4f   vol_max=%.4f   rel_err=%.4e\n",
           vd.seeds.N, vol_sum, vol_min, vol_max, vol_rel_err);

    write_vtk(&vd, OUTPUT_VTK);
    free_ncvx_vdiagram(&vd);
    return 0;
}
