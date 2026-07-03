/* Reproducer for the coplanar-input DT crash documented in ERROR.md.
 *
 * Trigger: 8 cube-corner mesh vertices (all on z=0 or z=1 planes) followed
 * by a 5×5×5 interior grid of 125 points, fed together to construct_dt_3d.
 * Expected: crash / assertion failure at point_id=12 (first interior point
 * inserted after all 8 coplanar corners occupy indices 4–11). */

#include <stdio.h>
#include <stdlib.h>
#include "delaunay.h"
#include "vor3d.h"
#include "dynarray.h"
#include "random.h"

static int randint_fn(void *ctx, int N)
{
    return (int)random_uniform_range_u64(ctx, (uint64_t)N);
}

static void write_vdiagram_vtk(const s_vdiagram *vd, const char *path)
{
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "write_vdiagram_vtk: cannot open %s\n", path); return; }

    int total_tris = 0;
    for (int i = 0; i < vd->seeds.N; i++)
        total_tris += vd->vcells[i].pieces[0].Nf;

    fprintf(f, "# vtk DataFile Version 2.0\n");
    fprintf(f, "Voronoi diagram\n");
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(f, "POINTS %d float\n", total_tris * 3);
    for (int i = 0; i < vd->seeds.N; i++) {
        const s_convh *ch = &vd->vcells[i].pieces[0];
        for (int fi = 0; fi < ch->Nf; fi++)
            for (int v = 0; v < 3; v++)
                fprintf(f, "%f %f %f\n",
                        ch->points.p[ch->faces[fi*3+v]].coords[0],
                        ch->points.p[ch->faces[fi*3+v]].coords[1],
                        ch->points.p[ch->faces[fi*3+v]].coords[2]);
    }

    fprintf(f, "CELLS %d %d\n", total_tris, total_tris * 4);
    for (int t = 0; t < total_tris; t++)
        fprintf(f, "3 %d %d %d\n", t*3, t*3+1, t*3+2);

    fprintf(f, "CELL_TYPES %d\n", total_tris);
    for (int t = 0; t < total_tris; t++)
        fprintf(f, "5\n");   /* VTK_TRIANGLE */

    fprintf(f, "CELL_DATA %d\n", total_tris);
    fprintf(f, "SCALARS cell_id int 1\n");
    fprintf(f, "LOOKUP_TABLE default\n");
    for (int i = 0; i < vd->seeds.N; i++)
        for (int fi = 0; fi < vd->vcells[i].pieces[0].Nf; fi++)
            fprintf(f, "%d\n", i);

    fclose(f);
    printf("Wrote %d triangles (%d cells) to %s\n", total_tris, vd->seeds.N, path);
}

static void write_bpoly_vtk(const s_bpoly *bp, const char *path)
{
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "write_bpoly_vtk: cannot open %s\n", path); return; }

    int Nf = bp->convh.Nf;
    fprintf(f, "# vtk DataFile Version 2.0\n");
    fprintf(f, "Bounding polytope\n");
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(f, "POINTS %d float\n", Nf * 3);
    for (int fi = 0; fi < Nf; fi++)
        for (int v = 0; v < 3; v++)
            fprintf(f, "%f %f %f\n",
                    bp->convh.points.p[bp->convh.faces[fi*3+v]].coords[0],
                    bp->convh.points.p[bp->convh.faces[fi*3+v]].coords[1],
                    bp->convh.points.p[bp->convh.faces[fi*3+v]].coords[2]);

    fprintf(f, "CELLS %d %d\n", Nf, Nf * 4);
    for (int t = 0; t < Nf; t++)
        fprintf(f, "3 %d %d %d\n", t*3, t*3+1, t*3+2);

    fprintf(f, "CELL_TYPES %d\n", Nf);
    for (int t = 0; t < Nf; t++)
        fprintf(f, "5\n");

    fclose(f);
    printf("Wrote %d triangles to %s\n", Nf, path);
}

/* Dump all DT points -- the real seeds and the mirrors added to bound the
 * boundary cells -- to a CSV with a header, for review. */
static void write_seeds_csv(const s_points *seeds, const s_dynarray *mirrors, const char *path)
{
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "write_seeds_csv: cannot open %s\n", path); return; }
    fprintf(f, "type,x,y,z\n");
    for (int i = 0; i < seeds->N; i++)
        fprintf(f, "seed,%.17g,%.17g,%.17g\n",
                seeds->p[i].coords[0], seeds->p[i].coords[1], seeds->p[i].coords[2]);
    const s_point *m = (const s_point *)mirrors->items;
    for (unsigned i = 0; i < mirrors->N; i++)
        fprintf(f, "mirror,%.17g,%.17g,%.17g\n",
                m[i].coords[0], m[i].coords[1], m[i].coords[2]);
    fclose(f);
    printf("Wrote %d seeds + %u mirrors to %s\n", seeds->N, mirrors->N, path);
}

int main(void)
{
    /* 8 cube corners — 4 on z=0, 4 on z=1 (two perfectly coplanar layers) */
    s_point mesh_verts[8] = {
        {{{0.0, 0.0, 0.0}}},
        {{{1.0, 0.0, 0.0}}},
        {{{0.0, 1.0, 0.0}}},
        {{{1.0, 1.0, 0.0}}},
        {{{0.0, 0.0, 1.0}}},
        {{{1.0, 0.0, 1.0}}},
        {{{0.0, 1.0, 1.0}}},
        {{{1.0, 1.0, 1.0}}},
    };

    /* 5×5×5 = 125 interior grid points, away from the faces */
    const int G = 5;
    s_point grid[125];
    int idx = 0;
    for (int ix = 1; ix <= G; ix++)
        for (int iy = 1; iy <= G; iy++)
            for (int iz = 1; iz <= G; iz++)
                grid[idx++] = (s_point){{{ix / (double)(G + 1),
                                          iy / (double)(G + 1),
                                          iz / (double)(G + 1)}}};

    /* Combine: mesh vertices first (will become DT point-ids 4–11),
     * interior grid second (first at point-id 12). */
    const int N = 8 + 125;
    s_point *pts = malloc((size_t)N * sizeof(s_point));
    if (!pts) { fputs("OOM\n", stderr); return 1; }
    for (int i = 0; i < 8;   i++) pts[i]     = mesh_verts[i];
    for (int i = 0; i < 125; i++) pts[8 + i] = grid[i];

    s_points points = {.N = N, .p = pts};

    printf("Calling construct_dt_3d with %d points (%d coplanar mesh verts + %d interior)...\n",
           N, 8, 125);
    fflush(stdout);

    s_scplx dt = construct_dt_3d(&points, NULL, false, 1e-12, NULL);

    if (!dt.head) {
        fputs("construct_dt_3d returned an empty complex (error or crash).\n", stderr);
        free(pts);
        return 1;
    }

    printf("construct_dt_3d succeeded — %d tetrahedra.\n", dt.N_ncells);

    s_bpoly bp = bpoly_from_points(&points, 1e-12);
    write_bpoly_vtk(&bp, "bpoly.vtk");
    s_random_context rctx = random_initialize(42);
    s_dynarray buff = dynarray_initialize(sizeof(s_point), 0);

    s_vdiagram vd = vor3d_in_bp(&points, &bp, 1e-3, 1e-12, 1e-12,
                                    randint_fn, &rctx, &buff, NULL);
    write_seeds_csv(&points, &buff, "seeds_and_mirrors.csv");
    dynarray_free(&buff);

    if (!vdiagram_is_valid(&vd)) {
        fputs("vor3d_in_bp returned an invalid diagram.\n", stderr);
        free_vdiagram(&vd);
        free_bpoly(&bp);
        free_complex(&dt);
        free(pts);
        return 1;
    }

    printf("vor3d_in_bp succeeded — %d Voronoi cells.\n", vd.seeds.N);

    double sum_vol = 0.0;
    for (int i = 0; i < vd.seeds.N; i++)
        sum_vol += vd.vcells[i].volume;
    double bp_vol = vd.bpoly.volume;
    printf("BP volume: %.10f\nSum of vcell volumes: %.10f\nRelative error: %.6e\n",
           bp_vol, sum_vol, (bp_vol - sum_vol) / bp_vol);

    write_vdiagram_vtk(&vd, "voronoi_coplanar.vtk");

    free_vdiagram(&vd);
    free_bpoly(&bp);
    free_complex(&dt);
    free(pts);
    return 0;
}
