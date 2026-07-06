#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include "vor3d.h"
#include "random.h"
#include "dynarray.h"

static double EPS_DEG = 1e-12;
static double TOL     = 1e-12;
static double VOL_REL_DIFF = 1e-2;   /* looser than convex: accounts for tet discretisation */

static s_random_context rctx;

static int randint_fn(void *ctx, int N)
{
    return (int)random_uniform_range_u64(ctx, (uint64_t)N);
}

static double randd01(void *ctx)
{
    return random_uniform_double(ctx);
}

/* -----------------------------------------------------------------------
 * Domain definitions
 * ----------------------------------------------------------------------- */

/* --- Cube: 2x2x2 centred at origin, volume = 8 --- */
// static s_trimesh make_cube(void)
// {
//     /* 8 vertices */
//     static const s_point pts[8] = {
//         {{{-1,-1,-1}}}, {{{ 1,-1,-1}}}, {{{ 1, 1,-1}}}, {{{-1, 1,-1}}},
//         {{{-1,-1, 1}}}, {{{ 1,-1, 1}}}, {{{ 1, 1, 1}}}, {{{-1, 1, 1}}},
//     };
//     /* 12 triangles (2 per face, outward winding — trimesh_from_arrays will correct if needed) */
//     static const int faces[12*3] = {
//         0,1,2, 0,2,3,   /* -z */
//         4,6,5, 4,7,6,   /* +z */
//         0,4,5, 0,5,1,   /* -y */
//         2,6,7, 2,7,3,   /* +y */
//         0,3,7, 0,7,4,   /* -x */
//         1,5,6, 1,6,2,   /* +x */
//     };
//     return trimesh_from_arrays(pts, 8, faces, 12, EPS_DEG);
// }

/* --- L-shape: 2x2x1 box minus a 1x1x1 corner, volume = 3 ---
 * Occupies x in [0,2], y in [0,2], z in [0,1].
 * Removed corner: x in [1,2], y in [1,2], z in [0,1]. */
static s_trimesh make_lshape(void)
{
    static const s_point verts[12] = {
        /* z=0 */
        {{{0,0,0}}}, {{{2,0,0}}}, {{{2,1,0}}}, {{{1,1,0}}}, {{{1,2,0}}}, {{{0,2,0}}},
        /* z=1 */
        {{{0,0,1}}}, {{{2,0,1}}}, {{{2,1,1}}}, {{{1,1,1}}}, {{{1,2,1}}}, {{{0,2,1}}},
    };
    /* Faces: bottom (z=0, inward normal = +z so winding CW from above),
     *        top (z=1, outward normal = +z),
     *        5 rectangular sides split into 2 triangles each. */
    static const int faces[20*3] = {
        /* bottom (z=0, normal points -z: CCW from below = CW from above) */
        0,5,4,  0,4,3,  0,3,2,  0,2,1,
        /* top (z=1, normal points +z) */
        6,7,8,  6,8,9,  6,9,10, 6,10,11,
        /* sides */
        0,1,7,  0,7,6,    /* y=0 */
        1,2,8,  1,8,7,    /* x=2, y in [0,1] */
        2,3,9,  2,9,8,    /* concave step: y=1, x in [1,2] */
        3,4,10, 3,10,9,   /* concave step: x=1, y in [1,2] */
        4,5,11, 4,11,10,  /* y=2 */
        5,0,6,  5,6,11,   /* x=0 */
    };
    return trimesh_from_arrays(verts, 12, faces, 20, EPS_DEG);
}

/* --- Plus-shape: two 3x1x1 boxes crossing at centre, volume = 5 ---
 * Box A: x in [-1.5,1.5], y in [-0.5,0.5], z in [-0.5,0.5]
 * Box B: x in [-0.5,0.5], y in [-1.5,1.5], z in [-0.5,0.5]
 * Union: volume = 3 + 3 - 1 = 5 */
static s_trimesh make_plus(void)
{
    /* Plus/cross shape: union of
     *   Box A: x in [-1.5,1.5], y in [-0.5,0.5], z in [-0.5,0.5]  (3x1x1)
     *   Box B: x in [-0.5,0.5], y in [-1.5,1.5], z in [-0.5,0.5]  (1x3x1)
     * Volume = 3 + 3 - 1 = 5.
     * Cross-section: 12-vertex non-convex polygon with 4 reflex corners
     * at v2, v5, v8, v11 (and their z=+0.5 counterparts v14,v17,v20,v23).
     * Vertices are listed CCW from above for z=-0.5. */
    static const s_point verts[24] = {
        /* z = -0.5 (v0..v11, CCW from above) */
        {{{-0.5,-1.5,-0.5}}}, {{{ 0.5,-1.5,-0.5}}},   /* v0,  v1  */
        {{{ 0.5,-0.5,-0.5}}}, {{{ 1.5,-0.5,-0.5}}},   /* v2,  v3  (v2 reflex) */
        {{{ 1.5, 0.5,-0.5}}}, {{{ 0.5, 0.5,-0.5}}},   /* v4,  v5  (v5 reflex) */
        {{{ 0.5, 1.5,-0.5}}}, {{{-0.5, 1.5,-0.5}}},   /* v6,  v7  */
        {{{-0.5, 0.5,-0.5}}}, {{{-1.5, 0.5,-0.5}}},   /* v8,  v9  (v8 reflex) */
        {{{-1.5,-0.5,-0.5}}}, {{{-0.5,-0.5,-0.5}}},   /* v10, v11 (v11 reflex) */
        /* z = +0.5 (v12..v23, same xy shifted by 12) */
        {{{-0.5,-1.5, 0.5}}}, {{{ 0.5,-1.5, 0.5}}},   /* v12, v13 */
        {{{ 0.5,-0.5, 0.5}}}, {{{ 1.5,-0.5, 0.5}}},   /* v14, v15 */
        {{{ 1.5, 0.5, 0.5}}}, {{{ 0.5, 0.5, 0.5}}},   /* v16, v17 */
        {{{ 0.5, 1.5, 0.5}}}, {{{-0.5, 1.5, 0.5}}},   /* v18, v19 */
        {{{-0.5, 0.5, 0.5}}}, {{{-1.5, 0.5, 0.5}}},   /* v20, v21 */
        {{{-1.5,-0.5, 0.5}}}, {{{-0.5,-0.5, 0.5}}},   /* v22, v23 */
    };
    /* 44 triangles: 10 bottom + 10 top + 24 sides (12 rectangles x 2).
     * Bottom/top decomposed into 5 axis-aligned rectangles x 2 triangles.
     * Bottom: CW from above (outward normal = -z).
     * Top:    CCW from above (outward normal = +z).
     * Sides:  (vi, v_{i+1}, v_{i+13}), (vi, v_{i+13}, v_{i+12}) = outward. */
    static const int faces[44*3] = {
        /* bottom (outward = -z, CW from above) */
        2,1,0,   11,2,0,     /* bottom arm:  v0,v1,v2,v11 */
        4,3,2,   5,4,2,      /* right arm:   v2,v3,v4,v5  */
        7,6,5,   8,7,5,      /* top arm:     v5,v6,v7,v8  */
        10,9,8,  11,10,8,    /* left arm:    v8,v9,v10,v11 */
        5,2,11,  8,5,11,     /* center:      v2,v5,v8,v11 */
        /* top (outward = +z, CCW from above) */
        12,13,14, 12,14,23,  /* bottom arm:  v12,v13,v14,v23 */
        14,15,16, 14,16,17,  /* right arm:   v14,v15,v16,v17 */
        17,18,19, 17,19,20,  /* top arm:     v17,v18,v19,v20 */
        20,21,22, 20,22,23,  /* left arm:    v20,v21,v22,v23 */
        23,14,17, 23,17,20,  /* center:      v14,v17,v20,v23 */
        /* 12 side rectangles: edge (vi -> v_{i+1 mod 12}) */
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


/* -----------------------------------------------------------------------
 * Test helpers
 * ----------------------------------------------------------------------- */

static void check_trimesh(const char *name, const s_trimesh *m, int expected_faces,
                           int *fail, int *total)
{
    (*total)++;
    int ok = 1;

    if (!trimesh_is_valid(m)) { printf("  [FAIL] %s: trimesh_is_valid\n", name); ok=0; }
    if (m->Nf != expected_faces) {
        printf("  [FAIL] %s: Nf=%d expected %d\n", name, m->Nf, expected_faces); ok=0;
    }
    for (int i = 0; i < m->Nf * 3 && ok; i++) {
        if (m->adjacency[i] == -1) {
            printf("  [FAIL] %s: boundary edge at face %d local %d\n",
                   name, i/3, i%3); ok=0;
        }
    }
    /* Outward normal sanity: divergence-theorem signed volume must be positive.
     * A per-face centroid dot-product check fails on non-convex meshes (concave
     * faces can have a centroid on the "wrong" side of the mesh centroid even
     * when the normal is correctly outward).  Signed volume is the right test. */
    {
        s_point c = point_average(&m->points);
        double V = 0.0;
        for (int i = 0; i < m->Nf; i++) {
            s_point v0 = subtract_points(m->points.p[m->faces[i*3+0]], c);
            s_point v1 = subtract_points(m->points.p[m->faces[i*3+1]], c);
            s_point v2 = subtract_points(m->points.p[m->faces[i*3+2]], c);
            V += dot_prod(v0, cross_prod(v1, v2));
        }
        V /= 6.0;
        if (V <= 0.0) {
            printf("  [FAIL] %s: signed volume %.6g <= 0 (normals point inward)\n",
                   name, V); ok=0;
        }
    }

    if (!ok) (*fail)++;
    else printf("  [OK] %s trimesh construction\n", name);
}

static void check_point_in_trimesh(const char *name, const s_trimesh *m,
                                    const s_point *inside_pts, int N_in,
                                    const s_point *outside_pts, int N_out,
                                    int *fail, int *total)
{
    (*total)++;
    int ok = 1;
    for (int i = 0; i < N_in; i++) {
        int r = point_in_trimesh(m, inside_pts[i], EPS_DEG, 20);
        if (r != 1) {
            printf("  [FAIL] %s: inside point %d returned %d\n", name, i, r); ok=0;
        }
    }
    for (int i = 0; i < N_out; i++) {
        int r = point_in_trimesh(m, outside_pts[i], EPS_DEG, 20);
        if (r != 0) {
            printf("  [FAIL] %s: outside point %d returned %d\n", name, i, r); ok=0;
        }
    }
    if (!ok) (*fail)++;
    else printf("  [OK] %s point_in_trimesh\n", name);
}

static void write_ncvx_vdiagram_vtk(const s_ncvx_vdiagram *vd, const char *path)
{
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "write_ncvx_vdiagram_vtk: cannot open %s\n", path); return; }

    int total_tris = 0;
    for (int i = 0; i < vd->seeds.N; i++)
        for (int p = 0; p < vd->vcells[i].N_pieces; p++)
            total_tris += vd->vcells[i].pieces[p].Nf;

    fprintf(f, "# vtk DataFile Version 2.0\n");
    fprintf(f, "Non-convex Voronoi diagram\n");
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(f, "POINTS %d float\n", total_tris * 3);
    for (int i = 0; i < vd->seeds.N; i++)
        for (int p = 0; p < vd->vcells[i].N_pieces; p++) {
            const s_convh *ch = &vd->vcells[i].pieces[p];
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
        for (int p = 0; p < vd->vcells[i].N_pieces; p++)
            for (int fi = 0; fi < vd->vcells[i].pieces[p].Nf; fi++)
                fprintf(f, "%d\n", i);

    fclose(f);
    printf("  Wrote %d triangles (%d cells) to %s\n", total_tris, vd->seeds.N, path);
}


static void check_domain_volume(const char *name, const s_ncvx_domain *domain,
                                 double analytical_vol, int *fail, int *total)
{
    (*total)++;
    double rel_err = fabs(domain->domain_volume - analytical_vol) / analytical_vol;
    printf("  %s domain volume: analytical=%.4f tet_sum=%.4f rel_err=%.2e\n",
           name, analytical_vol, domain->domain_volume, rel_err);
    /* 1% tolerance: interior tet mesh only approximates the domain volume */
    if (rel_err > 0.01) {
        printf("  [FAIL] %s: domain volume rel_err=%.2e > 1%%\n", name, rel_err);
        (*fail)++;
    } else {
        printf("  [OK] %s domain volume\n", name);
    }
}

static void check_ncvx_vdiagram(const char *name, const s_ncvx_vdiagram *vd,
                                  double analytical_vol, int *fail, int *total)
{
    (*total)++;
    int ok = 1;
    double sum = 0;
    for (int i = 0; i < vd->seeds.N; i++) {
        if (vd->vcells[i].N_pieces < 1) {
            printf("  [FAIL] %s: cell %d has 0 pieces\n", name, i); ok=0;
        }
        if (vd->vcells[i].volume <= 0) {
            printf("  [FAIL] %s: cell %d has volume <= 0\n", name, i); ok=0;
        }
        sum += vd->vcells[i].volume;
    }
    double rel_err = fabs(sum - analytical_vol) / analytical_vol;
    printf("  %s Voronoi volume: analytical=%.4f sum=%.4f rel_err=%.2e\n",
           name, analytical_vol, sum, rel_err);
    if (rel_err > VOL_REL_DIFF) {
        printf("  [FAIL] %s: Voronoi volume rel_err=%.2e\n", name, rel_err); ok=0;
    }
    if (!ok) (*fail)++;
    else printf("  [OK] %s Voronoi run\n", name);
}

/* Generate N random seeds strictly inside the trimesh using rejection sampling */
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

static void iterate_nonconvex(const char *name, int iterations, int N_seeds,
                               const s_ncvx_domain *domain, double analytical_vol,
                               int *fail, int *total)
{
    int local_fail = 0;
    for (int ii = 0; ii < iterations; ii++) {
        s_points seeds = random_seeds_inside(&domain->surface, N_seeds);
        if (!seeds.p) { local_fail++; continue; }

        s_dynarray buff = dynarray_initialize(sizeof(s_point), 0);
        s_ncvx_vdiagram vd = vor3d_in_ncvx_domain(&seeds, domain, VOL_REL_DIFF,
                                                   EPS_DEG, TOL, randint_fn, &rctx,
                                                   &buff, NULL, 0);
        dynarray_free(&buff);
        free_points(&seeds);

        if (vd.seeds.N == 0) { local_fail++; printf("  iter %d: EMPTY\n", ii); continue; }

        double sum = 0;
        for (int k = 0; k < vd.seeds.N; k++) sum += vd.vcells[k].volume;
        double rel_err = fabs(sum - analytical_vol) / analytical_vol;
        int failed = rel_err > VOL_REL_DIFF;
        if (failed) local_fail++;
        printf("  iter %2d: sum=%.6f rel_err=%.2e %s\n",
               ii, sum, rel_err, failed ? "[FAIL]" : "[OK]");

        free_ncvx_vdiagram(&vd);
    }
    printf("  %s stress (%d seeds): FAILED %d / %d\n",
           name, N_seeds, local_fail, iterations);
    (*total)++;
    if (local_fail > 0) (*fail)++;
}


/* -----------------------------------------------------------------------
 * Per-domain test block
 * ----------------------------------------------------------------------- */

static void test_domain(const char *name, s_trimesh mesh, int expected_faces,
                         double analytical_vol,
                         const s_point *inside_pts, int N_in,
                         const s_point *outside_pts, int N_out,
                         int stress_iters, int stress_seeds,
                         int *fail, int *total)
{
    printf("\n=== %s ===\n", name);

    /* 1. Trimesh construction */
    check_trimesh(name, &mesh, expected_faces, fail, total);

    /* 2. point_in_trimesh */
    check_point_in_trimesh(name, &mesh, inside_pts, N_in, outside_pts, N_out,
                            fail, total);

    /* 3. Domain preparation */
    s_ncvx_domain domain = ncvx_domain_from_trimesh(&mesh, EPS_DEG, TOL);
    fprintf(stderr, "DEBUG: ncvx_domain_from_trimesh returned, valid=%d\n",
            ncvx_domain_is_valid(&domain));
    if (!ncvx_domain_is_valid(&domain)) {
        printf("  [FAIL] %s: ncvx_domain_from_trimesh returned invalid domain\n", name);
        (*fail)++; (*total)++;
        free_trimesh(&mesh);
        return;
    }
    check_domain_volume(name, &domain, analytical_vol, fail, total);

    /* 4. Single Voronoi run */
    {
        s_points seeds = random_seeds_inside(&mesh, 8);
        s_dynarray buff = dynarray_initialize(sizeof(s_point), 0);
        s_ncvx_vdiagram vd = vor3d_in_ncvx_domain(&seeds, &domain, VOL_REL_DIFF,
                                               EPS_DEG, TOL, randint_fn, &rctx,
                                               &buff, NULL, 0);
        dynarray_free(&buff);
        free_points(&seeds);
        if (vd.seeds.N == 0) {
            printf("  [FAIL] %s: single Voronoi run returned empty diagram\n", name);
            (*fail)++; (*total)++;
        } else {
            check_ncvx_vdiagram(name, &vd, domain.domain_volume, fail, total);
            char vtk_path[256];
            snprintf(vtk_path, sizeof(vtk_path), "./%s_voronoi.vtk", name);
            write_ncvx_vdiagram_vtk(&vd, vtk_path);
        }
        free_ncvx_vdiagram(&vd);
    }

    /* 5. Stress iteration */
    iterate_nonconvex(name, stress_iters, stress_seeds, &domain,
                      domain.domain_volume, fail, total);

    free_ncvx_domain(&domain);
    free_trimesh(&mesh);
}


/* -----------------------------------------------------------------------
 * Degenerate seed checks
 * ----------------------------------------------------------------------- */

// static void test_degenerate(const s_ncvx_domain *domain, double domain_vol,
//                              int *fail, int *total)
// {
//     printf("\n=== Degenerate seed configurations ===\n");
//
//     /* Seed on the concave edge of the L-shape: (1, 1, 0.5) */
//     {
//         (*total)++;
//         s_point edge_seed = {{{1.0, 1.0, 0.5}}};
//         s_point other = {{{0.5, 0.5, 0.5}}};
//         s_point pts[2] = {edge_seed, other};
//         s_points seeds = { .N = 2, .p = pts };
//         s_dynarray buff = dynarray_initialize(sizeof(s_point), 0);
//         s_ncvx_vdiagram vd = vor3d_nonconvex(&seeds, domain, VOL_REL_DIFF,
//                                                EPS_DEG, TOL, randint_fn, &rctx,
//                                                &buff, NULL);
//         dynarray_free(&buff);
//         int ok = 1;
//         if (vd.seeds.N == 0) { printf("  [FAIL] degenerate: empty diagram\n"); ok=0; }
//         else {
//             double sum = 0;
//             for (int i = 0; i < vd.seeds.N; i++) {
//                 sum += vd.vcells[i].volume;
//                 if (vd.vcells[i].volume <= 0 || vd.vcells[i].N_pieces < 1)
//                     { printf("  [FAIL] degenerate: cell %d invalid\n", i); ok=0; }
//             }
//             double rel = fabs(sum - domain_vol) / domain_vol;
//             if (rel > VOL_REL_DIFF)
//                 { printf("  [FAIL] degenerate: volume rel_err=%.2e\n", rel); ok=0; }
//         }
//         free_ncvx_vdiagram(&vd);
//         if (!ok) (*fail)++;
//         else printf("  [OK] seed on concave edge\n");
//     }
//
//     /* Symmetric seeds: placed symmetrically inside the L so Voronoi boundaries
//      * may align with trimesh faces */
//     {
//         (*total)++;
//         s_point pts[4] = {
//             {{{0.5, 0.5, 0.5}}}, {{{1.5, 0.5, 0.5}}},
//             {{{0.5, 1.5, 0.5}}}, {{{0.5, 0.5, 0.5}}},  /* intentional symmetry */
//         };
//         /* use only 3 distinct */
//         pts[3] = (s_point){{{0.25, 0.25, 0.5}}};
//         s_points seeds = { .N = 4, .p = pts };
//         s_dynarray buff = dynarray_initialize(sizeof(s_point), 0);
//         s_ncvx_vdiagram vd = vor3d_nonconvex(&seeds, domain, VOL_REL_DIFF,
//                                                EPS_DEG, TOL, randint_fn, &rctx,
//                                                &buff, NULL);
//         dynarray_free(&buff);
//         int ok = 1;
//         if (vd.seeds.N == 0) { printf("  [FAIL] symmetric: empty diagram\n"); ok=0; }
//         else {
//             double sum = 0;
//             for (int i = 0; i < vd.seeds.N; i++) sum += vd.vcells[i].volume;
//             double rel = fabs(sum - domain_vol) / domain_vol;
//             if (rel > VOL_REL_DIFF)
//                 { printf("  [FAIL] symmetric: volume rel_err=%.2e\n", rel); ok=0; }
//         }
//         free_ncvx_vdiagram(&vd);
//         if (!ok) (*fail)++;
//         else printf("  [OK] symmetric seeds\n");
//     }
// }


/* -----------------------------------------------------------------------
 * Main
 * ----------------------------------------------------------------------- */

int main(void)
{
    rctx = random_initialize((uint64_t)time(NULL));

    int fail = 0, total = 0;

    /* --- Cube (convex reference) --- */
    // {
    //     static const s_point inside[] = { {{{0,0,0}}}, {{{0.5,0.5,0.5}}}, {{{-0.5,-0.5,0}}} };
    //     static const s_point outside[] = { {{{2,0,0}}}, {{{0,0,2}}}, {{{1.5,1.5,1.5}}} };
    //     test_domain("Cube", make_cube(), 12, 8.0,
    //                 inside, 3, outside, 3,
    //                 0.4, 10, 6, &fail, &total);
    // }

    /* --- L-shape --- */
    {
        static const s_point inside[] = {
            {{{0.5,0.5,0.5}}}, {{{1.5,0.5,0.5}}}, {{{0.5,1.5,0.5}}},
        };
        static const s_point outside[] = {
            {{{1.5,1.5,0.5}}},   /* in the removed corner */
            {{{3.0,0.0,0.5}}},   /* outside entirely */
            {{{0.5,0.5,2.0}}},   /* above */
        };
        s_trimesh lmesh = make_lshape();
        test_domain("L-shape", lmesh, 20, 3.0,
                    inside, 3, outside, 3,
                    10, 6, &fail, &total);

        /* Degenerate checks use a fresh L-shape domain */
        // s_trimesh lmesh2 = make_lshape();
        // s_ncvx_domain ldom = ncvx_domain_from_trimesh(&lmesh2, EPS_DEG, TOL);
        // if (ncvx_domain_is_valid(&ldom))
        //     test_degenerate(&ldom, ldom.domain_volume, &fail, &total);
        // free_ncvx_domain(&ldom);
        // free_trimesh(&lmesh2);
    }

    /* --- Plus-shape --- */
    {
        static const s_point inside[] = {
            {{{0,0,0}}}, {{{1.0,0,0}}}, {{{0,1.0,0}}},
        };
        static const s_point outside[] = {
            {{{1.0,1.0,0}}},    /* in a concave notch */
            {{{2.0,2.0,0}}},    /* outside */
            {{{0,0,1.0}}},      /* above */
        };
        test_domain("Plus-shape", make_plus(), 44, 5.0,
                    inside, 3, outside, 3,
                    10, 6, &fail, &total);
    }

    printf("\n============================\n");
    printf("TOTAL: FAILED %d / %d\n", fail, total);
    printf("============================\n");
    return fail > 0 ? 1 : 0;
}
