/* Medial axis extraction (medax.h).
 * For each test shape (jittered sphere, L-shape, plus, and the 5 lung lobes)
 * the inner medial axis is extracted together with its medial graph; the
 * example prints the point/edge counts and the deepest medial ball, and writes
 * <shape>_surface.vtk and <shape>_axis.vtk for viewing in ParaView.
 *
 * Run from the example/ build dir:  ./example_medax [r_scale]
 *   r_scale (default 1.0) multiplies every shape's sampling spacing r_sample
 *   (smaller = denser sampling, more medial points).
 *
 * The L-shape / plus builders are copied from example_nonconvex.c. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "trimesh.h"
#include "points.h"
#include "medax.h"
#include "random.h"

static double TOL     = 1e-9;

static s_random_context rctx;
static double randd01(void *ctx) { return random_uniform_double(ctx); }

/* ---------------- L-shape (from example_nonconvex.c) --------------------- */
/* NOTE: kept axis-aligned on purpose. The medial-axis inside/outside test
 * uses a winding-number oracle (not +x ray-casting), so axis-aligned meshes
 * work directly -- no rotation jitter needed. */
static s_trimesh make_lshape(void)
{
    static const s_point verts[12] = {
        {{{0,0,0}}}, {{{2,0,0}}}, {{{2,1,0}}}, {{{1,1,0}}}, {{{1,2,0}}}, {{{0,2,0}}},
        {{{0,0,1}}}, {{{2,0,1}}}, {{{2,1,1}}}, {{{1,1,1}}}, {{{1,2,1}}}, {{{0,2,1}}},
    };
    static const int faces[20*3] = {
        0,5,4,  0,4,3,  0,3,2,  0,2,1,
        6,7,8,  6,8,9,  6,9,10, 6,10,11,
        0,1,7,  0,7,6,   1,2,8,  1,8,7,   2,3,9,  2,9,8,
        3,4,10, 3,10,9,  4,5,11, 4,11,10, 5,0,6,  5,6,11,
    };
    return trimesh_from_arrays(verts, 12, faces, 20);
}

/* ---------------- Plus-shape (from example_nonconvex.c) ------------------ */
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
        2,1,0,   11,2,0,   4,3,2,   5,4,2,   7,6,5,   8,7,5,
        10,9,8,  11,10,8,  5,2,11,  8,5,11,
        12,13,14, 12,14,23, 14,15,16, 14,16,17, 17,18,19, 17,19,20,
        20,21,22, 20,22,23, 23,14,17, 23,17,20,
        0,1,13,  0,13,12,  1,2,14,  1,14,13,  2,3,15,  2,15,14,
        3,4,16,  3,16,15,  4,5,17,  4,17,16,  5,6,18,  5,18,17,
        6,7,19,  6,19,18,  7,8,20,  7,20,19,  8,9,21,  8,21,20,
        9,10,22, 9,22,21,  10,11,23, 10,23,22, 11,0,12, 11,12,23,
    };
    return trimesh_from_arrays(verts, 24, faces, 44);
}

/* ---------------- jittered icosphere ------------------------------------ */
typedef struct { int *a, *b, *m, n, cap; } s_midcache;

static int get_mid(s_midcache *c, int a, int b, s_point *V, int *nV)
{
    if (a > b) { int t = a; a = b; b = t; }
    for (int i = 0; i < c->n; i++)
        if (c->a[i] == a && c->b[i] == b) return c->m[i];
    s_point p = {{{ (V[a].x + V[b].x) * 0.5,
                    (V[a].y + V[b].y) * 0.5,
                    (V[a].z + V[b].z) * 0.5 }}};
    double L = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
    p.x /= L; p.y /= L; p.z /= L;                 /* project onto unit sphere */
    int idx = (*nV)++;
    V[idx] = p;
    if (c->n == c->cap) {
        c->cap = c->cap ? c->cap * 2 : 128;
        c->a = realloc(c->a, sizeof(int) * c->cap);
        c->b = realloc(c->b, sizeof(int) * c->cap);
        c->m = realloc(c->m, sizeof(int) * c->cap);
    }
    c->a[c->n] = a; c->b[c->n] = b; c->m[c->n] = idx; c->n++;
    return idx;
}

/* Icosahedron subdivided nsub times, radius R, with per-vertex outward radial
 * jitter (fraction of R) so the samples are NOT cospherical (the degenerate
 * case the paper excludes). */
static s_trimesh make_sphere(int nsub, double R, double jitter)
{
    const double t = (1.0 + sqrt(5.0)) / 2.0;
    s_point ico[12] = {
        {{{-1, t, 0}}}, {{{1, t, 0}}}, {{{-1,-t, 0}}}, {{{1,-t, 0}}},
        {{{0,-1, t}}}, {{{0, 1, t}}}, {{{0,-1,-t}}}, {{{0, 1,-t}}},
        {{{t, 0,-1}}}, {{{t, 0, 1}}}, {{{-t, 0,-1}}}, {{{-t, 0, 1}}},
    };
    int ico_f[20*3] = {
        0,11,5, 0,5,1, 0,1,7, 0,7,10, 0,10,11,
        1,5,9, 5,11,4, 11,10,2, 10,7,6, 7,1,8,
        3,9,4, 3,4,2, 3,2,6, 3,6,8, 3,8,9,
        4,9,5, 2,4,11, 6,2,10, 8,6,7, 9,8,1,
    };

    int maxV = 10 * (1 << (2 * nsub)) + 2;
    int maxF = 20 * (1 << (2 * nsub));
    s_point *V = malloc(sizeof(s_point) * maxV);
    int *F = malloc(sizeof(int) * 3 * maxF);
    int nV = 12, nF = 20;
    for (int i = 0; i < 12; i++) {
        double L = sqrt(ico[i].x*ico[i].x + ico[i].y*ico[i].y + ico[i].z*ico[i].z);
        V[i] = (s_point){{{ ico[i].x/L, ico[i].y/L, ico[i].z/L }}};
    }
    for (int i = 0; i < 20*3; i++) F[i] = ico_f[i];

    for (int s = 0; s < nsub; s++) {
        s_midcache c = {0};
        int *F2 = malloc(sizeof(int) * 3 * 4 * nF);
        int nF2 = 0;
        for (int f = 0; f < nF; f++) {
            int v0 = F[3*f], v1 = F[3*f+1], v2 = F[3*f+2];
            int a = get_mid(&c, v0, v1, V, &nV);
            int b = get_mid(&c, v1, v2, V, &nV);
            int d = get_mid(&c, v2, v0, V, &nV);
            int tri[4][3] = { {v0,a,d}, {v1,b,a}, {v2,d,b}, {a,b,d} };
            for (int k = 0; k < 4; k++) {
                F2[3*nF2+0] = tri[k][0];
                F2[3*nF2+1] = tri[k][1];
                F2[3*nF2+2] = tri[k][2];
                nF2++;
            }
        }
        free(c.a); free(c.b); free(c.m); free(F);
        F = F2; nF = nF2;
    }

    /* radial jitter + scale to radius R (vertices are currently unit) */
    for (int i = 0; i < nV; i++) {
        double r = R * (1.0 + jitter * randd01(&rctx));
        V[i].x *= r; V[i].y *= r; V[i].z *= r;
    }

    s_trimesh m = trimesh_from_arrays(V, nV, F, nF);
    free(V); free(F);
    return m;
}

/* Quick node/edge report from the CSR (no traversal needed). */
static void graph_stats(const s_medial_graph *g, int *nE)
{ *nE = g->adj_off[g->N] / 2; }

/* Write the mesh surface as VTK polydata (for ParaView). */
static void write_surface_vtk(const s_trimesh *m, const char *path)
{
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); return; }
    fprintf(f, "# vtk DataFile Version 3.0\nsurface\nASCII\nDATASET POLYDATA\n");
    fprintf(f, "POINTS %d double\n", m->points.N);
    for (int i = 0; i < m->points.N; i++)
        fprintf(f, "%.17g %.17g %.17g\n",
                m->points.p[i].x, m->points.p[i].y, m->points.p[i].z);
    fprintf(f, "POLYGONS %d %d\n", m->Nf, 4 * m->Nf);
    for (int i = 0; i < m->Nf; i++)
        fprintf(f, "3 %d %d %d\n", m->faces[3*i], m->faces[3*i+1], m->faces[3*i+2]);
    fclose(f);
}


/* Extract, report and write one shape. Frees the mesh. */
static void run_shape(const char *shape, s_trimesh mesh, double r_sample)
{
    if (!trimesh_is_valid(&mesh)) { printf("%-8s INVALID mesh\n", shape); return; }

    char path[256];
    snprintf(path, sizeof path, "%s_surface.vtk", shape);
    write_surface_vtk(&mesh, path);

    /* theta0 = rho0 = 0: the paper's default filters. */
    s_medax ma = medax_from_trimesh(&mesh, r_sample, 0, 0,
                                    &(s_vor3d_opts){ .tol_dup = TOL, .randd01 = randd01, .rctx = &rctx },
                                    /*want_graph=*/true);
    if (!ma.verts.p || !ma.graph) {
        printf("%-8s extraction FAILED\n", shape);
        free_trimesh(&mesh); return;
    }

    int nE = 0; graph_stats(ma.graph, &nE);
    int deep = 0;
    for (int i = 1; i < ma.verts.N; i++) if (ma.radius[i] > ma.radius[deep]) deep = i;
    s_point c = ma.verts.p[deep];
    printf("%-8s r_sample=%-9.3g points=%-7d edges=%-8d deepest ball: r=%.4g at (%.4g, %.4g, %.4g)\n",
           shape, r_sample, ma.verts.N, nE, ma.radius[deep], c.x, c.y, c.z);

    snprintf(path, sizeof path, "%s_axis.vtk", shape);
    medax_write_vtk(&ma, path);
    free_medax(&ma);
    free_trimesh(&mesh);
}

int main(int argc, char **argv)
{
    rctx = random_initialize(12345);
    double s = (argc > 1) ? atof(argv[1]) : 1.0;

    static const char *lobe_obj[5] = {
        "lobes/lung_upper_lobe_left.obj",  "lobes/lung_lower_lobe_left.obj",
        "lobes/lung_upper_lobe_right.obj", "lobes/lung_middle_lobe_right.obj",
        "lobes/lung_lower_lobe_right.obj",
    };
    static const char *lobe_name[5] = { "lobe_ul","lobe_ll","lobe_ur","lobe_mr","lobe_lr" };

    run_shape("sphere", make_sphere(3, 1.0, 0.02), 0.10*s);
    run_shape("lshape", make_lshape(),             0.08*s);
    run_shape("plus",   make_plus(),               0.08*s);
    for (int i = 0; i < 5; i++) {
        s_trimesh lobe = trimesh_from_obj(lobe_obj[i]);
        if (!trimesh_is_valid(&lobe)) {
            printf("%-8s could not load %s\n", lobe_name[i], lobe_obj[i]);
            continue;
        }
        s_point bmin, bmax;
        bounding_box_points(&lobe.points, &bmin, &bmax);
        double dx = bmax.x-bmin.x, dy = bmax.y-bmin.y, dz = bmax.z-bmin.z;
        run_shape(lobe_name[i], lobe, sqrt(dx*dx + dy*dy + dz*dz) / 40.0 * s);
    }
    return 0;
}
