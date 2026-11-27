/* Generates the README figures (docs/img): build a domain, call vor3d, draw.
 *
 *   ./readme_figures [out_dir [figure]]   (default out_dir: readme_figures)
 *   figure: gallery | power | cdt | orphans | medax | ashape | volsph (default: all)
 *
 * Non-convex domains are read from the .obj files in example/meshes; convex ones are built
 * from points. Drawing (shading, depth sorting, PNG/GIF) is in include/plot/vor3d_plot.h.
 * Cells are opaque without outlines: neighbouring cells never share a colour,
 * and faces between two drawn cells are skipped (they are hidden anyway).
 * Everything is deterministic (fixed RNG seeds). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include "vor3d.h"
#include "medax.h"
#include "alpha.h"
#include "volsph.h"
#include "delaunay.h"
#include "convh.h"
#include "random.h"
#include "dynarray.h"
#include "plot/vor3d_plot.h"

#ifndef MESH_DIR
#define MESH_DIR "../meshes"
#endif

static const double TOL     = 1e-12;
static const char  *OUT     = "readme_figures";
#define FIG_SIZE 2000   /* PNG side, pixels */

static s_random_context RNG;
static double rnd(void) { return random_uniform_double(&RNG); }
static double rnd01_cb(void *ctx) { return random_uniform_double(ctx); }
static s_point P3(double x, double y, double z) { return (s_point){{{x, y, z}}}; }


/* ------------------------------------------------------------------------
 * Domains, seeds, colours
 * ------------------------------------------------------------------------ */

static s_trimesh load_mesh(const char *name)
{
    char path[512];
    snprintf(path, sizeof path, "%s/%s.obj", MESH_DIR, name);
    s_trimesh T = trimesh_from_obj(path);
    if (!trimesh_is_valid(&T)) { fprintf(stderr, "cannot load %s\n", path); exit(1); }
    return T;
}

static s_vor3d_domain load_domain(const char *name)
{
    s_trimesh T = load_mesh(name);
    s_vor3d_domain d = vor3d_domain_from_trimesh(&T, &(s_vor3d_opts){ .tol_dup = TOL }, 1);
    if (!vor3d_domain_is_valid(&d)) { fprintf(stderr, "invalid domain %s\n", name); exit(1); }
    free_trimesh(&T);
    return d;
}

static s_vor3d_domain box_domain(s_point lo, s_point hi)
{
    s_point c[8];
    for (int i = 0; i < 8; i++)
        c[i] = P3(i & 1 ? hi.x : lo.x, i & 2 ? hi.y : lo.y, i & 4 ? hi.z : lo.z);
    return vor3d_domain_from_points(&(s_points){ 8, c });
}

/* Icosphere vertices on the sphere (c, r): a convex point set. */
static s_points sphere_points(s_point c, double r, int sub)
{
    return r3d_icosphere(c, r, sub, NULL, NULL);
}

static s_points random_seeds_in(const s_vor3d_domain *dom, int N)
{
    s_point lo = dom->bbox_min, hi = dom->bbox_max;
    s_points S = { 0, malloc(sizeof(s_point) * N) };
    while (S.N < N) {
        s_point p = P3(lo.x + (hi.x-lo.x)*rnd(), lo.y + (hi.y-lo.y)*rnd(), lo.z + (hi.z-lo.z)*rnd());
        if (vor3d_domain_contains(dom, p, TOL)) S.p[S.N++] = p;
    }
    return S;
}

static s_vdiagram diagram(const s_points *seeds, const double *w, const s_vor3d_domain *dom,
                          int *kept, bool merge_orphans)
{
    s_vdiagram vd = vor3d_diagram(seeds, w, dom, &(s_vor3d_opts){ .tol_dup = TOL }, kept, !dom->is_convex,
                                  merge_orphans);
    if (!vd.vcells) { fprintf(stderr, "vor3d_diagram failed\n"); exit(1); }
    printf("  diagram: %d cells, vol_rel_err %.2e\n", vd.seeds.N, vd.vol_rel_err);
    return vd;
}


/* ------------------------------------------------------------------------
 * Output
 * ------------------------------------------------------------------------ */

static void save_png(const r3d_scene *S, const char *name, r3d_view v)
{
    char path[512];
    snprintf(path, sizeof path, "%s/%s.png", OUT, name);
    r3d_png(S, path, v, FIG_SIZE);
    printf("  wrote %s\n", path);
}


/* ------------------------------------------------------------------------
 * Figures
 * ------------------------------------------------------------------------ */

/* Voronoi diagram of N_seeds random seeds in dom (takes ownership of dom). */
static void fig_diagram(const char *name, s_vor3d_domain dom, int N_seeds, r3d_view v, bool gif)
{
    printf("%s\n", name);
    s_points seeds = random_seeds_in(&dom, N_seeds);
    s_vdiagram vd = diagram(&seeds, NULL, &dom, NULL, false);
    unsigned *rgb = vor3d_plot_color_vcells(&vd);
    r3d_scene S = {0};
    vor3d_plot_add_vdiagram(&S, &vd, rgb, &dom);
    save_png(&S, name, v);
    if (gif) {
        char path[512];
        snprintf(path, sizeof path, "%s/%s.gif", OUT, name);
        r3d_spin_gif(&S, path, v, 60, 480, 20);
        printf("  wrote %s\n", path);
    }
    r3d_scene_free(&S);
    free(rgb);
    free_vdiagram(&vd);
    free_points(&seeds);
    free_vor3d_domain(&dom);
}

static void fig_gallery(void)
{
    r3d_view v = { 62, 32, 0.86 };
    fig_diagram("torus", load_domain("torus"), 110, (r3d_view){ 58, 30, 0.86 }, true);
    fig_diagram("cube", box_domain(P3(0, 0, 0), P3(1, 1, 1)), 70, v, false);
    s_points sp = sphere_points(P3(0, 0, 0), 1.0, 3);
    fig_diagram("sphere", vor3d_domain_from_points(&sp), 70, v, false);
    free_points(&sp);
    fig_diagram("cross", load_domain("cross"), 90, v, false);
    fig_diagram("trefoil", load_domain("trefoil"), 110, (r3d_view){ 20, 20, 0.86 }, false);
    fig_diagram("hollow_sphere", load_domain("bowl"), 80, (r3d_view){ 50, 30, 0.86 }, false);
}

/* Same seeds, unweighted vs power: thin slab seen from above, each seed drawn
 * with its circle of radius sqrt(w) (ignored by the unweighted diagram). A seed
 * keeps its colour in both panels; cells adjacent in either never share one. */
static void fig_power(void)
{
    printf("power\n");
    const int N = 40;
    const double H = 0.35;
    s_vor3d_domain d = box_domain(P3(0, 0, 0), P3(4, 4, H));
    s_points seeds = { N, malloc(sizeof(s_point) * N) };
    double *w = malloc(sizeof(double) * N);
    for (int i = 0; i < N; i++) {
        seeds.p[i] = P3(0.2 + 3.6*rnd(), 0.2 + 3.6*rnd(), H*(0.3 + 0.4*rnd()));
        double r = 0.55 * pow(rnd(), 1.5);
        w[i] = r*r;
    }
    s_vdiagram vd[2];
    int *kept[2], *seed_of[2];
    for (int pass = 0; pass < 2; pass++) {
        kept[pass] = malloc(sizeof(int) * N);
        vd[pass] = diagram(&seeds, pass ? w : NULL, &d, kept[pass], false);
        seed_of[pass] = malloc(sizeof(int) * vd[pass].seeds.N);
        for (int i = 0; i < N; i++) if (kept[pass][i] >= 0) seed_of[pass][kept[pass][i]] = i;
    }
    /* Greedy per-seed colouring on the union of both adjacencies. */
    int *col = malloc(sizeof(int) * N), used[VOR3D_PLOT_N_PALETTE] = {0};
    for (int i = 0; i < N; i++) {
        bool taken[VOR3D_PLOT_N_PALETTE] = {0};
        for (int pass = 0; pass < 2; pass++) {
            int ci = kept[pass][i];
            if (ci < 0) continue;
            for (int k = vd[pass].adj_offset[ci]; k < vd[pass].adj_offset[ci+1]; k++) {
                int j = seed_of[pass][vd[pass].adj_neighbor[k]];
                if (j < i) taken[col[j]] = true;
            }
        }
        int best = -1;
        for (int c = 0; c < VOR3D_PLOT_N_PALETTE; c++)
            if (!taken[c] && (best < 0 || used[c] < used[best])) best = c;
        if (best < 0) best = i % VOR3D_PLOT_N_PALETTE;
        col[i] = best;
        used[best]++;
    }
    for (int pass = 0; pass < 2; pass++) {
        unsigned *rgb = malloc(sizeof(unsigned) * vd[pass].seeds.N);
        for (int i = 0; i < N; i++) if (kept[pass][i] >= 0) rgb[kept[pass][i]] = vor3d_plot_palette(col[i]);
        r3d_scene S = { .maxlen = 0.06 };
        vor3d_plot_add_vdiagram(&S, &vd[pass], rgb, &d);
        for (int i = 0; i < N; i++) {
            double r = sqrt(w[i]);
            s_point c = P3(seeds.p[i].x, seeds.p[i].y, H + 0.03);
            if (r > 0.04) r3d_add_ring(&S, c, r - 0.011, r + 0.011, 0x262626);
            r3d_add_ring(&S, c, 0.0, 0.035, 0x262626);
        }
        save_png(&S, pass ? "power_b_laguerre" : "power_a_voronoi", (r3d_view){ 32, 20, 0.86 });
        r3d_scene_free(&S);
        free(rgb);
    }
    for (int pass = 0; pass < 2; pass++) {
        free_vdiagram(&vd[pass]);
        free(kept[pass]);
        free(seed_of[pass]);
    }
    free(col);
    free(w);
    free_points(&seeds);
    free_vor3d_domain(&d);
}

/* CDT figure: a tet is drawn if it is in the y <= 0 half; its faces shared
 * with another drawn tet are hidden. */
typedef struct { const s_tetmesh *M; int t; } s_hide_tet;

static bool tet_drawn(const s_tetmesh *M, int t)
{
    if (t < 0) return false;
    double y = 0.0;
    for (int k = 0; k < 4; k++) y += M->points.p[M->tets[4*t+k]].y;
    return y <= 0.0;
}

static bool hide_tet_face(const void *ctx, s_point fc, s_point n)
{
    (void)n;
    const s_hide_tet *h = ctx;
    const int *v = &h->M->tets[4*h->t];
    for (int k = 0; k < 4; k++) {
        if (!tet_drawn(h->M, h->M->neighbors[4*h->t+k])) continue;
        /* Shared face = the face opposite vertex k. */
        s_point g = P3(0, 0, 0);
        for (int j = 0; j < 4; j++) if (j != k) g = sum_points(g, h->M->points.p[v[j]]);
        if (distance(scale_point(g, 1.0/3.0), fc) < 1e-9) return true;
    }
    return false;
}

/* CDT of a coarse torus: the interior tets of one half, the other half's boundary as a ghost. */
static void fig_cdt(void)
{
    printf("cdt\n");
    s_trimesh T = load_mesh("torus_coarse");
    s_tetmesh M = tetrahedralize_trimesh(&T, &(s_vor3d_opts){ .tol_dup = TOL });
    r3d_scene S = {0};
    int N_drawn = 0;
    for (int t = 0; t < M.N_tets; t++) {
        if (!tet_drawn(&M, t)) continue;
        s_point v[4];
        for (int k = 0; k < 4; k++) v[k] = M.points.p[M.tets[4*t+k]];
        static const int F[12] = { 0,1,2,  0,1,3,  0,2,3,  1,2,3 };
        s_hide_tet H = { &M, t };
        /* open: no back-face culling (near-flat tets can have unreliable normals). */
        r3d_add_mesh(&S, v, 4, F, 4, (r3d_style){ .rgb = vor3d_plot_palette(((t+1) * 5) % VOR3D_PLOT_N_PALETTE),
                     .convex = true, .open = true, .hide = hide_tet_face, .hide_ctx = &H });
        N_drawn++;
    }
    int *ghost = malloc(sizeof(int) * 3 * T.Nf), Ng = 0;
    for (int f = 0; f < T.Nf; f++) {
        const int *t = T.faces + 3*f;
        if ((T.points.p[t[0]].y + T.points.p[t[1]].y + T.points.p[t[2]].y) / 3.0 <= 0.0) continue;
        memcpy(ghost + 3*Ng++, t, 3 * sizeof(int));
    }
    r3d_add_mesh(&S, T.points.p, T.points.N, ghost, Ng,
                 (r3d_style){ .rgb = VOR3D_PLOT_GHOST_RGB, .ghost = true });
    printf("  %d interior tets, %d drawn\n", M.N_tets, N_drawn);
    save_png(&S, "cdt_torus", (r3d_view){ 58, 200, 0.86 });
    r3d_scene_free(&S);
    free(ghost);
    free_trimesh(&T);
    free_tetmesh(&M);
}

/* Volume of the orphans of a cell (all surface components but the largest). */
static double orphan_volume(const s_vcell *c)
{
    if (c->N_surface < 2) return 0.0;
    double vmax = 0.0, sum = 0.0;
    for (int k = 0; k < c->N_surface; k++) {
        double v = fabs(volume_trimesh(&c->surface[k]));
        sum += v;
        vmax = fmax(vmax, v);
    }
    return sum - vmax;
}

/* Orphan merging on a U shape: raw (a cell split across the gap) vs merged.
 * Seeds: the largest orphan over 150 random tries. The owner of the orphan is
 * red, the cell that receives it is blue, the rest pale tints. */
static void fig_orphans(void)
{
    printf("orphans\n");
    s_vor3d_domain d = load_domain("ushape");
    s_points best = {0};
    double best_v = 0.0;
    for (int attempt = 0; attempt < 150; attempt++) {
        s_points seeds = random_seeds_in(&d, 8);
        s_vdiagram raw = vor3d_diagram(&seeds, NULL, &d, &(s_vor3d_opts){ .tol_dup = TOL }, NULL, true, false);
        double v = 0.0;
        for (int i = 0; raw.vcells && i < raw.seeds.N; i++) v = fmax(v, orphan_volume(&raw.vcells[i]));
        free_vdiagram(&raw);
        if (v > best_v) { best_v = v; free_points(&best); best = seeds; }
        else free_points(&seeds);
    }
    printf("  largest orphan volume %.3f (domain %.3f)\n", best_v, d.volume);
    s_vdiagram raw = diagram(&best, NULL, &d, NULL, false);
    s_vdiagram mg  = diagram(&best, NULL, &d, NULL, true);
    int owner = 0;
    for (int i = 0; i < raw.seeds.N; i++)
        if (orphan_volume(&raw.vcells[i]) > orphan_volume(&raw.vcells[owner])) owner = i;
    for (int pass = 0; pass < 2; pass++) {
        const s_vdiagram *vd = pass ? &mg : &raw;
        unsigned *col = vor3d_plot_color_vcells(vd);
        for (int i = 0; i < vd->seeds.N; i++) {
            bool receiver = mg.vcells[i].volume - raw.vcells[i].volume > 1e-9;
            if (i == owner) col[i] = 0xE15759;
            else if (receiver) col[i] = 0x4E79A7;
            else {   /* 60% towards light grey */
                unsigned pale = 0;
                for (int ch = 0; ch < 3; ch++)
                    pale |= (unsigned)(((col[i] >> (8*ch)) & 0xFF) * (1 - 0.6) + 0xE8 * 0.6) << (8*ch);
                col[i] = pale;
            }
        }
        r3d_scene S = {0};
        vor3d_plot_add_vdiagram(&S, vd, col, &d);
        save_png(&S, pass ? "orphans_b_merged" : "orphans_a_raw", (r3d_view){ 40, 20, 0.86 });
        r3d_scene_free(&S);
        free(col);
    }
    free_vdiagram(&raw);
    free_vdiagram(&mg);
    free_points(&best);
    free_vor3d_domain(&d);
}

static unsigned viridis(double t)
{
    static const double s[5][3] = {
        {0.267, 0.005, 0.329}, {0.229, 0.322, 0.546}, {0.128, 0.567, 0.551},
        {0.369, 0.789, 0.383}, {0.993, 0.906, 0.144},
    };
    t = fmin(fmax(t, 0.0), 1.0) * 4.0;
    int i = (int)t; if (i > 3) i = 3;
    double f = t - i, c[3];
    for (int k = 0; k < 3; k++) c[k] = s[i][k]*(1-f) + s[i+1][k]*f;
    return ((unsigned)(c[0]*255) << 16) | ((unsigned)(c[1]*255) << 8) | (unsigned)(c[2]*255);
}

/* Medial graph of an L: a sphere per node, a tube per edge, coloured by medial-ball radius. */
static void fig_medax(void)
{
    printf("medax\n");
    s_trimesh T = load_mesh("lshape");
    s_random_context r = random_initialize(7);
    s_medax ma = medax_from_trimesh(&T, 0.25, 0, 0, &(s_vor3d_opts){ .tol_dup = TOL, .randd01 = rnd01_cb, .rctx = &r }, true);
    if (!ma.radius || !ma.graph) { fprintf(stderr, "medax failed\n"); exit(1); }
    const s_medial_graph *G = ma.graph;
    printf("  %d medial points, %d graph edges\n", ma.verts.N, G->adj_off[G->N] / 2);
    double rmax = 0;
    for (int i = 0; i < ma.verts.N; i++) rmax = fmax(rmax, ma.radius[i]);
    r3d_scene S = {0};
    vor3d_plot_add_trimesh(&S, &T, (r3d_style){ .rgb = VOR3D_PLOT_GHOST_RGB, .ghost = true });
    for (int i = 0; i < G->N; i++) {
        r3d_add_sphere(&S, ma.verts.p[i], 0.024, viridis(ma.radius[i] / rmax));
        for (int k = G->adj_off[i]; k < G->adj_off[i+1]; k++) {
            int j = G->adj[k];
            if (j > i)
                r3d_add_tube(&S, ma.verts.p[i], ma.verts.p[j], 0.008,
                             viridis(0.5 * (ma.radius[i] + ma.radius[j]) / rmax));
        }
    }
    save_png(&S, "medial_axis", (r3d_view){ 44, 35, 0.86 });
    r3d_scene_free(&S);
    free_medax(&ma);
    free_trimesh(&T);
}

/* Alpha shape: random points inside the trefoil, then the reconstructed surface. */
static void fig_ashape(void)
{
    printf("ashape\n");
    s_vor3d_domain d = load_domain("trefoil");
    const int N = 6000;
    s_points pts = random_seeds_in(&d, N);
    double vol = d.volume, h = cbrt(vol / N);
    s_random_context r = random_initialize(11);
    s_alpha_shape as = alpha_shape_3d(&pts, NULL, pow(2.4*h, 2), &(s_vor3d_opts){ .tol_dup = 1e-12, .randd01 = rnd01_cb, .rctx = &r }, NULL, false);
    if (!alpha_shape_is_valid(&as)) { fprintf(stderr, "alpha shape failed\n"); exit(1); }
    printf("  alpha %.4g, %d faces, volume %.4f (domain %.4f)\n", as.alpha, as.surface.Nf, as.volume, vol);
    r3d_view v = { 20, 20, 0.86 };
    r3d_scene S = {0};
    for (int i = 0; i < N; i++) r3d_add_sphere(&S, pts.p[i], 0.032, vor3d_plot_palette(0));
    save_png(&S, "ashape_a_points", v);
    r3d_scene_free(&S);
    S = (r3d_scene){0};
    vor3d_plot_add_trimesh(&S, &as.surface, (r3d_style){ .rgb = vor3d_plot_palette(3) });
    save_png(&S, "ashape_b_surface", v);
    r3d_scene_free(&S);
    free_alpha_shape(&as);
    free_points(&pts);
    free_vor3d_domain(&d);
}

/* Union of balls split into per-ball contributions (ball INTERSECT power cell). */
static void fig_volsph(void)
{
    printf("volsph\n");
    const int N = 28;
    s_points c = { N, malloc(sizeof(s_point) * N) };
    double *rad = malloc(sizeof(double) * N), *w = malloc(sizeof(double) * N);
    /* Each ball after the first overlaps an earlier one, so none is isolated. */
    for (int i = 0; i < N; i++) {
        rad[i] = 0.35 + 0.35*rnd();
        w[i] = rad[i]*rad[i];
        bool overlaps;
        do {
            s_point p;
            do p = P3(2*rnd()-1, 2*rnd()-1, 2*rnd()-1); while (norm(p) > 1);
            c.p[i] = scale_point(p, 1.3);
            overlaps = i == 0;
            for (int j = 0; j < i && !overlaps; j++)
                overlaps = distance(c.p[i], c.p[j]) < rad[i] + rad[j];
        } while (!overlaps);
    }
    s_vor3d_domain d = box_domain(P3(-2.5, -2.5, -2.5), P3(2.5, 2.5, 2.5));
    int *kept = malloc(sizeof(int) * N);
    s_vdiagram vd = diagram(&c, w, &d, kept, false);

    double *contrib = malloc(sizeof(double) * N), Vcontrib = 0;
    double Vu = volume_union_spheres(&c, rad, &(s_vor3d_opts){ .tol_dup = 1e-12 }, contrib);
    for (int i = 0; i < N; i++) Vcontrib += contrib[i];

    /* Polyhedral piece of each ball: icosphere hull INTERSECT its power cell. */
    s_convh *piece = calloc(N, sizeof(s_convh));
    double Vpieces = 0;
    for (int i = 0; i < N; i++) {
        if (kept[i] < 0 || vd.vcells[kept[i]].N_pieces == 0) continue;
        s_points bp = sphere_points(c.p[i], rad[i], 3);
        s_convh ball;
        if (convhull_from_points(&bp, &ball) != 1) { free_points(&bp); continue; }
        int ret = intersection_convhulls(&ball, &vd.vcells[kept[i]].pieces[0], TOL, &piece[i]);
        if (ret != 1) {
            printf("  ball %d: intersection_convhulls returned %d (power cell volume %.4f)\n",
                   i, ret, vd.vcells[kept[i]].volume);
            piece[i] = convhull_NAN;
        }
        else
            Vpieces += volume_convhull(&piece[i]);
        free_convhull(&ball);
        free_points(&bp);
    }
    printf("  union volume %.5f, sum of contributions %.5f, polyhedral pieces %.5f\n",
           Vu, Vcontrib, Vpieces);

    unsigned *rgb = vor3d_plot_color_vcells(&vd);
    for (int pass = 0; pass < 2; pass++) {
        r3d_scene S = {0};
        for (int i = 0; i < N; i++) {
            if (!convhull_is_valid(&piece[i])) continue;
            s_point cm = convhull_volume_centroid(&piece[i]);
            s_point shift = pass ? scale_point(cm, 0.55) : P3(0, 0, 0);
            vor3d_plot_add_convh(&S, &piece[i], (r3d_style){ .shift = shift, .rgb = rgb[kept[i]] });
        }
        save_png(&S, pass ? "volsph_b_exploded" : "volsph_a_union", (r3d_view){ 62, 32, 0.86 });
        r3d_scene_free(&S);
    }
    for (int i = 0; i < N; i++) if (convhull_is_valid(&piece[i])) free_convhull(&piece[i]);
    free(piece);
    free(rgb);
    free(contrib);
    free_vdiagram(&vd);
    free(kept);
    free_vor3d_domain(&d);
    free(w); free(rad);
    free_points(&c);
}


int main(int argc, char **argv)
{
    if (argc > 1) OUT = argv[1];
    mkdir(OUT, 0755);
    const char *only = argc > 2 ? argv[2] : NULL;
#define RUN(name, fn) do { if (!only || !strcmp(only, name)) { RNG = random_initialize(2026); fn(); } } while (0)
    RUN("gallery", fig_gallery);
    RUN("power",   fig_power);
    RUN("cdt",     fig_cdt);
    RUN("orphans", fig_orphans);
    RUN("medax",   fig_medax);
    RUN("ashape",  fig_ashape);
    RUN("volsph",  fig_volsph);
    return 0;
}
