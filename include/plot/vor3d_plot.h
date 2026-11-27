/* vor3d_plot.h -- optional header-only plotting for vor3d: PNG or spinning GIF via gnuplot.
 *
 * Not part of libvor3d: nothing is compiled unless you include it.
 *     #include "plot/vor3d_plot.h"
 * Needs gnuplot (pngcairo) in PATH at run time; GIFs also need ffmpeg.
 *
 * API
 *   One call:
 *     vor3d_plot_vdiagram(vd, dom, path, view, size)   PNG of every cell; neighbours differ in colour
 *     vor3d_plot_vcell(vd, i, dom, path, view, size)   PNG of cell i inside the transparent domain
 *   Scenes:
 *     r3d_scene S = {0};                               empty scene
 *     r3d_add_mesh(&S, p, Np, tri, Nf, style)          triangles tri[3*f + k] into p[]
 *     r3d_add_sphere(&S, c, r, rgb)                    small lit sphere
 *     r3d_add_tube(&S, a, b, r, rgb)                   thin tube from a to b
 *     r3d_add_ring(&S, c, r0, r1, rgb)                 flat horizontal annulus (r0 = 0: disc)
 *     vor3d_plot_add_vcell(&S, vd, i, style)
 *     vor3d_plot_add_vdiagram(&S, vd, rgb, dom)        every cell, rgb[i] per cell
 *     vor3d_plot_add_trimesh(&S, mesh, style)          e.g. dom.surface, alpha shape as.surface
 *     vor3d_plot_add_convh(&S, convh, style)
 *     vor3d_plot_add_points(&S, pts, r, rgb)           a sphere of radius r per point
 *     vor3d_plot_add_tetmesh(&S, tm, shrink, style)    every tet, shrunk by shrink (e.g. 0.8)
 *     vor3d_plot_add_medax(&S, ma, r, rgb)             medial points (+ graph edges if built)
 *     vor3d_plot_add_alpha_complex(&S, ac, r, rgb)     outer tet faces, dangling faces/edges/vertices
 *     r3d_png(&S, path, view, size)                    square PNG, size x size pixels
 *     r3d_spin_gif(&S, path, view, N_frames, size, fps)
 *     r3d_scene_free(&S)
 *   Colours:
 *     vor3d_plot_color_vcells(vd)                      one colour per cell, neighbours never equal
 *     vor3d_plot_palette(i)                            colour i of the 12-colour cell palette
 *     VOR3D_PLOT_GHOST_RGB                             grey used for transparent context surfaces
 *   Helper:
 *     r3d_icosphere(c, r, sub, &tri, &Nf)              subdivided icosahedron on a sphere
 *
 * Arguments
 *   path   output file (.png, or .gif for r3d_spin_gif)
 *   view   (r3d_view){ rot_x, rot_z, fill }: camera angles in degrees, as gnuplot's "set view",
 *          and the fraction of the image the scene spans (0 = default 0.86, which fits it
 *          with a margin; larger zooms in, smaller zooms out)
 *   size   image side in pixels (images are square)
 *   rgb    a colour 0xRRGGBB (unsigned), e.g. 0x80B1D3
 *   r      radius of the drawn spheres (tubes are r/3 thick), in scene units
 *   style  (r3d_style){ .rgb = ... } plus optional fields, see r3d_style
 *
 * Colouring cells yourself: vor3d_plot_color_vcells returns a malloc'd rgb[vd->seeds.N]
 * (free it); change any entry, then pass it to vor3d_plot_add_vdiagram:
 *     unsigned *rgb = vor3d_plot_color_vcells(&vd);
 *     rgb[7] = 0xE15759;                               highlight cell 7 in red
 *     r3d_scene S = {0};
 *     vor3d_plot_add_vdiagram(&S, &vd, rgb, &dom);
 *     r3d_png(&S, "vd.png", (r3d_view){ 60, 30, 0.86 }, 1000);
 *     r3d_scene_free(&S);
 *     free(rgb);
 *
 * Output is deterministic: the same scene gives the same file. */
#ifndef VOR3D_PLOT_H
#define VOR3D_PLOT_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gplot.h"
#include "vor3d.h"
#include "delaunay.h"
#include "alpha.h"
#include "medax.h"

typedef struct r3d_view {
    double rot_x, rot_z;  /* gnuplot "set view rot_x, rot_z" (degrees) */
    double fill;          /* fraction of the image the scene spans; 0 = 0.86 (> 1 zooms in) */
} r3d_view;

typedef struct r3d_style {
    unsigned rgb;
    bool     ghost;       /* transparent context surface instead of opaque */
    s_point  shift;       /* translation (exploded views) */
    bool     convex;      /* orient normals away from the mesh centre (for unreliable winding) */
    bool     open;        /* not a closed outward surface: never back-face culled */
    /* Optional: faces with hide(hide_ctx, centroid, unit normal) true are skipped. */
    bool   (*hide)(const void *ctx, s_point centroid, s_point normal);
    const void *hide_ctx;
} r3d_style;

typedef struct r3d__polys {   /* internal: polygon soup */
    int nv, capv, np, capp;
    double *x, *y, *z;
    int *off;                 /* np+1 */
    unsigned *rgb;            /* np, unlit base colour */
    s_point *nrm;             /* np, unit normal (lit polygons) */
    unsigned char *flg;       /* np, R3D__LIT | R3D__CULL | R3D__SPLIT */
} r3d__polys;

typedef struct r3d_scene {
    double maxlen;            /* faces are split until edges are shorter; 0 = 0.04 x scene diameter */
    r3d__polys faces, ghost;  /* internal */
    bool has_bbox;            /* internal */
    s_point bmin, bmax;       /* internal */
} r3d_scene;


/* ------------------------------------------------------------------------
 * Internals
 * ------------------------------------------------------------------------ */

#define R3D__PI 3.14159265358979323846
enum { R3D__LIT = 1, R3D__CULL = 2, R3D__SPLIT = 4 };

static inline s_point r3d__p(double x, double y, double z) { return (s_point){{{ x, y, z }}}; }
static inline s_point r3d__unit(s_point a) { double n = norm(a); return n > 0 ? scale_point(a, 1.0/n) : a; }

static inline void r3d__vert(r3d__polys *P, s_point p)
{
    if (P->nv == P->capv) {
        P->capv = P->capv ? 2*P->capv : 1024;
        P->x = realloc(P->x, sizeof(double) * P->capv);
        P->y = realloc(P->y, sizeof(double) * P->capv);
        P->z = realloc(P->z, sizeof(double) * P->capv);
    }
    P->x[P->nv] = p.x; P->y[P->nv] = p.y; P->z[P->nv] = p.z;
    P->nv++;
}

static inline s_point r3d__vert_at(const r3d__polys *P, int k) { return r3d__p(P->x[k], P->y[k], P->z[k]); }

static inline void r3d__close(r3d__polys *P, unsigned rgb, s_point n, int flags)
{
    if (P->np + 1 >= P->capp) {
        P->capp = P->capp ? 2*P->capp : 1024;
        P->off = realloc(P->off, sizeof(int) * (P->capp + 1));
        P->rgb = realloc(P->rgb, sizeof(unsigned) * P->capp);
        P->nrm = realloc(P->nrm, sizeof(s_point) * P->capp);
        P->flg = realloc(P->flg, P->capp);
    }
    if (P->np == 0) P->off[0] = 0;
    P->rgb[P->np] = rgb;
    P->nrm[P->np] = n;
    P->flg[P->np] = flags;
    P->off[++P->np] = P->nv;
}

static inline void r3d__free(r3d__polys *P)
{
    free(P->x); free(P->y); free(P->z); free(P->off); free(P->rgb); free(P->nrm); free(P->flg);
    *P = (r3d__polys){0};
}

static inline void r3d__bbox(r3d_scene *S, s_point p)
{
    if (!S->has_bbox) { S->bmin = S->bmax = p; S->has_bbox = true; return; }
    for (int k = 0; k < 3; k++) {
        if (p.coords[k] < S->bmin.coords[k]) S->bmin.coords[k] = p.coords[k];
        if (p.coords[k] > S->bmax.coords[k]) S->bmax.coords[k] = p.coords[k];
    }
}

/* Emit a triangle, split 1->4 until every edge is below maxlen. */
static inline void r3d__emit_tri(r3d__polys *dst, s_point a, s_point b, s_point c, s_point n,
                                 unsigned rgb, int flags, double maxlen)
{
    double l = fmax(distance(a, b), fmax(distance(b, c), distance(c, a)));
    if (l > maxlen) {
        s_point ab = scale_point(sum_points(a, b), 0.5), bc = scale_point(sum_points(b, c), 0.5),
                ca = scale_point(sum_points(c, a), 0.5);
        r3d__emit_tri(dst, a, ab, ca, n, rgb, flags, maxlen);
        r3d__emit_tri(dst, ab, b, bc, n, rgb, flags, maxlen);
        r3d__emit_tri(dst, ca, bc, c, n, rgb, flags, maxlen);
        r3d__emit_tri(dst, ab, bc, ca, n, rgb, flags, maxlen);
        return;
    }
    r3d__vert(dst, a); r3d__vert(dst, b); r3d__vert(dst, c);
    r3d__close(dst, rgb, n, flags);
}

/* Copy of the opaque faces with the R3D__SPLIT triangles subdivided: gnuplot sorts
 * whole polygons by one depth (painter's algorithm), so big faces must be split. */
static inline r3d__polys r3d__split_faces(const r3d_scene *S)
{
    double maxlen = S->maxlen > 0 ? S->maxlen : 0.04 * distance(S->bmin, S->bmax);
    const r3d__polys *P = &S->faces;
    r3d__polys Q = {0};
    for (int i = 0; i < P->np; i++) {
        int o = P->off[i];
        if (P->flg[i] & R3D__SPLIT)
            r3d__emit_tri(&Q, r3d__vert_at(P, o), r3d__vert_at(P, o+1), r3d__vert_at(P, o+2), P->nrm[i],
                          P->rgb[i], P->flg[i] & ~R3D__SPLIT, maxlen);
        else {
            for (int k = o; k < P->off[i+1]; k++) r3d__vert(&Q, r3d__vert_at(P, k));
            r3d__close(&Q, P->rgb[i], P->nrm[i], P->flg[i]);
        }
    }
    return Q;
}

/* gnuplot's orthographic camera: direction towards the viewer, screen right / up. */
static inline void r3d__view_basis(double rot_x, double rot_z, s_point *dir, s_point *right, s_point *up)
{
    double rx = rot_x * R3D__PI / 180.0, rz = rot_z * R3D__PI / 180.0;
    *dir = r3d__p(sin(rx)*sin(rz), -sin(rx)*cos(rz), cos(rx));
    s_point z = r3d__p(0, 0, 1);
    *up = r3d__unit(subtract_points(z, scale_point(*dir, dot_prod(z, *dir))));
    if (norm(*up) < 0.5) *up = r3d__p(0, 1, 0);
    *right = cross_prod(*up, *dir);
}

/* Camera-relative key light (upper left), ambient, soft highlight; two-sided. */
static inline void r3d__shade(const r3d__polys *P, double rot_x, double rot_z, unsigned *out)
{
    s_point v, r, u;
    r3d__view_basis(rot_x, rot_z, &v, &r, &u);
    s_point l = r3d__unit(sum_points(sum_points(v, scale_point(u, 0.6)), scale_point(r, -0.45)));
    s_point h = r3d__unit(sum_points(l, v));
    for (int i = 0; i < P->np; i++) {
        unsigned c = P->rgb[i];
        if (!(P->flg[i] & R3D__LIT)) { out[i] = c; continue; }
        s_point n = P->nrm[i];
        if (dot_prod(n, v) < 0) n = scale_point(n, -1.0);
        double diff = fmax(0.0, dot_prod(n, l));
        double spec = 0.05 * pow(fmax(0.0, dot_prod(n, h)), 40.0);
        double k = 0.55 + 0.50 * diff;
        unsigned rgb = 0;
        for (int ch = 2; ch >= 0; ch--) {
            double x = ((c >> (8*ch)) & 0xFF) * k + 255.0 * spec;
            rgb |= (unsigned)fmin(255.0, fmax(0.0, x)) << (8*ch);
        }
        out[i] = rgb;
    }
}

static inline s_point r3d__bbox_centre(const r3d_scene *S)
{
    return scale_point(sum_points(S->bmin, S->bmax), 0.5);
}

/* Half-extent of the scene on screen, and the 3D point that projects to its centre. */
static inline double r3d__screen_extent(const r3d_scene *S, const r3d__polys *F, double rot_x, double rot_z,
                                        s_point *centre)
{
    s_point v, r, u;
    r3d__view_basis(rot_x, rot_z, &v, &r, &u);
    double umin = 1e300, umax = -1e300, vmin = 1e300, vmax = -1e300;
    const r3d__polys *L[2] = { F, &S->ghost };
    for (int k = 0; k < 2; k++)
        for (int i = 0; i < L[k]->nv; i++) {
            s_point p = r3d__vert_at(L[k], i);
            double a = dot_prod(p, r), b = dot_prod(p, u);
            umin = fmin(umin, a); umax = fmax(umax, a);
            vmin = fmin(vmin, b); vmax = fmax(vmax, b);
        }
    s_point c0 = r3d__bbox_centre(S);
    double du = 0.5*(umin + umax) - dot_prod(c0, r), dv = 0.5*(vmin + vmax) - dot_prod(c0, u);
    if (centre) *centre = sum_points(c0, sum_points(scale_point(r, du), scale_point(u, dv)));
    return 0.5 * fmax(umax - umin, vmax - vmin);
}

/* gnuplot clips to the ranges, so they are the bbox around c (half-width R);
 * the view scale then zooms so the on-screen extent E fills the requested
 * fraction of the canvas (one unit of range half-width spans ~0.252 of it). */
static inline double r3d__range_halfwidth(const r3d_scene *S, s_point c, double E)
{
    double R = E;
    for (int k = 0; k < 3; k++)
        R = fmax(R, fmax(S->bmax.coords[k] - c.coords[k], c.coords[k] - S->bmin.coords[k]));
    return R * 1.001;
}

static inline double r3d__view_scale(r3d_view v, double E, double R)
{
    double fill = v.fill > 0 ? v.fill : 0.86;
    return (fill / 2.0) / 0.252 * (R / E);
}

static inline void r3d__config(s_gplot *g, s_point c, double E, double R, r3d_view v)
{
    char cmd[512];
    gplot_config(g,
        "unset border", "unset tics", "unset key", "unset colorbox", "unset title",
        "set margins 0,0,0,0",
        "set view equal xyz",
        "set xyplane relative 0",
        "set pm3d depthorder",
        "set style fill solid 1.0 noborder");
    snprintf(cmd, sizeof cmd, "set xrange [%.17g:%.17g]", c.x - R, c.x + R); gplot_config(g, cmd);
    snprintf(cmd, sizeof cmd, "set yrange [%.17g:%.17g]", c.y - R, c.y + R); gplot_config(g, cmd);
    snprintf(cmd, sizeof cmd, "set zrange [%.17g:%.17g]", c.z - R, c.z + R); gplot_config(g, cmd);
    snprintf(cmd, sizeof cmd, "set view %g,%g,%g", v.rot_x, v.rot_z, r3d__view_scale(v, E, R));
    gplot_config(g, cmd);
}

/* F = the split opaque faces. Closed meshes are back-face culled: this also avoids
 * painter's-algorithm errors between big front faces and small faces behind them. */
static inline void r3d__draw(s_gplot *g, const r3d_scene *S, const r3d__polys *F, double rot_x, double rot_z)
{
    const r3d__polys *G = &S->ghost;
    if (G->np > 0)
        gplot_draw_polygons_3d(g, G->np, G->off, G->x, G->y, G->z, G->rgb,
                               .extra = "fs transparent solid 0.14 noborder");
    if (F->np > 0) {
        s_point v, r, u;
        r3d__view_basis(rot_x, rot_z, &v, &r, &u);
        unsigned *rgb = malloc(sizeof(unsigned) * F->np);
        r3d__shade(F, rot_x, rot_z, rgb);
        r3d__polys Q = {0};
        for (int i = 0; i < F->np; i++) {
            if ((F->flg[i] & R3D__CULL) && dot_prod(F->nrm[i], v) < 0.0) continue;
            for (int k = F->off[i]; k < F->off[i+1]; k++) r3d__vert(&Q, r3d__vert_at(F, k));
            r3d__close(&Q, rgb[i], F->nrm[i], 0);
        }
        gplot_draw_polygons_3d(g, Q.np, Q.off, Q.x, Q.y, Q.z, Q.rgb, .extra = NULL);
        r3d__free(&Q);
        free(rgb);
    }
}


/* ------------------------------------------------------------------------
 * Scenes
 * ------------------------------------------------------------------------ */

static inline void r3d_scene_free(r3d_scene *S)
{
    r3d__free(&S->faces);
    r3d__free(&S->ghost);
    *S = (r3d_scene){0};
}

/* Triangles tri[3*f + k] (indices into p[0..Np-1]). */
static inline void r3d_add_mesh(r3d_scene *S, const s_point *p, int Np, const int *tri, int Nf, r3d_style st)
{
    if (Nf <= 0) return;
    s_point *q = malloc(sizeof(s_point) * Np);
    for (int i = 0; i < Np; i++) q[i] = sum_points(p[i], st.shift);
    s_point cm = r3d__p(0, 0, 0);
    if (st.convex) {
        for (int f = 0; f < Nf; f++)
            for (int k = 0; k < 3; k++) cm = sum_points(cm, q[tri[3*f+k]]);
        cm = scale_point(cm, 1.0 / (3.0 * Nf));
    }
    r3d__polys *dst = st.ghost ? &S->ghost : &S->faces;
    int flags = st.ghost ? 0 : (R3D__LIT | R3D__SPLIT | (st.open ? 0 : R3D__CULL));
    for (int f = 0; f < Nf; f++) {
        s_point a = q[tri[3*f]], b = q[tri[3*f+1]], c = q[tri[3*f+2]];
        s_point n = r3d__unit(cross_prod(subtract_points(b, a), subtract_points(c, a)));
        s_point fc = scale_point(sum_points(sum_points(a, b), c), 1.0/3.0);
        if (st.convex && dot_prod(n, subtract_points(fc, cm)) < 0) n = scale_point(n, -1.0);
        if (st.hide && st.hide(st.hide_ctx, fc, n)) continue;
        r3d__vert(dst, a); r3d__vert(dst, b); r3d__vert(dst, c);
        r3d__close(dst, st.rgb, n, flags);
        r3d__bbox(S, a); r3d__bbox(S, b); r3d__bbox(S, c);
    }
    free(q);
}

/* Icosahedron subdivided sub times, projected to the sphere (c, r); shared vertices
 * merged, triangles wound outward. tri (malloc'd tri[3*Nf]) may be NULL. */
static inline s_points r3d_icosphere(s_point c, double r, int sub, int **tri, int *Nf)
{
    static const double t = 1.6180339887498949;
    static const double V[12][3] = {
        {-1, t, 0}, {1, t, 0}, {-1, -t, 0}, {1, -t, 0}, {0, -1, t}, {0, 1, t},
        {0, -1, -t}, {0, 1, -t}, {t, 0, -1}, {t, 0, 1}, {-t, 0, -1}, {-t, 0, 1},
    };
    static const int F[20][3] = {
        {0,11,5},{0,5,1},{0,1,7},{0,7,10},{0,10,11},{1,5,9},{5,11,4},{11,10,2},{10,7,6},{7,1,8},
        {3,9,4},{3,4,2},{3,2,6},{3,6,8},{3,8,9},{4,9,5},{2,4,11},{6,2,10},{8,6,7},{9,8,1},
    };
    int N_tri = 20;
    s_point *T = malloc(sizeof(s_point) * 3 * N_tri);
    for (int f = 0; f < 20; f++)
        for (int k = 0; k < 3; k++) T[3*f+k] = r3d__unit(r3d__p(V[F[f][k]][0], V[F[f][k]][1], V[F[f][k]][2]));
    for (int s = 0; s < sub; s++) {
        s_point *nt = malloc(sizeof(s_point) * 12 * N_tri);
        for (int f = 0; f < N_tri; f++) {
            s_point a = T[3*f], b = T[3*f+1], cc = T[3*f+2];
            s_point ab = r3d__unit(sum_points(a, b)), bc = r3d__unit(sum_points(b, cc)),
                    ca = r3d__unit(sum_points(cc, a));
            s_point q[12] = { a, ab, ca,  ab, b, bc,  ca, bc, cc,  ab, bc, ca };
            memcpy(nt + 12*f, q, sizeof q);
        }
        free(T);
        T = nt;
        N_tri *= 4;
    }
    s_point *P = malloc(sizeof(s_point) * 3 * N_tri);
    int *I = malloc(sizeof(int) * 3 * N_tri), n = 0;
    for (int f = 0; f < 3*N_tri; f++) {
        s_point p = sum_points(c, scale_point(T[f], r));
        int i = 0;
        while (i < n && fabs(P[i].x - p.x) + fabs(P[i].y - p.y) + fabs(P[i].z - p.z) >= 1e-12) i++;
        if (i == n) P[n++] = p;
        I[f] = i;
    }
    free(T);
    if (tri) { *tri = I; *Nf = N_tri; }
    else free(I);
    return (s_points){ n, realloc(P, sizeof(s_point) * n) };
}

/* Small lit sphere (icosphere, 1 subdivision = 80 faces) with radial normals. */
static inline void r3d_add_sphere(r3d_scene *S, s_point c, double r, unsigned rgb)
{
    static s_points U;
    static int *F, NF;
    if (!U.p) U = r3d_icosphere(r3d__p(0, 0, 0), 1.0, 1, &F, &NF);
    for (int f = 0; f < NF; f++) {
        s_point n = r3d__p(0, 0, 0);
        for (int k = 0; k < 3; k++) {
            s_point u = U.p[F[3*f+k]];
            n = sum_points(n, u);
            r3d__vert(&S->faces, sum_points(c, scale_point(u, r)));
        }
        r3d__close(&S->faces, rgb, r3d__unit(n), R3D__LIT | R3D__CULL);
    }
    r3d__bbox(S, sum_points(c, r3d__p(r, r, r)));
    r3d__bbox(S, subtract_points(c, r3d__p(r, r, r)));
}

/* Thin hexagonal tube from a to b (no caps). */
static inline void r3d_add_tube(r3d_scene *S, s_point a, s_point b, double r, unsigned rgb)
{
    s_point d = subtract_points(b, a);
    if (norm(d) < 1e-12) return;
    s_point t = r3d__unit(d);
    s_point e1 = r3d__unit(cross_prod(t, fabs(t.x) < 0.9 ? r3d__p(1, 0, 0) : r3d__p(0, 1, 0)));
    s_point e2 = cross_prod(t, e1);
    for (int k = 0; k < 6; k++) {
        double t0 = 2*R3D__PI*k/6, t1 = 2*R3D__PI*(k+1)/6, tm = 0.5*(t0 + t1);
        s_point u0 = sum_points(scale_point(e1, cos(t0)), scale_point(e2, sin(t0)));
        s_point u1 = sum_points(scale_point(e1, cos(t1)), scale_point(e2, sin(t1)));
        r3d__vert(&S->faces, sum_points(a, scale_point(u0, r)));
        r3d__vert(&S->faces, sum_points(b, scale_point(u0, r)));
        r3d__vert(&S->faces, sum_points(b, scale_point(u1, r)));
        r3d__vert(&S->faces, sum_points(a, scale_point(u1, r)));
        r3d__close(&S->faces, rgb, sum_points(scale_point(e1, cos(tm)), scale_point(e2, sin(tm))),
                   R3D__LIT | R3D__CULL);
    }
}

/* Flat unlit horizontal annulus [r0, r1] around c (r0 = 0: disc). */
static inline void r3d_add_ring(r3d_scene *S, s_point c, double r0, double r1, unsigned rgb)
{
    const int K = 72;
    for (int k = 0; k < K; k++) {
        double a0 = 2*R3D__PI*k/K, a1 = 2*R3D__PI*(k+1)/K;
        r3d__vert(&S->faces, sum_points(c, r3d__p(r1*cos(a0), r1*sin(a0), 0)));
        r3d__vert(&S->faces, sum_points(c, r3d__p(r1*cos(a1), r1*sin(a1), 0)));
        if (r0 > 0) r3d__vert(&S->faces, sum_points(c, r3d__p(r0*cos(a1), r0*sin(a1), 0)));
        r3d__vert(&S->faces, sum_points(c, r3d__p(r0*cos(a0), r0*sin(a0), 0)));
        r3d__close(&S->faces, rgb, r3d__p(0, 0, 1), 0);
    }
    r3d__bbox(S, sum_points(c, r3d__p(r1, r1, 0)));
    r3d__bbox(S, subtract_points(c, r3d__p(r1, r1, 0)));
}


/* ------------------------------------------------------------------------
 * Output
 * ------------------------------------------------------------------------ */

/* Square PNG, size x size pixels. */
static inline void r3d_png(const r3d_scene *S, const char *path, r3d_view v, int size)
{
    r3d__polys F = r3d__split_faces(S);
    s_point c;
    double E = r3d__screen_extent(S, &F, v.rot_x, v.rot_z, &c);
    s_gplot *g = gplot_start(GPLOT_PNG_3D, path, .size = {size, size});
    r3d__config(g, c, E, r3d__range_halfwidth(S, c, E), v);
    r3d__draw(g, S, &F, v.rot_x, v.rot_z);
    gplot_end(g);
    r3d__free(&F);
}

/* GIF of a full turn about the vertical axis through the scene centre (needs ffmpeg). */
static inline void r3d_spin_gif(const r3d_scene *S, const char *path, r3d_view v,
                                int N_frames, int size, int fps)
{
    char cmd[128];
    r3d__polys F = r3d__split_faces(S);
    s_point c = r3d__bbox_centre(S);
    double E = 0;
    for (int f = 0; f < 24; f++) {
        s_point cf;
        double e = r3d__screen_extent(S, &F, v.rot_x, v.rot_z + 15.0*f, &cf);
        E = fmax(E, e + norm(r3d__p(cf.x - c.x, cf.y - c.y, 0)));
    }
    long cores = sysconf(_SC_NPROCESSORS_ONLN);   /* one gnuplot process per core renders frames */
    int workers = cores < 1 ? 1 : cores > N_frames ? N_frames : (int)cores;
    s_gplot *g = gplot_start(GPLOT_GIF_3D, path, .size = {size, size}, .framerate = fps,
                             .video_parallel_threads = workers);
    double R = r3d__range_halfwidth(S, c, E);
    r3d__config(g, c, E, R, v);
    for (int f = 0; f < N_frames; f++) {
        double rz = fmod(v.rot_z + 360.0 * f / N_frames, 360.0);
        snprintf(cmd, sizeof cmd, "set view %g,%g,%g", v.rot_x, rz, r3d__view_scale(v, E, R));
        if (f == 0) gplot_config(g, cmd);
        else        gplot_next_frame(g, cmd);
        r3d__draw(g, S, &F, v.rot_x, rz);
    }
    gplot_end(g);
    r3d__free(&F);
}


/* ------------------------------------------------------------------------
 * vor3d objects
 * ------------------------------------------------------------------------ */

#define VOR3D_PLOT_N_PALETTE 12
#define VOR3D_PLOT_GHOST_RGB 0x9AA5B1u

/* Colour i of the cell palette (matplotlib Set3, reordered so that 0 = blue, 3 = orange). */
static inline unsigned vor3d_plot_palette(int i)
{
    static const unsigned P[VOR3D_PLOT_N_PALETTE] = {
        0x80B1D3, 0xB3DE69, 0xFB8072, 0xFDB462, 0xBEBADA, 0x8DD3C7,
        0xFCCDE5, 0xFFED6F, 0xBC80BD, 0xCCEBC5, 0xFFFFB3, 0xD9D9D9,
    };
    return P[((i % VOR3D_PLOT_N_PALETTE) + VOR3D_PLOT_N_PALETTE) % VOR3D_PLOT_N_PALETTE];
}

/* malloc'd rgb[vd->seeds.N]: greedy colouring on the cell adjacency, no two neighbours
 * share a colour (when the palette allows); among the free colours the least used. */
static inline unsigned *vor3d_plot_color_vcells(const s_vdiagram *vd)
{
    int N = vd->seeds.N, used[VOR3D_PLOT_N_PALETTE] = {0};
    int *col = malloc(sizeof(int) * (N > 0 ? N : 1));
    unsigned *rgb = malloc(sizeof(unsigned) * (N > 0 ? N : 1));
    for (int i = 0; i < N; i++) {
        bool taken[VOR3D_PLOT_N_PALETTE] = {0};
        for (int k = vd->adj_offset[i]; k < vd->adj_offset[i+1]; k++)
            if (vd->adj_neighbor[k] < i) taken[col[vd->adj_neighbor[k]]] = true;
        int best = -1;
        for (int c = 0; c < VOR3D_PLOT_N_PALETTE; c++) {
            int cc = (c + i) % VOR3D_PLOT_N_PALETTE;
            if (!taken[cc] && (best < 0 || used[cc] < used[best])) best = cc;
        }
        if (best < 0) best = i % VOR3D_PLOT_N_PALETTE;
        col[i] = best;
        used[best]++;
        rgb[i] = vor3d_plot_palette(best);
    }
    free(col);
    return rgb;
}

/* Normals are oriented outward whatever the face winding (st.convex forced on). */
static inline void vor3d_plot_add_convh(r3d_scene *S, const s_convh *C, r3d_style st)
{
    st.convex = true;
    r3d_add_mesh(S, C->points.p, C->points.N, C->faces, C->Nf, st);
}

static inline void vor3d_plot_add_trimesh(r3d_scene *S, const s_trimesh *M, r3d_style st)
{
    r3d_add_mesh(S, M->points.p, M->points.N, M->faces, M->Nf, st);
}

static inline void vor3d_plot_add_points(r3d_scene *S, const s_points *pts, double r, unsigned rgb)
{
    for (int i = 0; i < pts->N; i++) r3d_add_sphere(S, pts->p[i], r, rgb);
}

/* Internal: tet v[0..3] scaled about its centroid by shrink. */
static inline void vor3d_plot__add_tet(r3d_scene *S, const s_point v[4], double shrink, r3d_style st)
{
    static const int F[12] = { 0,1,2,  0,1,3,  0,2,3,  1,2,3 };
    s_point g = scale_point(sum_points(sum_points(v[0], v[1]), sum_points(v[2], v[3])), 0.25), w[4];
    for (int k = 0; k < 4; k++) w[k] = sum_points(g, scale_point(subtract_points(v[k], g), shrink));
    st.convex = true;   /* tet winding is not guaranteed */
    r3d_add_mesh(S, w, 4, F, 4, st);
}

static inline void vor3d_plot_add_tetmesh(r3d_scene *S, const s_tetmesh *tm, double shrink, r3d_style st)
{
    for (int t = 0; t < tm->N_tets; t++) {
        s_point v[4];
        for (int k = 0; k < 4; k++) v[k] = tm->points.p[tm->tets[4*t+k]];
        vor3d_plot__add_tet(S, v, shrink, st);
    }
}

/* A sphere per medial point; a tube per medial-graph edge if the graph was built. */
static inline void vor3d_plot_add_medax(r3d_scene *S, const s_medax *ma, double r, unsigned rgb)
{
    vor3d_plot_add_points(S, &ma->verts, r, rgb);
    const s_medial_graph *G = ma->graph;
    for (int i = 0; G && i < G->N; i++)
        for (int k = G->adj_off[i]; k < G->adj_off[i+1]; k++)
            if (G->adj[k] > i) r3d_add_tube(S, ma->verts.p[i], ma->verts.p[G->adj[k]], r / 3, rgb);
}

/* Internal: index of the ascending k-simplex key in the lexicographically sorted list, or -1. */
static inline int vor3d_plot__find(const int *list, int N, int k, const int *key)
{
    int lo = 0, hi = N - 1;
    while (lo <= hi) {
        int mid = (lo + hi) / 2, c = 0;
        for (int j = 0; j < k && !c; j++) c = (list[k*mid+j] > key[j]) - (list[k*mid+j] < key[j]);
        if (!c) return mid;
        if (c < 0) lo = mid + 1; else hi = mid - 1;
    }
    return -1;
}

/* The faces of K on the boundary of its tets, faces of no tet, edges of no face (tubes)
 * and vertices of no edge (spheres of radius r). Uses the sorted lists of s_alpha_complex. */
static inline void vor3d_plot_add_alpha_complex(r3d_scene *S, const s_alpha_complex *ac, double r, unsigned rgb)
{
    int *n_tets = calloc(ac->N_faces + 1, sizeof(int));   /* tets per face */
    bool *e_used = calloc(ac->N_edges + 1, sizeof(bool));  /* edge of some face */
    bool *v_used = calloc(ac->points.N + 1, sizeof(bool)); /* vertex of some edge */
    for (int t = 0; t < ac->N_tets; t++)
        for (int skip = 0; skip < 4; skip++) {
            int key[3], n = 0;
            for (int j = 0; j < 4; j++) if (j != skip) key[n++] = ac->tets[4*t+j];
            int f = vor3d_plot__find(ac->faces, ac->N_faces, 3, key);
            if (f >= 0) n_tets[f]++;
        }
    int *tri = malloc(sizeof(int) * 3 * (ac->N_faces + 1)), Nt = 0;
    for (int f = 0; f < ac->N_faces; f++) {
        const int *v = &ac->faces[3*f];
        for (int a = 0; a < 3; a++)
            for (int b = a + 1; b < 3; b++) {
                int e = vor3d_plot__find(ac->edges, ac->N_edges, 2, (int[2]){ v[a], v[b] });
                if (e >= 0) e_used[e] = true;
            }
        if (n_tets[f] < 2) { memcpy(tri + 3*Nt, v, 3 * sizeof(int)); Nt++; }
    }
    r3d_add_mesh(S, ac->points.p, ac->points.N, tri, Nt, (r3d_style){ .rgb = rgb, .open = true });
    for (int e = 0; e < ac->N_edges; e++) {
        const int *v = &ac->edges[2*e];
        v_used[v[0]] = v_used[v[1]] = true;
        if (!e_used[e]) r3d_add_tube(S, ac->points.p[v[0]], ac->points.p[v[1]], r / 3, rgb);
    }
    for (int i = 0; i < ac->points.N; i++)
        if (!v_used[i]) r3d_add_sphere(S, ac->points.p[i], r, rgb);
    free(tri); free(n_tets); free(e_used); free(v_used);
}

/* Cell i: its surface components if built, else its convex pieces. */
static inline void vor3d_plot_add_vcell(r3d_scene *S, const s_vdiagram *vd, int i, r3d_style st)
{
    const s_vcell *c = &vd->vcells[i];
    if (!vd->is_convex && c->N_surface > 0)
        for (int k = 0; k < c->N_surface; k++) vor3d_plot_add_trimesh(S, &c->surface[k], st);
    else
        for (int p = 0; p < c->N_pieces; p++) vor3d_plot_add_convh(S, &c->pieces[p], st);
}

/* Internal: a cell face is hidden if a point just outside it is still in the domain
 * (then another cell is drawn on its other side). */
typedef struct vor3d_plot__hide { const s_vor3d_domain *dom; double eps; } vor3d_plot__hide;

static inline bool vor3d_plot__hide_inner_face(const void *ctx, s_point fc, s_point n)
{
    const vor3d_plot__hide *h = ctx;
    return vor3d_domain_contains(h->dom, sum_points(fc, scale_point(n, h->eps)), 1e-12);
}

/* Every cell, rgb[i] per cell; faces between two cells are skipped (always hidden). */
static inline void vor3d_plot_add_vdiagram(r3d_scene *S, const s_vdiagram *vd, const unsigned *rgb,
                                          const s_vor3d_domain *dom)
{
    vor3d_plot__hide H = { dom, 1e-4 * distance(dom->bbox_min, dom->bbox_max) };
    for (int i = 0; i < vd->seeds.N; i++)
        vor3d_plot_add_vcell(S, vd, i, (r3d_style){ .rgb = rgb[i], .hide = vor3d_plot__hide_inner_face,
                                                    .hide_ctx = &H });
}

static inline void vor3d_plot_vdiagram(const s_vdiagram *vd, const s_vor3d_domain *dom,
                                      const char *path, r3d_view v, int size)
{
    unsigned *rgb = vor3d_plot_color_vcells(vd);
    r3d_scene S = {0};
    vor3d_plot_add_vdiagram(&S, vd, rgb, dom);
    r3d_png(&S, path, v, size);
    r3d_scene_free(&S);
    free(rgb);
}

static inline void vor3d_plot_vcell(const s_vdiagram *vd, int i, const s_vor3d_domain *dom,
                                   const char *path, r3d_view v, int size)
{
    r3d_scene S = {0};
    vor3d_plot_add_vcell(&S, vd, i, (r3d_style){ .rgb = vor3d_plot_palette(0) });
    vor3d_plot_add_trimesh(&S, &dom->surface, (r3d_style){ .rgb = VOR3D_PLOT_GHOST_RGB, .ghost = true });
    r3d_png(&S, path, v, size);
    r3d_scene_free(&S);
}

#endif
