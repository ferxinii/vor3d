/* Medial-axis center scaling study.
 * For each test shape (sphere, L-shape, plus, and the 5 lung lobes) the medial
 * axis is extracted at THREE sampling densities -- ~0.1x, 1x and ~10x the point
 * count (r_sample scaled by sqrt(10)) -- and for each we time medax_center (the
 * warm-started descent) and validate it against an OpenMP-parallel brute-force
 * all-pairs exact center. A summary table reports how the center search scales
 * with graph size and whether it finds the true minimum at every size.
 *
 * Each run also writes <shape>_surface.vtk and <shape>_<scale>_axis.vtk (the
 * medial axis at each of the three granularities) for viewing in ParaView.
 *
 * Run from the example/ build dir:  ./example_medax [base_scale] [novalidate]
 *   base_scale (default 1.0) multiplies every density's r_sample.
 *   novalidate (or 0) as a 2nd arg skips the O(N^2) exact check -- a fast run
 *     that still writes the three axis-VTK granularities.
 *
 * The L-shape / plus builders are copied from example_nonconvex.c. */

/* ---- Reference results (measured 2026-07, 10 OpenMP threads) --------------
 * medax_center = forced warm-started descent; exact = OpenMP brute-force
 * all-pairs argmin. The descent found the TRUE global minimum in ALL 24 cases
 * (same node, E_ratio = 1.000000). Times are wall-clock seconds.
 *
 *   shape    scale     nodes   extract    center     exact
 *   sphere   coarse      184     0.12s     0.003s     0.00s
 *   sphere   as-is       799     0.17s     0.050s     0.01s
 *   sphere   fine      26435     2.37s     3.151s     8.70s
 *   lshape   coarse     1412     0.03s     0.053s     0.02s
 *   lshape   as-is     15633     0.53s     0.902s     2.44s
 *   lshape   fine     162208    14.18s    13.922s   321.55s
 *   plus     coarse     2147     0.04s     0.077s     0.04s
 *   plus     as-is     23678     0.82s     1.507s     5.77s
 *   plus     fine     249884    22.16s    21.002s   786.19s
 *   lobe_ul  coarse    15278     1.79s     0.637s     2.22s
 *   lobe_ul  as-is     15203     2.12s     0.628s     2.22s
 *   lobe_ul  fine      79725     4.65s     4.658s    68.03s
 *   lobe_ll  coarse    16639     2.19s     1.003s     2.69s
 *   lobe_ll  as-is     16673     2.66s     1.244s     2.74s
 *   lobe_ll  fine      88239     4.93s    11.071s    86.19s
 *   lobe_ur  coarse    14179     1.61s     0.563s     1.92s
 *   lobe_ur  as-is     14065     1.78s     0.645s     1.90s
 *   lobe_ur  fine      76059     4.50s     7.150s    65.60s
 *   lobe_mr  coarse    15467     1.65s     0.908s     2.34s
 *   lobe_mr  as-is     15376     1.88s     0.788s     2.30s
 *   lobe_mr  fine     107649     6.15s    10.073s   139.07s
 *   lobe_lr  coarse    16309     2.11s     0.874s     2.61s
 *   lobe_lr  as-is     16160     2.39s     0.874s     2.57s
 *   lobe_lr  fine      81795     4.76s     6.804s    75.08s
 *
 * Takeaways:
 *  - Correctness: exact global minimum at every size, 184 -> 249884 nodes.
 *  - The center search scales ~N^1.2 (near-linear); brute force is O(N^2), so
 *    the descent's advantage grows with size: ~3.8x at 24k nodes, ~37x at 250k.
 *  - Crossover ~2-3k nodes: below that, exact is actually faster than the
 *    descent (multi-start + polish overhead) -- hence medax_center's exact_below
 *    default (~2500). See "plus coarse" above (descent 0.077s vs exact 0.04s).
 *  - Lobe coarse ~= as-is (~15k): the lobe meshes are already at the 1-sample-
 *    per-face floor, so r_sample cannot coarsen them further.
 *  - "sphere fine" is the per-node outlier: a near-degenerate tight cluster
 *    (nearly equal radii -> very flat energy) makes the confirmation work
 *    hardest, so its center time is high for its node count. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "trimesh.h"
#include "points.h"
#include "medax.h"
#include "volsph.h"
#include "random.h"

static double EPS_DEG = 1e-9;
static double TOL     = 1e-9;

static s_random_context rctx;
static double randd01(void *ctx) { return random_uniform_double(ctx); }

/* ---------------- L-shape (from example_nonconvex.c) --------------------- */
/* NOTE: kept axis-aligned on purpose. The medial-axis inside/outside test now
 * uses the exact CDT point-location oracle (not +x ray-casting), so axis-
 * aligned meshes work directly -- no rotation jitter needed. */
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
    return trimesh_from_arrays(verts, 12, faces, 20, EPS_DEG);
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
    return trimesh_from_arrays(verts, 24, faces, 44, EPS_DEG);
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

    s_trimesh m = trimesh_from_arrays(V, nV, F, nF, EPS_DEG);
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

/* ----- OpenMP brute-force ground truth (parallel all-pairs argmin) --------- */
/* The library keeps its exact solver static; for validation we reimplement the
 * O(N^2) all-pairs center here and parallelize the outer (per-source) loop with
 * OpenMP -- each thread runs its own Dijkstra, then we reduce to the argmin. */
typedef struct { double d; int v; } ex_hn;
static void ex_hup(ex_hn *h, int i)
{ while (i>0){ int p=(i-1)/2; if(h[p].d<=h[i].d)break; ex_hn t=h[p];h[p]=h[i];h[i]=t;i=p; } }
static void ex_hdown(ex_hn *h, int n, int i)
{ for(;;){ int l=2*i+1,r=l+1,m=i; if(l<n&&h[l].d<h[m].d)m=l; if(r<n&&h[r].d<h[m].d)m=r;
           if(m==i)break; ex_hn t=h[m];h[m]=h[i];h[i]=t;i=m; } }
static void ex_dijkstra(const s_medial_graph *g, int src, double *dist, ex_hn *h)
{
    int N=g->N;
    for(int i=0;i<N;i++) dist[i]=INFINITY;
    int hn=0; dist[src]=0.0; h[hn++]=(ex_hn){0.0,src};
    while(hn>0){ ex_hn top=h[0]; h[0]=h[--hn]; ex_hdown(h,hn,0);
        if(top.d>dist[top.v]) continue;
        for(int k=g->adj_off[top.v];k<g->adj_off[top.v+1];k++){
            int w=g->adj[k]; double nd=top.d+g->adj_w[k];
            if(nd<dist[w]){ dist[w]=nd; h[hn++]=(ex_hn){nd,w}; ex_hup(h,hn-1); } } }
}
/* connected components via BFS; returns largest component id, sets *largest */
static int ex_components(const s_medial_graph *g, int *comp, int *queue, int *largest)
{
    int N=g->N; for(int i=0;i<N;i++) comp[i]=-1;
    int nc=0,best=0,root=0;
    for(int s=0;s<N;s++){ if(comp[s]!=-1)continue; int qh=0,qt=0,sz=0; comp[s]=nc; queue[qt++]=s;
        while(qh<qt){ int u=queue[qh++]; sz++;
            for(int k=g->adj_off[u];k<g->adj_off[u+1];k++)
                if(comp[g->adj[k]]==-1){ comp[g->adj[k]]=nc; queue[qt++]=g->adj[k]; } }
        if(sz>best){best=sz;root=s;} nc++; }
    *largest=best; return comp[root];
}
/* parallel exact center over component lc: argmin_i sum_j w_j d(i,j)^2, where
 * w_j = mass[j] is the same governed-volume weight the library uses (so this
 * brute force validates the library's descent on the SAME energy). */
static int parallel_exact(const s_medial_graph *g, const double *mass,
                          const int *comp, int lc, double *out_E)
{
    int N=g->N, gbest=-1; double gE=INFINITY;
    #pragma omp parallel
    {
        double *dist=malloc(sizeof(double)*(size_t)N);
        ex_hn  *heap=malloc(sizeof(ex_hn)*(size_t)(g->adj_off[N]+8));
        int lbest=-1; double lE=INFINITY;
        if(dist&&heap){
            #pragma omp for schedule(dynamic,32)
            for(int i=0;i<N;i++){
                if(comp[i]!=lc) continue;
                ex_dijkstra(g,i,dist,heap);
                double e=0.0;
                for(int j=0;j<N;j++){ if(comp[j]!=lc||dist[j]==INFINITY)continue;
                    e+=mass[j]*dist[j]*dist[j]; }
                if(e<lE){ lE=e; lbest=i; }
            }
        }
        #pragma omp critical
        { if(lE<gE || (lE==gE && lbest>=0 && (gbest<0||lbest<gbest))){ gE=lE; gbest=lbest; } }
        free(dist); free(heap);
    }
    *out_E=gE; return gbest;
}

/* ---------------- scaling study: one measurement row -------------------- */
typedef struct {
    const char *shape;
    const char *scale;
    int    N, nE;
    double r_sample;
    double t_extract;   /* medial-axis extraction (wall) */
    double t_center;    /* medax_center, forced descent (wall) */
    double t_exact;     /* OpenMP brute-force all-pairs (wall) */
    int    same;        /* descent node == brute-force exact node? (-1 = not run) */
    double ratio;       /* E_descent / E_exact */
} s_row;

/* Time medax_center (forced descent) on the medial graph; if validate, also run
 * the OpenMP-parallel brute-force exact center and compare. Fills row's fields. */
static void measure(const s_medax *ma, s_row *row, int validate)
{
    const s_medial_graph *g = ma->graph;
    int N = g->N; graph_stats(g, &row->nE); row->N = N;

    double Ed, t1 = omp_get_wtime();
    int cd = medax_center(ma, 0, 3, 2.0, -1, &Ed);     /* forced descent (exact_below=0) */
    row->t_center = omp_get_wtime() - t1;

    if (validate) {
        int *comp = malloc(sizeof(int)*(size_t)N), *queue = malloc(sizeof(int)*(size_t)N);
        /* same governed-volume weight the library uses (see cen_weights): the
         * per-ball union contribution vol(B_j INTERSECT L_j). NULL buffer ok. */
        double *mass = malloc(sizeof(double)*(size_t)N);
        if (comp && queue && mass && cd >= 0) {
            volume_contribution_spheres(&ma->verts, ma->radius, 1e-10, 1e-9, NULL, mass);
            /* Use the RAW governed volume (no clamp): now that volsph emits no
             * negative per-ball contributions, these are exactly the weights the
             * library's cen_weights uses, so ratio == 1 is a genuine accuracy
             * check (a leak would resurface as ratio != 1, not be masked). */
            int largest, lc = ex_components(g, comp, queue, &largest);
            double Ex, t0 = omp_get_wtime();
            int cx = parallel_exact(g, mass, comp, lc, &Ex);
            row->t_exact = omp_get_wtime() - t0;
            row->same    = (cx == cd);
            row->ratio   = (Ex > 0) ? Ed / Ex : 1.0;
        }
        free(comp); free(queue); free(mass);
    }
    printf("  [%-6s %-8s] N=%-7d extract=%7.2fs center=%7.3fs exact=%8.2fs  "
           "same=%s ratio=%.6f\n", row->scale, row->shape, row->N, row->t_extract,
           row->t_center, row->t_exact, row->same==1?"yes":(row->same==0?"NO":"-"),
           row->ratio);
}

/* Extract the medial axis of one mesh at spacing r_sample, write its surface and
 * (per-scale) axis VTK for ParaView, then measure the center. Frees the mesh. */
static void run_shape(const char *shape, const char *scale, s_trimesh mesh,
                      double r_sample, s_row *row, int validate)
{
    *row = (s_row){ .shape = shape, .scale = scale, .r_sample = r_sample, .same = -1 };
    if (!trimesh_is_valid(&mesh)) { printf("  [%s %s] INVALID mesh\n", scale, shape); return; }

    char path[256];
    snprintf(path, sizeof path, "%s_surface.vtk", shape);
    write_surface_vtk(&mesh, path);

    double t = omp_get_wtime();
    s_medax ma = medax_from_trimesh(&mesh, r_sample, MEDAX_THETA0_DEFAULT,
                                    MEDAX_RHO0_DEFAULT, randd01, &rctx, TOL, EPS_DEG,
                                    /*build_medial_graph=*/true);
    row->t_extract = omp_get_wtime() - t;
    if (!ma.verts.p || !ma.graph) {
        printf("  [%s %s] extraction FAILED\n", scale, shape);
        free_trimesh(&mesh); return;
    }
    snprintf(path, sizeof path, "%s_%s_axis.vtk", shape, scale);   /* per-granularity */
    medax_write_vtk(&ma, path);

    measure(&ma, row, validate);
    free_medax(&ma);
    free_trimesh(&mesh);
}

int main(int argc, char **argv)
{
    rctx = random_initialize(12345);
    double base = (argc > 1) ? atof(argv[1]) : 1.0;
    /* pass "novalidate" (or 0) as a 2nd arg to skip the O(N^2) exact check --
     * a fast run that still writes the 3 granularity levels of axis VTK. */
    int validate = !(argc > 2 && (strcmp(argv[2], "novalidate") == 0 ||
                                  strcmp(argv[2], "0") == 0));

    /* Three sampling densities. Sample count ~ 1/r_sample^2, so scaling r_sample
     * by sqrt(10) changes the point (and graph) count by ~10x either way. */
    const double SMULT[3] = { 3.16227766, 1.0, 0.31622777 };  /* r_sample multiplier */
    const char  *SNAME[3] = { "coarse", "as-is", "fine" };    /* ~0.1x, 1x, ~10x points */

    static const char *lobe_obj[5] = {
        "lobes/lung_upper_lobe_left.obj",  "lobes/lung_lower_lobe_left.obj",
        "lobes/lung_upper_lobe_right.obj", "lobes/lung_middle_lobe_right.obj",
        "lobes/lung_lower_lobe_right.obj",
    };
    static const char *shape_name[8] =
        { "sphere","lshape","plus","lobe_ul","lobe_ll","lobe_ur","lobe_mr","lobe_lr" };

    s_row rows[3][8];   /* [scale][shape] */

    for (int sc = 0; sc < 3; sc++) {
        double s = base * SMULT[sc];
        printf("\n########## scale = %s (r_sample x %.3g, ~%s points) ##########\n",
               SNAME[sc], SMULT[sc], sc==0 ? "0.1x" : (sc==1 ? "1x" : "10x"));
        run_shape("sphere", SNAME[sc], make_sphere(3, 1.0, 0.02), 0.10*s, &rows[sc][0], validate);
        run_shape("lshape", SNAME[sc], make_lshape(),             0.08*s, &rows[sc][1], validate);
        run_shape("plus",   SNAME[sc], make_plus(),               0.08*s, &rows[sc][2], validate);
        for (int i = 0; i < 5; i++) {
            s_trimesh lobe = trimesh_from_obj(lobe_obj[i], EPS_DEG);
            if (!trimesh_is_valid(&lobe)) {
                printf("  [%s %s] could not load %s\n", SNAME[sc], shape_name[3+i], lobe_obj[i]);
                rows[sc][3+i] = (s_row){ .shape = shape_name[3+i], .scale = SNAME[sc], .same = -1 };
                continue;
            }
            s_point bmin, bmax;
            bounding_box_points(&lobe.points, &bmin, &bmax);
            double dx = bmax.x-bmin.x, dy = bmax.y-bmin.y, dz = bmax.z-bmin.z;
            double diag = sqrt(dx*dx + dy*dy + dz*dz);
            run_shape(shape_name[3+i], SNAME[sc], lobe, (diag/40.0)*s, &rows[sc][3+i], validate);
        }
    }

    /* summary table, grouped by shape so the 3 densities line up */
    printf("\n============================= SUMMARY =============================\n");
    printf("%-8s %-7s %9s %9s %9s %10s %5s %10s\n",
           "shape", "scale", "nodes", "extr(s)", "cntr(s)", "exact(s)", "same", "E_ratio");
    int all_exact = 1;
    for (int sh = 0; sh < 8; sh++) {
        for (int sc = 0; sc < 3; sc++) {
            s_row *r = &rows[sc][sh];
            printf("%-8s %-7s %9d %9.2f %9.3f %10.2f %5s %10.6f\n",
                   r->shape, r->scale, r->N, r->t_extract, r->t_center, r->t_exact,
                   r->same==1 ? "yes" : (r->same==0 ? "NO" : "-"), r->ratio);
            if (r->same == 0) all_exact = 0;
        }
    }
    printf("------------------------------------------------------------------\n");
    printf("descent found the true minimum in every measured case: %s\n",
           all_exact ? "YES" : "NO");
    return 0;
}
