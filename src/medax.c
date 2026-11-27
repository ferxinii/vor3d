/* Approximate medial axis extraction (Dey-Zhao MEDIAL).
 * See MEDIAL_AXIS_PLAN.md for the full derivation and pipeline.
 *
 * Step 1 (skeleton): public API compiles; free_medax / medax_write_vtk done.
 * Step 2: surface sampling (mesh vertices ++ area-proportional stratified face
 * scatter), per-sample virtual-umbrella records {normals, R}, and the Delaunay
 * triangulation of the samples, with the Nreal == Ns guard.
 * Step 3: enumerate unique Delaunay edges and filter them by the angle/ratio
 * conditions + inner test; the selected inner edges are stored in a hash set.
 * Step 4: emit one circumcenter + medial radius per DT tet that touches a
 * selected edge (its Voronoi vertex), with an on-by-default inside safety
 * filter.
 * Step 5: build the Voronoi 1-skeleton graph on those points -- an edge per
 * pair of face-adjacent emitted tets (dual Voronoi edge), stored as CSR
 * adjacency in the s_medax. Every edge is a certified-inside corridor.
 *
 * Inside/outside (Steps 3-4) is decided by point location in the flagged
 * winding-number test against the surface itself (formerly a convex-hull CDT;
 * replaced because CDT construction can fail on degenerate cut cells), a
 * oracle -- no ray-casting degeneracy, so axis-aligned meshes need no jitter. */

#include "medax.h"

#include "delaunay.h"
#include "scplx.h"
#include "./hlibs/dynarray.h"
#include "./hlibs/hash.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#define MEDAX_MARGIN 1e-6              /* barycentric interior margin (Step 2) */
static const double MEDAX_SQRT3 = 1.7320508075688772;

/* Dey-Zhao filter defaults (paper: theta0 = pi/8, rho0 = 8). The paper tunes
 * rho0 against umbrella-triangle circumradii ~ spacing/sqrt(3); our R is the
 * spacing itself, so the effective ratio threshold is rho0 / sqrt(3). */
#define MEDAX_THETA0_DEFAULT (3.14159265358979323846 / 8.0)
#define MEDAX_RHO0_DEFAULT   8.0
#define MEDAX_RHO_EFF(rho0)  ((rho0) / MEDAX_SQRT3)

/* ----- internal surface-sample bundle (Step 2 output) --------------------- */
/* samples: mesh vertices [0..Nv) ++ scattered face points [Nv..N); the scatter
 * range is Morton-sorted at generation time (see build_samples).
 * un_off/un: CSR of per-sample UNIT umbrella normals.
 * R_uniform: local size (== r_sample) for the uniform-scatter case;
 * R: per-sample local size, non-NULL only in the r_sample<=0 (vertices-only)
 *    fallback where the spacing is not uniform. */
typedef struct medax_samples {
    s_points  samples;
    int       Nv;
    int      *un_off;    /* size samples.N + 1 */
    s_point  *un;        /* size un_off[samples.N] */
    double    R_uniform;
    double   *R;         /* size samples.N, or NULL */
} s_medax_samples;

typedef struct scatpt { s_point p; int f; } s_scatpt;

static void free_samples(s_medax_samples *s)
{
    if (!s) return;
    free(s->samples.p);
    free(s->un_off);
    free(s->un);
    free(s->R);
    *s = (s_medax_samples){0};
}

/* ----- small geometry helpers -------------------------------------------- */
static inline s_point vsub(s_point a, s_point b)
{ return (s_point){{{ a.x - b.x, a.y - b.y, a.z - b.z }}}; }

static inline s_point vcross(s_point a, s_point b)
{
    return (s_point){{{ a.y * b.z - a.z * b.y,
                        a.z * b.x - a.x * b.z,
                        a.x * b.y - a.y * b.x }}};
}

static inline double vnorm(s_point a)
{ return sqrt(a.x * a.x + a.y * a.y + a.z * a.z); }

/* face corners (returns pointers into mesh->points) */
static inline void face_corners(const s_trimesh *m, int f,
                                s_point *A, s_point *B, s_point *C)
{
    *A = m->points.p[m->faces[3 * f + 0]];
    *B = m->points.p[m->faces[3 * f + 1]];
    *C = m->points.p[m->faces[3 * f + 2]];
}

static double face_area(const s_trimesh *m, int f)
{
    s_point A, B, C; face_corners(m, f, &A, &B, &C);
    return 0.5 * vnorm(vcross(vsub(B, A), vsub(C, A)));
}

/* unit outward normal from mesh->fnormals (guaranteed outward-oriented, but
 * stored unnormalized). Returns 0 for a degenerate (near-zero) normal. */
static int face_unit_normal(const s_trimesh *m, int f, double eps, s_point *out)
{
    s_point n = m->fnormals[f];
    double L = vnorm(n);
    if (L <= eps) return 0;
    *out = (s_point){{{ n.x / L, n.y / L, n.z / L }}};
    return 1;
}

/* triangle circumradius = a*b*c / (4*Area); used only for the vertices-only
 * (r_sample<=0) local-size fallback. Returns 0 for a degenerate face. */
static double face_circumradius(const s_trimesh *m, int f)
{
    s_point A, B, C; face_corners(m, f, &A, &B, &C);
    double a = vnorm(vsub(B, C)), b = vnorm(vsub(C, A)), c = vnorm(vsub(A, B));
    double area = face_area(m, f);
    if (area <= 0.0) return 0.0;
    return (a * b * c) / (4.0 * area);
}

/* ----- the surface scatter (Step 0.5) ------------------------------------ */
/* Stratified jittered sampling of one face: split into k^2 sub-triangles on a
 * barycentric lattice, one uniform jittered point per sub-triangle, each
 * clamped to a tiny interior margin so it cannot collide (within TOL) with a
 * mesh vertex or a scattered point on the adjacent face across a shared edge.
 * (alpha, beta) are the face-affine coords: p = A + alpha*(B-A) + beta*(C-A).
 * Pushes {p, f} into `out`. Returns 0 on allocation failure. */
static int scatter_face(const s_trimesh *m, int f, double r_sample,
                        double (*randd01)(void*), void *rctx, s_dynarray *out)
{
    double area = face_area(m, f);
    if (area <= 0.0) return 1;                 /* degenerate face: no scatter */

    double a_s = (MEDAX_SQRT3 / 4.0) * r_sample * r_sample;
    int k = (int)lround(sqrt(area / a_s));
    if (k < 1) k = 1;

    s_point A, B, C; face_corners(m, f, &A, &B, &C);
    s_point AB = vsub(B, A), AC = vsub(C, A);

    for (int i = 0; i < k; i++) {
        for (int j = 0; j < k - i; j++) {
            /* up-triangle (i,j) always exists; down-triangle when i+j<=k-2 */
            for (int down = 0; down <= 1; down++) {
                if (down && (i + j > k - 2)) continue;

                double r1 = randd01(rctx), r2 = randd01(rctx);
                double su = sqrt(r1);
                double w0 = 1.0 - su, w1 = su * (1.0 - r2), w2 = su * r2;

                double alpha, beta;
                if (!down) {                    /* corners (i,j)(i+1,j)(i,j+1) */
                    alpha = (i + w1) / k;
                    beta  = (j + w2) / k;
                } else {              /* corners (i+1,j)(i,j+1)(i+1,j+1) */
                    alpha = (i + 1.0 - w1) / k;
                    beta  = (j + 1.0 - w0) / k;
                }

                /* interior-margin clamp in barycentric coords, renormalized */
                double u = 1.0 - alpha - beta, v = alpha, w = beta;
                if (u < MEDAX_MARGIN || v < MEDAX_MARGIN || w < MEDAX_MARGIN) {
                    if (u < MEDAX_MARGIN) u = MEDAX_MARGIN;
                    if (v < MEDAX_MARGIN) v = MEDAX_MARGIN;
                    if (w < MEDAX_MARGIN) w = MEDAX_MARGIN;
                    double s = u + v + w; u /= s; v /= s; w /= s;
                    alpha = v; beta = w;
                }

                s_point p = {{{ A.x + alpha * AB.x + beta * AC.x,
                                A.y + alpha * AB.y + beta * AC.y,
                                A.z + alpha * AB.z + beta * AC.z }}};
                s_scatpt sp = { p, f };
                if (!dynarray_push(out, &sp)) return 0;
            }
        }
    }
    return 1;
}

/* ----- assemble samples + umbrella records + local sizes ----------------- */
static s_medax_samples build_samples(const s_trimesh *mesh, double r_sample,
                                     double (*randd01)(void*), void *rctx,
                                     double EPS_DEG)
{
    s_medax_samples out = {0};
    int Nv = mesh->points.N;
    int Nf = mesh->Nf;

    /* scatter phase (only when r_sample > 0) */
    s_dynarray scat = dynarray_initialize(sizeof(s_scatpt), 64);
    if (r_sample > 0.0) {
        if (!scat.items) return out;
        for (int f = 0; f < Nf; f++)
            if (!scatter_face(mesh, f, r_sample, randd01, rctx, &scat)) {
                dynarray_free(&scat);
                return out;
            }
    }
    int Nscat = (int)scat.N;
    int N = Nv + Nscat;

    /* Morton-sort the scatter round BEFORE anything derives from its order, so
     * the whole bundle (coords, umbrella CSR) is assembled in sorted order by
     * construction -- consecutive DT insertions are then spatially adjacent
     * and the walk hint locates them in O(1) (MEDAX_SPEEDUP_PLAN.md Fix 2).
     * Only the scatter round is sorted: the mesh vertices [0..Nv) stay first
     * as a coarse scaffold of the whole surface (BRIO rounds).  Sorting ALL
     * samples globally measurably REGRESSED the sample DT (0.36s -> 0.65s on
     * the lobe): clustered insertion leaves giant tets over unvisited regions
     * and later clusters pay large flip cascades there.  Deterministic (pure
     * function of coordinates); on allocation failure the scatter keeps its
     * face-by-face order, which is valid (just less local). */
    if (Nscat > 1) {
        s_scatpt *sp   = (s_scatpt *)scat.items;
        int      *perm = malloc(sizeof(int) * (size_t)Nscat);
        s_point  *crd  = malloc(sizeof(s_point) * (size_t)Nscat);
        s_scatpt *srt  = malloc(sizeof(s_scatpt) * (size_t)Nscat);
        if (perm && crd && srt) {
            for (int j = 0; j < Nscat; j++) crd[j] = sp[j].p;
            s_points cloud = { .N = Nscat, .p = crd };
            if (points_morton_permutation(&cloud, perm)) {
                for (int j = 0; j < Nscat; j++) srt[j] = sp[perm[j]];
                memcpy(sp, srt, sizeof(s_scatpt) * (size_t)Nscat);
            }
        }
        free(perm); free(crd); free(srt);
    }

    /* per-vertex incident-face count (valid faces only), for the umbrella CSR */
    int *deg = calloc((size_t)Nv, sizeof(int));
    out.un_off = malloc(sizeof(int) * (size_t)(N + 1));
    out.samples.p = malloc(sizeof(s_point) * (size_t)N);
    if (!deg || !out.un_off || !out.samples.p) goto fail;

    for (int f = 0; f < Nf; f++) {
        s_point nf;
        if (!face_unit_normal(mesh, f, EPS_DEG, &nf)) continue;
        for (int c = 0; c < 3; c++) deg[mesh->faces[3 * f + c]]++;
    }

    /* CSR offsets: vertices carry deg[v] normals, scattered points carry 1 */
    out.un_off[0] = 0;
    for (int v = 0; v < Nv; v++)   out.un_off[v + 1]     = out.un_off[v] + deg[v];
    for (int j = 0; j < Nscat; j++) out.un_off[Nv + j + 1] = out.un_off[Nv + j] + 1;
    int M = out.un_off[N];
    out.un = malloc(sizeof(s_point) * (size_t)(M > 0 ? M : 1));
    if (!out.un) goto fail;

    /* fill sample coordinates: mesh vertices, then scattered points */
    memcpy(out.samples.p, mesh->points.p, sizeof(s_point) * (size_t)Nv);
    for (int j = 0; j < Nscat; j++)
        out.samples.p[Nv + j] = ((s_scatpt*)scat.items)[j].p;
    out.samples.N = N;
    out.Nv = Nv;

    /* fill umbrella normals: per-vertex fan (running cursor), then scattered */
    int *cursor = malloc(sizeof(int) * (size_t)(Nv > 0 ? Nv : 1));
    if (!cursor) goto fail;
    for (int v = 0; v < Nv; v++) cursor[v] = out.un_off[v];
    for (int f = 0; f < Nf; f++) {
        s_point nf;
        if (!face_unit_normal(mesh, f, EPS_DEG, &nf)) continue;
        for (int c = 0; c < 3; c++) {
            int v = mesh->faces[3 * f + c];
            out.un[cursor[v]++] = nf;
        }
    }
    for (int j = 0; j < Nscat; j++) {
        int f = ((s_scatpt*)scat.items)[j].f;
        s_point nf;
        face_unit_normal(mesh, f, EPS_DEG, &nf);   /* scattered face is valid */
        out.un[out.un_off[Nv + j]] = nf;
    }
    free(cursor);

    /* local size R */
    if (r_sample > 0.0) {
        out.R_uniform = r_sample;
        out.R = NULL;
    } else {
        /* vertices-only fallback: R(v) = max incident-face circumradius */
        out.R = calloc((size_t)N, sizeof(double));
        if (!out.R) goto fail;
        for (int f = 0; f < Nf; f++) {
            double rc = face_circumradius(mesh, f);
            for (int c = 0; c < 3; c++) {
                int v = mesh->faces[3 * f + c];
                if (rc > out.R[v]) out.R[v] = rc;
            }
        }
    }

    free(deg);
    dynarray_free(&scat);
    return out;

fail:
    free(deg);
    dynarray_free(&scat);
    free_samples(&out);
    return out;
}

/* ----- inside/outside oracle: generalized winding number ------------------
 * Replaces the former domain-CDT walk (tetrahedralize_domain_flagged +
 * point location).  Rationale: the CDT oracle requires a successful global
 * construction, which fails (cleanly, but after minutes) on degenerate cell
 * surfaces -- e.g. Voronoi cells containing sub-ulp slit pockets, see
 * TRIMESH_REPAIR_PLAN.md.  The winding number is a pure per-query summation
 * over the faces: no construction, no failure mode, robust for any valid
 * closed trimesh; benchmarked at 20000/20000 agreement with the CDT oracle
 * wherever the latter builds, at comparable total cost for medax-scale query
 * counts (tests/bench_inside.c). */
/* grid: inside-classification cache over the winding number (exact agreement
 * by construction; surface-adjacent voxels fall back to the winding sum).
 * NULL = plain winding (build failure, or MEDAX_NOGRID=1 for A/B timing). */
static int mesh_inside(const s_trimesh *mesh, const s_trimesh_inside_grid *grid,
                       s_point p)
{
    if (grid) return point_in_trimesh_grid(grid, mesh, p);
    return point_in_trimesh_winding(mesh, p);
}

/* ----- Step 3: Delaunay edge enumeration + angle/ratio/inner filter ------- */
/* Edge set: key = packed (min<<32 | max) sample-id pair, value = 1 if the edge
 * is a selected INNER medial-axis edge, 0 if merely seen-but-rejected. */
static size_t edge_hash(const void *key)
{
    uint64_t x = *(const uint64_t*)key;              /* splitmix64 finalizer */
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33;
    return (size_t)x;
}
static bool edge_eq(const void *a, const void *b)
{ return *(const uint64_t*)a == *(const uint64_t*)b; }

static inline uint64_t edge_key(int a, int b)
{
    if (a > b) { int t = a; a = b; b = t; }
    return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b;
}

static inline double R_of(const s_medax_samples *S, int p)
{ return S->R ? S->R[p] : S->R_uniform; }

/* angle condition at sample p: the edge (unit dir e) makes an angle > theta0
 * with EVERY umbrella face plane, i.e. |dot(n_hat, e)| > sin(theta0) for all. */
static int angle_sel_at(const s_medax_samples *S, int p, s_point e, double sin0)
{
    int lo = S->un_off[p], hi = S->un_off[p + 1];
    if (hi <= lo) return 0;                    /* no umbrella normals */
    for (int k = lo; k < hi; k++) {
        s_point n = S->un[k];
        double d = fabs(n.x * e.x + n.y * e.y + n.z * e.z);
        if (d <= sin0) return 0;
    }
    return 1;
}

/* ratio condition at sample p: edge much longer than the local sample scale. */
static inline int ratio_sel_at(const s_medax_samples *S, int p, double len,
                               double rho_eff)
{ return len > rho_eff * R_of(S, p); }

typedef struct edge_counts {
    long unique, angle, ratio, candidate, selected;
} s_edge_counts;

/* Build the selected-inner-edge hash set from the DT. The inner test checks
 * the edge midpoint against the cell surface (winding number). Returns 0 on
 * error. */
static int build_selected_edges(const s_scplx *dt, const s_medax_samples *S,
                                const s_trimesh *mesh,
                                const s_trimesh_inside_grid *grid,
                                double theta0, double rho0,
                                s_hash_table *out, s_edge_counts *cnt)
{
    size_t expected = (size_t)7 * (size_t)S->samples.N + 16;   /* ~7 edges/vert */
    if (!hash_init(out, sizeof(uint64_t), sizeof(uint8_t), expected | 1,
                   expected, edge_hash, edge_eq, NULL))
        return 0;

    double sin0 = sin(theta0);
    double rho_eff = MEDAX_RHO_EFF(rho0);
    *cnt = (s_edge_counts){0};

    for (const s_ncell *c = dt->head; c; c = c->next) {
        for (int i = 0; i < 4; i++) for (int j = i + 1; j < 4; j++) {
            int a = c->vertex_id[i], b = c->vertex_id[j];
            uint64_t key = edge_key(a, b);
            if (hash_get(out, &key)) continue;         /* already evaluated */
            cnt->unique++;

            s_point P = S->samples.p[a], Q = S->samples.p[b];
            s_point d = vsub(Q, P);
            double len = vnorm(d);
            uint8_t selected = 0;

            if (len > 0.0) {
                s_point e = {{{ d.x / len, d.y / len, d.z / len }}};
                int ang = angle_sel_at(S, a, e, sin0) || angle_sel_at(S, b, e, sin0);
                int rat = ratio_sel_at(S, a, len, rho_eff) ||
                          ratio_sel_at(S, b, len, rho_eff);
                if (ang) cnt->angle++;
                if (rat) cnt->ratio++;
                if (ang || rat) {
                    cnt->candidate++;
                    s_point m = {{{ (P.x + Q.x) * 0.5,
                                    (P.y + Q.y) * 0.5,
                                    (P.z + Q.z) * 0.5 }}};
                    if (mesh_inside(mesh, grid, m) == 1) {
                        selected = 1;
                        cnt->selected++;
                    }
                }
            }
            if (hash_insert(out, &key, &selected) == 0) { hash_free(out); return 0; }
        }
    }
    return 1;
}

/* ----- Steps 4-5: emit Voronoi vertices + radii, and (opt) the medial graph - */
typedef struct emit_counts {
    long emitted, degenerate, outside, edges;
} s_emit_counts;

static inline uint64_t ptr_key(const void *p) { return (uint64_t)(uintptr_t)p; }

typedef struct med_edge { int i, j; double w; } s_med_edge;

/* Assemble the CSR s_medial_graph from an undirected edge list (each edge stored
 * in BOTH endpoints' slices). Returns NULL on allocation failure. */
static s_medial_graph *medgraph_from_edges(int N, const s_med_edge *edges, int E)
{
    s_medial_graph *g = malloc(sizeof(s_medial_graph));
    int    *adj_off = calloc((size_t)N + 1, sizeof(int));
    int    *adj     = (E > 0) ? malloc(sizeof(int)    * 2 * (size_t)E) : NULL;
    double *adj_w   = (E > 0) ? malloc(sizeof(double) * 2 * (size_t)E) : NULL;
    int    *cursor  = malloc(sizeof(int) * (size_t)(N > 0 ? N : 1));
    if (!g || !adj_off || (E > 0 && (!adj || !adj_w)) || !cursor) {
        free(g); free(adj_off); free(adj); free(adj_w); free(cursor);
        return NULL;
    }
    for (int e = 0; e < E; e++) { adj_off[edges[e].i + 1]++; adj_off[edges[e].j + 1]++; }
    for (int i = 0; i < N; i++)   adj_off[i + 1] += adj_off[i];
    for (int i = 0; i < N; i++)   cursor[i] = adj_off[i];
    for (int e = 0; e < E; e++) {
        int i = edges[e].i, j = edges[e].j; double w = edges[e].w;
        adj[cursor[i]] = j; adj_w[cursor[i]] = w; cursor[i]++;
        adj[cursor[j]] = i; adj_w[cursor[j]] = w; cursor[j]++;
    }
    free(cursor);
    g->N = N; g->adj_off = adj_off; g->adj = adj; g->adj_w = adj_w;
    return g;
}

static void free_medial_graph(s_medial_graph *g)
{
    if (!g) return;
    free(g->adj_off); free(g->adj); free(g->adj_w);
    free(g);
}

/* Step 5: join emitted tets that share a face (dual Voronoi edge) into the CSR
 * medial graph. tet2node maps a tet pointer to its emitted node index. Returns
 * NULL on error; sets *n_edges. */
static s_medial_graph *build_medgraph(const s_scplx *dt, s_hash_table *tet2node,
                                      const s_point *P, int N, long *n_edges)
{
    *n_edges = 0;
    s_dynarray edges = dynarray_initialize(sizeof(s_med_edge), 1024);
    if (!edges.items) return NULL;

    for (const s_ncell *c = dt->head; c; c = c->next) {
        uint64_t pk = ptr_key(c);
        int *ni = hash_get(tet2node, &pk);
        if (!ni) continue;
        for (int k = 0; k < 4; k++) {
            const s_ncell *nb = c->opposite[k];
            if (!nb) continue;
            uint64_t nk = ptr_key(nb);
            int *nj = hash_get(tet2node, &nk);
            if (!nj || *ni >= *nj) continue;      /* store each edge once (i<j) */
            s_med_edge e = { *ni, *nj, vnorm(vsub(P[*ni], P[*nj])) };
            if (!dynarray_push(&edges, &e)) { dynarray_free(&edges); return NULL; }
        }
    }

    s_medial_graph *g = medgraph_from_edges(N, (s_med_edge*)edges.items, (int)edges.N);
    dynarray_free(&edges);
    if (g) *n_edges = g->adj_off[N] / 2;
    return g;
}

/* A DT tet's circumcenter is a vertex of the dual Voronoi facet of each of its
 * 6 edges; emit it once iff any of those 6 edges is a selected inner edge (one
 * point per tet, auto-deduplicated). If build_graph, also build the medial
 * graph over the emitted points. */
static s_medax emit_medax(const s_scplx *dt, const s_trimesh *mesh,
                          const s_trimesh_inside_grid *grid,
                          s_hash_table *sel, double EPS_DEG, int inside_filter,
                          int build_graph, s_emit_counts *cnt)
{
    s_dynarray pts = dynarray_initialize(sizeof(s_point), 1024);
    s_dynarray rad = dynarray_initialize(sizeof(double), 1024);
    s_hash_table tet2node;                 /* tet pointer -> emitted node index */
    if (!pts.items || !rad.items) {
        dynarray_free(&pts); dynarray_free(&rad);
        return medax_NAN;
    }
    if (build_graph) {
        size_t exp = (size_t)(dt->N_ncells > 0 ? dt->N_ncells : 1);
        if (!hash_init(&tet2node, sizeof(uint64_t), sizeof(int), (2 * exp) | 1, exp,
                       edge_hash, edge_eq, NULL)) {
            dynarray_free(&pts); dynarray_free(&rad);
            return medax_NAN;
        }
    }
    *cnt = (s_emit_counts){0};

    /* Step 4: emit points (recording tet -> node index when building a graph). */
    for (const s_ncell *c = dt->head; c; c = c->next) {
        int emit = 0;
        for (int i = 0; i < 4 && !emit; i++)
            for (int j = i + 1; j < 4; j++) {
                uint64_t key = edge_key(c->vertex_id[i], c->vertex_id[j]);
                uint8_t *v = hash_get(sel, &key);
                if (v && *v) { emit = 1; break; }
            }
        if (!emit) continue;

        s_point verts[4]; extract_vertices_ncell(dt, c, verts);
        s_point cc;
        if (!circumcentre_tetrahedron(verts, EPS_DEG, &cc)) { cnt->degenerate++; continue; }
        if (inside_filter && mesh_inside(mesh, grid, cc) != 1) { cnt->outside++; continue; }

        double r = vnorm(vsub(cc, verts[0]));   /* all 4 verts equidistant */
        int node = (int)pts.N;
        int ok = dynarray_push(&pts, &cc) && dynarray_push(&rad, &r);
        if (ok && build_graph) {
            uint64_t pk = ptr_key(c);
            ok = (hash_insert(&tet2node, &pk, &node) != 0);
        }
        if (!ok) {
            dynarray_free(&pts); dynarray_free(&rad);
            if (build_graph) hash_free(&tet2node);
            return medax_NAN;
        }
        cnt->emitted++;
    }

    s_medax ma = medax_NAN;
    ma.verts.N = (int)pts.N;
    ma.verts.p = (s_point*)pts.items;      /* transfer ownership */
    ma.radius  = (double*)rad.items;

    if (build_graph) {
        long ne = 0;
        s_medial_graph *g = build_medgraph(dt, &tet2node, ma.verts.p, ma.verts.N, &ne);
        hash_free(&tet2node);
        if (!g) { free_medax(&ma); return medax_NAN; }
        cnt->edges = ne;
        ma.graph = g;
    }
    return ma;
}

/* Drop the records of samples the DT merged into an earlier one (kept[i] < 0),
 * in place; survivors keep their order, matching the DT's compacted points. */
static void compact_samples(s_medax_samples *S, const int *kept)
{
    int k = 0, o = 0, nv = 0;
    for (int i = 0; i < S->samples.N; i++) {
        if (kept[i] < 0) continue;
        int lo = S->un_off[i], n = S->un_off[i + 1] - lo;
        S->un_off[k] = o;
        for (int t = 0; t < n; t++) S->un[o++] = S->un[lo + t];
        S->samples.p[k] = S->samples.p[i];
        if (S->R) S->R[k] = S->R[i];
        if (i < S->Nv) nv++;
        k++;
    }
    S->un_off[k] = o;
    S->samples.N = k;
    S->Nv = nv;
}

/* ------------------------------------------------------------------------- */
s_medax medax_from_trimesh(const s_trimesh *mesh, double r_sample,
                           double theta0, double rho0,
                           double (*randd01)(void*), void *rctx,
                           double EPS_DEG, double TOL, bool want_graph)
{
    if (!trimesh_is_valid(mesh)) return medax_NAN;
    if (theta0 <= 0.0) theta0 = MEDAX_THETA0_DEFAULT;
    if (rho0   <= 0.0) rho0   = MEDAX_RHO0_DEFAULT;
    if (r_sample > 0.0 && !randd01) return medax_NAN;

    /* Step 2: samples + umbrella records */
    struct timespec _mt; clock_gettime(CLOCK_MONOTONIC, &_mt);
    double _t0 = _mt.tv_sec + 1e-9 * _mt.tv_nsec, _t1;
#define MEDAX_TICK(msg) do { \
        clock_gettime(CLOCK_MONOTONIC, &_mt); \
        _t1 = _mt.tv_sec + 1e-9 * _mt.tv_nsec; \
        if (getenv("MEDAX_PROFILE")) \
            fprintf(stderr, "[medax][t] %-10s %.3fs\n", msg, _t1 - _t0); \
        _t0 = _t1; \
    } while (0)

    s_medax_samples S = build_samples(mesh, r_sample, randd01, rctx, EPS_DEG);
    if (!S.samples.p) return medax_NAN;
    MEDAX_TICK("sampling");

    /* Step 2: Delaunay triangulation of the samples (lean, no mirrors). Samples
     * within TOL are merged (the first survives, so mesh vertices win). */
    int *kept = malloc(sizeof(int) * (size_t)S.samples.N);
    if (!kept) { free_samples(&S); return medax_NAN; }
    for (int i = 0; i < S.samples.N; i++) kept[i] = i;
    s_dt_builder b = dt_builder_begin(&S.samples, NULL, TOL, NULL, NULL, NULL);
    if (!b._stack) { free(kept); free_samples(&S); return medax_NAN; }
    int Nreal = 0;
    s_scplx dt = dt_builder_end(&b, false, &Nreal, kept, S.samples.N);
    if (Nreal != S.samples.N) compact_samples(&S, kept);
    free(kept);

    MEDAX_TICK("sample-DT");

    /* Inside/outside oracle: generalized winding number against the input
     * surface itself (see mesh_inside), cached in a classification grid --
     * O(1) per query away from the surface, exact winding fallback in
     * surface-adjacent voxels.  No construction step that can fail: if the
     * grid build fails (or MEDAX_NOGRID is set) every query just runs the
     * plain winding sum. */
    s_trimesh_inside_grid grid;
    const s_trimesh_inside_grid *gp = NULL;
    if (!getenv("MEDAX_NOGRID") &&
        trimesh_inside_grid_build(mesh, 128, &grid))
        gp = &grid;
    MEDAX_TICK("grid");

    /* Step 3: enumerate + filter Delaunay edges (angle/ratio/inner) */
    s_hash_table sel_edges;
    s_edge_counts cnt;
    if (!build_selected_edges(&dt, &S, mesh, gp, theta0, rho0, &sel_edges, &cnt)) {
        if (gp) trimesh_inside_grid_free(&grid);
        free_complex(&dt);
        free_samples(&S);
        return medax_NAN;
    }
    MEDAX_TICK("edge-sel");

    /* Step 4 (+5): emit circumcenters of tets touching a selected edge, and
     * optionally the medial graph over them. */
    s_emit_counts ec;
    s_medax ma = emit_medax(&dt, mesh, gp, &sel_edges, EPS_DEG, /*inside_filter=*/1,
                            want_graph, &ec);
    if (gp) trimesh_inside_grid_free(&grid);

    MEDAX_TICK("emit");
#undef MEDAX_TICK
    hash_free(&sel_edges);
    free_complex(&dt);
    free_samples(&S);
    return ma;
}

void free_medax(s_medax *ma)
{
    if (!ma) return;
    free(ma->verts.p);
    free(ma->radius);
    free_medial_graph(ma->graph);
    *ma = medax_NAN;
}

int medax_write_vtk(const s_medax *ma, const char *path)
{
    if (!ma || !path) return 1;

    FILE *f = fopen(path, "w");
    if (!f) return 1;

    int N = ma->verts.p ? ma->verts.N : 0;

    fprintf(f, "# vtk DataFile Version 3.0\n");
    fprintf(f, "medial axis (approx, Dey-Zhao)\n");
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET POLYDATA\n");

    fprintf(f, "POINTS %d double\n", N);
    for (int i = 0; i < N; i++) {
        s_point p = ma->verts.p[i];
        fprintf(f, "%.17g %.17g %.17g\n", p.x, p.y, p.z);
    }

    /* one VERTICES cell per point so ParaView renders them as a point cloud */
    fprintf(f, "VERTICES %d %d\n", N, 2 * N);
    for (int i = 0; i < N; i++) {
        fprintf(f, "1 %d\n", i);
    }

    /* medial graph edges as LINES cells (each undirected edge once, from i<j) */
    if (N > 0 && ma->graph && ma->graph->adj_off && ma->graph->adj) {
        const s_medial_graph *g = ma->graph;
        int nE = g->adj_off[N] / 2;
        fprintf(f, "LINES %d %d\n", nE, 3 * nE);
        for (int i = 0; i < N; i++)
            for (int k = g->adj_off[i]; k < g->adj_off[i + 1]; k++)
                if (g->adj[k] > i) fprintf(f, "2 %d %d\n", i, g->adj[k]);
    }

    if (ma->radius) {
        fprintf(f, "POINT_DATA %d\n", N);
        fprintf(f, "SCALARS radius double 1\n");
        fprintf(f, "LOOKUP_TABLE default\n");
        for (int i = 0; i < N; i++) {
            fprintf(f, "%.17g\n", ma->radius[i]);
        }
    }

    fclose(f);
    return 0;
}
