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

/* ----- performance note on medax_center (the canonical center) -------------
 * The center is argmin_i sum_j w_j d(i,j)^2 over the medial graph, where w_j is
 * the governed volume of ball j (its per-ball contribution vol(B_j INTERSECT
 * L_j) to the union of medial balls -- the Laguerre/power partition, computed by
 * volume_contribution_spheres; see MEDIAL_AXIS_BT.md sec.3). Each energy
 * evaluation needs a single-source distance field (one Dijkstra), so the cost is
 * ~ (number of candidate nodes evaluated) x (per-Dijkstra cost). The current
 * warm-started multi-start descent (expanding-ring polish + trajectory
 * memoization) is exact and runs ~0.7-0.9 s on a ~15 k-node lung lobe.
 *
 * TRIED and REVERTED: a coarse-to-fine (multigrid) scheme -- coarsen the graph
 * by heavy-edge matching, solve the small coarse graph exactly, map the result
 * back and refine on the fine graph. It was exact (validated against brute force)
 * but gave NO speedup: the fine refinement still has to run the H-hop polish
 * *confirmation* (~O(H^2) full-graph Dijkstras to certify the minimum), and that
 * floor dominates regardless of how good the warm start is. The plain descent
 * already sits on that same floor, so the coarse solve was pure overhead.
 *
 * The real remaining lever is to attack that floor directly: the confirmation
 * evaluates a ring of nodes whose Dijkstras are INDEPENDENT, so it parallelizes
 * (as the example's OpenMP brute-force check demonstrates ~10x). To do it safely
 * here: (a) give each thread its OWN scratch heap (cen_descent currently shares
 * one -- it would race); (b) this library is meant to run in parallel at a higher
 * level (e.g. many cells at once), so guard any inner "#pragma omp parallel" on
 * omp_get_active_level() / rely on nested parallelism being off, to avoid
 * oversubscribing cores. Prefer a single parallel level chosen by the caller. */

#include "medax.h"

#include "delaunay.h"
#include "scplx.h"
#include "dynarray.h"
#include "hash.h"
#include "volsph.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#define MEDAX_MARGIN 1e-6              /* barycentric interior margin (Step 2) */
static const double MEDAX_SQRT3 = 1.7320508075688772;

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

/* ------------------------------------------------------------------------- */
s_medax medax_from_trimesh(const s_trimesh *mesh, double r_sample,
                           double theta0, double rho0,
                           double (*randd01)(void*), void *rctx,
                           double TOL, double EPS_DEG, bool build_graph)
{
    if (!trimesh_is_valid(mesh)) return medax_NAN;
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

    /* Step 2: Delaunay triangulation of the samples (lean, no mirrors) */
    s_dt_builder b = dt_builder_begin(&S.samples, NULL, TOL, NULL, NULL);
    if (!b._stack) { free_samples(&S); return medax_NAN; }
    int Nreal = 0;
    s_scplx dt = dt_builder_end(&b, false, &Nreal, NULL, 0);
    if (Nreal != S.samples.N) {
        fprintf(stderr, "[medax] DT deduplicated points (Nreal=%d != Ns=%d); "
                        "aborting -- shrink margin/TOL interaction.\n",
                        Nreal, S.samples.N);
        free_complex(&dt);
        free_samples(&S);
        return medax_NAN;
    }

    fprintf(stderr, "[medax] Step 2 OK: Nv=%d Nscat=%d Ns=%d (DT built, Nreal==Ns).\n",
                    S.Nv, S.samples.N - S.Nv, S.samples.N);
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
    fprintf(stderr, "[medax] Step 3 OK: edges unique=%ld angle=%ld ratio=%ld "
                    "candidate=%ld selected(inner)=%ld.\n",
                    cnt.unique, cnt.angle, cnt.ratio, cnt.candidate, cnt.selected);
    MEDAX_TICK("edge-sel");

    /* Step 4 (+5): emit circumcenters of tets touching a selected edge, and
     * optionally the medial graph over them. */
    s_emit_counts ec;
    s_medax ma = emit_medax(&dt, mesh, gp, &sel_edges, EPS_DEG, /*inside_filter=*/1,
                            build_graph, &ec);
    if (gp) trimesh_inside_grid_free(&grid);
    if (build_graph)
        fprintf(stderr, "[medax] Step 4-5 OK: points=%ld edges=%ld (skipped "
                        "degenerate=%ld, outside=%ld).\n",
                        ec.emitted, ec.edges, ec.degenerate, ec.outside);
    else
        fprintf(stderr, "[medax] Step 4 OK: points=%ld (skipped degenerate=%ld, "
                        "outside=%ld).\n", ec.emitted, ec.degenerate, ec.outside);

    MEDAX_TICK("emit");
#undef MEDAX_TICK
    hash_free(&sel_edges);
    free_complex(&dt);
    free_samples(&S);
    return ma;
}

s_points medax_debug_sample_surface(const s_trimesh *mesh, double r_sample,
                                    double (*randd01)(void*), void *rctx,
                                    double EPS_DEG)
{
    if (!trimesh_is_valid(mesh)) return points_NAN;
    if (r_sample > 0.0 && !randd01) return points_NAN;

    s_medax_samples S = build_samples(mesh, r_sample, randd01, rctx, EPS_DEG);
    if (!S.samples.p) return points_NAN;

    s_points pts = S.samples;   /* transfer ownership of the coordinate array */
    S.samples.p = NULL;
    free_samples(&S);
    return pts;
}

/* ===== Canonical center: volume-weighted geodesic Frechet mean ============ */
/* The center minimizes f(y) = sum_j r_j^3 * d(y,j)^2 over medial nodes (the
 * geodesic Frechet / Karcher mean of the medial axis, weighted by ball volume).
 * Small components are solved exactly (all-pairs); large ones by warm-started
 * multi-start gradient descent -- one distance field per iteration instead of
 * all-pairs -- following:
 *
 *   [1] C. Mancinelli, E. Puppo, "Computing the Riemannian center of mass on
 *       meshes", Computer Aided Geometric Design 103 (2023) 102203.
 *       (Newton/gradient descent on the squared-distance energy; the gradient is
 *        the sum of log maps, so one distance field from the iterate gives it;
 *        warm start crucial; detects non-convexity via the Hessian.)
 *   [2] F. M. Rygaard, S. Hauberg, S. Markvorsen, "Simultaneous Optimization of
 *       Geodesics and Frechet Means" (GEORCE-FM), arXiv:2511.04301 (2025).
 *       (Global + local-quadratic convergence; stochastic/adaptive extension for
 *        large data; same convex-ball uniqueness caveat as Karcher 1977.)
 *
 * All static; the only public entry points are medax_center / medax_center_seeds.
 * Graph traversal (Dijkstra, energy, descent) is an internal detail here -- no
 * generic graph API is exposed. */

typedef struct cen_hn { double d; int v; } s_cen_hn;
static void cen_hup(s_cen_hn *h, int i)
{ while (i>0){int p=(i-1)/2; if(h[p].d<=h[i].d)break; s_cen_hn t=h[p];h[p]=h[i];h[i]=t;i=p;} }
static void cen_hdown(s_cen_hn *h, int n, int i)
{ for(;;){int l=2*i+1,r=2*i+2,m=i; if(l<n&&h[l].d<h[m].d)m=l; if(r<n&&h[r].d<h[m].d)m=r;
          if(m==i)break; s_cen_hn t=h[m];h[m]=h[i];h[i]=t;i=m;} }

/* Dijkstra over the CSR; fills dist[] (INFINITY if unreachable). If fh != NULL,
 * fh[j] = the first-hop neighbour of src toward j (the geodesic tangent).
 * h is a caller-owned scratch heap of size >= g->adj_off[N] + 8 (reused across
 * the many calls a single center computation makes, to avoid per-call malloc). */
static void cen_dijkstra(const s_medial_graph *g, int src, double *dist, int *fh,
                         s_cen_hn *h)
{
    int N = g->N;
    for (int i = 0; i < N; i++) { dist[i] = INFINITY; if (fh) fh[i] = -1; }
    int hn = 0; dist[src] = 0.0; h[hn++] = (s_cen_hn){ 0.0, src };
    while (hn > 0) {
        s_cen_hn top = h[0]; h[0] = h[--hn]; cen_hdown(h, hn, 0);
        if (top.d > dist[top.v]) continue;
        for (int k = g->adj_off[top.v]; k < g->adj_off[top.v+1]; k++) {
            int w = g->adj[k]; double nd = top.d + g->adj_w[k];
            if (nd < dist[w]) {
                dist[w] = nd; if (fh) fh[w] = (top.v == src) ? w : fh[top.v];
                h[hn++] = (s_cen_hn){ nd, w }; cen_hup(h, hn-1);
            }
        }
    }
}

/* f(node) = sum_j w_j d(node,j)^2 over component lc, from a filled dist[].
 * mass[j] = w_j = governed volume of ball j (see cen_weights). */
static double cen_energy(const double *dist, const double *mass,
                         const int *comp, int lc, int N)
{
    double E = 0.0;
    for (int j = 0; j < N; j++) {
        if (comp[j] != lc || dist[j] == INFINITY) continue;
        E += mass[j] * dist[j] * dist[j];
    }
    return E;
}

/* exact all-pairs argmin (O(N) Dijkstra); returns node, sets *E. */
static int cen_exact(const s_medial_graph *g, const double *mass,
                     const int *comp, int lc, double *dist, s_cen_hn *heap, double *E)
{
    int N = g->N, best = -1; double bE = INFINITY;
    for (int i = 0; i < N; i++) {
        if (comp[i] != lc) continue;
        cen_dijkstra(g, i, dist, NULL, heap);
        double e = cen_energy(dist, mass, comp, lc, N);
        if (e < bE) { bE = e; best = i; }
    }
    *E = bE; return best;
}

/* Straightest walk on the graph from y in ambient direction (dx,dy,dz) for up
 * to arc-length L; returns the landing node (>= one edge unless y is a sink). */
static int cen_walk(const s_medial_graph *g, const s_point *pos, int y,
                    double dx, double dy, double dz, double L)
{
    int cur = y; double rem = L;
    for (int steps = 0; steps < 1000000; steps++) {
        int off = g->adj_off[cur], deg = g->adj_off[cur+1] - off;
        int best = -1; double bestdot = 0.0, bestlen = 0.0;
        for (int t = 0; t < deg; t++) {
            int k = g->adj[off+t];
            double ex=pos[k].x-pos[cur].x, ey=pos[k].y-pos[cur].y, ez=pos[k].z-pos[cur].z;
            double el = sqrt(ex*ex+ey*ey+ez*ez); if (el == 0) continue;
            double dot = (ex*dx + ey*dy + ez*dz) / el;
            if (dot > bestdot) { bestdot = dot; best = k; bestlen = el; }
        }
        if (best < 0) break;                  /* no forward edge */
        cur = best; rem -= bestlen;
        if (rem <= 0.0) break;
    }
    return cur;
}

#define CEN_MAX_ITERS 500

/* Warm-started descent from `start`: gradient big-step (Karcher/Exp step in the
 * mean-log direction, with a backtracking line search), a best-neighbour
 * fallback for discretization stalls, and -- when that stalls at a 1-hop local
 * minimum -- an H-hop escape that expands ring by ring, stopping at the FIRST
 * ring holding a lower-energy node, jumps there and restarts the descent (so its
 * cost is proportional to the barrier distance, not H; H = polish_hops, <=0
 * disables it). O(k) Dijkstra. heap is caller scratch.
 *
 * resolved (shared across a multi-start, may be NULL): resolved[node] = the
 * descent result reachable from node, or -1. Because the descent is
 * deterministic, if it lands on an already-resolved node it MUST reach the same
 * result, so it short-circuits there (and the whole trajectory is then stamped).
 * This makes later seeds that merge onto an earlier path skip the expensive
 * final confirmation -- exact, no accuracy loss. Returns node, sets *E. */
static int cen_descent(const s_medial_graph *g, const s_point *pos,
                       const double *mass, const int *comp, int lc,
                       int start, int H, s_cen_hn *heap, int *resolved, double *E)
{
    int N = g->N;
    double *da = malloc(sizeof(double)*(size_t)N), *db = malloc(sizeof(double)*(size_t)N);
    int    *fa = malloc(sizeof(int)*(size_t)N),    *fb = malloc(sizeof(int)*(size_t)N);
    int    *seen = (H > 0) ? calloc((size_t)N, sizeof(int)) : NULL;
    int    *bq   = (H > 0) ? malloc(sizeof(int)*(size_t)N) : NULL;
    int    *traj = resolved ? malloc(sizeof(int)*(CEN_MAX_ITERS + 1)) : NULL;
    if (!da || !db || !fa || !fb || (H > 0 && (!seen || !bq)) || (resolved && !traj)) {
        free(da);free(db);free(fa);free(fb);free(seen);free(bq);free(traj);
        *E = INFINITY; return -1;
    }
    int stamp = 0, tn = 0;

    int y = start;
    cen_dijkstra(g, y, da, fa, heap);
    double fy = cen_energy(da, mass, comp, lc, N);

    for (int iter = 0; iter < CEN_MAX_ITERS; iter++) {
        /* merge short-circuit: this node already leads to a known result */
        if (resolved && resolved[y] >= 0) {
            int R = resolved[y];
            cen_dijkstra(g, R, db, NULL, heap);
            double ER = cen_energy(db, mass, comp, lc, N);
            for (int i = 0; i < tn; i++) resolved[traj[i]] = R;
            free(da);free(db);free(fa);free(fb);free(seen);free(bq);free(traj);
            *E = ER; return R;
        }
        if (resolved && tn <= CEN_MAX_ITERS) traj[tn++] = y;

        /* mean log vector V = (1/M) sum_j m_j d(y,j) * unit(pos[fh]-pos[y]) */
        double Vx=0, Vy=0, Vz=0, M=0;
        for (int j = 0; j < N; j++) {
            if (comp[j] != lc || fa[j] < 0 || da[j] == INFINITY) continue;
            double m = mass[j]; M += m;
            double ex=pos[fa[j]].x-pos[y].x, ey=pos[fa[j]].y-pos[y].y, ez=pos[fa[j]].z-pos[y].z;
            double el = sqrt(ex*ex+ey*ey+ez*ez); if (el == 0) continue;
            double w = m * da[j] / el; Vx += w*ex; Vy += w*ey; Vz += w*ez;
        }
        int moved = 0;
        if (M > 0) {
            Vx/=M; Vy/=M; Vz/=M;
            double Lstep = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
            if (Lstep > 0) {
                double dhx=Vx/Lstep, dhy=Vy/Lstep, dhz=Vz/Lstep, L = Lstep;
                for (int ls = 0; ls < 24; ls++) {           /* backtracking */
                    int y2 = cen_walk(g, pos, y, dhx, dhy, dhz, L);
                    if (y2 != y) {
                        cen_dijkstra(g, y2, db, fb, heap);
                        double f2 = cen_energy(db, mass, comp, lc, N);
                        if (f2 < fy) {
                            y = y2; fy = f2;
                            double *td=da; da=db; db=td; int *tf=fa; fa=fb; fb=tf;
                            moved = 1; break;
                        }
                    }
                    L *= 0.5; if (L < 1e-15) break;
                }
            }
        }
        if (moved) continue;
        /* fallback: pick the best of all <=4 neighbours */
        int off = g->adj_off[y], deg = g->adj_off[y+1]-off, bk = -1; double bf = fy;
        for (int t = 0; t < deg; t++) {
            int k = g->adj[off+t];
            cen_dijkstra(g, k, db, NULL, heap);
            double fk = cen_energy(db, mass, comp, lc, N);
            if (fk < bf) { bf = fk; bk = k; }
        }
        if (bk >= 0) { y = bk; fy = bf; cen_dijkstra(g, y, da, fa, heap); continue; }

        /* 1-hop local minimum: expanding-ring escape. Rings 2..H are evaluated
         * one at a time (ring 1 was just checked by the fallback); stop at the
         * first ring with a lower-energy node and jump to the best in it. */
        if (H <= 0) break;
        stamp++;
        seen[y] = stamp;
        int qt = 0;                            /* frontier = ring-1 nodes */
        for (int k = g->adj_off[y]; k < g->adj_off[y+1]; k++) {
            int w = g->adj[k];
            if (seen[w] != stamp) { seen[w] = stamp; bq[qt++] = w; }
        }
        int ring_start = 0, depth = 1, pbest = -1; double pf = fy;
        while (ring_start < qt && depth <= H) {
            int ring_end = qt;
            if (depth >= 2) {                  /* evaluate this ring */
                for (int i = ring_start; i < ring_end; i++) {
                    int u = bq[i]; if (comp[u] != lc) continue;
                    cen_dijkstra(g, u, db, NULL, heap);
                    double fu = cen_energy(db, mass, comp, lc, N);
                    if (fu < pf) { pf = fu; pbest = u; }
                }
                if (pbest >= 0) break;          /* nearest improving ring: stop */
            }
            if (depth < H)                      /* expand to next ring */
                for (int i = ring_start; i < ring_end; i++) {
                    int u = bq[i];
                    for (int k = g->adj_off[u]; k < g->adj_off[u+1]; k++) {
                        int w = g->adj[k];
                        if (seen[w] != stamp) { seen[w] = stamp; bq[qt++] = w; }
                    }
                }
            ring_start = ring_end; depth++;
        }
        if (pbest < 0) break;                  /* no lower node within H hops: done */
        y = pbest; fy = pf; cen_dijkstra(g, y, da, fa, heap);   /* jump + restart */
    }
    /* converged to y: stamp the whole trajectory as leading here */
    if (resolved) { for (int i = 0; i < tn; i++) resolved[traj[i]] = y; resolved[y] = y; }
    free(da); free(db); free(fa); free(fb); free(seen); free(bq); free(traj);
    *E = fy; return y;
}

typedef struct rpair { double r; int i; } s_rpair;
static int rpair_desc(const void *a, const void *b)
{
    double ra = ((const s_rpair*)a)->r, rb = ((const s_rpair*)b)->r;
    return (ra < rb) - (ra > rb);
}

/* Descent seeds: deep local-radius maxima with geodesic non-max suppression at
 * separation beta*r_seed (one seed per thick region). Fills seeds[<=Kmax],
 * returns the count. dmin/dist are scratch of size N. */
static int cen_seeds(const s_medial_graph *g, const double *radius,
                     const int *comp, int lc, int Kmax, double beta,
                     int *seeds, double *dmin, double *dist, s_cen_hn *heap)
{
    int N = g->N;
    s_rpair *cand = malloc(sizeof(s_rpair) * (size_t)N); int nc = 0;
    if (!cand) return 0;
    for (int i = 0; i < N; i++) {
        if (comp[i] != lc) continue;
        int off = g->adj_off[i], deg = g->adj_off[i+1]-off, ismax = 1;
        for (int t = 0; t < deg; t++) if (radius[g->adj[off+t]] > radius[i]) { ismax = 0; break; }
        if (ismax) { cand[nc].r = radius[i]; cand[nc].i = i; nc++; }
    }
    qsort(cand, (size_t)nc, sizeof(s_rpair), rpair_desc);   /* deepest first */

    for (int i = 0; i < N; i++) dmin[i] = INFINITY;
    int ns = 0;
    for (int c = 0; c < nc && ns < Kmax; c++) {
        int node = cand[c].i;
        if (dmin[node] <= beta * radius[node]) continue;    /* too close to a seed */
        seeds[ns++] = node;
        cen_dijkstra(g, node, dist, NULL, heap);
        for (int j = 0; j < N; j++) if (dist[j] < dmin[j]) dmin[j] = dist[j];
    }
    free(cand);
    return ns;
}

/* Label connected components into comp[] (BFS via queue[]); returns the largest
 * component's id and sets *out_largest to its size. */
static int cen_components(const s_medial_graph *g, int *comp, int *queue, int *out_largest)
{
    int N = g->N;
    for (int i = 0; i < N; i++) comp[i] = -1;
    int ncomp = 0, largest = 0, root = 0;
    for (int s = 0; s < N; s++) {
        if (comp[s] != -1) continue;
        int qh = 0, qt = 0, sz = 0; comp[s] = ncomp; queue[qt++] = s;
        while (qh < qt) {
            int u = queue[qh++]; sz++;
            for (int k = g->adj_off[u]; k < g->adj_off[u+1]; k++)
                if (comp[g->adj[k]] == -1) { comp[g->adj[k]] = ncomp; queue[qt++] = g->adj[k]; }
        }
        if (sz > largest) { largest = sz; root = s; }
        ncomp++;
    }
    *out_largest = largest;
    return comp[root];
}

/* Node masses w_j for the Frechet energy: the governed volume of each medial
 * ball -- its per-ball contribution vol(B_j INTERSECT L_j) to the union of all
 * medial balls (the Laguerre/power partition), so sum_j w_j == vol(Omega) with
 * ball overlaps counted once (MEDIAL_AXIS_BT.md sec.3). Computed over the FULL
 * ball set (verts/radius); the caller uses only the entries in its component.
 *
 * Returns a freshly malloc'd array of size ma->verts.N, or NULL on allocation
 * failure. If the union decomposition degenerates (non-finite or non-positive
 * total -- e.g. a pathological ball set), falls back to the ball-volume proxy
 * w_j = r_j^3 (the previous, overlap-double-counting weight), which is a safe
 * approximation of the same quantity. */
#define MEDAX_VOL_EPS    1e-10        /* volsph degeneracy epsilon */
#define MEDAX_VOL_TOLDUP 1e-9         /* volsph duplicate-center tolerance */
static double *cen_weights(const s_medax *ma)
{
    int N = ma->verts.N;
    double *w = malloc(sizeof(double) * (size_t)N);
    if (!w) return NULL;

    s_dynarray buff = dynarray_initialize(sizeof(s_ncell *), 64);
    volume_contribution_spheres(&ma->verts, ma->radius, MEDAX_VOL_EPS,
                                MEDAX_VOL_TOLDUP, &buff, w);
    dynarray_free(&buff);

    double sum = 0.0; int ok = 1;
    for (int j = 0; j < N; j++) {
        if (!isfinite(w[j])) { ok = 0; break; }
        if (w[j] < 0.0) w[j] = 0.0;      /* clamp tiny rounding negatives */
        sum += w[j];
    }
    if (!ok || !(sum > 0.0)) {           /* degenerate: fall back to r^3 */
        for (int j = 0; j < N; j++) {
            double r = ma->radius[j]; w[j] = r*r*r;
        }
    }
    return w;
}

int medax_center(const s_medax *ma, int exact_below, int K, double beta,
                 int polish_hops, double *out_energy)
{
    if (out_energy) *out_energy = INFINITY;
    if (!ma || !ma->graph || ma->verts.N == 0 || !ma->radius) return -1;
    const s_medial_graph *g = ma->graph;
    int N = g->N;
    if (N <= 0 || !g->adj_off) return -1;
    if (exact_below < 0) exact_below = 2500;
    if (K <= 0) K = 3;
    if (beta <= 0.0) beta = 2.0;
    if (polish_hops < 0) polish_hops = 8;

    int    *comp  = malloc(sizeof(int)    * (size_t)N);
    int    *queue = malloc(sizeof(int)    * (size_t)N);
    double *dist  = malloc(sizeof(double) * (size_t)N);
    double *mass  = cen_weights(ma);      /* node masses w_j = governed volume */
    s_cen_hn *heap = malloc(sizeof(s_cen_hn) * (size_t)(g->adj_off[N] + 8));
    if (!comp || !queue || !dist || !mass || !heap) {
        free(comp); free(queue); free(dist); free(mass); free(heap); return -1;
    }

    int largest;
    int lc = cen_components(g, comp, queue, &largest);

    int result = -1; double E = INFINITY;
    if (largest < exact_below) {
        result = cen_exact(g, mass, comp, lc, dist, heap, &E);
    } else {
        int    *seeds = malloc(sizeof(int)    * (size_t)K);
        double *dmin  = malloc(sizeof(double) * (size_t)N);
        int    *resolved = malloc(sizeof(int) * (size_t)N);   /* memoization (Fix 4) */
        if (resolved) for (int i = 0; i < N; i++) resolved[i] = -1;
        if (seeds && dmin && resolved) {
            /* seeds are deep-clearance maxima: selected on radius, not mass */
            int ns = cen_seeds(g, ma->radius, comp, lc, K, beta, seeds, dmin, dist, heap);
            for (int s = 0; s < ns; s++) {
                double Es;
                int c = cen_descent(g, ma->verts.p, mass, comp, lc,
                                    seeds[s], polish_hops, heap, resolved, &Es);
                if (c >= 0 && Es < E) { E = Es; result = c; }
            }
        }
        free(seeds); free(dmin); free(resolved);
    }
    free(comp); free(queue); free(dist); free(mass); free(heap);
    if (out_energy) *out_energy = E;
    return result;
}

int medax_center_seeds(const s_medax *ma, int K, double beta, int polish_hops,
                       int *seeds, int *results)
{
    if (!ma || !ma->graph || ma->verts.N == 0 || !ma->radius || !seeds || !results)
        return -1;
    const s_medial_graph *g = ma->graph;
    int N = g->N;
    if (N <= 0 || !g->adj_off) return -1;
    if (K <= 0) K = 3;
    if (beta <= 0.0) beta = 2.0;
    if (polish_hops < 0) polish_hops = 8;

    int    *comp  = malloc(sizeof(int)    * (size_t)N);
    int    *queue = malloc(sizeof(int)    * (size_t)N);
    double *dist  = malloc(sizeof(double) * (size_t)N);
    double *dmin  = malloc(sizeof(double) * (size_t)N);
    double *mass  = cen_weights(ma);      /* node masses w_j = governed volume */
    s_cen_hn *heap = malloc(sizeof(s_cen_hn) * (size_t)(g->adj_off[N] + 8));
    if (!comp || !queue || !dist || !dmin || !mass || !heap) {
        free(comp); free(queue); free(dist); free(dmin); free(mass); free(heap); return -1;
    }
    int    *resolved = malloc(sizeof(int) * (size_t)N);
    if (resolved) for (int i = 0; i < N; i++) resolved[i] = -1;
    int largest;
    int lc = cen_components(g, comp, queue, &largest);
    /* seeds are deep-clearance maxima: selected on radius, not mass */
    int ns = cen_seeds(g, ma->radius, comp, lc, K, beta, seeds, dmin, dist, heap);
    for (int s = 0; s < ns; s++) {
        double Es;
        results[s] = cen_descent(g, ma->verts.p, mass, comp, lc,
                                 seeds[s], polish_hops, heap, resolved, &Es);
    }
    free(comp); free(queue); free(dist); free(dmin); free(mass); free(heap); free(resolved);
    return ns;
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
