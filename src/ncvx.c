#include "vdiagram.h"
#include "vor3d.h"
#include "delaunay.h"
#include "scplx.h"
#include "dynarray.h"
#include "trimesh.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>


/* Full definition is private to this file; callers hold s_ncvx_tet only via pointer. */
struct ncvx_tet {
    s_point v[4];
    s_point aabb_min;
    s_point aabb_max;
};


/* -----------------------------------------------------------------------
 * Internal validity / free helpers
 * ----------------------------------------------------------------------- */

int ncvx_domain_is_valid(const s_ncvx_domain *d)
{
    return trimesh_is_valid(&d->surface) &&
           d->tets    != NULL &&
           d->N_tets  >  0   &&
           convhull_is_valid(&d->bpoly.convh);
}

void free_ncvx_domain(s_ncvx_domain *d)
{
    free_trimesh(&d->surface);
    free(d->tets);
    free_bpoly(&d->bpoly);
    memset(d, 0, sizeof(s_ncvx_domain));
}

void free_ncvx_vdiagram(s_ncvx_vdiagram *vd)
{
    int N = vd->seeds.N;   /* read before free_points zeroes the struct */
    free_points(&vd->seeds);
    free_ncvx_domain(&vd->domain);
    if (vd->vcells) {
        for (int i = 0; i < N; i++) free_vcell(&vd->vcells[i]);
        free(vd->vcells);
    }
    memset(vd, 0, sizeof(s_ncvx_vdiagram));
}


/* ---------------------------------------------------------------------- */

s_ncvx_domain ncvx_domain_from_trimesh(const s_trimesh *mesh,
                                        double EPS_DEG, double TOL)
{
    if (!trimesh_is_valid(mesh))
        return (s_ncvx_domain){0};

    /* Build CDT of the trimesh interior; returns only interior tetrahedra. */
    s_scplx dt = tetrahedralize_interior_trimesh(mesh, EPS_DEG, TOL);
    if (!dt.head) {
        fprintf(stderr, "ncvx_domain_from_trimesh: tetrahedralize_interior_trimesh failed\n");
        return (s_ncvx_domain){0};
    }

    /* Extract tets into the flat array used for Voronoi clipping. */
    s_dynarray tets = dynarray_initialize(sizeof(struct ncvx_tet), 256);
    if (!tets.items) { free_complex(&dt); return (s_ncvx_domain){0}; }

    double volume_sum = 0.0;

    for (s_ncell *nc = dt.head; nc; nc = nc->next) {
        s_point v[4];
        extract_vertices_ncell(&dt, nc, v);

        struct ncvx_tet tet;
        tet.v[0] = v[0]; tet.v[1] = v[1]; tet.v[2] = v[2]; tet.v[3] = v[3];

        tet.aabb_min.x = fmin(fmin(v[0].x, v[1].x), fmin(v[2].x, v[3].x));
        tet.aabb_min.y = fmin(fmin(v[0].y, v[1].y), fmin(v[2].y, v[3].y));
        tet.aabb_min.z = fmin(fmin(v[0].z, v[1].z), fmin(v[2].z, v[3].z));
        tet.aabb_max.x = fmax(fmax(v[0].x, v[1].x), fmax(v[2].x, v[3].x));
        tet.aabb_max.y = fmax(fmax(v[0].y, v[1].y), fmax(v[2].y, v[3].y));
        tet.aabb_max.z = fmax(fmax(v[0].z, v[1].z), fmax(v[2].z, v[3].z));

        volume_sum += fabs(signed_volume_tetra(v));

        if (!dynarray_push(&tets, &tet)) {
            dynarray_free(&tets); free_complex(&dt); return (s_ncvx_domain){0};
        }
    }

    free_complex(&dt);

    s_ncvx_domain domain = {0};

    domain.surface = copy_trimesh(mesh);
    if (!trimesh_is_valid(&domain.surface)) {
        dynarray_free(&tets); return (s_ncvx_domain){0};
    }

    domain.N_tets = (int)tets.N;
    if (domain.N_tets == 0) {
        fprintf(stderr, "ncvx_domain_from_trimesh: no interior tetrahedra found.\n");
        free_trimesh(&domain.surface); dynarray_free(&tets); return (s_ncvx_domain){0};
    }

    domain.tets = malloc(sizeof(struct ncvx_tet) * (size_t)domain.N_tets);
    if (!domain.tets) {
        free_trimesh(&domain.surface); dynarray_free(&tets); return (s_ncvx_domain){0};
    }
    memcpy(domain.tets, tets.items, sizeof(struct ncvx_tet) * (size_t)domain.N_tets);
    dynarray_free(&tets);

    domain.domain_volume = volume_sum;

    domain.bpoly = bpoly_from_points(&mesh->points, EPS_DEG);
    if (!convhull_is_valid(&domain.bpoly.convh)) {
        free_trimesh(&domain.surface);
        free(domain.tets);
        return (s_ncvx_domain){0};
    }

    domain.spatial_index = NULL;

    return domain;
}


/* ---------------------------------------------------------------------- */

static void candidate_tets(const s_ncvx_domain *domain,
                            s_point vcell_min, s_point vcell_max,
                            s_dynarray *out)
{
    for (int i = 0; i < domain->N_tets; i++) {
        const struct ncvx_tet *t = &domain->tets[i];
        if (vcell_max.x < t->aabb_min.x || t->aabb_max.x < vcell_min.x) continue;
        if (vcell_max.y < t->aabb_min.y || t->aabb_max.y < vcell_min.y) continue;
        if (vcell_max.z < t->aabb_min.z || t->aabb_max.z < vcell_min.z) continue;
        dynarray_push(out, &t);
    }
}

/* Build the outward-facing half-space plane for face fi of tet v[4]. */
static void tet_face_outward_plane(const s_point v[4], int fi, s_point plane[3])
{
    static const int idx[4][3] = {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
    s_point pa = v[idx[fi][0]], pb = v[idx[fi][1]], pc = v[idx[fi][2]], pd = v[fi];
    s_point n = cross_prod(subtract_points(pb, pa), subtract_points(pc, pa));
    if (dot_prod(subtract_points(pd, pa), n) > 0.0) {
        plane[0] = pa; plane[1] = pc; plane[2] = pb;  /* flip to point outward */
    } else {
        plane[0] = pa; plane[1] = pb; plane[2] = pc;
    }
}

static s_vcell clip_vcell_to_domain(const s_vcell *vcell,
                                     const s_ncvx_domain *domain,
                                     double EPS_DEG, double TOL)
{
    /* AABB of the vcell (it has exactly one piece on entry — it's a convex cell) */
    s_point vcell_min, vcell_max;
    bounding_box_points(&vcell->pieces[0].points, &vcell_min, &vcell_max);

    s_dynarray candidates = dynarray_initialize(sizeof(struct ncvx_tet*), 32);
    if (!candidates.items) return (s_vcell){0};
    candidate_tets(domain, vcell_min, vcell_max, &candidates);

    s_dynarray pieces = dynarray_initialize(sizeof(s_convh), 16);
    if (!pieces.items) { dynarray_free(&candidates); return (s_vcell){0}; }

    double volume = 0.0;

    for (unsigned int ci = 0; ci < candidates.N; ci++) {
        const struct ncvx_tet *t;
        dynarray_get_value(&candidates, ci, &t);

        s_convh piece = copy_convhull(&vcell->pieces[0]);
        bool ok = true;
        for (int fi = 0; fi < 4 && ok; fi++) {
            s_point plane[3];
            tet_face_outward_plane(t->v, fi, plane);
            s_convh clipped;
            int r = clip_convhull_halfspace(&piece, plane, EPS_DEG, TOL, &clipped);
            free_convhull(&piece);
            if (r != 1) { piece = convhull_NAN; ok = false; }
            else         { piece = clipped; }
        }

        if (!ok || !convhull_is_valid(&piece)) { free_convhull(&piece); continue; }

        volume += volume_convhull(&piece);
        if (!dynarray_push(&pieces, &piece)) {
            free_convhull(&piece);
            dynarray_free(&pieces);
            dynarray_free(&candidates);
            return (s_vcell){0};
        }
    }
    dynarray_free(&candidates);

    s_vcell result = {0};
    result.seed_id  = vcell->seed_id;
    result.N_pieces = (int)pieces.N;
    result.volume   = volume;

    if (pieces.N > 0) {
        result.pieces = malloc(sizeof(s_convh) * pieces.N);
        if (!result.pieces) {
            for (unsigned int i = 0; i < pieces.N; i++) {
                s_convh p; dynarray_get_value(&pieces, i, &p); free_convhull(&p);
            }
            dynarray_free(&pieces);
            return (s_vcell){0};
        }
        memcpy(result.pieces, pieces.items, sizeof(s_convh) * pieces.N);
    }
    dynarray_free(&pieces);

    return result;
}


/* ---------------------------------------------------------------------- */

s_ncvx_vdiagram vor3d_in_ncvx_domain(const s_points *seeds,
                                  const s_ncvx_domain *domain,
                                  double vol_max_rel_diff,
                                  double EPS_DEG, double TOL,
                                  int (*randint)(void*, int), void *rctx,
                                  s_dynarray *buff_points, int *out_kept_idx)
{
    if (!seeds || seeds->N <= 0 || !ncvx_domain_is_valid(domain))
        return (s_ncvx_vdiagram){0};

    /* Stage B: convex Voronoi on the convex hull of the trimesh */
    s_vdiagram vd = vor3d_in_bp(seeds, &domain->bpoly, vol_max_rel_diff,
                                 EPS_DEG, TOL, randint, rctx,
                                 buff_points, out_kept_idx);
    if (vd.seeds.N == 0) return (s_ncvx_vdiagram){0};

    /* Stage C: clip each convex cell to the domain interior */
    int N = vd.seeds.N;
    s_vcell *vcells = calloc((size_t)N, sizeof(s_vcell));
    if (!vcells) { free_vdiagram(&vd); return (s_ncvx_vdiagram){0}; }

    for (int i = 0; i < N; i++)
        vcells[i] = clip_vcell_to_domain(&vd.vcells[i], domain, EPS_DEG, TOL);

    /* Assemble result — copies seeds and domain */
    s_ncvx_vdiagram out;
    out.seeds = copy_points(&vd.seeds);
    out.domain = (s_ncvx_domain){
        .surface       = copy_trimesh(&domain->surface),
        .tets          = malloc(sizeof(struct ncvx_tet) * (size_t)domain->N_tets),
        .N_tets        = domain->N_tets,
        .domain_volume = domain->domain_volume,
        .bpoly         = bpoly_copy(&domain->bpoly),
        .spatial_index = NULL,
    };

    if (!out.seeds.p || !out.domain.tets ||
        !trimesh_is_valid(&out.domain.surface)) {
        free_points(&out.seeds);
        free_ncvx_domain(&out.domain);
        for (int i = 0; i < N; i++) free_vcell(&vcells[i]);
        free(vcells);
        free_vdiagram(&vd);
        return (s_ncvx_vdiagram){0};
    }
    memcpy(out.domain.tets, domain->tets,
           sizeof(struct ncvx_tet) * (size_t)domain->N_tets);

    out.vcells = vcells;
    free_vdiagram(&vd);
    return out;
}


s_ncvx_vdiagram vor3d_in_trimesh(const s_points *seeds,
                                   const s_trimesh *mesh,
                                   double vol_max_rel_diff,
                                   double EPS_DEG, double TOL,
                                   int (*randint)(void*, int), void *rctx,
                                   s_dynarray *buff_points, int *out_kept_idx)
{
    s_ncvx_domain domain = ncvx_domain_from_trimesh(mesh, EPS_DEG, TOL);
    if (!ncvx_domain_is_valid(&domain)) return (s_ncvx_vdiagram){0};
    s_ncvx_vdiagram out = vor3d_in_ncvx_domain(seeds, &domain, vol_max_rel_diff,
                                           EPS_DEG, TOL, randint, rctx,
                                           buff_points, out_kept_idx);
    free_ncvx_domain(&domain);  /* result owns its own copy */
    return out;
}
