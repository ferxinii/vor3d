/* Domain construction and queries for the unified vor3d API (vor3d.h).
 * The diagram entry point vor3d_diagram lives in ncvx.c next to both clip
 * kernels. */
#include "vor3d.h"
#include "vdiagram.h"
#include "delaunay.h"
#include "scplx.h"
#include "convh.h"
#include "trimesh.h"
#include "points.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>


bool vor3d_domain_is_valid(const s_vor3d_domain *d)
{
    if (!d || !trimesh_is_valid(&d->surface) || !convhull_is_valid(&d->bpoly.convh))
        return false;
    return d->is_convex || d->cdt.head != NULL;
}

void free_vor3d_domain(s_vor3d_domain *d)
{
    free_trimesh(&d->surface);
    free_complex(&d->cdt);
    free_bpoly(&d->bpoly);
    memset(d, 0, sizeof(s_vor3d_domain));
}


/* Convex domain from an owned bpoly: the surface is the hull's own triangles. */
static s_vor3d_domain domain_from_owned_bpoly(s_bpoly bp)
{
    if (!convhull_is_valid(&bp.convh)) { free_bpoly(&bp); return (s_vor3d_domain){0}; }
    s_trimesh surface = trimesh_from_arrays(bp.convh.points.p, bp.convh.points.N,
                                            bp.convh.faces, bp.convh.Nf, 0.0);
    if (!trimesh_is_valid(&surface)) {
        fprintf(stderr, "vor3d_domain: hull is not a closed manifold\n");
        free_bpoly(&bp);
        return (s_vor3d_domain){0};
    }
    return (s_vor3d_domain){
        .is_convex = true,
        .surface   = surface,
        .bpoly     = bp,
        .cdt       = (s_scplx){0},
        .volume    = bp.volume,
    };
}

s_vor3d_domain vor3d_domain_from_convh(const s_convh *convh)
{
    if (!convh || !convhull_is_valid(convh)) return (s_vor3d_domain){0};
    return domain_from_owned_bpoly(bpoly_from_convh(convh));
}

s_vor3d_domain vor3d_domain_from_points(const s_points *points, double EPS_DEG)
{
    if (!points || points->N < 4) return (s_vor3d_domain){0};
    return domain_from_owned_bpoly(bpoly_from_points(points, EPS_DEG));
}


s_vor3d_domain vor3d_domain_from_trimesh(const s_trimesh *mesh, double EPS_DEG,
                                         double TOL, int keep_exterior)
{
    if (!mesh || !trimesh_is_valid(mesh)) return (s_vor3d_domain){0};

    /* The hull is needed by both kernels: it IS the domain when convex, and
     * sizes the far auxiliary box / fast-rejects containment queries when not. */
    s_bpoly bp = bpoly_from_points(&mesh->points, EPS_DEG);
    if (!convhull_is_valid(&bp.convh)) {
        free_bpoly(&bp);
        return (s_vor3d_domain){0};
    }

    if (trimesh_is_convex(mesh)) {
        s_trimesh surface = copy_trimesh(mesh);
        if (!trimesh_is_valid(&surface)) { free_bpoly(&bp); return (s_vor3d_domain){0}; }
        return (s_vor3d_domain){
            .is_convex = true,
            .surface   = surface,
            .bpoly     = bp,
            .cdt       = (s_scplx){0},
            .volume    = bp.volume,
        };
    }

    /* Non-convex: CDT of the interior. keep_exterior=0 returns only interior
     * tets; keep_exterior=1 returns the full convex-hull complex (interior +
     * exterior pockets flagged via nc->interior), which the clip treats
     * identically and which also supports exact walk-based point location. */
    s_scplx cdt = keep_exterior ? tetrahedralize_domain_flagged(mesh, EPS_DEG, TOL)
                                : tetrahedralize_interior_trimesh(mesh, EPS_DEG, TOL);
    if (!cdt.head) {
        fprintf(stderr, "vor3d_domain_from_trimesh: tetrahedralization failed\n");
        free_bpoly(&bp);
        return (s_vor3d_domain){0};
    }

    double volume = 0.0;
    for (s_ncell *nc = cdt.head; nc; nc = nc->next) {
        if (!nc->interior) continue;              /* skip exterior pockets, if kept */
        s_point v[4];
        extract_vertices_ncell(&cdt, nc, v);
        volume += fabs(signed_volume_tetra(v));
    }

    s_trimesh surface = copy_trimesh(mesh);
    if (!trimesh_is_valid(&surface)) {
        free_complex(&cdt); free_bpoly(&bp);
        return (s_vor3d_domain){0};
    }
    return (s_vor3d_domain){
        .is_convex        = false,
        .surface          = surface,
        .bpoly            = bp,
        .cdt              = cdt,
        .cdt_has_exterior = keep_exterior != 0,
        .volume           = volume,
    };
}


bool vor3d_domain_contains(const s_vor3d_domain *d, s_point p, double EPS_DEG, double TOL)
{
    const s_bpoly *bp = &d->bpoly;
    for (int a = 0; a < 3; a++)
        if (p.coords[a] < bp->min.coords[a] || p.coords[a] > bp->max.coords[a]) return false;
    e_geom_test t = test_point_in_convhull(&bp->convh, p, EPS_DEG, TOL);
    if (t == TEST_OUT) return false;
    if (d->is_convex) return t == TEST_IN || t == TEST_BOUNDARY;

    /* Non-convex: a walk that lands in a tet is conclusive (exact orientation
     * tests). With exterior pockets kept the complex is the whole hull, so a
     * walk that steps off it is conclusive too. The lean (interior-only)
     * complex is not convex, so there a walk can leave through a concavity
     * while p is inside: fall back to the winding number. */
    const s_ncell *nc = in_ncell_walk(&d->cdt, p);
    if (nc) return nc->interior;
    if (d->cdt_has_exterior) return false;
    return point_in_trimesh_winding(&d->surface, p) == 1;
}
