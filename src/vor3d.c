/* Domain construction and queries for the unified vor3d API (vor3d.h).
 * The diagram entry point vor3d_diagram lives in ncvx.c next to both clip
 * kernels. */
#include "vor3d.h"
#include "dt.h"
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
    if (!d || !trimesh_is_valid(&d->surface) || !convhull_is_valid(&d->convh))
        return false;
    return d->is_convex || (d->cdt && d->cdt->head != NULL);
}

void free_vor3d_domain(s_vor3d_domain *d)
{
    free_trimesh(&d->surface);
    if (d->cdt) { free_complex(d->cdt); free(d->cdt); }
    free_convhull(&d->convh);
    memset(d, 0, sizeof(s_vor3d_domain));
}


/* Domain with only the hull fields set (takes ownership of convh). */
static s_vor3d_domain domain_with_convh(s_convh convh)
{
    s_vor3d_domain d = { .convh = convh };
    d.convh_vertex_avg = point_average(&convh.points);
    s_point mn = convh.points.p[0], mx = convh.points.p[0];
    for (int i = 1; i < convh.points.N; i++) {
        s_point q = convh.points.p[i];
        if (q.x < mn.x) mn.x = q.x;
        if (q.y < mn.y) mn.y = q.y;
        if (q.z < mn.z) mn.z = q.z;
        if (q.x > mx.x) mx.x = q.x;
        if (q.y > mx.y) mx.y = q.y;
        if (q.z > mx.z) mx.z = q.z;
    }
    d.bbox_min = mn;
    d.bbox_max = mx;
    d.volume = volume_convhull(&convh);
    return d;
}

/* Convex domain from an owned hull: the surface is the hull's own triangles. */
static s_vor3d_domain domain_from_owned_convh(s_convh convh)
{
    if (!convhull_is_valid(&convh)) { free_convhull(&convh); return (s_vor3d_domain){0}; }
    s_vor3d_domain d = domain_with_convh(convh);
    d.surface = trimesh_from_arrays(convh.points.p, convh.points.N, convh.faces, convh.Nf);
    if (!trimesh_is_valid(&d.surface)) {
        fprintf(stderr, "vor3d_domain: hull is not a closed manifold\n");
        free_vor3d_domain(&d);
        return (s_vor3d_domain){0};
    }
    d.is_convex = true;
    return d;
}

s_vor3d_domain vor3d_domain_from_convh(const s_convh *convh)
{
    if (!convh || !convhull_is_valid(convh)) return (s_vor3d_domain){0};
    return domain_from_owned_convh(copy_convhull(convh));
}

s_vor3d_domain vor3d_domain_from_points(const s_points *points)
{
    if (!points || points->N < 4) return (s_vor3d_domain){0};
    s_convh convh;
    if (convhull_from_points(points, &convh) != 1) return (s_vor3d_domain){0};
    return domain_from_owned_convh(convh);
}


s_vor3d_domain vor3d_domain_from_trimesh(const s_trimesh *mesh, const s_vor3d_opts *opts,
                                         int keep_exterior)
{
    const double TOL = opts ? opts->tol_dup : 0.0;
    if (!mesh || !trimesh_is_valid(mesh)) return (s_vor3d_domain){0};

    /* The hull is needed by both kernels: it IS the domain when convex, and
     * sizes the far auxiliary box / fast-rejects containment queries when not. */
    s_convh convh;
    if (convhull_from_points(&mesh->points, &convh) != 1) return (s_vor3d_domain){0};
    s_vor3d_domain d = domain_with_convh(convh);

    if (trimesh_is_convex(mesh)) {
        d.surface = copy_trimesh(mesh);
        if (!trimesh_is_valid(&d.surface)) { free_vor3d_domain(&d); return (s_vor3d_domain){0}; }
        d.is_convex = true;
        return d;
    }

    /* Non-convex: CDT of the interior. keep_exterior=0 returns only interior
     * tets; keep_exterior=1 returns the full convex-hull complex (interior +
     * exterior pockets flagged via nc->interior), which the clip treats
     * identically and which also supports exact walk-based point location. */
    d.cdt = malloc(sizeof *d.cdt);
    if (d.cdt) *d.cdt = keep_exterior ? tetrahedralize_domain_flagged(mesh, TOL)
                                      : tetrahedralize_interior_trimesh(mesh, TOL);
    if (!d.cdt || !d.cdt->head) {
        fprintf(stderr, "vor3d_domain_from_trimesh: tetrahedralization failed\n");
        free_vor3d_domain(&d);
        return (s_vor3d_domain){0};
    }
    d.cdt_has_exterior = keep_exterior != 0;

    double volume = 0.0;
    for (s_ncell *nc = d.cdt->head; nc; nc = nc->next) {
        if (!nc->interior) continue;              /* skip exterior pockets, if kept */
        s_point v[4];
        extract_vertices_ncell(d.cdt, nc, v);
        volume += fabs(signed_volume_tetra(v));
    }
    d.volume = volume;

    d.surface = copy_trimesh(mesh);
    if (!trimesh_is_valid(&d.surface)) { free_vor3d_domain(&d); return (s_vor3d_domain){0}; }
    return d;
}


bool vor3d_domain_contains(const s_vor3d_domain *d, s_point p, double tol_bdry)
{
    for (int a = 0; a < 3; a++)
        if (p.coords[a] < d->bbox_min.coords[a] || p.coords[a] > d->bbox_max.coords[a]) return false;
    e_geom_test t = test_point_in_convhull(&d->convh, p, tol_bdry);
    if (t == TEST_OUT) return false;
    if (d->is_convex) return t == TEST_IN || t == TEST_BOUNDARY;

    /* Non-convex: a walk that lands in a tet is conclusive (exact orientation
     * tests). With exterior pockets kept the complex is the whole hull, so a
     * walk that steps off it is conclusive too. The lean (interior-only)
     * complex is not convex, so there a walk can leave through a concavity
     * while p is inside: fall back to the winding number. */
    const s_ncell *nc = in_ncell_walk(d->cdt, p);
    if (nc) return nc->interior;
    if (d->cdt_has_exterior) return false;
    return point_in_trimesh_winding(&d->surface, p) == 1;
}
