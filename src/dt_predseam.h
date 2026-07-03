#ifndef VOR3D_DT_PREDSEAM_H
#define VOR3D_DT_PREDSEAM_H

/* Library-internal predicate seam for the DT builder (Phase 1 of CDT_REFACTOR).
 *
 * Every topological decision in the builder goes through orient3d / insphere.
 * These two helpers are the single point where that computation is selected:
 *   - default (s->exact_ids == 0): the CURRENT coordinate predicates, i.e. a
 *     verbatim reproduction of the calls the builder made before the seam
 *     existed.  This is the only path taken in Phase 1; nothing observable
 *     changes.
 *   - exact  (s->exact_ids == 1): the id-keyed indirect predicates from
 *     cdt_predicates (activated only for CDT builds, Phase 2).
 *
 * Both take scplx point INDICES (a,b,c,d,e); the exact path feeds them
 * straight to cdt_*, the default path dereferences s->points.p[].  Kept as
 * static inline in a private header so dt.c and scplx.c share one definition
 * (the point-location walk in scplx.c needs the same seam) without an extra
 * translation unit; unused-in-a-TU does not warn for static inline.
 *
 * This header is PRIVATE to the library -- do not install / expose it.
 *
 * SIGN CONVENTIONS (load-bearing; VERIFIED by tests/predseam_convention.c):
 *   dtp_orient   returns the sign of orient3d(a,b,c,d) in Shewchuk's
 *                test_orientation convention (the builder's existing one).
 *   dtp_insphere returns -1 when e is OUTSIDE circumsphere(a,b,c,d), +1 when
 *                INSIDE, 0 when ON -- i.e. test_insphere's convention, which
 *                are_locally_delaunay reads as "in1 == -1 => locally Delaunay".
 *
 * Two convention facts, both measured (not assumed) by the unit test, drive
 * the exact branch below:
 *   (1) genericPoint::orient3D == -Shewchuk orient3d (OPPOSITE sign).  So the
 *       exact dtp_orient NEGATES cdt_orient3d to match test_orientation.
 *   (2) test_insphere is orientation-dependent relative to the cdt insphere:
 *         test_insphere(a,b,c,d,e) == sign(cdt_orient3d(a,b,c,d)) * cdt_insphere(a,b,c,d,e)
 *       (test_insphere == +cdt_insphere for cdt-positive tets, -cdt_insphere
 *       for cdt-negative ones).  A plain negation would be correct ONLY if the
 *       builder guaranteed a fixed tet orientation; it does not, so the exact
 *       dtp_insphere folds in the orientation sign explicitly.  It uses the
 *       SoS (never-0) perturbed insphere, so on cocircular inputs it returns a
 *       consistent nonzero sign instead of 0.
 *
 * NOTE (Appendix A / determine_case): the exact determine_case mirror uses
 * cdt_orient3d DIRECTLY and returns a CASE enum from RELATIVE sign comparisons
 * only, which are invariant under the global sign flip in (1) -- so it is
 * correct regardless of the orient parity and does NOT need the negation here.
 * Sites that consume an orientation sign at face value (the walk, flip44 gate,
 * are_locally_delaunay flat-check, the p-on-edge test) DO go through this seam
 * and get the test_orientation convention.
 */

#include "scplx.h"          /* s_scplx, s_point (via points.h) */
#include "gtests.h"         /* test_orientation, test_insphere */
#include "cdt_predicates.h" /* cdt_orient3d, cdt_perturbed_insphere */

/* Translate a scplx-local vertex id to a cdt_predicates registry id.  Identity
 * for the global CDT DT (l2g_ids == NULL); for Phase B local cavity DTs it maps
 * local ids (sentinels 0..3 -> scratch registry slots, real 4+i -> global id).*/
static inline int dtp_g(const s_scplx *s, int id)
{
    return s->l2g_ids ? s->l2g_ids[id] : id;
}

/* orient3d sign of d vs plane(a,b,c).  Args are scplx point indices. */
static inline int dtp_orient(const s_scplx *s, int a, int b, int c, int d)
{
    if (s->exact_ids)   /* genericPoint::orient3D == -Shewchuk orient3d */
        return -cdt_orient3d(dtp_g(s,a), dtp_g(s,b), dtp_g(s,c), dtp_g(s,d));
    const s_point *p = s->points.p;
    return test_orientation((s_point[3]){ p[a], p[b], p[c] }, p[d]);
}

/* insphere of e vs circumsphere(a,b,c,d).  Args are scplx point indices.
 * -1 outside, +1 inside, 0 on (default path only; exact path never returns 0). */
static inline int dtp_insphere(const s_scplx *s, int a, int b, int c, int d, int e)
{
    if (s->exact_ids) {
        /* test_insphere == sign(cdt_orient3d(a,b,c,d)) * cdt_insphere(a,b,c,d,e).
         * Use the SoS (never-0) perturbed insphere; fold in orientation sign. */
        int ga=dtp_g(s,a), gb=dtp_g(s,b), gc=dtp_g(s,c), gd=dtp_g(s,d), ge=dtp_g(s,e);
        int so = cdt_orient3d(ga, gb, gc, gd);
        int r  = cdt_perturbed_insphere(ga, gb, gc, gd, ge);
        return (so < 0) ? -r : r;   /* so == 0 only on a flat tet (excluded upstream) */
    }
    const s_point *p = s->points.p;
    return test_insphere((s_point[4]){ p[a], p[b], p[c], p[d] }, p[e]);
}

/* Exact point-in-triangle (-1 OUT, 0 BND, +1 IN); ids are scplx-local. */
static inline int dtp_point_in_triangle(const s_scplx *s, int p, int a, int b, int c)
{
    return cdt_point_in_triangle(dtp_g(s,p), dtp_g(s,a), dtp_g(s,b), dtp_g(s,c));
}

/* Exact closed-segment cross of coplanar s1-s2 and p-q; ids are scplx-local. */
static inline int dtp_segments_cross(const s_scplx *s, int s1, int s2, int p, int q)
{
    return cdt_segments_cross(dtp_g(s,s1), dtp_g(s,s2), dtp_g(s,p), dtp_g(s,q));
}

#endif /* VOR3D_DT_PREDSEAM_H */
