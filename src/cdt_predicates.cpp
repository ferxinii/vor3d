#include "cdt_predicates.h"
#include "indirect_predicates.h"
#include <deque>
#include <vector>

/* explicitPoint3D objects live in a deque: stable references after push_back/resize. */
static std::deque<explicitPoint3D>       g_explicit;
/* g_lnc[id] is non-null iff point id is a Steiner; heap-allocated to hold
 * stable references into g_explicit. */
static std::vector<implicitPoint3D_LNC*> g_lnc;
/* g_pts[id] points to either the explicit or LNC object. */
static std::vector<genericPoint*>        g_pts;

static void ensure_capacity(int id)
{
    if ((int)g_explicit.size() <= id) g_explicit.resize(id + 1);
    if ((int)g_lnc.size()     <= id) g_lnc.resize(id + 1, nullptr);
    if ((int)g_pts.size()     <= id) g_pts.resize(id + 1, nullptr);
}

extern "C" {

void cdt_predicates_init()
{
    initFPU();
}

void cdt_predicates_clear()
{
    for (auto* p : g_lnc) delete p;
    g_explicit.clear();
    g_lnc.clear();
    g_pts.clear();
}

void cdt_point_set_explicit(int id, double x, double y, double z)
{
    ensure_capacity(id);
    g_explicit[id].set(x, y, z);
    if (g_lnc[id]) { delete g_lnc[id]; g_lnc[id] = nullptr; }
    g_pts[id] = &g_explicit[id];
}

void cdt_point_set_lnc(int id, int v1, int v2, double t)
{
    /* paper: t*V[v1] + (1-t)*V[v2]
     * library: (1-t_lib)*P + t_lib*Q  ->  P=V[v1], Q=V[v2], t_lib = 1.0-t */
    ensure_capacity(id);
    if (g_lnc[id]) delete g_lnc[id];
    g_lnc[id] = new implicitPoint3D_LNC(g_explicit[v1], g_explicit[v2], 1.0 - t);
    g_pts[id] = g_lnc[id];

    /* Also store this point's rounded coords as an explicit entry.  The library
     * LNC combines two EXPLICIT points, so a nested Steiner -- one placed on a
     * sub-segment whose endpoint is a prior Steiner -- must reference that
     * Steiner through g_explicit.  Compute with the SAME t*V[v1]+(1-t)*V[v2]
     * formula the builder uses for M (points.p), so g_explicit[id] == M exactly
     * and nested LNCs stay consistent with the rounded geometry.  (g_pts[id]
     * remains the exact LNC for direct use of this point.) */
    g_explicit[id].set(t*g_explicit[v1].X() + (1.0-t)*g_explicit[v2].X(),
                       t*g_explicit[v1].Y() + (1.0-t)*g_explicit[v2].Y(),
                       t*g_explicit[v1].Z() + (1.0-t)*g_explicit[v2].Z());
}

int cdt_orient3d(int a, int b, int c, int d)
{
    return genericPoint::orient3D(*g_pts[a], *g_pts[b], *g_pts[c], *g_pts[d]);
}

int cdt_insphere(int a, int b, int c, int d, int e)
{
    return genericPoint::inSphere(*g_pts[a], *g_pts[b], *g_pts[c], *g_pts[d], *g_pts[e]);
}

int cdt_perturbed_insphere(int i1, int i2, int i3, int i4, int i5)
{
    int r = genericPoint::inSphere(*g_pts[i1], *g_pts[i2], *g_pts[i3],
                                   *g_pts[i4], *g_pts[i5]);
    if (r != 0) return r;

    /* SoS: sort [i1..i5] ascending by insertion sort; track swap parity n. */
    int ids[5] = {i1, i2, i3, i4, i5};
    int n = 0;
    for (int i = 1; i < 5; i++) {
        int key = ids[i], j = i - 1;
        while (j >= 0 && ids[j] > key) { ids[j+1] = ids[j]; j--; n ^= 1; }
        ids[j+1] = key;
    }

    /* Try cofactor orient3D calls, skipping each row in turn.
     * Sign for skip k is (-1)^(k+n). */
    for (int skip = 0; skip < 5; skip++) {
        int sub[4], k = 0;
        for (int j = 0; j < 5; j++) if (j != skip) sub[k++] = ids[j];
        int o = genericPoint::orient3D(*g_pts[sub[0]], *g_pts[sub[1]],
                                       *g_pts[sub[2]], *g_pts[sub[3]]);
        if (o != 0) {
            return ((skip ^ n) & 1) ? -o : o;
        }
    }
    return 0; /* all five coincident -- degenerate input */
}

int cdt_point_in_tet(int p, int t0, int t1, int t2, int t3)
{
    /* Algorithm 3 from paper: p is inside iff all four orient3D tests are >= 0. */
    if (genericPoint::orient3D(*g_pts[t0], *g_pts[t1], *g_pts[t2], *g_pts[p])  < 0) return 0;
    if (genericPoint::orient3D(*g_pts[t0], *g_pts[t1], *g_pts[p],  *g_pts[t3]) < 0) return 0;
    if (genericPoint::orient3D(*g_pts[t0], *g_pts[p],  *g_pts[t2], *g_pts[t3]) < 0) return 0;
    if (genericPoint::orient3D(*g_pts[p],  *g_pts[t1], *g_pts[t2], *g_pts[t3]) < 0) return 0;
    return 1;
}

int cdt_segment_crosses_triangle(int s1, int s2, int v0, int v1, int v2)
{
    /* Algorithm 4 from paper: interior of s1-s2 crosses interior of <v0,v1,v2>. */
    int o1 = genericPoint::orient3D(*g_pts[v0], *g_pts[v1], *g_pts[v2], *g_pts[s1]);
    int o2 = genericPoint::orient3D(*g_pts[v0], *g_pts[v1], *g_pts[v2], *g_pts[s2]);
    if (o1 == o2) return 0; /* both on same side (or both on plane) */
    int a = genericPoint::orient3D(*g_pts[v0], *g_pts[v1], *g_pts[s1], *g_pts[s2]);
    int b = genericPoint::orient3D(*g_pts[v1], *g_pts[v2], *g_pts[s1], *g_pts[s2]);
    int c = genericPoint::orient3D(*g_pts[v2], *g_pts[v0], *g_pts[s1], *g_pts[s2]);
    if (a * b < 0 || b * c < 0) return 0;
    return 1;
}

int cdt_inner_seg_crosses_inner_tri(int s1, int s2, int v0, int v1, int v2)
{
    return genericPoint::innerSegmentCrossesInnerTriangle(
        *g_pts[s1], *g_pts[s2], *g_pts[v0], *g_pts[v1], *g_pts[v2]) ? 1 : 0;
}

int cdt_inner_segs_cross(int a, int b, int p, int q)
{
    return genericPoint::innerSegmentsCross(
        *g_pts[a], *g_pts[b], *g_pts[p], *g_pts[q]) ? 1 : 0;
}

int cdt_point_in_inner_tri(int p, int a, int b, int c)
{
    return genericPoint::pointInInnerTriangle(
        *g_pts[p], *g_pts[a], *g_pts[b], *g_pts[c]) ? 1 : 0;
}

int cdt_point_in_triangle(int p, int a, int b, int c)
{
    /* OUT if not coplanar with the triangle plane. */
    if (genericPoint::orient3D(*g_pts[a], *g_pts[b], *g_pts[c], *g_pts[p]) != 0) return -1;
    /* Coplanar: pointInTriangle returns in-or-on plus the three edge signs. */
    int o1, o2, o3;
    bool on = genericPoint::pointInTriangle(*g_pts[p], *g_pts[a], *g_pts[b], *g_pts[c], o1, o2, o3);
    if (!on) return -1;
    return (o1 == 0 || o2 == 0 || o3 == 0) ? 0 : 1;
}

int cdt_segments_cross(int s1, int s2, int p, int q)
{
    return genericPoint::segmentsCross(
        *g_pts[s1], *g_pts[s2], *g_pts[p], *g_pts[q]) ? 1 : 0;
}

int cdt_get_approx_coords(int id, double *x, double *y, double *z)
{
    if (id < 0 || id >= (int)g_pts.size() || !g_pts[id]) return 0;
    return g_pts[id]->getApproxXYZCoordinates(*x, *y, *z) ? 1 : 0;
}

} /* extern "C" */
