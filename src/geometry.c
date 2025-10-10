
#include "geometry.h"
#include "algebra.h"
#include "external/geometric_predicates/include/predicates.h"
#include "external/convhull_3d/convhull_3d.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>


int orientation(const s_point *p3, s_point q)
{
    double aux = orient3d(p3[0].coords, p3[1].coords, p3[2].coords, q.coords);
    if (aux > 0) return 1;
    else if (aux < 0) return -1;
    else return 0;
}


int in_sphere(const s_point *p4, s_point q)
{   

    int factor;
    if (orientation(p4, p4[3]) == 1) factor = 1;
    else factor = -1;

    double aux = insphere(p4[0].coords, p4[1].coords, p4[2].coords, p4[3].coords, q.coords);
    
    if (aux > 0) return factor;
    else if (aux < 0) return -factor;
    else return 0;
}


double dot_prod(s_point u, s_point v)
{
    return u.x*v.x + u.y*v.y + u.z*v.z;
}


s_point cross_prod(s_point u, s_point v)
{
    s_point out;
    out.x = u.y*v.z - u.z*v.y;
    out.y = u.z*v.x - u.x*v.z;
    out.z = u.x*v.y - u.y*v.x;
    return out;
}


double norm_squared(s_point v)
{
    return dot_prod(v, v);
}


double norm(s_point v)
{
    return sqrt(norm_squared(v));
}


double distance_squared(s_point a, s_point b)
{
    s_point diff = subtract_points(a, b);
    return norm_squared(diff);
}


double distance(s_point a, s_point b)
{
    return sqrt(distance_squared(a, b));
}


double max_distance(const s_point *p, int N, s_point query)
{   
    double maxd2 = 0;
    for (int ii=0; ii<N; ii++) {
        double d2 = distance_squared(p[ii], query);
        if (maxd2 < d2) maxd2 = d2;
    }
    return sqrt(maxd2);
}


s_point normalize_3d(s_point v)
{
    double n = norm(v);
    return (s_point){{{v.x/n, v.y/n, v.z/n}}};
}


s_point subtract_points(s_point u, s_point v)
{
    return (s_point){{{u.x - v.x, u.y - v.y, u.z - v.z}}};
}


s_point find_center_mass(const s_point *in, int N_points)
{
    s_point out = in[0];
    for (int ii=1; ii<N_points; ii++) {
        out.x += in[ii].x;
        out.y += in[ii].y;
        out.z += in[ii].z;
    }
    out.x /= N_points;
    out.y /= N_points;
    out.z /= N_points;
    return out;
}


int coord_with_largest_component_3d(s_point x)
{
    double *n = x.coords;
    int out = 2;
    if (fabs(n[0]) > fabs(n[1]) && fabs(n[0]) > fabs(n[2])) out = 0;
    else if (fabs(n[1]) > fabs(n[0]) && fabs(n[1]) > fabs(n[2])) out = 1;
    return out;
}


int segment_crosses_triangle_3d(const s_point *triangle, s_point a, s_point b)
{
    if (orientation(triangle, a) == orientation(triangle, b)) return 0;
    s_point aux[3];
    aux[2] = a; 

    aux[0] = triangle[0];   aux[1] = triangle[1];
    int s1 = orientation(aux, b);

    aux[0] = triangle[1];   aux[1] = triangle[2];
    int s2 = orientation(aux, b);

    aux[0] = triangle[2];   aux[1] = triangle[0];
    int s3 = orientation(aux, b);
    
    if (s1 == s2 && s2 == s3 && s3 == s1) return 1;
    else return 0;
}


int between_1d(double x, double a, double b, double eps)
{
    // returns 1 if x lies between [a,b]
    if (a > b) { double t=a; a=b; b=t; }
    return (x + eps >= a && x - eps <= b);
}


int segments_intersect_2d(const s_point *AB, const s_point *pd)
{
    const double EPS = 1e-9;
    double Ax = AB[0].x, Ay = AB[0].y;
    double Bx = AB[1].x, By = AB[1].y;
    double px = pd[0].x, py = pd[0].y;
    double dx = pd[1].x, dy = pd[1].y;

    // 1) Degenerate pd, treat as X
    if (fabs(px - dx) < EPS && fabs(py - dy) < EPS) {
        double X[2] = { px, py };
        if (fabs(Ax - Bx) < EPS && fabs(Ay - By) < EPS) {
            // If A, B also degenerate
            return (fabs(px - Ax) < EPS && fabs(py - Ay) < EPS);
        }
        if (orient2d(AB[0].coords, AB[1].coords, X) == 0 &&
            between_1d(px, Ax, Bx, EPS) &&
            between_1d(py, Ay, By, EPS)) {
            return 1;
        }
        return 0;
    }

    // 2) Degenerate AB, treat as Y
    if (fabs(Ax - Bx) < EPS && fabs(Ay - By) < EPS) {
        double Y[2] = { Ax, Ay };
        if (orient2d(pd[0].coords, pd[1].coords, Y) == 0 &&
            between_1d(Ax, px, dx, EPS) &&
            between_1d(Ay, py, dy, EPS)) {
            return 1;
        }
        return 0;
    }

    // 3) General case: two straddling tests
    // 3a) [A,B] vs {p,d}
    int o1 = orient2d(AB[0].coords, AB[1].coords, pd[0].coords);
    int o2 = orient2d(AB[0].coords, AB[1].coords, pd[1].coords);
    if (o1 != 0 && o2 != 0 && o1 == o2) return 0;

    // 3b) [p,d] vs {A,B}
    o1 = orient2d(pd[0].coords, pd[1].coords, AB[0].coords);
    o2 = orient2d(pd[0].coords, pd[1].coords, AB[1].coords);
    if (o1 != 0 && o2 != 0 && o1 == o2) return 0;

    return 1;
}


s_point closest_point_on_triangle(const s_point *triangle, s_point p)
{
    s_point A = triangle[0], B = triangle[1], C = triangle[2];
    s_point AB = {{{B.x-A.x, B.y-A.y, B.z-A.z}}};
    s_point AC = {{{C.x-A.x, C.y-A.y, C.z-A.z}}};
    s_point AP = {{{p.x-A.x, p.y-A.y, p.z-A.z}}};

    // In vertex A
    double d1 = dot_prod(AB, AP), d2 = dot_prod(AC, AP);
    if (d1 <= 0 && d2 <= 0) return A;

    // In vertex B
    s_point BP = {{{ p.x-B.x, p.y-B.y, p.z-B.z }}};
    double d3 = dot_prod(AB, BP), d4 = dot_prod(AC, BP);
    if (d3 >= 0 && d4 <= d3) return B;

    // In edge AB
    double vc = d1*d4 - d3*d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) {
        double v = d1 / (d1 - d3);
        return (s_point){{{A.x + v*AB.x, A.y + v*AB.y, A.z + v*AB.z}}};
    }

    // In vertex C
    s_point CP = {{{ p.x-C.x, p.y-C.y, p.z-C.z }}};
    double d5 = dot_prod(AB, CP), d6 = dot_prod(AC, CP);
    if (d6 >= 0 && d5 <= d6) return C;

    // In edge AC
    double vb = d5*d2 - d1*d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {
        double w = d2 / (d2 - d6);
        return (s_point){{{A.x + w*AC.x, A.y + w*AC.y, A.z + w*AC.z}}};
    }

    // In edge BC
    double va = d3*d6 - d5*d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        s_point BC = {{{ C.x-B.x, C.y-B.y, C.z-B.z }}};
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return (s_point){{{B.x + w*BC.x, B.y + w*BC.y, B.z + w*BC.z}}};
    }

    // Inside face region
    // Barycentric coordinates (u,v,w)
    double denom = va + vb + vc;
    assert(fabs(denom) > 1e-9  && "Triangle is degenerate?");  // Triangle is degenarate?
    double v = vb / denom;
    double w = vc / denom;
    return (s_point){{{A.x + AB.x*v + AC.x*w, A.y + AB.y*v + AC.y*w, A.z + AB.z*v + AC.z*w}}};
}


s_point closest_point_on_segment(const s_point *segment, s_point p)
{
    s_point AB = subtract_points(segment[1], segment[0]);
    s_point pA = subtract_points(p, segment[0]);
    // Project c onto ab, computing parameterized position d(t) = a + t*(b – a)
    double denom = dot_prod(AB, AB);
    assert(denom > 1e-9);
    double t = dot_prod(pA, AB) / denom;
    // If outside segment, clamp t (and therefore d) to the closest endpoint
    if (t < 0) t = 0;
    if (t > 1) t = 1;
    return (s_point){{{segment[0].x + t*AB.x, segment[0].y + t*AB.y, segment[0].z + t*AB.z}}};
}


int point_in_triangle_2d(const s_point *triangle, s_point p)
{
    int o1 = orient2d(triangle[0].coords, triangle[1].coords, p.coords);
    int o2 = orient2d(triangle[1].coords, triangle[2].coords, p.coords);
    int o3 = orient2d(triangle[2].coords, triangle[0].coords, p.coords);
    
    // Find reference sign (non-zero)
    int signs[3] = {o1, o2, o3};
    int ref_sign = 0;
    for (int ii=0; ii<3; ii++) {
        if (signs[ii] != 0) {
            ref_sign = signs[ii];
            break;
        }
    }
    assert(ref_sign != 0);

    for (int ii=0; ii<3; ii++) {
        if (signs[ii] != 0 && signs[ii] != ref_sign) return 0;
    }
    return 1;
}


int point_in_triangle_3d(const s_point *triangle, s_point p)
{
    s_point aux = closest_point_on_triangle(triangle, p);
    if (distance_squared(aux, p) < 1e-6) return 1;
    else return 0;
}



ch_vertex *malloc_points_to_chvertex(const s_point *points, int Np)
{
    // ch_vertex *out = malloc(Np * sizeof(ch_vertex));
    // memcpy(out, points, Np * sizeof(ch_vertex));
    // return out;
    if (Np <= 0) return NULL;
    ch_vertex *out = malloc((size_t)Np * sizeof(ch_vertex));
    if (!out) { puts("ERROR: malloc_points_to_chvertex"); return NULL; }

    for (int i = 0; i < Np; ++i) {
        out[i].x = (CH_FLOAT) points[i].x;
        out[i].y = (CH_FLOAT) points[i].y;
        out[i].z = (CH_FLOAT) points[i].z;
        /* OR: out[i].v[0] = (CH_FLOAT)points[i].coords[0]; ... */
    }
    return out;
}


s_point *extract_normals_from_ch(const ch_vertex *vertices, int *faces, int Nf, s_point ch_CM, int NORMALIZE)
{
    s_point *out = malloc(sizeof(s_point) * Nf);
    assert(out != NULL);

    for (int ii = 0; ii < Nf; ii++) {
        ch_vertex v0 = vertices[faces[ii * 3 + 0]];
        ch_vertex v1 = vertices[faces[ii * 3 + 1]];
        ch_vertex v2 = vertices[faces[ii * 3 + 2]];

        s_point d1 = {{{ v1.x-v0.x, v1.y-v0.y, v1.z-v0.z }}};
        s_point d2 = {{{ v2.x-v0.x, v2.y-v0.y, v2.z-v0.z }}};
        s_point n = cross_prod(d1, d2);
        if (NORMALIZE != 0) n = normalize_3d(n);

        // Orientation check
        s_point verts_face[3] = { (s_point){{{v0.x, v0.y, v0.z}}}, 
                                  (s_point){{{v1.x, v1.y, v1.z}}}, 
                                  (s_point){{{v2.x, v2.y, v2.z}}} };
        int o = orientation(verts_face, ch_CM);
        // if (o == 0) {
        //     printf("%f, %f, %f\n", verts_face[0].x, verts_face[0].y, verts_face[0].z);
        //     printf("%f, %f, %f\n", verts_face[1].x, verts_face[1].y, verts_face[1].z);
        //     printf("%f, %f, %f\n", verts_face[2].x, verts_face[2].y, verts_face[2].z);
        // }
        assert(o != 0);
        if (o < 0) {
            n.x = -n.x;  n.y = -n.y;  n.z = -n.z;
            int tmp = faces[ii * 3 + 1];
            faces[ii * 3 + 1] = faces[ii * 3 + 2];
            faces[ii * 3 + 2] = tmp;
        }

        out[ii] = n;
    }
    return out;
}


int is_inside_convhull(s_point query, const s_point *pch, const int *faces, int Nf)
{   
    // 1: inside, 0: outise, -1: in boundary
    if (Nf == 0) { puts("Error: is_inside_convhull: Nf == 0"); exit(1); }

    int prev_sign = 0;
    for (int f = 0; f < Nf; ++f) {
        s_point pf[3] = {pch[faces[3*f + 0]],
                         pch[faces[3*f + 1]],
                         pch[faces[3*f + 2]]};
        
        int sign = orientation(pf, query);

        if (sign == 0) {  // Point is coplanar
            if (point_in_triangle_3d(pf, query)) {
                return -1;
            } else continue;
        }

        // if we've already seen a non-zero sign, it must match
        if (prev_sign == 0) prev_sign = sign;
        else if (sign != prev_sign) return 0;  // outside!
    }

    assert(prev_sign != 0  && "Point is coplanar with all faces? Strange...");
    return 1;
}


int is_in_boundary_convhull(const int *faces, int Nf, int vid)
{
    return inarray(faces, Nf * 3, vid);
}


int mark_inside_convhull(const s_point *query, int Np, const s_point *pch, const int *faces, int Nf, int *mark)
{
    memset(mark, 0, sizeof(int) * Np);
    int count = 0;
    for (int ii=0; ii<Np; ii++) {
        if (is_inside_convhull(query[ii], pch, faces, Nf)) mark[ii] = 1;
    }
    return count;
}


s_point random_point_uniform_3d(s_point min, s_point max)
{
    double ux = rand() / ((double) RAND_MAX + 1.0);
    double uy = rand() / ((double) RAND_MAX + 1.0);
    double uz = rand() / ((double) RAND_MAX + 1.0);
    s_point out;
    out.x = min.x + (max.x - min.x) * ux;
    out.y = min.y + (max.y - min.y) * uy;
    out.z = min.z + (max.z - min.z) * uz;
    return out;
}


s_point random_point_inside_convhull(const s_point *pch, const int *faces, int Nf, 
                                     s_point min, s_point max)
{
    int MAX_IT = 10000;
    int it = 0;
    s_point out = random_point_uniform_3d(min, max);
    while (is_inside_convhull(out, pch, faces, Nf) != 1) {
        out = random_point_uniform_3d(min, max);
        assert(it < MAX_IT && "Reached maximum iters looking for point inside convhull.");
        it++;
    }
    return out;
}


double volume_tetrahedron_approx(s_point p1, s_point p2, s_point p3, s_point p4)
{
    return fabs(1.0/6.0 * orient3d(p1.coords, p2.coords, p3.coords, p4.coords));
}


double compute_volume_convhull(const s_point *points, const int *faces, const s_point *fnormals_UNNORMALIZED, int Nf)
{
    double vol = 0;
    for (int ii=0; ii<Nf; ii++) {
        double Nx = fnormals_UNNORMALIZED[ii].x;
        vol += Nx * (points[faces[ii*3 + 0]].x +
                     points[faces[ii*3 + 1]].x +
                     points[faces[ii*3 + 2]].x);
    }
    return vol / 6;
}


double compute_volume_convhull_from_points(const s_point *points, int Np)
{
    ch_vertex *ch_vertices = malloc_points_to_chvertex(points, Np);
    int *faces;
    int N_faces;
    convhull_3d_build(ch_vertices, Np, &faces, &N_faces);

    s_point CM = find_center_mass(points, Np);
    s_point *fnormals_UNNORMALIZED = extract_normals_from_ch(ch_vertices, faces, N_faces, CM, 0);
    
    double volume = compute_volume_convhull(points, faces, fnormals_UNNORMALIZED, N_faces);

    free(ch_vertices);
    free(faces);
    free(fnormals_UNNORMALIZED);
    return volume;
}


void convhull_from_points(const s_point *points, int Np, int **faces, s_point **fnormals, int *Nf)
{
    ch_vertex *ch_vertices = malloc_points_to_chvertex(points, Np);
    convhull_3d_build(ch_vertices, Np, faces, Nf);
    if (fnormals) {
        s_point CM = find_center_mass(points, Np);
        *fnormals = extract_normals_from_ch(ch_vertices, *faces, *Nf, CM, 1);
    }
    free(ch_vertices);
}


