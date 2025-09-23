
#include "geometry.h"
#include "algebra.h"
#include "predicates.h"
#include "convhull_3d.h"
#include "array_operations.c"
#include <stdio.h>
#include <string.h>
#include <math.h>


int orientation(double **p, double *q, int dim)
{
    if (dim == 2) {
        double aux = orient2d(p[0], p[1], q);
        if (aux > 0) return 1;
        else if (aux < 0) return -1;
        else return 0;
    }
    else if (dim == 3) {
        double aux = orient3d(p[0], p[1], p[2], q);
        if (aux > 0) return 1;
        else if (aux < 0) return -1;
        else return 0;
    }
    else {
        printf("orientation not implemented for this dimension.");
        exit(1);
    }
}


int in_sphere(double **p, double *q, int dim)
{   
    if (dim == 2) {
        int factor;
        if (orientation(p, p[dim], 2) == 1) factor = 1;
        else factor = -1;

        double aux = incircle(p[0], p[1], p[2], q);
        
        if (aux > 0) return factor;
        else if (aux < 0) return -factor;
        else return 0;
        
    } else if (dim == 3) {
        int factor;
        if (orientation(p, p[dim], 3) == 1) factor = 1;
        else factor = -1;

        double aux = insphere(p[0], p[1], p[2], p[3], q);
        
        if (aux > 0) return factor;
        else if (aux < 0) return -factor;
        else return 0;

    } else {
        printf("orientation not implemented for this dimension.");
        exit(1);
    }
}


int segment_crosses_triangle_3d(double **triangle, double *a, double *b)
{
    double STORAGE[3*3];
    double *aux[3] = {STORAGE, STORAGE + 3, STORAGE + 6};

    aux[2][0] = a[0];     aux[2][1] = a[1];     aux[2][2] = a[2];
    
    if (orientation(triangle, a, 3) == orientation(triangle, b, 3)) return 0;

    aux[0][0] = triangle[0][0];     aux[0][1] = triangle[0][1];     aux[0][2] = triangle[0][2];
    aux[1][0] = triangle[1][0];     aux[1][1] = triangle[1][1];     aux[1][2] = triangle[1][2];
    int s1 = orientation(aux, b, 3);

    aux[0][0] = triangle[1][0];     aux[0][1] = triangle[1][1];     aux[0][2] = triangle[1][2];
    aux[1][0] = triangle[2][0];     aux[1][1] = triangle[2][1];     aux[1][2] = triangle[2][2];
    int s2 = orientation(aux, b, 3);

    aux[0][0] = triangle[2][0];     aux[0][1] = triangle[2][1];     aux[0][2] = triangle[2][2];
    aux[1][0] = triangle[0][0];     aux[1][1] = triangle[0][1];     aux[1][2] = triangle[0][2];
    int s3 = orientation(aux, b, 3);
    
    if (s1 == s2 && s2 == s3 && s3 == s1) return 1;
    else return 0;
}


int between_1d(double x, double a, double b, double eps)
{
    // returns 1 if x lies between [a,b]
    if (a > b) { double t=a; a=b; b=t; }
    return (x + eps >= a && x - eps <= b);
}


int segments_intersect_2d(double *A, double *B, double *p, double *d)
{
    const double EPS = 1e-9;
    double Ax = A[0], Ay = A[1];
    double Bx = B[0], By = B[1];
    double px = p[0], py = p[1];
    double dx = d[0], dy = d[1];

    // 1) Degenerate pd, treat as X
    if (fabs(px - dx) < EPS && fabs(py - dy) < EPS) {
        double X[2] = { px, py };
        if (fabs(Ax - Bx) < EPS && fabs(Ay - By) < EPS) {
            // If A, B also degenerate
            return (fabs(px - Ax) < EPS && fabs(py - Ay) < EPS);
        }
        double *auxAB[2] = { A, B };
        if (orientation(auxAB, X, 2) == 0 &&
            between_1d(px, Ax, Bx, EPS) &&
            between_1d(py, Ay, By, EPS)) {
            return 1;
        }
        return 0;
    }

    // 2) Degenerate AB, treat as Y
    if (fabs(Ax - Bx) < EPS && fabs(Ay - By) < EPS) {
        double Y[2] = { Ax, Ay };
        double *auxPD[2] = { p, d };
        if (orient2d(auxPD[0], auxPD[1], Y) == 0 &&
            between_1d(Ax, px, dx, EPS) &&
            between_1d(Ay, py, dy, EPS)) {
            return 1;
        }
        return 0;
    }

    // 3) General case: two “straddling” tests
    // 3a) [A,B] vs {p,d}
    double *auxAB[2] = { A, B };
    int o1 = orient2d(auxAB[0], auxAB[1], p);
    int o2 = orient2d(auxAB[0], auxAB[1], d);
    if (o1 != 0 && o2 != 0 && o1 == o2) return 0;

    // 3b) [p,d] vs {A,B}
    double *auxPD[2] = { p, d };
    o1 = orient2d(auxPD[0], auxPD[1], A);
    o2 = orient2d(auxPD[0], auxPD[1], B);
    if (o1 != 0 && o2 != 0 && o1 == o2) return 0;

    return 1;
}


void find_center_mass(double **in, int N_points, int dim, double *out)
{
    for (int ii=0; ii<dim; ii++) {
        out[ii] = in[0][ii];
    }

    for (int ii=1; ii<N_points; ii++) {
        for (int jj=0; jj<dim; jj++) {
            out[jj] += in[ii][jj];
        }
    }

    for (int ii=0; ii<dim; ii++) {
        out[ii] /= N_points;
    }
}


int coord_with_largest_component_3d(double *n)
{
    int out = 2;
    if (fabs(n[0]) > fabs(n[1]) && fabs(n[0]) > fabs(n[2])) out = 0;
    else if (fabs(n[1]) > fabs(n[0]) && fabs(n[1]) > fabs(n[2])) out = 1;
    return out;
}


void cross_3d(const double *u, const double *v, double *out)
{
    out[0] = u[1]*v[2] - u[2]*v[1];
    out[1] = u[2]*v[0] - u[0]*v[2];
    out[2] = u[0]*v[1] - u[1]*v[0];
}


double dot_3d(const double *u, const double *v)
{
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}


void subtract_3d(const double *u, const double *v, double *result)
{
    result[0] = u[0] - v[0];
    result[1] = u[1] - v[1];
    result[2] = u[2] - v[2];
}


double distance_squared(const double *a, const double *b)
{
    double diff[3];
    subtract_3d(a, b, diff);
    return dot_3d(diff, diff);
}


void closest_point_on_triangle(const double *A, const double *B, const double *C, const double *p,
                               double *c_out) 
{
    double AB[3] = {B[0]-A[0], B[1]-A[1], B[2]-A[2]},
           AC[3] = {C[0]-A[0], C[1]-A[1], C[2]-A[2]},
           AP[3] = {p [0]-A[0], p [1]-A[1], p [2]-A[2]};

    double d1 = dot_3d(AB, AP), d2 = dot_3d(AC, AP);
    if (d1 <= 0 && d2 <= 0) { memcpy(c_out, A, 3*sizeof(double)); return; }

    double BP[3] = { p[0]-B[0], p[1]-B[1], p[2]-B[2] };
    double d3 = dot_3d(AB, BP), d4 = dot_3d(AC, BP);
    if (d3 >= 0 && d4 <= d3) { memcpy(c_out, B, 3*sizeof(double)); return; }

    double vc = d1*d4 - d3*d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) {
        double v = d1 / (d1 - d3);
        c_out[0] = A[0] + v * AB[0];
        c_out[1] = A[1] + v * AB[1];
        c_out[2] = A[2] + v * AB[2];
        return;
    }

    double CP[3] = { p[0]-C[0], p[1]-C[1], p[2]-C[2] };
    double d5 = dot_3d(AB, CP),
           d6 = dot_3d(AC, CP);
    if (d6 >= 0 && d5 <= d6) { memcpy(c_out, C, 3*sizeof(double)); return; }

    // Check edge AC region
    double vb = d5*d2 - d1*d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {
        double w = d2 / (d2 - d6);
        c_out[0] = A[0] + w * AC[0];
        c_out[1] = A[1] + w * AC[1];
        c_out[2] = A[2] + w * AC[2];
        return;
    }

    // Check edge BC region
    double va = d3*d6 - d5*d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        // projection onto BC
        double BC[3] = {C[0]-B[0], C[1]-B[1], C[2]-B[2]};
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        c_out[0] = B[0] + w * BC[0];
        c_out[1] = B[1] + w * BC[1];
        c_out[2] = B[2] + w * BC[2];
        return;
    }

    // Inside face region
    // Barycentric coordinates (u,v,w)
    double denom = va + vb + vc;
    assert(fabs(denom) > 1e-9);
    double v = vb / denom;
    double w = vc / denom;
    c_out[0] = A[0] + AB[0]*v + AC[0]*w;
    c_out[1] = A[1] + AB[1]*v + AC[1]*w;
    c_out[2] = A[2] + AB[2]*v + AC[2]*w;
}


void closest_point_on_segment(double *p, double *A, double *B, double *OUT)
{
    double ab[3] = { B[0]-A[0], B[1]-A[1], B[2]-B[2]};
    // Project c onto ab, computing parameterized position d(t) = a + t*(b – a)
    double pa[3]; subtract_3d(p, A, pa);
    double t = dot_3d(pa, ab) / dot_3d(ab, ab);
    // If outside segment, clamp t (and therefore d) to the closest endpoint
    if (t < 0) t = 0;
    if (t > 1) t = 1;
    // Compute projected position from the clamped t
    OUT[0] = A[0] + t * ab[0];
    OUT[1] = A[1] + t * ab[1];
    OUT[2] = A[2] + t * ab[2];
}


int point_in_triangle_2d(double *v1, double *v2, double *v3, double *p)
{
    int o1 = orient2d(v1, v2, p);
    int o2 = orient2d(v2, v3, p);
    int o3 = orient2d(v3, v1, p);
    
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


ch_vertex *convert_points_to_chvertex(double **points, int Np)
{
    ch_vertex *out = malloc(sizeof(ch_vertex) * Np);
    for (int ii=0; ii<Np; ii++) {
        out[ii].x = points[ii][0];
        out[ii].y = points[ii][1];
        out[ii].z = points[ii][2];
    }
    return out;
}


double **extract_normals_from_ch(ch_vertex *vertices, int *faces, int Nf, double *ch_CM)
{   // Normalized normals
    double STORAGE[3*3];
    double *vertices_face[3] = {STORAGE, STORAGE + 3, STORAGE + 6};

    double **out = malloc_matrix(Nf, 3);
    for (int ii=0; ii<Nf; ii++) {
        ch_vertex v0 = vertices[faces[ii*3+0]];
        ch_vertex v1 = vertices[faces[ii*3+1]];
        ch_vertex v2 = vertices[faces[ii*3+2]];

        double d1[3], d2[3];
        d1[0] = v1.x - v0.x;
        d1[1] = v1.y - v0.y;
        d1[2] = v1.z - v0.z;
        d2[0] = v2.x - v0.x;
        d2[1] = v2.y - v0.y;
        d2[2] = v2.z - v0.z;

        double n[3];
        cross_3d(d1, d2, n);
        normalize_inplace(n, 3);
        
        vertices_face[0][0] = v0.x;     vertices_face[0][1] = v0.y;     vertices_face[0][2] = v0.z;
        vertices_face[1][0] = v1.x;     vertices_face[1][1] = v1.y;     vertices_face[1][2] = v1.z;
        vertices_face[2][0] = v2.x;     vertices_face[2][1] = v2.y;     vertices_face[2][2] = v2.z;

        int o = orientation(vertices_face, ch_CM, 3);
        assert(o != 0);
        if (o < 0) {
            n[0] = -n[0];
            n[1] = -n[1];
            n[2] = -n[2];
            int tmp = faces[ii*3+1];
            faces[ii*3+1] = faces[ii*3+2];
            faces[ii*3+2] = tmp;
        }

        out[ii][0] = n[0];
        out[ii][1] = n[1];
        out[ii][2] = n[2];
    }
    return out;
}


double **extract_normals_from_ch_UNNORMALIZED(ch_vertex *vertices, int *faces, int Nf, double *ch_CM)
{ 
    double STORAGE[3*3];
    double *vertices_face[3] = {STORAGE, STORAGE + 3, STORAGE + 6};


    double **out = malloc_matrix(Nf, 3);
    for (int ii=0; ii<Nf; ii++) {
        ch_vertex v0 = vertices[faces[ii*3+0]];
        ch_vertex v1 = vertices[faces[ii*3+1]];
        ch_vertex v2 = vertices[faces[ii*3+2]];

        double d1[3], d2[3];
        d1[0] = v1.x - v0.x; d1[1] = v1.y - v0.y; d1[2] = v1.z - v0.z;
        d2[0] = v2.x - v0.x; d2[1] = v2.y - v0.y; d2[2] = v2.z - v0.z;

        double n[3];
        cross_3d(d1, d2, n);
        
        vertices_face[0][0] = v0.x; vertices_face[0][1] = v0.y; vertices_face[0][2] = v0.z;
        vertices_face[1][0] = v1.x; vertices_face[1][1] = v1.y; vertices_face[1][2] = v1.z;
        vertices_face[2][0] = v2.x; vertices_face[2][1] = v2.y; vertices_face[2][2] = v2.z;

        double fc[3], dir[3];
        find_center_mass(vertices_face, 3, 3, fc);
        subtract_3d(fc, ch_CM, dir);
        double dot = dot_3d(dir, n);
        if (dot < 0) {
            n[0] = -n[0];
            n[1] = -n[1];
            n[2] = -n[2];
        }

        out[ii][0] = n[0];
        out[ii][1] = n[1];
        out[ii][2] = n[2];
    }
    return out;
}


int is_inside_convhull(double *query, double **pch, int *faces, __attribute__((unused)) double **fnormals, int Nf)
{   
    int prev_sign = 0;
    for (int f = 0; f < Nf; ++f) {
        double *pf[3] = { pch[faces[3*f + 0]],
                          pch[faces[3*f + 1]],
                          pch[faces[3*f + 2]] };
        
        double dot = orient3d(pf[0], pf[1], pf[2], query);

        int sign = (dot > 0 ? +1 : -1);
        // if we've already seen a non-zero sign, it must match
        if (prev_sign != 0 && sign != prev_sign) return 0;  // we are on the "other" side

        prev_sign = sign;
    }

    if (Nf == 0) {
        puts("What? Nf == 0?");
        exit(1);
    }
    assert(prev_sign != 0);

    return 1;
}


int is_in_boundary_convhull(int *faces, int Nf, int vid)
{
    return inarray(faces, Nf * 3, vid);
}


void random_point_uniform_3d(double *min, double *max, double *out)
{
    double ux = rand() / ((double) RAND_MAX + 1.0);
    double uy = rand() / ((double) RAND_MAX + 1.0);
    double uz = rand() / ((double) RAND_MAX + 1.0);

    out[0] = min[0] + (max[0] - min[0]) * ux;
    out[1] = min[1] + (max[1] - min[1]) * uy;
    out[2] = min[2] + (max[2] - min[2]) * uz;
}


void random_point_inside_convhull(double **pch, int *faces, double **fnormals, int Nf, 
                                  double *min, double *max, double *out)
{
    int MAX_IT = 1000;
    int it = 0;
    random_point_uniform_3d(min, max, out);
    while (!is_inside_convhull(out, pch, faces, fnormals, Nf)) {
        random_point_uniform_3d(min, max, out);
        assert(it < MAX_IT && "Reached maximum iters looking for point inside convhull.");
        it++;
    }
}



int mark_inside_convhull(double **points, int Np, double **pch, int *faces, double **fnormals, int Nf, int *mark)
{
    memset(mark, 0, sizeof(int) * Np);
    int count = 0;
    for (int ii=0; ii<Np; ii++) {
        if (is_inside_convhull(points[ii], pch, faces, fnormals, Nf)) mark[ii] = 1;
    }
    return count;
}


double volume_tetrahedron_approx(double *p1, double *p2, double *p3, double *p4)
{
    return fabs(1.0/6.0 * orient3d(p1, p2, p3, p4));
}


double compute_volume_convhull(double **points, int *faces, double **fnormals, int Nf)
{
    double vol = 0;
    for (int ii=0; ii<Nf; ii++) {
        double Nx = fnormals[ii][0];
        vol += Nx * (points[faces[ii*3 + 0]][0] +
                     points[faces[ii*3 + 1]][0] +
                     points[faces[ii*3 + 2]][0]);
    }
    return vol / 6;
}


double compute_volume_convhull_from_points(double **points, int Np)
{
    ch_vertex *ch_vertices = convert_points_to_chvertex(points, Np);
    int *faces;
    int N_faces;
    convhull_3d_build(ch_vertices, Np, &faces, &N_faces);

    double CM[3];
    find_center_mass(points, Np, 3, CM);
    double **fnormals = extract_normals_from_ch_UNNORMALIZED(ch_vertices, faces, N_faces, CM);
    
    double volume = compute_volume_convhull(points, faces, fnormals, N_faces);

    free(ch_vertices);
    free(faces);
    free_matrix(fnormals, N_faces);
    return volume;
}


void extract_faces_convhull_from_points(double **points, int Np, int **faces, double ***fnormals, int *Nf)
{
    ch_vertex *ch_vertices = convert_points_to_chvertex(points, Np);
    convhull_3d_build(ch_vertices, Np, faces, Nf);

    if (fnormals) {
        double CM[3];
        find_center_mass(points, Np, 3, CM);
        *fnormals = extract_normals_from_ch(ch_vertices, *faces, *Nf, CM);
    }

    free(ch_vertices);
}


void segment_convex_hull_intersection(const double *p0, const double *p1, double **pch, const int *faces, 
                                     double **fnormals, int Nf, int *ind1, double *i1, int *ind2, double *i2)
{       // ind1 indicates if i1 is a valid intersection, simil 2
    double tmin = 0, tmax = 1;

    double p1_p0[3];
    subtract_3d(p1, p0, p1_p0);

    for (int ii=0; ii<Nf; ii++) {
        double *n = fnormals[ii];
        double den = dot_3d(n, p1_p0);
        
        assert(fabs(den) > 1e-12 && "Perpendicular normal with segment.");

        double *q = pch[faces[ii*3]];
        double q_p0[3];     
        subtract_3d(q, p0, q_p0);
        double t = dot_3d(n, q_p0) / den;

        if (den < 0) {
            tmin = fmax(tmin, t);
        } else if (den >0) {
            tmax = fmin(tmax, t);
        }
    }
    
    *ind1 = 0;
    *ind2 = 0;

    if (tmin < 0) {
        printf("WAGNING! tmin < 0\n");
        tmin = 0;
    }
    if (tmax > 0) {
        printf("WAGNING! tmax > 0\n");
        tmax = 1;
    }
    if (tmin < tmax) {
        *ind1 = 1;
        i1[0] = p0[0] + tmin * p1_p0[0];  
        i1[1] = p0[1] + tmin * p1_p0[1];  
        i1[2] = p0[2] + tmin * p1_p0[2];

        i2[0] = p0[0] + tmax * p1_p0[0];  
        i2[1] = p0[1] + tmax * p1_p0[1];  
        i2[2] = p0[2] + tmax * p1_p0[2];
    }
}


