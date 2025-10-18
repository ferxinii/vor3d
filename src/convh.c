
#include "convh.h"
#include "algebra.h"
#include "geometry.h"
#include "convhull_3d.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>


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


