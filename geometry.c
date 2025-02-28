#ifndef GEOMETRY_C
#define GEOMETRY_C

#include "algebra.c"
#include "predicates.h"


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
    // static double **M = NULL;
    // static int dim_prev = 0;
    // if (dim != dim_prev) {
    //     if (M) free(M);
    //     M = malloc_matrix(dim, dim);
    //     dim_prev = dim;
    // }
    //
    // for (int ii=0; ii<dim; ii++) {
    //     for (int jj=0; jj<dim; jj++) {
    //         M[ii][jj] = p[ii][jj] - q[jj];
    //     }
    // }
    // 
    // const double TOL = 1e-6;
    //
    // double det = determinant(M, dim);
    // if (fabs(det) < TOL) {
    //     puts("Determinant inside tolerance.");
    //     return 0;
    // }
    // if (det > 0) {
    //     return 1;
    // } else if (det < 0) {
    //     return -1;
    // } else {
    //     return 0;
    // }
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

        // else {
        //     puts("WHAT TO DO? Simply return output of either?");
        //     exit(1);
        // }
    } else {
        printf("orientation not implemented for this dimension.");
        exit(1);
    }
    // static double **p_aux = NULL;
    // static int prev_dim = 0;
    // if (dim != prev_dim) {
    //     if (p_aux) free_matrix(p_aux, prev_dim+1);
    //     p_aux = malloc_matrix(dim+1, dim+1);
    //     prev_dim = dim;
    // }
    // copy_matrix(p, p_aux, dim+1, dim);
    //
    // // check if p is oriented, if not, orient p_aux!
    // int oriented = orientation(p, p[dim], dim);  
    // if (oriented == 0) {
    //     puts("COPLANAR POINTS! Is this legal?\n");
    //     return 0;
    // }
    // if (oriented == -1) {
    //     double *row_aux = p_aux[0];
    //     p_aux[0] = p_aux[1];
    //     p_aux[1] = row_aux;
    //     if (oriented == orientation(p_aux, p_aux[dim], dim) ) {
    //         printf("WHAT?, %d, %d\n", oriented, orientation(p_aux, p_aux[dim], dim));
    //         print_matrix(p, dim+1, dim);
    //         print_matrix(p_aux, dim+1, dim);
    //     }
    // }
    // assert(orientation(p_aux, p_aux[dim], dim) == 1 && "points are not oriented");
    //
    // double norm_q_2 = norm_squared(q, dim);
    // double **M = malloc_matrix(dim+2, dim+2);
    // for (int ii=0; ii<dim+1; ii++) {
    //     for (int jj=0; jj<dim; jj++) {
    //         M[ii][jj] = p_aux[ii][jj] - q[jj];
    //     }
    //     M[ii][dim] = norm_squared(p_aux[ii], dim) - norm_q_2;
    // }
    //
    // const double TOL = 1e-6;
    // double det = determinant(M, dim+1);
    // if (fabs(det) < TOL) {
    //     puts("Determinant inside tolerance.");
    //     return 0;
    // }
    // if (det > 0) {
    //     return 1;
    // } else if (det < 0) {
    //     return -1;
    // } else {
    //     return 0;
    // }
}


// --------------------------- SEGMENT TRIANGLE INTERSECTION --------------------------------

// Compute the cross product: result = u x v
void cross(const double u[3], const double v[3], double result[3]) {
    result[0] = u[1]*v[2] - u[2]*v[1];
    result[1] = u[2]*v[0] - u[0]*v[2];
    result[2] = u[0]*v[1] - u[1]*v[0];
}

// Compute the dot product of two vectors
double dot(const double u[3], const double v[3]) {
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2];
}

// Compute the vector difference: result = u - v
void subtract(const double u[3], const double v[3], double result[3]) {
    result[0] = u[0] - v[0];
    result[1] = u[1] - v[1];
    result[2] = u[2] - v[2];
}

/*
 * Returns 1 if the segment from p0 to p1 intersects the triangle defined by vertices v0, v1, v2;
 * otherwise, returns 0.
 *
 * This implementation is based on the Möller–Trumbore algorithm for ray–triangle intersection.
 * The segment is parameterized such that t = 0 corresponds to p0 and t = 1 to p1.
 */
int segment_crosses_triangle_3d(double **triangle, double *p0, const double *p1) 
{
    const double EPSILON = 1e-9;

    double *v0 = triangle[0];
    double *v1 = triangle[1];
    double *v2 = triangle[2];

    double d[3];
    subtract(p1, p0, d); // Direction of the segment

    double e1[3], e2[3];
    subtract(v1, v0, e1);
    subtract(v2, v0, e2);

    double h[3];
    cross(d, e2, h);
    double a = dot(e1, h);
    
    if (fabs(a) < EPSILON)
        return 0;  // The segment is parallel to the triangle plane or triangle is degenerate

    double f = 1.0 / a;
    double s[3];
    subtract(p0, v0, s);
    // u and v are the barycentric coordinates
    double u = f * dot(s, h);
    if (u < 0.0 || u > 1.0)
        return 0;  // The intersection lies outside of the triangle

    double q[3];
    cross(s, e1, q);
    double v = f * dot(d, q);
    if (v < 0.0 || u + v > 1.0)
        return 0;  // The intersection lies outside of the triangle

    double t = f * dot(e2, q);
    if (t < 0.0 || t > 1.0)
        return 0;  // The intersection point is not on the segment

    return 1;
}


// int segment_crosses_triangle_3d(double **triangle, double *a, double *b)
// {
//     static double **aux = NULL;  // multithreaded will not work!
//     if (!aux) aux = malloc_matrix(3, 3);
//     aux[2][0] = a[0];     aux[2][1] = a[1];     aux[2][2] = a[2];
//     
//     if (orientation(triangle, a, 3) == orientation(triangle, b, 3)) return 0;
//
//     aux[0][0] = triangle[0][0];     aux[0][1] = triangle[0][1];     aux[0][2] = triangle[0][2];
//     aux[1][0] = triangle[1][0];     aux[1][1] = triangle[1][1];     aux[1][2] = triangle[1][2];
//     int s1 = orientation(aux, b, 3);
//
//     aux[0][0] = triangle[1][0];     aux[0][1] = triangle[1][1];     aux[0][2] = triangle[1][2];
//     aux[1][0] = triangle[2][0];     aux[1][1] = triangle[2][1];     aux[1][2] = triangle[2][2];
//     int s2 = orientation(aux, b, 3);
//
//     aux[0][0] = triangle[2][0];     aux[0][1] = triangle[2][1];     aux[0][2] = triangle[2][2];
//     aux[1][0] = triangle[0][0];     aux[1][1] = triangle[0][1];     aux[1][2] = triangle[0][2];
//     int s3 = orientation(aux, b, 3);
//     
//     if (s1 == s2 && s2 == s3 && s3 == s1) return 1;
//     else return 0;
// }


#endif
