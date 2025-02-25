#ifndef GEOMETRY_C
#define GEOMETRY_C

#include "algebra.c"


int orientation(double **p, double *q, int dim)
{
    double **M = malloc_matrix(dim+1, dim+1);
    for (int ii=0; ii<dim; ii++) {
        for (int jj=0; jj<dim; jj++) {
            M[ii][jj] = p[ii][jj] - q[jj];
        }
    }
    double det = determinant(M, dim);
    if (det > 0) {
        return 1;
    } else if (det < 0) {
        return -1;
    } else {
        return 0;
    }
}


int insphere(double **p, double *q, int dim)
{   
    static double **p_aux = NULL;
    static int prev_dim = 0;
    if (dim != prev_dim) {
        if (p_aux) free_matrix(p_aux, prev_dim+1);
        p_aux = malloc_matrix(dim+1, dim+1);
        prev_dim = dim;
    }
    copy_matrix(p, p_aux, dim+1, dim);

    // check if p is oriented, if not, orient p_aux!
    int oriented = orientation(p, p[dim], dim);  
    if (oriented == -1) {
        double *row_aux = p_aux[0];
        p_aux[0] = p_aux[1];
        p_aux[1] = row_aux;
    }
    assert(orientation(p_aux, p_aux[dim], dim) == 1 && "points are not oriented");

    double norm_q_2 = norm_squared(q, dim);
    double **M = malloc_matrix(dim+2, dim+2);
    for (int ii=0; ii<dim+1; ii++) {
        for (int jj=0; jj<dim; jj++) {
            M[ii][jj] = p_aux[ii][jj] - q[jj];
        }
        M[ii][dim] = norm_squared(p_aux[ii], dim) - norm_q_2;
    }
    double det = determinant(M, dim+1);
    if (det > 0) {
        return 1;
    } else if (det < 0) {
        return -1;
    } else {
        return 0;
    }
}


int segment_crosses_triangle_3d(double **triangle, double *a, double *b)
{
    static double **aux = NULL;  // multithreaded will not work!
    if (!aux) aux = malloc_matrix(3, 3);
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


#endif
