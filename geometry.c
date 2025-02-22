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
    // check if p is oriented
    int oriented = orientation(p, p[dim], dim);  
    assert(oriented == 1 && "points are not oriented");

    double norm_q_2 = norm_squared(q, dim);
    double **M = malloc_matrix(dim+2, dim+2);
    for (int ii=0; ii<dim+1; ii++) {
        for (int jj=0; jj<dim; jj++) {
            M[ii][jj] = p[ii][jj] - q[jj];
        }
        M[ii][dim] = norm_squared(p[ii], dim) - norm_q_2;
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

#endif
