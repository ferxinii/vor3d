
#include "algebra.h"
#include "geometry.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>


int id_where_equal_int(const int *arr, int N, int entry) 
{
    for (int ii=0; ii<N; ii++) {
        if (arr[ii] == entry) return ii;
    }
    assert(1 == 0 && "Could not find id."); 
    exit(1);
}


int inarray(const int *arr1, int N, int a)
{
    for (int ii=0; ii<N; ii++) {
        if (arr1[ii] == a) return 1;
    }
    return 0;
}


double **malloc_matrix(int N1, int N2)
{
    double **array = malloc(sizeof(double*) * N1);
    for (int ii=0; ii<N1; ii++) {
        array[ii] = malloc(sizeof(double) * N2);
    }
    return array;
}


double **realloc_matrix(double **matrix, int old_rows, int new_rows, int ncols)
{
    if (new_rows < old_rows) {
        for (int ii = new_rows; ii < old_rows; ii++) {
            free(matrix[ii]);
        }
    }

    double **new_matrix = realloc(matrix, sizeof(double*) * new_rows);

    for (int ii=old_rows; ii<new_rows; ii++) {
        new_matrix[ii] = malloc(sizeof(double) * ncols);
    }
    return new_matrix;
}


int **malloc_matrix_int(int N1, int N2)
{
    int **array = malloc(sizeof(int*) * N1);
    for (int ii=0; ii<N1; ii++) {
        array[ii] = malloc(sizeof(int) * N2);
    }
    return array;
}


int **realloc_matrix_int(int **matrix, int old_rows, int new_rows, int ncols)
{
    if (new_rows < old_rows) {
        for (int ii = new_rows; ii < old_rows; ii++) {
            free(matrix[ii]);
        }
    }

    int **new_matrix = realloc(matrix, sizeof(int*) * new_rows);

    for (int ii=old_rows; ii<new_rows; ii++) {
        new_matrix[ii] = malloc(ncols * sizeof(int));
    }
    return new_matrix;
}


void free_matrix(double **array, int N1)
{
    for (int ii=0; ii<N1; ii++) {
        free(array[ii]);
    }
    free(array);
}


void free_matrix_int(int **array, int N1)
{
    for (int ii=0; ii<N1; ii++) {
        free(array[ii]);
    }
    free(array);
}


void copy_matrix(double **in, double **out, int N1, int N2)
{
    for (int ii=0; ii<N1; ii++) {
        for (int jj=0; jj<N2; jj++) {
            out[ii][jj] = in[ii][jj];
        }
    }
}


void copy_matrix_int(int **in, int **out, int N1, int N2)
{
    for (int ii=0; ii<N1; ii++) {
        for (int jj=0; jj<N2; jj++) {
            out[ii][jj] = in[ii][jj];
        }
    }
}


void print_matrix(double **array, int N1, int N2)
{
    for (int ii=0; ii<N1; ii++) {
        for (int jj=0; jj<N2; jj++) {
            printf("%f ", array[ii][jj]);
        }
        printf("\n");
    }
}


double norm_squared(const double *v, int dim)
{
    double out = 0;
    for (int ii=0; ii<dim; ii++) {
        out += v[ii] * v[ii];
    }
    return out;
}


double norm_difference(const double *a, const double *b, int dim)
{
    double out = 0;
    for (int ii=0; ii<dim; ii++) {
        out += (a[ii] - b[ii]) * (a[ii] - b[ii]);
    }
    return sqrt(out);
}


double norm_difference_squared(const double *a, const double *b, int dim)
{
    double out = 0;
    for (int ii=0; ii<dim; ii++) {
        out += (a[ii] - b[ii]) * (a[ii] - b[ii]);
    }
    return out;
}


double max_distance(double **p, int N, int dim, double *q)
{   
    double maxd = 0;
    for (int ii=0; ii<N; ii++) {
        double d = norm_difference_squared(p[ii], q, dim);
        if (maxd < d) maxd = d;
    }
    return sqrt(maxd);
}


s_point normalize_3d(s_point v)
{
    double norm = sqrt(norm_squared(v.coords, 3));
    return (s_point) {v.x/norm, v.y/norm, v.z/norm};
}

