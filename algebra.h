#ifndef ALGEBRA_H
#define ALGEBRA_H

double **malloc_matrix(int N1, int N2);
double **realloc_matrix(double **matrix, int old_rows, int new_rows, int ncols);
int **malloc_matrix_int(int N1, int N2);
int **realloc_matrix_int(int **matrix, int old_rows, int new_rows, int ncols);
void free_matrix(double **array, int N1);
void free_matrix_int(int **array, int N1);
void copy_matrix(double **in, double **out, int N1, int N2);
void copy_matrix_int(int **in, int **out, int N1, int N2);
void print_matrix(double **array, int N1, int N2);
double norm_squared(const double *v, int dim);
double norm_difference(const double *a, const double *b, int dim);
double max_distance(double **p, int N, int dim, double *q);
void normalize_inplace(double *v, int dim);

#endif
