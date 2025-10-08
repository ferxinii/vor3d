#ifndef VOR3D_ALGEBRA_H
#define VOR3D_ALGEBRA_H

int id_where_equal_int(const int *arr, int N, int entry);
int inarray(const int *arr1, int N, int a);

double **malloc_matrix(int N1, int N2);
double **realloc_matrix(double **matrix, int old_rows, int new_rows, int ncols);
int **malloc_matrix_int(int N1, int N2);
int **realloc_matrix_int(int **matrix, int old_rows, int new_rows, int ncols);
void free_matrix(double **array, int N1);
void free_matrix_int(int **array, int N1);
void copy_matrix(double **in, double **out, int N1, int N2);
void copy_matrix_int(int **in, int **out, int N1, int N2);
void print_matrix(double **array, int N1, int N2);

#endif
