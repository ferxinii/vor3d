#ifndef ALGEBRA_C
#define ALGEBRA_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>


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
    return sqrt(out);
}


double norm_difference(const double *a, const double *b, int dim)
{
    double out = 0;
    for (int ii=0; ii<dim; ii++) {
        out += (a[ii] - b[ii]) * (a[ii] - b[ii]);
    }
    return sqrt(out);
}


double max_distance(double **p, int N, int dim, double *q)
{   
    double maxd = 0;
    for (int ii=0; ii<N; ii++) {
        double d = norm_difference(p[ii], q, dim);
        if (maxd < d) maxd = d;
    }
    return maxd;
}


void normalize_inplace(double *v, int dim)
{
    double norm = sqrt(norm_squared(v, dim));
    for (int ii=0; ii<dim; ii++) {
        v[ii] = v[ii] / norm;
    }
}


int solve3x3(const double A[3][3], const double b[3], double x[3])
{
    double TOL = 1e-9;
    // Create an augmented matrix M of size 3x4
    double M[3][4];
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            M[i][j] = A[i][j];
        }
        M[i][3] = b[i];
    }
    
    // Forward elimination with partial pivoting
    for (int i=0; i<3; i++) {
        // Find pivot: the row with the maximum absolute value in column i
        int maxRow = i;
        for (int k=i+1; k<3; k++) {
            if (fabs(M[k][i]) > fabs(M[maxRow][i])) {
                maxRow = k;
            }
        }
        // If the pivot element is nearly 0, the matrix is singular
        if (fabs(M[maxRow][i]) < TOL) {
            return 0; // no unique solution
        }
        // Swap the current row with the pivot row if necessary
        if (maxRow != i) {
            for (int j=i; j<4; j++) {
                double tmp = M[i][j];
                M[i][j] = M[maxRow][j];
                M[maxRow][j] = tmp;
            }
        }
        // Eliminate entries below the pivot
        for (int k=i+1; k<3; k++) {
            double factor = M[k][i] / M[i][i];
            for (int j=i; j<3; j++) {
                M[k][j] -= factor * M[i][j];
            }
        }
    }
    
    // Back substitution to solve for x
    for (int i=2; i>=0; i--) {
        if (fabs(M[i][i]) < TOL) {
            return 0; // Should not occur here because we checked during elimination
        }
        x[i] = M[i][3] / M[i][i];
        for (int k = i - 1; k >= 0; k--) {
            M[k][3] -= M[k][i] * x[i];
        }
    }
    return 1; // solution found
}


// int PLU_inplace(double **A, int N, double tol, int *P, int *num_permut) 
// {
//     *num_permut = 0;
//     for (int ii=0; ii<N; ii++) {
//         P[ii] = ii;  // Unit permutation vector
//     }
//
//     for (int ii=0; ii<N; ii++) {
//         double maxA = 0;
//         int idmax = ii;
//
//         for (int jj=ii; jj<N; jj++) {  // Find largest value in row
//             double absA = fabs(A[jj][ii]);
//             if (absA > maxA) {
//                 maxA = absA;
//                 idmax = jj;
//             }
//         }
//         if (maxA < tol) return 1;  // Failure, A degenerate
//         if (idmax != ii) {
//             // Pivoting P
//             int aux = P[ii];
//             P[ii] = P[idmax];
//             P[idmax] = aux;
//
//             // Pivoting rows of A
//             double *row = A[ii];
//             A[ii] = A[idmax];
//             A[idmax] = row;
//
//             (*num_permut)++;
//         }
//         
//         // Update A in-place
//         for (int jj=ii+1; jj<N; jj++) {
//             A[jj][ii] /= A[ii][ii];
//
//             for (int kk=ii+1; kk<N; kk++)
//                 A[jj][kk] -= A[jj][ii] * A[ii][kk];
//         }
//     }
//     return 0; 
// }
//
//
// double determinant_PLU(double **A, int N) {
//     int P[N];  // Permutation vector
//     int num_permut;  // number of permutations
//     double tol = 1e-12;
//     int out = PLU_inplace(A, N, tol, P, &num_permut);
//     if (out == 1) {  // Matrix is singular
//         return 0;
//     }
//
//     double det = A[0][0];
//     for (int ii=1; ii<N; ii++) {
//         det *= A[ii][ii];
//     }
//     return num_permut % 2 == 0 ? det : -det;
// }
//
//
// double determinant_3x3(double **A)
// {
//     double det;
//     det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
//         - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
//         + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
//         return det;
// }
//
//
// double determinant(double **A, int N) {
//     switch(N) {  // In small dimensions, no need for PLU!
//         case 1:
//             return **A;
//         case 2:
//             return A[0][0] * A[1][1] - A[1][0] * A[0][1];
//         case 3:
//             return determinant_3x3(A);
//         default:
//             return determinant_PLU(A, N);  // A is modified!
//     }
// }

#endif
