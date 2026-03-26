/* 
 * Header-only library to solve systems of equations.
 *
 * Copyright (c) 2026 Fernando Muñoz
 * MIT License. See bottom of file.
 */

#ifndef HLIBS_LINALG_H
#define HLIBS_LINALG_H

#include <math.h>
#include <stdlib.h>

/* INTERFACE */
static inline int solve_2x2_cramer(const double A[2][2], const double b[2], double x[2], double pivot_tol);

static inline int solve_3x3_ppivot(const double A[3][3], const double b[3], double x[3], double pivot_tol);
static inline int solve_3x3_ppivot_inplace(double A[3][3], double b[3], double x[3], double pivot_tol);

static inline int solve_NxN_ppivot(int N, const double *A, const double *b, double *x, double pivot_tol);
static inline int solve_NxN_ppivot_inplace(int N, double *A, double *b, double *x, double pivot_tol);



/* IMPLEMENTATION */
static inline int solve_2x2_cramer(const double A[2][2], const double b[2], double x[2], double pivot_tol)
{
    double det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
    if (fabs(det) <= pivot_tol) return 0;
    x[0] = (b[0]*A[1][1] - b[1]*A[0][1]) / det;
    x[1] = (A[0][0]*b[1] - A[1][0]*b[0]) / det;
    return 2;
}


static inline int solve_3x3_ppivot_inplace(double A[3][3], double b[3], double x[3], double pivot_tol)
{
    for (int col = 0; col < 3; col++) {
        /* find pivot */
        int piv = col;
        double best = fabs(A[col][col]);
        for (int r = col + 1; r < 3; r++) {
            double val = fabs(A[r][col]);
            if (val > best) { best = val; piv = r; }
        }
        if (best <= pivot_tol) return col;
        if (piv != col) {
            for (int c = col; c < 3; c++) {
                double t = A[col][c]; A[col][c] = A[piv][c]; A[piv][c] = t;
            }
            double t = b[col]; b[col] = b[piv]; b[piv] = t;
        }
        /* eliminate */
        for (int r = col + 1; r < 3; r++) {
            double fac = A[r][col] / A[col][col];
            for (int c = col; c < 3; c++) A[r][c] -= fac * A[col][c];
            b[r] -= fac * b[col];
        }
    }
    /* back-substitute */
    for (int i = 2; i >= 0; i--) {
        double s = b[i];
        for (int j = i + 1; j < 3; j++) s -= A[i][j] * x[j];
        x[i] = s / A[i][i];
    }
    return 3;
}

static inline int solve_3x3_ppivot(const double A[3][3], const double b[3], double x[3], double pivot_tol)
{
    double Ac[3][3] = {{A[0][0], A[0][1], A[0][2]},
                       {A[1][0], A[1][1], A[1][2]},
                       {A[2][0], A[2][1], A[2][2]}};
    double bc[3] = {b[0], b[1], b[2]};
    return solve_3x3_ppivot_inplace(Ac, bc, x, pivot_tol);
}


static inline int solve_NxN_ppivot_inplace(int N, double *A, double *b, double *x, double pivot_tol)
{
    for (int col = 0; col < N; col++) {
        /* find pivot */
        int piv = col;
        double best = fabs(A[col*N+col]);
        for (int r = col + 1; r < N; r++) {
            double val = fabs(A[r*N+col]);
            if (val > best) { best = val; piv = r; }
        }
        if (best <= pivot_tol) return col;
        if (piv != col) {
            for (int c = col; c < N; c++) {
                double t = A[col*N+c]; A[col*N+c] = A[piv*N+c]; A[piv*N+c] = t;
            }
            double t = b[col]; b[col] = b[piv]; b[piv] = t;
        }
        /* eliminate */
        for (int r = col + 1; r < N; r++) {
            double fac = A[r*N+col] / A[col*N+col];
            for (int c = col; c < N; c++) A[r*N+c] -= fac * A[col*N+c];
            b[r] -= fac * b[col];
        }
    }
    /* back-substitute */
    for (int i = N - 1; i >= 0; i--) {
        double s = b[i];
        for (int j = i + 1; j < N; j++) s -= A[i*N+j] * x[j];
        x[i] = s / A[i*N+i];
    }
    return N;
}

static inline int solve_NxN_ppivot(int N, const double *A, const double *b, double *x, double pivot_tol)
{
    #define STACK_MAX 8
    double stack_A[STACK_MAX * STACK_MAX];
    double stack_b[STACK_MAX];
    double *Ac = (N <= STACK_MAX) ? stack_A : (double *)malloc(N * N * sizeof(double));
    double *bc = (N <= STACK_MAX) ? stack_b : (double *)malloc(N * sizeof(double));
    if (!Ac || !bc) { 
        if (N > STACK_MAX) { free(Ac); free(bc); }
        return -1; 
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) Ac[i*N+j] = A[i*N+j];
        bc[i] = b[i];
    }
    int rank = solve_NxN_ppivot_inplace(N, Ac, bc, x, pivot_tol);
    if (N > STACK_MAX) { free(Ac); free(bc); }
    #undef STACK_MAX
    return rank;
}

#endif


/* MIT License.
 *
 * Copyright (c) 2026 Fernando Muñoz.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

