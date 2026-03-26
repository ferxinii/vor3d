/* 
 * Header-only Linear Programming helpers.
 * 
 * Copyright (c) 2026 Fernando Muñoz
 * MIT license. See bottom of file.
 */

#ifndef HLIBS_LINEAR_PROGRAMMING_H
#define HLIBS_LINEAR_PROGRAMMING_H

#include <stdbool.h>
#include <math.h>


/* INTERFACE */
typedef struct LPconstraint2D {  
    double a, b, c;  /* Constraint: a*x + b*y <= c */
} s_LPconstraint2D;

/* constraints[N] is mutated in-place. x_guess and out_x can be NULL */
static inline bool LP_is_feasible_2D(int (*randint)(void* rctx, int), void* rctx, 
                                     int N, s_LPconstraint2D constraints[N], double TOL,
                                     double x_guess[2], double out_x[2]);



/* IMPLEMENTATION */
static inline bool constraint_valid(s_LPconstraint2D con, double TOL2)
{
    double n2 = con.a*con.a + con.b*con.b;
    if (n2 <= TOL2) return false;
    else return true;
}

static inline bool point_satisfies_constraint(s_LPconstraint2D con, const double p[2], double TOL) 
{   /* a*x + b*y <= c + TOL */
    double val = con.a * p[0] + con.b * p[1] - con.c;
    return (val <= TOL);
}

static inline void line_from_constraint(s_LPconstraint2D con, double out_p[2], double out_v[2]) 
{   /* Returns point on line p and direction v */
    out_v[0] = -con.b;  out_v[1] = con.a; 
    /* To find p, set either x=0 or y=0 */
    if (fabs(con.b) > fabs(con.a)) { out_p[0] = 0;  out_p[1] = con.c/con.b; }
    else { out_p[0] = con.c/con.a; out_p[1] = 0; }
}


static inline void shuffle_constraints_2D(int (*randint)(void*, int), void* rctx, int N, s_LPconstraint2D constraints[N])
{
	for (int i = N - 1; i > 0; --i) {
		int j = randint(rctx, i+1);
		s_LPconstraint2D tmp = constraints[i];
		constraints[i] = constraints[j];
		constraints[j] = tmp;
	}
}

static inline bool LP_is_feasible_2D(int (*randint)(void* rctx, int), void* rctx,
                                     int N, s_LPconstraint2D constraints[N], double TOL,
                                     double x_guess[2], double out_x[2]) 
{   /* Seidel algorithm */
    const double TOL2 = TOL * TOL;
    shuffle_constraints_2D(randint, rctx, N, constraints);
    double x[2]; if (x_guess) {x[0]=x_guess[0]; x[1]=x_guess[1];} else {x[0]=0; x[1]=0;}; 

    for (int ii=0; ii<N; ii++) {  
        if (!constraint_valid(constraints[ii], TOL2)) {
            if (constraints[ii].c < -TOL) return false;
            else continue;
        }
        if (point_satisfies_constraint(constraints[ii], x, TOL)) continue;

        /* Find new x on boundary of ii satisfiying all previous constraints */
        double p[2], v[2];
        line_from_constraint(constraints[ii], p, v);
        double norm_v = sqrt(v[0]*v[0] + v[1]*v[1]);

        double t_min = -INFINITY;
        double t_max = +INFINITY;
        double TOL_t = (norm_v > 0.0) ? (TOL / norm_v) : TOL;
        for (int jj=0; jj<ii; jj++) {
            double alpha = constraints[jj].a * v[0] + constraints[jj].b * v[1];
            double beta  = constraints[jj].c - (constraints[jj].a * p[0] + constraints[jj].b * p[1]);
            beta += TOL;  /* To be consistent with inclusion test */

            if (fabs(alpha) <= TOL * norm_v) {  /* alpha ~ 0 */
                if (beta < -TOL) return false;  /* Impossible */
                else continue;  /* jj places no restriction on this line */
            }

            double t_bound = beta / alpha;
            if (alpha > 0.0) {  /* t <= t_bound */
                if (t_bound < t_max) t_max = t_bound;
            } else {  /* alpha < 0 -> t >= t_bound */
                if (t_bound > t_min) t_min = t_bound;
            }

            /* check emptiness with t-tolerance mapped from coordinate TOL */
            if (t_min > t_max + TOL_t) return false;
        } 

        double t;
        if (isinf(t_min) && isinf(t_max)) t = 0.0;
        else if (isinf(t_min)) t = t_max - 1.0;
        else if (isinf(t_max)) t = t_min + 1.0;
        else t = 0.5 * (t_min + t_max);

        x[0] = p[0] + t * v[0];
        x[1] = p[1] + t * v[1];
    } 
    if (out_x) { out_x[0] = x[0];  out_x[1] = x[1]; }
    return true;
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

