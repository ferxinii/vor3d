/* 
 * Header-only PRNG implementation.
 * Allows sampling from various distributions in a robust, unbiased, and fast
 * manner. Support for multiple threads (each has a non-overlapping context 
 * derived from the same seed). Internally uses XOSHIRO256** PRNG. 
 * 
 * Copyright (c) 2026 Fernando Muñoz
 * MIT license. See bottom of file.
 *
 * TODO: Random sample from a given cdf
 */

#ifndef HLIBS_RANDOM_H
#define HLIBS_RANDOM_H
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>



typedef struct random_context {  /* Ignore contents, implementation details */
    uint64_t XOSHIRO256_state[4];      

    /* Normal distribution algorithm generates two samples each time, so we cache them */
    bool stored_standard_normal;
    double next_standard_normal;
} s_random_context;



/* INTERFACE */
static inline s_random_context random_initialize(uint64_t seed);
static inline void random_initialize_threads(uint64_t seed, int N_MPI, int N_OMP, s_random_context out[N_MPI*N_OMP]);
static inline uint64_t random_uniform_range_u64(s_random_context *ctx, uint64_t N);  /* [0,N) */
static inline double random_uniform_double(s_random_context *ctx);  /* [0,1) */
static inline double random_normal(s_random_context *ctx, double mean, double std);
static inline int random_poisson(s_random_context *ctx, double lambda);
static inline void random_shuffle(s_random_context *ctx, int N, int out[N]);
static inline void random_pdf_to_cdf(int N, const double pdf[N], double cdf[N]);  /* Can be used in-place */
static inline int random_sample_cdf(s_random_context *ctx, int N, const double cdf[N]);  /* No need to be normalised */




/* IMPLEMENTATION */
/* Internal PRNG: xoshiro256** */
static inline uint64_t XOSHIRO256_splitmix64(uint64_t *x) {   /* Modifies x */
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static inline uint64_t XOSHIRO256_uint64_rotl(const uint64_t x, int k) 
{   /* Rotate left: move bits to the left k spots and wrap around */
    return (x << k) | (x >> (64 - k));
}

static inline uint64_t XOSHIRO256_next(s_random_context *ctx) 
{
    uint64_t *s = ctx->XOSHIRO256_state;

    const uint64_t result = XOSHIRO256_uint64_rotl(s[1] * 5, 7) * 9;
    const uint64_t t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];

    s[2] ^= t;
    s[3] = XOSHIRO256_uint64_rotl(s[3], 45);
    return result;
}

/* The following are used to jump ahead 2^128 or 2^192 calls to next. 
 * Use when giving context to different threads to generate different 
 * starting points on each, but with the same seed.
 * */
static inline void XOSHIRO256_jump(s_random_context *ctx) {
    static const uint64_t JUMP[4] = { 0x180ec6d33cfd0abaULL,
                                      0xd5a61266f0c9392cULL,
                                      0xa9582618e03fc9aaULL,
                                      0x39abdc4529b1661cULL };
    uint64_t t[4] = {0,0,0,0};
    for(int i=0; i<4; i++) {
        for(int b=0; b<64; b++) if (JUMP[i] & (1ULL<<b)) {
            for(int j=0; j<4; j++) t[j] ^= ctx->XOSHIRO256_state[j];
        }
        XOSHIRO256_next(ctx);
    }

    for (int j=0; j<4; j++) ctx->XOSHIRO256_state[j] = t[j];
}

static inline void XOSHIRO256_long_jump(s_random_context *ctx)
{
    static const uint64_t LONG_JUMP[4] = { 0x76e15d3efefdcbbfULL,
                                           0xc5004e441c522fb3ULL,
                                           0x77710069854ee241ULL,
                                           0x39109bb02acbe635ULL };
    uint64_t t[4] = {0, 0, 0, 0};
    for (int i=0; i<4; i++) {
        for (int b=0; b<64; b++) if (LONG_JUMP[i] & (1ULL << b)) {
            for (int j=0; j<4; j++) t[j] ^= ctx->XOSHIRO256_state[j];
        }
        XOSHIRO256_next(ctx);
    }

    for (int j=0; j<4; j++) ctx->XOSHIRO256_state[j] = t[j];
}




static inline s_random_context random_initialize(uint64_t seed)
{
    s_random_context out;

    /* Initialize xoshiro256 state */
    uint64_t x = seed;  /* gets modified */
    for (int i = 0; i < 4; ++i) {
        out.XOSHIRO256_state[i] = XOSHIRO256_splitmix64(&x);
    }

    out.stored_standard_normal = false;
    out.next_standard_normal = 0.0;
    return out;
}

static inline void random_initialize_threads(uint64_t seed, int N_MPI, int N_OMP, s_random_context out[N_MPI*N_OMP])
{
    assert(N_MPI > 0 && N_OMP > 0);

    out[0] = random_initialize(seed);  /* Main context */
    
    for (int i=0; i<N_MPI; i++) {
        /* Separate MPI processes by long jumps (2^192) */
        if (i != 0) {
            out[i*N_OMP] = out[(i-1)*N_OMP];
            XOSHIRO256_long_jump(&out[i*N_OMP]);
        }
        for (int j=1; j<N_OMP; j++) {
            /* Separate OMP processed by "short" jump (2^128) */
            out[i*N_OMP + j] = out[i*N_OMP + (j-1)];
            XOSHIRO256_jump(&out[i*N_OMP + j]);
        }
    }
}

static inline uint64_t random_uniform_u64(s_random_context *ctx)
{   /* Main interface to XOSHIROT256. [0, 2^64) */
    return XOSHIRO256_next(ctx);
}


static inline uint64_t random_uniform_range_u64(s_random_context *ctx, uint64_t N)
{   /* Lemire algorith. Random int in [0, N) */
    if (N <= 1) return 0;  /* If N==1, the only possible output is 0 */
    
    uint64_t t = -N % N ;
    uint64_t x = random_uniform_u64(ctx);
    __uint128_t m = (__uint128_t)x * (__uint128_t)N;
    uint64_t l = (uint64_t)m;
    if (l < N) {
        while (l < t) {
            x = random_uniform_u64(ctx);
            m = (__uint128_t)x * (__uint128_t)N;
            l = ( uint64_t )m;
        }
    }
    return (uint64_t)(m >> 64);
}


static inline double random_uniform_double(s_random_context *ctx) 
{   /* [0, 1) */
    uint64_t x = random_uniform_u64(ctx);
    const uint64_t top53 = x >> 11;  /* top 53 bits, uniform */
    return (double)top53 * (1.0 / 9007199254740992.0);  /* divide by 2^53 */
}


static inline double random_normal(s_random_context *ctx, double mean, double std)
{   /* Box-Muller algorithm */
    if (ctx->stored_standard_normal) {
        ctx->stored_standard_normal = false;
        return mean + std * ctx->next_standard_normal;
    }

    double U, V;
    do {  /* Get U strictly > 0 */
        U = random_uniform_double(ctx); 
    } while (U <= 0.0);
    V = random_uniform_double(ctx);

    double mag = sqrt(-2.0 * log(U));
    ctx->next_standard_normal = mag * sin(2.0 * M_PI * V);
    ctx->stored_standard_normal = true;
    return mean + std * mag * cos(2.0 * M_PI * V);
}


static inline double log1pexp(double y) 
{   /* Stable log(1+exp(y)) */
    if (y > 0) return y + log1p(exp(-y));
    return log1p(exp(y));
}

static inline int random_poisson_KNUTH(s_random_context *ctx, double lambda)
{   /* Knuth algorithm */
    double L = exp(-lambda);
    double p = 1.0;
    int k = 0;
    do {
        k++;
        p *= random_uniform_double(ctx);
    } while (p > L);
    return k - 1;
}

static inline int random_poisson_ATKINSON(s_random_context *ctx, double lambda)
{   /* Algorithm PA by Atkinson */
    const double c = 0.767 - 3.36 / lambda;
    const double beta = M_PI / sqrt(3.0 * lambda);
    const double alpha = beta * lambda;
    const double k = log(c) - lambda - log(beta);
    const double log_lambda = log(lambda);
    for (;;) {
        double U = random_uniform_double(ctx);
        if (U <= 0.0 || U >= 1.0) continue;

        double x = (alpha - log((1.0 - U) / U)) / beta;
        int N = floor(x + 0.5);
        if (N < 0) continue;

        double V = random_uniform_double(ctx);
        if (V <= 0.0) continue;

        double y = alpha - beta * x;

        double lhs = y + log(V) - 2.0 * log1pexp(y);  /* y + log(v/(1+exp(y))^2) */
        double rhs = k + N * log_lambda - lgamma((double)N+1.0);
        if (lhs <= rhs) return N;
    }
}

static inline int random_poisson(s_random_context *ctx, double lambda) 
{
    if (!(lambda >= 0.0)) return 0;  

    if (lambda <= 30.0) return random_poisson_KNUTH(ctx, lambda);
    else return random_poisson_ATKINSON(ctx, lambda);
}


static inline void random_shuffle(s_random_context *ctx, int N, int out[N])
{  /* Fisher-Yates algorithm */
    for (int i=0; i<N; i++) out[i] = i;  /* Initialize */

    for (int i = N-1; i>0; i--) {
		int j = random_uniform_range_u64(ctx, i+1);
	    int tmp = out[i];
		out[i] = out[j];
		out[j] = tmp;
	}
}


static inline void random_pdf_to_cdf(int N, const double pdf[N], double cdf[N])
{   /* Can be used in-place */
    cdf[0] = pdf[0];
    for (int i = 1; i < N; i++) {
        cdf[i] = cdf[i-1] + pdf[i];
    }
}


static inline int random_sample_cdf(s_random_context *ctx, int N, const double cdf[N])
{   /* Using binary search, no need to be normalised. */
    double r = random_uniform_double(ctx) * cdf[N-1];  /* Uniform in [0, F(N-1)) */

    int lo = 0, hi = N - 1;
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (r < cdf[mid]) hi = mid;
        else lo = mid + 1;
    }

    return lo;
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

