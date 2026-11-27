#ifndef VOR3D_OPTS_H
#define VOR3D_OPTS_H

/* Per-call options. opts == NULL and &(s_vor3d_opts){0} both mean all defaults.
 * Per-thread: rctx is mutated during a call. */
typedef struct vor3d_opts {
    double      tol_dup;           /* points closer than this (absolute length) are
                                      duplicates; 0 = exact duplicates only */
    double    (*randd01)(void *);  /* RNG in [0,1), randomizes point-location walks
                                      (start + face order); NULL = fixed order */
    void       *rctx;              /* passed to randd01 */
} s_vor3d_opts;

#endif
