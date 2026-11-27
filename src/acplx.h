#ifndef VOR3D_ACPLX_H
#define VOR3D_ACPLX_H

/* Library-internal alpha-complex marking engine (src/acplx.c).
 *
 * Works in place on a (weighted) Delaunay triangulation s_scplx.  Power
 * distance |x-p|^2 - w; alpha is a SQUARED (power) radius.  A simplex is in the
 * alpha complex K iff its orthoball has squared radius <= alpha and is
 * empty/attached-free (Edelsbrunner).  The public API (include/alpha.h) wraps
 * this; volsph.c uses it directly because it needs per-simplex lookup by DT
 * vertex id over ALL DT tets. */

#include <stdbool.h>
#include "scplx.h"
#include "./hlibs/hash.h"
#include "./hlibs/dynarray.h"

/* Tet predicate: true iff the tet's orthoball has squared radius <= alpha
 * (tets of a regular triangulation are always conflict-free).  With
 * has_big_tetra, sentinel-touching tets (any vertex id < 4) are never in. */
bool acplx_tet_in(const s_scplx *scplx, const s_ncell *nc, double alpha, bool has_big_tetra);

/* Mark the alpha complex K of scplx.
 *   - nc->mask_alpha = true for every tet in K (left untouched otherwise, so
 *     the caller must clear it first if the complex was marked before).
 *   - out_faces: key int[3] (sorted DT vertex ids), value bool: true iff the
 *     face is in K.  Every face of the DT is present.
 *   - out_edges: key int[2] (sorted), value bool, likewise for every DT edge.
 *   - out_vertices (length scplx->points.N, caller-allocated): true iff the
 *     vertex is in K.
 * buff_ncellPTR: caller-provided scratch dynarray of s_ncell* (vertex star walk).
 * Returns 1 on success, 0 on allocation failure; on 0 nothing is left
 * allocated in out_faces / out_edges and out_vertices is unspecified.  On 1 the
 * caller frees both tables with hash_free. */
int alpha_complex_mark(s_scplx *scplx, bool has_big_tetra, double alpha, s_dynarray *buff_ncellPTR,
                       s_hash_table *out_faces, s_hash_table *out_edges,
                       bool *out_vertices);

#endif
