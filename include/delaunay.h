#ifndef VOR3D_DELAUNAY
#define VOR3D_DELAUNAY

#include "scplx.h"

s_scplx construct_dt_3d(const s_points *points, const double *weights, bool keep_big_tetra, double TOL_duplicates);  /* Copy of points inside s_scplx */

int N_ncells_not_big_tetra(s_scplx *scplx);  /* Only makes sense if scplx has kept big_tetra */


void extract_alpha_complex(s_scplx *scplx, bool has_big_tetra, double alpha, 
                           double EPS_DEGEN, s_dynarray *buff_ncellPTR,
                           s_hash_table *out_faces, s_hash_table *out_edges, 
                           bool *out_vertices);

#endif
