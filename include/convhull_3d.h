/*
 Copyright (c) 2017-2021 Leo McCormack
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
*/
/*
 * Filename:
 *     convhull_3d.h
 * Description:
 *     A header only C implementation of the 3-D quickhull algorithm.
 *     The code is largely derived from the "computational-geometry-toolbox"
 *     by George Papazafeiropoulos (c) 2014, originally distributed under
 *     the BSD (2-clause) license.
 *     To include this implementation in a project, simply add this:
 *         #define CONVHULL_3D_ENABLE
 *         #include "convhull_3d.h"
 *     By default, the algorithm uses double floating point precision. To
 *     use single precision (less accurate but quicker), also add this:
 *         #define CONVHULL_3D_USE_SINGLE_PRECISION
 *     If your project has CBLAS linked, then you can also speed things up
 *     a tad by adding this:
 *         #define CONVHULL_3D_USE_CBLAS
 *     The code is C++ compiler safe.
 *     Reference: "The Quickhull Algorithm for Convex Hull, C. Bradford
 *                 Barber, David P. Dobkin and Hannu Huhdanpaa, Geometry
 *                 Center Technical Report GCG53, July 30, 1993"
 * Dependencies:
 *     cblas (optional for speed ups, especially for very large meshes)
 *     (Available in e.g. Apple Accelerate Framework, or Intel MKL)
 * Author, date created:
 *     Leo McCormack, 02.10.2017
 */


#ifndef CONVHULL_3D_H
#define CONVHULL_3D_H

typedef double CH_FLOAT;

typedef struct _ch_vertex {
    union {
        CH_FLOAT v[3];
        struct{
             CH_FLOAT x, y, z;
        };
    };
} ch_vertex;
typedef ch_vertex ch_vec3;

/* builds the 3-D convexhull */
void convhull_3d_build(/* input arguments */
                       ch_vertex* const in_vertices,            /* vector of input vertices; nVert x 1 */
                       const int nVert,                         /* number of vertices */
                       /* output arguments */
                       int** out_faces,                         /* & of empty int*, output face indices; flat: nOut_faces x 3 */
                       int* nOut_faces);                        /* & of int, number of output face indices */
    
/* exports the vertices, face indices, and face normals, as an 'obj' file, ready for GPU (for 3d convexhulls only) */
void convhull_3d_export_obj(/* input arguments */
                            ch_vertex* const vertices,          /* vector of input vertices; nVert x 1 */
                            const int nVert,                    /* number of vertices */
                            int* const faces,                   /* face indices; flat: nFaces x 3 */
                            const int nFaces,                   /* number of faces in hull */
                            const int keepOnlyUsedVerticesFLAG, /* 0: exports in_vertices, 1: exports only used vertices  */
                            char* const obj_filename);          /* obj filename, WITHOUT extension */
    
/* exports the vertices, face indices, and face normals, as an 'm' file, for MatLab verification (for 3d convexhulls only) */
void convhull_3d_export_m(/* input arguments */
                          ch_vertex* const vertices,            /* vector of input vertices; nVert x 1 */
                          const int nVert,                      /* number of vertices */
                          int* const faces,                     /* face indices; flat: nFaces x 3 */
                          const int nFaces,                     /* number of faces in hull */
                          char* const m_filename);              /* m filename, WITHOUT extension */
    
/* reads an 'obj' file and extracts only the vertices (for 3d convexhulls only) */
void extract_vertices_from_obj_file(/* input arguments */
                                    char* const obj_filename,       /* obj filename, WITHOUT extension */
                                    /* output arguments */
                                    ch_vertex** out_vertices,       /* & of empty ch_vertex*, output vertices; out_nVert x 1 */
                                    int* out_nVert);                /* & of int, number of vertices */

/**** NEW! ****/

/* builds the N-Dimensional convexhull of a grid of points */
void convhull_nd_build(/* input arguments */
                       CH_FLOAT* const in_points,               /* Matrix of points in 'd' dimensions; FLAT: nPoints x d */
                       const int nPoints,                       /* number of points */
                       const int d,                             /* Number of dimensions */
                       /* output arguments */
                       int** out_faces,                         /* (&) output face indices; FLAT: nOut_faces x d */
                       CH_FLOAT** out_cf,                       /* (&) contains the coefficients of the planes (set to NULL if not wanted); FLAT: nOut_faces x d */
                       CH_FLOAT** out_df,                       /* (&) contains the constant terms of the planes (set to NULL if not wanted); nOut_faces x 1 */
                       int* nOut_faces);                        /* (&) number of output face indices */

/* Computes the Delaunay triangulation (mesh) of an arrangement of points in N-dimensional space */
void delaunay_nd_mesh(/* input Arguments */
                      const float* points,                      /* The input points; FLAT: nPoints x nd */
                      const int nPoints,                        /* Number of points */
                      const int nd,                             /* The number of dimensions */
                      /* output Arguments */
                      int** Mesh,                               /* (&) the indices defining the Delaunay triangulation of the points; FLAT: nMesh x (nd+1) */
                      int* nMesh);                              /* (&) Number of triangulations */

/**** CUSTOM ALLOCATOR VERSIONS ****/

/* builds the 3-D convexhull */
void convhull_3d_build_alloc(/* input arguments */
                             ch_vertex* const in_vertices,            /* vector of input vertices; nVert x 1 */
                             const int nVert,                         /* number of vertices */
                             /* output arguments */
                             int** out_faces,                         /* & of empty int*, output face indices; flat: nOut_faces x 3 */
                             int* nOut_faces,                         /* & of int, number of output face indices */
                             void* allocator);                        /* & of an allocator */

/* builds the N-Dimensional convexhull of a grid of points */
void convhull_nd_build_alloc(/* input arguments */
                             CH_FLOAT* const in_points,               /* Matrix of points in 'd' dimensions; FLAT: nPoints x d */
                             const int nPoints,                       /* number of points */
                             const int d,                             /* Number of dimensions */
                             /* output arguments */
                             int** out_faces,                         /* (&) output face indices; FLAT: nOut_faces x d */
                             CH_FLOAT** out_cf,                       /* (&) contains the coefficients of the planes (set to NULL if not wanted); FLAT: nOut_faces x d */
                             CH_FLOAT** out_df,                       /* (&) contains the constant terms of the planes (set to NULL if not wanted); nOut_faces x 1 */
                             int* nOut_faces,                         /* (&) number of output face indices */
                             void* allocator);                        /* & of an allocator */

/* Computes the Delaunay triangulation (mesh) of an arrangement of points in N-dimensional space */
void delaunay_nd_mesh_alloc(/* input Arguments */
                            const float* points,                      /* The input points; FLAT: nPoints x nd */
                            const int nPoints,                        /* Number of points */
                            const int nd,                             /* The number of dimensions */
                            /* output Arguments */
                            int** Mesh,                               /* (&) the indices defining the Delaunay triangulation of the points; FLAT: nMesh x (nd+1) */
                            int* nMesh,                               /* (&) Number of triangulations */
                            void* allocator);                         /* & of an allocator */

/* reads an 'obj' file and extracts only the vertices (for 3d convexhulls only) */
void extract_vertices_from_obj_file_alloc(/* input arguments */
                                          char* const obj_filename,       /* obj filename, WITHOUT extension */
                                          /* output arguments */
                                          ch_vertex** out_vertices,       /* & of empty ch_vertex*, output vertices; out_nVert x 1 */
                                          int* out_nVert,                 /* & of int, number of vertices */
                                          void* allocator);               /* & of an allocator */


#endif 

