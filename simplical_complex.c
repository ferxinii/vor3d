#ifndef SIMPLICAL_COMPLEX_C
#define SIMPLICAL_COMPLEX_C

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "geometry.c"
#include "array_operations.c"


typedef struct setup {
    int dim;
    int N_points;
    double **points;  // (n_points + dim + 1) x dim
                      // The last (dim+1) points correspond to the big n_cell
    int N_ncells;
    struct ncell *head;  // Linked list of ncells
} s_setup;


typedef struct ncell {
    int *vertex_id;
    struct ncell **opposite;
    struct ncell *next;  // Linked list of cells
    struct ncell *prev;
    int mark;  // Used to mark particular ncells
    int count;
    double volume;  // DEBUGGING
} s_ncell;


s_ncell *malloc_ncell(const s_setup *setup)
{
    s_ncell *out = malloc(sizeof(s_ncell));
    out->vertex_id = malloc(sizeof(int) * (setup->dim + 1));
    out->opposite = malloc(sizeof(s_ncell*) * (setup->dim + 1));
    out->next = NULL;
    out->prev = NULL;
    out->mark = 0;
    assert(out->vertex_id && out->opposite && "Could not malloc ncell!");
    return out;
}


void free_ncell(s_ncell *ncell)
{
    free(ncell->vertex_id);
    free(ncell->opposite);
    free(ncell);
}


void free_complex(s_setup *setup)
{   
    s_ncell *current = setup->head;
    s_ncell *next = current->next;
    free_ncell(current);
    for (int ii=1; ii<setup->N_ncells; ii++) {
        current = next;
        next = current->next;
        free_ncell(current);
    }

    free_matrix(setup->points, setup->N_points);
    free(setup);
}


void print_ncell(const s_setup *setup, const s_ncell *ncell)
{
    printf("%p, ( ", (void*)ncell);
    for (int ii=0; ii<setup->dim+1; ii++) {
        printf("%d ", ncell->vertex_id[ii]);
    }
    printf(") \n");
}


void print_ncells(const s_setup *setup)
{
    puts("NCELLS");
    s_ncell *current = setup->head;
    int ii = 0;
    while (current) {
        printf("%d  |  p : %p  |  prev : %p  |  marked : %d  |  vertex_ids :", ii, (void*)current, (void*)current->prev, current->mark);
        for (int jj=0; jj<setup->dim+1; jj++) {
            printf(" %d", current->vertex_id[jj]);
        }
        printf("  |  opposite :");
        for (int jj=0; jj<setup->dim+1; jj++) {
            printf(" %p", (void*)current->opposite[jj]);
        }
        printf("\n");
        ii++;
        current = current->next;
    }
    puts("");
}


void write_ncell3d_file(s_setup *setup, s_ncell *ncell, FILE *file)
{
    fprintf(file, "%f %f %f\n", setup->points[ncell->vertex_id[0]][0],
                                setup->points[ncell->vertex_id[0]][1],
                                setup->points[ncell->vertex_id[0]][2]);
    fprintf(file, "%f %f %f\n", setup->points[ncell->vertex_id[1]][0],
                                setup->points[ncell->vertex_id[1]][1],
                                setup->points[ncell->vertex_id[1]][2]);
    fprintf(file, "%f %f %f\n\n", setup->points[ncell->vertex_id[2]][0],
                                  setup->points[ncell->vertex_id[2]][1],
                                  setup->points[ncell->vertex_id[2]][2]);

    fprintf(file, "%f %f %f\n", setup->points[ncell->vertex_id[0]][0],
                                setup->points[ncell->vertex_id[0]][1],
                                setup->points[ncell->vertex_id[0]][2]);
    fprintf(file, "%f %f %f\n", setup->points[ncell->vertex_id[1]][0],
                                setup->points[ncell->vertex_id[1]][1],
                                setup->points[ncell->vertex_id[1]][2]);
    fprintf(file, "%f %f %f\n\n", setup->points[ncell->vertex_id[3]][0],
                                  setup->points[ncell->vertex_id[3]][1],
                                  setup->points[ncell->vertex_id[3]][2]);

    fprintf(file, "%f %f %f\n", setup->points[ncell->vertex_id[3]][0],
                                setup->points[ncell->vertex_id[3]][1],
                                setup->points[ncell->vertex_id[3]][2]);
    fprintf(file, "%f %f %f\n", setup->points[ncell->vertex_id[1]][0],
                                setup->points[ncell->vertex_id[1]][1],
                                setup->points[ncell->vertex_id[1]][2]);
    fprintf(file, "%f %f %f\n\n", setup->points[ncell->vertex_id[2]][0],
                                  setup->points[ncell->vertex_id[2]][1],
                                  setup->points[ncell->vertex_id[2]][2]);

    fprintf(file, "%f %f %f\n", setup->points[ncell->vertex_id[0]][0],
                                setup->points[ncell->vertex_id[0]][1],
                                setup->points[ncell->vertex_id[0]][2]);
    fprintf(file, "%f %f %f\n", setup->points[ncell->vertex_id[3]][0],
                                setup->points[ncell->vertex_id[3]][1],
                                setup->points[ncell->vertex_id[3]][2]);
    fprintf(file, "%f %f %f\n\n\n", setup->points[ncell->vertex_id[2]][0],
                                  setup->points[ncell->vertex_id[2]][1],
                                  setup->points[ncell->vertex_id[2]][2]);
}


void write_dt3d_file(s_setup *setup, FILE *file)
{
    s_ncell *current = setup->head;
    while (current) {
        write_ncell3d_file(setup, current, file);
        fprintf(file, "\n\n");
        current = current->next;
    }
}



void initialize_ncells_counter(const s_setup *setup)
{
    s_ncell *current = setup->head;
    int ii = 0;
    while (current) {
        current->count = ii;
        current = current->next;
        ii++;
    }
}


void initialize_ncells_mark(const s_setup *setup)  // TODO CHECK IF THIS WORKS?
{  
    s_ncell *current = setup->head;
    for (int ii=0; ii<setup->N_ncells; ii++) {
        current->mark = 0;
        current = current->next;
    }
}


void print_marked(const s_setup *setup)
{
    puts("Marked node ids:");
    s_ncell *current = setup->head;
    for (int ii=0; ii<setup->N_ncells; ii++) {
        if (current->mark == 1) {
            printf("    %d\n", ii);
        }
        current = current->next;
    }
}


int count_marked(const s_setup *setup) 
{   
    int count = 0;
    s_ncell *current = setup->head;
    for (int ii=0; ii<setup->N_ncells; ii++) {
        if (current->mark == 1) {
            count++;
        }
        current = current->next;
    }
    return count;
}


void extract_vertices_ncell(const s_setup *setup, const s_ncell *ncell, double **out)
{
    for (int ii=0; ii<setup->dim+1; ii++) {
        for (int jj=0; jj<setup->dim; jj++) {
            out[ii][jj] = setup->points[ncell->vertex_id[ii]][jj];
        }
    }
}


void extract_ids_face(const s_setup *setup, const s_ncell *ncell, const int *v_localid, int dim_face, int *out)
{   // v_localid[dim-dim_face], out[dim_face+1]
    int kk = 0;
    for (int ii=0; ii<setup->dim+1; ii++) {
        if (!inarray(v_localid, setup->dim - dim_face, ii)) {
            out[kk] = ncell->vertex_id[ii];
            kk++;
        }
    }
}


void extract_vertices_face(const s_setup *setup, const s_ncell *ncell, const int *v_localid, int dim_face, double **out)
{   // v_localid[dim-dim_face], out[dim_face+1][dim]
    int face_vertex_id[dim_face+1];
    extract_ids_face(setup, ncell, v_localid, dim_face, face_vertex_id);
    for (int ii=0; ii<dim_face+1; ii++) {
        for (int jj=0; jj<setup->dim; jj++) {
            out[ii][jj] = setup->points[face_vertex_id[ii]][jj];
        }
    }
}


void extract_face_center_and_normal(const s_setup *setup, const s_ncell *ncell, int face_localid, double *fc, double *n)
{
    static double **face_v = NULL;
    if (!face_v) face_v = malloc_matrix(3, 3);
    static double **ncell_v = NULL;
    if (!ncell_v) ncell_v = malloc_matrix(4, 3);

    extract_vertices_face(setup, ncell, &face_localid, 2, face_v);
    extract_vertices_ncell(setup, ncell, ncell_v);

    double d1[3], d2[3];
    subtract_3d(face_v[1], face_v[0], d1);
    subtract_3d(face_v[2], face_v[0], d2);
    cross_3d(d1, d2, n);

    double cc[3], v[3];
    find_center_mass(face_v, 3, 3, fc);
    find_center_mass(ncell_v, 4, 3, cc);
    subtract_3d(fc, cc, v);

    double dir_aux = dot_3d(v, n);
    assert(dir_aux != 0 && "Vectors are perpendicular?");
    if (dir_aux < 0) {
        n[0] = -n[0];
        n[1] = -n[1];
        n[2] = -n[2];
    } 
}


void face_localid_of_adjacent_ncell(const s_setup *setup, const s_ncell *ncell, const int *v_localid,
                                         int dim_face, int id_adjacent, int *out_v_localid)
{
    s_ncell *adjacent = ncell->opposite[id_adjacent];

    int vertex_id_face[dim_face+1];
    extract_ids_face(setup, ncell, v_localid, dim_face, vertex_id_face);

    int kk = 0;
    for (int ii=0; ii<setup->dim+1; ii++) {
        if (!inarray(vertex_id_face, dim_face+1, adjacent->vertex_id[ii])) {
            out_v_localid[kk] = ii;
            kk++;
        }
    }
}


s_ncell *next_ncell_ridge_cycle(const s_setup *setup, const s_ncell *ncell, int v_localid_main, int v_localid_2, 
                                int *new_v_localid_main, int *new_v_localid_2)
{
    s_ncell *next = ncell->opposite[v_localid_main];
    *new_v_localid_main = id_where_equal_int(next->vertex_id, setup->dim+1, ncell->vertex_id[v_localid_2]);
    for (int ii=0; ii<setup->dim+1; ii++) {
        if (next->opposite[ii] == ncell) {
            *new_v_localid_2 = ii;
            return next;
        }
    }
    assert(1 == 0 && "Could not find next ncell around ridge."); 
    exit(1);
}


int count_cycle_ridge(const s_setup *setup, const s_ncell *ncell, int v_localid_main, int v_localid_2)  // local 0:dim
{   
    int counter = 1;
    int maxit = 10000;
    const s_ncell *current = ncell;
    while (counter < maxit) {
        int new_v_localid_main, new_v_localid_2;
        if (current->opposite[v_localid_main] == NULL) {
            return -1;
        }

        s_ncell *next = next_ncell_ridge_cycle(setup, current, v_localid_main, v_localid_2, 
                                              &new_v_localid_main, &new_v_localid_2);
        if (next == ncell) return counter;
        counter++;
        current = next;
        v_localid_main = new_v_localid_main;
        v_localid_2 = new_v_localid_2;
    }
    assert(1==0 && "Reached maximum iterations in ridge cycle?");
    exit(1);
}


void mark_ncells_incident_face_STEP(const s_setup *setup, const s_ncell *ncell, const int *v_localid, int dim_face)
{
    for (int ii=0; ii<setup->dim+1; ii++) {
        if (inarray(v_localid, setup->dim - dim_face, ii)) {
            s_ncell *adjacent_ncell = ncell->opposite[ii];  // This ncell should already share the face!
            if (adjacent_ncell && adjacent_ncell->mark == 0) {
                adjacent_ncell->mark = 1;
                
                int new_v_localid[setup->dim - dim_face];
                face_localid_of_adjacent_ncell(setup, ncell, v_localid, dim_face, ii, new_v_localid);

                // Recursion
                mark_ncells_incident_face_STEP(setup, adjacent_ncell, new_v_localid, dim_face);
            }
        }
    }
}


void mark_ncells_incident_face(const s_setup *setup, s_ncell *ncell, const int *v_localid, int dim_face)
{
    initialize_ncells_mark(setup);
    ncell->mark = 1;
    
    // Recursion:
    mark_ncells_incident_face_STEP(setup, ncell, v_localid, dim_face);
}


int are_locally_delaunay_strict(const s_setup *setup, const s_ncell *ncell, int id_opposite)
{   
    // Create array for coords1 in static memory, CANNOT BE MULTI-THREADED! FIXME
    static int prev_dim = 0;
    static double **coords1 = NULL, **coords2 = NULL;
    if (setup->dim != prev_dim) {
        if (coords1) free_matrix(coords1, prev_dim + 1);
        if (coords2) free_matrix(coords2, prev_dim + 1);
        prev_dim = setup->dim;
        coords1 = malloc_matrix(setup->dim + 1, setup->dim);
        coords2 = malloc_matrix(setup->dim + 1, setup->dim);
    }

    extract_vertices_ncell(setup, ncell, coords1);
    extract_vertices_ncell(setup, ncell->opposite[id_opposite], coords2);
            
    // Extract vertex_id of opposite's cell face
    int opp_face_localid;
    face_localid_of_adjacent_ncell(setup, ncell, &id_opposite, setup->dim-1, id_opposite, &opp_face_localid);
    int opp_face_vertex_id = (ncell->opposite[id_opposite])->vertex_id[opp_face_localid];
    
    int in1 = in_sphere(coords1, setup->points[opp_face_vertex_id], setup->dim);
    // in_sphere(coords2, setup->points[ncell->vertex_id[id_opposite]], setup->dim) == -1)
    // if (in1 == 0) puts("DEBUG: ARE_LOCALLY_DELAUNAY_STRICT IS CIRCUMSCRIBED");
    if (in1 == -1) return 1;
    else return 0;
}


// int pd_intersects_face(double *s1, double *s2, double *p, double *d, int drop_coord)
// {
//     int i1, i2;
//     if (drop_coord == 0) {
//         i1 = 1;     i2 = 2;
//     } else if (drop_coord == 1) {
//         i1 = 2;     i2 = 0;
//     } else {
//         i1 = 0;     i2 = 1;
//     }
//
//     double A[2], B[2], paux[2], daux[2];
//     A[0] = s1[i1];    A[1] = s1[i1];
//     B[0] = s2[i1];    B[1] = s2[i2];
//     paux[0] = p[i1];  paux[1] = p[i2];
//     daux[0] = d[i1];  daux[1] = d[i2];
//     
//     // printf("A: (%f, %f), B: (%f, %f), p: (%f, %f), d: (%f, %f)\n", A[0], A[1], B[0], B[1], p[0], p[1], d[0], d[1]);
//     assert(!(p[0] == d[0] && p[1] == d[1]));
//     assert(!(A[0] == B[0] && A[1] == B[1]));
//     return segments_intersect_2d(A, B, paux, daux);
// }
//
//
// int point_in_face(double **vertices_face, double *p, double *d) 
// {
//     double n[3], d1[3], d2[3];
//     d1[0] = vertices_face[1][0] - vertices_face[0][0];
//     d1[1] = vertices_face[1][1] - vertices_face[0][1];
//     d1[2] = vertices_face[1][2] - vertices_face[0][2];
//     d2[0] = vertices_face[2][0] - vertices_face[0][0];
//     d2[1] = vertices_face[2][1] - vertices_face[0][1];
//     d2[2] = vertices_face[2][2] - vertices_face[0][2];
//     cross_3d(d1, d2, n);
//
//     int drop_coord = 2;
//     if (fabs(n[0]) > fabs(n[1]) && fabs(n[0]) > fabs(n[2])) drop_coord = 0;
//     else if (fabs(n[1]) > fabs(n[0]) && fabs(n[1]) > fabs(n[2])) drop_coord = 1;
//
//     double *aux[3];
//     aux[0] = p; 
//
//     aux[1] = vertices_face[0];
//     aux[2] = vertices_face[1];
//     if (orientation(aux, d, 3) == 0 &&
//         pd_intersects_face(vertices_face[0], vertices_face[1], p, d, drop_coord)) return 1;
//
//     aux[1] = vertices_face[1];
//     aux[2] = vertices_face[2];
//     if (orientation(aux, d, 3) == 0 &&
//         pd_intersects_face(vertices_face[1], vertices_face[2], p, d, drop_coord)) return 1;
//
//     aux[1] = vertices_face[0];
//     aux[2] = vertices_face[2];
//     if (orientation(aux, d, 3) == 0 &&
//         pd_intersects_face(vertices_face[0], vertices_face[2], p, d, drop_coord)) return 1;
//     
//     return 0;
// }


int point_in_triangle_2d(double *v1, double *v2, double *v3, double *p)
{
    int o1 = orient2d(v1, v2, p);
    int o2 = orient2d(v2, v3, p);
    int o3 = orient2d(v3, v1, p);
    
    // Find reference sign (non-zero)
    int signs[3] = {o1, o2, o3};
    int ref_sign = 0;
    for (int ii=0; ii<3; ii++) {
        if (signs[ii] != 0) {
            ref_sign = signs[ii];
            break;
        }
    }

    for (int ii=0; ii<4; ii++) {
        if (signs[ii] != 0 && signs[ii] != ref_sign) return 0;
    }
    return 1;
}


int point_in_face(double **vertices_face, double *p)
{
    double n[3], d1[3], d2[3];
    d1[0] = vertices_face[1][0] - vertices_face[0][0];
    d1[1] = vertices_face[1][1] - vertices_face[0][1];
    d1[2] = vertices_face[1][2] - vertices_face[0][2];
    d2[0] = vertices_face[2][0] - vertices_face[0][0];
    d2[1] = vertices_face[2][1] - vertices_face[0][1];
    d2[2] = vertices_face[2][2] - vertices_face[0][2];
    cross_3d(d1, d2, n);

    int drop_coord = 2;
    if (fabs(n[0]) > fabs(n[1]) && fabs(n[0]) > fabs(n[2])) drop_coord = 0;
    else if (fabs(n[1]) > fabs(n[0]) && fabs(n[1]) > fabs(n[2])) drop_coord = 1;

    if (orientation(vertices_face, p, 3) != 0) return 0;
    

    int i1, i2;
    if (drop_coord == 0) {
        i1 = 1;     i2 = 2;
    } else if (drop_coord == 1) {
        i1 = 2;     i2 = 0;
    } else {
        i1 = 0;     i2 = 1;
    }
    double v1[2] = {vertices_face[0][i1], vertices_face[0][i2]},
           v2[2] = {vertices_face[1][i1], vertices_face[1][i2]},
           v3[2] = {vertices_face[2][i1], vertices_face[2][i2]},
           paux[2] = {p[i1], p[i2]};

    return point_in_triangle_2d(v1, v2, v3, paux);
}


int point_in_tetra(s_setup *setup, double *x, s_ncell *ncell)
{
    assert(setup->dim == 3 && "next test not implemented, just in 3D.");

    static double **facet_vertices = NULL;
    if (!facet_vertices) facet_vertices = malloc_matrix(3, 3);

    for (int ii=0; ii<4; ii++) {
        double *opposite_vertex = setup->points[ncell->vertex_id[ii]];
        extract_vertices_face(setup, ncell, &ii, setup->dim-1, facet_vertices);

        int o1 = orientation(facet_vertices, opposite_vertex, 3);
        int o2 = orientation(facet_vertices, x, 3);
        if (o1 == 0) {
            if (point_in_face(facet_vertices, x)) return 1;
            else continue;
        }
        if (o2 != 0 && o1 * o2 < 0) return 0;
    }
    return 1;
}


s_ncell *bruteforce_find_ncell_containing(s_setup *setup, double *p)
{
    s_ncell *current = setup->head;
    while (current) {
        if (point_in_tetra(setup, p, current)) return current;
        current = current->next;
    }
    puts("DID NOT FIND CONTAINER NCELL!");
    exit(1);
}


s_ncell *in_ncell_walk_NEW(s_setup *setup, double *p)
{
    // return bruteforce_find_ncell_containing(setup, p);

    s_ncell *current = setup->head;
    assert(setup->N_ncells >= 1 && "N_ncells < 1");
    int randi = (rand() % setup->N_ncells);
    for (int ii=0; ii<randi; ii++) {  // Select random ncell to start
        current = current->next;
    }

    // Create array for facet_vertices in static memory, CANNOT BE MULTI-THREADED! FIXME
    static int prev_dim = 0;
    static double **facet_vertices = NULL;
    if (setup->dim != prev_dim) {
        if (facet_vertices) free_matrix(facet_vertices, prev_dim);
        prev_dim = setup->dim;
        facet_vertices = malloc_matrix(setup->dim, setup->dim);
    }

    int steps = 0;
    int max_steps = 1000;
    while (steps++ < max_steps) {
        int moved = 0;
        for (int ii=0; ii<setup->dim+1; ii++) {
            double *opposite_vertex = setup->points[current->vertex_id[ii]];
            
            s_ncell *next = current->opposite[ii];
            if (next) {
                extract_vertices_face(setup, current, &ii, setup->dim-1, facet_vertices);

                int o1 = orientation(facet_vertices, opposite_vertex, setup->dim);
                int o2 = orientation(facet_vertices, p, setup->dim);

                if (o1 == 0 && o2 == 0) {  // p in plane of flat tetrahedron
                    puts("P in plane of flat tetrahedron.");
                    exit(1);
                } else if (o1 == 0) { // Tetrahedron is flat
                    printf("\n\nDEBUG WALK: tetrahedron is flat. current = %p, next = %p\n\n", (void*)current, (void*)next);
                    exit(1);
                } else if (o2 == 0) {  // p in face's plane
                    printf("DEBUG WALK: p in plane of face. o1 = %d, o2 = %d. INTETRA = %d\n", o1, o2, point_in_tetra(setup, p, current));
                    if (point_in_tetra(setup, p, current)) return current;
                    else continue;
                } else if (o1 != o2) {
                    current = next;
                    moved = 1;
                    break;
                }
            }
        }
        if (!moved) break;
    }
    assert(steps < max_steps);
    assert(point_in_tetra(setup, p, current)); // DEBUG TODO FIXME, REMOVE, TESTING
    return current;
}


s_ncell *in_ncell_walk(s_setup *setup, double *p)  // Should make sure that p is inside the convull of all points (inside an n-cell)
{
    s_ncell *current = setup->head;
    assert(setup->N_ncells >= 1 && "N_ncells < 1");
    int randi = (rand() % setup->N_ncells);
    for (int ii=0; ii<randi; ii++) {  // Select random ncell to start
        current = current->next;
    }

    // Create array for facet_vertices in static memory, CANNOT BE MULTI-THREADED! FIXME
    static int prev_dim = 0;
    static double **facet_vertices = NULL;
    if (setup->dim != prev_dim) {
        if (facet_vertices) free_matrix(facet_vertices, prev_dim);
        prev_dim = setup->dim;
        facet_vertices = malloc_matrix(setup->dim, setup->dim);
    }

    STEP:
    for (int ii=0; ii<setup->dim+1; ii++) {
        double *opposite_vertex = setup->points[current->vertex_id[ii]];

        s_ncell *next = current->opposite[ii];
        if (next) {
            extract_vertices_face(setup, current, &ii, setup->dim-1, facet_vertices);

            int o1 = orientation(facet_vertices, opposite_vertex, setup->dim);
            int o2 = orientation(facet_vertices, p, setup->dim);
            assert(o1 != 0);
            if (o2 == 0) {
                if (point_in_tetra(setup, p, current)) return current;
                // if (point_in_tetra(setup, p, next)) return next;
            } else if (o1 != o2) {
                current = next;
                goto STEP;
            }
        }
    }
    
    return current;
}



int is_delaunay_3d(const s_setup *setup)
{
    static double **vertices_ncell = NULL;
    if (!vertices_ncell) vertices_ncell = malloc_matrix(4, 3);

    s_ncell *current = setup->head;
    while (current) {
        extract_vertices_ncell(setup, current, vertices_ncell);
        s_ncell *query = current->next;
        while (query) {
            for (int ii=0; ii<4; ii++) {
                if (!inarray(current->vertex_id, 4, query->vertex_id[ii])) {
                    int in_sph = in_sphere(vertices_ncell, setup->points[query->vertex_id[ii]], 3);
                    if (in_sph == 1) {
                        printf("CONFLICT: (%d, %d, %d, %d) and %d, in_sphere: %d\n", 
                                current->vertex_id[0], current->vertex_id[1], 
                                current->vertex_id[2], current->vertex_id[3], 
                                query->vertex_id[ii], in_sphere(vertices_ncell, 
                                setup->points[query->vertex_id[ii]], 3));
                        return 0;
                    }
                }
            }
            query = query->next;
        }
        current = current->next;
    }
    return 1;
}


void add_ncell_volume_3d(s_setup *setup, s_ncell *ncell)
{   // THIS IS JUST FOR DEBUGGING!!
    
    ncell->volume = fabs(1.0/6.0 * 
                         orient3d(setup->points[ncell->vertex_id[0]], 
                                  setup->points[ncell->vertex_id[1]],
                                  setup->points[ncell->vertex_id[2]],
                                  setup->points[ncell->vertex_id[3]]));
    // static double **vertices = NULL;
    // if (!vertices) vertices = malloc_matrix(4, 3);
    //
    // extract_vertices_ncell(setup, ncell, vertices);
    //
    // // ncell->volume = compute_volume_convhull_from_points(vertices, 4);
    // double aux_vol = compute_volume_convhull_from_points(vertices, 4);
    // printf("DEBUG VOLUME NCELL: %f, %f\n", ncell->volume, aux_vol);
    // if (ncell->volume < 0) {
    //     printf("DEBUG ADD_NCELL_VOLUME: vol = %f, vertices = (%d, %d, %d, %d)\n", 
    //             ncell->volume, ncell->vertex_id[0], ncell->vertex_id[1], 
    //             ncell->vertex_id[2], ncell->vertex_id[3]);
    // }
}


double compute_volume_complex(s_setup *setup)
{
    return compute_volume_convhull_from_points(setup->points, setup->N_points);
}




// ----------------------------------------------------------------------------------------------
// --------------------------------------- PLOTS ------------------------------------------------
// ----------------------------------------------------------------------------------------------


void plot_add_ncell(FILE *pipe, s_setup *setup, s_ncell *ncell, char *config)
{
    static double **face_vertices = NULL;
    if (!face_vertices) face_vertices = malloc_matrix(3, 3);
    for (int ii=0; ii<4; ii++) {
        extract_vertices_face(setup, ncell, &ii, 2, face_vertices);
        fprintf(pipe, "\"<echo \'");
        fprintf(pipe, "%f %f %f\\n", face_vertices[0][0], face_vertices[0][1], face_vertices[0][2]);
        fprintf(pipe, "%f %f %f\\n", face_vertices[1][0], face_vertices[1][1], face_vertices[1][2]);
        fprintf(pipe, "%f %f %f'\"", face_vertices[2][0], face_vertices[2][1], face_vertices[2][2]);
        fprintf(pipe, "%s, ", config);
    }
}


void plot_ncell_3d(s_setup *setup, s_ncell *ncell, char *f_name, double *ranges)
{
    FILE *pipe = popen("gnuplot -persistent 2>&1", "w");
    fprintf(pipe, "set terminal pngcairo enhanced font 'Arial,18' size 1080,1080 enhanced \n");
    fprintf(pipe, "set pm3d depthorder\n");
    fprintf(pipe, "set pm3d border lc 'black' lw 0.5\n");
    fprintf(pipe, "set view 100, 60, \n");
    fprintf(pipe, "set xyplane at 0\n");
    if (ranges) {
        fprintf(pipe, "set xrange [%f:%f]\n", ranges[0], ranges[1]);
        fprintf(pipe, "set yrange [%f:%f]\n", ranges[2], ranges[3]);
        fprintf(pipe, "set zrange [%f:%f]\n", ranges[4], ranges[5]);
    }
    fprintf(pipe, "set xlabel 'x'\n");
    fprintf(pipe, "set ylabel 'y'\n");
    fprintf(pipe, "set zlabel 'z'\n");
    fflush(pipe);

    fprintf(pipe, "set output '%s.png'\n", f_name);
    fprintf(pipe, "splot ");
    plot_add_ncell(pipe, setup, ncell, "w polygons fs transparent solid 0.2 fc rgb '#000090' notitle");

    fprintf(pipe, "\n");
    pclose(pipe);
}


void plot_dt_3d(s_setup *setup, char *f_name, double *ranges, int max_files)
{
    char colors[][20] = { "#000090", "#000fff", "#0090ff", "#0fffee", 
        "#90ff70", "#ffee00", "#ff7000", "#ee0000", "#7f0000" };
    char buff[1024];

    FILE *pipe = popen("gnuplot -persistent 2>&1", "w");
    fprintf(pipe, "set terminal pngcairo enhanced font 'Arial,18' size 1080,1080 enhanced \n");
    fprintf(pipe, "set pm3d depthorder\n");
    fprintf(pipe, "set pm3d border lc 'black' lw 0.5\n");
    fprintf(pipe, "set view 100, 60, \n");
    fprintf(pipe, "set xyplane at 0\n");
    if (ranges) {
        fprintf(pipe, "set xrange [%f:%f]\n", ranges[0], ranges[1]);
        fprintf(pipe, "set yrange [%f:%f]\n", ranges[2], ranges[3]);
        fprintf(pipe, "set zrange [%f:%f]\n", ranges[4], ranges[5]);
    }
    fprintf(pipe, "set xlabel 'x'\n");
    fprintf(pipe, "set ylabel 'y'\n");
    fprintf(pipe, "set zlabel 'z'\n");
    fflush(pipe);


    // PLOT ALL CELLS
    fprintf(pipe, "set output '%s_v1.png'\n", f_name);
    fprintf(pipe, "splot ");

    int it = 0;
    s_ncell *current = setup->head;
    while (current) {
        snprintf(buff, 1024, "w polygons fs transparent solid 0.2 fc rgb '%s' notitle", colors[it%8]);
        plot_add_ncell(pipe, setup, current, buff);
        current = current->next;
        it++;
    }
    fprintf(pipe, "\n");

    fprintf(pipe, "set output '%s_v2.png'\n", f_name);
    fprintf(pipe, "set view 100, 90, 1.5\n");
    fprintf(pipe, "replot\n");

    fprintf(pipe, "set output '%s_v3.png'\n", f_name);
    fprintf(pipe, "set view 100, 180, 1.5\n");
    fprintf(pipe, "replot\n");

    fprintf(pipe, "set output '%s_v4.png'\n", f_name);
    fprintf(pipe, "set view 100, 270, 1.5\n");
    fprintf(pipe, "replot\n");



    // A NEW PLOT FOR EACH CELL (UNTIL MAX_FILES)
    it = 0;
    current = setup->head;
    while (current) {
        if (max_files != 0 && it > max_files) break;

        fprintf(pipe, "set output '%s_%d.png'\n", f_name, it);
        fprintf(pipe, "splot ");
        
        snprintf(buff, 1024, "w polygons fs transparent solid 0.2 fc rgb '%s' notitle", colors[it%8]);
        plot_add_ncell(pipe, setup, current, buff);

        for (int jj=0; jj<setup->N_points; jj++) {
            fprintf(pipe, "\"<echo \'");
            fprintf(pipe, "%f %f %f\\n", setup->points[jj][0], setup->points[jj][1], setup->points[jj][2]);
            fprintf(pipe, "'\" pt 7 lc rgb 'black' notitle, ");
        }
        fprintf(pipe, "\n");

        it++;
        current = current->next;
    }



    pclose(pipe);
}

#endif
