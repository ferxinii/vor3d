
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
    double *CM;  // Center of Mass
    int N_ncells;
    struct ncell *head;  // Linked list of ncells
} s_setup;


typedef struct ncell {
    int *vertex_id;
    struct ncell **opposite;
    struct ncell *next;  // Linked list of cells
    int mark;  // Used in flood-fill algorithms, i.e. to find the star of an (n-3)-cell.
} s_ncell;


s_ncell *malloc_ncell(const s_setup *setup)
{
    s_ncell *out = malloc(sizeof(s_ncell));
    out->vertex_id = malloc(sizeof(int) * (setup->dim + 1));
    out->opposite = malloc(sizeof(s_ncell*) * (setup->dim + 1));
    out->next = NULL;
    out->mark = 0;
    return out;
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


void extract_vertices_ncell(const s_setup *setup, const s_ncell *ncell, double **out)
{
    for (int ii=0; ii<setup->dim+1; ii++) {
        for (int jj=0; jj<setup->dim; jj++) {
            out[ii][jj] = setup->points[ncell->vertex_id[ii]][jj];
        }
    }
}


// void extract_ids_facet(const s_setup *setup, const s_ncell *ncell1, int v1_localid, int *out)
// {  
//     int kk=0;
//     for (int ii=0; ii<setup->dim+1; ii++) {
//         if (ii != v1_localid) {
//             out[kk] = ncell1->vertex_id[ii];
//             kk++;
//         }
//     }
// }


// void extract_vertices_facet(const s_setup *setup, const s_ncell *ncell1, int v1_localid, double **out)
// {
//     int facet_vertex_id[setup->dim];
//     extract_ids_facet(setup, ncell1, v1_localid, facet_vertex_id);
//     for (int ii=0; ii<setup->dim; ii++) {
//         for (int jj=0; jj<setup->dim; jj++) {
//             out[ii][jj] = setup->points[facet_vertex_id[ii]][jj];
//         }
//     }
// }


// void extract_ids_ridge(const s_setup *setup, const s_ncell *ncell1, int v1_localid, int v2_localid, int *out)
// {
//     // TODO   
// }


// void extract_ids_n3cell(const s_setup *setup, const s_ncell *ncell, int v1_localid, int v2_localid, int v3_localid, int *out)
// {   // TODO TEST !!
//     assert(v2_localid != v1_localid && v2_localid != v3_localid && v1_localid != v3_localid && "Must be different ids.");
//     int kk = 0;
//     for (int ii=0; ii<setup->dim+1; ii++) {
//         if (ii != v1_localid && ii != v2_localid && ii != v3_localid) {
//             out[kk] = ncell->vertex_id[ii];
//             kk++;
//         }
//     }
// }


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


void find_face_localid_of_adjacent_ncell(const s_setup *setup, const s_ncell *ncell, const int *v_localid,
                                         int dim_face, int id_adjacent, int *out_v_localid)
{
    s_ncell *adjacent = ncell->opposite[id_adjacent];

    int N_v = setup->dim - dim_face;
    int vertex_id[dim_face+1];
    extract_ids_face(setup, ncell, v_localid, dim_face, vertex_id);

    int kk = 0;
    for (int ii=0; ii<setup->dim+1; ii++) {
        if (!inarray(vertex_id, N_v, adjacent->vertex_id[ii])) {
            out_v_localid[kk] = ii;
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


void adjacent_localid_of_face(const s_setup *setup, const s_ncell *ncell, const int *v_localid,
                              int dim_face, const s_ncell *adjacent, int *out_v_localid)
{
    int vertex_id_face[setup->dim - dim_face];
    extract_ids_face(setup, ncell, v_localid, dim_face, out_v_localid);
    
    int kk = 0;
    for (int ii=0; ii<setup->dim+1; ii++) {
        if (!inarray(vertex_id_face, setup->dim - dim_face, adjacent->vertex_id[ii])) {
            out_v_localid[kk] = ii;
            kk++;
        }
    }
}


void mark_ncells_incident_face_STEP(const s_setup *setup, const s_ncell *ncell, const int *v_localid, int dim_face)
{
    for (int ii=0; ii<setup->dim+1; ii++) {
        if (inarray(v_localid, setup->dim - dim_face, ii)) {
            s_ncell *adjacent_ncell = ncell->opposite[ii];  // This ncell should already share the face!
            if (adjacent_ncell->mark == 0) {
                adjacent_ncell->mark = 1;
                
                int new_v_localid[setup->dim - dim_face];
                adjacent_localid_of_face(setup, ncell, v_localid, dim_face, adjacent_ncell, new_v_localid);

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


void find_center_mass(double **in, int N_points, int dim, double *out)
{
    for (int ii=0; ii<dim; ii++) {
        out[ii] = in[0][ii];
    }

    for (int ii=1; ii<N_points; ii++) {
        for (int jj=0; jj<dim; jj++) {
            out[jj] += in[ii][jj];
        }
    }

    for (int ii=0; ii<dim; ii++) {
        out[ii] /= N_points;
    }
}


s_setup *initialize_setup(int dim, double **points, int N_points)
{
    // setup->points is EXTENDED FOR THE EXTRA NODES OF BIG_NCELL
    double **setup_points = malloc_matrix(N_points + dim + 1, dim);
    for (int ii=0; ii<N_points; ii++) {
        for (int jj=0; jj<dim; jj++) {
            setup_points[ii][jj] = points[ii][jj];
        }
    }
    
    double *CM = malloc(sizeof(double) * dim);
    find_center_mass(points, N_points, dim, CM);
    double maxd = max_distance(points, N_points, dim, CM);

    // Build the vertices of a regular simplex in R^(dim+1) with circumsphere of radius s centered at origin
    double s = 1.5 * maxd * dim * sqrt((dim+1.0)/dim);  // Scale so that inradius = 1.5 * maxd, original norm = sqrt(dim/(dim+1))
    double **V = malloc_matrix(dim+1, dim+1); 
    for (int ii = 0; ii < dim+1; ii++) {
        for (int jj = 0; jj < dim+1; jj++) {
            V[ii][jj] = ((ii == jj) ? 1.0 : 0.0) - 1.0/(dim+1);
            V[ii][jj] *= s;
        }
    }
    // Project in Rn by removing the last coordinate, and add CM to center around points
    for (int ii=0; ii<dim+1; ii++) {
        for (int jj=0; jj<dim; jj++) {
            setup_points[N_points + ii][jj] = CM[jj] + V[ii][jj];
        }
    }
    free_matrix(V, dim+1);
    
    for (int ii=0; ii<N_points; ii++) {  // DEBUG 
        assert(insphere(&setup_points[N_points], setup_points[ii], dim) == 1 && "big_ncell does not enclose all points.");
    }

    s_setup *setup = malloc(sizeof(s_setup));
    setup->dim = dim;
    setup->N_points = N_points + dim + 1;
    setup->points = setup_points;
    setup->CM = CM;

    s_ncell *big_ncell = malloc_ncell(setup);
    for (int ii=0; ii<setup->dim+1; ii++) {
        big_ncell->vertex_id[ii] = N_points + ii;
        big_ncell->opposite[ii] = NULL;
    }
    setup->head = big_ncell;
    setup->N_ncells = 1;
    
    return setup;
}


int are_locally_delaunay(const s_setup *setup, const s_ncell *ncell1, const s_ncell *ncell2)
{   // TODO Check if incident?
    int v1, v2;
    find_unique_pair(ncell1->vertex_id, ncell2->vertex_id, setup->dim + 1, &v1, &v2);
    
    // Create array for coords1 in static memory, CANNOT BE MULTI-THREADED! FIXME
    static int prev_dim = 0;
    static double **coords1 = NULL;
    if (setup->dim != prev_dim) {
        prev_dim = setup->dim;
        if (coords1) free_matrix(coords1, setup->dim + 1);
        coords1 = malloc_matrix(setup->dim + 1, setup->dim);
    }

    extract_vertices_ncell(setup, ncell1, coords1);
    int out = insphere(coords1, setup->points[v2], setup->dim);
    if (out == -1) {
        return 1;
    } else {
        return 0;
    }
}


s_ncell *in_ncell_walk(s_setup *setup, double *p)  // Should make sure that p is inside the convull of all points (inside an n-cell)
{
    s_ncell *current = setup->head;
    assert(setup->N_ncells >= 1 && "N_ncells < 1");
    for (int ii=0; ii<(rand() % setup->N_ncells); ii++) {  // Select random ncell to start
        current = current->next;
    }

    // Create array for facet_vertices in static memory, CANNOT BE MULTI-THREADED! FIXME
    static int prev_dim = 0;
    static double **facet_vertices = NULL;
    if (setup->dim != prev_dim) {
        prev_dim = setup->dim;
        if (facet_vertices) free_matrix(facet_vertices, setup->dim);
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
            assert(o1 != 0 && o2 != 0 && "Could not walk. p inside a facet?");
            if (o1 != o2) {
                current = next;
                goto STEP;
            }
        }
    }
    
    return current;
}
