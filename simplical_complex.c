
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "geometry.c"


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
    struct facet **facet;
    struct ncell *next;  // Linked list of cells
} s_ncell;


typedef struct facet {
    int *vertex_id;
    struct ncell *incident_ncell[2];
} s_facet;


s_ncell *malloc_ncell(const s_setup *setup)
{
    s_ncell *out = malloc(sizeof(s_ncell));
    out->vertex_id = malloc(sizeof(int) * (setup->dim + 1));
    out->facet = malloc(sizeof(s_facet*) * (setup->dim + 1));
    out->next = NULL;
    return out;
}


s_facet *malloc_facet(const s_setup *setup)
{
    s_facet *out = malloc(sizeof(s_facet));
    out->vertex_id = malloc(sizeof(int) * setup->dim);
    return out;
}


s_ncell *traverse_facet(const s_setup *setup, const s_ncell *ncell, const int facet_id)
{
    assert((facet_id <= setup->dim + 1) && "facet_id > dim+1"); 

    s_facet *facet = ncell->facet[facet_id];
    assert(facet && "facet is NULL");

    s_ncell *ncell_1 = facet->incident_ncell[0]; 
    s_ncell *ncell_2 = facet->incident_ncell[1];
    if (ncell_1 != ncell && ncell_2 == ncell) {
        return ncell_1;
    } 
    if (ncell_1 == ncell && ncell_2 != ncell) {
        return ncell_2;
    }
    printf("Could not traverse facet, no incident ncell corresponds to the query ncell.");
    exit(1);
}


int id_where_equal_int(const int *arr, int N, int entry) 
{
    for (int ii=0; ii<N; ii++) {
        if (arr[ii] == entry) return ii;
    }
    assert(1 == 0 && "Could not find id."); 
    exit(1);
}


int count_cycle_ridge(const s_setup *setup, s_ncell *ncell, int v_localid_main, int v_localid_2)  // local 0:dim
{  // Traverses around the ridge towards facet pointed by v_localid_main
    // double ridge[setup->dim - 1];
    // int kk = 0;
    // for (int ii=0; ii<setup->dim-1; ii++) {
    //     if (ii != v_localid_main && ii != v_localid_2) {
    //         ridge[kk] = ncell->vertex_id[ii];
    //         kk++;
    //     }
    // }
    int counter = 1;
    s_ncell *current_ncell = ncell;
    s_ncell *next_ncell = traverse_facet(setup, ncell, v_localid_main);
    while (next_ncell != ncell) {
        counter++;
        // Find new ridge (look for s')
        v_localid_main = id_where_equal_int(next_ncell->vertex_id, setup->dim+1, current_ncell->vertex_id[v_localid_2]);
        for (int ii=0; ii<setup->dim-1; ii++) {
            if (traverse_facet(setup, next_ncell, ii) == ncell) {
                v_localid_2 = ii;
                break;
            }
        }
        current_ncell = next_ncell;
        next_ncell = traverse_facet(setup, next_ncell, v_localid_main);
    }

    return counter;
}


void extract_vertices_ncell(const s_setup *setup, const s_ncell *ncell, double **out)
{
    for (int ii=0; ii<setup->dim+1; ii++) {
        for (int jj=0; jj<setup->dim; jj++) {
            out[ii][jj] = setup->points[ncell->vertex_id[ii]][jj];
        }
    }
}


void extract_vertices_facet(const s_setup *setup, const s_facet *facet, double **out)
{
    for (int ii=0; ii<setup->dim; ii++) {
        for (int jj=0; jj<setup->dim; jj++) {
            out[ii][jj] = setup->points[facet->vertex_id[ii]][jj];
        }
    }
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
        big_ncell->facet[ii] = malloc_facet(setup);
        big_ncell->facet[ii]->incident_ncell[0] = big_ncell;
        big_ncell->facet[ii]->incident_ncell[1] = NULL;
        int kk = 0;  //kk is necessary for correct indexing
        for (int jj=0; jj<setup->dim+1; jj++) {
            if (ii != jj) {
                big_ncell->facet[ii]->vertex_id[kk] = N_points + jj;
                kk++;
            }
        }
    }
    setup->head = big_ncell;
    setup->N_ncells = 1;
    
    return setup;
}


int inarray(const int *arr1, int N, int a)
{
    for (int ii=0; ii<N; ii++) {
        if (arr1[ii] == a) return 1;
    }
    return 0;
}


void find_unique_pair(const int *arr1, const int *arr2, int N, int *unique1, int *unique2)
{
    int xorAll = 0; // Compute XOR for all elements in both lists.
    for (int ii=0; ii<N; ii++) {
        xorAll ^= arr1[ii];
        xorAll ^= arr2[ii];
    }
    int diffBit = xorAll & -xorAll; // Find rightmost set bit

    *unique1 = 0; *unique2 = 0;
    for (int ii=0; ii<N; ii++) {  // Partition elements based on the set bit.
        if (arr1[ii] & diffBit) {
            *unique1 ^= arr1[ii];
        } else {
            *unique2 ^= arr1[ii];
        }
        if (arr2[ii] & diffBit) {
            *unique1 ^= arr2[ii];
        } else {
            *unique2 ^= arr2[ii];
        }
    }

    if (!inarray(arr1, N, *unique1)) {  // Correct output, flip 1 and 2 if necessary
        int aux = *unique1;
        *unique1 = *unique2;
        *unique2 = aux;
    }
}


int is_locally_delaunay(const s_setup *setup, const s_facet *facet)
{
    s_ncell *ncell1 = facet->incident_ncell[0];
    s_ncell *ncell2 = facet->incident_ncell[1];
    
    int v1, v2;
    find_unique_pair(ncell1->vertex_id, ncell2->vertex_id, setup->dim + 1, &v1, &v2);
    // printf("v1: %d, v2: %d\n", v1, v2);
    
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
    // printf("oriented: %d, insphere: %d\n", orientation(coords1, coords1[2], 2), out);
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
    for (int ii=0; ii<rand()%setup->N_ncells; ii++) {  // Select random ncell to start
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

        s_facet *test_facet = current->facet[ii];
        s_ncell *next = traverse_facet(setup, current, ii);
        if (next) {
            extract_vertices_facet(setup, test_facet, facet_vertices);

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
