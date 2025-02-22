

#include <stdio.h>
#include <time.h>
// #include "algebra.c"
#include "geometry.c"
#include "simplical_complex.c"

int main() {
    srand(time(NULL));

    // FIND UNIQUE PAIR, seems OK -----------------------------------
    int arr1[] = {1, 3, 7, 9, 5};
    int arr2[] = {7, 5, 2, 1, 9};
    int u1, u2;
    find_unique_pair(arr1, arr2, 5, &u1, &u2);
    printf("u1: %d, u2: %d\n\n", u1, u2);


    // INSPHERE, seems OK -------------------------------------------
    double **p = malloc_matrix(4, 3);
    p[0][0] = 0;    p[0][1] = 1;    p[0][2] = 0;
    p[1][0] = -1;   p[1][1] = 0;    p[1][2] = 0;
    p[2][0] = 0;    p[2][1] = 0;    p[2][2] = 1;
    p[3][0] = 0;    p[3][1] = -1;   p[3][2] = 0;

    double q[3] = {0, 0, 0.2};

    int result = orientation(p, p[3], 3);
    printf("orientation: %d\n", result);

    result = insphere(p, q, 3);
    printf("insphere: %d\n\n", result);


    // ISDELAUNAY, seems OK -----------------------------------------
    double **points = malloc_matrix(5, 2);
    points[0][0] = 1;    points[0][1] = 0;
    points[1][0] = 0;    points[1][1] = 1;
    points[2][0] = -3;   points[2][1] = 0;
    points[3][0] = 0;    points[3][1] = -1;
    points[4][0] = 3;    points[4][1] = -3;
    s_setup setup = {.dim = 2, .N_points = 5, .points = points, .N_ncells = 4};

    s_ncell *ncell1 = malloc_ncell(&setup);
    s_ncell *ncell2 = malloc_ncell(&setup);
    s_ncell *ncell3 = malloc_ncell(&setup);
    s_ncell *ncell4 = malloc_ncell(&setup);
    s_facet *facet1 = malloc_facet(&setup);
    s_facet *facet2 = malloc_facet(&setup);
    s_facet *facet3 = malloc_facet(&setup);
    s_facet *facet4 = malloc_facet(&setup);
    s_facet *facetO_1 = malloc_facet(&setup);
    s_facet *facetO_2 = malloc_facet(&setup);
    s_facet *facetO_3 = malloc_facet(&setup);
    s_facet *facetO_4 = malloc_facet(&setup);
    
    setup.head = ncell2;
    ncell2->next = ncell1;
    ncell1->next = ncell3;
    ncell3->next = ncell4;

    ncell1->vertex_id[0] = 0;
    ncell1->vertex_id[1] = 1;
    ncell1->vertex_id[2] = 3;
    ncell1->facet[0] = facet1;
    ncell1->facet[1] = facet2;
    ncell1->facet[2] = facetO_1;

    ncell2->vertex_id[0] = 1;
    ncell2->vertex_id[1] = 2;
    ncell2->vertex_id[2] = 3;
    ncell2->facet[0] = facet4;
    ncell2->facet[1] = facet1;
    ncell2->facet[2] = facetO_2;

    ncell3->vertex_id[0] = 0;
    ncell3->vertex_id[1] = 3;
    ncell3->vertex_id[2] = 4;
    ncell3->facet[0] = facet3;
    ncell3->facet[1] = facetO_3;
    ncell3->facet[2] = facet2;

    ncell4->vertex_id[0] = 2;
    ncell4->vertex_id[1] = 4;
    ncell4->vertex_id[2] = 3;
    ncell4->facet[0] = facet3;
    ncell4->facet[1] = facet4;
    ncell4->facet[2] = facetO_4;

    facet1->vertex_id[0] = 1;
    facet1->vertex_id[1] = 3;
    facet1->incident_ncell[0] = ncell1;
    facet1->incident_ncell[1] = ncell2;

    facet2->vertex_id[0] = 0;
    facet2->vertex_id[1] = 3;
    facet2->incident_ncell[0] = ncell1;
    facet2->incident_ncell[1] = ncell3;

    facet3->vertex_id[0] = 3;
    facet3->vertex_id[1] = 4;
    facet3->incident_ncell[0] = ncell3;
    facet3->incident_ncell[1] = ncell4;

    facet4->vertex_id[0] = 2;
    facet4->vertex_id[1] = 4;
    facet4->incident_ncell[0] = ncell2;
    facet4->incident_ncell[1] = ncell4;

    facetO_1->vertex_id[0] = 0;
    facetO_1->vertex_id[1] = 1;
    facetO_1->incident_ncell[0] = ncell1;
    facetO_1->incident_ncell[1] = NULL;

    facetO_2->vertex_id[0] = 1;
    facetO_2->vertex_id[1] = 2;
    facetO_2->incident_ncell[0] = ncell2;
    facetO_2->incident_ncell[1] = NULL;

    facetO_3->vertex_id[0] = 0;
    facetO_3->vertex_id[1] = 4;
    facetO_3->incident_ncell[0] = ncell3;
    facetO_3->incident_ncell[1] = NULL;

    facetO_4->vertex_id[0] = 2;
    facetO_4->vertex_id[1] = 4;
    facetO_4->incident_ncell[0] = ncell4;
    facetO_4->incident_ncell[1] = NULL;

    result = is_locally_delaunay(&setup, facet1);
    printf("is_locally_delaunay(f1): %d\n", result);
    result = is_locally_delaunay(&setup, facet2);
    printf("is_locally_delaunay(f2): %d\n\n", result);


    // WALKING, seems OK --------------------------------------------
    double query[2] = {2, -2};
    s_ncell *selected = in_ncell_walk(&setup, query);
    printf("inside: %d, %d, %d\n\n", selected->vertex_id[0], selected->vertex_id[1], selected->vertex_id[2]);

    query[0] = -0.2;   query[1] = 0;
    selected = in_ncell_walk(&setup, query);
    printf("inside: %d, %d, %d\n\n", selected->vertex_id[0], selected->vertex_id[1], selected->vertex_id[2]);


    // SETUP COMPLEX
    s_setup *setup_test = initialize_setup(2, points, 5);


    // CYCLIC ORDERING AROUND RIDGE
    int count = count_cycle_ridge(&setup, ncell1, 0, 1);
    printf("count of ncells around ridge (0,0) : %d\n\n", count);

}
