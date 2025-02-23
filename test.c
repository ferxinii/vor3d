

#include <stdio.h>
#include <time.h>
// #include "algebra.c"
#include "geometry.c"
#include "simplical_complex.c"

int main(void) {
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

    setup.head = ncell1;
    ncell1->next = ncell2;
    ncell2->next = ncell3;
    ncell3->next = ncell4;

    ncell1->vertex_id[0] = 0;
    ncell1->vertex_id[1] = 1;
    ncell1->vertex_id[2] = 3;
    ncell1->opposite[0] = ncell2;
    ncell1->opposite[1] = ncell3;
    ncell1->opposite[2] = NULL;

    ncell2->vertex_id[0] = 1;
    ncell2->vertex_id[1] = 2;
    ncell2->vertex_id[2] = 3;
    ncell2->opposite[0] = ncell4;
    ncell2->opposite[1] = ncell1;
    ncell2->opposite[2] = NULL;

    ncell3->vertex_id[0] = 0;
    ncell3->vertex_id[1] = 3;
    ncell3->vertex_id[2] = 4;
    ncell3->opposite[0] = ncell4;
    ncell3->opposite[1] = NULL;
    ncell3->opposite[2] = ncell1;

    ncell4->vertex_id[0] = 2;
    ncell4->vertex_id[1] = 4;
    ncell4->vertex_id[2] = 3;
    ncell4->opposite[0] = ncell3;
    ncell4->opposite[1] = ncell2;
    ncell4->opposite[2] = NULL;

    result = are_locally_delaunay(&setup, ncell1, ncell2);
    printf("is_locally_delaunay(f1): %d\n", result);
    result = are_locally_delaunay(&setup, ncell1, ncell3);
    printf("is_locally_delaunay(f2): %d\n", result);
    result = are_locally_delaunay(&setup, ncell3, ncell4);
    printf("is_locally_delaunay(f3): %d\n", result);
    result = are_locally_delaunay(&setup, ncell2, ncell4);
    printf("is_locally_delaunay(f4): %d\n\n", result);


    // WALKING, seems OK --------------------------------------------
    double query[2] = {2, -2};
    s_ncell *selected = in_ncell_walk(&setup, query);
    printf("inside: %d, %d, %d\n\n", selected->vertex_id[0], selected->vertex_id[1], selected->vertex_id[2]);

    query[0] = -0.2;   query[1] = 0;
    selected = in_ncell_walk(&setup, query);
    printf("inside: %d, %d, %d\n\n", selected->vertex_id[0], selected->vertex_id[1], selected->vertex_id[2]);


    // SETUP COMPLEX
    s_setup *setup_test = initialize_setup(2, points, 5);
    (void) setup_test;


    // CYCLIC ORDERING AROUND RIDGE
    int count = count_cycle_ridge(&setup, ncell1, 0, 1);
    printf("count of ncells around ridge (0,0) : %d\n\n", count);

}
