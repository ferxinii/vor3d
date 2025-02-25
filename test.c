

#include <stdio.h>
#include <time.h>
// #include "algebra.c"
#include "geometry.c"
#include "simplical_complex.c"
#include "dt_3d_incremental.c"

int main(void) {
    srand(time(NULL));

    // FIND UNIQUE PAIR, seems OK  ----------------------------------
    int arr1[] = {1, 3, 7, 9, 5};
    int arr2[] = {7, 5, 2, 1, 9};
    int u1, u2;
    find_unique_pair(arr1, arr2, 5, &u1, &u2);
    printf("u1: %d, u2: %d\n\n", u1, u2);


    // INSPHERE, seems OK  ------------------------------------------
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


    // ISDELAUNAY, seems OK  ----------------------------------------
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

    ncell1->vertex_id[0] = 0;       ncell1->vertex_id[1] = 1;       ncell1->vertex_id[2] = 3;
    ncell1->opposite[0] = ncell2;   ncell1->opposite[1] = ncell3;   ncell1->opposite[2] = NULL;

    ncell2->vertex_id[0] = 1;       ncell2->vertex_id[1] = 2;       ncell2->vertex_id[2] = 3;
    ncell2->opposite[0] = ncell4;   ncell2->opposite[1] = ncell1;   ncell2->opposite[2] = NULL;

    ncell3->vertex_id[0] = 0;       ncell3->vertex_id[1] = 3;       ncell3->vertex_id[2] = 4;
    ncell3->opposite[0] = ncell4;   ncell3->opposite[1] = NULL;     ncell3->opposite[2] = ncell1;

    ncell4->vertex_id[0] = 2;       ncell4->vertex_id[1] = 4;       ncell4->vertex_id[2] = 3;
    ncell4->opposite[0] = ncell3;   ncell4->opposite[1] = ncell2;   ncell4->opposite[2] = NULL;

    result = are_locally_delaunay(&setup, ncell1, 0);
    printf("is_locally_delaunay(f1): %d\n", result);
    result = are_locally_delaunay(&setup, ncell1, 1);
    printf("is_locally_delaunay(f2): %d\n", result);
    result = are_locally_delaunay(&setup, ncell3, 0);
    printf("is_locally_delaunay(f3): %d\n", result);
    result = are_locally_delaunay(&setup, ncell4, 1);
    printf("is_locally_delaunay(f4): %d\n\n", result);


    // WALKING, seems OK  -------------------------------------------
    double query[2] = {2, -2};
    s_ncell *selected = in_ncell_walk(&setup, query);
    printf("inside: %d, %d, %d\n\n", selected->vertex_id[0], selected->vertex_id[1], selected->vertex_id[2]);

    query[0] = -0.2;   query[1] = 0;
    selected = in_ncell_walk(&setup, query);
    printf("inside: %d, %d, %d\n\n", selected->vertex_id[0], selected->vertex_id[1], selected->vertex_id[2]);


    // SETUP COMPLEX  -----------------------------------------------
    s_setup *setup_test = initialize_setup(points, 5, 2);
    (void) setup_test;


    // CYCLIC ORDERING AROUND RIDGE  --------------------------------
    int count = count_cycle_ridge(&setup, ncell1, 0, 1);
    printf("count of ncells around ridge (0,0) : %d\n\n", count);
    

    // MARKING NCELLS SHARING FACE, seems OK  -----------------------
    int v_localid[2] = {0, 1};
    initialize_ncells_mark(&setup);
    print_marked(&setup);
    mark_ncells_incident_face(&setup, ncell1, v_localid, 0);
    count = count_marked(&setup);
    printf("count with general function : %d\n\n", count);
    print_marked(&setup);

    int v_localid_2[1] = {0};
    initialize_ncells_mark(&setup);
    print_marked(&setup);
    mark_ncells_incident_face(&setup, ncell1, v_localid_2, 1);
    count = count_marked(&setup);
    printf("count with general function : %d\n\n", count);
    print_marked(&setup);

    initialize_ncells_mark(&setup);
    print_marked(&setup);
    mark_ncells_incident_face(&setup, ncell2, v_localid_2, 1);
    count = count_marked(&setup);
    printf("count with general function : %d\n\n", count);
    print_marked(&setup);


    // INCREMENTAL ALGORITHM
    double **p_dt = malloc_matrix(2, 3);
    p_dt[0][0] = -1;    p[0][1] = 0;    p[0][2] = 0;
    p_dt[1][0] = 1;     p[1][1] = 0;    p[1][2] = 0;

    s_stack *stack = stack_create();
    s_setup *setup_dt = initialize_setup(p_dt, 2, 3);
    print_ncells(setup_dt);

    
    // FLIP14, seems OK
    flip14(setup_dt, setup_dt->head, 0, stack);
    print_ncells(setup_dt);

    stack_print(stack);


    // INSERT_ONE_POINT
    s_stack *stack2 = stack_create();
    s_setup *setup2 = initialize_setup(p_dt, 2, 3);
    printf("Before insertion, N_ncells = %d, dim = %d\n", setup2->N_ncells, setup2->dim);
    puts("\nInserting 0:");
    insert_one_point(setup2, 0, stack2);
    printf("After insertion, N_ncells = %d\n", setup2->N_ncells);
    print_ncells(setup_dt);

    puts("\nInserting 1:");
    insert_one_point(setup2, 1, stack2);
    printf("After insertion, N_ncells = %d\n", setup2->N_ncells);
    print_ncells(setup_dt);


    // SEGMENT_CROSSES_TRIANGLE_3D
    double **triangle = malloc_matrix(3, 3);
    triangle[0][0] = 1;     triangle[0][1] = 0;     triangle[0][2] = 0;
    triangle[1][0] = 0;     triangle[1][1] = 1;     triangle[1][2] = 0;
    triangle[2][0] = 0;     triangle[2][1] = 0;     triangle[2][2] = 1;
    double a[3] = {0.1, 0.1, 3};
    double b[3] = {0.1, 0.1, -3};
    
    printf("crosses_triangle: %d\n", segment_crosses_triangle_3d(triangle, a, b));


    // FLIP23
    double **p2 = malloc_matrix(5, 3);
    p2[0][0] = 1;    p2[0][1] = 0;    p2[0][2] = 0;
    p2[1][0] = -1;   p2[1][1] = 1;    p2[1][2] = 0;
    p2[2][0] = 0;    p2[2][1] = -1;   p2[2][2] = 0;
    p2[3][0] = 0;    p2[3][1] = 0;    p2[3][2] = 1;
    p2[4][0] = 0;    p2[4][1] = 0;    p2[4][2] = -1;

    s_setup s2 = {.dim = 3, .N_ncells = 2, .points = p2};
    s_ncell *nc1 = malloc_ncell(&s2);
    s_ncell *nc2 = malloc_ncell(&s2);
    s2.head = nc1;
    nc1->next = nc2;

    nc1->vertex_id[0] = 0;      nc1->vertex_id[1] = 1;
    nc1->vertex_id[2] = 2;      nc1->vertex_id[3] = 3;
    nc1->opposite[0] = NULL;    nc1->opposite[1] = NULL;
    nc1->opposite[2] = NULL;    nc1->opposite[3] = nc2;

    nc2->vertex_id[0] = 0;      nc2->vertex_id[1] = 1;
    nc2->vertex_id[2] = 2;      nc2->vertex_id[3] = 4;
    nc2->opposite[0] = NULL;    nc2->opposite[1] = NULL;
    nc2->opposite[2] = NULL;    nc2->opposite[3] = nc1;

    print_ncells(&s2);

    flip23(&s2, stack, nc1, 3, 3);
    puts("\nFLIP23:");
    print_ncells(&s2);

    flip32(&s2, stack, nc1, 0, 1, 2); 
    puts("\nFLIP32:");
    print_ncells(&s2);


    // CONSTRUCT_DT_3D
    s_setup *dt_setup = construct_dt_3d(p2, 5);
    print_ncells(dt_setup);

    
}
