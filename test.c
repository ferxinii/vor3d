// TODO : 
// [ ] IS DELAUNAY FUNCTION ONLY ON REAL NCELLS??
// [ ] Check if points are in general position
// [ ] Function to check if a point is inside a simplex, using barycentric coordinates?
// [ ] Plotting 3D, to check and produce visual output?
// [ ] Better robustness in geometric predicates? Tolerance?
// [ ] Generation of points
// [ ] Extracting the Voronoi Diagram

#include <stdio.h>
#include <time.h>
#include "geometry.c"
#include "simplical_complex.c"
#include "dt_3d_incremental.c"


void random_points_3d(int N, double **out) 
{
    for (int ii = 0; ii < N; ii++) {
        out[ii][0] = (double) rand() / (double) RAND_MAX;
        out[ii][1] = (double) rand() / (double) RAND_MAX;
        out[ii][2] = (double) rand() / (double) RAND_MAX;
    }
}


int main(void) {
    srand(time(NULL));

    // FIND UNIQUE PAIR, seems OK  ----------------------------------
    // int arr1[] = {1, 3, 7, 9, 5};
    // int arr2[] = {7, 5, 2, 1, 9};
    // int u1, u2;
    // find_unique_pair(arr1, arr2, 5, &u1, &u2);
    // printf("u1: %d, u2: %d\n\n", u1, u2);


    // in_sphere, seems OK  ------------------------------------------
    double **p = malloc_matrix(4, 3);
    p[0][0] = 0;    p[0][1] = 1;    p[0][2] = 0;
    p[1][0] = -1;   p[1][1] = 0;    p[1][2] = 0;
    p[2][0] = 0;    p[2][1] = 0;    p[2][2] = 1;
    p[3][0] = 0;    p[3][1] = -1;   p[3][2] = 0;

    double q[3] = {0, 0, 0.2};

    int result = orientation(p, p[3], 3);
    printf("orientation: %d\n", result);

    result = in_sphere(p, q, 3);
    printf("in_sphere: %d\n\n", result);


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
    p_dt[0][0] = -1;    p_dt[0][1] = 0;      p_dt[0][2] = 0;
    p_dt[1][0] = 1;     p_dt[1][1] = 0.5;    p_dt[1][2] = 0;

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
    print_matrix(setup2->points, 6, 3);

    // plot_ncell_3d(setup2, setup2->head, "flip14/before");

    printf("Before insertion, N_ncells = %d, dim = %d\n", setup2->N_ncells, setup2->dim);
    puts("\nInserting 0:");
    insert_one_point(setup2, 0, stack2);
    printf("After insertion, N_ncells = %d\n", setup2->N_ncells);
    print_ncells(setup2);

    // plot_ncell_3d(setup2, setup2->head, "flip14/after1_1");
    // plot_ncell_3d(setup2, setup2->head->next, "flip14/after1_2");
    // plot_ncell_3d(setup2, setup2->head->next->next, "flip14/after1_3");
    // plot_ncell_3d(setup2, setup2->head->next->next->next, "flip14/after1_4");

    puts("\nInserting 1:");
    insert_one_point(setup2, 1, stack2);
    printf("After insertion, N_ncells = %d\n", setup2->N_ncells);
    print_ncells(setup2);

    // plot_ncell_3d(setup2, setup2->head, "flip14/after2_1");
    // plot_ncell_3d(setup2, setup2->head->next, "flip14/after2_2");
    // plot_ncell_3d(setup2, setup2->head->next->next, "flip14/after2_3");
    // plot_ncell_3d(setup2, setup2->head->next->next->next, "flip14/after2_4");
    // plot_ncell_3d(setup2, setup2->head->next->next->next->next, "flip14/after2_5");
    // plot_ncell_3d(setup2, setup2->head->next->next->next->next->next, "flip14/after2_6");
    // plot_ncell_3d(setup2, setup2->head->next->next->next->next->next->next, "flip14/after2_7");


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

    s_setup s2 = {.dim = 3, .N_ncells = 2, .N_points = 5, .points = p2};
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
    // plot_ncell_3d(&s2, s2.head, "flip23/before1");
    // plot_ncell_3d(&s2, s2.head->next, "flip23/before2");

    flip23(&s2, stack, nc1, 3, 3);
    puts("\nFLIP23:");
    print_ncells(&s2);
    // plot_ncell_3d(&s2, s2.head, "flip23/after1");
    // plot_ncell_3d(&s2, s2.head->next, "flip23/after2");
    // plot_ncell_3d(&s2, s2.head->next->next, "flip23/after3");

    flip32(&s2, stack, nc1, 0, 1, 2); 
    puts("\nFLIP32:");
    print_ncells(&s2);

    
    // FLIP44, seems ok!
    double **p3 = malloc_matrix(6, 3);
    p3[0][0] = 1;    p3[0][1] = 0;    p3[0][2] = 0;
    p3[1][0] = 0;    p3[1][1] = 1;    p3[1][2] = 0;
    p3[2][0] = -1;   p3[2][1] = 0;    p3[2][2] = 0;
    p3[3][0] = 0;    p3[3][1] = -1;   p3[3][2] = 0;
    p3[4][0] = 0;    p3[4][1] = 0;    p3[4][2] = 1;
    p3[5][0] = 0;    p3[5][1] = 0;    p3[5][2] = -1;

    s_setup s3 = {.dim = 3, .N_ncells = 4, .N_points = 6, .points = p3};
    nc1 = malloc_ncell(&s3);
    nc2 = malloc_ncell(&s3);
    s_ncell *nc3 = malloc_ncell(&s3);
    s_ncell *nc4 = malloc_ncell(&s3);

    s3.head = nc1;
    nc1->next = nc2;

    nc2->next = nc3;
    nc2->prev = nc1;

    nc3->next = nc4;
    nc3->prev = nc2;

    nc4->prev = nc3;

    nc1->vertex_id[0] = 0;      nc1->vertex_id[1] = 1;
    nc1->vertex_id[2] = 4;      nc1->vertex_id[3] = 3;
    nc1->opposite[0] = nc2;     nc1->opposite[1] = NULL;
    nc1->opposite[2] = nc4;     nc1->opposite[3] = NULL;

    nc2->vertex_id[0] = 1;      nc2->vertex_id[1] = 2;
    nc2->vertex_id[2] = 3;      nc2->vertex_id[3] = 4;
    nc2->opposite[0] = NULL;    nc2->opposite[1] = nc1;
    nc2->opposite[2] = NULL;    nc2->opposite[3] = nc3;

    nc3->vertex_id[0] = 1;      nc3->vertex_id[1] = 2;
    nc3->vertex_id[2] = 3;      nc3->vertex_id[3] = 5;
    nc3->opposite[0] = NULL;    nc3->opposite[1] = nc4;
    nc3->opposite[2] = NULL;    nc3->opposite[3] = nc2;

    nc4->vertex_id[0] = 1;      nc4->vertex_id[1] = 3;
    nc4->vertex_id[2] = 5;      nc4->vertex_id[3] = 0;
    nc4->opposite[0] = NULL;    nc4->opposite[1] = NULL;
    nc4->opposite[2] = nc1;     nc4->opposite[3] = nc3;

    puts("\nFLIP44\nbefore:");
    print_ncells(&s3);
    // double ranges3[6] = {-2, 2, -2, 2, -2, 2};
    // plot_ncell_3d(&s3, s3.head, "flip44/1", ranges3);
    // plot_ncell_3d(&s3, s3.head->next, "flip44/2", ranges3);
    // plot_ncell_3d(&s3, s3.head->next->next, "flip44/3", ranges3);
    // plot_ncell_3d(&s3, s3.head->next->next->next, "flip44/4", ranges3);
    // plot_dt_3d(&s3, "flip44/before", ranges3);
    
    flip44(&s3, stack, nc1, 0, 2);
    puts("after:");
    print_ncells(&s3);
    // plot_dt_3d(&s3, "flip44/after", ranges3);
    // plot_ncell_3d(&s3, s3.head, "flip44/after1", ranges3);
    // plot_ncell_3d(&s3, s3.head->next, "flip44/after2", ranges3);
    // plot_ncell_3d(&s3, s3.head->next->next, "flip44/after3", ranges3);
    // plot_ncell_3d(&s3, s3.head->next->next->next, "flip44/after4", ranges3);



    // CONSTRUCT_DT_3D
    puts("\n-------------- TESTING DT ---------------");
    print_matrix(p2, 5, 3);
    printf("IS IN GENERAL POSITION: %d\n", are_in_general_position_3d(p2, 5));
    s_setup *dt_setup = construct_dt_3d(p2, 5);
    print_ncells(dt_setup);
    printf("IS DELAUNAY: %d\n", is_delaunay_3d(dt_setup));
    printf("N POINTS : %d\n", dt_setup->N_points);
    
    // double ranges[6] = {-1.5, 1.5, -1.5, 1.5, -1.5, 1.5};
    // plot_ncell_3d(dt_setup, dt_setup->head, "toy/test", ranges);
    // plot_ncell_3d(dt_setup, dt_setup->head->next, "toy/test_1", ranges);
    // plot_ncell_3d(dt_setup, dt_setup->head->next->next, "toy/test_2", ranges);
    // plot_dt_3d(dt_setup,  "toy/test_3", ranges);
    //

    // TODO!! STILL CANNOT TRIANGULATE SETS THAT ARE NOT IN GENERAL POSITION
    // puts("\n-------------- TESTING DT 2 ---------------");
    // print_matrix(p3, 5, 3);
    // printf("IS IN GENERAL POSITION: %d\n", are_in_general_position_3d(p3, 5));
    // dt_setup = construct_dt_3d(p3, 5);
    // print_ncells(dt_setup);
    // printf("IS DELAUNAY: %d\n", is_delaunay_3d(dt_setup));
    // printf("N POINTS : %d\n", dt_setup->N_points);
    // 
    // double ranges[6] = {-1.5, 1.5, -1.5, 1.5, -1.5, 1.5};
    // plot_ncell_3d(dt_setup, dt_setup->head, "toy3/test", ranges);
    // plot_ncell_3d(dt_setup, dt_setup->head->next, "toy3/test_1", ranges);
    // plot_dt_3d(dt_setup,  "toy3/test_3", ranges);


    // TEST AGAIN WITH MORE POINTS
    puts("\n-------------------------------");
    puts("TESTING AGAIN WITH MORE POINTS");
    double **p2_2 = malloc_matrix(6, 3);
    p2_2[0][0] = 1;    p2_2[0][1] = 0;    p2_2[0][2] = 0;
    p2_2[1][0] = -1;   p2_2[1][1] = 1;    p2_2[1][2] = 0;
    p2_2[2][0] = 0;    p2_2[2][1] = -1;   p2_2[2][2] = 0;
    p2_2[3][0] = 0;    p2_2[3][1] = 0;    p2_2[3][2] = 1;
    p2_2[4][0] = 0;    p2_2[4][1] = 0;    p2_2[4][2] = -1;
    p2_2[5][0] = 2;    p2_2[5][1] = 2;    p2_2[5][2] = 2;

    printf("\n\nIS IN GENERAL POSITION: %d\n", are_in_general_position_3d(p2_2, 6));

    s_setup *dt_setup_2 = construct_dt_3d(p2_2, 6);
    print_ncells(dt_setup_2);
    printf("IS DELAUNAY: %d\n", is_delaunay_3d(dt_setup_2));

    // double ranges_2[6] = {-2.5, 2.5, -2.5, 2.5, -2.5, 2.5};
    // plot_dt_3d(dt_setup_2,  "toy2/test", ranges_2);
    // plot_ncell_3d(dt_setup_2, dt_setup_2->head, "toy2/test_1", ranges_2);
    // plot_ncell_3d(dt_setup_2, dt_setup_2->head->next,  "toy2/test_2", ranges_2);
    // plot_ncell_3d(dt_setup_2, dt_setup_2->head->next->next, "toy2/test_3", ranges_2);
    // plot_ncell_3d(dt_setup_2, dt_setup_2->head->next->next->next, "toy2/test_4", ranges_2);
    // plot_ncell_3d(dt_setup_2, dt_setup_2->head->next->next->next->next, "toy2/test_5", ranges_2);
    // plot_ncell_3d(dt_setup_2, dt_setup_2->head->next->next->next->next->next, "toy2/test_6", ranges_2);
    
    puts("--------------------------------------\n");

    // ARE IN GENERAL POSITION, BY CREATING A CO_SPHERICAL SET OF POINTS
    double **p4 = malloc_matrix(5, 3);
    p4[0][0] = 1;    p4[0][1] = 0;    p4[0][2] = 0;
    p4[1][0] = 0;    p4[1][1] = 1;    p4[1][2] = 0;
    p4[2][0] = 0;    p4[2][1] = -1;   p4[2][2] = 0;
    p4[3][0] = 0;    p4[3][1] = 0;    p4[3][2] = 1;
    p4[4][0] = 0;    p4[4][1] = -1.0/sqrt(2);    p4[4][2] = -1.0/sqrt(2);
    printf("\nIN GENERAL POSITION TEST = %d (should be 0)\n", are_in_general_position_3d(p3, 5));


    // TESTING
    int N = 10;
    double **pp = malloc_matrix(N, 3);
    random_points_3d(N, pp);
    printf("\nRANDOM POINTS, IN GENERAL POSITION: %d\n", are_in_general_position_3d(pp, N));
    s_setup *ss = construct_dt_3d(pp, N);
    printf("RANDOM SET OF POINTS, IS DELAUNAY: %d\n", is_delaunay_3d(ss));
    double ranges2[6] = {0, 1, 0, 1, 0, 1};
    plot_dt_3d(ss, "random_example/out", ranges2);
    
}
