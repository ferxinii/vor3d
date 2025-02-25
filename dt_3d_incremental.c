// TODO Make plots to visualize what is happening? Or simply print the values... TEST!!

#include <assert.h>
#include "simplical_complex.c"


#define MAX_STACK_SIZE 10000

typedef struct stack {
    s_ncell *entry[MAX_STACK_SIZE];
    int size;
} s_stack;


s_stack *stack_create(void)
{
    s_stack *stack = malloc(sizeof(s_stack));
    stack->size = 0;
    return stack;
}


void stack_push(s_stack *stack, s_ncell *ncell)
{
    assert(stack->size < MAX_STACK_SIZE && "Reached stack limit.");
    stack->entry[stack->size] = ncell;
    stack->size++;
}


s_ncell *stack_pop(s_stack *stack)
{   
    if (stack->size == 0) return NULL;
    stack->size--;
    s_ncell *ncell = stack->entry[stack->size];
    return ncell;
}


void stack_print(s_stack *stack)
{
    puts("STACK:");
    for (int ii=0; ii<stack->size; ii++) {
        printf("%d, %p\n", ii, (void*)stack->entry[ii]);
    }
    puts("");
}


void flip14(s_setup *setup, s_ncell *container_ncell, int point_id, s_stack *stack)
{   
    setup->N_ncells += 3;

    s_ncell *nc2 = malloc_ncell(setup);
    s_ncell *nc3 = malloc_ncell(setup);
    s_ncell *nc4 = malloc_ncell(setup);
    
    
    // Update linked-list of ncells
    nc4->next = container_ncell->next;
    container_ncell->next = nc2;
    nc2->next = nc3;
    nc3->next = nc4;
    
    s_ncell *opposite_aux[4] = {container_ncell->opposite[0],
                                container_ncell->opposite[1],
                                container_ncell->opposite[2],
                                container_ncell->opposite[3]};
    int v_aux[4] = {container_ncell->vertex_id[0], 
                    container_ncell->vertex_id[1],
                    container_ncell->vertex_id[2],
                    container_ncell->vertex_id[3]};
    
    container_ncell->vertex_id[0] = point_id;
    container_ncell->opposite[1] = nc2;
    container_ncell->opposite[2] = nc3;
    container_ncell->opposite[3] = nc4;


    int opp_id, aux;
    nc2->vertex_id[0] = v_aux[0];
    nc2->vertex_id[1] = point_id;
    nc2->vertex_id[2] = v_aux[2];
    nc2->vertex_id[3] = v_aux[3];
    nc2->opposite[0] = container_ncell;
    nc2->opposite[1] = opposite_aux[1];
    nc2->opposite[2] = nc3;
    nc2->opposite[3] = nc4;

    if (opposite_aux[1]) {  // This is untested? TODO
        aux = 1;
        face_localid_of_adjacent_ncell(setup, nc2, &aux, 2, 1, &opp_id); 
        opposite_aux[1]->opposite[opp_id] = nc2;
    }

    nc3->vertex_id[0] = v_aux[0];
    nc3->vertex_id[1] = v_aux[1];
    nc3->vertex_id[2] = point_id;
    nc3->vertex_id[3] = v_aux[3];
    nc3->opposite[0] = container_ncell;
    nc3->opposite[1] = nc2;
    nc3->opposite[2] = opposite_aux[2];
    nc3->opposite[3] = nc4;
    if (opposite_aux[2]) {
        aux = 2;
        face_localid_of_adjacent_ncell(setup, nc3, &aux, 2, 2, &opp_id); 
        (nc3->opposite[2])->opposite[opp_id] = nc3;
    }

    nc4->vertex_id[0] = v_aux[0];
    nc4->vertex_id[1] = v_aux[1];
    nc4->vertex_id[2] = v_aux[2];
    nc4->vertex_id[3] = point_id;
    nc4->opposite[0] = container_ncell;
    nc4->opposite[1] = nc2;
    nc4->opposite[2] = nc3;
    nc4->opposite[3] = opposite_aux[3];
    if (opposite_aux[3]) {
        aux = 3;
        face_localid_of_adjacent_ncell(setup, nc4, &aux, 2, 3, &opp_id); 
        (nc4->opposite[3])->opposite[opp_id] = nc4;
    }

    stack_push(stack, container_ncell);
    stack_push(stack, nc2);
    stack_push(stack, nc3);
    stack_push(stack, nc4);
}


void flip23(s_setup *setup, s_stack *stack, s_ncell *nc1, int opp_cell_id, int opp_face_localid)
{   
    setup->N_ncells += 1;

    s_ncell *nc2 = nc1->opposite[opp_cell_id];
    s_ncell *nc3 = malloc_ncell(setup);

    // Update linked-list of ncells
    nc3->next = nc2->next;
    nc2->next = nc3;

    int face_vertex_id[3];
    extract_ids_face(setup, nc1, &opp_cell_id, 2, face_vertex_id);
    int a = face_vertex_id[0];
    int b = face_vertex_id[1];
    int c = face_vertex_id[2];
    int d = nc2->vertex_id[opp_face_localid];
    int p = nc1->vertex_id[opp_cell_id];

    int nc1_vertex_id_old[4] = {nc1->vertex_id[0], nc1->vertex_id[1], nc1->vertex_id[2], nc1->vertex_id[3]};
    int nc2_vertex_id_old[4] = {nc2->vertex_id[0], nc2->vertex_id[1], nc2->vertex_id[2], nc2->vertex_id[3]};
    s_ncell *opp_1_old[4] = {nc1->opposite[0], nc1->opposite[1], nc1->opposite[2], nc1->opposite[3]};
    s_ncell *opp_2_old[4] = {nc2->opposite[0], nc2->opposite[1], nc2->opposite[2], nc2->opposite[3]};
    
    s_ncell *nc_aux;
    int aux, aux2, v_localid_opp;
    
    // NC1
    nc1->vertex_id[id_where_equal_int(nc1_vertex_id_old, 4, c)] = d;
    nc1->opposite[id_where_equal_int(nc1_vertex_id_old, 4, a)] = nc2;
    nc1->opposite[id_where_equal_int(nc1_vertex_id_old, 4, b)] = nc3;
    nc1->opposite[id_where_equal_int(nc1_vertex_id_old, 4, p)] = opp_2_old[id_where_equal_int(nc2_vertex_id_old, 4, c)];

    aux = id_where_equal_int(nc2_vertex_id_old, 4, c); 
    nc_aux = opp_2_old[aux];
    if (nc_aux) {
        aux2 = id_where_equal_int(nc1_vertex_id_old, 4, p);
        face_localid_of_adjacent_ncell(setup, nc1, &aux2, 2, aux2, &v_localid_opp);
        face_localid_of_adjacent_ncell(setup, nc1, &aux, 2, aux, &v_localid_opp);
        nc_aux->opposite[v_localid_opp] = nc1;
    }


    
    // NC2
    nc2->vertex_id[id_where_equal_int(nc2_vertex_id_old, 4, a)] = p;
    nc2->opposite[id_where_equal_int(nc2_vertex_id_old, 4, b)] = nc3;
    nc2->opposite[id_where_equal_int(nc2_vertex_id_old, 4, c)] = nc1;
    nc2->opposite[id_where_equal_int(nc2_vertex_id_old, 4, d)] = opp_1_old[id_where_equal_int(nc1_vertex_id_old, 4, a)];

    aux = id_where_equal_int(nc1_vertex_id_old, 4, a);
    nc_aux = opp_1_old[aux];
    if (nc_aux) {
        face_localid_of_adjacent_ncell(setup, nc2, &aux, 2, aux, &v_localid_opp);
        nc_aux->opposite[v_localid_opp] = nc2;
    }

    // NC3
    nc3->vertex_id[0] = p;  nc3->vertex_id[1] = c;  nc3->vertex_id[2] = d;  nc3->vertex_id[3] = a;
    nc3->opposite[0] = opp_2_old[id_where_equal_int(nc2_vertex_id_old, 4, b)];
    nc3->opposite[1] = nc1;
    nc3->opposite[2] = opp_1_old[id_where_equal_int(nc1_vertex_id_old, 4, b)];
    nc3->opposite[3] = nc2;

    aux = id_where_equal_int(nc1_vertex_id_old, 4, b);
    nc_aux = opp_1_old[aux];
    if (nc_aux) {
        face_localid_of_adjacent_ncell(setup, nc3, &aux, 2, aux, &v_localid_opp);
        nc_aux->opposite[v_localid_opp] = nc3;
    }

    aux = id_where_equal_int(nc2_vertex_id_old, 4, b);
    nc_aux = opp_2_old[aux];
    if (nc_aux) {
        face_localid_of_adjacent_ncell(setup, nc3, &aux, 2, aux, &v_localid_opp);
        nc_aux->opposite[v_localid_opp] = nc3;
    }

    stack_push(stack, nc1);
    stack_push(stack, nc2);
    stack_push(stack, nc3);
}


int AUX_can_perform_flip_32(const s_setup *setup, const s_ncell *ncell, int opp_cell_id, int *ridge_id_2)
{
    if (setup->N_ncells < 3) return 0;

    int face_vertex_id[3];
    extract_ids_face(setup, ncell, &opp_cell_id, 2, face_vertex_id);

    int num_cells_ridge;
    
    num_cells_ridge = count_cycle_ridge(setup, ncell, opp_cell_id, id_where_equal_int(ncell->vertex_id, 4, face_vertex_id[0]));
    if (num_cells_ridge == 3) {
        *ridge_id_2 = id_where_equal_int(face_vertex_id, 3, face_vertex_id[0]);
        return 1;
    }

    num_cells_ridge = count_cycle_ridge(setup, ncell, opp_cell_id, id_where_equal_int(ncell->vertex_id, 4, face_vertex_id[1]));
    if (num_cells_ridge == 3) {
        *ridge_id_2 = id_where_equal_int(face_vertex_id, 3, face_vertex_id[1]);
        return 1;
    }

    num_cells_ridge = count_cycle_ridge(setup, ncell, opp_cell_id, id_where_equal_int(ncell->vertex_id, 4, face_vertex_id[2]));
    if (num_cells_ridge == 3) {
        *ridge_id_2 = id_where_equal_int(face_vertex_id, 3, face_vertex_id[2]);
        return 1;
    }
    
    return 0;
}


void flip32(s_setup *setup, s_stack *stack, s_ncell *nc1, int opp_cell_id, int ridge_id_2, int opp_face_localid)
{
    setup->N_ncells -= 1;

    int v2_main, v2_2, v3_main, v3_2;
    s_ncell *nc2 = next_ncell_ridge_cycle(setup, nc1, opp_cell_id, ridge_id_2, &v2_main, &v2_2);
    s_ncell *nc3 = next_ncell_ridge_cycle(setup, nc2, v2_main, v2_2, &v3_main, &v3_2);
    
    nc2->next = nc3->next;

    int localids_ridge[2] = {opp_cell_id, ridge_id_2};
    int vertices_ridge[2];
    extract_ids_face(setup, nc1, localids_ridge, 1, vertices_ridge);

    int p = nc1->vertex_id[opp_cell_id];
    int a = nc1->vertex_id[ridge_id_2];
    int b = vertices_ridge[0];
    int c = vertices_ridge[1];
    int d = nc2->vertex_id[opp_face_localid];

    int nc1_vertex_id_old[4] = {nc1->vertex_id[0], nc1->vertex_id[1], nc1->vertex_id[2], nc1->vertex_id[3]};
    int nc2_vertex_id_old[4] = {nc2->vertex_id[0], nc2->vertex_id[1], nc2->vertex_id[2], nc2->vertex_id[3]};
    int nc3_vertex_id_old[4] = {nc3->vertex_id[0], nc3->vertex_id[1], nc3->vertex_id[2], nc3->vertex_id[3]};
    s_ncell *opp_1_old[4] = {nc1->opposite[0], nc1->opposite[1], nc1->opposite[2], nc1->opposite[3]};
    s_ncell *opp_2_old[4] = {nc2->opposite[0], nc2->opposite[1], nc2->opposite[2], nc2->opposite[3]};
    s_ncell *opp_3_old[4] = {nc3->opposite[0], nc3->opposite[1], nc3->opposite[2], nc3->opposite[3]};
    
    s_ncell *nc_aux;
    int aux, v_localid_opp;

    // NC1
    nc1->vertex_id[id_where_equal_int(nc1_vertex_id_old, 4, c)] = d;
    nc1->opposite[id_where_equal_int(nc1_vertex_id_old, 4, b)] = nc2;
    nc1->opposite[id_where_equal_int(nc1_vertex_id_old, 4, p)] = opp_2_old[id_where_equal_int(nc2_vertex_id_old, 4, c)];
    nc1->opposite[id_where_equal_int(nc1_vertex_id_old, 4, a)] = opp_3_old[id_where_equal_int(nc3_vertex_id_old, 4, c)];
    
    aux = id_where_equal_int(nc2_vertex_id_old, 4, c);
    nc_aux = opp_2_old[aux];
    if (nc_aux) {
        face_localid_of_adjacent_ncell(setup, nc1, &aux, 2, aux, &v_localid_opp);
        nc_aux->opposite[v_localid_opp] = nc1;
    }

    aux = id_where_equal_int(nc3_vertex_id_old, 4, c);
    nc_aux = opp_3_old[aux];
    if (nc_aux) {
        face_localid_of_adjacent_ncell(setup, nc1, &aux, 2, aux, &v_localid_opp);
        nc_aux->opposite[v_localid_opp] = nc1;
    }

    // NC2
    nc2->vertex_id[id_where_equal_int(nc2_vertex_id_old, 4, b)] = p; 
    nc2->opposite[id_where_equal_int(nc2_vertex_id_old, 4, c)] = nc1;
    nc2->opposite[id_where_equal_int(nc2_vertex_id_old, 4, d)] = opp_1_old[id_where_equal_int(nc1_vertex_id_old, 4, b)];
    nc2->opposite[id_where_equal_int(nc2_vertex_id_old, 4, a)] = opp_3_old[id_where_equal_int(nc3_vertex_id_old, 4, b)];

    aux = id_where_equal_int(nc1_vertex_id_old, 4, b);
    nc_aux = opp_1_old[aux];
    if (nc_aux) {
        face_localid_of_adjacent_ncell(setup, nc2, &aux, 2, aux, &v_localid_opp);
        nc_aux->opposite[v_localid_opp] = nc2;
    }

    aux = id_where_equal_int(nc3_vertex_id_old, 4, b);
    nc_aux = opp_3_old[aux];
    if (nc_aux) {
        face_localid_of_adjacent_ncell(setup, nc2, &aux, 2, aux, &v_localid_opp);
        nc_aux->opposite[v_localid_opp] = nc2;
    }

    free_ncell(nc3);

    stack_push(stack, nc1);
    stack_push(stack, nc2);
}


void flip_tetrahedra(s_setup *setup, s_stack *stack, s_ncell *ncell, int opp_cell_id)
{
    print_ncell(setup, ncell);
    printf("Opp_cell_id: %d\n", opp_cell_id);
    print_ncell(setup, ncell->opposite[opp_cell_id]);

    static double **coords_face = NULL;
    if (!coords_face) {
        coords_face = malloc_matrix(3, 3);
    }
    extract_vertices_face(setup, ncell, &opp_cell_id, 2, coords_face);

    // Extract id of vertex in opposite cell corresponding to the face
    int opp_face_localid;
    face_localid_of_adjacent_ncell(setup, ncell, &opp_cell_id, 2, opp_cell_id, &opp_face_localid);
    int opp_face_vertex_id = (ncell->opposite[opp_cell_id])->vertex_id[opp_face_localid];

    printf("Opp_face_localid: %d\n", opp_face_localid);

    double *p = setup->points[ncell->vertex_id[opp_cell_id]];
    double *d = setup->points[opp_face_vertex_id];

    if (segment_crosses_triangle_3d(coords_face, p, d) == 1) {  // Convex polygon!
        puts("flip23");
        flip23(setup, stack, ncell, opp_cell_id, opp_face_localid);
        print_ncell(setup, stack->entry[stack->size-1]);
        print_ncell(setup, stack->entry[stack->size-2]);
        print_ncell(setup, stack->entry[stack->size-3]);
    } else {
        int ridge_id_2;
        if (AUX_can_perform_flip_32(setup, ncell, opp_cell_id, &ridge_id_2)) {
            puts("flip32");
            flip32(setup, stack, ncell, opp_cell_id, ridge_id_2, opp_face_localid);
        }
    }

    return;
}


void insert_one_point(s_setup *setup, int point_id, s_stack *stack)
{
    double *point = setup->points[point_id];
    s_ncell *container_ncell = in_ncell_walk(setup, point);

    // Insert p in container_ncell with a flip14
    puts("flip14");
    flip14(setup, container_ncell, point_id, stack);

    while (stack->size > 0) {
        s_ncell *current = stack_pop(stack);

        int opp_cell_id = id_where_equal_int(current->vertex_id, 4, point_id);
        if (current->opposite[opp_cell_id]) {
            if (are_locally_delaunay(setup, current, opp_cell_id) != 1) {
                flip_tetrahedra(setup, stack, current, opp_cell_id);
            }
        }
    }
}


s_setup *construct_dt_3d(double **points, int N_points)
{
    s_stack *stack = stack_create();
    s_setup *setup = initialize_setup(points, N_points, 3);

    for (int ii=0; ii<N_points; ii++) {
        insert_one_point(setup, ii, stack);
    }

    return setup;
}
