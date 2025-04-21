// TODO Make plots to visualize what is happening? Or simply print the values... TEST!!

#include <assert.h>
#include "simplical_complex.c"


// ---------------------------------------------------------------------------------------
// ----------------------------------- STACK ---------------------------------------------
// ---------------------------------------------------------------------------------------

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


void stack_free(s_stack *stack)
{
    free(stack->entry);
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


void stack_remove_ncell(s_stack *stack, s_ncell *ncell) {
    int newSize = 0;
    // Iterate over all entries in the stack.
    for (int ii = 0; ii < stack->size; ii++) {
        // Only copy entries that are not the target.
        if (stack->entry[ii] != ncell) {
            stack->entry[newSize++] = stack->entry[ii];
        }
    }
    stack->size = newSize;
}


void stack_print(s_stack *stack)
{
    puts("STACK:");
    for (int ii=0; ii<stack->size; ii++) {
        printf("%d, %p\n", ii, (void*)stack->entry[ii]);
    }
    puts("");
}


// ---------------------------------------------------------------------------------------
// ----------------------------------- FLIPS ---------------------------------------------
// ---------------------------------------------------------------------------------------

void flip14(s_setup *setup, s_ncell *container_ncell, int point_id, s_stack *stack)
{   
    setup->N_ncells += 3;

    s_ncell *nc2 = malloc_ncell(setup);
    s_ncell *nc3 = malloc_ncell(setup);
    s_ncell *nc4 = malloc_ncell(setup);
    
    // Update linked-list of ncells
    // ---- NC1 ----------------------------
    // ---- NC1 --- NC2 --- NC3 --- NC4 ----
    nc4->next = container_ncell->next;
    if (nc4->next) (nc4->next)->prev = nc4;

    container_ncell->next = nc2;
    nc2->prev = container_ncell;

    nc2->next = nc3;
    nc3->prev = nc2;

    nc3->next = nc4;
    nc4->prev = nc3;
    

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


void flip23(s_setup *setup, s_stack *stack, s_ncell *nc1, int opp_cell_id, int opp_face_localid, s_ncell **OUT_PTRS)
{   
    setup->N_ncells += 1;

    s_ncell *nc2 = nc1->opposite[opp_cell_id];
    s_ncell *nc3 = malloc_ncell(setup);

    // Update linked-list of ncells
    // ---- NC1 ----- NC2 ---------------- 
    // ---- NC1 ----- NC2 ----- NC3 ------
    nc3->next = nc2->next;
    nc3->prev = nc2;
    if (nc3->next) nc3->next->prev = nc3;

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
        aux2 = id_where_equal_int(nc1->vertex_id, 4, p);
        face_localid_of_adjacent_ncell(setup, nc1, &aux2, 2, aux2, &v_localid_opp);
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
        aux2 = id_where_equal_int(nc2->vertex_id, 4, d);
        face_localid_of_adjacent_ncell(setup, nc2, &aux2, 2, aux2, &v_localid_opp);
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
        aux2 = id_where_equal_int(nc3->vertex_id, 4, d);
        face_localid_of_adjacent_ncell(setup, nc3, &aux2, 2, aux2, &v_localid_opp);
        nc_aux->opposite[v_localid_opp] = nc3;
    }

    aux = id_where_equal_int(nc2_vertex_id_old, 4, b);
    nc_aux = opp_2_old[aux];
    if (nc_aux) {
        aux2 = id_where_equal_int(nc3->vertex_id, 4, p);
        face_localid_of_adjacent_ncell(setup, nc3, &aux2, 2, aux2, &v_localid_opp);
        nc_aux->opposite[v_localid_opp] = nc3;
    }
    
    stack_push(stack, nc1);
    stack_push(stack, nc2);
    stack_push(stack, nc3);
    if (OUT_PTRS) {
        OUT_PTRS[0] = nc1;
        OUT_PTRS[1] = nc2;
        OUT_PTRS[2] = nc3;
    }
}


int can_perform_flip32(const s_setup *setup, const s_ncell *ncell, int opp_cell_id, int *ridge_id_2)
{
    if (setup->N_ncells < 3) return 0;

    int face_vertex_id[3];
    extract_ids_face(setup, ncell, &opp_cell_id, 2, face_vertex_id);
    
    int opp_face_localid;
    face_localid_of_adjacent_ncell(setup, ncell, &opp_cell_id, 2, opp_cell_id, &opp_face_localid);
    
    s_ncell *opp_ncell = ncell->opposite[opp_cell_id];
    for (int ii=0; ii<4; ii++) {
        s_ncell *opp_opp = opp_ncell->opposite[ii];
        if (opp_opp && 
            opp_opp != ncell &&
            inarray(opp_opp->vertex_id, 4, ncell->vertex_id[opp_cell_id]) &&
            inarray(opp_opp->vertex_id, 4, opp_ncell->vertex_id[opp_face_localid])) {
                *ridge_id_2 = id_where_equal_int(ncell->vertex_id, 4, opp_ncell->vertex_id[ii]);
                return 1;
            }
    }
    return 0;
}


void flip32(s_setup *setup, s_stack *stack, s_ncell *nc1, int opp_cell_id, int ridge_id_2, int opp_face_localid, s_ncell **OUT_PTRS)
{
    setup->N_ncells -= 1;

    int v2_main, v2_2, v3_main, v3_2;
    s_ncell *nc2 = next_ncell_ridge_cycle(setup, nc1, opp_cell_id, ridge_id_2, &v2_main, &v2_2);
    s_ncell *nc3 = next_ncell_ridge_cycle(setup, nc2, v2_main, v2_2, &v3_main, &v3_2);

    // ------ NC3->PREV ----- NC3 ------- NC3->NEXT ------
    s_ncell *nc3_next = nc3->next;
    if (nc3->next) nc3->next->prev = nc3->prev;
    if (nc3->prev) nc3->prev->next = nc3_next;
    else setup->head = nc3->next;


    int localids_ridge[2] = {opp_cell_id, ridge_id_2};
    int vertices_ridge[2];
    extract_ids_face(setup, nc1, localids_ridge, 1, vertices_ridge);

    int p = nc1->vertex_id[opp_cell_id];
    int a = nc1->vertex_id[ridge_id_2];
    int b = vertices_ridge[0];
    int c = vertices_ridge[1];
    int d = nc2->vertex_id[opp_face_localid];

    // WE NEED TO MAKE SURE THEY ARE ORIENTED ??? TODO NOT SURE ABOUT THIS, IS IT NECESSARY?
    static double *face_vertices[3];
    face_vertices[0] = setup->points[a];
    face_vertices[1] = setup->points[b];
    face_vertices[2] = setup->points[c];
    if ( orientation(face_vertices, setup->points[nc1->vertex_id[opp_cell_id]], 3) == -1 ) {
        int temp = b;
        b = c;
        c = temp;
    }
    face_vertices[0] = setup->points[a];
    face_vertices[1] = setup->points[b];
    face_vertices[2] = setup->points[c];
    assert(orientation(face_vertices, setup->points[nc1->vertex_id[opp_cell_id]], 3) == 1);

    int nc1_vertex_id_old[4] = {nc1->vertex_id[0], nc1->vertex_id[1], nc1->vertex_id[2], nc1->vertex_id[3]};
    int nc2_vertex_id_old[4] = {nc2->vertex_id[0], nc2->vertex_id[1], nc2->vertex_id[2], nc2->vertex_id[3]};
    int nc3_vertex_id_old[4] = {nc3->vertex_id[0], nc3->vertex_id[1], nc3->vertex_id[2], nc3->vertex_id[3]};
    s_ncell *opp_1_old[4] = {nc1->opposite[0], nc1->opposite[1], nc1->opposite[2], nc1->opposite[3]};
    s_ncell *opp_2_old[4] = {nc2->opposite[0], nc2->opposite[1], nc2->opposite[2], nc2->opposite[3]};
    s_ncell *opp_3_old[4] = {nc3->opposite[0], nc3->opposite[1], nc3->opposite[2], nc3->opposite[3]};
    
    s_ncell *nc_aux;
    int aux, aux2, v_localid_opp;

    // NC1
    nc1->vertex_id[id_where_equal_int(nc1_vertex_id_old, 4, c)] = d;
    nc1->opposite[id_where_equal_int(nc1_vertex_id_old, 4, b)] = nc2;
    nc1->opposite[id_where_equal_int(nc1_vertex_id_old, 4, p)] = opp_2_old[id_where_equal_int(nc2_vertex_id_old, 4, c)];
    nc1->opposite[id_where_equal_int(nc1_vertex_id_old, 4, a)] = opp_3_old[id_where_equal_int(nc3_vertex_id_old, 4, c)];
    
    aux = id_where_equal_int(nc2_vertex_id_old, 4, c);
    nc_aux = opp_2_old[aux];
    if (nc_aux) {
        aux2 = id_where_equal_int(nc1->vertex_id, 4, p);
        assert(aux2 == opp_cell_id);
        face_localid_of_adjacent_ncell(setup, nc1, &aux2, 2, aux2, &v_localid_opp);
        nc_aux->opposite[v_localid_opp] = nc1;
    }

    aux = id_where_equal_int(nc3_vertex_id_old, 4, c);
    nc_aux = opp_3_old[aux];
    if (nc_aux) {
        aux2 = id_where_equal_int(nc1->vertex_id, 4, a);
        face_localid_of_adjacent_ncell(setup, nc1, &aux2, 2, aux2, &v_localid_opp);
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
        aux2 = id_where_equal_int(nc2->vertex_id, 4, d);
        face_localid_of_adjacent_ncell(setup, nc2, &aux2, 2, aux2, &v_localid_opp);
        nc_aux->opposite[v_localid_opp] = nc2;
    }

    aux = id_where_equal_int(nc3_vertex_id_old, 4, b);
    nc_aux = opp_3_old[aux];
    if (nc_aux) {
        aux2 = id_where_equal_int(nc2->vertex_id, 4, a);
        face_localid_of_adjacent_ncell(setup, nc2, &aux2, 2, aux2, &v_localid_opp);
        nc_aux->opposite[v_localid_opp] = nc2;
    }
    
    stack_remove_ncell(stack, nc3);

    free_ncell(nc3);

    stack_push(stack, nc1);
    stack_push(stack, nc2);
    if (OUT_PTRS) {
        OUT_PTRS[0] = nc1;
        OUT_PTRS[1] = nc2;
    }
}


int can_perform_flip44(const s_setup *setup, const s_ncell *ncell, double **vertices_face, int opp_cell_id, int *ridge_id_2)
{
    if (setup->N_ncells < 4) return 0;
    
    s_ncell *opp = ncell->opposite[opp_cell_id];
    int opp_face_localid;
    face_localid_of_adjacent_ncell(setup, ncell, &opp_cell_id, 2, opp_cell_id, &opp_face_localid);


    // DETERMINE RIDGE 
    int face_vertex_id[4];
    extract_ids_face(setup, ncell, &opp_cell_id, 2, face_vertex_id);

    double *aux[3];
    aux[0] = setup->points[ncell->vertex_id[opp_cell_id]]; 

    aux[1] = vertices_face[0];
    aux[2] = vertices_face[1];
    if (orientation(aux, setup->points[opp->vertex_id[opp_face_localid]], 3) == 0) {
        *ridge_id_2 = id_where_equal_int(ncell->vertex_id, 4, face_vertex_id[2]);
    }
    aux[1] = vertices_face[1];
    aux[2] = vertices_face[2];
    if (orientation(aux, setup->points[opp->vertex_id[opp_face_localid]], 3) == 0) {
        *ridge_id_2 = id_where_equal_int(ncell->vertex_id, 4, face_vertex_id[0]);
    }
    aux[1] = vertices_face[0];
    aux[2] = vertices_face[2];
    if (orientation(aux, setup->points[opp->vertex_id[opp_face_localid]], 3) == 0) {
        *ridge_id_2 = id_where_equal_int(ncell->vertex_id, 4, face_vertex_id[1]);
    }
    

    int nc2_id1, nc2_id2;
    s_ncell *nc2 = next_ncell_ridge_cycle(setup, ncell, opp_cell_id, *ridge_id_2, &nc2_id1, &nc2_id2);

    int nc3_id1, nc3_id2;
    s_ncell *nc3 = next_ncell_ridge_cycle(setup, nc2, nc2_id1, nc2_id2, &nc3_id1, &nc3_id2);

    int nc4_id1, nc4_id2;
    s_ncell *nc4 = next_ncell_ridge_cycle(setup, nc3, nc3_id1, nc3_id2, &nc4_id1, &nc4_id2);

    if (nc4->opposite[nc4_id1] == ncell) return 1;
    else return 0;
}


void debug_print_44_before(s_setup *setup, s_ncell *ncell, int opp_cell_id, int ridge_id_2, char *file)
{
    int nc2_id1, nc2_id2;
    s_ncell *nc2 = next_ncell_ridge_cycle(setup, ncell, opp_cell_id, ridge_id_2, &nc2_id1, &nc2_id2);

    int nc3_id1, nc3_id2;
    s_ncell *nc3 = next_ncell_ridge_cycle(setup, nc2, nc2_id1, nc2_id2, &nc3_id1, &nc3_id2);

    int nc4_id1, nc4_id2;
    s_ncell *nc4 = next_ncell_ridge_cycle(setup, nc3, nc3_id1, nc3_id2, &nc4_id1, &nc4_id2);

    printf("FLIP44 BEFORE: \n \
            (%d, %d), (%d, %d, %d, %d)\n \
            (%d, %d), (%d, %d, %d, %d)\n \
            (%d, %d), (%d, %d, %d, %d)\n \
            (%d, %d), (%d, %d, %d, %d)\n \n", 
            ncell->vertex_id[opp_cell_id], ncell->vertex_id[ridge_id_2], 
             ncell->vertex_id[0], ncell->vertex_id[1], ncell->vertex_id[2], ncell->vertex_id[3], 
            nc2->vertex_id[nc2_id1], nc2->vertex_id[nc2_id2], 
            nc2->vertex_id[0], nc2->vertex_id[1], nc2->vertex_id[2], nc2->vertex_id[3],
            nc3->vertex_id[nc3_id1],  nc3->vertex_id[nc3_id2],
            nc3->vertex_id[0], nc3->vertex_id[1], nc3->vertex_id[2], nc3->vertex_id[3],
            nc4->vertex_id[nc4_id1],  nc4->vertex_id[nc4_id2],
            nc4->vertex_id[0], nc4->vertex_id[1], nc4->vertex_id[2], nc4->vertex_id[3]);

    FILE *file_debug = fopen(file, "a");
    write_ncell3d_file((s_setup*)setup, (s_ncell*)ncell, file_debug);
    write_ncell3d_file((s_setup*)setup, (s_ncell*)nc2, file_debug);
    write_ncell3d_file((s_setup*)setup, (s_ncell*)nc3, file_debug);
    write_ncell3d_file((s_setup*)setup, (s_ncell*)nc4, file_debug);
    fprintf(file_debug, "\n\n");
    fclose(file_debug);
}


void debug_print_44_after(s_setup *setup, s_ncell **PTRS, char *file)
{
    s_ncell *nc1 = PTRS[0];
    s_ncell *nc2 = PTRS[1];
    s_ncell *nc3 = PTRS[2];
    s_ncell *nc4 = PTRS[3];
    printf("FLIP44 AFTER: \n \
            (%d, %d, %d, %d)\n \
            (%d, %d, %d, %d)\n \
            (%d, %d, %d, %d)\n \
            (%d, %d, %d, %d)\n \n", 
            nc1->vertex_id[0], nc1->vertex_id[1], nc1->vertex_id[2], nc1->vertex_id[3], 
            nc2->vertex_id[0], nc2->vertex_id[1], nc2->vertex_id[2], nc2->vertex_id[3],
            nc3->vertex_id[0], nc3->vertex_id[1], nc3->vertex_id[2], nc3->vertex_id[3],
            nc4->vertex_id[0], nc4->vertex_id[1], nc4->vertex_id[2], nc4->vertex_id[3]);

    FILE *file_debug = fopen(file, "a");
    write_ncell3d_file((s_setup*)setup, (s_ncell*)nc1, file_debug);
    write_ncell3d_file((s_setup*)setup, (s_ncell*)nc2, file_debug);
    write_ncell3d_file((s_setup*)setup, (s_ncell*)nc3, file_debug);
    write_ncell3d_file((s_setup*)setup, (s_ncell*)nc4, file_debug);
    fprintf(file_debug, "\n\n");
    fclose(file_debug);
}



void flip44(s_setup *setup, s_stack *stack, s_ncell *ncell, int id_ridge_1, int id_ridge_2, s_ncell **OUT_PTRS) 
{
    // FILE *f = fopen("flip44_pre.txt", "w");
    // fclose(f);
    // debug_print_44_before(setup, ncell, id_ridge_1, id_ridge_2, "flip44_pre.txt");

    // FIRST, A FLIP23
    int opp_face_localid;
    face_localid_of_adjacent_ncell(setup, ncell, &id_ridge_1, 2, id_ridge_1, &opp_face_localid);
    int opp_face_vertexid = ncell->opposite[id_ridge_1]->vertex_id[opp_face_localid];  // Store before flip23!

    int id_a, id_c;
    for (int ii=0; ; ii++) {
        if (ii != id_ridge_1 && ii != id_ridge_2) {
            id_a = ii;
            break;
        }
    }
    for (int ii=0; ; ii++) {
        if (ii != id_ridge_1 && ii != id_ridge_2 && ii != id_a) {
            id_c = ii;
            break;
        }
    }
    int p = ncell->vertex_id[id_ridge_1];
    int a = ncell->vertex_id[id_a];
    int c = ncell->vertex_id[id_c];
    int d = opp_face_vertexid;

    s_ncell *FLIP23_PTRS[3];
    flip23(setup, stack, ncell, id_ridge_1, opp_face_localid, FLIP23_PTRS);  // TOWARDS NC2
    s_ncell *nc5;
    if (inarray(FLIP23_PTRS[0]->vertex_id, 4, a) && inarray(FLIP23_PTRS[0]->vertex_id, 4, c) &&
        inarray(FLIP23_PTRS[0]->vertex_id, 4, d) && inarray(FLIP23_PTRS[0]->vertex_id, 4, p)) {
        nc5 = FLIP23_PTRS[0];
        stack_push(stack, FLIP23_PTRS[1]);
        stack_push(stack, FLIP23_PTRS[2]);
        if (OUT_PTRS) {
            OUT_PTRS[0] = FLIP23_PTRS[1];
            OUT_PTRS[1] = FLIP23_PTRS[2];
        }
    } else if (inarray(FLIP23_PTRS[1]->vertex_id, 4, a) && inarray(FLIP23_PTRS[1]->vertex_id, 4, c) &&
        inarray(FLIP23_PTRS[1]->vertex_id, 4, d) && inarray(FLIP23_PTRS[1]->vertex_id, 4, p)) {
        nc5 = FLIP23_PTRS[1];
        stack_push(stack, FLIP23_PTRS[0]);
        stack_push(stack, FLIP23_PTRS[2]);
        if (OUT_PTRS) {
            OUT_PTRS[0] = FLIP23_PTRS[0];
            OUT_PTRS[1] = FLIP23_PTRS[2];
        }
    } else if (inarray(FLIP23_PTRS[2]->vertex_id, 4, a) && inarray(FLIP23_PTRS[2]->vertex_id, 4, c) &&
        inarray(FLIP23_PTRS[2]->vertex_id, 4, d) && inarray(FLIP23_PTRS[2]->vertex_id, 4, p)) {
        nc5 = FLIP23_PTRS[2];
        stack_push(stack, FLIP23_PTRS[0]);
        stack_push(stack, FLIP23_PTRS[1]);
        if (OUT_PTRS) {
            OUT_PTRS[0] = FLIP23_PTRS[0];
            OUT_PTRS[1] = FLIP23_PTRS[1];
        }
    } else { 
        puts("SHOULD NOT BE HERE! FLIP44 COULD NOT FIND NC5");
        exit(1);
    }

    // NEXT, A FLIP32
    int nc5_p = id_where_equal_int(nc5->vertex_id, 4, ncell->vertex_id[id_ridge_1]);
    
    s_ncell *nc3 = nc5->opposite[id_where_equal_int(nc5->vertex_id, 4, ncell->vertex_id[id_ridge_1])];
    int nc3_id1 = id_where_equal_int(nc3->vertex_id, 4, opp_face_vertexid);
    int nc3_id2;
    face_localid_of_adjacent_ncell(setup, nc5, &nc5_p, 2, nc5_p, &nc3_id2);
    int nc3_opp_face_localid;
    face_localid_of_adjacent_ncell(setup, nc3, &nc3_id1, 2, nc3_id1, &nc3_opp_face_localid);

    s_ncell *FLIP32_PTRS[2];
    flip32(setup, stack, nc3, nc3_id1, nc3_id2, nc3_opp_face_localid, FLIP32_PTRS);

    if (OUT_PTRS) {
        OUT_PTRS[2] = FLIP32_PTRS[0];
        OUT_PTRS[3] = FLIP32_PTRS[1];
    }

    stack_push(stack, FLIP32_PTRS[0]);
    stack_push(stack, FLIP32_PTRS[1]);

    // f = fopen("flip44_after.txt", "w");
    // fclose(f);
    // debug_print_44_after(setup, OUT_PTRS, "flip44_after.txt");
}


// ---------------------------------------------------------------------------------------
// ----------------------------------- CASES ---------------------------------------------
// ---------------------------------------------------------------------------------------


int is_case_3_AUX(double *s1, double *s2, double *p, double *d, int drop_coord)
{
    double ax, ay, bx, by, px, py, dx, dy;
    if (drop_coord == 0) {
      ax = s1[1]; ay = s1[2];
      bx = s2[1]; by = s2[2];
      px = p[1];  py = p[2];
      dx = d[1];  dy = d[2];
    } else if (drop_coord == 1) {
      ax = s1[2]; ay = s1[0];
      bx = s2[2]; by = s2[0];
      px = p[2];  py = p[0];
      dx = d[2];  dy = d[0];
    } else {
      ax = s1[0]; ay = s1[1];
      bx = s2[0]; by = s2[1];
      px = p[0]; py = p[1];
      dx = d[0]; dy = d[1];
    }

    assert(px != dx);
    assert(py != dy);

    double A[2] = {ax, ay},
           B[2] = {bx, by},
           paux[2] = {px, py}, 
           daux[2] = {dx, dy};
    return segments_intersect_2d(A, B, paux, daux);
}


int is_case_3(double **vertices_face, double *p, double *d)
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

    
    double *aux[3];
    aux[0] = p; 

    aux[1] = vertices_face[0];
    aux[2] = vertices_face[1];
    if (orientation(aux, d, 3) == 0 &&
        is_case_3_AUX(vertices_face[0], vertices_face[1], p, d, drop_coord)) return 1;

    aux[1] = vertices_face[1];
    aux[2] = vertices_face[2];
    if (orientation(aux, d, 3) == 0 &&
        is_case_3_AUX(vertices_face[1], vertices_face[2], p, d, drop_coord)) return 1;

    aux[1] = vertices_face[0];
    aux[2] = vertices_face[2];
    if (orientation(aux, d, 3) == 0 &&
        is_case_3_AUX(vertices_face[0], vertices_face[2], p, d, drop_coord)) return 1;
    
    return 0;


}


int is_case_1(double **vertices_face, double *p, double *d) 
{
    if (segment_crosses_triangle_3d(vertices_face, p, d) == 1) return 1;
    return 0;
}


int is_case_2(double **vertices_face, double *p, double *d)
{
    if (segment_crosses_triangle_3d(vertices_face, p, d) == 0) return 1;
    return 0;
}


int is_case_4(double **vertices_face, double *p, double *d)
{
    if (orientation(vertices_face, p, 3) == 0) {  // abcp live in the same plane
        puts("DEBUG IS CASE 4: abcp in same plane!");
        exit(4);
        double *aux[3];
        aux[0] = p; 

        aux[1] = vertices_face[0];
        aux[2] = vertices_face[1];
        if (orientation(aux, d, 3) == 0) return 1;

        aux[1] = vertices_face[1];
        aux[2] = vertices_face[2];
        if (orientation(aux, d, 3) == 0) return 1;

        aux[1] = vertices_face[0];
        aux[2] = vertices_face[2];
        if (orientation(aux, d, 3) == 0) return 1;
    }
    return 0;
}


int determine_case(double **vertices_face, double *p, double *d) 
{
    if (is_case_4(vertices_face, p, d)) return 4;
    else if (is_case_3(vertices_face, p, d)) return 3;
    else if (is_case_2(vertices_face, p, d)) return 2;
    else if (is_case_1(vertices_face, p, d)) return 1;
    puts("Should never reach this!");
    exit(1);
}


// ---------------------------------------------------------------------------------------
// --------------------------------- ALGORITHM -------------------------------------------
// ---------------------------------------------------------------------------------------

s_setup *initialize_setup(double **points, int N_points, int dim)
{
    // setup->points is EXTENDED FOR THE EXTRA NODES OF BIG_NCELL
    double **setup_points = malloc_matrix(N_points + dim + 1, dim);
    for (int ii=0; ii<N_points; ii++) {
        for (int jj=0; jj<dim; jj++) {
            setup_points[ii][jj] = points[ii][jj];
        }
    }
    
    double CM[dim]; /*  = malloc(sizeof(double) * dim); */
    find_center_mass(points, N_points, dim, CM);
    double maxd = max_distance(points, N_points, dim, CM);

    // Build the vertices of a regular simplex in R^(dim+1) with circumsphere of radius s centered at origin
    double s = 3 * maxd * dim * sqrt((dim+1.0)/dim);  // Scale so that inradius = 1.5 * maxd, original norm = sqrt(dim/(dim+1))
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
            double aux = 2.0 * rand() / RAND_MAX - 1;  // ADD SOME NOISE TO AVOID COLINEARITIES
            setup_points[N_points + ii][jj] = CM[jj] + V[ii][jj] + 0.001 * aux * s ; 
        }
    }
    free_matrix(V, dim+1);
    

    // CHECKS???
    // if (dim == 3) assert(are_in_general_position_3d(setup_points, N_points+dim+1) == 1 && "setup points are not in general position.");
    // THIS WAS WAY TOO SLOW !!
    if (N_points > 4) {
        for (int ii=0; ii<N_points; ii++) {  // DEBUG 
            assert(in_sphere(&setup_points[N_points], setup_points[ii], dim) == 1 && "big_ncell does not enclose all points.");
        }
    }


    s_setup *setup = malloc(sizeof(s_setup));
    setup->dim = dim;
    setup->N_points = N_points + dim + 1;
    setup->points = setup_points;

    s_ncell *big_ncell = malloc_ncell(setup);
    for (int ii=0; ii<setup->dim+1; ii++) {
        big_ncell->vertex_id[ii] = N_points + ii;
        big_ncell->opposite[ii] = NULL;
    }
    setup->head = big_ncell;
    setup->N_ncells = 1;
    
    return setup;
}


void flip_tetrahedra(s_setup *setup, s_stack *stack, s_ncell *ncell, int opp_cell_id)
{
    static double **coords_face = NULL;
    if (!coords_face) {
        coords_face = malloc_matrix(3, 3);
    }
    extract_vertices_face(setup, ncell, &opp_cell_id, 2, coords_face);

    // Extract id of vertex in opposite cell corresponding to the face
    int opp_face_localid;
    face_localid_of_adjacent_ncell(setup, ncell, &opp_cell_id, 2, opp_cell_id, &opp_face_localid);
    int opp_face_vertex_id = (ncell->opposite[opp_cell_id])->vertex_id[opp_face_localid];

    double *p = setup->points[ncell->vertex_id[opp_cell_id]];
    double *d = setup->points[opp_face_vertex_id];
    
    switch(determine_case(coords_face, p, d)) {
        int ridge_id_2;
        case 1:
            flip23(setup, stack, ncell, opp_cell_id, opp_face_localid, NULL);
            break;
        case 2:
            if (can_perform_flip32(setup, ncell, opp_cell_id, &ridge_id_2)) {
                flip32(setup, stack, ncell, opp_cell_id, ridge_id_2, opp_face_localid, NULL);
            }
            break;
        case 3:
            // printf("CASE 3: a: (%f, %f, %f)\nb: (%f, %f, %f)\nc: (%f, %f, %f)\nd: (%f, %f, %f)\np: (%f, %f, %f)\n", 
            //         coords_face[0][0], coords_face[0][1], coords_face[0][2],
            //         coords_face[1][0], coords_face[1][1], coords_face[1][2],
            //         coords_face[2][0], coords_face[2][1], coords_face[2][2],
            //         d[0], d[1], d[2], p[0], p[1], p[2]);
            if (can_perform_flip44(setup, ncell, coords_face, opp_cell_id, &ridge_id_2)) {
                // printf("executed! opp_cell_id = %d, ridge_id_2 = %d, ncell:", opp_cell_id, ridge_id_2);
                // print_ncell(setup, ncell);
                s_ncell *FLIP44_PTRS[4];
                flip44(setup, stack, ncell, opp_cell_id, ridge_id_2, FLIP44_PTRS);
            }
            break;
        case 4:
            puts("TODO CASE 4...");
            exit(1);
    }

    return;
}


void insert_one_point(s_setup *setup, int point_id, s_stack *stack)
{
    double *point = setup->points[point_id];
    s_ncell *container_ncell = in_ncell_walk(setup, point);

    // Insert p in container_ncell with a flip14
    flip14(setup, container_ncell, point_id, stack);

    while (stack->size > 0) {
        s_ncell *current = stack_pop(stack);

        if (current) {  // UNSURE IF THIS IS RIGHT... TODO
            int opp_cell_id = id_where_equal_int(current->vertex_id, 4, point_id);
            if (current->opposite[opp_cell_id]) {
                if (are_locally_delaunay_strict(setup, current, opp_cell_id) != 1) {
                    flip_tetrahedra(setup, stack, current, opp_cell_id);
                }
            }
        }
    }
}


void remove_big_tetra(s_setup *setup)
{
    s_ncell *current = setup->head;
    while (current) {
        s_ncell *next = current->next;
        for (int ii=0; ii<4; ii++) {
            if (current->vertex_id[ii] >= setup->N_points - 4) {  // This checks if a vertex is part of BIG TETRA
                if (current->next) (current->next)->prev = current->prev;
                if (current->prev) (current->prev)->next = next;
                else setup->head = current->next;

                // Update opposite's opposite to NULL
                for (int jj=0; jj<4; jj++) {
                    if (current->opposite[jj]) {
                        for (int kk=0; kk<4; kk++) {
                            if (current->opposite[jj]->opposite[kk] == current) {
                                current->opposite[jj]->opposite[kk] = NULL;
                                break;
                            }
                        }
                    }
                }

                free_ncell(current);
                setup->N_ncells--;
                break;
            }
        }
        current = next;
    }
    setup->points = realloc_matrix(setup->points, setup->N_points, setup->N_points-4, 3);
    setup->N_points -= 4;
}


int count_valid_ncells_reduced_triangulation(const s_setup *setup)
{
    int kk = 0;
    s_ncell *current = setup->head;
    while (current) {
        s_ncell *next = current->next;
        for (int ii=0; ii<4; ii++) {
            if (current->vertex_id[ii] >= setup->N_points - 4) {  // This checks if a vertex is part of BIG TETRA
                kk++;
                break;
            }
        }
        current = next;
    }
    return setup->N_ncells - kk;
}


s_setup *construct_dt_3d(double **points, int N_points)
{
    s_stack *stack = stack_create();
    s_setup *setup = initialize_setup(points, N_points, 3);
    
     assert(is_delaunay_3d(setup) == 1 && "Setup is not delaunay");  // DEBUG
    for (int ii=0; ii<N_points; ii++) {
        insert_one_point(setup, ii, stack);
        if (is_delaunay_3d(setup) != 1) {  // TODO, DEBUG
            printf("ERROR? NOT DELAUNAY AFTER INSERTING: %d\n", ii);
            exit(1);
        }
    }
    
    stack_free(stack);  
    remove_big_tetra(setup);

    // DEBUG: COMPUTE VOLUMES!
    s_ncell *current = setup->head;
    while (current) {
        add_ncell_volume_3d(setup, current);
        current = current->next;
    }

    return setup;
}
