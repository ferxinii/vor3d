
#include "delaunay.h"
#include "scplx.h"
#include "array.h"
#include "geometry.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

// ---------------------------------------------------------------------------------------
// ----------------------------------- STACK ---------------------------------------------
// ---------------------------------------------------------------------------------------

#define MAX_STACK_SIZE 10000

typedef struct dstack {
    s_ncell *entry[MAX_STACK_SIZE];
    int size;
} s_dstack;

static s_dstack *stack_create(void);
static void stack_free(s_dstack *stack);
static void stack_push(s_dstack *stack, s_ncell *ncell);
static s_ncell *stack_pop(s_dstack *stack);
static void stack_remove_ncell(s_dstack *stack, s_ncell *ncell);


static s_dstack *stack_create(void)
{
    s_dstack *stack = malloc(sizeof(s_dstack));
    stack->size = 0;
    return stack;
}


static void stack_free(s_dstack *stack)
{
    free(stack->entry);
}


static void stack_push(s_dstack *stack, s_ncell *ncell)
{
    assert(stack->size < MAX_STACK_SIZE && "Reached stack limit.");
    for (int ii=0; ii<stack->size; ii++) {
        if (stack->entry[ii] == ncell) return;
    }
    stack->entry[stack->size] = ncell;
    stack->size++;
}


static s_ncell *stack_pop(s_dstack *stack)
{   
    if (stack->size == 0) return NULL;
    stack->size--;
    s_ncell *ncell = stack->entry[stack->size];
    return ncell;
}


static void stack_remove_ncell(s_dstack *stack, s_ncell *ncell) {
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



// ---------------------------------------------------------------------------------------
// ----------------------------------- FLIPS ---------------------------------------------
// ---------------------------------------------------------------------------------------

static void flip14(s_scplx *setup, s_ncell *container_ncell, int point_id, s_dstack *stack);
static void flip23(s_scplx *setup, s_dstack *stack, s_ncell *nc1, int opp_cell_id, int opp_face_localid, s_ncell **OUT_PTRS);
static int can_perform_flip32(const s_scplx *setup, const s_ncell *ncell, int opp_cell_id, int *ridge_id_2);
static void flip32(s_scplx *setup, s_dstack *stack, s_dstack *stack_blocked, s_ncell *nc1, int opp_cell_id, int ridge_id_2, int opp_face_localid, s_ncell **OUT_PTRS);
static int can_perform_flip44(const s_scplx *setup, const s_ncell *ncell, s_point *vertices_face, int opp_cell_id, int *ridge_id_2);
static void flip44(s_scplx *setup, s_dstack *stack, s_dstack *stack_blocked, s_ncell *ncell, int id_ridge_1, int id_ridge_2, s_ncell **OUT_PTRS);


static void flip14(s_scplx *setup, s_ncell *container_ncell, int point_id, s_dstack *stack)
{   
    setup->N_ncells += 3;

    s_ncell *nc2 = malloc_ncell();
    s_ncell *nc3 = malloc_ncell();
    s_ncell *nc4 = malloc_ncell();
    
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
        face_localid_of_adjacent_ncell(nc2, 2, &aux, 1, &opp_id); 
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
        face_localid_of_adjacent_ncell(nc3, 2, &aux, 2, &opp_id); 
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
        face_localid_of_adjacent_ncell(nc4, 2, &aux, 3, &opp_id); 
        (nc4->opposite[3])->opposite[opp_id] = nc4;
    }

    stack_push(stack, container_ncell);
    stack_push(stack, nc2);
    stack_push(stack, nc3);
    stack_push(stack, nc4);
}


static void flip23(s_scplx *setup, s_dstack *stack, s_ncell *nc1, int opp_cell_id, int opp_face_localid, s_ncell **OUT_PTRS)
{   
    setup->N_ncells += 1;

    s_ncell *nc2 = nc1->opposite[opp_cell_id];
    s_ncell *nc3 = malloc_ncell();

    // Update linked-list of ncells
    // ---- NC1 ----- NC2 ---------------- 
    // ---- NC1 ----- NC2 ----- NC3 ------
    nc3->next = nc2->next;
    nc3->prev = nc2;
    if (nc3->next) nc3->next->prev = nc3;

    nc2->next = nc3;

    int face_vertex_id[3];
    extract_ids_face(nc1, 2, &opp_cell_id, face_vertex_id);

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
        face_localid_of_adjacent_ncell(nc1, 2, &aux2, aux2, &v_localid_opp);
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
        face_localid_of_adjacent_ncell(nc2, 2, &aux2, aux2, &v_localid_opp);
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
        face_localid_of_adjacent_ncell(nc3, 2, &aux2, aux2, &v_localid_opp);
        nc_aux->opposite[v_localid_opp] = nc3;
    }

    aux = id_where_equal_int(nc2_vertex_id_old, 4, b);
    nc_aux = opp_2_old[aux];
    if (nc_aux) {
        aux2 = id_where_equal_int(nc3->vertex_id, 4, p);
        face_localid_of_adjacent_ncell(nc3, 2, &aux2, aux2, &v_localid_opp);
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


static int can_perform_flip32(const s_scplx *setup, const s_ncell *ncell, int opp_cell_id, int *ridge_id_2)
{
    if (setup->N_ncells < 3) return 0;

    int face_vertex_id[3];
    extract_ids_face(ncell, 2, &opp_cell_id, face_vertex_id);
    
    int opp_face_localid;
    face_localid_of_adjacent_ncell(ncell, 2, &opp_cell_id, opp_cell_id, &opp_face_localid);
    
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


static void flip32(s_scplx *setup, s_dstack *stack, s_dstack *stack_blocked, s_ncell *nc1, int opp_cell_id, int ridge_id_2, int opp_face_localid, s_ncell **OUT_PTRS)
{
    setup->N_ncells -= 1;

    int v2_main, v2_2, v3_main, v3_2;
    s_ncell *nc2 = next_ncell_ridge_cycle(nc1, opp_cell_id, ridge_id_2, &v2_main, &v2_2);
    s_ncell *nc3 = next_ncell_ridge_cycle(nc2, v2_main, v2_2, &v3_main, &v3_2);

    // ------ NC3->PREV ----- NC3 ------- NC3->NEXT ------
    s_ncell *nc3_next = nc3->next;
    if (nc3->next) nc3->next->prev = nc3->prev;
    if (nc3->prev) nc3->prev->next = nc3_next;
    else setup->head = nc3->next;


    int localids_ridge[2] = {opp_cell_id, ridge_id_2};
    int vertices_ridge[2];
    extract_ids_face(nc1, 1, localids_ridge, vertices_ridge);

    int p = nc1->vertex_id[opp_cell_id];
    int a = nc1->vertex_id[ridge_id_2];
    int b = vertices_ridge[0];
    int c = vertices_ridge[1];
    int d = nc2->vertex_id[opp_face_localid];

    // WE NEED TO MAKE SURE THEY ARE ORIENTED ??? TODO NOT SURE ABOUT THIS, IS IT NECESSARY?
    // s_point face_vertices[3];
    // face_vertices[0] = setup->points[a];
    // face_vertices[1] = setup->points[b];
    // face_vertices[2] = setup->points[c];
    // if ( orientation(face_vertices, setup->points[nc1->vertex_id[opp_cell_id]], 3) == -1 ) {
    //     int temp = b;
    //     b = c;
    //     c = temp;
    // }
    // face_vertices[0] = setup->points[a];
    // face_vertices[1] = setup->points[b];
    // face_vertices[2] = setup->points[c];
    // assert(orientation(face_vertices, setup->points[nc1->vertex_id[opp_cell_id]], 3) == 1);

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
        face_localid_of_adjacent_ncell(nc1, 2, &aux2, aux2, &v_localid_opp);
        nc_aux->opposite[v_localid_opp] = nc1;
    }

    aux = id_where_equal_int(nc3_vertex_id_old, 4, c);
    nc_aux = opp_3_old[aux];
    if (nc_aux) {
        aux2 = id_where_equal_int(nc1->vertex_id, 4, a);
        face_localid_of_adjacent_ncell(nc1, 2, &aux2, aux2, &v_localid_opp);
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
        face_localid_of_adjacent_ncell(nc2, 2, &aux2, aux2, &v_localid_opp);
        nc_aux->opposite[v_localid_opp] = nc2;
    }

    aux = id_where_equal_int(nc3_vertex_id_old, 4, b);
    nc_aux = opp_3_old[aux];
    if (nc_aux) {
        aux2 = id_where_equal_int(nc2->vertex_id, 4, a);
        face_localid_of_adjacent_ncell(nc2, 2, &aux2, aux2, &v_localid_opp);
        nc_aux->opposite[v_localid_opp] = nc2;
    }
    
    stack_remove_ncell(stack, nc3);
    if (stack_blocked) {
        stack_remove_ncell(stack_blocked, nc3);
    }

    free_ncell(nc3);

    stack_push(stack, nc1);
    stack_push(stack, nc2);
    if (OUT_PTRS) {
        OUT_PTRS[0] = nc1;
        OUT_PTRS[1] = nc2;
    }
}


static int can_perform_flip44(const s_scplx *setup, const s_ncell *ncell, s_point *vertices_face, int opp_cell_id, int *ridge_id_2)
{
    if (setup->N_ncells < 4) return 0;
    
    s_ncell *opp = ncell->opposite[opp_cell_id];
    int opp_face_localid;
    face_localid_of_adjacent_ncell(ncell, 2, &opp_cell_id, opp_cell_id, &opp_face_localid);


    // DETERMINE RIDGE 
    int face_vertex_id[4];
    extract_ids_face(ncell, 2, &opp_cell_id, face_vertex_id);

    s_point aux[3];
    aux[0] = setup->points.p[ncell->vertex_id[opp_cell_id]]; 

    aux[1] = vertices_face[0];
    aux[2] = vertices_face[1];
    if (orientation(aux, setup->points.p[opp->vertex_id[opp_face_localid]]) == 0) {
        *ridge_id_2 = id_where_equal_int(ncell->vertex_id, 4, face_vertex_id[2]);
    }
    aux[1] = vertices_face[1];
    aux[2] = vertices_face[2];
    if (orientation(aux, setup->points.p[opp->vertex_id[opp_face_localid]]) == 0) {
        *ridge_id_2 = id_where_equal_int(ncell->vertex_id, 4, face_vertex_id[0]);
    }
    aux[1] = vertices_face[0];
    aux[2] = vertices_face[2];
    if (orientation(aux, setup->points.p[opp->vertex_id[opp_face_localid]]) == 0) {
        *ridge_id_2 = id_where_equal_int(ncell->vertex_id, 4, face_vertex_id[1]);
    }
    

    int nc2_id1, nc2_id2;
    s_ncell *nc2 = next_ncell_ridge_cycle(ncell, opp_cell_id, *ridge_id_2, &nc2_id1, &nc2_id2);

    int nc3_id1, nc3_id2;
    s_ncell *nc3 = next_ncell_ridge_cycle(nc2, nc2_id1, nc2_id2, &nc3_id1, &nc3_id2);

    int nc4_id1, nc4_id2;
    s_ncell *nc4 = next_ncell_ridge_cycle(nc3, nc3_id1, nc3_id2, &nc4_id1, &nc4_id2);

    if (nc4->opposite[nc4_id1] == ncell) return 1;
    else return 0;
}


static void flip44(s_scplx *setup, s_dstack *stack, s_dstack *stack_blocked, s_ncell *ncell, int id_ridge_1, int id_ridge_2, s_ncell **OUT_PTRS) 
{

    // FIRST, A FLIP23
    int opp_face_localid;
    face_localid_of_adjacent_ncell(ncell, 2, &id_ridge_1, id_ridge_1, &opp_face_localid);
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
    face_localid_of_adjacent_ncell(nc5, 2, &nc5_p, nc5_p, &nc3_id2);
    int nc3_opp_face_localid;
    face_localid_of_adjacent_ncell(nc3, 2, &nc3_id1, nc3_id1, &nc3_opp_face_localid);

    s_ncell *FLIP32_PTRS[2];
    flip32(setup, stack, stack_blocked, nc3, nc3_id1, nc3_id2, nc3_opp_face_localid, FLIP32_PTRS);

    if (OUT_PTRS) {
        OUT_PTRS[2] = FLIP32_PTRS[0];
        OUT_PTRS[3] = FLIP32_PTRS[1];
    }
    stack_push(stack, FLIP32_PTRS[0]);
    stack_push(stack, FLIP32_PTRS[1]);
}



// ---------------------------------------------------------------------------------------
// ----------------------------------- CASES ---------------------------------------------
// ---------------------------------------------------------------------------------------

static int pd_intersects_edge(s_point s1, s_point s2, s_point p, s_point d, int drop_coord);
static int is_case_3(const s_point vertices_face[3], s_point p, s_point d);
static int is_case_1(const s_point vertices_face[3], s_point p, s_point d);
static int is_case_2(const s_point vertices_face[3], s_point p, s_point d);
static int is_case_4(const s_point vertices_face[3], s_point p, s_point d);
static int determine_case(const s_point vertices_face[3], s_point p, s_point d);


static int pd_intersects_edge(s_point s1, s_point s2, s_point p, s_point d, int drop_coord)
{
    int i1, i2;
    if (drop_coord == 0)      { i1 = 1; i2 = 2; } 
    else if (drop_coord == 1) { i1 = 2; i2 = 0; } 
    else                      { i1 = 0; i2 = 1; }

    s_point A, B, paux, daux;
    A.x = s1.coords[i1];    A.y = s1.coords[i2];
    B.x = s2.coords[i1];    B.y = s2.coords[i2];
    paux.x = p.coords[i1];  paux.x = p.coords[i2];
    daux.y = d.coords[i1];  daux.y = d.coords[i2];
    
    s_point AB[2] = {A, B},     pd[2] = {paux, daux};
    return segments_intersect_2d(AB, pd);
}


static int is_case_3(const s_point vertices_face[3], s_point p, s_point d)
{
    s_point d1 = subtract_points(vertices_face[1], vertices_face[0]);
    s_point d2 = subtract_points(vertices_face[2], vertices_face[0]);
    s_point n = cross_prod(d1, d2);
    int drop_coord = coord_with_largest_component_3d(n);

    s_point aux[3];
    aux[0] = p; 

    aux[1] = vertices_face[0];
    aux[2] = vertices_face[1];
    if (orientation(aux, d) == 0 &&
        pd_intersects_edge(vertices_face[0], vertices_face[1], p, d, drop_coord)) return 1;

    aux[1] = vertices_face[1];
    aux[2] = vertices_face[2];
    if (orientation(aux, d) == 0 &&
        pd_intersects_edge(vertices_face[1], vertices_face[2], p, d, drop_coord)) return 1;

    aux[1] = vertices_face[0];
    aux[2] = vertices_face[2];
    if (orientation(aux, d) == 0 &&
        pd_intersects_edge(vertices_face[0], vertices_face[2], p, d, drop_coord)) return 1;
    
    return 0;
}


static int is_case_1(const s_point vertices_face[3], s_point p, s_point d) 
{
    if (segment_crosses_triangle_3d(vertices_face, p, d) == 1) return 1;
    return 0;
}


static int is_case_2(const s_point vertices_face[3], s_point p, s_point d)
{
    if (segment_crosses_triangle_3d(vertices_face, p, d) == 0) return 1;
    return 0;
}


static int is_case_4(const s_point vertices_face[3], s_point p, s_point d)
{
    if (orientation(vertices_face, p) == 0) {  // abcp live in the same plane
        // puts("DEBUG ISCASE4: ABCP IN SAME PLANE!");
        s_point aux[3];
        aux[0] = p; 

        aux[1] = vertices_face[0];
        aux[2] = vertices_face[1];
        if (orientation(aux, d) == 0) return 1;

        aux[1] = vertices_face[1];
        aux[2] = vertices_face[2];
        if (orientation(aux, d) == 0) return 1;

        aux[1] = vertices_face[0];
        aux[2] = vertices_face[2];
        if (orientation(aux, d) == 0) return 1;
    }
    // puts("DEBUG ISCASE4: BUT NOT IN FACE...");
    return 0;
}


static int determine_case(const s_point vertices_face[3], s_point p, s_point d) 
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

static s_scplx initialize_setup(const s_points *points);
static int flip_tetrahedra(s_scplx *setup, s_dstack *stack, s_dstack *stack_blocked, s_ncell *ncell, int opp_cell_id);
static void remove_point_setup(s_scplx *setup, int point_id);
static int insert_one_point(s_scplx *setup, int point_id, s_dstack *stack, s_dstack *stack_blocked);
static void remove_big_tetra(s_scplx *setup);


static s_scplx initialize_setup(const s_points *points)
{
    // setup->points is EXTENDED FOR THE EXTRA NODES OF BIG_NCELL, PUT AT THE BEGINNING!
    s_points setup_points = { .N = points->N + 4, 
                              .p = malloc(sizeof(s_point) * (points->N + 4)) };
    for (int ii=0; ii<points->N; ii++) {
        setup_points.p[ii+4] = points->p[ii];
    }
    
    s_point CM = point_average(points);
    double maxd = max_distance(points, CM);

    // Build the vertices of a regular simplex in R^(dim+1) with circumsphere of radius s centered at origin
    int dim = 3;
    double s = 3 * maxd * dim * sqrt((dim+1.0)/dim);
    int nvert = dim + 1;
    double shift = 1.0 / (dim + 1);

    for (int ii = 0; ii < nvert; ++ii) {
        // compute only the first 3 coordinates (we drop the last one)
        double v0 = ((ii == 0) ? 1.0 : 0.0) - shift;
        double v1 = ((ii == 1) ? 1.0 : 0.0) - shift;
        double v2 = ((ii == 2) ? 1.0 : 0.0) - shift;
        // note: v3 = ((ii==3)?1.0:0.0) - shift exists but we drop it

        double jitter = 0.001 * s;
        setup_points.p[ii].x = CM.x + v0 * s + jitter;
        setup_points.p[ii].y = CM.y + v1 * s + jitter;
        setup_points.p[ii].z = CM.z + v2 * s + jitter;
    }

    
    s_scplx setup;
    setup.points.N = points->N + dim + 1;
    setup.points = setup_points;

    s_ncell *big_ncell = malloc_ncell();
    for (int ii=0; ii<4; ii++) {
        big_ncell->vertex_id[ii] = ii;
        big_ncell->opposite[ii] = NULL;
    }
    setup.head = big_ncell;
    setup.N_ncells = 1;
    
    return setup;
}


static int flip_tetrahedra(s_scplx *setup, s_dstack *stack, s_dstack *stack_blocked, s_ncell *ncell, int opp_cell_id)
{
    s_point coords_face[3];
    extract_vertices_face(setup, ncell, 2, &opp_cell_id, coords_face);

    // Extract id of vertex in opposite cell corresponding to the face
    int opp_face_localid;
    face_localid_of_adjacent_ncell(ncell, 2, &opp_cell_id, opp_cell_id, &opp_face_localid);
    int opp_face_vertex_id = (ncell->opposite[opp_cell_id])->vertex_id[opp_face_localid];

    s_point p = setup->points.p[ncell->vertex_id[opp_cell_id]];
    s_point d = setup->points.p[opp_face_vertex_id];
    
    switch(determine_case(coords_face, p, d)) {
        int ridge_id_2;
        case 1:
            flip23(setup, stack, ncell, opp_cell_id, opp_face_localid, NULL);
            return 1;
        case 2:
            if (can_perform_flip32(setup, ncell, opp_cell_id, &ridge_id_2)) {
                flip32(setup, stack, stack_blocked, ncell, opp_cell_id, ridge_id_2, opp_face_localid, NULL);
                return 1;
            } else {if (stack_blocked) stack_push(stack_blocked, ncell);}
            break;
        case 3:
            if (can_perform_flip44(setup, ncell, coords_face, opp_cell_id, &ridge_id_2)) {
                s_ncell *FLIP44_PTRS[4];
                flip44(setup, stack, stack_blocked, ncell, opp_cell_id, ridge_id_2, FLIP44_PTRS);
                return 1;
            } else {if (stack_blocked) stack_push(stack_blocked, ncell);}
            break;
        case 4:
            puts("CASE 4... UNSURE, UNTESTED");
            flip23(setup, stack, ncell, opp_cell_id, opp_face_localid, NULL);
            return 1;
    }

    return 0;
}


static void remove_point_setup(s_scplx *setup, int point_id)
{
    if (point_id < setup->points.N-1) {
        for (int ii=point_id; ii<setup->points.N-1; ii++) {
            setup->points.p[ii] = setup->points.p[ii+1];
        }
    }
    setup->points.p = realloc(setup->points.p, sizeof(s_point) * setup->points.N-1);
    setup->points.N--;
}


static int insert_one_point(s_scplx *setup, int point_id, s_dstack *stack, s_dstack *stack_blocked)
{
    s_point point = setup->points.p[point_id];
    s_ncell *container_ncell = in_ncell_walk(setup, point);

    double EPS2 = 1e-12;
    if (distance_squared(setup->points.p[container_ncell->vertex_id[0]], point) < EPS2 || 
        distance_squared(setup->points.p[container_ncell->vertex_id[1]], point) < EPS2 ||
        distance_squared(setup->points.p[container_ncell->vertex_id[2]], point) < EPS2 ||
        distance_squared(setup->points.p[container_ncell->vertex_id[3]], point) < EPS2) {
        // puts("insert_one_point: POINT EXISTS!");
        remove_point_setup(setup, point_id);
        return 0;
    }

    // Insert p in container_ncell with a flip14
    flip14(setup, container_ncell, point_id, stack);

    stack_blocked->size = 0;
    while (stack->size > 0) {
        s_ncell *current = stack_pop(stack);

        if (current) {  // UNSURE IF THIS IS RIGHT... TODO
            int opp_cell_id = id_where_equal_int(current->vertex_id, 4, point_id);
            if (current->opposite[opp_cell_id]) {
                if (!are_locally_delaunay_strict(setup, current, opp_cell_id)) {
                    flip_tetrahedra(setup, stack, stack_blocked, current, opp_cell_id);
                }
            }
        }
    }
    return 1;
}


static void remove_big_tetra(s_scplx *setup)
{
    s_ncell *current = setup->head;
    while (current) {
        s_ncell *next = current->next;
        for (int ii=0; ii<4; ii++) {
            if (current->vertex_id[ii] < 4) { // This checks if a vertex is part of BIG TETRA
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
    for (int ii=0; ii<setup->points.N-4; ii++) {
        setup->points.p[ii] = setup->points.p[ii+4];
    }
    setup->points.p = realloc(setup->points.p, sizeof(s_point) * (setup->points.N-4));
    setup->points.N -= 4;

    // Reindex all remaining tetrahedra so their vertex_id points at correct coords
    for (s_ncell *c=setup->head; c; c=c->next) {
        for (int kk=0; kk<4; kk++) {
            c->vertex_id[kk] -= 4;
        }
    }
}




// ---------------------------------------------------------------------------------------
// -------------------------------------- MAIN -------------------------------------------
// ---------------------------------------------------------------------------------------

s_scplx construct_dt_3d(const s_points *points)
{
    s_dstack *stack = stack_create();
    s_dstack *stack_blocked = stack_create();
    s_scplx setup = initialize_setup(points);
    
    int ii = 4;  // First 4 are big tetra, which already is inserted!
    while (ii < setup.points.N) {
        if (insert_one_point(&setup, ii, stack, stack_blocked)) {
            // SHOULD ALWAYS BE DELAUNAY! (only locally?)
            ii++;
        } else continue; //printf("Not inserted! ii=%d\n", ii);
        // if (!is_delaunay_3d(&setup)) printf("Not delaunay!! ii=%d\n", ii);
    }
    
    stack_free(stack);  
    remove_big_tetra(&setup);

    return setup;
}

