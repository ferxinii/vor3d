/* Algorithm from article: 
 * Hugo Ledoux. Computing the 3D Voronoi Diagram Robustly: An Easy Explanation. 
 * But considering the possibility of weights (regular triangulation) */
#include "delaunay.h"
#include "scplx.h"
#include "points.h"
#include "gtests.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>


bool are_locally_delaunay(const s_scplx *scplx, const s_ncell *ncell, int id_opposite,
                          e_delaunay_test_type type)
{   
    s_point coords1[4], coords2[4];
    extract_vertices_ncell(scplx, ncell, coords1);
    extract_vertices_ncell(scplx, ncell->opposite[id_opposite], coords2);

    if (test_orientation(coords1, coords1[3]) == 0 || 
        test_orientation(coords2, coords2[3]) == 0) {
        return false;  /* If any is flat, return 0 */
    }

    int opp_face_localid;
    face_localid_of_adjacent_ncell(ncell, 2, &id_opposite, id_opposite, &opp_face_localid);
    int opp_face_vertex_id = (ncell->opposite[id_opposite])->vertex_id[opp_face_localid];
    
    int in1;
    if (scplx->weights) {
        double weights1[4]; extract_weights_ncell(scplx, ncell, weights1);
        in1 = test_orthosphere(4, coords1, weights1, 
                scplx->points.p[opp_face_vertex_id], scplx->weights[opp_face_vertex_id]);
    } else {
        in1 = test_insphere(coords1, scplx->points.p[opp_face_vertex_id]);
    }

    switch (type) {
        case DELAUNAY_TEST_STRICT:
            if (in1 == -1) return true;
            else return false;
        case DELAUNAY_TEST_NONSTRICT:
            if (in1 != 1) return true;
            else return false;
    }
}

bool is_delaunay_3d(const s_scplx *scplx, e_delaunay_test_type type)
{
    s_ncell *current = scplx->head;
    while (current) {
        s_point vertices_ncell[4];
        extract_vertices_ncell(scplx, current, vertices_ncell);
        if (test_orientation(vertices_ncell, vertices_ncell[3]) == 0) {
            fprintf(stderr, "Flat tetra!! vids: %d %d %d %d\n",
                current->vertex_id[0], current->vertex_id[1],
                current->vertex_id[2], current->vertex_id[3]);
            return false;
        }
        for (int ii=0; ii<4; ii++) {
            if (current->opposite[ii] &&
                !are_locally_delaunay(scplx, current, ii, type)) {
                    fprintf(stderr, "Insphere failure: ncell vids %d %d %d %d, opposite vids %d %d %d %d\n",
                        current->vertex_id[0], current->vertex_id[1],
                        current->vertex_id[2], current->vertex_id[3],
                        current->opposite[ii]->vertex_id[0],
                        current->opposite[ii]->vertex_id[1],
                        current->opposite[ii]->vertex_id[2],
                        current->opposite[ii]->vertex_id[3]);
                return false;
            }
        }
        current = current->next;
    }
    return true;
}


static bool p_locally_redundant_in_ncell(const s_scplx *scplx, const s_ncell *nc, int p_id)
{   /* In unweighted triangulation, all points are NON-redundant */
    if (!scplx->weights) return false;

    /* Never redundant in a tetrahedron containing sentinel vertices */
    for (int i = 0; i < 4; i++)
        if (nc->vertex_id[i] < 4) return false;

    s_point v[4]; double w[4];
    extract_vertices_and_weights_ncell(scplx, nc, v, w);
    return (test_orthosphere(4, v, w, scplx->points.p[p_id], scplx->weights[p_id]) != -1);



    // s_point v[4];  double w[4];
    // extract_vertices_and_weights_ncell(scplx, nc, v, w);
    // return (test_orthosphere(4, v, w, scplx->points.p[p_id], scplx->weights[p_id]) == 1);
    //
    //
    // /* Count how many and which vertices of nc are part of big tetra */
    // int N_bigtetra = 0, bigtetra_ids[4] = {0};
    // for (int i=0; i<4; i++) 
    //     if (nc->vertex_id[i] < 4) bigtetra_ids[N_bigtetra++] = i;
    //
    // /* Select the vertices that are NOT part of big tetra */
    // s_point v_real[4];  double w_real[4];
    // for (int i = 0, k = 0; i < 4; i++) {
    //     bool is_bigtetra = false;
    //     for (int j = 0; j < N_bigtetra; j++)
    //         if (bigtetra_ids[j] == i) { is_bigtetra = true; break; }
    //     if (!is_bigtetra) {
    //         v_real[k] = v[i];
    //         w_real[k] = w[i];
    //         k++;
    //     }
    // }
    //
    // switch (N_bigtetra) {
    //     case 0:
    //         return (test_orthosphere(4, v, w, scplx->points.p[p_id], scplx->weights[p_id]) == 1);
    //     case 1:
    //         // return (test_orthosphere(3, v_real, w_real, scplx->points.p[p_id], scplx->weights[p_id]) == 1);
    //         return false;
    //     case 2:
    //         return (test_orthosphere(2, v_real, w_real, scplx->points.p[p_id], scplx->weights[p_id]) == 1);
    //     case 3:
    //         return (test_orthosphere(1, v_real, w_real, scplx->points.p[p_id], scplx->weights[p_id]) == 1);
    //     case 4:
    //         return false;
    //     default:
    //         fprintf(stderr, "p_locally_redundant_in_ncell: error in switch!\n");
    //         exit(1);
    // }
}



// ----------------------------------- HELPERS ---------------------------------------------
//
static int id_where_equal_int(const int *arr, int N, int entry) 
{
    for (int ii=0; ii<N; ii++) if (arr[ii] == entry) return ii;
    fprintf(stderr, "id_where_equal_int: Could not find id.\n"); 
    assert(1==0);
    return -10;
}

static int inarray(const int *arr1, int N, int a)
{
    for (int ii=0; ii<N; ii++) if (arr1[ii] == a) return 1;
    return 0;
}

// ----------------------------------- STACK ---------------------------------------------

#define INIT_STACK_SIZE 100

typedef struct dstack {
    s_ncell **entry;
    int size;
    int capacity;
} s_dstack;

static s_dstack stack_create(void)
{
    s_dstack stack;
    stack.size = 0;
    stack.capacity = INIT_STACK_SIZE;
    stack.entry = malloc(stack.capacity * sizeof(s_ncell *));
    return stack;
}

static void stack_free(s_dstack *stack)
{
    free(stack->entry);
}

static int stack_push(s_dstack *stack, s_ncell *ncell)
{
    if (ncell->in_stack) return 1;  /* Avoid duplicates */

    if (stack->size == stack->capacity) { /* Expand if needed */
        stack->capacity *= 2;
        s_ncell **tmp = realloc(stack->entry, stack->capacity * sizeof(s_ncell *));
        if (!tmp) { fprintf(stderr, "delaunay.c, stack_push: realloc failed increasing stack\n"); return 0; }
        stack->entry = tmp;
    }

    stack->entry[stack->size++] = ncell;
    ncell->in_stack = true;
    return 1;
}

static s_ncell *stack_pop(s_dstack *stack)
{   
    while (stack->size > 0) {
        s_ncell *ncell = stack->entry[--stack->size];
        if (ncell->in_stack) {
            ncell->in_stack = false;
            return ncell;
        }
    }
    return NULL;
}

static void stack_remove_ncell(s_dstack *stack, s_ncell *ncell) {
    int newSize = 0;
    for (int ii = 0; ii < stack->size; ii++)
        if (stack->entry[ii] != ncell)  /* Only copy entries that are not the target. */
            stack->entry[newSize++] = stack->entry[ii];
    stack->size = newSize;
    ncell->in_stack = false;
}



// ----------------------------------- FLIPS ---------------------------------------------
static int flip_tetrahedra(s_scplx *scplx, s_dstack *stack, s_ncell *ncell, int opp_cell_id, bool *ignored);

static inline void vertex_ids_ncell(s_ncell *ncell, int out[4])
{
    for (int ii=0; ii<4; ii++) 
        out[ii] = ncell->vertex_id[ii];
}

static inline void opposite_pointers_ncell(s_ncell *ncell, s_ncell *out[4])
{
    for (int ii=0; ii<4; ii++)
        out[ii] = ncell->opposite[ii];
}

static inline void set_ncell_vids(s_ncell *ncell, int v1, int v2, int v3, int v4)
{
    ncell->vertex_id[0] = v1;   ncell->vertex_id[1] = v2;
    ncell->vertex_id[2] = v3;   ncell->vertex_id[3] = v4;
}

static inline void set_ncell_opposite_pointers(s_ncell *ncell, s_ncell *o1, s_ncell *o2, s_ncell *o3, s_ncell *o4)
{
    ncell->opposite[0] = o1;    ncell->opposite[1] = o2;
    ncell->opposite[2] = o3;    ncell->opposite[3] = o4;
}

static int flip14(s_scplx *scplx, s_ncell *nc1, int point_id, s_dstack *stack)
{   
    scplx->N_ncells += 3;
    s_ncell *nc2 = malloc_ncell(), *nc3 = malloc_ncell(), *nc4 = malloc_ncell();
    if (!nc2 || !nc3 || !nc4) return 0;

    s_ncell *opp_aux[4]; opposite_pointers_ncell(nc1, opp_aux);
    int v_aux[4]; vertex_ids_ncell(nc1, v_aux);
    
    /* Update linked-list of ncells */
    /* ---- NC1 ---------------------------- */
    /* ---- NC1 --- NC2 --- NC3 --- NC4 ---- */
    nc4->next = nc1->next; if (nc4->next) (nc4->next)->prev = nc4;  /* Last migh be null! */
    nc3->next = nc4; nc4->prev = nc3;
    nc2->next = nc3; nc3->prev = nc2;
    nc1->next = nc2; nc2->prev = nc1;

    /* Update nc1 */
    nc1->vertex_id[0] = point_id;
    nc1->opposite[1] = nc2;  nc1->opposite[2] = nc3;  nc1->opposite[3] = nc4;

    /* scplx nc2 */
    set_ncell_vids(nc2,    v_aux[0], point_id, v_aux[2], v_aux[3]);
    set_ncell_opposite_pointers(nc2,    nc1, opp_aux[1], nc3, nc4);
    if (nc2->opposite[1]) {
        int opp_id; face_localid_of_adjacent_ncell(nc2, 2, &(int){1}, 1, &opp_id); 
        nc2->opposite[1]->opposite[opp_id] = nc2;
    }

    /* scplx nc3 */
    set_ncell_vids(nc3,    v_aux[0], v_aux[1], point_id, v_aux[3]);
    set_ncell_opposite_pointers(nc3,    nc1, nc2, opp_aux[2], nc4);
    if (nc3->opposite[2]) {
        int opp_id; face_localid_of_adjacent_ncell(nc3, 2, &(int){2}, 2, &opp_id);
        nc3->opposite[2]->opposite[opp_id] = nc3;
    }

    /* scplx nc4 */
    set_ncell_vids(nc4,    v_aux[0], v_aux[1], v_aux[2], point_id);
    set_ncell_opposite_pointers(nc4,    nc1, nc2, nc3, opp_aux[3]);
    if (nc4->opposite[3]) {
        int opp_id; face_localid_of_adjacent_ncell(nc4, 2, &(int){3}, 3, &opp_id);
        nc4->opposite[3]->opposite[opp_id] = nc4;
    }

    /* Push to stack */
    if (stack) { 
        if (!stack_push(stack, nc1)) return 0;
        if (!stack_push(stack, nc2)) return 0;
        if (!stack_push(stack, nc3)) return 0;
        if (!stack_push(stack, nc4)) return 0;
    }
    return 1;
}


static inline void map_vid_lid(s_ncell *ncell, int v1, int *l1, int v2, int *l2, int v3, int *l3, int v4, int *l4)
{   /* Map vertex ids to localids */
    *l1 = id_where_equal_int(ncell->vertex_id, 4, v1);
    *l2 = id_where_equal_int(ncell->vertex_id, 4, v2);
    *l3 = id_where_equal_int(ncell->vertex_id, 4, v3);
    *l4 = id_where_equal_int(ncell->vertex_id, 4, v4);
}   

static int flip23(s_scplx *scplx, s_dstack *stack, s_ncell *nc1, int opp_cell_id, int opp_face_localid, s_ncell *OUT_PTRS[3])
{   
    scplx->N_ncells += 1;
    s_ncell *nc2 = nc1->opposite[opp_cell_id], *nc3 = malloc_ncell();
    if (!nc2) return 0;

    s_ncell *nc1_opp_old[4]; opposite_pointers_ncell(nc1, nc1_opp_old);
    s_ncell *nc2_opp_old[4]; opposite_pointers_ncell(nc2, nc2_opp_old);

    /* Update linked-list of ncells */
    /* ---- NC1 ----- NC2 ---------------- */
    /* ---- NC1 ----- NC2 ----- NC3 ------ */
    nc3->next = nc2->next; if (nc3->next) nc3->next->prev = nc3;
    nc2->next = nc3; nc3->prev = nc2;

    /* scplx important indices */
    int face_vid[3]; extract_ids_face(nc1, 2, &opp_cell_id, face_vid);
    int a = face_vid[0];
    int b = face_vid[1];
    int c = face_vid[2];
    int d = nc2->vertex_id[opp_face_localid];
    int p = nc1->vertex_id[opp_cell_id];

    int nc1_id_a, nc1_id_b, nc1_id_c, nc1_id_p;
    map_vid_lid(nc1, a, &nc1_id_a, b, &nc1_id_b, c, &nc1_id_c, p, &nc1_id_p);
    int nc2_id_a, nc2_id_b, nc2_id_c, nc2_id_d;
    map_vid_lid(nc2, a, &nc2_id_a, b, &nc2_id_b, c, &nc2_id_c, d, &nc2_id_d);

     /* Update nc1 */
    nc1->vertex_id[nc1_id_c] = d;
    nc1->opposite[nc1_id_a] = nc2;
    nc1->opposite[nc1_id_b] = nc3;
    nc1->opposite[nc1_id_p] = nc2_opp_old[nc2_id_c];
    if (nc2_opp_old[nc2_id_c]) {
        int opp_aux; face_localid_of_adjacent_ncell(nc1, 2, &nc1_id_p, nc1_id_p, &opp_aux);
        nc2_opp_old[nc2_id_c]->opposite[opp_aux] = nc1;
    }

    /* Update nc2 */
    nc2->vertex_id[nc2_id_a] = p;
    nc2->opposite[nc2_id_b] = nc3;
    nc2->opposite[nc2_id_c] = nc1;
    nc2->opposite[nc2_id_d] = nc1_opp_old[nc1_id_a];
    if (nc1_opp_old[nc1_id_a]) {
        int opp_aux; face_localid_of_adjacent_ncell(nc2, 2, &nc2_id_d, nc2_id_d, &opp_aux);
        nc1_opp_old[nc1_id_a]->opposite[opp_aux] = nc2;
    }

    /* Build nc3 */
    set_ncell_vids(nc3, p, c, d, a);
    set_ncell_opposite_pointers(nc3,  nc2_opp_old[nc2_id_b], nc1, nc1_opp_old[nc1_id_b], nc2);
    if (nc1_opp_old[nc1_id_b]) {
        int opp_aux; face_localid_of_adjacent_ncell(nc3, 2, &(int){2}, 2, &opp_aux);
        nc1_opp_old[nc1_id_b]->opposite[opp_aux] = nc3;
    }
    if (nc2_opp_old[nc2_id_b]) {
        int opp_aux; face_localid_of_adjacent_ncell(nc3, 2, &(int){0}, 0, &opp_aux);
        nc2_opp_old[nc2_id_b]->opposite[opp_aux] = nc3;
    }

    if (stack) { 
        if (!stack_push(stack, nc1)) return 0;
        if (!stack_push(stack, nc2)) return 0;
        if (!stack_push(stack, nc3)) return 0;
    }
    if (OUT_PTRS) { OUT_PTRS[0] = nc1; OUT_PTRS[1] = nc2; OUT_PTRS[2] = nc3; }
    return 1;
}

static int can_perform_flip32(const s_ncell *ncell, int opp_cell_id, int *ridge_id_2)
{   /* Checks if tetra abpd exists */
    int face_vid[3]; extract_ids_face(ncell, 2, &opp_cell_id, face_vid);
    int opp_face_lid; face_localid_of_adjacent_ncell(ncell, 2, &opp_cell_id, opp_cell_id, &opp_face_lid);
    
    s_ncell *opp_ncell = ncell->opposite[opp_cell_id];
    for (int ii=0; ii<4; ii++) {
        s_ncell *opp_opp = opp_ncell->opposite[ii];
        if (opp_opp && opp_opp != ncell &&
            inarray(opp_opp->vertex_id, 4, ncell->vertex_id[opp_cell_id]) &&
            inarray(opp_opp->vertex_id, 4, opp_ncell->vertex_id[opp_face_lid])) {
                *ridge_id_2 = id_where_equal_int(ncell->vertex_id, 4, opp_ncell->vertex_id[ii]);
                return 1;
        }
    }
    return 0;
}


static int flip32(s_scplx *scplx, s_dstack *stack, s_ncell *nc1, int opp_cell_id, int ridge_id_2, int opp_face_localid, s_ncell *OUT_PTRS[2])
{
    scplx->N_ncells -= 1;
    s_ncell *nc2, *nc3;
    {
        int v2_main, v2_2, v3_main, v3_2;
        nc2 = next_ncell_ridge_cycle(nc1, opp_cell_id, ridge_id_2, &v2_main, &v2_2);
        nc3 = next_ncell_ridge_cycle(nc2, v2_main, v2_2, &v3_main, &v3_2);
    }
    s_ncell *nc1_opp_old[4]; opposite_pointers_ncell(nc1, nc1_opp_old);
    s_ncell *nc2_opp_old[4]; opposite_pointers_ncell(nc2, nc2_opp_old);
    s_ncell *nc3_opp_old[4]; opposite_pointers_ncell(nc3, nc3_opp_old);

    /* Remove nc3 from linked list */
    s_ncell *nc3_next = nc3->next;
    if (nc3->next) nc3->next->prev = nc3->prev;
    if (nc3->prev) nc3->prev->next = nc3_next;
    else scplx->head = nc3->next;

    /* scplx important indices */
    int lid_ridge[2] = {opp_cell_id, ridge_id_2};
    int vid_ridge[2]; extract_ids_face(nc1, 1, lid_ridge, vid_ridge);
    int p = nc1->vertex_id[opp_cell_id];
    int a = nc1->vertex_id[ridge_id_2];
    int b = vid_ridge[0];
    int c = vid_ridge[1];
    int d = nc2->vertex_id[opp_face_localid];

    int nc1_id_a, nc1_id_b, nc1_id_c, nc1_id_p;
    map_vid_lid(nc1, a, &nc1_id_a, b, &nc1_id_b, c, &nc1_id_c, p, &nc1_id_p);
    int nc2_id_a, nc2_id_b, nc2_id_c, nc2_id_d;
    map_vid_lid(nc2, a, &nc2_id_a, b, &nc2_id_b, c, &nc2_id_c, d, &nc2_id_d);
    int nc3_id_b, nc3_id_c, nc3_id_d, nc3_id_p;
    map_vid_lid(nc3, b, &nc3_id_b, c, &nc3_id_c, d, &nc3_id_d, p, &nc3_id_p);

     /* Update nc1 */
    nc1->vertex_id[nc1_id_c] = d;
    nc1->opposite[nc1_id_b] = nc2;
    nc1->opposite[nc1_id_p] = nc2_opp_old[nc2_id_c];
    nc1->opposite[nc1_id_a] = nc3_opp_old[nc3_id_c];
    if (nc2_opp_old[nc2_id_c]) {
        int opp_aux; face_localid_of_adjacent_ncell(nc1, 2, &nc1_id_p, nc1_id_p, &opp_aux);
        nc2_opp_old[nc2_id_c]->opposite[opp_aux] = nc1;
    }
    if (nc3_opp_old[nc3_id_c]) {
        int opp_aux; face_localid_of_adjacent_ncell(nc1, 2, &nc1_id_a, nc1_id_a, &opp_aux);
        nc3_opp_old[nc3_id_c]->opposite[opp_aux] = nc1;
    }

    /* Update nc2 */
    nc2->vertex_id[nc2_id_b] = p; 
    nc2->opposite[nc2_id_c] = nc1;
    nc2->opposite[nc2_id_d] = nc1_opp_old[nc1_id_b];
    nc2->opposite[nc2_id_a] = nc3_opp_old[nc3_id_b];
    if (nc1_opp_old[nc1_id_b]) {
        int opp_aux; face_localid_of_adjacent_ncell(nc2, 2, &nc2_id_d, nc2_id_d, &opp_aux);
        nc1_opp_old[nc1_id_b]->opposite[opp_aux] = nc2;
    }
    if (nc3_opp_old[nc3_id_b]) {
        int opp_aux; face_localid_of_adjacent_ncell(nc2, 2, &nc2_id_a, nc2_id_a, &opp_aux);
        nc3_opp_old[nc3_id_b]->opposite[opp_aux] = nc2;
    }

    if (stack) stack_remove_ncell(stack, nc3);
    free_ncell(nc3);
    if (stack) { 
        if (!stack_push(stack, nc1)) return 0;
        if (!stack_push(stack, nc2)) return 0;
    }
    if (OUT_PTRS) { OUT_PTRS[0] = nc1; OUT_PTRS[1] = nc2; }
    return 1;
}


static int can_perform_flip44(const s_scplx *scplx, const s_ncell *ncell, int opp_cell_id, int *ridge_id_2)
{   /* In general, in config44 no need for coplanarity. But since this flip is only done on a degenerate case, 
       it makes the test simpler. ridge_id_2 is the other vertex NOT belonging to the ridge. */
    /* Determine ridge along which we are in config44 */
    int opp_face_localid; face_localid_of_adjacent_ncell(ncell, 2, &opp_cell_id, opp_cell_id, &opp_face_localid);
    int face_vid[3]; extract_ids_face(ncell, 2, &opp_cell_id, face_vid);
    s_point face_points[3] = {scplx->points.p[face_vid[0]], scplx->points.p[face_vid[1]], scplx->points.p[face_vid[2]]};
    s_point point_p = scplx->points.p[ncell->vertex_id[opp_cell_id]];
    s_point point_d = scplx->points.p[ncell->opposite[opp_cell_id]->vertex_id[opp_face_localid]];

    int o0 = test_orientation((s_point[3]){point_p, face_points[0], face_points[1]}, point_d);
    int o1 = test_orientation((s_point[3]){point_p, face_points[1], face_points[2]}, point_d);
    int o2 = test_orientation((s_point[3]){point_p, face_points[2], face_points[0]}, point_d);
    assert((o0 == 0) + (o1 == 0) + (o2 == 0) == 1 ||
           (o0 == 0) + (o1 == 0) + (o2 == 0) == 2);

    int AUX_ridge_id_2[2];  int k=0;
    if (o0 == 0) { AUX_ridge_id_2[k++] = id_where_equal_int(ncell->vertex_id, 4, face_vid[2]); } 
    if (o1 == 0) { AUX_ridge_id_2[k++] = id_where_equal_int(ncell->vertex_id, 4, face_vid[0]); } 
    if (o2 == 0) { AUX_ridge_id_2[k++] = id_where_equal_int(ncell->vertex_id, 4, face_vid[1]); } 

    for (int i=0; i<k; i++) {
        *ridge_id_2 = AUX_ridge_id_2[i];

        int nc2_id1, nc2_id2;
        s_ncell *nc2 = next_ncell_ridge_cycle(ncell, opp_cell_id, *ridge_id_2, &nc2_id1, &nc2_id2);
        if (!nc2) continue;

        int nc3_id1, nc3_id2;
        s_ncell *nc3 = next_ncell_ridge_cycle(nc2, nc2_id1, nc2_id2, &nc3_id1, &nc3_id2);
        if (!nc3) continue;

        int nc4_id1, nc4_id2;
        s_ncell *nc4 = next_ncell_ridge_cycle(nc3, nc3_id1, nc3_id2, &nc4_id1, &nc4_id2);
        if (!nc4) continue;

        if (nc4->opposite[nc4_id1] == ncell) return 1;
    }
    // printf("NOT FOUND RIDGE\n");
    return 0;
}

static int flip44(s_scplx *scplx, s_dstack *stack, s_ncell *ncell, int id_ridge_1, int id_ridge_2, s_ncell *OUT_PTRS[4], bool *ignored)
{
    /* 1) flip23 */
    /* Store some indices before flip */
    int opp_face_lid; face_localid_of_adjacent_ncell(ncell, 2, &id_ridge_1, id_ridge_1, &opp_face_lid);
    int opp_face_vid = ncell->opposite[id_ridge_1]->vertex_id[opp_face_lid];  
    int lid_ridge[2] = {id_ridge_1, id_ridge_2};  
    int vid_ridge[2]; extract_ids_face(ncell, 1, lid_ridge, vid_ridge);
    int p = ncell->vertex_id[id_ridge_1];
    int a = vid_ridge[0];
    int c = vid_ridge[1];
    int d = opp_face_vid;

    s_ncell *FLIP23_PTRS[3];
    if (!flip23(scplx, NULL, ncell, id_ridge_1, opp_face_lid, FLIP23_PTRS)) return 0;

    /* Find which ncell added shares ridge. Currently, no way to predict this. */
    s_ncell *nc5 = NULL;
    for (int ii=0; ii<3; ii++) {
        if (inarray(FLIP23_PTRS[ii]->vertex_id, 4, a) && inarray(FLIP23_PTRS[ii]->vertex_id, 4, c) &&
            inarray(FLIP23_PTRS[ii]->vertex_id, 4, d) && inarray(FLIP23_PTRS[ii]->vertex_id, 4, p)) {
            nc5 = FLIP23_PTRS[ii];
            if (stack) {
                if (!stack_push(stack, FLIP23_PTRS[(ii+1)%3])) return 0;
                if (!stack_push(stack, FLIP23_PTRS[(ii+2)%3])) return 0;
            }
            if (OUT_PTRS) { OUT_PTRS[0] = FLIP23_PTRS[(ii+1)%3]; OUT_PTRS[1] = FLIP23_PTRS[(ii+2)%3]; }
            break;
        }
    }
    assert(nc5 != NULL && "Could not perform flip44...");

    /* 2) flip32 */
    int nc5_p = id_where_equal_int(nc5->vertex_id, 4, ncell->vertex_id[id_ridge_1]);
    s_ncell *nc3 = nc5->opposite[id_where_equal_int(nc5->vertex_id, 4, ncell->vertex_id[id_ridge_1])];
    int nc3_id1 = id_where_equal_int(nc3->vertex_id, 4, opp_face_vid);
    int nc3_id2; face_localid_of_adjacent_ncell(nc5, 2, &nc5_p, nc5_p, &nc3_id2);
    int nc3_opp_face_lid; face_localid_of_adjacent_ncell(nc3, 2, &nc3_id1, nc3_id1, &nc3_opp_face_lid);

    s_ncell *FLIP32_PTRS[2];
    flip32(scplx, NULL, nc3, nc3_id1, nc3_id2, nc3_opp_face_lid, FLIP32_PTRS);
    if(stack) {
        if (!stack_push(stack, FLIP32_PTRS[0])) return 0;
        if (!stack_push(stack, FLIP32_PTRS[1])) return 0;

        /* The Bowyer-Watson loop pops each cell and checks exactly one face:
         * opposite[local_p], i.e. the face opposite the newly inserted point p.
         * This is correct for flip14/flip23/flip32, where every face shared between
         * two cells that both contain p is guaranteed locally Delaunay by construction
         * -- those are "interior" star faces that never need re-examination.
         *
         * flip44 breaks that invariant.  It calls flip23 then flip32 on a
         * sub-configuration determined by the CASE_FLAT ridge, not the outer
         * boundary of the star.  The two cells produced by flip32 share a face
         * that CONTAINS p (so it is opposite some old vertex, not opposite p).
         * That face was never an outer boundary candidate at any point in the
         * sequence, so no flip ever pushed it for checking.  When the BW loop
         * pops these two cells it checks their outer faces (opposite p), which
         * are fine -- and the interior shared face goes unexamined.
         *
         * Fix: check that shared face explicitly here, right after flip32 returns.
         * Use NONSTRICT so that coplanar inputs (where the tested point lies exactly
         * on the circumsphere) are accepted without flipping -- otherwise the two
         * valid diagonalisations of a degenerate square cycle forever. */
        for (int fi = 0; fi < 4; fi++) {
            if (FLIP32_PTRS[0]->opposite[fi] == FLIP32_PTRS[1]) {
                if (!are_locally_delaunay(scplx, FLIP32_PTRS[0], fi, DELAUNAY_TEST_NONSTRICT)) {
                    stack_remove_ncell(stack, FLIP32_PTRS[0]);
                    stack_remove_ncell(stack, FLIP32_PTRS[1]);
                    if (flip_tetrahedra(scplx, stack, FLIP32_PTRS[0], fi, ignored) == -1) return -1;
                }
                break;
            }
        }
    }
    if (OUT_PTRS) { OUT_PTRS[2] = FLIP32_PTRS[0]; OUT_PTRS[3] = FLIP32_PTRS[1]; }
    return 1;
}


static int can_perform_flip41(const s_ncell *ncell, int opp_cell_id, int *redundant_localid)
{
    int face_vid[3]; extract_ids_face(ncell, 2, &opp_cell_id, face_vid);
    int opp_face_lid; face_localid_of_adjacent_ncell(ncell, 2, &opp_cell_id, opp_cell_id, &opp_face_lid);
    s_ncell *opp_ncell = ncell->opposite[opp_cell_id];
    int p_vid = ncell->vertex_id[opp_cell_id];

    /* For each edge of shared face, check if degree 3 tetrahedron exists */
    int degree3_edges[3][2], n_degree3 = 0;
    for (int i=0; i<3; i++) {
        int va = face_vid[i], vb = face_vid[(i+1)%3];
        for (int j=0; j<4; j++) {
            s_ncell *nb = opp_ncell->opposite[j];
            if (nb && nb != ncell &&
                inarray(nb->vertex_id, 4, p_vid) &&
                inarray(nb->vertex_id, 4, va)    &&
                inarray(nb->vertex_id, 4, vb)) {
                degree3_edges[n_degree3][0] = va;
                degree3_edges[n_degree3][1] = vb;
                n_degree3++;
                break;
            }
        }
    }

    if (n_degree3 != 2) return 0;

    /* Redundant vertex is common to both degree-3 edges */
    int redundant_vid = -1;
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            if (degree3_edges[0][i] == degree3_edges[1][j])
                redundant_vid = degree3_edges[0][i];

    assert(redundant_vid != -1);
    *redundant_localid = id_where_equal_int(ncell->vertex_id, 4, redundant_vid);
    return 1;
}


static int flip41(s_scplx *scplx, s_dstack *stack, s_ncell *ncell, int r_localid, bool *ignored)
{
    int r_vid = ncell->vertex_id[r_localid];
    scplx->N_ncells -= 3;
    ignored[r_vid] = true;

    /* Determine the 4 ncells in the star of redundant_vid */
    s_ncell *star[4]; 
    star[0] = ncell;
    int k = 1;
    for (int i = 0; i < 4; i++) {
        if (i == r_localid) continue;
        s_ncell *nb = ncell->opposite[i];
        if (nb && inarray(nb->vertex_id, 4, r_vid))
            star[k++] = nb;
    }
    assert(k == 4 && "flip41: redundant vertex does not have exactly 4 tetrahedra in its star");

    /* Precompute old opposite pointers before modifying anything */
    s_ncell *nc1_opp_old[4]; opposite_pointers_ncell(star[0], nc1_opp_old);
    s_ncell *nc2_opp_old[4]; opposite_pointers_ncell(star[1], nc2_opp_old);
    s_ncell *nc3_opp_old[4]; opposite_pointers_ncell(star[2], nc3_opp_old);
    s_ncell *nc4_opp_old[4]; opposite_pointers_ncell(star[3], nc4_opp_old);

    /* Local id of redundant vertex in each ncell */
    int nc1_r = id_where_equal_int(star[0]->vertex_id, 4, r_vid);
    int nc2_r = id_where_equal_int(star[1]->vertex_id, 4, r_vid);
    int nc3_r = id_where_equal_int(star[2]->vertex_id, 4, r_vid);
    int nc4_r = id_where_equal_int(star[3]->vertex_id, 4, r_vid);

    /* 4 vertices that will be kept */
    int a = ncell->vertex_id[(r_localid+1)%4];  
    int b = ncell->vertex_id[(r_localid+2)%4];  
    int c = ncell->vertex_id[(r_localid+3)%4];
    int d = -1;
    for (int i=0; i<4; i++) {
        int aux = star[1]->vertex_id[i];
        if (aux != a && aux != b && aux != c && aux != r_vid) d = aux;
    }
    assert(d != -1);

    /* Update nc1 to be the surviving tetrahedron */
    set_ncell_vids(star[0], a, b, c, d);

    /* Correct opposite pointer correspondence. 
     * out_i is opposite to whichever of {a,b,c,d} is absent from star[i]. */
    int abcd[4] = {a, b, c, d};
    s_ncell *out_raw[4]  = {nc1_opp_old[nc1_r], nc2_opp_old[nc2_r],
                            nc3_opp_old[nc3_r], nc4_opp_old[nc4_r]};
    s_ncell *out[4]      = {NULL, NULL, NULL, NULL};
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            if (!inarray(star[i]->vertex_id, 4, abcd[j]))
                { out[j] = out_raw[i]; break; }

    set_ncell_opposite_pointers(star[0], out[0], out[1], out[2], out[3]);

    /* Update outer neighbors. local ids now correctly match out[] */
    for (int j = 0; j < 4; j++) {
        if (out[j]) {
            int opp_aux; face_localid_of_adjacent_ncell(star[0], 2, &j, j, &opp_aux);
            out[j]->opposite[opp_aux] = star[0];
        }
    }

    /* Remove nc2, nc3, nc4 from linked list */
    if (star[1]->next) star[1]->next->prev = star[1]->prev;
    if (star[1]->prev) star[1]->prev->next = star[1]->next;
    else scplx->head = star[1]->next;

    if (star[2]->next) star[2]->next->prev = star[2]->prev;
    if (star[2]->prev) star[2]->prev->next = star[2]->next;
    else scplx->head = star[2]->next;

    if (star[3]->next) star[3]->next->prev = star[3]->prev;
    if (star[3]->prev) star[3]->prev->next = star[3]->next;
    else scplx->head = star[3]->next;

    if (stack) {
        stack_remove_ncell(stack, star[1]);
        stack_remove_ncell(stack, star[2]);
        stack_remove_ncell(stack, star[3]);
        if (!stack_push(stack, star[0])) return 0;
    }

    free_ncell(star[1]);
    free_ncell(star[2]);
    free_ncell(star[3]);
    return 1;
}



// --------------------------------- MAIN ALGORITHM -------------------------------------------

static void regular_tetrahedron(s_point centre, double inradius, s_point out[4])
{   /* Direct 3D construction. */
    double R = 3.0 * inradius;  /* circumradius */

    /* One vertex on top */
    out[0].x = centre.x;
    out[0].y = centre.y;
    out[0].z = centre.z + R;

    /* Three vertices on bottom circle at height centre.z - inradius */
    double r_circle = sqrt(R*R - inradius*inradius);  /* = R*sqrt(8/9) = 2*sqrt(2)*inradius */
    for (int ii = 0; ii < 3; ii++) {
        double angle = 2.0 * M_PI * ii / 3.0;
        out[ii+1].x = centre.x + r_circle * cos(angle);
        out[ii+1].y = centre.y + r_circle * sin(angle);
        out[ii+1].z = centre.z - inradius;
    }
}

static void perturb_big_tetra(s_point out[4], double inradius)
{   /* Small random perturbation relative to inradius to break degeneracies */
    double eps = inradius * 1e-6;
    for (int ii = 0; ii < 4; ii++) {
        out[ii].x += eps * ((double)rand()/RAND_MAX - 0.5);
        out[ii].y += eps * ((double)rand()/RAND_MAX - 0.5);
        out[ii].z += eps * ((double)rand()/RAND_MAX - 0.5);
    }
}

static int initialize_scplx(const s_points *points, const double *weights, s_scplx *out,
                             const s_point *bb_min_hint, const s_point *bb_max_hint)
{
    /* scplx->points is extended for the extra nodes of big_ncell, put at the beginning. */
    s_points scplx_points = { .N = points->N + 4,
                              .p = malloc(sizeof(s_point) * (points->N + 4)) };
    if (!scplx_points.p) return 0;
    for (int ii=0; ii<points->N; ii++)
        scplx_points.p[ii+4] = points->p[ii];

    /* Initialize the regular tetrahedron big enough to contain all points. */
    s_point bmin, bmax;
    if (points->N > 0) {
        bounding_box_points(points, &bmin, &bmax);
    } else if (bb_min_hint && bb_max_hint) {
        bmin = *bb_min_hint; bmax = *bb_max_hint;
    } else {
        bmin = (s_point){.x=0,.y=0,.z=0}; bmax = (s_point){.x=0,.y=0,.z=0};
    }
    /* Expand bounding box with optional hint (needed when seeds are few/coincident). */
    if (bb_min_hint) {
        bmin.x = fmin(bmin.x, bb_min_hint->x);
        bmin.y = fmin(bmin.y, bb_min_hint->y);
        bmin.z = fmin(bmin.z, bb_min_hint->z);
    }
    if (bb_max_hint) {
        bmax.x = fmax(bmax.x, bb_max_hint->x);
        bmax.y = fmax(bmax.y, bb_max_hint->y);
        bmax.z = fmax(bmax.z, bb_max_hint->z);
    }
    double bdiag = distance(bmin, bmax);
    s_point centre = scale_point(sum_points(bmax, bmin), 0.5);
    /* Add max radius so big tet contains all balls, not just centers. */
    double maxr = 0.0;
    if (weights) {
        for (int ii = 0; ii < points->N; ii++)
            maxr = fmax(maxr, sqrt(weights[ii]));
    }
    /* inradius min: bdiag/2 + 2*maxr, but the bigger it makes less errors. */
    double inradius = (2*bdiag + 8*maxr);    // TODO TEST WHY TIGHTENING RESULTS IN FAILS
    // double inradius = 1.1 * (bdiag/2 + 2*maxr);  
    regular_tetrahedron(centre, inradius, scplx_points.p);
    perturb_big_tetra(scplx_points.p, inradius);
    
    out->points.N = points->N + 4;
    out->points = scplx_points;
    s_ncell *big_ncell = malloc_ncell();
    if (!big_ncell) { free_points(&scplx_points); return 0; }
    for (int ii=0; ii<4; ii++) {
        big_ncell->vertex_id[ii] = ii;
        big_ncell->opposite[ii] = NULL;
    }
    out->head = big_ncell;
    out->N_ncells = 1;
    if (weights) {
        out->weights = malloc(sizeof(double) * scplx_points.N);
        if (!out->weights) {free_ncell(big_ncell); free_points(&scplx_points); return 0; }
        for (int ii=0; ii<4; ii++) out->weights[ii] = 0;
        for (int ii=4, jj=0; ii<points->N+4; ii++) out->weights[ii] = weights[jj++];
    } else out->weights = NULL;

    return 1;
}


typedef enum type_union_tetra {
    CASE_CONVEX,
    CASE_NON_CONVEX,
    CASE_FLAT,
    CASE_P_IN_EDGE,
    CASE_ERROR
} e_type_union_tetra;

static e_type_union_tetra determine_case(const s_point vertices_face[3], s_point p, s_point d) 
{   /* two ncells sharing face abc, with opposite vertices p and d */
    e_geom_test test_p_in_face = test_point_in_triangle_3D(vertices_face, p, 0, 0);
    if (test_p_in_face == TEST_BOUNDARY) return CASE_P_IN_EDGE;

    /* Checks intersection of segment pd with face */
    switch (test_segment_triangle_intersect_3D((s_point[2]){p,d}, vertices_face, 0, 0)) {
        case INTERSECT_DEGENERATE:  
            /* either pd has endpoint inside abc, or pd intersects edge/vertex */
            if (test_p_in_face == TEST_IN) return CASE_CONVEX;
            else return CASE_FLAT;
        case INTERSECT_EMPTY: return CASE_NON_CONVEX;
        case INTERSECT_NONDEGENERATE: return CASE_CONVEX;
        case INTERSECT_ERROR: 
            printf("ERROR IN DETERMINE CASE! Degenerate triangle?\n");
            return CASE_ERROR;
    }
}

static int flip_tetrahedra(s_scplx *scplx, s_dstack *stack, s_ncell *ncell, int opp_cell_id, bool *ignored)
{   /* -1 ERROR, 0 NOT FLIPPED, 1 FLIPPED */
    if (!ncell->opposite[opp_cell_id]) return 0;

    s_point coords_face[3];
    extract_vertices_face(scplx, ncell, 2, &opp_cell_id, coords_face);

    int opp_face_localid; 
    face_localid_of_adjacent_ncell(ncell, 2, &opp_cell_id, opp_cell_id, &opp_face_localid);
    int opp_face_vertex_id = (ncell->opposite[opp_cell_id])->vertex_id[opp_face_localid];

    s_point p = scplx->points.p[ncell->vertex_id[opp_cell_id]];
    s_point d = scplx->points.p[opp_face_vertex_id];
    
    switch (determine_case(coords_face, p, d)) {
        int ridge_id_2;
        int redundant_localid;
        case CASE_ERROR: return 0;
        case CASE_CONVEX:
            if (!flip23(scplx, stack, ncell, opp_cell_id, opp_face_localid, NULL)) return -1;
            else return 1;
        case CASE_NON_CONVEX:
            if (scplx->weights &&
                can_perform_flip41(ncell, opp_cell_id, &redundant_localid)) {
                if (!flip41(scplx, stack, ncell, redundant_localid, ignored)) return -1;
                return 1;
            } else if (can_perform_flip32(ncell, opp_cell_id, &ridge_id_2)) {
                if (!flip32(scplx, stack, ncell, opp_cell_id, ridge_id_2, opp_face_localid, NULL)) return -1;
                return 1;
            } else return 0;
        case CASE_FLAT:
            if (can_perform_flip44(scplx, ncell, opp_cell_id, &ridge_id_2)) {
                if (flip44(scplx, stack, ncell, opp_cell_id, ridge_id_2, NULL, ignored) == -1) return -1;
                return 1;
            } else return 0;
        case CASE_P_IN_EDGE:
            if (!flip23(scplx, stack, ncell, opp_cell_id, opp_face_localid, NULL)) return -1;
            else return 1;
    }
    return 0;
}

static bool point_close_to_ncell_vertex(s_scplx *scplx, s_ncell *ncell, s_point point, double TOL)
{
    const double TOL2 = TOL*TOL;
    if (distance_squared(scplx->points.p[ncell->vertex_id[0]], point) <= TOL2 || 
        distance_squared(scplx->points.p[ncell->vertex_id[1]], point) <= TOL2 ||
        distance_squared(scplx->points.p[ncell->vertex_id[2]], point) <= TOL2 ||
        distance_squared(scplx->points.p[ncell->vertex_id[3]], point) <= TOL2)
        return true;
    else return false;
}

static int insert_one_point(s_scplx *scplx, int point_id, s_dstack *stack, double TOL_dup, bool *ignored)
{   /* -1: ERROR, 0: not inserted, 1: inserted */
    s_point point = scplx->points.p[point_id];
    s_ncell *container_ncell = in_ncell_walk(scplx, point);

    if (point_close_to_ncell_vertex(scplx, container_ncell, point, TOL_dup) ||
        p_locally_redundant_in_ncell(scplx, container_ncell, point_id)) {
        ignored[point_id] = true;
        return 0;
    }

    if (!flip14(scplx, container_ncell, point_id, stack)) return -1;

    /* Bowyer-Watson flip loop.  Each entry on the stack is a cell containing
     * point_id.  We check only opposite[local_p] -- the single face of that
     * cell that is opposite point_id (an outer boundary face of the star of p).
     * The invariant is: interior faces of the star (shared between two cells
     * that both contain p) are always locally Delaunay after flip14/flip23/flip32.
     * flip44 breaks this invariant and must therefore check its own interior
     * faces explicitly before returning (see the comment there). */
    while (stack->size > 0) {
        s_ncell *current = stack_pop(stack);
        int opp_cell_id = id_where_equal_int(current->vertex_id, 4, point_id);
        if (!current->opposite[opp_cell_id]) continue;
        if (!are_locally_delaunay(scplx, current, opp_cell_id, DELAUNAY_TEST_NONSTRICT)) {
            if (flip_tetrahedra(scplx, stack, current, opp_cell_id, ignored) == -1) return -1;
        }
    }
    return 1;
}

/*
 * kept_idx, if non-NULL, is composed in place: for each i in [0, N_kept_idx),
 * if kept_idx[i] is a valid index into the pre-compaction "real seed" range
 * (i.e. was produced by an earlier filtering stage as an index there), it is
 * rewritten to that seed's final compacted index, or -1 if the seed was
 * dropped here as a near-duplicate. Reuses the remap table already built
 * for compaction below -- no extra allocation needed.
 */
static void remove_ignored_points(s_scplx *scplx, bool *ignored, bool keep_big_tetra,
                                  int *kept_idx, int N_kept_idx)
{
    if (!keep_big_tetra) {
        /* Mark first 4 indices as ignored */
        for (int i = 0; i < 4; i++) ignored[i] = true;
        
        /* Remove any ncell referencing big tetra vertices */
        s_ncell *current = scplx->head;
        while (current) {
            s_ncell *next = current->next;
            for (int ii=0; ii<4; ii++) if (current->vertex_id[ii] < 4) { 
                if (current->next) (current->next)->prev = current->prev;
                if (current->prev) (current->prev)->next = next;
                else scplx->head = current->next;

                /* Update opposite's opposite to NULL */
                for (int jj=0; jj<4; jj++) if (current->opposite[jj]) {
                    for (int kk=0; kk<4; kk++) if (current->opposite[jj]->opposite[kk] == current) {
                        current->opposite[jj]->opposite[kk] = NULL;
                        break;
                    }
                }

                free_ncell(current);
                scplx->N_ncells--;
                break;
            }
            current = next;
        }
    }

    /* Compact points / weights, and build remap table */
    int *remap = malloc(sizeof(int) * scplx->points.N);
    int k = 0;
    for (int i = 0; i < scplx->points.N; i++) {
        if (ignored[i]) { remap[i] = -1; continue; }
        scplx->points.p[k] = scplx->points.p[i];
        if (scplx->weights) scplx->weights[k] = scplx->weights[i];
        remap[i] = k++;
    }
    scplx->points.N = k;
    scplx->points.p = realloc(scplx->points.p, sizeof(s_point) * k);
    if (scplx->weights) scplx->weights = realloc(scplx->weights, sizeof(double) * k);

    /* Update vertex ids */
    for (s_ncell *c = scplx->head; c; c = c->next) {
        for (int i = 0; i < 4; i++) {
            assert(remap[c->vertex_id[i]] != -1 && "Ignored point still referenced by tetra.");
            c->vertex_id[i] = remap[c->vertex_id[i]];
        }
    }

    if (kept_idx)
        for (int i = 0; i < N_kept_idx; i++)
            if (kept_idx[i] >= 0)
                kept_idx[i] = remap[4 + kept_idx[i]];

    free(remap);
}




/* BUILDER */
s_dt_builder dt_builder_begin(const s_points *seeds, const double *weights, double TOL_dup,
                              const s_point *bb_min_hint, const s_point *bb_max_hint)
{
    s_dt_builder b = {0};

    bool *ignored = calloc(seeds->N + 4, sizeof(bool));
    if (!ignored) return b;

    s_dstack *stack = malloc(sizeof(s_dstack));
    if (!stack) { free(ignored); return b; }
    *stack = stack_create();
    if (!stack->entry) { free(ignored); free(stack); return b; }

    if (!initialize_scplx(seeds, weights, &b.dt, bb_min_hint, bb_max_hint)) {
        free(ignored); stack_free(stack); free(stack); return b;
    }

    for (int ii = 4; ii < b.dt.points.N; ii++) {
        if (insert_one_point(&b.dt, ii, stack, TOL_dup, ignored) == -1) {
            free(ignored); stack_free(stack); free(stack);
            free_complex(&b.dt);
            return (s_dt_builder){0};
        }
    }

    b._ignored = ignored;
    b._stack = stack;
    b._N_original = seeds->N;
    return b;
}


bool dt_builder_extend(s_dt_builder *b, const s_points *new_points, double TOL_dup)
{
    int N_old = b->dt.points.N;
    int N_add = new_points->N;
    int N_new = N_old + N_add;

    s_point *tmp_p = realloc(b->dt.points.p, N_new * sizeof(s_point));
    if (!tmp_p) return false;
    b->dt.points.p = tmp_p;
    memcpy(&b->dt.points.p[N_old], new_points->p, N_add * sizeof(s_point));
    b->dt.points.N = N_new;

    bool *tmp_ig = realloc(b->_ignored, N_new * sizeof(bool));
    if (!tmp_ig) return false;
    b->_ignored = tmp_ig;
    memset(&b->_ignored[N_old], 0, N_add * sizeof(bool));

    if (b->dt.weights) {
        double *tmp_w = realloc(b->dt.weights, N_new * sizeof(double));
        if (!tmp_w) return false;
        b->dt.weights = tmp_w;
        for (int i = N_old; i < N_new; i++) b->dt.weights[i] = 0.0;
    }

    s_dstack *stack = (s_dstack *)b->_stack;
    for (int i = N_old; i < N_new; i++)
        if (insert_one_point(&b->dt, i, stack, TOL_dup, b->_ignored) == -1)
            return false;

    return true;
}


s_scplx dt_builder_end(s_dt_builder *b, bool keep_big_tetra, int *out_Nreal,
                       int *kept_idx, int N_kept_idx)
{
    // if (!is_delaunay_3d(&b->dt, DELAUNAY_TEST_NONSTRICT)) {
    //     fprintf(stderr, "WARNING: DT is not Delaunay!\n");
    //     write_points_to_csv("error.csv", "w", &b->dt.points);
    //     // s_ncell *c = b->dt.head;
    //     // while (c) {
    //     //     s_point v[4]; extract_vertices_ncell(&b->dt, c, v);
    //     //     if (test_orientation(v, v[3]) == 0) {
    //     //         fprintf(stderr, "Flat tetra coords:\n");
    //     //         fprintf(stderr, "  v0(%d): %.6f %.6f %.6f\n", c->vertex_id[0], v[0].x, v[0].y, v[0].z);
    //     //         fprintf(stderr, "  v1(%d): %.6f %.6f %.6f\n", c->vertex_id[1], v[1].x, v[1].y, v[1].z);
    //     //         fprintf(stderr, "  v2(%d): %.6f %.6f %.6f\n", c->vertex_id[2], v[2].x, v[2].y, v[2].z);
    //     //         fprintf(stderr, "  v3(%d): %.6f %.6f %.6f\n", c->vertex_id[3], v[3].x, v[3].y, v[3].z);
    //     //     }
    //     //     c = c->next;
    //     // }
    //     exit(1);
    // }

    for (s_ncell *c = b->dt.head; c; c = c->next) {
        s_point v[4]; extract_vertices_ncell(&b->dt, c, v);
        if (test_orientation(v, v[3]) == 0)
            fprintf(stderr, "Flat tetra in final result (%p): %d %d %d %d\n",
                (void *)c,
                c->vertex_id[0], c->vertex_id[1],
                c->vertex_id[2], c->vertex_id[3]);
    }

    if (out_Nreal != NULL) {
        int surviving = 0;
        for (int i = 4; i < 4 + b->_N_original; i++)
            if (!b->_ignored[i]) surviving++;
        *out_Nreal = surviving;
    }

    remove_ignored_points(&b->dt, b->_ignored, keep_big_tetra, kept_idx, N_kept_idx);

    s_dstack *stack = (s_dstack *)b->_stack;
    stack_free(stack);
    free(stack);
    free(b->_ignored);
    b->_stack = NULL;
    b->_ignored = NULL;

    return b->dt;
}


s_scplx construct_dt_3d(const s_points *points, const double *weights,
                        bool keep_big_tetra, double TOL_duplicates, int *out_Nreal)
{
    s_dt_builder b = dt_builder_begin(points, weights, TOL_duplicates, NULL, NULL);
    if (!b._stack) {
        fprintf(stderr, "construct_dt_3d: Error.\n");
        return (s_scplx){0};
    }
    /* Preserve the old behaviour: *out_Nreal on entry specifies how many of the
     * first points->N entries are "real" seeds (the rest are mirrors added by
     * extend_sites_mirroring).  Store that count so dt_builder_end uses it. */
    if (out_Nreal != NULL && *out_Nreal <= points->N)
        b._N_original = *out_Nreal;
    return dt_builder_end(&b, keep_big_tetra, out_Nreal, NULL, 0);
}




int N_ncells_not_big_tetra(s_scplx *scplx)
{
    int N_remove = 0;

    s_ncell *current = scplx->head;
    while (current) {
        for (int ii=0; ii<4; ii++) if (current->vertex_id[ii] < 4) {
            N_remove++;
            break;
        }
        current = current->next;
    }
    return scplx->N_ncells - N_remove;
}

