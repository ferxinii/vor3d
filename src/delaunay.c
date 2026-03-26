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


static int id_where_equal_int(const int *arr, int N, int entry) 
{
    for (int ii=0; ii<N; ii++) if (arr[ii] == entry) return ii;
    fprintf(stderr, "id_where_equal_int: Could not find id.\n"); 
    assert(1==0);
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

static void stack_push(s_dstack *stack, s_ncell *ncell)
{
    for (int ii=0; ii<stack->size; ii++) /* Avoid duplicates */
        if (stack->entry[ii] == ncell) return;

    if (stack->size == stack->capacity) { /* Expand if needed */
        stack->capacity *= 2;
        s_ncell **tmp = realloc(stack->entry, stack->capacity * sizeof(s_ncell *));
        assert(tmp && "realloc failed when increasing stack memory");
        stack->entry = tmp;
    }

    stack->entry[stack->size++] = ncell;
}

static s_ncell *stack_pop(s_dstack *stack)
{   
    if (stack->size == 0) return NULL;
    stack->size--;
    return stack->entry[stack->size];
}

static void stack_remove_ncell(s_dstack *stack, s_ncell *ncell) {
    int newSize = 0;
    for (int ii = 0; ii < stack->size; ii++)
        if (stack->entry[ii] != ncell)  /* Only copy entries that are not the target. */
            stack->entry[newSize++] = stack->entry[ii];
    stack->size = newSize;
}



// ----------------------------------- FLIPS ---------------------------------------------

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
    if (stack) { stack_push(stack, nc1); stack_push(stack, nc2); stack_push(stack, nc3); stack_push(stack, nc4); }
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
    nc3->next = nc2->next; if (nc3->next) nc3->next->prev = nc3;  /* Last migh be null! */
    nc2->next = nc3; nc3->prev = nc2;

    /* scplx important indices */
    int face_vid[3]; extract_ids_face(nc1, 2, &opp_cell_id, face_vid);
    int a = face_vid[0];
    int b = face_vid[1];
    int c = face_vid[2];
    int d = nc2->vertex_id[opp_face_localid];
    int p = nc1->vertex_id[opp_cell_id];

    int nc1_id_a, nc1_id_b, nc1_id_c, nc1_id_p;  /* Map from global vertex ids to local ids */
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

    /* scplx nc3 */
    set_ncell_vids(nc3,    p, c, d, a);
    set_ncell_opposite_pointers(nc3,  nc2_opp_old[nc2_id_b], nc1, nc1_opp_old[nc1_id_b], nc2);
    if (nc1_opp_old[nc1_id_b]) {
        int opp_aux; face_localid_of_adjacent_ncell(nc3, 2, &(int){2}, 2, &opp_aux);
        nc1_opp_old[nc1_id_b]->opposite[opp_aux] = nc3;
    }
    if (nc2_opp_old[nc2_id_b]) {
        int opp_aux; face_localid_of_adjacent_ncell(nc3, 2, &(int){0}, 0, &opp_aux);
        nc2_opp_old[nc2_id_b]->opposite[opp_aux] = nc3;
    }
    
    if (stack) { stack_push(stack, nc1); stack_push(stack, nc2); stack_push(stack, nc3); }
    if (OUT_PTRS) { OUT_PTRS[0] = nc1; OUT_PTRS[1] = nc2; OUT_PTRS[2] = nc3; }
    return 1;
}


static int can_perform_flip32(const s_scplx *scplx, const s_ncell *ncell, int opp_cell_id, int *ridge_id_2)
{   /* Checks if tetra abpd exists */
    if (scplx->N_ncells < 3) return 0;

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

static void flip32(s_scplx *scplx, s_dstack *stack, s_ncell *nc1, int opp_cell_id, int ridge_id_2, int opp_face_localid, s_ncell *OUT_PTRS[2])
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

    /* We remove nc3 */
    /* ----- NC3->PREV ----- NC3 ----- NC3->NEXT ----- */
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

    int nc1_id_a, nc1_id_b, nc1_id_c, nc1_id_p;  /* Map from global vertex ids to local ids */
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
    if (stack) { stack_push(stack, nc1); stack_push(stack, nc2); }
    if (OUT_PTRS) { OUT_PTRS[0] = nc1; OUT_PTRS[1] = nc2; }
}


static int can_perform_flip44(const s_scplx *scplx, const s_ncell *ncell, int opp_cell_id, int *ridge_id_2)
{   /* In general, in config44 no need for coplanarity. But since this flip is only done on a degenerate case, 
       it makes the test simpler. */
    if (scplx->N_ncells < 4) return 0;

    /* Determine ridge */
    int opp_face_localid; face_localid_of_adjacent_ncell(ncell, 2, &opp_cell_id, opp_cell_id, &opp_face_localid);
    int face_vid[3]; extract_ids_face(ncell, 2, &opp_cell_id, face_vid);
    s_point face_points[3] = {scplx->points.p[face_vid[0]], scplx->points.p[face_vid[1]], scplx->points.p[face_vid[2]]};
    s_point point_p = scplx->points.p[ncell->vertex_id[opp_cell_id]];
    s_point point_d = scplx->points.p[ncell->opposite[opp_cell_id]->vertex_id[opp_face_localid]];

    if (orientation_robust((s_point[3]){point_p, face_points[0], face_points[1]}, point_d) == 0) {
        *ridge_id_2 = id_where_equal_int(ncell->vertex_id, 4, face_vid[2]);
    } else if (orientation_robust((s_point[3]){point_p, face_points[1], face_points[2]}, point_d) == 0) {
        *ridge_id_2 = id_where_equal_int(ncell->vertex_id, 4, face_vid[0]);
    } else if (orientation_robust((s_point[3]){point_p, face_points[2], face_points[0]}, point_d) == 0) {
        *ridge_id_2 = id_where_equal_int(ncell->vertex_id, 4, face_vid[1]);
    } else {
        assert(1 == 0 && "Should never reach this! Edge triangle intersect already gave a degenerate intersection...");
    }

    /* Check that in config44 by traversing ridge */
    int nc2_id1, nc2_id2;
    s_ncell *nc2 = next_ncell_ridge_cycle(ncell, opp_cell_id, *ridge_id_2, &nc2_id1, &nc2_id2);
    if (!nc2) return 0;
    int nc3_id1, nc3_id2;
    s_ncell *nc3 = next_ncell_ridge_cycle(nc2, nc2_id1, nc2_id2, &nc3_id1, &nc3_id2);
    if (!nc3) return 0;
    int nc4_id1, nc4_id2;
    s_ncell *nc4 = next_ncell_ridge_cycle(nc3, nc3_id1, nc3_id2, &nc4_id1, &nc4_id2);
    if (!nc4) return 0;
    return (nc4->opposite[nc4_id1] == ncell);

}

static int flip44(s_scplx *scplx, s_dstack *stack, s_ncell *ncell, int id_ridge_1, int id_ridge_2, s_ncell *OUT_PTRS[4]) 
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
    if (!flip23(scplx, NULL, ncell, id_ridge_1, opp_face_lid, FLIP23_PTRS)) return 0;  /* TOWARDS NC2 */

    /* Find which ncell added shares ridge. Currently, no way to predict this. */
    s_ncell *nc5;
    int debug_found = 0;
    for (int ii=0; ii<3; ii++) {
        if (inarray(FLIP23_PTRS[ii]->vertex_id, 4, a) && inarray(FLIP23_PTRS[ii]->vertex_id, 4, c) &&
            inarray(FLIP23_PTRS[ii]->vertex_id, 4, d) && inarray(FLIP23_PTRS[ii]->vertex_id, 4, p)) {
            nc5 = FLIP23_PTRS[ii];
            if (stack) { stack_push(stack, FLIP23_PTRS[(ii+1)%3]); stack_push(stack, FLIP23_PTRS[(ii+2)%3]); }
            if (OUT_PTRS) { OUT_PTRS[0] = FLIP23_PTRS[(ii+1)%3]; OUT_PTRS[1] = FLIP23_PTRS[(ii+2)%3]; }
            debug_found = 1;
            break;
        }
    }
    assert(debug_found == 1 && "Could not perform flip44...");

    /* 2) flip32 */
    int nc5_p = id_where_equal_int(nc5->vertex_id, 4, ncell->vertex_id[id_ridge_1]);
    s_ncell *nc3 = nc5->opposite[id_where_equal_int(nc5->vertex_id, 4, ncell->vertex_id[id_ridge_1])];
    int nc3_id1 = id_where_equal_int(nc3->vertex_id, 4, opp_face_vid);
    int nc3_id2; face_localid_of_adjacent_ncell(nc5, 2, &nc5_p, nc5_p, &nc3_id2);
    int nc3_opp_face_lid; face_localid_of_adjacent_ncell(nc3, 2, &nc3_id1, nc3_id1, &nc3_opp_face_lid);

    s_ncell *FLIP32_PTRS[2];
    flip32(scplx, NULL, nc3, nc3_id1, nc3_id2, nc3_opp_face_lid, FLIP32_PTRS);
    if(stack) { stack_push(stack, FLIP32_PTRS[0]); stack_push(stack, FLIP32_PTRS[1]); }
    if (OUT_PTRS) { OUT_PTRS[2] = FLIP32_PTRS[0]; OUT_PTRS[3] = FLIP32_PTRS[1]; }
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

static int initialize_scplx(const s_points *points, const double *weights, s_scplx *out)
{
    /* scplx->points is extended for the extra nodes of big_ncell, put at the beginning. */
    s_points scplx_points = { .N = points->N + 4, 
                              .p = malloc(sizeof(s_point) * (points->N + 4)) };
    if (!scplx_points.p) return 0;
    for (int ii=0; ii<points->N; ii++)
        scplx_points.p[ii+4] = points->p[ii];
    
    /* Initialize the regular tetrahedron big enough to contain all points. */
    s_point bmin, bmax; bounding_box_points(points, &bmin, &bmax);
    double bdiag = distance(bmin, bmax);
    s_point centre = scale_point(sum_points(bmax, bmin), 0.5);
    /* Add max radius so big tet contains all balls, not just centers. */
    double maxr = 0.0;
    if (weights) {
        for (int ii = 0; ii < points->N; ii++)
            maxr = fmax(maxr, sqrt(weights[ii]));
    }
    /* inradius min: bdiag/2 + 2*maxr, but the bigger it makes less errors. */
    double inradius = (2*bdiag + 8*maxr);  
    // double inradius = 1.1 * (bdiag/2 + 2*maxr);  
    regular_tetrahedron(centre, inradius, scplx_points.p);
    perturb_big_tetra(scplx_points.p, inradius);
    
    for (int i=0; i<points->N; i++) {
        if (test_point_in_tetrahedron(scplx_points.p, points->p[i], 1e-12, 1e-12) != TEST_IN) {
            printf("OUTSIDE\n");
        }
    }

    out->points.N = points->N + 4;
    out->points = scplx_points;
    s_ncell *big_ncell = malloc_ncell();
    if (!big_ncell) return 0;
    for (int ii=0; ii<4; ii++) {
        big_ncell->vertex_id[ii] = ii;
        big_ncell->opposite[ii] = NULL;
    }
    out->head = big_ncell;
    out->N_ncells = 1;
    if (weights) {
        out->weights = malloc(sizeof(double) * scplx_points.N);
        if (!out->weights) return 0;
        for (int ii=0; ii<4; ii++) out->weights[ii] = 0;
        for (int ii=4, jj=0; ii<points->N+4; ii++) out->weights[ii] = weights[jj++];
    } else out->weights = NULL;

    return 1;
}

static int determine_case(const s_point vertices_face[3], s_point p, s_point d) 
{
    if (test_point_in_triangle_3D(vertices_face, p, 0, 0) == TEST_BOUNDARY) return 4;
    e_intersect_type type = test_segment_triangle_intersect_3D((s_point[2]){p,d}, vertices_face, 0, 0);
    if (type == INTERSECT_DEGENERATE) return 3;
    if (type == INTERSECT_EMPTY) return 2;
    if (type == INTERSECT_NONDEGENERATE) return 1;
    fprintf(stderr, "determine_case: Should never reach this!");
    assert(1 == 0);
}

static int flip_tetrahedra(s_scplx *scplx, s_dstack *stack, s_ncell *ncell, int opp_cell_id)
{   /* -1 ERROR, 0 NOT FLIPPED, 1 FLIPPED */
    s_point coords_face[3];
    extract_vertices_face(scplx, ncell, 2, &opp_cell_id, coords_face);

    // Extract id of vertex in opposite cell corresponding to the face
    int opp_face_localid;
    face_localid_of_adjacent_ncell(ncell, 2, &opp_cell_id, opp_cell_id, &opp_face_localid);
    int opp_face_vertex_id = (ncell->opposite[opp_cell_id])->vertex_id[opp_face_localid];

    s_point p = scplx->points.p[ncell->vertex_id[opp_cell_id]];
    s_point d = scplx->points.p[opp_face_vertex_id];
    
    switch(determine_case(coords_face, p, d)) {
        int ridge_id_2;
        case 1:
            if (!flip23(scplx, stack, ncell, opp_cell_id, opp_face_localid, NULL)) return -1;
            else return 1;
        case 2:
            if (can_perform_flip32(scplx, ncell, opp_cell_id, &ridge_id_2)) {
                flip32(scplx, stack, ncell, opp_cell_id, ridge_id_2, opp_face_localid, NULL);
                return 1;
            } else return 0; 
        case 3:
            if (can_perform_flip44(scplx, ncell, opp_cell_id, &ridge_id_2)) {
                // printf("CASE 3\n");
                flip44(scplx, stack, ncell, opp_cell_id, ridge_id_2, NULL);
                return 1;
            } else return 0;
        case 4:
            fprintf(stderr, "DEBUG delaunay.c: flip_tetrahedra: CASE 4... UNSURE, UNTESTED\n");
            if (!flip23(scplx, stack, ncell, opp_cell_id, opp_face_localid, NULL)) return -1;
            else return 1;
    }
    return 0;
}

static void remove_point_scplx(s_scplx *scplx, int point_id)
{
    if (point_id < scplx->points.N-1)
        for (int ii=point_id; ii<scplx->points.N-1; ii++)
            scplx->points.p[ii] = scplx->points.p[ii+1];

    scplx->points.p = realloc(scplx->points.p, sizeof(s_point) * (scplx->points.N-1));
    scplx->points.N--;

    if (scplx->weights) {
        if (point_id < scplx->points.N)  // N already decremented
            for (int ii=point_id; ii<scplx->points.N; ii++)
                scplx->weights[ii] = scplx->weights[ii+1];
        scplx->weights = realloc(scplx->weights, sizeof(double) * scplx->points.N);
    }
}

static int point_close_to_ncell_vertex(s_scplx *scplx, s_ncell *ncell, s_point point, double TOL)
{
    const double TOL2 = TOL*TOL;
    if (distance_squared(scplx->points.p[ncell->vertex_id[0]], point) <= TOL2 || 
        distance_squared(scplx->points.p[ncell->vertex_id[1]], point) <= TOL2 ||
        distance_squared(scplx->points.p[ncell->vertex_id[2]], point) <= TOL2 ||
        distance_squared(scplx->points.p[ncell->vertex_id[3]], point) <= TOL2)
        return 1;
    else return 0;
}

static int insert_one_point(s_scplx *scplx, int point_id, s_dstack *stack, double TOL_dup)
{   /* -1: ERROR, 0: not inserted, 1: inserted */
    s_point point = scplx->points.p[point_id];
    s_ncell *container_ncell = in_ncell_walk(scplx, point);

    if (point_close_to_ncell_vertex(scplx, container_ncell, point, TOL_dup)) {
        remove_point_scplx(scplx, point_id);
        return 0;
    }

    if (!flip14(scplx, container_ncell, point_id, stack)) return -1;

    while (stack->size > 0) {
        s_ncell *current = stack_pop(stack);
        int opp_cell_id = id_where_equal_int(current->vertex_id, 4, point_id);
        if (current->opposite[opp_cell_id] && !are_locally_delaunay(scplx, current, opp_cell_id, DELAUNAY_TEST_STRICT))
            if (flip_tetrahedra(scplx, stack, current, opp_cell_id) == -1) return -1;
    }
    return 1;
}

static void remove_big_tetra(s_scplx *scplx)
{
    s_ncell *current = scplx->head;
    while (current) {
        s_ncell *next = current->next;
        for (int ii=0; ii<4; ii++) if (current->vertex_id[ii] < 4) {  // This checks if a vertex is part of BIG TETRA
            if (current->next) (current->next)->prev = current->prev;
            if (current->prev) (current->prev)->next = next;
            else scplx->head = current->next;

            // Update opposite's opposite to NULL
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

    /* Update points and their ids */
    for (int ii=0; ii<scplx->points.N-4; ii++)
        scplx->points.p[ii] = scplx->points.p[ii+4];

    scplx->points.p = realloc(scplx->points.p, sizeof(s_point) * (scplx->points.N-4));
    scplx->points.N -= 4;

    for (s_ncell *c=scplx->head; c; c=c->next) 
        for (int kk=0; kk<4; kk++)
            c->vertex_id[kk] -= 4;

    /* Weights */
    if (scplx->weights) {
        for (int ii=0; ii<scplx->points.N; ii++)   // N already decremented by 4
            scplx->weights[ii] = scplx->weights[ii+4];
        scplx->weights = realloc(scplx->weights, sizeof(double) * scplx->points.N);
    }
}


s_scplx construct_dt_3d(const s_points *points, const double *weights, bool keep_big_tetra, double TOL_duplicates)
{
    s_dstack stack = stack_create();
    s_scplx scplx; if (!initialize_scplx(points, weights, &scplx)) goto error;

    int ii = 4;  /* First 4 are big tetra, which already is inserted */
    while (ii < scplx.points.N) {
        int o = insert_one_point(&scplx, ii, &stack, TOL_duplicates);
        if (o == 1) ii++;
        else if (o == -1) goto error;
    }
    
    if (!keep_big_tetra) remove_big_tetra(&scplx);

    if (!is_delaunay_3d(&scplx, DELAUNAY_TEST_NONSTRICT))
        fprintf(stderr, "WARNING: DT is not Delaunay!\n");

    stack_free(&stack);  
    return scplx;

    error:
        fprintf(stderr, "construct_dt_3d: Error.\n");
        return (s_scplx){0};
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

