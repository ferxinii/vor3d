#ifndef VOR3D_DELAUNAY
#define VOR3D_DELAUNAY

#include "simplical_complex.h"


#define MAX_STACK_SIZE 10000


typedef struct stack {
    s_ncell *entry[MAX_STACK_SIZE];
    int size;
} s_stack;


s_stack *stack_create(void);
void stack_free(s_stack *stack);
void stack_push(s_stack *stack, s_ncell *ncell);
s_ncell *stack_pop(s_stack *stack);
s_ncell *stack_peek(s_stack *stack);
void stack_shuffle(s_stack *stack);
void stack_remove_ncell(s_stack *stack, s_ncell *ncell);
void stack_remove_entry(s_stack *stack, int id);
void stack_print(s_stack *stack);


void flip14(s_scplx *setup, s_ncell *container_ncell, int point_id, s_stack *stack);
void flip23(s_scplx *setup, s_stack *stack, s_ncell *nc1, int opp_cell_id, int opp_face_localid, s_ncell **OUT_PTRS);
int can_perform_flip32(const s_scplx *setup, const s_ncell *ncell, int opp_cell_id, int *ridge_id_2);
void flip32(s_scplx *setup, s_stack *stack, s_stack *stack_blocked, s_ncell *nc1, int opp_cell_id, int ridge_id_2, int opp_face_localid, s_ncell **OUT_PTRS);
int can_perform_flip44(const s_scplx *setup, const s_ncell *ncell, s_point *vertices_face, int opp_cell_id, int *ridge_id_2);
void flip44(s_scplx *setup, s_stack *stack, s_stack *stack_blocked, s_ncell *ncell, int id_ridge_1, int id_ridge_2, s_ncell **OUT_PTRS);


int pd_intersects_edge(s_point s1, s_point s2, s_point p, s_point d, int drop_coord);
int is_case_3(const s_point vertices_face[3], s_point p, s_point d);
int is_case_1(const s_point vertices_face[3], s_point p, s_point d);
int is_case_2(const s_point vertices_face[3], s_point p, s_point d);
int is_case_4(const s_point vertices_face[3], s_point p, s_point d);
int determine_case(const s_point vertices_face[3], s_point p, s_point d);


s_scplx *initialize_setup(const s_point *points, int N_points, int dim);
int flip_tetrahedra(s_scplx *setup, s_stack *stack, s_stack *stack_blocked, s_ncell *ncell, int opp_cell_id);
void remove_point_setup(s_scplx *setup, int point_id);
int insert_one_point(s_scplx *setup, int point_id, s_stack *stack, s_stack *stack_blocked);
void remove_big_tetra(s_scplx *setup);
int count_valid_ncells_reduced_triangulation(const s_scplx *setup);


s_scplx *construct_dt_3d(const s_point *points, int N_points);  // Copy of points inside s_scplx

#endif
