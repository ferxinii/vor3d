#ifndef DT_3D_INCREMENTAL_H
#define DT_3D_INCREMENTAL_H

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


void flip14(s_setup *setup, s_ncell *container_ncell, int point_id, s_stack *stack);
void flip23(s_setup *setup, s_stack *stack, s_ncell *nc1, int opp_cell_id, int opp_face_localid, s_ncell **OUT_PTRS);
int can_perform_flip32(const s_setup *setup, const s_ncell *ncell, int opp_cell_id, int *ridge_id_2);
void flip32(s_setup *setup, s_stack *stack, s_stack *stack_blocked, s_ncell *nc1, int opp_cell_id, int ridge_id_2, int opp_face_localid, s_ncell **OUT_PTRS);
int can_perform_flip44(const s_setup *setup, const s_ncell *ncell, double **vertices_face, int opp_cell_id, int *ridge_id_2);
void flip44(s_setup *setup, s_stack *stack, s_stack *stack_blocked, s_ncell *ncell, int id_ridge_1, int id_ridge_2, s_ncell **OUT_PTRS);


int pd_intersects_edge(double *s1, double *s2, double *p, double *d, int drop_coord);
int is_case_3(double **vertices_face, double *p, double *d);
int is_case_1(double **vertices_face, double *p, double *d);
int is_case_2(double **vertices_face, double *p, double *d);
int is_case_4(double **vertices_face, double *p, double *d);
int determine_case(double **vertices_face, double *p, double *d);


s_setup *initialize_setup(double **points, int N_points, int dim);
int flip_tetrahedra(s_setup *setup, s_stack *stack, s_stack *stack_blocked, s_ncell *ncell, int opp_cell_id);
void remove_point_setup(s_setup *setup, int point_id);
int insert_one_point(s_setup *setup, int point_id, s_stack *stack, s_stack *stack_blocked);
void remove_big_tetra(s_setup *setup);
int count_valid_ncells_reduced_triangulation(const s_setup *setup);


s_setup *construct_dt_3d(double **points, int N_points);

#endif
