#ifndef VOR3D_MIRRORING_H
#define VOR3D_MIRRORING_H

#include "points.h"
#include "bpoly.h"

typedef struct extruding_vertex {
    s_point vertex;
    int vcell;
} s_extruding_vertex;


typedef struct vertex_list {
    s_extruding_vertex *list;
    int N;
    int Nmax;
} s_vertex_list;

s_vertex_list initialize_vertex_list(int Nmax);
int increase_memory_vertex_list_if_needed(s_vertex_list *v_list, int N_needed);
void free_vertex_list(s_vertex_list *v_list);


int extend_sites_mirroring_initial(const s_bpoly *bp, double EPS_degenerate, double TOL, s_points *inout_seeds);
int extend_sites_mirroring_extruding(const s_bpoly *bp, double EPS_degenerate, double TOL, const s_vertex_list *extruding, int Nreal, s_points *inout_seeds);

#endif
