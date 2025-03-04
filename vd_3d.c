
#include "simplical_complex.c"

typedef struct poly {
    int Nv;
    int *vertex_id;
    struct poly **opposite;
    struct poly *next;  // Linked list of cells
    struct poly *prev;
    int mark;  // Used to mark particular ncells
} s_poly;


s_poly *malloc_poly(const s_setup *setup, int N_vertices)
{
    s_poly *out = malloc(sizeof(s_poly));
    out->Nv = N_vertices;
    out->vertex_id = malloc(sizeof(int) * N_vertices);
    out->opposite = malloc(sizeof(s_poly*) * N_vertices);
    out->next = NULL;
    out->prev = NULL;
    out->mark = 0;
    return out;
}


void extract_vertex(const s_setup *setup, const s_ncell *ncell, double *out)
{
    static double **vertices_ncell = NULL;
    if (!vertices_ncell) vertices_ncell = malloc_matrix(4, 3);
    
    extract_vertices_ncell(setup, ncell, vertices_ncell);

    find_center_mass(vertices_ncell, 4, 3, out);
}


void extract_polyhedron
