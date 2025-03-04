#ifndef SIMPLICAL_COMPLEX_C
#define SIMPLICAL_COMPLEX_C

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "geometry.c"
#include "array_operations.c"


typedef struct setup {
    int dim;
    int N_points;
    double **points;  // (n_points + dim + 1) x dim
                      // The last (dim+1) points correspond to the big n_cell
    double *CM;  // Center of Mass
    int N_ncells;
    struct ncell *head;  // Linked list of ncells
} s_setup;


typedef struct ncell {
    int *vertex_id;
    struct ncell **opposite;
    struct ncell *next;  // Linked list of cells
    struct ncell *prev;
    int mark;  // Used in flood-fill algorithms, i.e. to find the star of an (n-3)-cell.
} s_ncell;


s_ncell *malloc_ncell(const s_setup *setup)
{
    s_ncell *out = malloc(sizeof(s_ncell));
    out->vertex_id = malloc(sizeof(int) * (setup->dim + 1));
    out->opposite = malloc(sizeof(s_ncell*) * (setup->dim + 1));
    out->next = NULL;
    out->prev = NULL;
    out->mark = 0;
    assert(out->vertex_id && out->opposite && "Could not malloc ncell!");
    return out;
}


void free_ncell(s_ncell *ncell)
{
    free(ncell->vertex_id);
    free(ncell->opposite);
    free(ncell);
}


void print_ncell(const s_setup *setup, const s_ncell *ncell)
{
    printf("%p, ( ", (void*)ncell);
    for (int ii=0; ii<setup->dim+1; ii++) {
        printf("%d ", ncell->vertex_id[ii]);
    }
    printf(") \n");
}


void print_ncells(const s_setup *setup)
{
    puts("NCELLS");
    s_ncell *current = setup->head;
    int ii = 0;
    while (current) {
        printf("%d  |  p : %p  |  prev : %p  |  marked : %d  |  vertex_ids :", ii, (void*)current, (void*)current->prev, current->mark);
        for (int jj=0; jj<setup->dim+1; jj++) {
            printf(" %d", current->vertex_id[jj]);
        }
        printf("  |  opposite :");
        for (int jj=0; jj<setup->dim+1; jj++) {
            printf(" %p", (void*)current->opposite[jj]);
        }
        printf("\n");
        ii++;
        current = current->next;
    }
    puts("");
}


void initialize_ncells_mark(const s_setup *setup)  // TODO CHECK IF THIS WORKS?
{  
    s_ncell *current = setup->head;
    for (int ii=0; ii<setup->N_ncells; ii++) {
        current->mark = 0;
        current = current->next;
    }
}


void print_marked(const s_setup *setup)
{
    puts("Marked node ids:");
    s_ncell *current = setup->head;
    for (int ii=0; ii<setup->N_ncells; ii++) {
        if (current->mark == 1) {
            printf("    %d\n", ii);
        }
        current = current->next;
    }
}


int count_marked(const s_setup *setup) 
{   
    int count = 0;
    s_ncell *current = setup->head;
    for (int ii=0; ii<setup->N_ncells; ii++) {
        if (current->mark == 1) {
            count++;
        }
        current = current->next;
    }
    return count;
}


void extract_vertices_ncell(const s_setup *setup, const s_ncell *ncell, double **out)
{
    for (int ii=0; ii<setup->dim+1; ii++) {
        for (int jj=0; jj<setup->dim; jj++) {
            out[ii][jj] = setup->points[ncell->vertex_id[ii]][jj];
        }
    }
}


void extract_ids_face(const s_setup *setup, const s_ncell *ncell, const int *v_localid, int dim_face, int *out)
{   // v_localid[dim-dim_face], out[dim_face+1]
    int kk = 0;
    for (int ii=0; ii<setup->dim+1; ii++) {
        if (!inarray(v_localid, setup->dim - dim_face, ii)) {
            out[kk] = ncell->vertex_id[ii];
            kk++;
        }
    }
}


void extract_vertices_face(const s_setup *setup, const s_ncell *ncell, const int *v_localid, int dim_face, double **out)
{   // v_localid[dim-dim_face], out[dim_face+1][dim]
    int face_vertex_id[dim_face+1];
    extract_ids_face(setup, ncell, v_localid, dim_face, face_vertex_id);
    for (int ii=0; ii<dim_face+1; ii++) {
        for (int jj=0; jj<setup->dim; jj++) {
            out[ii][jj] = setup->points[face_vertex_id[ii]][jj];
        }
    }
}


void face_localid_of_adjacent_ncell(const s_setup *setup, const s_ncell *ncell, const int *v_localid,
                                         int dim_face, int id_adjacent, int *out_v_localid)
{
    s_ncell *adjacent = ncell->opposite[id_adjacent];

    int vertex_id_face[dim_face+1];
    extract_ids_face(setup, ncell, v_localid, dim_face, vertex_id_face);

    int kk = 0;
    for (int ii=0; ii<setup->dim+1; ii++) {
        if (!inarray(vertex_id_face, dim_face+1, adjacent->vertex_id[ii])) {
            out_v_localid[kk] = ii;
            kk++;
        }
    }
}


s_ncell *next_ncell_ridge_cycle(const s_setup *setup, const s_ncell *ncell, int v_localid_main, int v_localid_2, 
                                int *new_v_localid_main, int *new_v_localid_2)
{
    s_ncell *next = ncell->opposite[v_localid_main];
    *new_v_localid_main = id_where_equal_int(next->vertex_id, setup->dim+1, ncell->vertex_id[v_localid_2]);
    for (int ii=0; ii<setup->dim+1; ii++) {
        if (next->opposite[ii] == ncell) {
            *new_v_localid_2 = ii;
            return next;
        }
    }
    assert(1 == 0 && "Could not find next ncell around ridge."); 
    exit(1);
}


int count_cycle_ridge(const s_setup *setup, const s_ncell *ncell, int v_localid_main, int v_localid_2)  // local 0:dim
{   
    int counter = 1;
    int maxit = 10000;
    const s_ncell *current = ncell;
    while (counter < maxit) {
        int new_v_localid_main, new_v_localid_2;
        if (!current->opposite[v_localid_main]) return 0;  // THIS IS TO CHECK THERE IS INDEED A NEXT CELL; TODO HANDLE THIS BETTER??

        s_ncell *next = next_ncell_ridge_cycle(setup, current, v_localid_main, v_localid_2, 
                                              &new_v_localid_main, &new_v_localid_2);
        if (next == ncell) return counter;
        counter++;
        current = next;
        v_localid_main = new_v_localid_main;
        v_localid_2 = new_v_localid_2;
    }
    assert(1==0 && "Reached maximum iterations in ridge cycle?");
    exit(1);
}


void mark_ncells_incident_face_STEP(const s_setup *setup, const s_ncell *ncell, const int *v_localid, int dim_face)
{
    for (int ii=0; ii<setup->dim+1; ii++) {
        if (inarray(v_localid, setup->dim - dim_face, ii)) {
            s_ncell *adjacent_ncell = ncell->opposite[ii];  // This ncell should already share the face!
            if (adjacent_ncell && adjacent_ncell->mark == 0) {
                adjacent_ncell->mark = 1;
                
                int new_v_localid[setup->dim - dim_face];
                face_localid_of_adjacent_ncell(setup, ncell, v_localid, dim_face, ii, new_v_localid);

                // Recursion
                mark_ncells_incident_face_STEP(setup, adjacent_ncell, new_v_localid, dim_face);
            }
        }
    }
}


void mark_ncells_incident_face(const s_setup *setup, s_ncell *ncell, const int *v_localid, int dim_face)
{
    initialize_ncells_mark(setup);
    ncell->mark = 1;
    
    // Recursion:
    mark_ncells_incident_face_STEP(setup, ncell, v_localid, dim_face);
}


void find_center_mass(double **in, int N_points, int dim, double *out)
{
    for (int ii=0; ii<dim; ii++) {
        out[ii] = in[0][ii];
    }

    for (int ii=1; ii<N_points; ii++) {
        for (int jj=0; jj<dim; jj++) {
            out[jj] += in[ii][jj];
        }
    }

    for (int ii=0; ii<dim; ii++) {
        out[ii] /= N_points;
    }
}


int are_locally_delaunay(const s_setup *setup, const s_ncell *ncell, int id_opposite)
{   
    // Create array for coords1 in static memory, CANNOT BE MULTI-THREADED! FIXME
    static int prev_dim = 0;
    static double **coords1 = NULL, **coords2 = NULL;
    if (setup->dim != prev_dim) {
        if (coords1) free_matrix(coords1, prev_dim + 1);
        if (coords2) free_matrix(coords2, prev_dim + 1);
        prev_dim = setup->dim;
        coords1 = malloc_matrix(setup->dim + 1, setup->dim);
        coords2 = malloc_matrix(setup->dim + 1, setup->dim);
    }

    extract_vertices_ncell(setup, ncell, coords1);
    extract_vertices_ncell(setup, ncell->opposite[id_opposite], coords2);
            
    // Extract vertex_id of opposite's cell face
    int opp_face_localid;
    face_localid_of_adjacent_ncell(setup, ncell, &id_opposite, setup->dim-1, id_opposite, &opp_face_localid);
    int opp_face_vertex_id = (ncell->opposite[id_opposite])->vertex_id[opp_face_localid];

    if (in_sphere(coords1, setup->points[opp_face_vertex_id], setup->dim) == -1 &&
        in_sphere(coords2, setup->points[ncell->vertex_id[id_opposite]], setup->dim) == -1) {
        return 1;
    } else {
        return 0;
    }
}


s_ncell *in_ncell_walk(s_setup *setup, double *p)  // Should make sure that p is inside the convull of all points (inside an n-cell)
{
    s_ncell *current = setup->head;
    assert(setup->N_ncells >= 1 && "N_ncells < 1");
    int randi = (rand() % setup->N_ncells);
    for (int ii=0; ii<randi; ii++) {  // Select random ncell to start
        current = current->next;
    }

    // Create array for facet_vertices in static memory, CANNOT BE MULTI-THREADED! FIXME
    static int prev_dim = 0;
    static double **facet_vertices = NULL;
    if (setup->dim != prev_dim) {
        if (facet_vertices) free_matrix(facet_vertices, prev_dim);
        prev_dim = setup->dim;
        facet_vertices = malloc_matrix(setup->dim, setup->dim);
    }

    STEP:
    for (int ii=0; ii<setup->dim+1; ii++) {
        double *opposite_vertex = setup->points[current->vertex_id[ii]];

        s_ncell *next = current->opposite[ii];
        if (next) {
            extract_vertices_face(setup, current, &ii, setup->dim-1, facet_vertices);

            int o1 = orientation(facet_vertices, opposite_vertex, setup->dim);
            int o2 = orientation(facet_vertices, p, setup->dim);
            assert(o1 != 0 && o2 != 0 && "Could not walk. p inside a facet?");
            if (o1 != o2) {
                current = next;
                goto STEP;
            }
        }
    }
    
    return current;
}


int is_delaunay_3d(const s_setup *setup)
{
    static double **vertices_ncell = NULL;
    if (!vertices_ncell) vertices_ncell = malloc_matrix(4, 3);

    s_ncell *current = setup->head;
    while (current) {
        extract_vertices_ncell(setup, current, vertices_ncell);
        s_ncell *query = current->next;
        while (query) {
            for (int ii=0; ii<4; ii++) {
                if (!inarray(current->vertex_id, 4, query->vertex_id[ii]) && 
                    in_sphere(vertices_ncell, setup->points[query->vertex_id[ii]], 3) != -1) {
                    printf("CONFLICT: (%d, %d, %d, %d) and %d\n", current->vertex_id[0], current->vertex_id[1], 
                                           current->vertex_id[2], current->vertex_id[3], query->vertex_id[ii]);
                    return 0;
                }
            }
            query = query->next;
        }
        current = current->next;
    }
    return 1;
}


// ----------------------------------------------------------------------------------------------
// --------------------------------------- PLOTS ------------------------------------------------
// ----------------------------------------------------------------------------------------------

void plot_ncell_3d(s_setup *setup, s_ncell *ncell, char *f_name, double *ranges)
{
    FILE *pipe = popen("gnuplot -persistent 2>&1", "w");
    fprintf(pipe, "set terminal pngcairo enhanced font 'Arial,18' size 1080,1080 enhanced \n");
    fprintf(pipe, "set output '%s.png'\n", f_name);
    fprintf(pipe, "set pm3d depthorder\n");
    fprintf(pipe, "set view 100, 60, \n");
    fprintf(pipe, "set xyplane at 0\n");
    fprintf(pipe, "set xrange [%f:%f]\n", ranges[0], ranges[1]);
    fprintf(pipe, "set yrange [%f:%f]\n", ranges[2], ranges[3]);
    fprintf(pipe, "set zrange [%f:%f]\n", ranges[4], ranges[5]);
    fprintf(pipe, "set xlabel 'x'\n");
    fprintf(pipe, "set ylabel 'y'\n");
    fprintf(pipe, "set zlabel 'z'\n");
    fflush(pipe);
    fprintf(pipe, "splot ");

    static double **face_vertices = NULL;
    if (!face_vertices) face_vertices = malloc_matrix(3, 3);
    for (int ii=0; ii<4; ii++) {
        extract_vertices_face(setup, ncell, &ii, 2, face_vertices);
        fprintf(pipe, "\"<echo \'");
        fprintf(pipe, "%f %f %f\\n", face_vertices[0][0], face_vertices[0][1], face_vertices[0][2]);
        fprintf(pipe, "%f %f %f\\n", face_vertices[1][0], face_vertices[1][1], face_vertices[1][2]);
        fprintf(pipe, "%f %f %f'\"", face_vertices[2][0], face_vertices[2][1], face_vertices[2][2]);
        fprintf(pipe, "w polygons fs transparent solid 0.1 notitle, ");
    }
    fprintf(pipe, "\"<echo \'");
    for (int ii=0; ii<setup->N_points; ii++) {
        fprintf(pipe, "%f %f %f\\n", setup->points[ii][0], setup->points[ii][1], setup->points[ii][2]);
    }
    fprintf(pipe, "'\" pt 7 lc rgb 'black' notitle, ");


    fflush(pipe);
    fprintf(pipe, "\n");
    pclose(pipe);
}


void plot_dt_3d(s_setup *setup, char *f_name, double *ranges)
{
    FILE *pipe = popen("gnuplot -persistent 2>&1", "w");
    fprintf(pipe, "set terminal pngcairo enhanced font 'Arial,18' size 1080,1080 enhanced \n");
    fprintf(pipe, "set output '%s.png'\n", f_name);
    fprintf(pipe, "set pm3d depthorder\n");
    fprintf(pipe, "set view 100, 60, \n");
    fprintf(pipe, "set xyplane at 0\n");
    fprintf(pipe, "set xrange [%f:%f]\n", ranges[0], ranges[1]);
    fprintf(pipe, "set yrange [%f:%f]\n", ranges[2], ranges[3]);
    fprintf(pipe, "set zrange [%f:%f]\n", ranges[4], ranges[5]);
    fprintf(pipe, "set xlabel 'x'\n");
    fprintf(pipe, "set ylabel 'y'\n");
    fprintf(pipe, "set zlabel 'z'\n");
    fflush(pipe);
    fprintf(pipe, "splot ");
    static double **face_vertices = NULL;
    if (!face_vertices) face_vertices = malloc_matrix(3, 3);

    s_ncell *current = setup->head;
    while (current) {
        for (int ii=0; ii<4; ii++) {
            extract_vertices_face(setup, current, &ii, 2, face_vertices);
            fprintf(pipe, "\"<echo \'");
            fprintf(pipe, "%f %f %f\\n", face_vertices[0][0], face_vertices[0][1], face_vertices[0][2]);
            fprintf(pipe, "%f %f %f\\n", face_vertices[1][0], face_vertices[1][1], face_vertices[1][2]);
            fprintf(pipe, "%f %f %f'\"", face_vertices[2][0], face_vertices[2][1], face_vertices[2][2]);
            fprintf(pipe, "w polygons fs transparent solid 0.1 notitle, ");
        }
        for (int ii=0; ii<4; ii++) {
            extract_vertices_face(setup, current, &ii, 2, face_vertices);
            fprintf(pipe, "\"<echo \'");
            fprintf(pipe, "%f %f %f\\n", face_vertices[0][0], face_vertices[0][1], face_vertices[0][2]);
            fprintf(pipe, "%f %f %f\\n", face_vertices[1][0], face_vertices[1][1], face_vertices[1][2]);
            fprintf(pipe, "%f %f %f'\"", face_vertices[2][0], face_vertices[2][1], face_vertices[2][2]);
            fprintf(pipe, "pt 7 lc rgb 'black' notitle, ");
        }
        current = current->next;
    }

    fflush(pipe);
    fprintf(pipe, "\n");
    pclose(pipe);
}

#endif
