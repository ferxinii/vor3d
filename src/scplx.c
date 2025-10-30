
#include "scplx.h"
#include "geometry.h"
#include "algebra.h"
#include "gnuplotc.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>


s_ncell *malloc_ncell(void)
{
    s_ncell *out = malloc(sizeof(s_ncell));
    assert(out && "Could not malloc ncell");
    memset(out, 0, sizeof(s_ncell));
    return out;
}


void free_ncell(s_ncell *ncell)
{
    free(ncell);
}


void free_complex(s_scplx *setup)
{
    free_points(&setup->points);

    s_ncell *current = setup->head;
    while (current) {
        s_ncell *next = current->next;
        free_ncell(current);
        current = next;
    }

    memset(setup, 0, sizeof(s_scplx));
}


void print_ncell(const s_ncell *ncell)
{
    printf("%p, ( ", (void*)ncell);
    for (int ii=0; ii<4; ii++) {
        printf("%d ", ncell->vertex_id[ii]);
    }
    printf(") \n");
}


void print_scomplex(const s_scplx *setup)
{
    puts("NCELLS");
    s_ncell *current = setup->head;
    int ii = 0;
    while (current) {
        printf("%d  |  p : %p  |  prev : %p  |  marked : %d  |  vertex_ids :", ii, (void*)current, (void*)current->prev, current->mark);
        for (int jj=0; jj<4; jj++) {
            printf(" %d", current->vertex_id[jj]);
        }
        printf("  |  opposite :");
        for (int jj=0; jj<4; jj++) {
            printf(" %p", (void*)current->opposite[jj]);
        }
        printf("\n");
        ii++;
        current = current->next;
    }
    puts("");
}


void write_ncell3d_file(const s_scplx *setup, const s_ncell *ncell, FILE *file)
{
    fprintf(file, "%f %f %f\n", setup->points.p[ncell->vertex_id[0]].x,
                                setup->points.p[ncell->vertex_id[0]].y,
                                setup->points.p[ncell->vertex_id[0]].z);
    fprintf(file, "%f %f %f\n", setup->points.p[ncell->vertex_id[1]].x,
                                setup->points.p[ncell->vertex_id[1]].y,
                                setup->points.p[ncell->vertex_id[1]].z);
    fprintf(file, "%f %f %f\n\n", setup->points.p[ncell->vertex_id[2]].x,
                                  setup->points.p[ncell->vertex_id[2]].y,
                                  setup->points.p[ncell->vertex_id[2]].z);

    fprintf(file, "%f %f %f\n", setup->points.p[ncell->vertex_id[0]].x,
                                setup->points.p[ncell->vertex_id[0]].y,
                                setup->points.p[ncell->vertex_id[0]].z);
    fprintf(file, "%f %f %f\n", setup->points.p[ncell->vertex_id[1]].x,
                                setup->points.p[ncell->vertex_id[1]].y,
                                setup->points.p[ncell->vertex_id[1]].z);
    fprintf(file, "%f %f %f\n\n", setup->points.p[ncell->vertex_id[3]].x,
                                  setup->points.p[ncell->vertex_id[3]].y,
                                  setup->points.p[ncell->vertex_id[3]].z);

    fprintf(file, "%f %f %f\n", setup->points.p[ncell->vertex_id[3]].x,
                                setup->points.p[ncell->vertex_id[3]].y,
                                setup->points.p[ncell->vertex_id[3]].z);
    fprintf(file, "%f %f %f\n", setup->points.p[ncell->vertex_id[1]].x,
                                setup->points.p[ncell->vertex_id[1]].y,
                                setup->points.p[ncell->vertex_id[1]].z);
    fprintf(file, "%f %f %f\n\n", setup->points.p[ncell->vertex_id[2]].x,
                                  setup->points.p[ncell->vertex_id[2]].y,
                                  setup->points.p[ncell->vertex_id[2]].z);

    fprintf(file, "%f %f %f\n", setup->points.p[ncell->vertex_id[0]].x,
                                setup->points.p[ncell->vertex_id[0]].y,
                                setup->points.p[ncell->vertex_id[0]].z);
    fprintf(file, "%f %f %f\n", setup->points.p[ncell->vertex_id[3]].x,
                                setup->points.p[ncell->vertex_id[3]].y,
                                setup->points.p[ncell->vertex_id[3]].z);
    fprintf(file, "%f %f %f\n\n\n", setup->points.p[ncell->vertex_id[2]].x,
                                    setup->points.p[ncell->vertex_id[2]].y,
                                    setup->points.p[ncell->vertex_id[2]].z);
}


void write_scomplex_file(const s_scplx *setup, FILE *file)
{
    s_ncell *current = setup->head;
    while (current) {
        write_ncell3d_file(setup, current, file);
        fprintf(file, "\n\n");
        current = current->next;
    }
}


void initialize_ncells_counter(const s_scplx *setup)
{
    s_ncell *current = setup->head;
    int ii = 0;
    while (current) {
        current->count = ii;
        current = current->next;
        ii++;
    }
}


void initialize_ncells_mark(const s_scplx *setup)
{  
    s_ncell *current = setup->head;
    for (int ii=0; ii<setup->N_ncells; ii++) {
        current->mark = 0;
        current = current->next;
    }
}


void print_marked(const s_scplx *setup)
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


int count_marked(const s_scplx *setup) 
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


void extract_vertices_ncell(const s_scplx *setup, const s_ncell *ncell, s_point out[4])
{
    for (int ii=0; ii<4; ii++) {
        out[ii] = setup->points.p[ncell->vertex_id[ii]];
    }
}


void extract_ids_face(const s_ncell *ncell, int dim_face, const int *v_localid, int *out)
{   // v_localid[dim-dim_face], out[dim_face+1]
    if (dim_face < 0 || dim_face >= 3) {
        fprintf(stderr, "dim_face must be >= 0 && < 3\n");
        exit(1);
    }

    int kk = 0;
    for (int ii=0; ii<4; ii++) {
        if (!inarray(v_localid, 3 - dim_face, ii)) {
            out[kk] = ncell->vertex_id[ii];
            kk++;
        }
    }
}


void extract_vertices_face(const s_scplx *setup, const s_ncell *ncell, int dim_face, const int v_localid[3-dim_face], s_point out[dim_face+1])
{   // v_localid[dim-dim_face], out[dim_face+1]
    if (dim_face < 0 || dim_face >= 3) {
        fprintf(stderr, "dim_face must be >= 0 && < 3\n");
        exit(1);
    }

    int face_vertex_id[dim_face+1];
    extract_ids_face(ncell, dim_face, v_localid, face_vertex_id);
    for (int ii=0; ii<dim_face+1; ii++) {
        out[ii] = setup->points.p[face_vertex_id[ii]];
    }

}


void face_localid_of_adjacent_ncell(const s_ncell *ncell, int dim_face, const int v_localid[3-dim_face], int id_adjacent, int out_v_localid[3-dim_face])
{
    if (dim_face < 0 || dim_face >= 3) {
        fprintf(stderr, "dim_face must be >= 0 && < 3\n");
        exit(1);
    }

    s_ncell *adjacent = ncell->opposite[id_adjacent];

    int vertex_id_face[dim_face+1];
    extract_ids_face(ncell, dim_face, v_localid, vertex_id_face);

    int kk = 0;
    for (int ii=0; ii<4; ii++) {
        if (!inarray(vertex_id_face, dim_face+1, adjacent->vertex_id[ii])) {
            out_v_localid[kk] = ii;
            kk++;
        }
    }
}


s_ncell *next_ncell_ridge_cycle(const s_ncell *ncell, int v_localid_main, int v_localid_2, int *new_v_localid_main, int *new_v_localid_2)
{
    s_ncell *next = ncell->opposite[v_localid_main];
    *new_v_localid_main = id_where_equal_int(next->vertex_id, 4, ncell->vertex_id[v_localid_2]);
    for (int ii=0; ii<4; ii++) {
        if (next->opposite[ii] == ncell) {
            *new_v_localid_2 = ii;
            return next;
        }
    }
    fprintf(stderr, "Could not find next ncell around ridge.\n"); 
    exit(1);
}


void mark_ncells_incident_face_STEP(const s_ncell *ncell, int dim_face, const int v_localid[3-dim_face])
{
    for (int ii=0; ii<4; ii++) {
        if (inarray(v_localid, 3-dim_face, ii)) {
            s_ncell *adjacent_ncell = ncell->opposite[ii];  // This ncell should already share the face!
            if (adjacent_ncell && adjacent_ncell->mark == 0) {
                adjacent_ncell->mark = 1;
                
                int new_v_localid[3-dim_face];
                face_localid_of_adjacent_ncell(ncell, dim_face, v_localid, ii, new_v_localid);

                // Recursion
                mark_ncells_incident_face_STEP(adjacent_ncell, dim_face, new_v_localid);
            }
        }
    }
}


void mark_ncells_incident_face(const s_scplx *setup, s_ncell *ncell, int dim_face, const int v_localid[3-dim_face])
{
    if (dim_face < 0 || dim_face >= 3) {
        fprintf(stderr, "dim_face must be >= 0 && < 3\n");
        exit(1);
    }

    initialize_ncells_mark(setup);
    ncell->mark = 1;
    
    // Recursion:
    mark_ncells_incident_face_STEP(ncell, dim_face, v_localid);
}


int are_locally_delaunay_strict(const s_scplx *setup, const s_ncell *ncell, int id_opposite)
{   // I.E., only return true if the point is INSIDE circumscr., not on it.
    s_point coords1[4];
    s_point coords2[4];
    extract_vertices_ncell(setup, ncell, coords1);
    extract_vertices_ncell(setup, ncell->opposite[id_opposite], coords2);

    assert(!( orientation(coords1, coords1[3]) == 0 && 
              orientation(coords2, coords2[3]) == 0 ) && 
           "Cannot check delaunayness, both ncells are degenerate.");

    if (orientation(coords1, coords1[3]) == 0 || 
        orientation(coords2, coords2[3]) == 0) return 0;

    // Extract vertex_id of opposite's cell face
    int opp_face_localid;
    face_localid_of_adjacent_ncell(ncell, 2, &id_opposite, id_opposite, &opp_face_localid);
    int opp_face_vertex_id = (ncell->opposite[id_opposite])->vertex_id[opp_face_localid];
    
    int in1 = in_sphere(coords1, setup->points.p[opp_face_vertex_id]);
    if (in1 == -1) return 1;
    else return 0;
}


int are_locally_delaunay_nonstrict(const s_scplx *setup, const s_ncell *ncell, int id_opposite)
{   // Points on circumscribed sphere are VALID
    s_point coords1[4];
    s_point coords2[4];
    extract_vertices_ncell(setup, ncell, coords1);
    extract_vertices_ncell(setup, ncell->opposite[id_opposite], coords2);

    assert(!(orientation(coords1, coords1[3]) == 0 && 
             orientation(coords2, coords2[3]) == 0) &&
           "Cannot check delaunayness, both ncells are degenerate.");

    if (orientation(coords1, coords1[3]) == 0 || 
        orientation(coords2, coords2[3]) == 0) return 0;

    // Extract vertex_id of opposite's cell face
    int opp_face_localid;
    face_localid_of_adjacent_ncell(ncell, 2, &id_opposite, id_opposite, &opp_face_localid);
    int opp_face_vertex_id = (ncell->opposite[id_opposite])->vertex_id[opp_face_localid];
    
    int in1 = in_sphere(coords1, setup->points.p[opp_face_vertex_id]);
    if (in1 != 1) return 1;
    else return 0;
}


int in_ncell(const s_scplx *setup, const s_ncell *ncell, s_point query)
{   // TODO Assumes consisten ordering of vertices?
    s_point v0 = setup->points.p[ncell->vertex_id[0]];
    s_point v1 = setup->points.p[ncell->vertex_id[1]];
    s_point v2 = setup->points.p[ncell->vertex_id[2]];
    s_point v3 = setup->points.p[ncell->vertex_id[3]];
    
    s_point tetra[4] = {v0, v1, v2, v3};
    return in_tetrahedron(tetra, query);
 }


s_ncell *bruteforce_find_ncell_containing(const s_scplx *setup, s_point p)
{
    s_ncell *current = setup->head;
    while (current) {
        if (in_ncell(setup, current, p)) return current;
        current = current->next;
    }
    fprintf(stderr, "did not find container ncell.\n");
    exit(1);
}


s_ncell *in_ncell_walk(const s_scplx *setup, s_point p)
{
    s_ncell *current = setup->head;
    assert(setup->N_ncells >= 1 && "N_ncells < 1");
    int randi = (rand() % setup->N_ncells);
    for (int ii=0; ii<randi; ii++) {  // Select random ncell to start
        current = current->next;
    }

    s_point facet_vertices[3];
    s_ncell *prev = current;
    STEP:
    for (int ii=0; ii<4; ii++) {
        s_point opposite_vertex = setup->points.p[current->vertex_id[ii]];

        s_ncell *next = current->opposite[ii];
        if (next) {
            extract_vertices_face(setup, current, 2, &ii, facet_vertices);

            int o1 = orientation(facet_vertices, opposite_vertex);
            int o2 = orientation(facet_vertices, p);

            if (o1 == 0 && next != prev) {
                prev = current;
                current = next;
                goto STEP;
            } else if (o1 == 0) continue;

            if (o2 == 0) {
                if (in_ncell(setup, current, p)) { return current; }
                else if (next != prev) {
                    prev = current;
                    current = next;
                    goto STEP;
                } else continue;
            } else if (o1 != o2) {
                prev = current;
                current = next;
                goto STEP;
            }
        }
    }
    
    return current;
}


int is_delaunay_3d(const s_scplx *setup)
{
    s_point vertices_ncell[4];
    s_ncell *current = setup->head;
    while (current) {
        extract_vertices_ncell(setup, current, vertices_ncell);
        if (orientation(vertices_ncell, vertices_ncell[3]) == 0) {
            puts("Flat tetra!!");
            return 0;
        }
        for (int ii=0; ii<4; ii++) {
            if (current->opposite[ii] &&
                !are_locally_delaunay_nonstrict(setup, current, ii)) return 0;
        }
        current = current->next;
    }
    return 1;
}





// -----------------------------------------------------------------------------------------
// --------------------------------------- PLOTS -------------------------------------------
// -----------------------------------------------------------------------------------------


void plot_add_ncell(s_gnuplot *interface, const s_scplx *setup, const s_ncell *ncell, char *config)
{
    s_point face_vertices[3];
    for (int ii=0; ii<4; ii++) {
        extract_vertices_face(setup, ncell, 2, &ii, face_vertices);
        draw_solid_triangle_3d(interface, face_vertices[0].coords, face_vertices[1].coords,
                               face_vertices[2].coords, config);
    }
}


void plot_ncell_3d(const s_scplx *setup, const s_ncell *ncell, char *f_name, s_point ranges[2])
{
    s_gnuplot *interface = gnuplot_start(PNG_3D, f_name, (int[2]){1080, 1080}, 18);
    gnuplot_config(interface, "set pm3d depthorder",
                              "set pm3d border lc 'black' lw 0.5",
                              "set view 100, 60",
                              "set xyplane at 0");
    if (ranges) {
        char buff[1024];
        snprintf(buff, 1024, "set xrange [%f:%f]\n set yrange [%f:%f]\n set zrange [%f:%f]", 
                 ranges[0].x, ranges[1].x, ranges[0].y, ranges[1].y, ranges[0].z, ranges[1].z); 
        gnuplot_config(interface, buff);
    }

    plot_add_ncell(interface, setup, ncell, "fs transparent solid 0.2 fc rgb '#000090'");
    gnuplot_end(interface);
}


void plot_all_ncells_3d(s_scplx *setup, char *f_name, s_point ranges[2], char *view_command)
{
     char colors[][20] = { "#000090", "#000fff", "#0090ff", "#0fffee", 
        "#90ff70", "#ffee00", "#ff7000", "#ee0000", "#7f0000" };
    char buff[1024];

    s_gnuplot *interface = gnuplot_start(PNG_3D, f_name, (int[2]){1080, 1080}, 18);
    gnuplot_config(interface, "set pm3d depthorder",
                              "set pm3d border lc 'black' lw 0.5",
                              "set xyplane at 0",
                              view_command);
    if (ranges) {
        snprintf(buff, 1024, "set xrange [%f:%f]\n set yrange [%f:%f]\n set zrange [%f:%f]", 
                 ranges[0].x, ranges[1].x, ranges[0].y, ranges[1].y, ranges[0].z, ranges[1].z); 
    }
    if (ranges) gnuplot_config(interface, buff);

    int it = 0;
    s_ncell *current = setup->head;
    while (current) {
        snprintf(buff, 1024, "fs transparent solid 0.2 fc rgb '%s'", colors[it%8]);
        plot_add_ncell(interface, setup, current, buff);
        current = current->next;
        it++;
    }
    gnuplot_end(interface);
}


void plot_dt_differentviews(s_scplx *setup, char *f_name, s_point ranges[2])
{
    char final_name[512];
    snprintf(final_name, 512, "%s_v1.png", f_name);
    plot_all_ncells_3d(setup, final_name, ranges, "set view 100, 60, 1.5");

    snprintf(final_name, 512, "%s_v2.png", f_name);
    plot_all_ncells_3d(setup, final_name, ranges, "set view 100, 90, 1.5");

    snprintf(final_name, 512, "%s_v3.png", f_name);
    plot_all_ncells_3d(setup, final_name, ranges, "set view 100, 180, 1.5");

    snprintf(final_name, 512, "%s_v4.png", f_name);
    plot_all_ncells_3d(setup, final_name, ranges, "set view 100, 270, 1.5");

}

