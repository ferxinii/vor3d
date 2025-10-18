
#include "simplical_complex.h"
#include "geometry.h"
#include "algebra.h"
#include "convh.h"
#include "gnuplotc.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>


s_ncell *malloc_ncell(const s_scplx *setup)
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


void print_ncell(const s_scplx *setup, const s_ncell *ncell)
{
    printf("%p, ( ", (void*)ncell);
    for (int ii=0; ii<setup->dim+1; ii++) {
        printf("%d ", ncell->vertex_id[ii]);
    }
    printf(") \n");
}


void print_ncells(const s_scplx *setup)
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


void write_ncell3d_file(s_scplx *setup, s_ncell *ncell, FILE *file)
{
    fprintf(file, "%f %f %f\n", setup->points[ncell->vertex_id[0]].x,
                                setup->points[ncell->vertex_id[0]].y,
                                setup->points[ncell->vertex_id[0]].z);
    fprintf(file, "%f %f %f\n", setup->points[ncell->vertex_id[1]].x,
                                setup->points[ncell->vertex_id[1]].y,
                                setup->points[ncell->vertex_id[1]].z);
    fprintf(file, "%f %f %f\n\n", setup->points[ncell->vertex_id[2]].x,
                                  setup->points[ncell->vertex_id[2]].y,
                                  setup->points[ncell->vertex_id[2]].z);

    fprintf(file, "%f %f %f\n", setup->points[ncell->vertex_id[0]].x,
                                setup->points[ncell->vertex_id[0]].y,
                                setup->points[ncell->vertex_id[0]].z);
    fprintf(file, "%f %f %f\n", setup->points[ncell->vertex_id[1]].x,
                                setup->points[ncell->vertex_id[1]].y,
                                setup->points[ncell->vertex_id[1]].z);
    fprintf(file, "%f %f %f\n\n", setup->points[ncell->vertex_id[3]].x,
                                  setup->points[ncell->vertex_id[3]].y,
                                  setup->points[ncell->vertex_id[3]].z);

    fprintf(file, "%f %f %f\n", setup->points[ncell->vertex_id[3]].x,
                                setup->points[ncell->vertex_id[3]].y,
                                setup->points[ncell->vertex_id[3]].z);
    fprintf(file, "%f %f %f\n", setup->points[ncell->vertex_id[1]].x,
                                setup->points[ncell->vertex_id[1]].y,
                                setup->points[ncell->vertex_id[1]].z);
    fprintf(file, "%f %f %f\n\n", setup->points[ncell->vertex_id[2]].x,
                                  setup->points[ncell->vertex_id[2]].y,
                                  setup->points[ncell->vertex_id[2]].z);

    fprintf(file, "%f %f %f\n", setup->points[ncell->vertex_id[0]].x,
                                setup->points[ncell->vertex_id[0]].y,
                                setup->points[ncell->vertex_id[0]].z);
    fprintf(file, "%f %f %f\n", setup->points[ncell->vertex_id[3]].x,
                                setup->points[ncell->vertex_id[3]].y,
                                setup->points[ncell->vertex_id[3]].z);
    fprintf(file, "%f %f %f\n\n\n", setup->points[ncell->vertex_id[2]].x,
                                    setup->points[ncell->vertex_id[2]].y,
                                    setup->points[ncell->vertex_id[2]].z);
}


void write_dt3d_file(s_scplx *setup, FILE *file)
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


void initialize_ncells_mark(const s_scplx *setup)  // TODO CHECK IF THIS WORKS?
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


void extract_vertices_ncell(const s_scplx *setup, const s_ncell *ncell, s_point *out)
{
    for (int ii=0; ii<setup->dim+1; ii++) {
        out[ii] = setup->points[ncell->vertex_id[ii]];
    }
}


void extract_ids_face(const s_scplx *setup, const s_ncell *ncell, const int *v_localid, int dim_face, int *out)
{   // v_localid[dim-dim_face], out[dim_face+1]
    int kk = 0;
    for (int ii=0; ii<setup->dim+1; ii++) {
        if (!inarray(v_localid, setup->dim - dim_face, ii)) {
            out[kk] = ncell->vertex_id[ii];
            kk++;
        }
    }
}


void extract_vertices_face(const s_scplx *setup, const s_ncell *ncell, const int *v_localid, int dim_face, s_point *out)
{   // v_localid[dim-dim_face], out[dim_face+1][dim]
    int face_vertex_id[dim_face+1];
    extract_ids_face(setup, ncell, v_localid, dim_face, face_vertex_id);
    for (int ii=0; ii<dim_face+1; ii++) {
        out[ii] = setup->points[face_vertex_id[ii]];
    }
}


void extract_face_center_and_normal(const s_scplx *setup, const s_ncell *ncell, int face_localid, s_point *fc, s_point *n)
{
    assert(setup->dim == 3 && "Currently only supports 3D");
    s_point face_v[3];
    s_point ncell_v[3];

    extract_vertices_face(setup, ncell, &face_localid, 2, face_v);
    extract_vertices_ncell(setup, ncell, ncell_v);

    s_point d1 = subtract_points(face_v[1], face_v[0]);
    s_point d2 = subtract_points(face_v[2], face_v[0]);
    *n = cross_prod(d1, d2);

    *fc = find_center_mass(face_v, 3);
    s_point cc = find_center_mass(ncell_v, 4);
    s_point v = subtract_points(*fc, cc);

    double dir_aux = dot_prod(v, *n);
    assert(dir_aux != 0 && "Vectors are perpendicular?");
    if (dir_aux < 0) {
        (*n).x = -(*n).x;
        (*n).y = -(*n).y;
        (*n).z = -(*n).z;
    } 
}


void face_localid_of_adjacent_ncell(const s_scplx *setup, const s_ncell *ncell, const int *v_localid,
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


s_ncell *next_ncell_ridge_cycle(const s_scplx *setup, const s_ncell *ncell, int v_localid_main, int v_localid_2, 
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


int count_cycle_ridge(const s_scplx *setup, const s_ncell *ncell, int v_localid_main, int v_localid_2)  // local 0:dim
{   
    int counter = 1;
    int maxit = 10000;
    const s_ncell *current = ncell;
    while (counter < maxit) {
        int new_v_localid_main, new_v_localid_2;
        if (current->opposite[v_localid_main] == NULL) {
            return -1;
        }

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


void mark_ncells_incident_face_STEP(const s_scplx *setup, const s_ncell *ncell, const int *v_localid, int dim_face)
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


void mark_ncells_incident_face(const s_scplx *setup, s_ncell *ncell, const int *v_localid, int dim_face)
{
    initialize_ncells_mark(setup);
    ncell->mark = 1;
    
    // Recursion:
    mark_ncells_incident_face_STEP(setup, ncell, v_localid, dim_face);
}


int are_locally_delaunay_strict(const s_scplx *setup, const s_ncell *ncell, int id_opposite)
{   // I.E., only return true if the point is INSIDE circumscr., not on it.
    assert(setup->dim == 3 && "Currently only supports 3D");
    s_point coords1[4];
    s_point coords2[4];
    extract_vertices_ncell(setup, ncell, coords1);
    extract_vertices_ncell(setup, ncell->opposite[id_opposite], coords2);

    assert(!(orientation(coords1, coords1[3]) == 0 && 
             orientation(coords2, coords2[3]) == 0));

    if (orientation(coords1, coords1[3]) == 0 || 
        orientation(coords2, coords2[3]) == 0) return 0;

    // Extract vertex_id of opposite's cell face
    int opp_face_localid;
    face_localid_of_adjacent_ncell(setup, ncell, &id_opposite, setup->dim-1, id_opposite, &opp_face_localid);
    int opp_face_vertex_id = (ncell->opposite[id_opposite])->vertex_id[opp_face_localid];
    
    int in1 = in_sphere(coords1, setup->points[opp_face_vertex_id]);
    if (in1 == -1) return 1;
    else return 0;
}


int are_locally_delaunay_nonstrict(const s_scplx *setup, const s_ncell *ncell, int id_opposite)
{   // Points on circumscribed sphere are VALID
    assert(setup->dim == 3 && "Currently only supports 3D");
    s_point coords1[4];
    s_point coords2[4];
    extract_vertices_ncell(setup, ncell, coords1);
    extract_vertices_ncell(setup, ncell->opposite[id_opposite], coords2);

    assert(!(orientation(coords1, coords1[3]) == 0 && 
             orientation(coords2, coords2[3]) == 0));

    if (orientation(coords1, coords1[3]) == 0 || 
        orientation(coords2, coords2[3]) == 0) return 0;

    // Extract vertex_id of opposite's cell face
    int opp_face_localid;
    face_localid_of_adjacent_ncell(setup, ncell, &id_opposite, setup->dim-1, id_opposite, &opp_face_localid);
    int opp_face_vertex_id = (ncell->opposite[id_opposite])->vertex_id[opp_face_localid];
    
    int in1 = in_sphere(coords1, setup->points[opp_face_vertex_id]);
    if (in1 != 1) return 1;
    else return 0;
}


int point_in_face(s_point vertices_face[3], s_point p)
{
    s_point d1 = subtract_points(vertices_face[1], vertices_face[0]);
    s_point d2 = subtract_points(vertices_face[2], vertices_face[0]);
    s_point n = cross_prod(d1, d2);

    int drop_coord = coord_with_largest_component_3d(n);

    if (orientation(vertices_face, p) != 0) return 0;

    int i1, i2;
    if (drop_coord == 0) {      i1 = 1;     i2 = 2; }
    else if (drop_coord == 1) { i1 = 2;     i2 = 0; }
    else {                      i1 = 0;     i2 = 1; }
    
    s_point v1 = {{{vertices_face[0].coords[i1], vertices_face[0].coords[i2], 0}}};
    s_point v2 = {{{vertices_face[1].coords[i1], vertices_face[1].coords[i2], 0}}};
    s_point v3 = {{{vertices_face[2].coords[i1], vertices_face[2].coords[i2], 0}}};
    s_point triangle[3] = {v1, v2, v3};
    s_point paux = {{{p.coords[i1], p.coords[i2], 0}}};

    return point_in_triangle_2d(triangle, paux);
}


int point_in_tetra(const s_scplx *setup, s_point x, const s_ncell *nc)
{   // TODO Assumes consisten ordering of vertices?
    assert(setup->dim == 3 && "Only supports 3D");
    s_point v0 = setup->points[nc->vertex_id[0]];
    s_point v1 = setup->points[nc->vertex_id[1]];
    s_point v2 = setup->points[nc->vertex_id[2]];
    s_point v3 = setup->points[nc->vertex_id[3]];

    // compute signed volumes (orientation)
    s_point facet_vertices[3];

    facet_vertices[0] = v1; facet_vertices[1] = v2; facet_vertices[2] = v3; 
    int s0 = orientation(facet_vertices, x);  // face opposite v0

    facet_vertices[0] = v0; facet_vertices[1] = v3; facet_vertices[2] = v2; 
    int s1 = orientation(facet_vertices, x);  // face opposite v0
                                                 //
    facet_vertices[0] = v0; facet_vertices[1] = v1; facet_vertices[2] = v3; 
    int s2 = orientation(facet_vertices, x);  // face opposite v0

    facet_vertices[0] = v0; facet_vertices[1] = v2; facet_vertices[2] = v1; 
    int s3 = orientation(facet_vertices, x);  // face opposite v0

    // find a nonzero reference sign
    int ref = 0;
    if (s0 != 0) ref = s0;
    else if (s1 != 0) ref = s1;
    else if (s2 != 0) ref = s2;
    else if (s3 != 0) ref = s3;
    if (ref == 0) {
        puts("DEGENERATE!");
        exit(1);
    }

    // if any nonzero sign disagrees, x is outside
    if ((s0 && s0 != ref) || (s1 && s1 != ref) ||
        (s2 && s2 != ref) || (s3 && s3 != ref)   )
        return 0;
    
    return 1;
}


s_ncell *bruteforce_find_ncell_containing(const s_scplx *setup, s_point p)
{
    s_ncell *current = setup->head;
    while (current) {
        if (point_in_tetra(setup, p, current)) return current;
        current = current->next;
    }
    puts("DID NOT FIND CONTAINER NCELL!");
    exit(1);
}


s_ncell *in_ncell_walk(const s_scplx *setup, s_point p)
{
    assert(setup->dim == 3 && "Only supports 3D");

    s_ncell *current = setup->head;
    assert(setup->N_ncells >= 1 && "N_ncells < 1");
    int randi = (rand() % setup->N_ncells);
    for (int ii=0; ii<randi; ii++) {  // Select random ncell to start
        current = current->next;
    }

    s_point facet_vertices[3];
    s_ncell *prev = current;
    STEP:
    for (int ii=0; ii<setup->dim+1; ii++) {
        s_point opposite_vertex = setup->points[current->vertex_id[ii]];

        s_ncell *next = current->opposite[ii];
        if (next) {
            extract_vertices_face(setup, current, &ii, setup->dim-1, facet_vertices);

            int o1 = orientation(facet_vertices, opposite_vertex);
            int o2 = orientation(facet_vertices, p);

            if (o1 == 0 && next != prev) {
                prev = current;
                current = next;
                goto STEP;
            } else if (o1 == 0) continue;

            if (o2 == 0) {
                if (point_in_tetra(setup, p, current)) { return current; }
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


void add_ncell_volume_3d(const s_scplx *setup, s_ncell *ncell)
{   // THIS IS JUST FOR DEBUGGING!!
    ncell->volume = volume_tetrahedron_approx(
                    setup->points[ncell->vertex_id[0]], 
                    setup->points[ncell->vertex_id[1]],
                    setup->points[ncell->vertex_id[2]],
                    setup->points[ncell->vertex_id[3]]);
}


double compute_volume_complex(s_scplx *setup)
{
    return compute_volume_convhull_from_points(setup->points, setup->N_points);
}




// -----------------------------------------------------------------------------------------
// --------------------------------------- PLOTS -------------------------------------------
// -----------------------------------------------------------------------------------------


void plot_add_ncell(s_gnuplot *interface, const s_scplx *setup, const s_ncell *ncell, char *config)
{
    s_point face_vertices[3];
    for (int ii=0; ii<4; ii++) {
        extract_vertices_face(setup, ncell, &ii, 2, face_vertices);
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


    // A NEW PLOT FOR EACH CELL (UNTIL MAX_FILES)
    // int it = 0;
    // current = setup->head;
    // while (current) {
    //     if (max_files != 0 && it > max_files) break;
    //
    //     fprintf(pipe, "set output '%s_%d.png'\n", f_name, it);
    //     fprintf(pipe, "splot ");
    //     
    //     snprintf(buff, 1024, "w polygons fs transparent solid 0.2 fc rgb '%s' notitle", colors[it%8]);
    //     plot_add_ncell(pipe, setup, current, buff);
    //
    //     for (int jj=0; jj<setup->N_points; jj++) {
    //         fprintf(pipe, "\"<echo \'");
    //         fprintf(pipe, "%f %f %f\\n", setup->points[jj][0], setup->points[jj][1], setup->points[jj][2]);
    //         fprintf(pipe, "'\" pt 7 lc rgb 'black' notitle, ");
    //     }
    //     fprintf(pipe, "\n");
    //
    //     it++;
    //     current = current->next;
    // }
    // pclose(pipe);
}

