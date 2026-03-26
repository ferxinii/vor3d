
#include "scplx.h"
#include "points.h"
#include "gtests.h"
#include "gnuplotc.h"
#include "dynarray.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <float.h>

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


void free_complex(s_scplx *scplx)
{
    free_points(&scplx->points);
    if (scplx->weights) free(scplx->weights);

    s_ncell *current = scplx->head;
    while (current) {
        s_ncell *next = current->next;
        free_ncell(current);
        current = next;
    }

    memset(scplx, 0, sizeof(s_scplx));
}


void print_ncell(const s_ncell *ncell)
{
    printf("%p, ( ", (void*)ncell);
    for (int ii=0; ii<4; ii++) {
        printf("%d ", ncell->vertex_id[ii]);
    }
    printf(") \n");
}


void print_scomplex(const s_scplx *scplx)
{
    puts("NCELLS");
    s_ncell *current = scplx->head;
    int ii = 0;
    while (current) {
        printf("%d  |  p : %p  |  prev : %p  |  marked : %d  |  vertex_ids :", ii, (void*)current, (void*)current->prev, current->mark_token);
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


// void write_ncell3d_file(const s_scplx *scplx, const s_ncell *ncell, FILE *file)
// {
//     fprintf(file, "%f %f %f\n", scplx->points.p[ncell->vertex_id[0]].x,
//                                 scplx->points.p[ncell->vertex_id[0]].y,
//                                 scplx->points.p[ncell->vertex_id[0]].z);
//     fprintf(file, "%f %f %f\n", scplx->points.p[ncell->vertex_id[1]].x,
//                                 scplx->points.p[ncell->vertex_id[1]].y,
//                                 scplx->points.p[ncell->vertex_id[1]].z);
//     fprintf(file, "%f %f %f\n\n", scplx->points.p[ncell->vertex_id[2]].x,
//                                   scplx->points.p[ncell->vertex_id[2]].y,
//                                   scplx->points.p[ncell->vertex_id[2]].z);
//
//     fprintf(file, "%f %f %f\n", scplx->points.p[ncell->vertex_id[0]].x,
//                                 scplx->points.p[ncell->vertex_id[0]].y,
//                                 scplx->points.p[ncell->vertex_id[0]].z);
//     fprintf(file, "%f %f %f\n", scplx->points.p[ncell->vertex_id[1]].x,
//                                 scplx->points.p[ncell->vertex_id[1]].y,
//                                 scplx->points.p[ncell->vertex_id[1]].z);
//     fprintf(file, "%f %f %f\n\n", scplx->points.p[ncell->vertex_id[3]].x,
//                                   scplx->points.p[ncell->vertex_id[3]].y,
//                                   scplx->points.p[ncell->vertex_id[3]].z);
//
//     fprintf(file, "%f %f %f\n", scplx->points.p[ncell->vertex_id[3]].x,
//                                 scplx->points.p[ncell->vertex_id[3]].y,
//                                 scplx->points.p[ncell->vertex_id[3]].z);
//     fprintf(file, "%f %f %f\n", scplx->points.p[ncell->vertex_id[1]].x,
//                                 scplx->points.p[ncell->vertex_id[1]].y,
//                                 scplx->points.p[ncell->vertex_id[1]].z);
//     fprintf(file, "%f %f %f\n\n", scplx->points.p[ncell->vertex_id[2]].x,
//                                   scplx->points.p[ncell->vertex_id[2]].y,
//                                   scplx->points.p[ncell->vertex_id[2]].z);
//
//     fprintf(file, "%f %f %f\n", scplx->points.p[ncell->vertex_id[0]].x,
//                                 scplx->points.p[ncell->vertex_id[0]].y,
//                                 scplx->points.p[ncell->vertex_id[0]].z);
//     fprintf(file, "%f %f %f\n", scplx->points.p[ncell->vertex_id[3]].x,
//                                 scplx->points.p[ncell->vertex_id[3]].y,
//                                 scplx->points.p[ncell->vertex_id[3]].z);
//     fprintf(file, "%f %f %f\n\n\n", scplx->points.p[ncell->vertex_id[2]].x,
//                                     scplx->points.p[ncell->vertex_id[2]].y,
//                                     scplx->points.p[ncell->vertex_id[2]].z);
// }
// void write_scomplex_file(const s_scplx *scplx, FILE *file)
// {
//     s_ncell *current = scplx->head;
//     while (current) {
//         write_ncell3d_file(scplx, current, file);
//         fprintf(file, "\n\n");
//         current = current->next;
//     }
// }



void extract_vertices_ncell(const s_scplx *scplx, const s_ncell *ncell, s_point out[4])
{
    for (int ii=0; ii<4; ii++)
        out[ii] = scplx->points.p[ncell->vertex_id[ii]];
}


void extract_weights_ncell(const s_scplx *scplx, const s_ncell *ncell, double out[4])
{
    if (!scplx->weights) { out[0]=0; out[1]=0; out[2]=0; out[3]=0; return; }
    for (int ii=0; ii<4; ii++)
        out[ii] = scplx->weights[ncell->vertex_id[ii]];
}



void extract_vertices_and_weights_ncell(const s_scplx *scplx, const s_ncell *ncell, s_point p_out[4], double w_out[4])
{
    for (int ii=0; ii<4; ii++) {
        int vid = ncell->vertex_id[ii];
        p_out[ii] = scplx->points.p[vid];
        w_out[ii] = scplx->weights[vid];
    }
}


void extract_ids_face(const s_ncell *ncell, int dim_face, const int *v_localid, int *out)
{   // v_localid[dim-dim_face], out[dim_face+1]
    if (dim_face < 0 || dim_face >= 3) {
        fprintf(stderr, "scplx.c: dim_face must be >= 0 && < 3\n");
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


void extract_vertices_face(const s_scplx *scplx, const s_ncell *ncell, int dim_face, 
                           const int v_localid[3-dim_face], s_point out[dim_face+1])
{   // v_localid[dim-dim_face], out[dim_face+1]
    if (dim_face < 0 || dim_face >= 3) {
        fprintf(stderr, "scplx.c: dim_face must be >= 0 && < 3\n");
        exit(1);
    }

    int face_vertex_id[dim_face+1];
    extract_ids_face(ncell, dim_face, v_localid, face_vertex_id);
    for (int ii=0; ii<dim_face+1; ii++) {
        out[ii] = scplx->points.p[face_vertex_id[ii]];
    }
}


void extract_weights_face(const s_scplx *scplx, const s_ncell *ncell, int dim_face,
                          const int v_localid[3-dim_face], double out[dim_face+1])
{
    if (dim_face < 0 || dim_face >= 3) {
        fprintf(stderr, "scplx.c: dim_face must be >= 0 && < 3\n");
        exit(1);
    }

    if (!scplx->weights) {
        for (int i=0; i<dim_face+1; i++) out[i] = 0;
        return;
    }

    int ids[dim_face+1];
    extract_ids_face(ncell, dim_face, v_localid, ids);
    for (int i=0; i<dim_face+1; i++)
        out[i] = scplx->weights ? scplx->weights[ids[i]] : 0.0;
}


void extract_vertices_and_weights_face(const s_scplx *scplx, const s_ncell *ncell, int dim_face,
                                       const int v_localid[3-dim_face],
                                       s_point pts_out[dim_face+1], double wts_out[dim_face+1])
{
    if (dim_face < 0 || dim_face >= 3) {
        fprintf(stderr, "scplx.c: dim_face must be >= 0 && < 3\n");
        exit(1);
    }

    int ids[dim_face+1];
    extract_ids_face(ncell, dim_face, v_localid, ids);
    for (int i = 0; i < dim_face+1; i++) {
        pts_out[i] = scplx->points.p[ids[i]];
        wts_out[i] = scplx->weights ? scplx->weights[ids[i]] : 0.0;
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


s_ncell *next_ncell_ridge_cycle(const s_ncell *ncell, int v_localid_main, int v_localid_2,
                                int *new_v_localid_main, int *new_v_localid_2)
{
    s_ncell *next = ncell->opposite[v_localid_main];
    *new_v_localid_main = id_where_equal_int(next->vertex_id, 4, ncell->vertex_id[v_localid_2]);
    for (int ii=0; ii<4; ii++) {
        if (next->opposite[ii] == ncell) {
            *new_v_localid_2 = ii;
            return next;
        }
    }
    return NULL;
}


void walk_ridge_cycle_and_check_ncells(const s_ncell *start, const int omitted[2],
                                       bool (*check)(void *ctx, const s_ncell *), void *ctx)
{
    const s_ncell *cur = start;
    int cur_main = omitted[0], cur_2 = omitted[1];
    bool hit_boundary = false;
    do {
        if (!check(ctx, cur)) return;
        if (!cur->opposite[cur_main]) { hit_boundary = true; break; }
        int next_main, next_2;
        const s_ncell *next = next_ncell_ridge_cycle(cur, cur_main, cur_2,
                                                     &next_main, &next_2);
        if (next == start || !next) break;
        cur = next; cur_main = next_main; cur_2 = next_2;
    } while (cur != start);

    if (hit_boundary) {
        cur = start; cur_main = omitted[1]; cur_2 = omitted[0];
        if (!cur->opposite[cur_main]) return;  // boundary on both sides
        int next_main, next_2;
        cur = next_ncell_ridge_cycle(cur, cur_main, cur_2, &next_main, &next_2);
        cur_main = next_main; cur_2 = next_2;
        while (cur && cur != start) {
            if (!check(ctx, cur)) return;
            if (!cur->opposite[cur_main]) break;
            const s_ncell *next = next_ncell_ridge_cycle(cur, cur_main, cur_2,
                                                          &next_main, &next_2);
            if (next == start || !next) break;
            cur = next; cur_main = next_main; cur_2 = next_2;
        }
    }
}



static int new_mark_stamp(s_scplx *scplx)
{
    scplx->mark_stamp++;
    return scplx->mark_stamp;
}

static int mark_ncells_incident_face_STEP(int mark_stamp, s_ncell *ncell, int dim_face, const int v_localid[3-dim_face], s_dynarray *out_ncell_ptrs)
{
    for (int ii=0; ii<4; ii++) if (inarray(v_localid, 3-dim_face, ii)) {
        s_ncell *adjacent_ncell = ncell->opposite[ii];  
        if (adjacent_ncell && adjacent_ncell->mark_token != mark_stamp) {
            adjacent_ncell->mark_token = mark_stamp;
            if (!dynarray_push(out_ncell_ptrs, &adjacent_ncell)) return 0;
            
            int new_v_localid[3-dim_face];
            face_localid_of_adjacent_ncell(ncell, dim_face, v_localid, ii, new_v_localid);

            /* Recursion: */
            if (!mark_ncells_incident_face_STEP(mark_stamp, adjacent_ncell, dim_face, new_v_localid, out_ncell_ptrs)) return 0;
        }
    }
    return 1;
}

int ncells_incident_face(s_scplx *scplx, s_ncell *ncell, int dim_face, const int v_localid[3-dim_face], s_dynarray *out_ncell_ptrs)
{   /* 0 ERROR, 1 OK */
    if (dim_face < 0 || dim_face >= 3) {
        fprintf(stderr, "dim_face must be >= 0 && < 3\n");
        exit(1);
    }

    int mark_stamp = new_mark_stamp(scplx);
    ncell->mark_token = mark_stamp;

    out_ncell_ptrs->N = 0;
    if (!dynarray_push(out_ncell_ptrs, &ncell)) return 0;
    
    /* Recursion, Flood-Fill */
    if (!mark_ncells_incident_face_STEP(mark_stamp, ncell, dim_face, v_localid, out_ncell_ptrs)) return 0;
    return 1;
}


int are_locally_delaunay(const s_scplx *scplx, const s_ncell *ncell, int id_opposite, e_delaunay_test_type type)
{   // Points on circumscribed sphere are VALID
    s_point coords1[4];
    s_point coords2[4];
    extract_vertices_ncell(scplx, ncell, coords1);
    extract_vertices_ncell(scplx, ncell->opposite[id_opposite], coords2);

    // assert(!(orientation_robust(coords1, coords1[3]) == 0 && 
    //          orientation_robust(coords2, coords2[3]) == 0) &&
    //        "Cannot check delaunayness, both ncells are degenerate.");

    if (orientation_robust(coords1, coords1[3]) == 0 || 
        orientation_robust(coords2, coords2[3]) == 0) return 1;  /* TODO unsure of what to return! */

    // Extract vertex_id of opposite's cell face
    int opp_face_localid;
    face_localid_of_adjacent_ncell(ncell, 2, &id_opposite, id_opposite, &opp_face_localid);
    int opp_face_vertex_id = (ncell->opposite[id_opposite])->vertex_id[opp_face_localid];
    
    // Perform insphere test: weighted or unweighted.
    int in1;
    if (scplx->weights) {
        double weights1[4]; extract_weights_ncell(scplx, ncell, weights1);
        in1 = insphere_weighted_robust(coords1, weights1, 
                scplx->points.p[opp_face_vertex_id], scplx->weights[opp_face_vertex_id]);
    } else {
        in1 = insphere_robust(coords1, scplx->points.p[opp_face_vertex_id]);
    }

    switch (type) {
        case DELAUNAY_TEST_STRICT:
            if (in1 == -1) return 1;
            else return 0;     
        case DELAUNAY_TEST_NONSTRICT:
            if (in1 != 1) return 1;
            else return 0;
    }
}


int is_delaunay_3d(const s_scplx *scplx, e_delaunay_test_type type)
{
    s_point vertices_ncell[4];
    s_ncell *current = scplx->head;
    while (current) {
        extract_vertices_ncell(scplx, current, vertices_ncell);
        if (orientation_robust(vertices_ncell, vertices_ncell[3]) == 0) {
            puts("Flat tetra!!");
            return 0;
        }
        for (int ii=0; ii<4; ii++) {
            if (current->opposite[ii] &&
                !are_locally_delaunay(scplx, current, ii, type)) return 0;
        }
        current = current->next;
    }
    return 1;
}

int test_point_in_ncell(const s_scplx *scplx, const s_ncell *ncell, s_point query)
{   // TODO Assumes consistent ordering of vertices?
    s_point v0 = scplx->points.p[ncell->vertex_id[0]];
    s_point v1 = scplx->points.p[ncell->vertex_id[1]];
    s_point v2 = scplx->points.p[ncell->vertex_id[2]];
    s_point v3 = scplx->points.p[ncell->vertex_id[3]];
    s_point tetra[4] = {v0, v1, v2, v3};
    return test_point_in_tetrahedron(tetra, query, 0, 0);
 }


s_ncell *bruteforce_find_ncell_containing(const s_scplx *scplx, s_point p)
{
    s_ncell *current = scplx->head;
    while (current) {
        e_geom_test test = test_point_in_ncell(scplx, current, p);
        if (test == TEST_IN || test == TEST_BOUNDARY) return current;
        if (test == TEST_ERROR) fprintf(stderr, "Warning: test error\n");
        // if (test == TEST_DEGENERATE);  /* These are expected! Just skip them */
        current = current->next;
    }
    fprintf(stderr, "Warning: Did not find container ncell.\n");
    fprintf(stderr, "p: %g, %g, %g\n", p.x, p.y, p.z);
    exit(1);
}


static s_point ncell_centroid(const s_scplx *scplx, const s_ncell *ncell) {
    s_point vertices[4];
    extract_vertices_ncell(scplx, ncell, vertices);
    s_point out = { .x = (vertices[0].x + vertices[1].x + vertices[2].x + vertices[3].x) / 4.0,
                    .y = (vertices[0].y + vertices[1].y + vertices[2].y + vertices[3].y) / 4.0, 
                    .z = (vertices[0].z + vertices[1].z + vertices[2].z + vertices[3].z) / 4.0 };
    return out;
}


static void random_order_04(int out[4])
{
    for (int ii=0; ii<4; ii++) out[ii] = ii;

    for (int ii=3; ii>0; ii--) {
        int jj = rand() % (ii + 1);
        int tmp = out[ii];
        out[ii] = out[jj];
        out[jj] = tmp;
    }
}

s_ncell *in_ncell_walk(const s_scplx *scplx, s_point p)
{
    // return bruteforce_find_ncell_containing(scplx, p);

    s_ncell *current = scplx->head;
    assert(scplx->N_ncells >= 1 && "N_ncells < 1");
    int randi = (rand() % scplx->N_ncells);  /* Select random ncell to start */
    for (int ii=0; ii<randi; ii++) 
        current = current->next;   

    s_point facet_vertices[3];
    s_ncell *prev = current;

    int steps = 0;
    STEP:
    steps++;
    if (steps > 3 * scplx->N_ncells) {
        fprintf(stderr, "Warning: in_ncell_walk seems to be in a loop. Now bruteforcing to find inside which ncell.\n");
        return bruteforce_find_ncell_containing(scplx, p);
    }
    int order[4]; random_order_04(order);
    for (int kk=0; kk<4; kk++) {  /* Visit faces in random order to prevent loops ? */
        int ii = order[kk];
        s_point opposite_vertex = scplx->points.p[current->vertex_id[ii]];

        s_ncell *next = current->opposite[ii];
        if (!next) continue;

        extract_vertices_face(scplx, current, 2, &ii, facet_vertices);

        int o1 = orientation_robust(facet_vertices, opposite_vertex);
        int o2 = orientation_robust(facet_vertices, p);
        
        if (o1 == 0) {  /* Tetrahedron is degenerate */
            if (next != prev) {  /* We come from a different adjacent one, so walk towards */
                prev = current;
                current = next;
                goto STEP;
            } 
            else continue; 
        }

        if (o2 == 0) {  /* Query is coplanar with face */
            e_geom_test test = test_point_in_triangle_3D(facet_vertices, p, 0, 0);
            if (test == TEST_IN || test == TEST_BOUNDARY) { return current; }
            if (next != prev) {
                // TIE-BREAKING: move to the neighbor whose centroid is closer to p.
                double d_curr = distance_squared(ncell_centroid(scplx, current), p);
                double d_next = distance_squared(ncell_centroid(scplx, next), p);
                if (d_next < d_curr) {
                    prev = current;
                    current = next;
                    goto STEP;
                }
            } 
            continue;
        } else if (o1 != o2) {  // Regular step
            prev = current;
            current = next;
            goto STEP;
        }
    }
    
    return current; 
}


double power_distance_point_vertex(const s_scplx *scplx, int vid, s_point p)
{
    if (!scplx->weights) return distance_squared(scplx->points.p[vid], p);
    else return distance_squared(scplx->points.p[vid], p) - scplx->weights[vid];
}




/* PLOTS */
void plot_add_ncell(s_gnuplot *interface, const s_scplx *scplx, const s_ncell *ncell, char *config)
{
    s_point face_vertices[3];
    for (int ii=0; ii<4; ii++) {
        extract_vertices_face(scplx, ncell, 2, &ii, face_vertices);
        draw_solid_triangle_3d(interface, face_vertices[0].coords, face_vertices[1].coords,
                               face_vertices[2].coords, config);
    }
}


void plot_ncell_3d(const s_scplx *scplx, const s_ncell *ncell, char *f_name, s_point ranges[2])
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

    plot_add_ncell(interface, scplx, ncell, "fs transparent solid 0.2 fc rgb '#000090'");
    gnuplot_end(interface);
}


void plot_all_ncells_3d(const s_scplx *scplx, char *f_name, s_point ranges[2], char *view_command)
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
    s_ncell *current = scplx->head;
    while (current) {
        snprintf(buff, 1024, "fs transparent solid 0.2 fc rgb '%s'", colors[it%8]);
        plot_add_ncell(interface, scplx, current, buff);
        current = current->next;
        it++;
    }
    gnuplot_end(interface);
}


void plot_dt_3d_differentviews(const s_scplx *scplx, char *f_name, s_point ranges[2])
{
    char final_name[512];
    snprintf(final_name, 512, "%s_v1.png", f_name);
    plot_all_ncells_3d(scplx, final_name, ranges, "set view 100, 60, 1.5");

    snprintf(final_name, 512, "%s_v2.png", f_name);
    plot_all_ncells_3d(scplx, final_name, ranges, "set view 100, 90, 1.5");

    snprintf(final_name, 512, "%s_v3.png", f_name);
    plot_all_ncells_3d(scplx, final_name, ranges, "set view 100, 180, 1.5");

    snprintf(final_name, 512, "%s_v4.png", f_name);
    plot_all_ncells_3d(scplx, final_name, ranges, "set view 100, 270, 1.5");

}

