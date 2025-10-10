// [ ] TODO Improve malloc of vertices, or check if i have reached VCELL_MAX_VERTICES to increase size as needed
#include "vdiagram.h"
#include "simplical_complex.h"
#include "bpoly.h"
#include "algebra.h"
#include "geometry.h"
#include "bpoly.h"
#include "external/gnuplotC/gnuplotc.h"
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>


// DEBUG:
void dump_setup_points(const s_scplx *setup, const int *ids, int n)
{
    for (int i = 0; i < n; ++i) {
        int id = ids[i];
        if (id < 0 || id >= setup->N_points) {
            fprintf(stderr, "dump: bad id %d\n", id);
            continue;
        }
        s_point p = setup->points[id];
        fprintf(stderr, "setup->points[%d] = (%g, %g, %g)\n", id, p.x, p.y, p.z);
    }
}


void free_vcell(s_vcell *vcell)
{
    free(vcell->vertices);
    free_matrix_int(vcell->origin_vertices, vcell->Nv_capacity);
    if (vcell->faces) {
        free(vcell->faces);
        free(vcell->fnormals);
    }
    free(vcell);
}


void free_vdiagram(s_vdiagram *vdiagram)
{
    for (int ii=0; ii<vdiagram->N; ii++) {
        // printf("Freing vcell %d\n", ii);
        if (vdiagram->vcells[ii]) free_vcell(vdiagram->vcells[ii]);
    }
    // puts("Freeing vcells");
    free(vdiagram->vcells);
    // puts("Freeing seeds");
    free(vdiagram->seeds);
    // puts("Freeing bpoly");
    free_bpoly((s_bpoly *)vdiagram->bpoly);
    // puts("Freeing vdiagram");
    free(vdiagram);
}


void write_vcell_file(const s_vcell *vcell, FILE *file)
{
    for (int ii=0; ii<vcell->Nf; ii++) {
        fprintf(file, "%f %f %f\n", vcell->vertices[vcell->faces[ii*3 + 0]].x,
                                    vcell->vertices[vcell->faces[ii*3 + 0]].y,
                                    vcell->vertices[vcell->faces[ii*3 + 0]].z);
        fprintf(file, "%f %f %f\n", vcell->vertices[vcell->faces[ii*3 + 1]].x,
                                    vcell->vertices[vcell->faces[ii*3 + 1]].y,
                                    vcell->vertices[vcell->faces[ii*3 + 1]].z);
        fprintf(file, "%f %f %f\n\n", vcell->vertices[vcell->faces[ii*3 + 2]].x,
                                      vcell->vertices[vcell->faces[ii*3 + 2]].y,
                                      vcell->vertices[vcell->faces[ii*3 + 2]].z);
    }
}


void write_vd_file(const s_vdiagram *vd, FILE *file)
{
    for (int ii=0; ii<vd->N; ii++) {
        write_vcell_file(vd->vcells[ii], file);
        fprintf(file, "\n\n");
    }
}


s_vdiagram *malloc_vdiagram(const s_scplx *setup, int Nreal)
{
    s_vdiagram *out = malloc(sizeof(s_vdiagram));
    out->seeds = malloc(sizeof(s_point) * setup->N_points);
    memcpy(out->seeds, setup->points, setup->N_points * sizeof(s_point));

    out->N = Nreal;
    out->vcells = malloc(sizeof(s_vcell*) * Nreal);
    for (int ii=0; ii<Nreal; ii++) out->vcells[ii] = NULL;
    out->bpoly = NULL;
    return out;
}


void print_vcell(const s_vcell *vcell)
{
    puts("VCELL");
    printf("Seed: %d, Nv = %d, vol = %f\n", vcell->seed_id, vcell->Nv, vcell->volume);
    printf("%f, %f, %f; %d, %d, %d, %d\n", vcell->vertices[0].x, vcell->vertices[0].y,
                            vcell->vertices[0].z, vcell->origin_vertices[0][3],
                            vcell->origin_vertices[0][0], vcell->origin_vertices[0][1], 
                            vcell->origin_vertices[0][2]);
    for (int ii=1; ii<vcell->Nv; ii++) {
        printf("%f, %f, %f; %d, %d, %d, %d\n", vcell->vertices[ii].x, vcell->vertices[ii].y, 
                            vcell->vertices[ii].z, vcell->origin_vertices[ii][3],
                            vcell->origin_vertices[ii][0], vcell->origin_vertices[ii][1],
                            vcell->origin_vertices[ii][2]);
    }
    printf("\n");
}


void print_vdiagram(const s_vdiagram *vdiagram)
{
    puts("------- VORONOI DIAGRAM ------");
    for (int ii=0; ii<vdiagram->N; ii++) {
        if (vdiagram->vcells[ii]) print_vcell(vdiagram->vcells[ii]);
    }
    puts("------------------------------");
}


s_vcell *malloc_vcell(int seed_id)
{
    s_vcell *out = malloc(sizeof(s_vcell));
    out->seed_id = seed_id;
    out->Nv = 0;
    out->Nv_capacity = VCELL_BLOCK_VERTICES;
    out->vertices = malloc(sizeof(s_point) * VCELL_BLOCK_VERTICES);
    out->origin_vertices = malloc_matrix_int(VCELL_BLOCK_VERTICES, 4);
    return out;
}


void increase_num_vertices_if_needed(s_vcell *vcell)
{
   // Check if we need to increase the capacity.
    if (vcell->Nv + 1 >= vcell->Nv_capacity) {
        int new_capacity = vcell->Nv_capacity + VCELL_BLOCK_VERTICES;

        vcell->vertices = realloc(vcell->vertices, sizeof(s_point) * new_capacity);
        vcell->origin_vertices = realloc_matrix_int(vcell->origin_vertices, vcell->Nv, new_capacity, 4);
        // printf("DEBUG: Increased ncell capacity: old=%d, new=%d\n", vcell->Nv_capacity, new_capacity);
        vcell->Nv_capacity = new_capacity;
    }
}


int add_vvertex_from_ncell(const s_scplx *setup, const s_ncell *ncell, s_vcell *vcell)
{   // Returns the index of the vertex
    increase_num_vertices_if_needed(vcell);

    s_point v[4];
    extract_vertices_ncell(setup, ncell, v);

    // int ids[] = { ncell->vertex_id[0], ncell->vertex_id[1],
    //           ncell->vertex_id[2], ncell->vertex_id[3] };
    // dump_setup_points(setup, ids, 4);
    // // Debug: print the original four points from the ncell
    // fprintf(stderr, "DEBUG: add_vvertex_from_ncell: ncell->count=%d, vertex_ids = [%d,%d,%d,%d]\n",
    //         ncell->count,
    //         ncell->vertex_id[0], ncell->vertex_id[1],
    //         ncell->vertex_id[2], ncell->vertex_id[3]);
    // for (int i=0;i<4;i++) {
    //     fprintf(stderr, "  v[%d] = (%g, %g, %g)\n", i, v[i].x, v[i].y, v[i].z);
    // }


    s_point a = subtract_points(v[1], v[0]);
    s_point b = subtract_points(v[2], v[0]);
    s_point c = subtract_points(v[3], v[0]);

    double a2 = norm_squared(a);
    double b2 = norm_squared(b);
    double c2 = norm_squared(c);
    
    s_point bc = cross_prod(b, c);
    s_point ca = cross_prod(c, a);
    s_point ab = cross_prod(a, b);

    double f = 2 * dot_prod(ab, c);
    if (fabs(f) < 1e-9) {
        // printf("circumcenter: NEARLY SINGULAR!\n");
        return -1;
    }
    f = 1.0 / f;

    s_point circumcenter = {{{
        f * ( a2 * bc.x + b2 * ca.x + c2 * ab.x) + v[0].x,
        f * ( a2 * bc.y + b2 * ca.y + c2 * ab.y) + v[0].y,
        f * ( a2 * bc.z + b2 * ca.z + c2 * ab.z) + v[0].z,
    }}};

    for (int ii=0; ii<vcell->Nv; ii++) {
        if (distance_squared(circumcenter, vcell->vertices[ii]) < 1e-5) {
            return -1;
        }
    }

    vcell->vertices[vcell->Nv] = circumcenter;

    vcell->origin_vertices[vcell->Nv][3] = ncell->count;
    vcell->Nv++;

    return vcell->Nv - 1;
}


void compute_vcell_volume(s_vcell *vcell)
{
    double vol = 0;
    for (int ii=0; ii<vcell->Nf; ii++) {
        int i0 = 0;
        int i1 = 1;
        int i2 = 2;

        double Nx = (vcell->vertices[vcell->faces[ii*3 + 1]].coords[i1] -
                     vcell->vertices[vcell->faces[ii*3 + 0]].coords[i1]) *
                    (vcell->vertices[vcell->faces[ii*3 + 2]].coords[i2] -
                     vcell->vertices[vcell->faces[ii*3 + 0]].coords[i2]) 
                    -
                    (vcell->vertices[vcell->faces[ii*3 + 1]].coords[i2] -
                     vcell->vertices[vcell->faces[ii*3 + 0]].coords[i2]) *
                    (vcell->vertices[vcell->faces[ii*3 + 2]].coords[i1] -
                     vcell->vertices[vcell->faces[ii*3 + 0]].coords[i1]);

        vol += Nx * (vcell->vertices[vcell->faces[ii*3 + 0]].coords[i0] +
                     vcell->vertices[vcell->faces[ii*3 + 1]].coords[i0] +
                     vcell->vertices[vcell->faces[ii*3 + 2]].coords[i0]);
    }
    vcell->volume = vol / 6;
}


int add_convex_hull_vcell(s_vcell *vcell)
{
    convhull_from_points(vcell->vertices, vcell->Nv, &vcell->faces, &vcell->fnormals, &vcell->Nf);
    if(!vcell->faces) return 0;
    return 1;
}


int bounded_extraction(const s_scplx *setup, s_vcell *vcell)
{
    s_ncell *current = setup->head;
    while (current) {
        if (current->mark == 1) add_vvertex_from_ncell(setup, current, vcell);
        current = current->next;
    }
    // printf("DEBUG: NV=%d\n", vcell->Nv);
    int out = add_convex_hull_vcell(vcell);
    return out;
}


s_vcell *extract_voronoi_cell(const s_scplx *setup, int vertex_id, const s_bpoly *bp)
{
    // Find an ncell with this vertex
    s_ncell *ncell = setup->head;
    while (!inarray(ncell->vertex_id, 4, vertex_id)) {
        ncell = ncell->next;
    }
    // Find "complementary" indices to localid
    int v_localid = id_where_equal_int(ncell->vertex_id, 4, vertex_id);
    int v_localid_COMP[3];
    int kk=0;
    for (int ii=0; ii<4; ii++) {
        if (ii != v_localid) {
            v_localid_COMP[kk] = ii;
            kk++;
        }
    }
    // Mark incident ncells to this point
    initialize_ncells_mark(setup);
    mark_ncells_incident_face(setup, ncell, v_localid_COMP, 0);
    // printf("DEBUG: N=%d\n", count_marked(setup));


    s_vcell *vcell = malloc_vcell(vertex_id);
    int out = bounded_extraction(setup, vcell);
    if (out == 0) {
        puts("extract_voronoi_cell: ERROR extracting convhull of cell!");
        // print_vcell(vcell);
        free_vcell(vcell);
        return NULL;
    }
    
    // Snap into bp if cell spikes outward, KNOWN PROBLEM FIXME TODO
    for (int ii=0; ii<vcell->Nv; ii++) {
        if (!is_inside_convhull(vcell->vertices[ii], bp->points, bp->faces, bp->Nf)) {  
            // (If -1 (on face) dont do anything)
            vcell->vertices[ii] = find_closest_point_on_bp(bp, vcell->vertices[ii]);
            // if (!is_inside_convhull(new, bp->points, bp->faces, bp->fnormals, bp->Nf)) {
            //     puts("WARNING! SNAP IS STILL OUTSIDE?");
            // }
        }
    }

    // Compute volume
    compute_vcell_volume(vcell);
    // printf("VOL: %f\n", vcell->volume);
    return vcell;
}



s_vdiagram *voronoi_from_delaunay_3d(const s_scplx *setup, const s_bpoly *bpoly, int Nreal)
{
    initialize_ncells_counter(setup);
    
    s_vdiagram *vdiagram = malloc_vdiagram(setup, Nreal);
    vdiagram->bpoly = bpoly;
    
    for (int ii=0; ii<Nreal; ii++) {
        // print_vdiagram(vdiagram);
        for (int tries = 0; tries < 2; tries ++) { 
            vdiagram->vcells[ii] = extract_voronoi_cell(setup, ii, bpoly);
            if (vdiagram->vcells[ii] == NULL) continue;
            else break;
        }
        if (vdiagram->vcells[ii] == NULL || vdiagram->vcells[ii]->volume <= 0) {
            // printf("%p, %f\n", (void*)vdiagram->vcells[ii], vdiagram->vcells[ii]->volume);
            // print_vcell(vdiagram->vcells[ii]);
            puts("ERROR: Could not construct vdiagram. Exiting voronoi_from_delaunay_3d");
            free_vdiagram(vdiagram);
            return NULL;
        }
    }

    return vdiagram;
}


int find_inside_which_vcell(const s_vdiagram *vd, s_point x)
{
    for (int ii=0; ii<vd->N; ii++) {
        if (is_inside_convhull(x, vd->vcells[ii]->vertices, 
                               vd->vcells[ii]->faces, vd->vcells[ii]->Nf)) {
            return ii;
        }
    }
    return -1;
}


// ----------------------------------------------------------------------------------------------
// --------------------------------------- PLOTS ------------------------------------------------
// ----------------------------------------------------------------------------------------------

void randomize_colors(int N, char **colors) 
{
    for (int i = N - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        char *tmp = colors[i];
        colors[i] = colors[j];
        colors[j] = tmp;
    }
}


void plot_add_vcell(s_gnuplot *interface, const s_vcell *vcell, char *config)
{
    for (int ii=0; ii<vcell->Nf; ii++) {
        draw_solid_triangle_3d(interface, vcell->vertices[vcell->faces[ii*3]].coords, 
                vcell->vertices[vcell->faces[ii*3+1]].coords, 
                vcell->vertices[vcell->faces[ii*3+2]].coords, config);
        }
}


void plot_add_bpoly(s_gnuplot *interface, const s_bpoly *bpoly, char *config)
{
    for (int ii=0; ii<bpoly->Nf; ii++) {
        draw_solid_triangle_3d(interface, bpoly->points[bpoly->faces[ii*3]].coords, 
                bpoly->points[bpoly->faces[ii*3+1]].coords, 
                bpoly->points[bpoly->faces[ii*3+2]].coords, config);
    }
}


void plot_vcell(const s_vdiagram *vdiag, const s_vcell *vcell, char *f_name, const s_point ranges[2])
{
    int size[2] = {1080, 1080};
    s_gnuplot *interface = gnuplot_start(PNG_3D, f_name, size, 18);
    gnuplot_config(interface,  "set pm3d depthorder",
                               "set pm3d border lc 'black' lw 0.1",
                               "set view 100, 10", 
                               "set xyplane at 0");
    if (ranges) {
        char buff[1024];
        snprintf(buff, 1024, "set xrange [%f:%f]\n set yrange [%f:%f]\n set zrange [%f:%f]", 
                 ranges[0].x, ranges[1].x, ranges[0].y, ranges[1].y, ranges[0].z, ranges[1].z); 
        gnuplot_config(interface, buff);
    }

    plot_add_vcell(interface, vcell, "fs transparent solid 0.6");
    plot_add_bpoly(interface, vdiag->bpoly, "fs transparent solid 0.1");
    gnuplot_end(interface);
}


void plot_all_vcells(const s_vdiagram *vdiagram, char *f_name, const s_point ranges[2], char *view_command)
{
    char buff[1024];
    char *colors[] = { "#000090", "#ee0000", "#7f0000", "#0090ff", "#0fffee", 
        "#90ff70", "#ffee00", "#000fff",  "#ff7000" };
    randomize_colors(9, colors);
    
    int size[2] = {2160, 2160};
    s_gnuplot *interface = gnuplot_start(PNG_3D, f_name, size, 18);
    gnuplot_config(interface,  "set pm3d depthorder",
                               "set pm3d border lc 'black' lw 0.01",
                               "unset border",
                               "unset xtics",
                               "unset ytics",
                               "unset ztics",
                               view_command);
    if (ranges) {
        snprintf(buff, 1024, "set xrange [%f:%f]\n set yrange [%f:%f]\n set zrange [%f:%f]", 
                 ranges[0].x, ranges[1].x, ranges[0].y, ranges[1].y, ranges[0].z, ranges[1].z); 
        gnuplot_config(interface, buff);
    }


    for (int ii=0; ii<vdiagram->N; ii++) {
        s_vcell *vcell = vdiagram->vcells[ii];
        snprintf(buff, 1024, "fs transparent solid 0.2 fc rgb '%s'", colors[ii%8]);
        plot_add_vcell(interface, vcell, buff);  
    }
    gnuplot_end(interface);
}


void plot_vdiagram_differentviews(const s_vdiagram *vdiagram, char *f_name, const s_point ranges[2])
{
    char final_name[1024];
    snprintf(final_name, 1024, "%s_v1.png", f_name);
    plot_all_vcells(vdiagram, final_name, ranges, "set view 100, 60, 1.5"); 

    snprintf(final_name, 1024, "%s_v2.png", f_name);
    plot_all_vcells(vdiagram, final_name, ranges, "set view 100, 90, 1.5");

    snprintf(final_name, 1024, "%s_v3.png", f_name);
    plot_all_vcells(vdiagram, final_name, ranges, "set view 100, 180, 1.5");

    snprintf(final_name, 1024, "%s_v4.png", f_name);
    plot_all_vcells(vdiagram, final_name, ranges, "set view 100, 270, 1.5");
}



void clear_volumes_file(char *fname)
{
    FILE *f = fopen(fname, "w");
    fclose(f);
}


void append_volumes_to_file(s_vdiagram *vdiagram, char *fname, int id)
{
    FILE *f = fopen(fname, "a");
    assert(f && "Could not open file to write volumes to");

    for (int ii=0; ii<vdiagram->N; ii++) {
        fprintf(f, "%d, %f, %f, %f, %f\n", id, vdiagram->seeds[vdiagram->vcells[ii]->seed_id].x,
                                               vdiagram->seeds[vdiagram->vcells[ii]->seed_id].y, 
                                               vdiagram->seeds[vdiagram->vcells[ii]->seed_id].z, 
                                               vdiagram->vcells[ii]->volume);   
    }
    fclose(f);
}

