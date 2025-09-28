// [ ] TODO Improve malloc of vertices, or check if i have reached VCELL_MAX_VERTICES to increase size as needed
#include "vdiagram.h"
#include "simplical_complex.h"
#include "bpoly.h"
#include "algebra.h"
#include "geometry.h"
#include "bpoly.h"
#include "gnuplotc.h"
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>


void free_vcell(s_vcell *vcell)
{
    free_matrix(vcell->vertices, vcell->Nv_capacity);
    free_matrix_int(vcell->origin_vertices, vcell->Nv_capacity);
    free(vcell->faces);
    free_matrix(vcell->fnormals, vcell->Nf);
    free(vcell);
}


void free_vdiagram(s_vdiagram *vdiagram)
{
    for (int ii=0; ii<vdiagram->N_vcells; ii++) {
        if (vdiagram->vcells[ii]) free_vcell(vdiagram->vcells[ii]);
    }
    free(vdiagram->vcells);
    free_matrix(vdiagram->seeds, vdiagram->N_vcells);
    free_bpoly((s_bound_poly *)vdiagram->bpoly);
    free(vdiagram);
}


void write_vcell_file(s_vcell *vcell, FILE *file)
{
    for (int ii=0; ii<vcell->Nf; ii++) {
        fprintf(file, "%f %f %f\n", vcell->vertices[vcell->faces[ii*3 + 0]][0],
                                     vcell->vertices[vcell->faces[ii*3 + 0]][1],
                                     vcell->vertices[vcell->faces[ii*3 + 0]][2]);
        fprintf(file, "%f %f %f\n", vcell->vertices[vcell->faces[ii*3 + 1]][0],
                                     vcell->vertices[vcell->faces[ii*3 + 1]][1],
                                     vcell->vertices[vcell->faces[ii*3 + 1]][2]);
        fprintf(file, "%f %f %f\n\n", vcell->vertices[vcell->faces[ii*3 + 2]][0],
                                      vcell->vertices[vcell->faces[ii*3 + 2]][1],
                                      vcell->vertices[vcell->faces[ii*3 + 2]][2]);
    }
}


void write_vd_file(s_vdiagram *vd, FILE *file)
{
    for (int ii=0; ii<vd->N_vcells; ii++) {
        write_vcell_file(vd->vcells[ii], file);
        fprintf(file, "\n\n");
    }
}


s_vdiagram *malloc_vdiagram(const s_setup *setup, int Nreal)
{
    s_vdiagram *out = malloc(sizeof(s_vdiagram));
    out->seeds = malloc_matrix(setup->N_points, 3);
    copy_matrix(setup->points, out->seeds, setup->N_points, 3);

    out->N_vcells = Nreal;
    out->vcells = malloc(sizeof(s_vcell*) * Nreal);
    for (int ii=0; ii<Nreal; ii++) out->vcells[ii] = NULL;
    out->bpoly = NULL;
    return out;
}


void print_vcell(const s_vcell *vcell)
{
    puts("VCELL");
    printf("Seed: %d, Nv = %d, vol = %f\n", vcell->seed_id, vcell->Nv, vcell->volume);
    printf("%f, %f, %f; %d, %d, %d, %d\n", vcell->vertices[0][0], vcell->vertices[0][1],
                            vcell->vertices[0][2], vcell->origin_vertices[0][3],
                            vcell->origin_vertices[0][0], vcell->origin_vertices[0][1], 
                            vcell->origin_vertices[0][2]);
    for (int ii=1; ii<vcell->Nv; ii++) {
        printf("%f, %f, %f; %d, %d, %d, %d\n", vcell->vertices[ii][0], vcell->vertices[ii][1], vcell->vertices[ii][2], 
                               vcell->origin_vertices[ii][3], vcell->origin_vertices[ii][0],
                               vcell->origin_vertices[ii][1], vcell->origin_vertices[ii][2]);
    }

    printf("\n");
}


void print_vdiagram(const s_vdiagram *vdiagram)
{
    puts("------- VORONOI DIAGRAM ------");
    for (int ii=0; ii<vdiagram->N_vcells; ii++) {
        print_vcell(vdiagram->vcells[ii]);
    }
    puts("------------------------------");
}


s_vcell *malloc_vcell(int seed_id)
{
    s_vcell *out = malloc(sizeof(s_vcell));
    out->seed_id = seed_id;
    out->Nv = 0;
    out->Nv_capacity = VCELL_BLOCK_VERTICES;
    out->vertices = malloc_matrix(VCELL_BLOCK_VERTICES, 3);
    out->origin_vertices = malloc_matrix_int(VCELL_BLOCK_VERTICES, 4);
    return out;
}


void increase_num_vertices_if_needed(s_vcell *vcell)
{
   // Check if we need to increase the capacity.
    if (vcell->Nv + 1 >= vcell->Nv_capacity) {
        int new_capacity = vcell->Nv_capacity + VCELL_BLOCK_VERTICES;

        vcell->vertices = realloc_matrix(vcell->vertices, vcell->Nv, new_capacity, 3);
        
        vcell->origin_vertices = realloc_matrix_int(vcell->origin_vertices, vcell->Nv, new_capacity, 4);
        
        // printf("DEBUG: Increased ncell capacity: old=%d, new=%d\n", vcell->Nv_capacity, new_capacity);
        vcell->Nv_capacity = new_capacity;
    }
}


int add_vvertex_from_ncell(const s_setup *setup, const s_ncell *ncell, s_vcell *vcell)
{   // Returns the index of the vertex
    increase_num_vertices_if_needed(vcell);

    static double **v = NULL;
    if (!v) v = malloc_matrix(4, 3);

    extract_vertices_ncell(setup, ncell, v);

    double a[3], b[3], c[3];
    subtract_3d(v[1], v[0], a);
    subtract_3d(v[2], v[0], b);
    subtract_3d(v[3], v[0], c);

    double a2 = norm_squared(a, 3);
    double b2 = norm_squared(b, 3);
    double c2 = norm_squared(c, 3);
    
    double bc[3], ca[3], ab[3];
    cross_3d(b, c, bc);
    cross_3d(c, a, ca);
    cross_3d(a, b, ab);

    double f = 2 * dot_3d(ab, c);
    if (fabs(f) < 1e-9) {
        // printf("circumcenter: NEARLY SINGULAR!\n");
        // printf("ab: (%f, %f, %f), c: (%f, %f, %f), denom: %.16f\n", ab[0], ab[1], ab[2], c[0], c[1], c[2],f);
        // print_matrix(v, 4, 3);
        return -1;
    }
    assert(fabs(f) > 1e-9 && "NEARLY SINGULAR?");
    f = 1.0 / f;

    double circumcenter[3] = {
        f * ( a2 * bc[0] + b2 * ca[0] + c2 * ab[0]) + v[0][0],
        f * ( a2 * bc[1] + b2 * ca[1] + c2 * ab[1]) + v[0][1],
        f * ( a2 * bc[2] + b2 * ca[2] + c2 * ab[2]) + v[0][2],
    };

    for (int ii=0; ii<vcell->Nv; ii++) {
        if (norm_difference(circumcenter, vcell->vertices[ii], 3) < 1e-3) {
            return -1;
        }
    }

    vcell->vertices[vcell->Nv][0] = circumcenter[0];
    vcell->vertices[vcell->Nv][1] = circumcenter[1];
    vcell->vertices[vcell->Nv][2] = circumcenter[2];

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

        double Nx = (vcell->vertices[vcell->faces[ii*3 + 1]][i1] -
                     vcell->vertices[vcell->faces[ii*3 + 0]][i1]) *
                    (vcell->vertices[vcell->faces[ii*3 + 2]][i2] -
                     vcell->vertices[vcell->faces[ii*3 + 0]][i2]) 
                    -
                    (vcell->vertices[vcell->faces[ii*3 + 1]][i2] -
                     vcell->vertices[vcell->faces[ii*3 + 0]][i2]) *
                    (vcell->vertices[vcell->faces[ii*3 + 2]][i1] -
                     vcell->vertices[vcell->faces[ii*3 + 0]][i1]);

        vol += Nx * (vcell->vertices[vcell->faces[ii*3 + 0]][i0] +
                     vcell->vertices[vcell->faces[ii*3 + 1]][i0] +
                     vcell->vertices[vcell->faces[ii*3 + 2]][i0]);
    }
    vcell->volume = vol / 6;
}


int add_convex_hull_vcell(s_vcell *vcell)
{
    double CM[3];
    find_center_mass(vcell->vertices, vcell->Nv, 3, CM);
    int *faces, N_faces;
    ch_vertex *vertices_aux = convert_points_to_chvertex(vcell->vertices, vcell->Nv);
    convhull_3d_build(vertices_aux, vcell->Nv, &faces, &N_faces);
    if(!faces) return 0;
    vcell->faces = faces;
    vcell->Nf = N_faces;
    double **fnormals = extract_normals_from_ch(vertices_aux, faces, N_faces, CM);
    vcell->fnormals = fnormals;
    free(vertices_aux);
    return 1;
}


int bounded_extraction(const s_setup *setup, s_vcell *vcell)
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


s_vcell *extract_voronoi_cell(const s_setup *setup, int vertex_id, s_bound_poly *bp)
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
        puts("ERROR EXTRACTING VCELL!");
        // print_vcell(vcell);
        free_vcell(vcell);
        return NULL;
    }
    
    // Snap into bp if cell spikes outward, KNOWN PROBLEM FIXME TODO
    for (int ii=0; ii<vcell->Nv; ii++) {
        if (!is_inside_convhull(vcell->vertices[ii], bp->points, bp->faces, bp->fnormals, bp->Nf)) {
            double new[3];
            find_closest_point_on_bp(bp, vcell->vertices[ii], new);
            vcell->vertices[ii][0] = new[0];
            vcell->vertices[ii][1] = new[1];
            vcell->vertices[ii][2] = new[2];
            // if (!is_inside_convhull(new, bp->points, bp->faces, bp->fnormals, bp->Nf)) {
            //     puts("WARNING! SNAP IS STILL OUTSIDE?");
            // }
        }
    }

    // Compute volume
    compute_vcell_volume(vcell);
    return vcell;
}



s_vdiagram *voronoi_from_delaunay_3d(const s_setup *setup, s_bound_poly *bpoly, int Nreal)
{
    initialize_ncells_counter(setup);
    
    s_vdiagram *vdiagram = malloc_vdiagram(setup, Nreal);
    vdiagram->bpoly = bpoly;
    
    for (int ii=0; ii<Nreal; ii++) {
        for (int tries = 0; tries < 2; tries ++) { 
            vdiagram->vcells[ii] = extract_voronoi_cell(setup, ii, bpoly);
            if (vdiagram->vcells[ii] == NULL) continue;
            else break;
        }
        if (vdiagram->vcells[ii] == NULL || vdiagram->vcells[ii]->volume <= 0) {
            puts("OBS: Could not construct vdiagram.");
            // free_vdiagram(vdiagram);
            return NULL;
        }
    }

    return vdiagram;
}


int find_inside_which_vcell(s_vdiagram *vd, double *x)
{
    for (int ii=0; ii<vd->N_vcells; ii++) {
        if (is_inside_convhull(x, vd->vcells[ii]->vertices, vd->vcells[ii]->faces, 
            vd->vcells[ii]->fnormals, vd->vcells[ii]->Nf)) {
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


void plot_add_vcell(t_gnuplot *interface, const s_vcell *vcell, char *config)
{
    for (int ii=0; ii<vcell->Nf; ii++) {
        draw_solid_triangle_3d(interface, vcell->vertices[vcell->faces[ii*3]], 
                vcell->vertices[vcell->faces[ii*3+1]], vcell->vertices[vcell->faces[ii*3+2]], config);
        }
}


void plot_add_bpoly(t_gnuplot *interface, const s_bound_poly *bpoly, char *config)
{
    for (int ii=0; ii<bpoly->Nf; ii++) {
        draw_solid_triangle_3d(interface, bpoly->points[bpoly->faces[ii*3]], 
                bpoly->points[bpoly->faces[ii*3+1]], bpoly->points[bpoly->faces[ii*3+2]], config);
    }
}


void plot_vcell(s_vdiagram *vdiag, s_vcell *vcell, char *f_name, double *ranges)
{
    int size[2] = {1080, 1080};
    t_gnuplot *interface = gnuplot_start(PNG_3D, f_name, size, 18);
    gnuplot_config(interface,  "set pm3d depthorder",
                               "set pm3d border lc 'black' lw 0.1",
                               "set view 100, 10", 
                               "set xyplane at 0",
                               "set xlabel 'x'",
                               "set ylabel 'y'",
                               "set zlabel 'z'");
    if (ranges) {
        char buff[1024];
        snprintf(buff, 1024, "set xrange [%f:%f]\n set yrange [%f:%f]\n set zrange [%f:%f]", 
                 ranges[0], ranges[1], ranges[2], ranges[3], ranges[4], ranges[5]); 
        gnuplot_config(interface, buff);
    }

    plot_add_vcell(interface, vcell, "fs transparent solid 0.6");
    plot_add_bpoly(interface, vdiag->bpoly, "fs transparent solid 0.1");
    gnuplot_end(interface);
}


void plot_all_vcells(s_vdiagram *vdiagram, char *f_name, double *ranges, char *view_command)
{
    char buff[1024];
    char *colors[] = { "#000090", "#ee0000", "#7f0000", "#0090ff", "#0fffee", 
        "#90ff70", "#ffee00", "#000fff",  "#ff7000" };
    randomize_colors(9, colors);
    
    int size[2] = {2160, 2160};
    t_gnuplot *interface = gnuplot_start(PNG_3D, f_name, size, 18);
    gnuplot_config(interface,  "set pm3d depthorder",
                               "set pm3d border lc 'black' lw 0.01",
                               "unset border",
                               "unset xtics",
                               "unset ytics",
                               "unset ztics",
                               "set xlabel 'x'",
                               "set ylabel 'y'",
                               "set zlabel 'z'",
                               view_command);
    if (ranges) {
        snprintf(buff, 1024, "set xrange [%f:%f]\n set yrange [%f:%f]\n set zrange [%f:%f]", 
                 ranges[0], ranges[1], ranges[2], ranges[3], ranges[4], ranges[5]); 
        gnuplot_config(interface, buff);
    }


    for (int ii=0; ii<vdiagram->N_vcells; ii++) {
        s_vcell *vcell = vdiagram->vcells[ii];
        snprintf(buff, 1024, "fs transparent solid 0.2 fc rgb '%s'", colors[ii%8]);
        plot_add_vcell(interface, vcell, buff);  
    }
    gnuplot_end(interface);
}


void plot_vdiagram_differentviews(s_vdiagram *vdiagram, char *f_name, double *ranges)
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

    for (int ii=0; ii<vdiagram->N_vcells; ii++) {
        fprintf(f, "%d, %f, %f, %f, %f\n", id, vdiagram->seeds[vdiagram->vcells[ii]->seed_id][0],
                                               vdiagram->seeds[vdiagram->vcells[ii]->seed_id][1], 
                                               vdiagram->seeds[vdiagram->vcells[ii]->seed_id][2], 
                                               vdiagram->vcells[ii]->volume);   
    }
    fclose(f);
}

