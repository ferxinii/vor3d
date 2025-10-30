// [ ] TODO Improve malloc of vertices, or check if i have reached VCELL_MAX_VERTICES to increase size as needed
#include "geometry.h"
#include "convh.h"
#include "scplx.h"
#include "bpoly.h"
#include "algebra.h"
#include "geometry.h"
#include "vdiagram.h"
#include "convh.h"
#include "bpoly.h"
#include "gnuplotc.h"
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>


typedef struct vertex_buff {
    int max;
    int N;
    s_point *v;
} s_vbuff;


static s_vbuff init_vbuff()
{
    s_vbuff vbuff;
    vbuff.max = VCELL_BLOCK_VERTICES;
    vbuff.N = 0;
    vbuff.v = malloc(sizeof(s_point) * VCELL_BLOCK_VERTICES);
    return vbuff;
}


static void free_vbuff(s_vbuff *vbuff)
{
    free(vbuff->v);
    memset(vbuff, 0, sizeof(s_vbuff));
}


static void increase_vbuff_if_needed(s_vbuff *vbuff)
{
    if (vbuff->N + 1 >= vbuff->max) {
        vbuff->max += VCELL_BLOCK_VERTICES;
        vbuff->v = realloc(vbuff->v, sizeof(s_point) * vbuff->max);
    }
}


static int add_vvertex_from_ncell(const s_scplx *setup, const s_ncell *ncell, s_vbuff *vbuff)
{   // Returns the index of the vertex
    increase_vbuff_if_needed(vbuff);

    s_point v[4];
    extract_vertices_ncell(setup, ncell, v);

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
        return EXIT_FAILURE;
    }
    f = 1.0 / f;

    s_point circumcenter = {{{
        f * ( a2 * bc.x + b2 * ca.x + c2 * ab.x) + v[0].x,
        f * ( a2 * bc.y + b2 * ca.y + c2 * ab.y) + v[0].y,
        f * ( a2 * bc.z + b2 * ca.z + c2 * ab.z) + v[0].z,
    }}};

    for (int ii=0; ii<vbuff->N; ii++) {
        if (distance_squared(circumcenter, vbuff->v[ii]) < 1e-6) {
            return EXIT_FAILURE;
        }
    }

    vbuff->v[vbuff->N++] = circumcenter;

    return EXIT_SUCCESS;
}


static s_convhull convhull_from_vbuff(s_vbuff *vbuff)
{
    s_points points = {.N = vbuff->N, .p = vbuff->v};
    // Makes copy of points inside.
    return convhull_from_points(&points);
}


static s_vdiagram malloc_vdiagram(const s_scplx *setup, int Nreal)
{
    assert(Nreal <= setup->points.N && "Nreal needs to be <= than setup->N"); 

    s_vdiagram out;

    // setup->points may have more points... Need to copy manually.
    out.seeds = (s_points) { .N = Nreal,
                             .p = malloc(sizeof(s_point) * Nreal) };
    memcpy(out.seeds.p, setup->points.p, sizeof(s_point) * Nreal);
    
    out.vcells = malloc(sizeof(s_vcell) * Nreal);
    memset(out.vcells, 0, sizeof(s_vcell) * Nreal);

    memset(&out.bpoly, 0, sizeof(s_bpoly));
    return out;
}


void free_vcell(s_vcell *vcell)
{
    free_convhull(&vcell->convh);
    memset(vcell, 0, sizeof(s_vcell));
}


void free_vdiagram(s_vdiagram *vdiagram)
{
    for (int ii=0; ii<vdiagram->seeds.N; ii++) {
        if (vdiagram->vcells[ii].seed_id != 0) free_vcell(&vdiagram->vcells[ii]);
    }
    free(vdiagram->vcells);
    free_points(&vdiagram->seeds);
    free_bpoly(&vdiagram->bpoly);
    memset(vdiagram, 0, sizeof(s_vdiagram));
}


void write_vcell_file(const s_vcell *vcell, FILE *file)
{
    for (int ii=0; ii<vcell->convh.Nf; ii++) {
        s_point face[3];
        convh_get_face(&vcell->convh, ii, face);
        fprintf(file, "%f %f %f\n", face[0].x, face[0].y, face[0].z);
        fprintf(file, "%f %f %f\n", face[1].x, face[1].y, face[1].z);
        fprintf(file, "%f %f %f\n\n", face[2].x, face[2].y, face[2].z);
    }
}


void write_vd_file(const s_vdiagram *vd, FILE *file)
{
    for (int ii=0; ii<vd->seeds.N; ii++) {
        write_vcell_file(&vd->vcells[ii], file);
        fprintf(file, "\n\n");
    }
}


void print_vcell(const s_vcell *vcell)
{
    puts("VCELL");
    printf("Seed: %d, Nv = %d, vol = %f\n", vcell->seed_id, vcell->convh.points.N, vcell->volume);
    printf("%f, %f, %f\n", vcell->convh.points.p[0].x, vcell->convh.points.p[0].y, vcell->convh.points.p[0].z);
    for (int ii=1; ii<vcell->convh.points.N; ii++) {
        printf("%f, %f, %f\n", vcell->convh.points.p[ii].x, vcell->convh.points.p[ii].y, vcell->convh.points.p[ii].z);
    }
    printf("\n");
}


void print_vdiagram(const s_vdiagram *vdiagram)
{
    puts("------- VORONOI DIAGRAM ------");
    for (int ii=0; ii<vdiagram->seeds.N; ii++) {
        print_vcell(&vdiagram->vcells[ii]);
    }
    puts("------------------------------");
}


static s_convhull convhull_from_marked_ncells(const s_scplx *setup, const s_bpoly *bp)
{
    s_vbuff vbuff = init_vbuff();

    s_ncell *current = setup->head;
    while (current) {
        if (current->mark == 1) add_vvertex_from_ncell(setup, current, &vbuff);
        current = current->next;
    }

    // Snap into bp if cell spikes outward, KNOWN PROBLEM FIXME TODO
    for (int ii=0; ii<vbuff.N; ii++) {
        if (is_inside_convhull(&bp->convh, vbuff.v[ii]) == 0) {
            s_point tmp;
            double distance = find_closest_point_on_bp(bp, vbuff.v[ii], &tmp);
            if (distance > 1e-12) vbuff.v[ii] = tmp;
        }
    }

    s_convhull out = convhull_from_vbuff(&vbuff);

    free_vbuff(&vbuff);
    return out;
}


static s_vcell extract_voronoi_cell(const s_scplx *setup, int vertex_id, const s_bpoly *bp)
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
    mark_ncells_incident_face(setup, ncell, 0, v_localid_COMP);


    s_vcell out;
    out.seed_id = vertex_id;
    out.convh = convhull_from_marked_ncells(setup, bp);
    if (out.convh.Nf == 0) {
        fprintf(stderr, "extract_voronoi_cell: Error extracting convhull of cell\n");
        memset(&out, 0, sizeof(s_vcell));
        return out;
    }

    out.volume = volume_convhull(&out.convh);
    return out;
}



s_vdiagram voronoi_from_delaunay_3d(const s_scplx *setup, const s_bpoly *bpoly, int Nreal)
{   // copy of bpoly
    initialize_ncells_counter(setup);
    
    s_vdiagram out = malloc_vdiagram(setup, Nreal);
    out.bpoly = bpoly_copy(bpoly);
    
    for (int ii=0; ii<Nreal; ii++) {
        for (int tries = 0; tries < 2; tries ++) { 
            out.vcells[ii] = extract_voronoi_cell(setup, ii, bpoly);
            if (out.vcells[ii].convh.Nf == 0) continue;
            else break;
        }
        if (out.vcells[ii].convh.Nf == 0 || out.vcells[ii].volume <= 0) {
            // printf("%p, %f\n", (void*)vdiagram->vcells[ii], vdiagram->vcells[ii]->volume);
            // print_vcell(vdiagram->vcells[ii]);
            puts("ERROR: Could not construct vdiagram. Exiting voronoi_from_delaunay_3d");
            free_vdiagram(&out);
            return out;
        }
    }

    return out;
}


int find_inside_which_vcell(const s_vdiagram *vd, s_point x)
{
    for (int ii=0; ii<vd->seeds.N; ii++) {
        if (is_inside_convhull(&vd->vcells[ii].convh, x)) return ii;
    }
    return -1;
}


// ----------------------------------------------------------------------------------------------
// --------------------------------------- PLOTS ------------------------------------------------
// ----------------------------------------------------------------------------------------------

static void randomize_colors(int N, char **colors) 
{
    for (int i = N - 1; i > 0; i--) {
        int j = rand() % (i + 1);
        char *tmp = colors[i];
        colors[i] = colors[j];
        colors[j] = tmp;
    }
}


static void plot_add_vcell(s_gnuplot *interface, const s_vcell *vcell, char *config)
{
    for (int ii=0; ii<vcell->convh.Nf; ii++) {
        s_point face[3];
        convh_get_face(&vcell->convh, ii, face);
        draw_solid_triangle_3d(interface, face[0].coords, face[1].coords, face[2].coords, config);

        }
}


static void plot_add_bpoly(s_gnuplot *interface, const s_bpoly *bpoly, char *config)
{
    for (int ii=0; ii<bpoly->convh.Nf; ii++) {
        s_point face[3];
        convh_get_face(&bpoly->convh, ii, face);
        draw_solid_triangle_3d(interface, face[0].coords, face[1].coords, face[2].coords, config);
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
    plot_add_bpoly(interface, &vdiag->bpoly, "fs transparent solid 0.1");
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


    for (int ii=0; ii<vdiagram->seeds.N; ii++) {
        snprintf(buff, 1024, "fs transparent solid 0.2 fc rgb '%s'", colors[ii%8]);
        plot_add_vcell(interface, &vdiagram->vcells[ii], buff);  
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

    for (int ii=0; ii<vdiagram->seeds.N; ii++) {
        fprintf(f, "%d, %f, %f, %f, %f\n", id, vdiagram->seeds.p[vdiagram->vcells[ii].seed_id].x,
                                               vdiagram->seeds.p[vdiagram->vcells[ii].seed_id].y, 
                                               vdiagram->seeds.p[vdiagram->vcells[ii].seed_id].z, 
                                               vdiagram->vcells[ii].volume);   
    }
    fclose(f);
}

