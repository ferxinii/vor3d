// [ ] TODO Improve malloc of vertices, or check if i have reached VCELL_MAX_VERTICES to increase size as needed
#include "points.h"
#include "convh.h"
#include "scplx.h"
#include "vdiagram.h"
#include "volsph.h"
#include "convh.h"
#include "./hlibs/dynarray.h"
#include "gnuplotc.h"
#include "./hlibs/random.h"
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>



void free_vcell(s_vcell *c)
{
    for (int i = 0; i < c->N_pieces; i++) free_convhull(&c->pieces[i]);
    free(c->pieces);
    for (int i = 0; i < c->N_surface; i++) free_trimesh(&c->surface[i]);
    free(c->surface);
    memset(c, 0, sizeof(s_vcell));
}


/* Memory of vdiagram, printing, and other utilities. */
int vdiagram_is_valid(const s_vdiagram *vd)
{
    if (vd->vcells == NULL) return 0;
    else return 1;
}

void free_vdiagram(s_vdiagram *vd)
{
    if (vd->vcells) {
        for (int i = 0; i < vd->seeds.N; i++) free_vcell(&vd->vcells[i]);
        free(vd->vcells);
    }
    free_points(&vd->seeds);
    free(vd->adj_offset);
    free(vd->adj_neighbor);
    memset(vd, 0, sizeof(s_vdiagram));
}

void print_vcell(const s_vcell *vcell)
{
    puts("VCELL");
    printf("Seed: %d, N_pieces = %d, vol = %f\n", vcell->seed_id, vcell->N_pieces, vcell->volume);
    if (vcell->N_pieces > 0) {
        printf("Piece 0: Nv = %d\n", vcell->pieces[0].points.N);
        printf("%f, %f, %f\n", vcell->pieces[0].points.p[0].x,
               vcell->pieces[0].points.p[0].y, vcell->pieces[0].points.p[0].z);
        for (int ii=1; ii<vcell->pieces[0].points.N; ii++) {
            printf("%f, %f, %f\n", vcell->pieces[0].points.p[ii].x,
                   vcell->pieces[0].points.p[ii].y, vcell->pieces[0].points.p[ii].z);
        }
    }
    printf("\n");
}

void print_vdiagram(const s_vdiagram *vdiagram)
{
    puts("------- VORONOI DIAGRAM ------");
    for (int ii=0; ii<vdiagram->seeds.N; ii++)
        print_vcell(&vdiagram->vcells[ii]);
    puts("------------------------------");
}

int find_inside_which_vcell(const s_vdiagram *vd, s_point x, double EPS_degenerate, double TOL)
{
    for (int ii=0; ii<vd->seeds.N; ii++) {
        for (int p = 0; p < vd->vcells[ii].N_pieces; p++) {
            e_geom_test test = test_point_in_convhull(&vd->vcells[ii].pieces[p], x, EPS_degenerate, TOL);
            if (test == TEST_DEGENERATE || test == TEST_ERROR) fprintf(stderr, "find_inside_which_vcell: TEST_DEGENERATE or TEST_ERROR. Continuing.\n");
            if (test == TEST_IN || test == TEST_BOUNDARY) return ii;
        }
    }
    return -1;
}


// --------------------------------------- PLOTS ------------------------------------------------
static void randomize_colors(int N, char **colors)
{
    /* Plot-only cosmetic shuffle. Uses a locally-seeded PRNG (not the process-global
     * rand()) so it neither depends on nor perturbs any shared RNG state. */
    s_random_context rng = random_initialize(0x9E3779B97F4A7C15ULL);
    for (int i = N - 1; i > 0; i--) {
        int j = (int)random_uniform_range_u64(&rng, (uint64_t)(i + 1));
        char *tmp = colors[i];
        colors[i] = colors[j];
        colors[j] = tmp;
    }
}

static void plot_add_vcell(s_gnuplot *interface, const s_vcell *vcell, char *config)
{
    for (int p = 0; p < vcell->N_pieces; p++) {
        const s_convh *ch = &vcell->pieces[p];
        for (int ii=0; ii<ch->Nf; ii++) {
            s_point face[3];
            convh_get_face(ch, ii, face);
            draw_solid_triangle_3d(interface, face[0].coords, face[1].coords, face[2].coords, config);
        }
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


void plot_vcell(const s_vor3d_domain *domain, const s_vcell *vcell, char *f_name, const s_point ranges[2])
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
    plot_add_bpoly(interface, &domain->bpoly, "fs transparent solid 0.1");
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
