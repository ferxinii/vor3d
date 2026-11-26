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


/* Serialize / deserialize
 * Format (v2): seeds | bpoly | per-cell: seed_id, volume, N_pieces, pieces[]
 * N_pieces is always 1 for convex diagrams. */
static size_t size_serialize_vdiagram(const s_vdiagram *vd)
{
    int N_seeds = vd->seeds.N;

    size_t size = 0;
    /* seeds */
    size += sizeof(int);
    size += sizeof(s_point) * N_seeds;

    /* bpoly */
    size_t size_ch = 0;
    serialize_convhull(&vd->bpoly.convh, NULL, &size_ch, NULL);
    size += size_ch;
    size += sizeof(s_point);    /* CM */
    size += sizeof(s_point);    /* min */
    size += sizeof(s_point);    /* max */
    size += sizeof(double);     /* volume */

    /* vcells */
    for (int ii=0; ii<N_seeds; ii++) {
        size += sizeof(int);    /* seed_id */
        size += sizeof(double); /* volume */
        size += sizeof(int);    /* N_pieces */
        for (int p = 0; p < vd->vcells[ii].N_pieces; p++) {
            size_t size_p = 0;
            serialize_convhull(&vd->vcells[ii].pieces[p], NULL, &size_p, NULL);
            size += size_p;
        }
    }

    return size;
}


int serialize_vdiagram(const s_vdiagram *vd, uint8_t *buff_write, size_t *size, uint8_t **out)
{
    size_t s = size_serialize_vdiagram(vd);
    if (size) *size = s;

    if (!buff_write && !out) return 1;
    if (buff_write && out) {
        fprintf(stderr, "ERROR serialize_vdiagram: either provide buffer or let function malloc, but not both.\n");
        return 0;
    }

    uint8_t *destination = NULL;
    if (buff_write) {
        destination = buff_write;
    } else {
        destination = malloc(s);
        if (!destination) return 0;
        *out = destination;
    }

    uint8_t *p = destination;

    /* seeds */
    memcpy(p, &vd->seeds.N, sizeof(int));
    p += sizeof(int);
    memcpy(p, vd->seeds.p, sizeof(s_point) * vd->seeds.N);
    p += sizeof(s_point) * vd->seeds.N;

    /* bpoly */
    size_t size_ch;
    serialize_convhull(&vd->bpoly.convh, p, &size_ch, NULL);
    p += size_ch;
    memcpy(p, &vd->bpoly.CM, sizeof(s_point));      p += sizeof(s_point);
    memcpy(p, &vd->bpoly.min, sizeof(s_point));     p += sizeof(s_point);
    memcpy(p, &vd->bpoly.max, sizeof(s_point));     p += sizeof(s_point);
    memcpy(p, &vd->bpoly.volume, sizeof(double));   p += sizeof(double);

    /* vcells */
    for (int ii=0; ii<vd->seeds.N; ii++) {
        memcpy(p, &vd->vcells[ii].seed_id,  sizeof(int));    p += sizeof(int);
        memcpy(p, &vd->vcells[ii].volume,   sizeof(double)); p += sizeof(double);
        memcpy(p, &vd->vcells[ii].N_pieces, sizeof(int));    p += sizeof(int);
        for (int pi=0; pi<vd->vcells[ii].N_pieces; pi++) {
            serialize_convhull(&vd->vcells[ii].pieces[pi], p, &size_ch, NULL);
            p += size_ch;
        }
    }

    return 1;
}

int deserialize_vdiagram(const uint8_t *data, s_vdiagram *out, size_t *bytes_read)
{
    if (!out || !data) return 0;

    s_vdiagram vd;
    const uint8_t *p = data;

    /* points */
    memcpy(&vd.seeds.N, p, sizeof(int));
    p += sizeof(int);
    assert(vd.seeds.N != 0);
    vd.seeds.p = malloc(sizeof(s_point) * vd.seeds.N);
    if (!vd.seeds.p) return 0;
    memcpy(vd.seeds.p, p, sizeof(s_point) * vd.seeds.N);
    p += sizeof(s_point) * vd.seeds.N;

    /* bpoly */
    size_t read;
    if (!deserialize_convhull(p, &vd.bpoly.convh, &read))
        { free(vd.seeds.p); return 0; };
    p += read;
    memcpy(&vd.bpoly.CM, p, sizeof(s_point));      p += sizeof(s_point);
    memcpy(&vd.bpoly.min, p, sizeof(s_point));     p += sizeof(s_point);
    memcpy(&vd.bpoly.max, p, sizeof(s_point));     p += sizeof(s_point);
    memcpy(&vd.bpoly.volume, p, sizeof(double));   p += sizeof(double);

    /* vcells */
    vd.vcells = malloc(sizeof(s_vcell) * vd.seeds.N);
    if (!vd.vcells) { free(vd.seeds.p); free_bpoly(&vd.bpoly); return 0; }
    int cells_read = 0;
    for (int ii=0; ii<vd.seeds.N; ii++) {
        memcpy(&vd.vcells[ii].seed_id,  p, sizeof(int));    p += sizeof(int);
        memcpy(&vd.vcells[ii].volume,   p, sizeof(double)); p += sizeof(double);
        memcpy(&vd.vcells[ii].N_pieces, p, sizeof(int));    p += sizeof(int);
        vd.vcells[ii].surface = NULL;          /* not serialized (yet); see plan */
        vd.vcells[ii].N_surface = 0;

        vd.vcells[ii].pieces = malloc(sizeof(s_convh) * vd.vcells[ii].N_pieces);
        if (!vd.vcells[ii].pieces) goto error;

        int pieces_read = 0;
        for (int pi=0; pi<vd.vcells[ii].N_pieces; pi++) {
            if (!deserialize_convhull(p, &vd.vcells[ii].pieces[pi], &read)) {
                for (int jj=0; jj<pieces_read; jj++) free_convhull(&vd.vcells[ii].pieces[jj]);
                free(vd.vcells[ii].pieces);
                goto error;
            }
            p += read;
            pieces_read++;
        }
        cells_read++;
    }

    *out = vd;
    if (bytes_read) *bytes_read = (size_t)(p - data);
    return 1;

error:
    for (int jj=0; jj<cells_read; jj++) free_vcell(&vd.vcells[jj]);
    free(vd.vcells);
    free(vd.seeds.p);
    free_bpoly(&vd.bpoly);
    return 0;
}

int write_serialized_vdiagram(const char *file, const uint8_t *data, size_t size)
{
	FILE *f = fopen(file, "wb");
	if (!f) return 0;

	/* write 8-byte length (little-endian) so the reader knows the size */
	uint64_t sz = (uint64_t)size;
	uint8_t header[8];
	for (int i = 0; i < 8; ++i) header[i] = (uint8_t)((sz >> (i * 8)) & 0xFF);

	if (fwrite(header, 1, 8, f) != 8) { fclose(f); return 0; }
	if (fwrite(data, 1, size, f) != size) { fclose(f); return 0; }
	fclose(f);
	return 1;
}

int read_serialized_vdiagram(const char *file, uint8_t **outbuf, size_t *outsize)
{
	FILE *f = fopen(file, "rb");
	if (!f) return 0;

	uint8_t header[8];
	if (fread(header, 1, 8, f) != 8) { fclose(f); return 0; }

	uint64_t sz = 0;
	for (int i = 0; i < 8; ++i) sz |= ((uint64_t)header[i]) << (i * 8);

	uint8_t *buf = NULL;
	if (sz > 0) {
		buf = (uint8_t*)malloc((size_t)sz);
		if (!buf) { fclose(f); return 0; }
		if (fread(buf, 1, (size_t)sz, f) != (size_t)sz)
            { free(buf); fclose(f); return 0; }
	}

	fclose(f);
	*outbuf = buf;
	*outsize = (size_t)sz;
	return 0;
}



/* Memory of vdiagram, printing, and other utilities. */
int vdiagram_is_valid(const s_vdiagram *vd)
{
    if (vd->vcells == NULL) return 0;
    else return 1;
}

void free_vdiagram(s_vdiagram *vdiagram)
{
    if (vdiagram->vcells) {
        for (int ii=0; ii<vdiagram->seeds.N; ii++)
            free_vcell(&vdiagram->vcells[ii]);
        free(vdiagram->vcells);
    }
    free_points(&vdiagram->seeds);
    free_bpoly(&vdiagram->bpoly);
    memset(vdiagram, 0, sizeof(s_vdiagram));
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
    /* s_vdiagram is always convex (N_pieces == 1) */
    for (int ii=0; ii<vd->seeds.N; ii++) {
        e_geom_test test = test_point_in_convhull(&vd->vcells[ii].pieces[0], x, EPS_degenerate, TOL);
        if (test == TEST_DEGENERATE || test == TEST_ERROR) fprintf(stderr, "find_inside_which_vcell: TEST_DEGENERATE or TEST_ERROR. Continuing.\n");
        if (test == TEST_IN || test == TEST_BOUNDARY) return ii;
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
