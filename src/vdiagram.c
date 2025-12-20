// [ ] TODO Improve malloc of vertices, or check if i have reached VCELL_MAX_VERTICES to increase size as needed
#include "points.h"
#include "convh.h"
#include "scplx.h"
#include "bpoly.h"
#include "vdiagram.h"
#include "convh.h"
#include "bpoly.h"
#include "lists.h"
#include "gnuplotc.h"
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>


static int id_where_equal(const int *arr1, int N, int a)
{
    for (int ii=0; ii<N; ii++) if (arr1[ii] == a) return ii;
    return -1;
}


/* Serialize / deserialize */
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
    size += sizeof(double);     /* dmax */
    size += sizeof(s_point);    /* CM */
    size += sizeof(s_point);    /* min */
    size += sizeof(s_point);    /* max */
    size += sizeof(double);     /* volume */

    /* vcells */
    size += sizeof(int) * N_seeds;      /* seed_id */
    size += sizeof(double) * N_seeds;   /* volume */
    for (int ii=0; ii<N_seeds; ii++) {
        serialize_convhull(&vd->vcells[ii].convh, NULL, &size_ch, NULL);
        size += size_ch;
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
    memcpy(p, &vd->bpoly.dmax, sizeof(double));     p += sizeof(double);     
    memcpy(p, &vd->bpoly.CM, sizeof(s_point));      p += sizeof(s_point);    
    memcpy(p, &vd->bpoly.min, sizeof(s_point));     p += sizeof(s_point);    
    memcpy(p, &vd->bpoly.max, sizeof(s_point));     p += sizeof(s_point);    
    memcpy(p, &vd->bpoly.volume, sizeof(double));   p += sizeof(double);     
    
    /* vcells */
    for (int ii=0; ii<vd->seeds.N; ii++) {
        memcpy(p, &vd->vcells[ii].seed_id, sizeof(int));
        p += sizeof(int);
        memcpy(p, &vd->vcells[ii].volume, sizeof(double));
        p += sizeof(double);

        serialize_convhull(&vd->vcells[ii].convh, p, &size_ch, NULL);
        p += size_ch;
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
    memcpy(&vd.bpoly.dmax, p, sizeof(double));     p += sizeof(double);     
    memcpy(&vd.bpoly.CM, p, sizeof(s_point));      p += sizeof(s_point);    
    memcpy(&vd.bpoly.min, p, sizeof(s_point));     p += sizeof(s_point);    
    memcpy(&vd.bpoly.max, p, sizeof(s_point));     p += sizeof(s_point);    
    memcpy(&vd.bpoly.volume, p, sizeof(double));   p += sizeof(double);     
    
    /* vcells */
    vd.vcells = malloc(sizeof(s_vcell) * vd.seeds.N);
    if (!vd.vcells) { free(vd.seeds.p); free_bpoly(&vd.bpoly); return 0; }
    for (int ii=0; ii<vd.seeds.N; ii++) {
        memcpy(&vd.vcells[ii].seed_id, p, sizeof(int));
        p += sizeof(int);
        memcpy(&vd.vcells[ii].volume, p, sizeof(double));
        p += sizeof(double);

        if (!deserialize_convhull(p, &vd.vcells[ii].convh, &read)) { 
            for (int jj=0; jj<ii; jj++) free_convhull(&vd.vcells[jj].convh);
            free(vd.vcells);
            free(vd.seeds.p);
            free_bpoly(&vd.bpoly);
            return 0;
        }
        p += read;
    }

    *out = vd;
    if (bytes_read) *bytes_read = (size_t)(p - data);
    return 1;
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

static s_vdiagram malloc_vdiagram(const s_scplx *setup, int Nreal)
{
    assert(Nreal <= setup->points.N && "Nreal needs to be <= than setup->N"); 

    s_vdiagram out;
    memset(&out.bpoly, 0, sizeof(s_bpoly));

    // setup->points may have more points... Need to copy manually.
    out.seeds = (s_points) { .N = Nreal, .p = malloc(sizeof(s_point)*Nreal) };
    memcpy(out.seeds.p, setup->points.p, sizeof(s_point) * Nreal);
    
    out.vcells = malloc(sizeof(s_vcell) * Nreal);
    memset(out.vcells, 0, sizeof(s_vcell) * Nreal);
    return out;
}

void free_vdiagram(s_vdiagram *vdiagram)
{
    for (int ii=0; ii<vdiagram->seeds.N; ii++)
        if (vdiagram->vcells[ii].convh.Nf != 0) free_convhull(&vdiagram->vcells[ii].convh);
    free(vdiagram->vcells);
    free_points(&vdiagram->seeds);
    free_bpoly(&vdiagram->bpoly);
    memset(vdiagram, 0, sizeof(s_vdiagram));
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

int find_inside_which_vcell(const s_vdiagram *vd, s_point x, double EPS_degenerate, double TOL)
{
    for (int ii=0; ii<vd->seeds.N; ii++) {
        e_geom_test test = test_point_in_convhull(&vd->vcells[ii].convh, x, EPS_degenerate, TOL);
        if (test == TEST_DEGENERATE || test == TEST_ERROR) fprintf(stderr, "find_inside_which_vcell: TEST_DEGENERATE or TEST_ERROR. Continuing.\n");
        if (test == TEST_IN || test == TEST_BOUNDARY) return ii;
    }
    return -1;
}



/* Main functions to construct vdiagram */
static int add_vvertex_from_ncell(const s_scplx *setup, const s_ncell *ncell, double EPS_degenerate, double TOL,  s_list *vertices)
{   /* 1: added, 0: not added, -1: error */
    s_point v[4];
    extract_vertices_ncell(setup, ncell, v);
    assert(orientation_robust(v, v[3]) != 0);

    // s_point a = subtract_points(v[1], v[0]);
    // s_point b = subtract_points(v[2], v[0]);
    // s_point c = subtract_points(v[3], v[0]);
    //
    // double a2 = norm_squared(a);
    // double b2 = norm_squared(b);
    // double c2 = norm_squared(c);
    // 
    // s_point bc = cross_prod(b, c);
    // s_point ca = cross_prod(c, a);
    // s_point ab = cross_prod(a, b);
    //
    // double f = 2 * dot_prod(ab, c);
    // if (fabs(f) < EPS_degenerate) goto error_circumsphere;
    // f = 1.0 / f;
    //
    // s_point circumsphere = {{{
    //     f * ( a2 * bc.x + b2 * ca.x + c2 * ab.x) + v[0].x,
    //     f * ( a2 * bc.y + b2 * ca.y + c2 * ab.y) + v[0].y,
    //     f * ( a2 * bc.z + b2 * ca.z + c2 * ab.z) + v[0].z,
    // }}};
    //

    s_point circumcentre;
    if (!circumcentre_from_points(v, EPS_degenerate, &circumcentre)) goto error_circumsphere;

    /* Check if vertex already exists */
    const double TOL2 = TOL*TOL;
    for (unsigned ii=0; ii<vertices->N; ii++) {  
        s_point p; 
        if (!list_get_value(vertices, ii, &p)) goto error_list;
        if (distance_squared(circumcentre, p) < TOL2) return 0;
    }

    if (!list_push(vertices, &circumcentre)) goto error_list;
    return 1;

    error_circumsphere:
        fprintf(stderr, "Error in add_vvertex_from_ncell: Could not compute circumsphere.\n");
        /* printf("DEBUG: EPS = %g, f = %f\n", EPS_degenerate, f); */
        // printf("%d: %g, %g, %g\n", ncell->vertex_id[0], v[0].x, v[0].y, v[0].z);
        // printf("%d: %g, %g, %g\n", ncell->vertex_id[1], v[1].x, v[1].y, v[1].z);
        // printf("%d: %g, %g, %g\n", ncell->vertex_id[2], v[2].x, v[2].y, v[2].z);
        // printf("%d: %g, %g, %g\n", ncell->vertex_id[3], v[3].x, v[3].y, v[3].z);
        // printf("vol: %g\n", signed_volume_tetra(v));
        // printf("orient3d: %d\n", orientation_robust(v, v[3]));
        return 0;
    error_list:
        fprintf(stderr, "Error in add_vvertex_from_ncell. List problem.\n");
        return -1;
}

static s_convh convhull_from_marked_ncells(const s_scplx *setup, const s_bpoly *bp, double EPS_degenerate, double TOL, s_list *buff_v)
{
    s_list *vertices = buff_v;

    vertices->N = 0;
    s_ncell *current = setup->head;
    for (int ii=0; ii<setup->N_ncells; ii++, current = current->next) {
        if (current->mark) 
            if (add_vvertex_from_ncell(setup, current, EPS_degenerate, TOL, vertices) == -1) goto error;
    }

    /* Snap point back into bounding polyhedron if outside. 
     * This is a limitation of mirroring approach used to bound vdiagram? Or something deeper? */
    for (unsigned ii=0; ii<vertices->N; ii++) {
        s_point p;
        if (!list_get_value(vertices, ii, &p)) goto error;
        if (test_point_in_convhull(&bp->convh, p, EPS_degenerate, TOL) == TEST_OUT) {
            s_point tmp;
            double distance = find_closest_point_on_bp(bp, p, EPS_degenerate, &tmp);
            if (distance > TOL) {
                // printf("DEBUG: Snap! %g\n", distance);
                if (!list_change_entry(vertices, ii, &tmp)) goto error;
            }
        }
    }

    s_points points = {.N = vertices->N, .p = vertices->items};  /* No need to free, external buffer */
    s_convh ch;
    if (convhull_from_points(&points, EPS_degenerate, TOL, &ch) != 1) goto error;
    else return ch;

    error:
        fprintf(stderr, "Error in convhull_from_marked_ncells.\n");
        return convhull_NAN;
}

static s_vcell extract_voronoi_cell(const s_scplx *setup, const s_bpoly *bp, int seed_id, double EPS_degenerate, double TOL, s_list *buff_v)
{
    /* Find an ncell with this vertex */
    s_ncell *ncell = setup->head;
    int v_localid = id_where_equal(ncell->vertex_id, 4, seed_id);
    while (v_localid == -1) { 
        ncell = ncell->next;
        assert(ncell);
        v_localid = id_where_equal(ncell->vertex_id, 4, seed_id);
    }

    /* Find "complementary" indices to localid */
    int v_localid_COMP[3];
    for (int ii=0, kk=0; ii<4; ii++) 
        if (ii != v_localid) v_localid_COMP[kk++] = ii;

    /* Mark incident ncells to this point */
    mark_ncells_incident_face(setup, ncell, 0, v_localid_COMP);

    s_vcell out = {0};
    out.seed_id = seed_id;
    out.convh = convhull_from_marked_ncells(setup, bp, EPS_degenerate, TOL, buff_v);
    if (!convhull_is_valid(&out.convh)) {
        fprintf(stderr, "extract_voronoi_cell: Error extracting convhull of cell.\n");
        return (s_vcell){0};
    }
    out.volume = volume_convhull(&out.convh);
    return out;
}


s_vdiagram voronoi_from_delaunay_3d(const s_scplx *setup, const s_bpoly *bpoly, int Nreal, double EPS_degenerate, double TOL)
{   /* copy of bpoly inside */
    initialize_ncells_counter(setup);

    s_vdiagram out = malloc_vdiagram(setup, Nreal);
    out.bpoly = bpoly_copy(bpoly);
    
    s_list buff_v = list_initialize(sizeof(s_point), 10);  /* Place to store vertices of each cell before convhull */
    for (int ii=0; ii<Nreal; ii++) {
        out.vcells[ii] = extract_voronoi_cell(setup, bpoly, ii, EPS_degenerate, TOL, &buff_v);
        if (!convhull_is_valid(&out.vcells[ii].convh)) {
            fprintf(stderr, "voronoi_from_delaunay_3d: Could not construct vdiagram.\n");
            goto error;
        }
    }

    free_list(&buff_v);
    return out;

    error:
        free_vdiagram(&out);
        free_list(&buff_v);
        return out;
}





// --------------------------------------- PLOTS ------------------------------------------------
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

