// [ ] TODO Improve malloc of vertices, or check if i have reached VCELL_MAX_VERTICES to increase size as needed
#ifndef VD_3D_C
#define VD_3D_C

#include "simplical_complex.c"
#include <float.h>
#define CONVHULL_3D_ENABLE
#include "convhull_3d.h"
#include "bpoly.c"

#define VCELL_BLOCK_VERTICES 1000


typedef struct vdiagram {
    int N_vcells;
    struct vcell **vcells;  // Array of pointers to the cells, N_vcells x 1
    const struct bound_poly *bpoly;
    double **seeds;
} s_vdiagram;


typedef struct vcell {
    int seed_id;
    // struct vcell *next;     // linked list of vcells
    int Nv;
    int Nv_capacity;
    double **vertices;  // Nv x 3
    int **origin_vertices;  // Nv x 4, LAST column indicates the dual ncell if POSITIVE,
                            // if -1: comes from circumcenter delaunay. The rest of 
                            //        the columns indicate the delaunay indices of the 
                            //        face whose normal was extended, 
                            // if -2: ARTIFICIAL extension, first column is id of face of
                            //        convhull of setup points
                            // if -3: it is from the bounding polyhedron, and the first 
                            //        index is the point id 
                            // if -4: it comes from some coords, supposedly from the ray 
                            //        intersection with bounding convex hull, and the first
                            //        columns correspond to face's vertex_id
    int Nf;
    int *faces;
    double **fnormals;
    double volume;
} s_vcell;


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
        free_vcell(vdiagram->vcells[ii]);
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


s_vdiagram *malloc_vdiagram(const s_setup *setup)
{
    s_vdiagram *out = malloc(sizeof(s_vdiagram));
    out->seeds = malloc_matrix(setup->N_points, 3);
    copy_matrix(setup->points, out->seeds, setup->N_points, 3);

    out->N_vcells = setup->N_points;
    out->vcells = malloc(sizeof(s_vcell*) * setup->N_points);
    out->bpoly = NULL;
    return out;
}


void print_vcell(const s_vcell *vcell)
{
    puts("VCELL");
    printf("Seed: %d\n", vcell->seed_id);
    printf("Nv = %d\n", vcell->Nv);
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

    // First check if the vertex already exists
    for (int ii=0; ii<vcell->Nv; ii++) {
        if (vcell->origin_vertices[ii][3] == ncell->count) {
            puts("DEBUG ADD_VVERTEX_FROM_NCELL: ALREADY EXISTS!");
            return ii;
        }
    }

    static double **vertices_ncell = NULL;
    if (!vertices_ncell) vertices_ncell = malloc_matrix(4, 3);

    extract_vertices_ncell(setup, ncell, vertices_ncell);

    double CM[3];
    find_center_mass(vertices_ncell, 4, 3, CM);

    vcell->vertices[vcell->Nv][0] = CM[0];
    vcell->vertices[vcell->Nv][1] = CM[1];
    vcell->vertices[vcell->Nv][2] = CM[2];

    vcell->origin_vertices[vcell->Nv][3] = ncell->count;
    vcell->Nv++;

    return vcell->Nv - 1;
}


int add_vvertex_from_coords(double *coords, int *dt_face_vid, s_vcell *vcell)
{   // This function assumes that the vertex does not already exist!
    increase_num_vertices_if_needed(vcell);

    vcell->vertices[vcell->Nv][0] = coords[0];
    vcell->vertices[vcell->Nv][1] = coords[1];
    vcell->vertices[vcell->Nv][2] = coords[2];

    vcell->origin_vertices[vcell->Nv][3] = -4;
    vcell->origin_vertices[vcell->Nv][0] = dt_face_vid[0];
    vcell->origin_vertices[vcell->Nv][1] = dt_face_vid[1];
    vcell->origin_vertices[vcell->Nv][2] = dt_face_vid[2];

    vcell->Nv++;

    return vcell->Nv - 1;
}


int add_vvertex_from_segment(double *i1, int bp_vid_1, int bp_vid_2, int face_id, s_vcell *vcell)
{   // i1 has the coordinates of the intersection!
    increase_num_vertices_if_needed(vcell);

    vcell->vertices[vcell->Nv][0] = i1[0];
    vcell->vertices[vcell->Nv][1] = i1[1];
    vcell->vertices[vcell->Nv][2] = i1[2];
    vcell->origin_vertices[vcell->Nv][3] = -4;
    vcell->origin_vertices[vcell->Nv][0] = bp_vid_1;
    vcell->origin_vertices[vcell->Nv][1] = bp_vid_2;
    vcell->origin_vertices[vcell->Nv][2] = face_id;

    vcell->Nv++;

    return vcell->Nv - 1;
}


int segment_vvertex_already_exists(const s_vcell *vcell, int p1, int p2)
{
    for (int ii=0; ii<vcell->Nv; ii++) {
        if (vcell->origin_vertices[ii][3] == -4 &&
            inarray(vcell->origin_vertices[ii], 2, p1) && 
            inarray(vcell->origin_vertices[ii], 2, p2)) {
            return 1;
        }
    }
    return 0;
}


int add_vvertex_as_extension(double *coords, int *dt_face_vid, s_vcell *vcell)
{   // This function assumes that the vertex does not already exist!
    increase_num_vertices_if_needed(vcell);

    vcell->vertices[vcell->Nv][0] = coords[0];
    vcell->vertices[vcell->Nv][1] = coords[1];
    vcell->vertices[vcell->Nv][2] = coords[2];

    vcell->origin_vertices[vcell->Nv][3] = -2;
    vcell->origin_vertices[vcell->Nv][0] = dt_face_vid[0];
    vcell->origin_vertices[vcell->Nv][1] = dt_face_vid[1];
    vcell->origin_vertices[vcell->Nv][2] = dt_face_vid[2];

    vcell->Nv++;

    return vcell->Nv - 1;
}


int add_vvertex_from_bpoly(const s_bound_poly *bp, int id, s_vcell *vcell)
{
    increase_num_vertices_if_needed(vcell);

    vcell->vertices[vcell->Nv][0] = bp->points[id][0];
    vcell->vertices[vcell->Nv][1] = bp->points[id][1];
    vcell->vertices[vcell->Nv][2] = bp->points[id][2];

    vcell->origin_vertices[vcell->Nv][3] = -3;
    vcell->origin_vertices[vcell->Nv][0] = id;

    vcell->Nv++;

    return vcell->Nv - 1;
}


void remove_artificial_extended_vertices(s_vcell *vcell)
{
    int kk = 0;
    for (int ii=0; ii<vcell->Nv; ii++) {
        if (vcell->origin_vertices[ii][3] != -2) {
            copy_matrix(&vcell->vertices[ii], &vcell->vertices[ii-kk], 1, 3);
            copy_matrix_int(&vcell->origin_vertices[ii], &vcell->origin_vertices[ii-kk], 1, 4);
        } else {
            kk++;
        }
    }
    vcell->Nv = vcell->Nv - kk;
}


void correctly_bound_with_bp_vertices(const s_vdiagram *vdiagram, s_vcell *vcell)
{   
    const s_bound_poly *bp = vdiagram->bpoly;
    for (int ii=0; ii<bp->Np; ii++) {
        if (is_inside_convhull(bp->points[ii], vcell->vertices, vcell->faces, vcell->fnormals, vcell->Nf)) {
            add_vvertex_from_bpoly(bp, ii, vcell);
        }
    }

    for (int ii=0; ii<bp->Nf; ii++) {
        // CHECK SEGMENT CONVHULL INTERSECTION
        double *p0, *p1, i1[3], i2[3]; 
        int ind1, ind2;
        if (!segment_vvertex_already_exists(vcell, bp->faces[ii*3 + 0], bp->faces[ii*3 + 1])) {
            p0 = bp->points[bp->faces[ii*3 + 0]];
            p1 = bp->points[bp->faces[ii*3 + 1]];
            segment_convex_hull_intersection(p0, p1, vcell->vertices, vcell->faces, 
                                                 vcell->fnormals, vcell->Nf, &ind1, i1, &ind2, i2);
            if (ind1) add_vvertex_from_segment(i1, bp->faces[ii*3+0], bp->faces[ii*3+1], ii, vcell);
            if (ind2) add_vvertex_from_segment(i2, bp->faces[ii*3+0], bp->faces[ii*3+1], ii, vcell);
        }

        if (!segment_vvertex_already_exists(vcell, bp->faces[ii*3 + 0], bp->faces[ii*3 + 2])) {
            p0 = bp->points[bp->faces[ii*3 + 0]];
            p1 = bp->points[bp->faces[ii*3 + 2]];
            segment_convex_hull_intersection(p0, p1, vcell->vertices, vcell->faces, 
                                                 vcell->fnormals, vcell->Nf, &ind1, i1, &ind2, i2);
            if (ind1) add_vvertex_from_segment(i1, bp->faces[ii*3+0], bp->faces[ii*3+2], ii, vcell);
            if (ind2) add_vvertex_from_segment(i2, bp->faces[ii*3+0], bp->faces[ii*3+2], ii, vcell);
        }

        if (!segment_vvertex_already_exists(vcell, bp->faces[ii*3 + 2], bp->faces[ii*3 + 1])) {
            p0 = bp->points[bp->faces[ii*3 + 2]];
            p1 = bp->points[bp->faces[ii*3 + 1]];
            segment_convex_hull_intersection(p0, p1, vcell->vertices, vcell->faces, 
                                                 vcell->fnormals, vcell->Nf, &ind1, i1, &ind2, i2);
            if (ind1) add_vvertex_from_segment(i1, bp->faces[ii*3+2], bp->faces[ii*3+1], ii, vcell);
            if (ind2) add_vvertex_from_segment(i2, bp->faces[ii*3+2], bp->faces[ii*3+1], ii, vcell);
        }
    }


    remove_artificial_extended_vertices(vcell);
}


void compute_vcell_volume(s_vcell *vcell)
{   // TODO??
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


void add_vvertices_extended_normals(const s_setup *setup, int *faces_ch_setup, 
        double **fnormals, int Nf_ch_setup, int vertex_id, const s_vdiagram *vd, s_vcell *vcell)
{
    for (int ii = 0; ii<Nf_ch_setup; ii++) {
        if (inarray(&faces_ch_setup[ii*3], 3, vertex_id)) {
            // printf("DEBUG EXTENDED: INARRAY! %d, %d, %d; %d\n", faces_ch_setup[ii*3], faces_ch_setup[ii*3 + 1], faces_ch_setup[ii*3+2], vertex_id);
            
            double o[3];
            o[0] = (setup->points[faces_ch_setup[ii*3 + 0]][0] + 
                    setup->points[faces_ch_setup[ii*3 + 1]][0] + 
                    setup->points[faces_ch_setup[ii*3 + 2]][0]) / 3;
            o[1] = (setup->points[faces_ch_setup[ii*3 + 0]][1] + 
                    setup->points[faces_ch_setup[ii*3 + 1]][1] + 
                    setup->points[faces_ch_setup[ii*3 + 2]][1]) / 3;
            o[2] = (setup->points[faces_ch_setup[ii*3 + 0]][2] + 
                    setup->points[faces_ch_setup[ii*3 + 1]][2] + 
                    setup->points[faces_ch_setup[ii*3 + 2]][2]) / 3;

            double s = vd->bpoly->dmax;

            double v_coords[3] = {o[0] + s * fnormals[ii][0],
                                  o[1] + s * fnormals[ii][1],
                                  o[2] + s * fnormals[ii][2]};
            add_vvertex_as_extension(v_coords, &faces_ch_setup[ii*3], vcell);

            double intersection[3];
            inside_ray_convhull_intersection(vd->bpoly->points, vd->bpoly->faces, vd->bpoly->fnormals, vd->bpoly->Nf, o, fnormals[ii], intersection);
            add_vvertex_from_coords(intersection, &faces_ch_setup[ii*3], vcell);
        }
    }
}


void add_convex_hull_vcell(s_vcell *vcell)
{
    double CM[3];
    find_center_mass(vcell->vertices, vcell->Nv, 3, CM);
    int *faces, N_faces;
    ch_vertex *vertices_aux = convert_points_to_chvertex(vcell->vertices, vcell->Nv);
    convhull_3d_build(vertices_aux, vcell->Nv, &faces, &N_faces);
    vcell->faces = faces;
    vcell->Nf = N_faces;
    double **fnormals = extract_normals_from_ch(vertices_aux, faces, N_faces, CM);
    vcell->fnormals = fnormals;
    free(vertices_aux);
}


void bounded_extraction(const s_setup *setup, s_vcell *vcell)
{
    s_ncell *current = setup->head;
    while (current) {
        if (current->mark == 1) add_vvertex_from_ncell(setup, current, vcell);
        current = current->next;
    }
    
    add_convex_hull_vcell(vcell);
}


void unbounded_extraction(const s_setup *setup, const s_vdiagram *vdiagram, s_vcell *vcell, 
        int *faces_ch_setup, double **fnormals_ch_setup, int Nf_ch_setup, int vertex_id)
{
    s_ncell *current = setup->head;
    while (current) {
        if (current->mark == 1) {
            add_vvertex_from_ncell(setup, current, vcell);
        }
        current = current->next;
    }
    add_vvertices_extended_normals(setup, faces_ch_setup, fnormals_ch_setup, Nf_ch_setup, vertex_id, vdiagram, vcell);

    // First build convex hull of this extended vcell
    add_convex_hull_vcell(vcell);
    
    // Fix bounding box with bounding polyhedra and redo triangulation
    correctly_bound_with_bp_vertices(vdiagram, vcell);
    free_matrix(vcell->fnormals, vcell->Nf);
    add_convex_hull_vcell(vcell);
}


s_vcell *extract_voronoi_cell(const s_vdiagram *vdiagram, const s_setup *setup, 
                              int *faces_ch_setup, double **fnormals_ch_setup,
                              int Nf_ch_setup, int vertex_id)
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


    s_vcell *vcell = malloc_vcell(vertex_id);
    if (is_in_boundary_convhull(faces_ch_setup, Nf_ch_setup, vertex_id)) {
        unbounded_extraction(setup, vdiagram, vcell, faces_ch_setup, fnormals_ch_setup, Nf_ch_setup, vertex_id);
    } else {
        bounded_extraction(setup, vcell);
    }

    // Compute volume
    compute_vcell_volume(vcell);
    return vcell;
}



s_vdiagram *voronoi_from_delaunay_3d(const s_setup *setup, s_bound_poly *bpoly)
{
    initialize_ncells_counter(setup);
    
    // Extract faces convhull setup
    int *faces_setup, Nf_setup;
    double **fnormals_setup;
    extract_faces_convhull_from_points(setup->points, setup->N_points, 
                                       &faces_setup, &fnormals_setup, &Nf_setup);
    
    s_vdiagram *vdiagram = malloc_vdiagram(setup);
    vdiagram->bpoly = bpoly;
    
    for (int ii=0; ii<setup->N_points; ii++) {
        vdiagram->vcells[ii] = extract_voronoi_cell(vdiagram, setup, faces_setup,
                                fnormals_setup, Nf_setup, ii);
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

void plot_add_vcell(FILE *pipe, const s_vcell *vcell, char *config)
{
    for (int ii=0; ii<vcell->Nf; ii++) {
            fprintf(pipe, "\"<echo \'");
            fprintf(pipe, "%f %f %f\\n", vcell->vertices[vcell->faces[ii*3 + 0]][0],
                                         vcell->vertices[vcell->faces[ii*3 + 0]][1],
                                         vcell->vertices[vcell->faces[ii*3 + 0]][2]);
            fprintf(pipe, "%f %f %f\\n", vcell->vertices[vcell->faces[ii*3 + 1]][0],
                                         vcell->vertices[vcell->faces[ii*3 + 1]][1],
                                         vcell->vertices[vcell->faces[ii*3 + 1]][2]);
            fprintf(pipe, "%f %f %f\'\"", vcell->vertices[vcell->faces[ii*3 + 2]][0],
                                          vcell->vertices[vcell->faces[ii*3 + 2]][1],
                                          vcell->vertices[vcell->faces[ii*3 + 2]][2]);
            fprintf(pipe, "%s, ", config);
        }
}


void plot_add_bpoly(FILE *pipe, const s_bound_poly *bpoly, char *config)
{
    for (int ii=0; ii<bpoly->Nf; ii++) {
        fprintf(pipe, "\"<echo \'");
        fprintf(pipe, "%f %f %f\\n", bpoly->points[bpoly->faces[ii*3]][0],
                                     bpoly->points[bpoly->faces[ii*3]][1], 
                                     bpoly->points[bpoly->faces[ii*3]][2]);
        fprintf(pipe, "%f %f %f\\n", bpoly->points[bpoly->faces[ii*3+1]][0],
                                     bpoly->points[bpoly->faces[ii*3+1]][1], 
                                     bpoly->points[bpoly->faces[ii*3+1]][2]);
        fprintf(pipe, "%f %f %f'\"", bpoly->points[bpoly->faces[ii*3+2]][0],
                                     bpoly->points[bpoly->faces[ii*3+2]][1], 
                                     bpoly->points[bpoly->faces[ii*3+2]][2]);
        fprintf(pipe, "%s, ", config);
    }
}


void plot_vcell(s_vdiagram *vdiag, s_vcell *vcell, char *f_name, double *ranges)
{
    FILE *pipe = popen("gnuplot -persistent 2>&1", "w");
    fprintf(pipe, "set terminal pngcairo enhanced font 'Arial,18' size 1080,1080 enhanced \n");
    fprintf(pipe, "set output '%s.png'\n", f_name);
    fprintf(pipe, "set pm3d depthorder\n");
    fprintf(pipe, "set view 100, 10, \n");
    fprintf(pipe, "set xyplane at 0\n");
    fprintf(pipe, "set xrange [%f:%f]\n", ranges[0], ranges[1]);
    fprintf(pipe, "set yrange [%f:%f]\n", ranges[2], ranges[3]);
    fprintf(pipe, "set zrange [%f:%f]\n", ranges[4], ranges[5]);
    fprintf(pipe, "set xlabel 'x'\n");
    fprintf(pipe, "set ylabel 'y'\n");
    fprintf(pipe, "set zlabel 'z'\n");
    fflush(pipe);

    fprintf(pipe, "set output '%s.png'\n", f_name);
    fprintf(pipe, "splot ");
    plot_add_vcell(pipe, vcell, "w polygons fs transparent solid 0.6 notitle");
    plot_add_bpoly(pipe, vdiag->bpoly, "w polygonf transparent solid 0.1 notitle");
    fprintf(pipe, "\n");
    pclose(pipe);
}


void plot_vdiagram(s_vdiagram *vdiagram, char *f_name, double *ranges, int max_files, double **aux_points, int *N_aux)
{
    char colors[][20] = { "#000090", "#000fff", "#0090ff", "#0fffee", 
        "#90ff70", "#ffee00", "#ff7000", "#ee0000", "#7f0000" };
    char buff[1024];

    FILE *pipe = popen("gnuplot -persistent 2>&1", "w");
    fprintf(pipe, "set terminal pngcairo enhanced font 'Arial,18' size 1080,1080 enhanced \n");
    fprintf(pipe, "set pm3d depth\n");
    fprintf(pipe, "set pm3d border lc 'black' lw 0.5\n");
    fprintf(pipe, "set view 100, 10, 1.5\n");
    fprintf(pipe, "unset border\n");
    fprintf(pipe, "unset xtics\n");
    fprintf(pipe, "unset ytics\n");
    fprintf(pipe, "unset ztics\n");
    fprintf(pipe, "set xrange [%f:%f]\n", ranges[0], ranges[1]);
    fprintf(pipe, "set yrange [%f:%f]\n", ranges[2], ranges[3]);
    fprintf(pipe, "set zrange [%f:%f]\n", ranges[4], ranges[5]);
    fflush(pipe);

    // PLOT POINTS
    if (aux_points) {
        fprintf(pipe, "set output '%s_aux_0.png'\n", f_name);
        fprintf(pipe, "splot ");

        for (int jj=0; jj<vdiagram->N_vcells; jj++) {
            fprintf(pipe, "\"<echo \'");
            fprintf(pipe, "%f %f %f\\n", vdiagram->seeds[jj][0], vdiagram->seeds[jj][1], vdiagram->seeds[jj][2]);
            fprintf(pipe, "'\" pt 7 lc rgb 'black' notitle, ");
        }

        fprintf(pipe, "\"<echo \'");
        for (int jj=0; jj<*N_aux; jj++) {
            fprintf(pipe, "%f %f %f\\n", aux_points[jj][0], aux_points[jj][1], aux_points[jj][2]);
        }
        fprintf(pipe, "'\" pt 3 lc rgb 'red' notitle, ");

        plot_add_bpoly(pipe, vdiagram->bpoly, "w polygons fs transparent solid 0.01 fc 'black' notitle");
        fprintf(pipe, "\n");

        fprintf(pipe, "set output '%s_aux_1.png'\n", f_name);
        fprintf(pipe, "set view 100, 90, 1.5\n");
        fprintf(pipe, "replot\n");

        fprintf(pipe, "set output '%s_aux_2.png'\n", f_name);
        fprintf(pipe, "set view 100, 180, 1.5\n");
        fprintf(pipe, "replot\n");

        fprintf(pipe, "set output '%s_aux_3.png'\n", f_name);
        fprintf(pipe, "set view 100, 270, 1.5\n");
        fprintf(pipe, "replot\n");
    }
    


    // PLOT WITH ALL CELLS
    fprintf(pipe, "set output '%s_v1.png'\n", f_name);
    fprintf(pipe, "splot ");
    for (int ii=0; ii<vdiagram->N_vcells; ii++) {
        s_vcell *vcell = vdiagram->vcells[ii];
        snprintf(buff, 1024, "w polygons fs transparent solid 0.2 fc rgb '%s' notitle", colors[ii%8]);
        plot_add_vcell(pipe, vcell, buff);  
    }       
    // plot_add_bpoly(pipe, vdiagram->bpoly, "w polygons fs transparent solid 0.01 fc 'black' notitle");

    if (aux_points) {
        fprintf(pipe, "\"<echo \'");
        for (int jj=0; jj<*N_aux; jj++) {
            fprintf(pipe, "%f %f %f\\n", aux_points[jj][0], aux_points[jj][1], aux_points[jj][2]);
        }
        fprintf(pipe, "'\" pt 3 lc rgb 'red' notitle, ");
    }

    fprintf(pipe, "\n");

    fprintf(pipe, "set output '%s_v2.png'\n", f_name);
    fprintf(pipe, "set view 100, 90, 1.5\n");
    fprintf(pipe, "replot\n");

    fprintf(pipe, "set output '%s_v3.png'\n", f_name);
    fprintf(pipe, "set view 100, 180, 1.5\n");
    fprintf(pipe, "replot\n");

    fprintf(pipe, "set output '%s_v4.png'\n", f_name);
    fprintf(pipe, "set view 100, 270, 1.5\n");
    fprintf(pipe, "replot\n");


    // A NEW PLOT FOR EACH CELL (UNTIL MAX_FILES)
    for (int ii=0; ii<vdiagram->N_vcells; ii++) {
        if (max_files != 0 && ii > max_files) break;

        s_vcell *vcell = vdiagram->vcells[ii];
        fprintf(pipe, "set output '%s_%d.png'\n", f_name, ii);
        fprintf(pipe, "splot ");

        snprintf(buff, 1024, "w polygons fs transparent solid 0.2 fc rgb '%s' notitle", colors[ii%8]);
        plot_add_vcell(pipe, vcell, buff);

        plot_add_bpoly(pipe, vdiagram->bpoly, "w polygons fs transparent solid 0.01 fc 'black' notitle");

        for (int jj=0; jj<vdiagram->N_vcells; jj++) {
            fprintf(pipe, "\"<echo \'");
            fprintf(pipe, "%f %f %f\\n", vdiagram->seeds[jj][0], vdiagram->seeds[jj][1], vdiagram->seeds[jj][2]);
            fprintf(pipe, "'\" pt 7 lc rgb 'black' notitle, ");
        }

        fprintf(pipe, "\"<echo \'");
        fprintf(pipe, "%f %f %f\\n", vdiagram->seeds[ii][0], vdiagram->seeds[ii][1], vdiagram->seeds[ii][2]);
        fprintf(pipe, "'\" pt 6 ps 2 lc rgb 'black' notitle, ");

        if (aux_points) {
            fprintf(pipe, "\"<echo \'");
            for (int jj=0; jj<*N_aux; jj++) {
                fprintf(pipe, "%f %f %f\\n", aux_points[jj][0], aux_points[jj][1], aux_points[jj][2]);
            }
            fprintf(pipe, "'\" pt 3 lc rgb 'red' notitle, ");
    }

        fprintf(pipe, "\n");
    }

    pclose(pipe);
}


void append_volumes_to_file(s_vdiagram *vdiagram, char *fname)
{
    FILE *f = fopen(fname, "a");
    assert(f && "Could not open file to write volumes to");

    for (int ii=0; ii<vdiagram->N_vcells; ii++) {
        fprintf(f, "%f, %f, %f, %f\n", vdiagram->seeds[vdiagram->vcells[ii]->seed_id][0],
                                       vdiagram->seeds[vdiagram->vcells[ii]->seed_id][1], 
                                       vdiagram->seeds[vdiagram->vcells[ii]->seed_id][2], 
                                       vdiagram->vcells[ii]->volume);   
    }
    


    fclose(f);
}


#endif
