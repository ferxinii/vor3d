// TODO Improve malloc of vertices, or check if i have reached VCELL_MAX_VERTICES to increase size as needed

#define CONVHULL_3D_ENABLE

#include "simplical_complex.c"
#include <float.h>
#include "convhull_3d.h"

#define VCELL_MAX_VERTICES 1000

typedef struct bound_poly {
    int Np;
    double **points;
    int Nf;
    int **faces;
} s_bound_poly;


typedef struct vdiagram {
    int N_vcells;
    struct vcell *head, *tail;  // linked list of vcells
    const struct bound_poly *bpoly;
} s_vdiagram;


typedef struct vcell {
    double seed[3];
    struct vcell *next;     // linked list of vcells
    int Nv;
    double **vertices;  // Nv x 3
    int **origin_vertices;  // Nv x 4, LAST column indicates the dual ncell if POSITIVE,
                            // if -1, the rest of the columns indicate the delaunay indices
                            // of the face whose normal was extended
    int Nf;
    int *faces;
} s_vcell;


void add_noise_to_bp(s_bound_poly *bpoly)
{   // ADD SOME NOISE TO AVOID COLINEARITIES
    for (int ii=0; ii<bpoly->Np; ii++) {
        double aux = 2.0 * rand() / RAND_MAX - 1;  
        bpoly->points[ii][0] += 0.0001 * aux;
    }
}


s_vdiagram *malloc_vdiagram(void)
{
    s_vdiagram *out = malloc(sizeof(s_vdiagram));
    out->N_vcells = 0;
    out->bpoly = NULL;
    out->head = NULL;
    out->tail = NULL;
    return out;
}


void initialize_ncells_counter(const s_setup *setup)
{
    s_ncell *current = setup->head;
    int ii = 0;
    while (current) {
        current->count = ii;
        current = current->next;
        ii++;
    }
}


void print_vcell(const s_vcell *vcell)
{
    puts("VCELL");
    printf("Seed: %f, %f, %f\n", vcell->seed[0], vcell->seed[1], vcell->seed[2]);
    printf("Nv = %d\n", vcell->Nv);
    printf("%f, %f, %f; %d\n", vcell->vertices[0][0], vcell->vertices[0][1], vcell->vertices[0][2], 
                               vcell->origin_vertices[0][3]);
    for (int ii=1; ii<vcell->Nv; ii++) {
        printf("%f, %f, %f; %d\n", vcell->vertices[ii][0], vcell->vertices[ii][1], vcell->vertices[ii][2], 
                               vcell->origin_vertices[ii][3]);
    }

    printf("Nf = %d\n", vcell->Nf);
    printf("%d, %d, %d\n", vcell->faces[0], vcell->faces[1], vcell->faces[2]);
    for (int ii=1; ii<vcell->Nv; ii++) {
        printf("%d, %d, %d\n", vcell->faces[ii*3 + 0], vcell->faces[ii*3 + 1], vcell->faces[ii*3 + 2]);
    }
    printf("\n");
}


void print_vdiagram(const s_vdiagram *vdiagram)
{
    s_vcell *vcell = vdiagram->head;
    print_vcell(vcell);
    for (int ii=1; ii<vdiagram->N_vcells; ii++) {
        vcell = vcell->next;
        print_vcell(vcell);
    }
}


s_vcell *malloc_vcell(void)
{
    s_vcell *out = malloc(sizeof(s_vcell));
    out->Nv = 0;
    out->vertices = malloc_matrix(VCELL_MAX_VERTICES, 3);
    out->origin_vertices = malloc_matrix_int(VCELL_MAX_VERTICES, 4);
    out->next = NULL;
    return out;
}


int add_vvertex_from_ncell(const s_setup *setup, const s_ncell *ncell, s_vcell *vcell)
{   // Returns the index of the vertex
    // First check if the vertex already exists
    for (int ii=0; ii<vcell->Nv; ii++) {
        if (vcell->origin_vertices[ii][3] == ncell->count) {
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


int add_vvertex_from_coords(double x, double y, double z, int *dt_face_vid, s_vcell *vcell)
{   // This function assumes that the vertex does not already exist!
    vcell->vertices[vcell->Nv][0] = x;
    vcell->vertices[vcell->Nv][1] = y;
    vcell->vertices[vcell->Nv][2] = z;

    vcell->origin_vertices[vcell->Nv][3] = -1;
    vcell->origin_vertices[vcell->Nv][0] = dt_face_vid[0];
    vcell->origin_vertices[vcell->Nv][1] = dt_face_vid[1];
    vcell->origin_vertices[vcell->Nv][2] = dt_face_vid[2];

    vcell->Nv++;

    return vcell->Nv - 1;
}


void extract_vertices_face_bpoly(const s_bound_poly *bpoly, int *face, double **out)
{
    out[0][0] = bpoly->points[face[0]][0];
    out[0][1] = bpoly->points[face[0]][1];
    out[0][2] = bpoly->points[face[0]][2];

    out[1][0] = bpoly->points[face[1]][0];
    out[1][1] = bpoly->points[face[1]][1];
    out[1][2] = bpoly->points[face[1]][2];

    out[2][0] = bpoly->points[face[2]][0];
    out[2][1] = bpoly->points[face[2]][1];
    out[2][2] = bpoly->points[face[2]][2];
}


void find_intersection_with_bounding_poly(const s_bound_poly *bpoly, const double *origin, const double *dir, double *intersection)
{
    static double **face_v = NULL;
    if (!face_v) face_v = malloc_matrix(3, 3);
    double dmin = DBL_MAX;

    for (int ii=0; ii<bpoly->Nf; ii++) {
        extract_vertices_face_bpoly(bpoly, bpoly->faces[ii], face_v);
        double intersection_aux[3];
        int indicator = ray_triangle_intersection_3d(face_v, origin, dir, intersection_aux);
        if (indicator == 1) {
            double d = norm_difference(intersection_aux, origin, 3);
            if (d < dmin) {
                dmin = d;
                intersection[0] = intersection_aux[0];
                intersection[1] = intersection_aux[1];
                intersection[2] = intersection_aux[2];
            }
        }
    }
    assert(dmin != DBL_MAX && "Did not find any intersection with bounding poly.");
}


int face_extension_already_exists(const s_vcell *vcell, int *vertex_id)
{
    for (int ii=0; ii<vcell->Nv; ii++) {
        if (vcell->origin_vertices[ii][3] == -1 &&
        inarray(vcell->origin_vertices[ii], 3, vertex_id[0]) &&
        inarray(vcell->origin_vertices[ii], 3, vertex_id[1]) &&
        inarray(vcell->origin_vertices[ii], 3, vertex_id[2])) {
            return 1;
        }
    }
    return 0;
}


void extract_vface_BOUNDING(const s_vdiagram *vdiagram, const s_setup *setup, const s_ncell *ncell, int v_localid_vertex, int v_localid_2, s_vcell *vcell)
{   
    // FIRST, FIND BOTH ENDS OF THE UNCLOSED CYCLE
    const s_ncell *current = ncell;
    int v1_A = v_localid_vertex;
    int v2_A = v_localid_2;
    while (current->opposite[v1_A] != NULL) {
        int new_v1_A, new_v2_A;
        s_ncell *next = next_ncell_ridge_cycle(setup, current, v1_A, v2_A, &new_v1_A, &new_v2_A);
        add_vvertex_from_ncell(setup, current, vcell);

        current = next;
        v1_A = new_v1_A;
        v2_A = new_v2_A;
    }
    const s_ncell *ncell_A = current;

    static double **face_v = NULL;
    if (!face_v) face_v = malloc_matrix(3, 3);
    static double **ncell_v = NULL;
    if (!ncell_v) ncell_v = malloc_matrix(4, 3);

    if (!face_extension_already_exists(vcell, ncell_A->vertex_id)) {
        // DETERMINE EXTENDING DIRECTIONS, AND ADD CROSSING POINT
        extract_vertices_face(setup, ncell_A, &v1_A, 2, face_v);
        extract_vertices_ncell(setup, ncell_A, ncell_v);
        double d1[3], d2[3], n_A[3];
        subtract_3d(face_v[1], face_v[0], d1);
        subtract_3d(face_v[2], face_v[0], d2);
        cross_3d(d1, d2, n_A);
        double fc_A[3], cc[3], v[3];
        find_center_mass(face_v, 3, 3, fc_A);
        find_center_mass(ncell_v, 4, 3, cc);
        subtract_3d(fc_A, cc, v);
        double dir_aux = dot_3d(v, n_A);
        assert(dir_aux != 0 && "Vectors are perpendicular?");
        if (dir_aux < 0) {
            n_A[0] = -n_A[0];
            n_A[1] = -n_A[1];
            n_A[2] = -n_A[2];
        } 
    
        double intersection_A[3];
        find_intersection_with_bounding_poly(vdiagram->bpoly, fc_A, n_A, intersection_A);
        printf("INTERSECTION: %f, %f, %f\n", intersection_A[0], intersection_A[1], intersection_A[2]);
        add_vvertex_from_coords(intersection_A[0], intersection_A[1], intersection_A[2], 
                                ncell_A->vertex_id, vcell);
    }
    
    current = ncell;
    int v1_B = v_localid_2;
    int v2_B = v_localid_vertex;
    while (current->opposite[v1_B] != NULL) {
        int new_v1_B, new_v2_B;
        s_ncell *next = next_ncell_ridge_cycle(setup, current, v1_B, v2_B, &new_v1_B, &new_v2_B);
        add_vvertex_from_ncell(setup, current, vcell);

        current = next;
        v1_B = new_v1_B;
        v2_B = new_v2_B;
    }
    const s_ncell *ncell_B = current;
    
    if (!face_extension_already_exists(vcell, ncell_B->vertex_id)) {
        extract_vertices_face(setup, ncell_B, &v1_B, 2, face_v);
        extract_vertices_ncell(setup, ncell_B, ncell_v);
        double d1[3], d2[3], n_B[3];
        subtract_3d(face_v[1], face_v[0], d1);
        subtract_3d(face_v[2], face_v[0], d2);
        cross_3d(d1, d2, n_B);
        double fc_B[3], cc[3], v[3];
        find_center_mass(face_v, 3, 3, fc_B);
        find_center_mass(ncell_v, 4, 3, cc);
        subtract_3d(fc_B, cc, v);
        double dir_aux = dot_3d(v, n_B);
        assert(dir_aux != 0 && "Vectors are perpendicular?");
        if (dir_aux < 0) {
            n_B[0] = -n_B[0];
            n_B[1] = -n_B[1];
            n_B[2] = -n_B[2];
        }

        double intersection_B[3];
        find_intersection_with_bounding_poly(vdiagram->bpoly, fc_B, n_B, intersection_B);
        printf("INTERSECTION: %f, %f, %f\n", intersection_B[0], intersection_B[1], intersection_B[2]);
        
        add_vvertex_from_coords(intersection_B[0], intersection_B[1], intersection_B[2], 
                                ncell_B->vertex_id, vcell);
    }
}


void extract_vface(const s_vdiagram *vdiagram, const s_setup *setup, const s_ncell *ncell, int v_localid_vertex, int v_localid_2, s_vcell *vcell)
{
    // The ridge actually corresponds to the complementary indices
    int v_ridge_1 = 0;
    while (v_ridge_1 == v_localid_vertex || v_ridge_1 == v_localid_2) {
        v_ridge_1++;
    }
    int v_ridge_2 = 0;
    while (v_ridge_2 == v_localid_vertex || v_ridge_2 == v_localid_2 || v_ridge_2 == v_ridge_1) {
        v_ridge_2++;
    }

    int Nv = count_cycle_ridge(setup, ncell, v_ridge_1, v_ridge_2);
    printf("Nv = %d\n", Nv);
    if (Nv == -1) {
        extract_vface_BOUNDING(vdiagram, setup, ncell, v_ridge_1, v_ridge_2, vcell);
        printf("EXTRACTED VFACE BOUNDING, %p\n", (void*) ncell);
    } else {
        printf("Nv = %d\n", Nv);

        int new_v_ridge_1, new_v_ridge_2;
        const s_ncell *current_ncell = ncell;
        for (int ii=1; ii<Nv; ii++) {
            printf("ii=%d\n", ii);
            s_ncell *next_ncell = next_ncell_ridge_cycle(setup, current_ncell, v_ridge_1, v_ridge_2, &new_v_ridge_1, &new_v_ridge_2);
            
            add_vvertex_from_ncell(setup, next_ncell, vcell);

            v_ridge_1 = new_v_ridge_1;
            v_ridge_2 = new_v_ridge_2;
            current_ncell = next_ncell;
        }
    }
}



void extract_vfaces_from_ncell_and_vertex(const s_vdiagram *vdiagram, const s_setup *setup, s_vcell *vcell, const s_ncell *ncell, int vertex_id)
{   // The idea is to cycle around all ridges sharing vertex_id
    int vertex_localid = id_where_equal_int(ncell->vertex_id, 4, vertex_id);
    for (int v_localid_2 = 0; v_localid_2 < 4; v_localid_2++) {
        if (v_localid_2 != vertex_localid) {
            extract_vface(vdiagram, setup, ncell, vertex_localid, v_localid_2, vcell);
        }
    }
}


ch_vertex *convert_points_to_chvertex(double **points, int Np)
{
    ch_vertex *out = malloc(sizeof(ch_vertex) * Np);
    for (int ii=0; ii<Np; ii++) {
        out[ii].x = points[ii][0];
        out[ii].y = points[ii][1];
        out[ii].z = points[ii][2];
    }
    return out;
}


s_vcell *extract_voronoi_cell(const s_vdiagram *vdiagram, const s_setup *setup, int vertex_id)
{   
    s_vcell *vcell = malloc_vcell();

    // Find an ncell with this vertex
    s_ncell *ncell = setup->head;
    while (ncell) {
        if (inarray(ncell->vertex_id, 4, vertex_id)) {
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

            // Mark ncells incident to point (dim 0)
            initialize_ncells_mark(setup);
            mark_ncells_incident_face(setup, ncell, v_localid_COMP, 0);

            vcell->seed[0] = setup->points[vertex_id][0];
            vcell->seed[1] = setup->points[vertex_id][1];
            vcell->seed[2] = setup->points[vertex_id][2];
            
            // Extract vertices of voronoi cell
            s_ncell *current = setup->head;
            while (current) {
                if (current->mark == 1) {
                    add_vvertex_from_ncell(setup, current, vcell);
                    extract_vfaces_from_ncell_and_vertex(vdiagram, setup, vcell, current, vertex_id);
                }
                current = current->next;
            }

            int *faces, N_faces;
            ch_vertex *vertices_aux = convert_points_to_chvertex(vcell->vertices, vcell->Nv);
            convhull_3d_build(vertices_aux, vcell->Nv, &faces, &N_faces);
            vcell->faces = faces;
            vcell->Nf = N_faces;

        }
        ncell = ncell->next;
    }

    return vcell;
}


s_vdiagram *voronoi_from_delaunay_3d(const s_setup *setup, s_bound_poly *bpoly)
{   
    initialize_ncells_counter(setup);

    s_vdiagram *vdiagram = malloc_vdiagram();
    add_noise_to_bp(bpoly);
    vdiagram->bpoly = bpoly;
    
    for (int ii=0; ii<setup->N_points; ii++) {
        s_vcell *vcell = extract_voronoi_cell(vdiagram, setup, ii);
        if (!vdiagram->head) {
            vdiagram->head = vcell;
            vdiagram->tail = vcell;
        } else {
            vdiagram->tail->next = vcell;
            vdiagram->tail = vcell;
        }
        vdiagram->N_vcells++;
    }

    return vdiagram;
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

    for (int ii=0; ii<vcell->Nf; ii++) {
        fprintf(pipe, "set output '%s_%d.png'\n", f_name, ii);
        fprintf(pipe, "splot ");
        fprintf(pipe, "\"<echo \'");
        fprintf(pipe, "%f %f %f\\n", vcell->vertices[vcell->faces[ii*3 + 0]][0],
                                     vcell->vertices[vcell->faces[ii*3 + 0]][1],
                                     vcell->vertices[vcell->faces[ii*3 + 0]][2]);
        fprintf(pipe, "%f %f %f\\n", vcell->vertices[vcell->faces[ii*3 + 1]][0],
                                     vcell->vertices[vcell->faces[ii*3 + 1]][1],
                                     vcell->vertices[vcell->faces[ii*3 + 1]][2]);
        fprintf(pipe, "%f %f %f'\"", vcell->vertices[vcell->faces[ii*3 + 2]][0],
                                     vcell->vertices[vcell->faces[ii*3 + 2]][1],
                                     vcell->vertices[vcell->faces[ii*3 + 2]][2]);
        fprintf(pipe, "w polygons fs transparent solid 0.6 notitle, ");

        for (int ii=0; ii<vdiag->bpoly->Nf; ii++) {
            fprintf(pipe, "\"<echo \'");
            fprintf(pipe, "%f %f %f\\n", vdiag->bpoly->points[vdiag->bpoly->faces[ii][0]][0],
                                         vdiag->bpoly->points[vdiag->bpoly->faces[ii][0]][1], 
                                         vdiag->bpoly->points[vdiag->bpoly->faces[ii][0]][2]);
            fprintf(pipe, "%f %f %f\\n", vdiag->bpoly->points[vdiag->bpoly->faces[ii][1]][0],
                                         vdiag->bpoly->points[vdiag->bpoly->faces[ii][1]][1], 
                                         vdiag->bpoly->points[vdiag->bpoly->faces[ii][1]][2]);
            fprintf(pipe, "%f %f %f'\"", vdiag->bpoly->points[vdiag->bpoly->faces[ii][2]][0],
                                          vdiag->bpoly->points[vdiag->bpoly->faces[ii][2]][1], 
                                          vdiag->bpoly->points[vdiag->bpoly->faces[ii][2]][2]);
            fprintf(pipe, "w polygons fs transparent solid 0.2 notitle, ");
        }
        fprintf(pipe, "\"<echo \'");
        fprintf(pipe, "%f %f %f\'\"", vcell->seed[0], vcell->seed[1], vcell->seed[2]);
        fprintf(pipe, "pt 7 lc rgb 'black' notitle");
        fprintf(pipe, "\n");
    }
    
    fprintf(pipe, "set output '%s.png'\n", f_name);
    fprintf(pipe, "splot ");
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
    fprintf(pipe, "w polygons fs transparent solid 0.6 notitle, ");
    }
    fprintf(pipe, "\"<echo \'");
    fprintf(pipe, "%f %f %f'\"", vcell->seed[0], vcell->seed[1], vcell->seed[2]);
    fprintf(pipe, "pt 7 lc rgb 'black' notitle, ");
    for (int ii=0; ii<vdiag->bpoly->Nf; ii++) {
            fprintf(pipe, "\"<echo \'");
            fprintf(pipe, "%f %f %f\\n", vdiag->bpoly->points[vdiag->bpoly->faces[ii][0]][0],
                                         vdiag->bpoly->points[vdiag->bpoly->faces[ii][0]][1], 
                                         vdiag->bpoly->points[vdiag->bpoly->faces[ii][0]][2]);
            fprintf(pipe, "%f %f %f\\n", vdiag->bpoly->points[vdiag->bpoly->faces[ii][1]][0],
                                         vdiag->bpoly->points[vdiag->bpoly->faces[ii][1]][1], 
                                         vdiag->bpoly->points[vdiag->bpoly->faces[ii][1]][2]);
            fprintf(pipe, "%f %f %f'\"", vdiag->bpoly->points[vdiag->bpoly->faces[ii][2]][0],
                                          vdiag->bpoly->points[vdiag->bpoly->faces[ii][2]][1], 
                                          vdiag->bpoly->points[vdiag->bpoly->faces[ii][2]][2]);
            fprintf(pipe, "w polygons fs transparent solid 0.1 notitle, ");
        }
    fprintf(pipe, "\n");
    
    pclose(pipe);
}


void plot_vdiagram(s_vdiagram *vdiagram, char *f_name, double *ranges)
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


    fprintf(pipe, "splot ");
    s_vcell *vcell = vdiagram->head;
    for (int jj=0; jj<vdiagram->N_vcells; jj++) {
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
        fprintf(pipe, "w polygons fs transparent solid 0.2 notitle, ");
        }
        fprintf(pipe, "\"<echo \'");
        fprintf(pipe, "%f %f %f'\"", vcell->seed[0], vcell->seed[1], vcell->seed[2]);
        fprintf(pipe, "pt 7 lc rgb 'black' notitle, ");

        vcell = vcell->next;
    }       
    fprintf(pipe, "\n");

    vcell = vdiagram->head;
    for (int jj=0; jj<vdiagram->N_vcells; jj++) {
        fprintf(pipe, "set output '%s_%d.png'\n", f_name, jj);
        for (int ii=0; ii<vcell->Nf; ii++) {
            fprintf(pipe, "splot ");
            fprintf(pipe, "\"<echo \'");
            fprintf(pipe, "%f %f %f\\n", vcell->vertices[vcell->faces[ii*3 + 0]][0],
                                         vcell->vertices[vcell->faces[ii*3 + 0]][1],
                                         vcell->vertices[vcell->faces[ii*3 + 0]][2]);
            fprintf(pipe, "%f %f %f\\n", vcell->vertices[vcell->faces[ii*3 + 1]][0],
                                         vcell->vertices[vcell->faces[ii*3 + 1]][1],
                                         vcell->vertices[vcell->faces[ii*3 + 1]][2]);
            fprintf(pipe, "%f %f %f'\"", vcell->vertices[vcell->faces[ii*3 + 2]][0],
                                         vcell->vertices[vcell->faces[ii*3 + 2]][1],
                                         vcell->vertices[vcell->faces[ii*3 + 2]][2]);
            fprintf(pipe, "w polygons fs transparent solid 0.6 notitle, ");

            fprintf(pipe, "\"<echo \'");
            fprintf(pipe, "%f %f %f\'\"", vcell->seed[0], vcell->seed[1], vcell->seed[2]);
            fprintf(pipe, "pt 7 lc rgb 'black' notitle");
            fprintf(pipe, "\n");
        }
        vcell = vcell->next;
    }

    pclose(pipe);
}
