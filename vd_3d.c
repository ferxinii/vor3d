
// #define CONVHULL_3D_ENABLE

#include "simplical_complex.c"
#include <float.h>
// #include "convhull_3d.h"

#define VCELL_MAX_VERTICES 1000

typedef struct bound_poly {
    int Np;
    double **points;
    int Nf;
    int **faces;
} s_bound_poly;


typedef struct vdiagram {
    int N_vcells;
    struct vcell *head;
    struct bound_poly *bpoly;
} s_vdiagram;


typedef struct vcell {
    int Nv;
    double **vertices;  // Nv x 4, where the last column indicates either the dual ncell (POSITIVE) 
                        // or the corresponding bound_poly index (NEGATIVE)
    //int Nf;
    // struct vface *head;  // linked list of faces
    // struct vface *tail;
} s_vcell;


// typedef struct vface {
//     int Nv;
//     int *vertex_id;
//     struct vface *next;
//     int bounded;
//     int id_ridge_1, id_ridge_2;  // These correspond to the ridges of the vt dual to the face
// } s_vface;


s_vcell *malloc_vcell(void)
{
    s_vcell *out = malloc(sizeof(s_vcell));
    out->Nv = 0;
    out->vertices = malloc_matrix(VCELL_MAX_VERTICES, 3);
    // out->head = NULL;
    // out->tail = NULL;
    return out;
}


// s_vface *malloc_vface(int N_vertices)
// {
//     s_vface *out = malloc(sizeof(s_vface));
//     out->Nv = N_vertices;
//     out->vertex_id = malloc(sizeof(int) * N_vertices);
//     out->next = NULL;
//     return out;
// }


// void print_vface(const s_vface *vface)
// {
//     printf("VFACE:   Nv = %d;    ", vface->Nv);
//     printf("%d", vface->vertex_id[0]);
//     for (int ii=1; ii<vface->Nv; ii++) {
//         printf(", %d", vface->vertex_id[ii]);
//     }
//     printf("\n");
// }


void print_vcell(const s_vcell *vcell)
{
    puts("VCELL");
    printf("Nv = %d\n", vcell->Nv);
    printf("%f, %f, %f\n", vcell->vertices[0][0], vcell->vertices[0][1], vcell->vertices[0][2]);
    for (int ii=1; ii<vcell->Nv; ii++) {
        printf("%f, %f, %f\n", vcell->vertices[ii][0], vcell->vertices[ii][1], vcell->vertices[ii][2]);
    }
    printf("\n");

    // s_vface *current = vcell->head;
    // while (current) {
    //     print_vface(current);
    //     current = current->next;
    // }
}


int add_vvertex_from_ncell(const s_setup *setup, const s_ncell *ncell, s_vcell *vcell)
{   // Returns the index of the vertex
    // Check if the vertex already exists
    for (int ii=0; vcell->Nv; ii++) {
        if ((int)vcell->vertices[ii][3] == ncell->count) return ii;
    } 

    static double **vertices_ncell = NULL;
    if (!vertices_ncell) vertices_ncell = malloc_matrix(4, 3);

    extract_vertices_ncell(setup, ncell, vertices_ncell);

    double CM[3];
    find_center_mass(vertices_ncell, 4, 3, CM);

    vcell->vertices[vcell->Nv][0] = CM[0];
    vcell->vertices[vcell->Nv][1] = CM[1];
    vcell->vertices[vcell->Nv][2] = CM[2];
    vcell->vertices[vcell->Nv][3] = ncell->count;
    vcell->Nv++;

    return vcell->Nv - 1;
}


int add_vvertex_from_coords(double x, double y, double z, s_vcell *vcell)
{
    vcell->vertices[vcell->Nv][0] = x;
    vcell->vertices[vcell->Nv][1] = y;
    vcell->vertices[vcell->Nv][2] = z;
    vcell->vertices[vcell->Nv][3] = -1;
    vcell->Nv++;

    return vcell->Nv - 1;
}



// void extract_vertex(const s_setup *setup, const s_ncell *ncell, double *out)
// {
//     static double **vertices_ncell = NULL;
//     if (!vertices_ncell) vertices_ncell = malloc_matrix(4, 3);
//     
//     extract_vertices_ncell(setup, ncell, vertices_ncell);
//
//     find_center_mass(vertices_ncell, 4, 3, out);
// }


// void extract_vertices_vface(const s_vcell *vcell, const s_vface *vface, double **out)
// {
//     for (int ii=0; ii<vface->Nv; ii++) {
//         out[ii][0] = vcell->vertices[vface->vertex_id[ii]][0];
//         out[ii][1] = vcell->vertices[vface->vertex_id[ii]][1];
//         out[ii][2] = vcell->vertices[vface->vertex_id[ii]][2];
//     }
// }


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
            puts("DEBUG: INTERSECTION!!");
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
    

    // DETERMINE EXTENDING DIRECTIONS
    static double **face_v = NULL, **ncell_v = NULL;
    if (!face_v) face_v = malloc_matrix(3, 3);
    if (!ncell_v) ncell_v = malloc_matrix(4, 3);
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

    extract_vertices_face(setup, ncell_B, &v1_B, 2, face_v);
    extract_vertices_ncell(setup, ncell_B, ncell_v);
    double n_B[3];
    subtract_3d(face_v[1], face_v[0], d1);
    subtract_3d(face_v[2], face_v[0], d2);
    cross_3d(d1, d2, n_B);
    double fc_B[3];
    find_center_mass(face_v, 3, 3, fc_B);
    find_center_mass(ncell_v, 4, 3, cc);
    subtract_3d(fc_B, cc, v);
    dir_aux = dot_3d(v, n_B);
    assert(dir_aux != 0 && "Vectors are perpendicular?");
    if (dir_aux < 0) {
        n_B[0] = -n_B[0];
        n_B[1] = -n_B[1];
        n_B[2] = -n_B[2];
    }


    // FIND INTERSECTION POINT WITH BOUNDING POLYHEDRON
    double intersection_A[3];
    find_intersection_with_bounding_poly(vdiagram->bpoly, fc_A, n_A, intersection_A);
    printf("INTERSECTION: %f, %f, %f\n", intersection_A[0], intersection_A[1], intersection_A[2]);

    double intersection_B[3];
    find_intersection_with_bounding_poly(vdiagram->bpoly, fc_B, n_B, intersection_B);
    printf("INTERSECTION: %f, %f, %f\n", intersection_B[0], intersection_B[1], intersection_B[2]);
    
    // ADD INTERSECTION TO CELL VERTICES
    add_vvertex_from_coords(intersection_A[0], intersection_A[1], intersection_A[2], vcell);
    add_vvertex_from_coords(intersection_B[0], intersection_B[1], intersection_B[2], vcell);
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
    if (Nv == -1) {
        extract_vface_BOUNDING(vdiagram, setup, ncell, v_localid_vertex, v_localid_2, vcell);
        printf("We must deal with this uncomplete ridge...");

        exit(1);
    } else {
        printf("Nv = %d\n", Nv);

        // s_vface *vface = malloc_vface(Nv);
        // vface->id_ridge_1 = ncell->vertex_id[v_localid_vertex];
        // vface->id_ridge_1 = ncell->vertex_id[v_localid_2];
        // vface->bounded = 1;
        //
        // vface->vertex_id[0] = ncell->count;


        int new_v_ridge_1, new_v_ridge_2;
        const s_ncell *current_ncell = ncell;
        for (int ii=1; ii<Nv; ii++) {
            s_ncell *next_ncell = next_ncell_ridge_cycle(setup, current_ncell, v_ridge_1, v_ridge_2, &new_v_ridge_1, &new_v_ridge_2);
            
            add_vvertex_from_ncell(setup, next_ncell, vcell);
            // vface->vertex_id[ii] = id;

            v_ridge_1 = new_v_ridge_1;
            v_ridge_2 = new_v_ridge_2;
            current_ncell = next_ncell;
        }
    }
}


s_vdiagram *initialize_vdiagram(const s_setup *setup)
{
    s_vdiagram *out = malloc(sizeof(s_vdiagram));
    out->N_vcells = setup->N_points;
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


// int face_exists(s_vcell *vcell, int id_ridge_1, int id_ridge_2) 
// {
//     s_vface *current = vcell->head;
//     while (current) {
//         if ( ((current->id_ridge_1 == id_ridge_1) && (current->id_ridge_2 == id_ridge_2)) ||
//              ((current->id_ridge_1 == id_ridge_2) && (current->id_ridge_2 == id_ridge_1)) ) {
//             return 1;
//         }
//         current = current->next;
//     }
//     return 0;
// }


void extract_vfaces_from_ncell_and_vertex(const s_vdiagram *vdiagram, const s_setup *setup, s_vcell *vcell, const s_ncell *ncell, int vertex_id)
{
    int v_localid_1 = id_where_equal_int(ncell->vertex_id, 4, vertex_id);
    for (int v_localid_2 = 0; v_localid_2 < 4; v_localid_2++) {
        // if (v_localid_2 != v_localid_1 && !face_exists(vcell, vertex_id, ncell->vertex_id[v_localid_2])) {
        extract_vface(vdiagram, setup, ncell, v_localid_1, v_localid_2, vcell);
            
            // if (!vcell->head) {
            //     vcell->head = vface;
            //     vcell->tail = vface;
            // } else {
            //     vcell->tail->next = vface;
            //     vcell->tail = vface;
            // }
        // }
    }
}


s_vcell *extract_voronoi_cell(const s_vdiagram *vdiagram, const s_setup *setup, int vertex_id)
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

    // Mark ncells incident to point (dim 0)
    initialize_ncells_mark(setup);
    initialize_ncells_counter(setup);
    mark_ncells_incident_face(setup, ncell, v_localid_COMP, 0);

    // int Nv = count_marked(setup);
    s_vcell *vcell = malloc_vcell();
    
    // Extract vertices of voronoi cell
    s_ncell *current = setup->head;
    // kk = 0;
    while (current) {
        if (current->mark == 1) {
            add_vvertex_from_ncell(setup, current, vcell);
            extract_vfaces_from_ncell_and_vertex(vdiagram, setup, vcell, current, vertex_id);
            // kk++;
        }
        current = current->next;
    }

    return vcell;
}


// TRIANGULATION OF THE FACES

// double g_QSORT_CM[3];
// double g_QSORT_ref[3];
// double g_QSORT_normal[3];
// const s_vdiagram *g_QSORT_vdiagram;
//
//
// void INIT_compare_angles(const s_vdiagram *vdiagram, const double *CM, const double *ref, const double *normal)
// {
//     g_QSORT_vdiagram = vdiagram;
//
//     g_QSORT_CM[0] = CM[0];
//     g_QSORT_CM[1] = CM[1];
//     g_QSORT_CM[2] = CM[2];
//
//     g_QSORT_ref[0] = ref[0];
//     g_QSORT_ref[1] = ref[1];
//     g_QSORT_ref[2] = ref[2];
//
//     g_QSORT_normal[0] = normal[0];
//     g_QSORT_normal[1] = normal[1];
//     g_QSORT_normal[2] = normal[2];
// }
//
//
// int FUN_compare_angles(const void *a, const void *b)
// {
//     const int *id1 = (const int *)a;
//     const int *id2 = (const int *)b;
//
//     double *v1 = g_QSORT_vdiagram->vertices[*id1];
//     double *v2 = g_QSORT_vdiagram->vertices[*id2];
//
//     double d1[3], d2[3];
//     subtract_3d(v1, g_QSORT_CM, d1);
//     subtract_3d(v2, g_QSORT_CM, d2);
//
//     double angle1 = compute_signed_angle(g_QSORT_ref, d1, g_QSORT_normal);
//     double angle2 = compute_signed_angle(g_QSORT_ref, d2, g_QSORT_normal);
//     
//     if (angle1 < angle2)
//         return -1;
//     else if (angle1 > angle2)
//         return 1;
//     else
//         return 0;
// }
//
//
// void order_vertices_face(const s_vdiagram *vdiagram, s_vface *vface)
// {   
//     static double **vertices_face = NULL;
//     static int prev_Nv = 0;
//     if (prev_Nv != vface->Nv) {
//         if (vertices_face) free_matrix(vertices_face, prev_Nv);
//         vertices_face = malloc_matrix(vface->Nv, 3);
//     }
//     
//     extract_vertices_vface(vdiagram, vface, vertices_face);
//     double CM[3];
//     find_center_mass(vertices_face, vface->Nv, 3, CM);
//     
//     double *v1 = vdiagram->vertices[vface->vertex_id[0]];
//     double *v2 = vdiagram->vertices[vface->vertex_id[1]];
//     double d1[3], d2[3];  // Difference vector wrt CM
//     subtract_3d(v1, CM, d1);
//     subtract_3d(v2, CM, d2);
//
//     double normal[3];
//     cross_3d(d1, d2, normal);
//     
//     INIT_compare_angles(vdiagram, CM, d1, normal);
//     qsort(vface->vertex_id, vface->Nv, sizeof(int), &FUN_compare_angles);
// }


// void triangulate_vface(const s_vdiagram *vdiagram, s_vface *vface, int **out)
// {   // OUT MUST BE PREALLOCATED AS int[vface->Nv - 2][3]
//     // order_vertices_face(vdiagram, vface);
//     (void)vdiagram;
//     for (int ii = 0; ii<vface->Nv-2; ii++) {
//         out[ii][0] = vface->vertex_id[0];
//         out[ii][1] = vface->vertex_id[ii+1];
//         out[ii][2] = vface->vertex_id[ii+2];
//     }
// }


// void triangulate_vcell(const s_vdiagram *vdiagram, s_vcell *vcell, int ***out, int *N_triangles)
// {
//     *N_triangles = 0;
//     s_vface *current = vcell->head;
//     while (current) {
//         *N_triangles += current->Nv - 2;
//         current = current->next;
//     }
//
//     *out = malloc_matrix_int(*N_triangles, 3);
//     current = vcell->head;
//     int kk = 0;
//     while (current) {
//         triangulate_vface(vdiagram, current, &((*out)[kk]));
//         kk += current->Nv - 2;
//         current = current->next;
//     }
// }


// void plot_vcell(const s_vdiagram *vdiagram, s_vcell *vcell, char *f_name, double *ranges)
// {
//     int **triangulation, N_triangles;
//     triangulate_vcell(vdiagram, vcell, &triangulation, &N_triangles);
//
//
//     FILE *pipe = popen("gnuplot -persistent 2>&1", "w");
//     fprintf(pipe, "set terminal pngcairo enhanced font 'Arial,18' size 1080,1080 enhanced \n");
//     fprintf(pipe, "set output '%s.png'\n", f_name);
//     fprintf(pipe, "set pm3d depthorder\n");
//     fprintf(pipe, "set view 100, 10, \n");
//     fprintf(pipe, "set xyplane at 0\n");
//     fprintf(pipe, "set xrange [%f:%f]\n", ranges[0], ranges[1]);
//     fprintf(pipe, "set yrange [%f:%f]\n", ranges[2], ranges[3]);
//     fprintf(pipe, "set zrange [%f:%f]\n", ranges[4], ranges[5]);
//     fprintf(pipe, "set xlabel 'x'\n");
//     fprintf(pipe, "set ylabel 'y'\n");
//     fprintf(pipe, "set zlabel 'z'\n");
//     fflush(pipe);
//
//     for (int ii=0; ii<N_triangles; ii++) {
//         fprintf(pipe, "set output '%s_%d.png'\n", f_name, ii);
//         fprintf(pipe, "splot ");
//         fprintf(pipe, "\"<echo \'");
//         fprintf(pipe, "%f %f %f\\n", vdiagram->vertices[triangulation[ii][0]][0], 
//                                      vdiagram->vertices[triangulation[ii][0]][1], 
//                                      vdiagram->vertices[triangulation[ii][0]][2]);
//         fprintf(pipe, "%f %f %f\\n", vdiagram->vertices[triangulation[ii][1]][0], 
//                                      vdiagram->vertices[triangulation[ii][1]][1], 
//                                      vdiagram->vertices[triangulation[ii][1]][2]);
//         fprintf(pipe, "%f %f %f'\"", vdiagram->vertices[triangulation[ii][2]][0], 
//                                      vdiagram->vertices[triangulation[ii][2]][1], 
//                                      vdiagram->vertices[triangulation[ii][2]][2]);
//         fprintf(pipe, "w polygons fs transparent solid 0.6 notitle");
//         fprintf(pipe, "\n");
//     }
//
//     fprintf(pipe, "set output '%s.png'\n", f_name);
//         fprintf(pipe, "splot ");
//     for (int ii=0; ii<N_triangles; ii++) {
//         fprintf(pipe, "\"<echo \'");
//         fprintf(pipe, "%f %f %f\\n", vdiagram->vertices[triangulation[ii][0]][0], 
//                                      vdiagram->vertices[triangulation[ii][0]][1], 
//                                      vdiagram->vertices[triangulation[ii][0]][2]);
//         fprintf(pipe, "%f %f %f\\n", vdiagram->vertices[triangulation[ii][1]][0], 
//                                      vdiagram->vertices[triangulation[ii][1]][1], 
//                                      vdiagram->vertices[triangulation[ii][1]][2]);
//         fprintf(pipe, "%f %f %f'\"", vdiagram->vertices[triangulation[ii][2]][0], 
//                                      vdiagram->vertices[triangulation[ii][2]][1], 
//                                      vdiagram->vertices[triangulation[ii][2]][2]);
//         fprintf(pipe, "w polygons fs transparent solid 0.6 notitle, ");
//     }
//     fprintf(pipe, "\n");
//     
//     pclose(pipe);
// }
