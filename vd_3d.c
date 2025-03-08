// [ ] TODO Improve malloc of vertices, or check if i have reached VCELL_MAX_VERTICES to increase size as needed
// [ ] IMPROVE RAY INTERSECTION WITH BOUNDING BOX...
// [ ] Plot each cell a different color



#include "simplical_complex.c"
#include <float.h>
#define CONVHULL_3D_ENABLE
#include "convhull_3d.h"

#define VCELL_MAX_VERTICES 1000

typedef struct bound_poly {
    int Np;
    double **points;
    int Nf;
    int *faces;  // Its flat! Nf x 3
    double **fnormals;
    double dmax;  // Max distance between two pairs of points
} s_bound_poly;


typedef struct vdiagram {
    int N_vcells;
    // struct vcell *head, *tail;  // linked list of vcells
    struct vcell **vcells;  // N_vcells x 1
    const struct bound_poly *bpoly;
    double **seeds;
} s_vdiagram;


typedef struct vcell {
    int seed_id;
    // struct vcell *next;     // linked list of vcells
    int Nv;
    double **vertices;  // Nv x 3
    int **origin_vertices;  // Nv x 4, LAST column indicates the dual ncell if POSITIVE,
                            // if -1: real extension. The rest of the columns indicate the delaunay indices
                            //        of the face whose normal was extended, 
                            // if -2: ARTIFICIAL extension
                            // if -3: it is from the bounding polyhedron, and the first index is the point id
                            // if -4: it comes from a segment, and the first to columns correspond to the polygon vertex ids 
                            //        where the segment comes from, and the third column to the face id
    int Nf;
    int *faces;
    double **fnormals;
} s_vcell;


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


void add_noise_to_bp(s_bound_poly *bpoly)
{   // ADD SOME NOISE TO AVOID COLINEARITIES... Necessary??
    for (int ii=0; ii<bpoly->Np; ii++) {
        for (int jj=0; jj<3; jj++) {
            double aux = 2.0 * rand() / RAND_MAX - 1;  
            bpoly->points[ii][jj] += 0.001 * aux;
        }
    }
}


void extract_dmax_bp(s_bound_poly *bpoly)
{
    double dmax = 0; 
    for (int ii=0; ii<bpoly->Np-1; ii++) {
        for (int jj=ii+1; jj<bpoly->Np; jj++) {
            double d = norm_difference(bpoly->points[ii], bpoly->points[jj], 3);
            if (d > dmax) 
                dmax = d;
        }
    }
    bpoly->dmax = dmax;
}


double **extract_normals_from_ch(ch_vertex *vertices, int Nv, int *faces, int Nf, double *ch_CM)
{   // Normalized normals
    static double **vertices_face = NULL;
    if (!vertices_face) vertices_face = malloc_matrix(3, 3);

    double **out = malloc_matrix(Nf, 3);
    for (int ii=0; ii<Nf; ii++) {
        ch_vertex v0 = vertices[faces[ii*3+0]];
        ch_vertex v1 = vertices[faces[ii*3+1]];
        ch_vertex v2 = vertices[faces[ii*3+2]];

        double d1[3], d2[3];
        d1[0] = v1.x - v0.x;
        d1[1] = v1.y - v0.y;
        d1[2] = v1.z - v0.z;
        d2[0] = v2.x - v0.x;
        d2[1] = v2.y - v0.y;
        d2[2] = v2.z - v0.z;

        double n[3];
        cross_3d(d1, d2, n);
        normalize_inplace(n, 3);
        
        vertices_face[0][0] = v0.x;     vertices_face[0][1] = v0.y;     vertices_face[0][2] = v0.z;
        vertices_face[1][0] = v1.x;     vertices_face[1][1] = v1.y;     vertices_face[1][2] = v1.z;
        vertices_face[2][0] = v2.x;     vertices_face[2][1] = v2.y;     vertices_face[2][2] = v2.z;
        double fc[3], dir[3];
        find_center_mass(vertices_face, 3, 3, fc);
        subtract_3d(fc, ch_CM, dir);
        double dot = dot_3d(dir, n);
        if (dot < 0) {
            n[0] = -n[0];
            n[1] = -n[1];
            n[2] = -n[2];
            puts("DEBUG: Normal flipped! Is this normal?");
        }

        out[ii][0] = n[0];
        out[ii][1] = n[1];
        out[ii][2] = n[2];
    }
    return out;
}


void extract_convhull_bp(s_bound_poly *bpoly)
{   
    ch_vertex *pch = convert_points_to_chvertex(bpoly->points, bpoly->Np);

    convhull_3d_build(pch, bpoly->Np, &bpoly->faces, &bpoly->Nf);
    
    double CM[3];
    find_center_mass(bpoly->points, bpoly->Np, 3, CM);
    double **normals = extract_normals_from_ch(pch, bpoly->Np, bpoly->faces, bpoly->Nf, CM);
    bpoly->fnormals = normals;
    free(pch);
}


int is_inside_convhull(double *query, double **pch, int Np, int *faces, double **fnormals, int Nf)
{   
    double *pf[3], fc[3];
    pf[0] = pch[faces[0]];
    pf[1] = pch[faces[1]];
    pf[2] = pch[faces[2]];
    find_center_mass(pf, 3, 3, fc);
    
    double d[3];
    subtract_3d(fc, query, d);
    double dot = dot_3d(d, fnormals[0]);
    assert(dot != 0 && "parallel normal...");

    int prev_sign = dot > 0 ? 1 : -1;

    for (int ii=1; ii<Nf; ii++) {
        pf[0] = pch[faces[ii*3 + 0]];
        pf[1] = pch[faces[ii*3 + 1]];
        pf[2] = pch[faces[ii*3 + 2]];
        find_center_mass(pf, 3, 3, fc);
        
        subtract_3d(fc, query, d);
        dot = dot_3d(d, fnormals[ii]);
        assert(dot != 0 && "parallel normal...");
        int sign = dot > 0 ? 1 : -1;
        if (prev_sign * sign < 0) {// Check change in sign
            return 0;
        }

        prev_sign = sign;
    }
    return 1;
}


void segment_convex_hull_intersection(const double *p0, const double *p1, double **pch, int Np, const int *faces, 
                                     double **fnormals, int Nf, int *ind1, double *i1, int *ind2, double *i2)
{       // ind1 indicates if i1 is a valid intersection, simil 2
    double tmin = 0, tmax = 1;

    double p1_p0[3];
    subtract_3d(p1, p0, p1_p0);

    for (int ii=0; ii<Nf; ii++) {
        double *n = fnormals[ii];
        double den = dot_3d(n, p1_p0);
        assert(fabs(den) > 1e-6 && "Perpendicular normal with segment.");

        double *q = pch[faces[ii*3]];
        double q_p0[3];     
        subtract_3d(q, p0, q_p0);
        double t = dot_3d(n, q_p0) / den;

        if (den < 0) {
            tmin = fmax(tmin, t);
        } else if (den >0) {
            tmax = fmin(tmax, t);
        }
    }
    
    *ind1 = 0;
    *ind2 = 0;
    int intersection = 0;
    double TOL = 1e-6;
    if (tmin > TOL && tmin < tmax) {
        *ind1 = 1;
        i1[0] = p0[0] + tmin * p1_p0[0];  
        i1[1] = p0[1] + tmin * p1_p0[1];  
        i1[2] = p0[2] + tmin * p1_p0[2];  
        intersection = 1;
    }
    if (tmax < 1-TOL && tmin < tmax) {
        *ind2 = 1;
        i2[0] = p0[0] + tmax * p1_p0[0];  
        i2[1] = p0[1] + tmax * p1_p0[1];  
        i2[2] = p0[2] + tmax * p1_p0[2];  
        intersection = 1;
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

    printf("Nf = %d\n", vcell->Nf);
    printf("%d, %d, %d\n", vcell->faces[0], vcell->faces[1], vcell->faces[2]);
    for (int ii=1; ii<vcell->Nf; ii++) {
        printf("%d, %d, %d\n", vcell->faces[ii*3 + 0], vcell->faces[ii*3 + 1], vcell->faces[ii*3 + 2]);
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
    out->vertices = malloc_matrix(VCELL_MAX_VERTICES, 3);
    out->origin_vertices = malloc_matrix_int(VCELL_MAX_VERTICES, 4);
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


int add_vvertex_from_coords(double *coords, int *dt_face_vid, s_vcell *vcell)
{   // This function assumes that the vertex does not already exist!
    vcell->vertices[vcell->Nv][0] = coords[0];
    vcell->vertices[vcell->Nv][1] = coords[1];
    vcell->vertices[vcell->Nv][2] = coords[2];

    vcell->origin_vertices[vcell->Nv][3] = -1;
    vcell->origin_vertices[vcell->Nv][0] = dt_face_vid[0];
    vcell->origin_vertices[vcell->Nv][1] = dt_face_vid[1];
    vcell->origin_vertices[vcell->Nv][2] = dt_face_vid[2];

    vcell->Nv++;

    return vcell->Nv - 1;
}


int add_vvertex_as_extension(double *coords, int *dt_face_vid, s_vcell *vcell)
{   // This function assumes that the vertex does not already exist!
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
    vcell->vertices[vcell->Nv][0] = bp->points[id][0];
    vcell->vertices[vcell->Nv][1] = bp->points[id][1];
    vcell->vertices[vcell->Nv][2] = bp->points[id][2];

    vcell->origin_vertices[vcell->Nv][3] = -3;
    vcell->origin_vertices[vcell->Nv][0] = id;

    vcell->Nv++;

    return vcell->Nv - 1;
}


int add_vvertex_from_segment(double *i1, int bp_vid_1, int bp_vid_2, int face_id, s_vcell *vcell)
{
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


void find_intersection_with_bounding_poly(const s_bound_poly *bpoly, const double *origin, double *dir, double *intersection, double *extension)
{   // in extension, we extend the intersection point in the ray direction a distance of 2*bpoly->dmax
    normalize_inplace(dir, 3);

    static double **face_v = NULL;
    if (!face_v) face_v = malloc_matrix(3, 3);
    double dmin = DBL_MAX;
    int intersected = 0;

    for (int ii=0; ii<bpoly->Nf; ii++) {
        extract_vertices_face_bpoly(bpoly, &bpoly->faces[ii*3], face_v);
        double intersection_aux[3];
        int indicator = ray_triangle_intersection_3d(face_v, origin, dir, intersection_aux);
        if (indicator == 1) {
            double d = norm_difference(intersection_aux, origin, 3);
            if (d < dmin) {
                intersected = 1;
                dmin = d;
                intersection[0] = intersection_aux[0];
                intersection[1] = intersection_aux[1];
                intersection[2] = intersection_aux[2];
            }
        }
    }
    
    assert(bpoly->dmax > 0.1 && "Is bpoly->dmax initialized?");
    double s = 2;
    extension[0] = intersection[0] + s * bpoly->dmax * dir[0];
    extension[1] = intersection[1] + s * bpoly->dmax * dir[1];
    extension[2] = intersection[2] + s * bpoly->dmax * dir[2];

    assert(intersected != 0 && "Did not find any intersection with bounding poly.");
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


void extract_face_center_and_normal(const s_setup *setup, const s_ncell *ncell, int face_localid, double *fc, double *n)
{
    static double **face_v = NULL;
    if (!face_v) face_v = malloc_matrix(3, 3);
    static double **ncell_v = NULL;
    if (!ncell_v) ncell_v = malloc_matrix(4, 3);

    extract_vertices_face(setup, ncell, &face_localid, 2, face_v);
    extract_vertices_ncell(setup, ncell, ncell_v);

    double d1[3], d2[3];
    subtract_3d(face_v[1], face_v[0], d1);
    subtract_3d(face_v[2], face_v[0], d2);
    cross_3d(d1, d2, n);

    double cc[3], v[3];
    find_center_mass(face_v, 3, 3, fc);
    find_center_mass(ncell_v, 4, 3, cc);
    subtract_3d(fc, cc, v);

    double dir_aux = dot_3d(v, n);
    assert(dir_aux != 0 && "Vectors are perpendicular?");
    if (dir_aux < 0) {
        n[0] = -n[0];
        n[1] = -n[1];
        n[2] = -n[2];
    } 
}


void extract_vface_BOUNDING(const s_vdiagram *vdiagram, const s_setup *setup, const s_ncell *ncell, int v_localid_vertex, int v_localid_2, s_vcell *vcell)
{   
    // DIRECTION 1 OF UNCLOSED CYCLE 
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
    
    int face_id[3];
    extract_ids_face(setup, ncell_A, &v1_A, 2, face_id);
    if (!face_extension_already_exists(vcell, face_id)) {
        double fc_A[3], n_A[3];
        extract_face_center_and_normal(setup, ncell_A, v1_A, fc_A, n_A);

        double intersection_A[3], extension_A[3];
        find_intersection_with_bounding_poly(vdiagram->bpoly, fc_A, n_A, intersection_A, extension_A);
        // printf("INTERSECTION: %f, %f, %f\n", intersection_A[0], intersection_A[1], intersection_A[2]);
        // printf("EXTENSION: %f, %f, %f\n", extension_A[0], extension_A[1], extension_A[2]);

        add_vvertex_from_coords(intersection_A, face_id, vcell);
        add_vvertex_as_extension(extension_A, face_id, vcell);
    }
    
    // DIRECTION 2 OF UNCLOSED CYCLE 
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

    extract_ids_face(setup, ncell_B, &v1_B, 2, face_id);
    // printf("DEBUG: ncell_B=%d, v1=%d, v2=%d, face=%d,%d,%d\n", ncell_B->count, ncell_B->vertex_id[v1_B], ncell_B->vertex_id[v2_B], face_id[0], face_id[1], face_id[2]);
    if (!face_extension_already_exists(vcell, face_id)) {
        double fc_B[3], n_B[3];
        extract_face_center_and_normal(setup, ncell_B, v1_B, fc_B, n_B);

        double intersection_B[3], extension_B[3];
        find_intersection_with_bounding_poly(vdiagram->bpoly, fc_B, n_B, intersection_B, extension_B);
        // printf("INTERSECTION: %f, %f, %f\n", intersection_B[0], intersection_B[1], intersection_B[2]);
        // printf("EXTENSION: %f, %f, %f\n", extension_B[0], extension_B[1], extension_B[2]);

        
        add_vvertex_from_coords(intersection_B, face_id, vcell);
        add_vvertex_as_extension(extension_B, face_id, vcell);
    }
}


void extract_vface(const s_vdiagram *vdiagram, const s_setup *setup, const s_ncell *ncell, int v_localid_vertex, int v_localid_2, s_vcell *vcell)
{
    // The ridge actually is identified with the complementary indices
    int v_ridge_1 = 0;
    while (v_ridge_1 == v_localid_vertex || v_ridge_1 == v_localid_2) {
        v_ridge_1++;
    }
    int v_ridge_2 = 0;
    while (v_ridge_2 == v_localid_vertex || v_ridge_2 == v_localid_2 || v_ridge_2 == v_ridge_1) {
        v_ridge_2++;
    }
    // printf("vertex: %d, ridge: %d, v1: %d, v2: %d\n", ncell->vertex_id[v_localid_vertex], ncell->vertex_id[v_localid_2], 
                                                    // ncell->vertex_id[v_ridge_1], ncell->vertex_id[v_ridge_2]);


    int Nv = count_cycle_ridge(setup, ncell, v_ridge_1, v_ridge_2);
    if (Nv == -1) {
        extract_vface_BOUNDING(vdiagram, setup, ncell, v_ridge_1, v_ridge_2, vcell);
    } else {
        int new_v_ridge_1, new_v_ridge_2;
        const s_ncell *current_ncell = ncell;
        for (int ii=1; ii<Nv; ii++) {
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


void remove_artificial_extended_vertices(s_vcell *vcell)
{
    int kk = 0;
    for (int ii=0; ii<vcell->Nv; ii++) {
        if (vcell->origin_vertices[ii][3] != -2) {
            vcell->vertices[ii - kk] = vcell->vertices[ii];
            vcell->origin_vertices[ii - kk] = vcell->origin_vertices[ii];
        } else {
            kk++;
        }
    }
    vcell->Nv = vcell->Nv - kk;
}


void correctly_bound_with_bp_vertices(const s_vdiagram *vdiagram, s_vcell *vcell)
{   
    int old_Nv = vcell->Nv;  // Just in case I add more points in the meantime
    const s_bound_poly *bp = vdiagram->bpoly;
    for (int ii=0; ii<bp->Np; ii++) {
        if (is_inside_convhull(bp->points[ii], vcell->vertices, old_Nv, vcell->faces, vcell->fnormals, vcell->Nf)) {
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
            segment_convex_hull_intersection(p0, p1, vcell->vertices, old_Nv, vcell->faces, 
                                                 vcell->fnormals, vcell->Nf, &ind1, i1, &ind2, i2);
            if (ind1) add_vvertex_from_segment(i1, bp->faces[ii*3+0], bp->faces[ii*3+1], ii, vcell);
            if (ind2) add_vvertex_from_segment(i2, bp->faces[ii*3+0], bp->faces[ii*3+1], ii, vcell);
        }

        if (!segment_vvertex_already_exists(vcell, bp->faces[ii*3 + 0], bp->faces[ii*3 + 2])) {
            p0 = bp->points[bp->faces[ii*3 + 0]];
            p1 = bp->points[bp->faces[ii*3 + 2]];
            segment_convex_hull_intersection(p0, p1, vcell->vertices, old_Nv, vcell->faces, 
                                                 vcell->fnormals, vcell->Nf, &ind1, i1, &ind2, i2);
            if (ind1) add_vvertex_from_segment(i1, bp->faces[ii*3+0], bp->faces[ii*3+2], ii, vcell);
            if (ind2) add_vvertex_from_segment(i2, bp->faces[ii*3+0], bp->faces[ii*3+2], ii, vcell);
        }

        if (!segment_vvertex_already_exists(vcell, bp->faces[ii*3 + 2], bp->faces[ii*3 + 1])) {
            p0 = bp->points[bp->faces[ii*3 + 2]];
            p1 = bp->points[bp->faces[ii*3 + 1]];
            segment_convex_hull_intersection(p0, p1, vcell->vertices, old_Nv, vcell->faces, 
                                                 vcell->fnormals, vcell->Nf, &ind1, i1, &ind2, i2);
            if (ind1) add_vvertex_from_segment(i1, bp->faces[ii*3+2], bp->faces[ii*3+1], ii, vcell);
            if (ind2) add_vvertex_from_segment(i2, bp->faces[ii*3+2], bp->faces[ii*3+1], ii, vcell);
        }
    }

    remove_artificial_extended_vertices(vcell);
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
    mark_ncells_incident_face(setup, ncell, v_localid_COMP, 0);

    s_vcell *vcell = malloc_vcell(vertex_id);
    
    // Extract vertices of voronoi cell
    s_ncell *current = setup->head;
    while (current) {
        if (current->mark == 1) {
            add_vvertex_from_ncell(setup, current, vcell);
            extract_vfaces_from_ncell_and_vertex(vdiagram, setup, vcell, current, vertex_id);
        }
        current = current->next;
    }


    double CM[3];
    find_center_mass(vcell->vertices, vcell->Nv, 3, CM);
    int *faces, N_faces;
    ch_vertex *vertices_aux = convert_points_to_chvertex(vcell->vertices, vcell->Nv);
    convhull_3d_build(vertices_aux, vcell->Nv, &faces, &N_faces);
    vcell->faces = faces;
    vcell->Nf = N_faces;
    double **fnormals = extract_normals_from_ch(vertices_aux, vcell->Nv, faces, N_faces, CM);
    vcell->fnormals = fnormals;
    free(vertices_aux);

    correctly_bound_with_bp_vertices(vdiagram, vcell);
    find_center_mass(vcell->vertices, vcell->Nv, 3, CM);
    vertices_aux = convert_points_to_chvertex(vcell->vertices, vcell->Nv);
    convhull_3d_build(vertices_aux, vcell->Nv, &faces, &N_faces);
    vcell->faces = faces;
    vcell->Nf = N_faces;
    double **fnormals_2 = extract_normals_from_ch(vertices_aux, vcell->Nv, faces, N_faces, CM);
    vcell->fnormals = fnormals_2;  // TODO free previous normals!
    free(vertices_aux);

    return vcell;
}


s_vdiagram *voronoi_from_delaunay_3d(const s_setup *setup, s_bound_poly *bpoly)
{
    initialize_ncells_counter(setup);

    s_vdiagram *vdiagram = malloc_vdiagram(setup);
    add_noise_to_bp(bpoly);
    vdiagram->bpoly = bpoly;
    
    for (int ii=0; ii<setup->N_points; ii++) {
        vdiagram->vcells[ii] = extract_voronoi_cell(vdiagram, setup, ii);
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

        // for (int ii=0; ii<vdiag->bpoly->Nf; ii++) {
        //     fprintf(pipe, "\"<echo \'");
        //     fprintf(pipe, "%f %f %f\\n", vdiag->bpoly->points[vdiag->bpoly->faces[ii][0]][0],
        //                                  vdiag->bpoly->points[vdiag->bpoly->faces[ii][0]][1], 
        //                                  vdiag->bpoly->points[vdiag->bpoly->faces[ii][0]][2]);
        //     fprintf(pipe, "%f %f %f\\n", vdiag->bpoly->points[vdiag->bpoly->faces[ii][1]][0],
        //                                  vdiag->bpoly->points[vdiag->bpoly->faces[ii][1]][1], 
        //                                  vdiag->bpoly->points[vdiag->bpoly->faces[ii][1]][2]);
        //     fprintf(pipe, "%f %f %f'\"", vdiag->bpoly->points[vdiag->bpoly->faces[ii][2]][0],
        //                                   vdiag->bpoly->points[vdiag->bpoly->faces[ii][2]][1], 
        //                                   vdiag->bpoly->points[vdiag->bpoly->faces[ii][2]][2]);
        //     fprintf(pipe, "w polygons fs transparent solid 0.2 notitle, ");
        // }
        fprintf(pipe, "\"<echo \'");
        fprintf(pipe, "%f %f %f\'\"", vdiag->seeds[vcell->seed_id][0], vdiag->seeds[vcell->seed_id][1], vdiag->seeds[vcell->seed_id][2]);
        fprintf(pipe, "pt 7 lc rgb 'black' notitle, ");
    }
    
    for (int ii=0; ii<vdiag->bpoly->Np; ii++) {
        fprintf(pipe, "\"<echo \'");
        fprintf(pipe, "%f %f %f\'\"", vdiag->bpoly->points[ii][0], vdiag->bpoly->points[ii][1], vdiag->bpoly->points[ii][2]);
        fprintf(pipe, "pt 7 lc rgb 'black' notitle, ");    
    }
    fprintf(pipe, "\n");
    
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
    fprintf(pipe, "%f %f %f'\"", vdiag->seeds[vcell->seed_id][0], vdiag->seeds[vcell->seed_id][1], vdiag->seeds[vcell->seed_id][2]);
    fprintf(pipe, "pt 7 lc rgb 'black' notitle, ");
    for (int ii=0; ii<vdiag->bpoly->Np; ii++) {
        fprintf(pipe, "\"<echo \'");
        fprintf(pipe, "%f %f %f\'\"", vdiag->bpoly->points[ii][0], vdiag->bpoly->points[ii][1], vdiag->bpoly->points[ii][2]);
        fprintf(pipe, "pt 7 lc rgb 'black' notitle, ");    
    }
    // for (int ii=0; ii<vdiag->bpoly->Nf; ii++) {
    //         fprintf(pipe, "\"<echo \'");
    //         fprintf(pipe, "%f %f %f\\n", vdiag->bpoly->points[vdiag->bpoly->faces[ii][0]][0],
    //                                      vdiag->bpoly->points[vdiag->bpoly->faces[ii][0]][1], 
    //                                      vdiag->bpoly->points[vdiag->bpoly->faces[ii][0]][2]);
    //         fprintf(pipe, "%f %f %f\\n", vdiag->bpoly->points[vdiag->bpoly->faces[ii][1]][0],
    //                                      vdiag->bpoly->points[vdiag->bpoly->faces[ii][1]][1], 
    //                                      vdiag->bpoly->points[vdiag->bpoly->faces[ii][1]][2]);
    //         fprintf(pipe, "%f %f %f'\"", vdiag->bpoly->points[vdiag->bpoly->faces[ii][2]][0],
    //                                       vdiag->bpoly->points[vdiag->bpoly->faces[ii][2]][1], 
    //                                       vdiag->bpoly->points[vdiag->bpoly->faces[ii][2]][2]);
    //         fprintf(pipe, "w polygons fs transparent solid 0.1 notitle, ");
    //     }
    fprintf(pipe, "\n");
    
    pclose(pipe);
}


void plot_vdiagram(s_vdiagram *vdiagram, char *f_name, double *ranges)
{
    char colors[][20] = {
        "#000090",
        "#000fff",
        "#0090ff",
        "#0fffee",
        "#90ff70",
        "#ffee00",
        "#ff7000",
        "#ee0000",
        "#7f0000"
    };

    FILE *pipe = popen("gnuplot -persistent 2>&1", "w");
    fprintf(pipe, "set terminal pngcairo enhanced font 'Arial,18' size 1080,1080 enhanced \n");
    fprintf(pipe, "set output '%s.png'\n", f_name);
    fprintf(pipe, "set pm3d depth\n");
    fprintf(pipe, "set pm3d border lc 'black' lw 0.5\n");
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
    for (int jj=0; jj<vdiagram->N_vcells; jj++) {
        s_vcell *vcell = vdiagram->vcells[jj];
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
            fprintf(pipe, "w polygons fs transparent solid 0.2 fc rgb '%s' notitle, ", colors[jj%8]);

        }
        for (int ii=0; ii<vdiagram->bpoly->Nf; ii++) {
            fprintf(pipe, "\"<echo \'");
            fprintf(pipe, "%f %f %f\\n", vdiagram->bpoly->points[vdiagram->bpoly->faces[ii*3]][0],
                                         vdiagram->bpoly->points[vdiagram->bpoly->faces[ii*3]][1], 
                                         vdiagram->bpoly->points[vdiagram->bpoly->faces[ii*3]][2]);
            fprintf(pipe, "%f %f %f\\n", vdiagram->bpoly->points[vdiagram->bpoly->faces[ii*3+1]][0],
                                         vdiagram->bpoly->points[vdiagram->bpoly->faces[ii*3+1]][1], 
                                         vdiagram->bpoly->points[vdiagram->bpoly->faces[ii*3+1]][2]);
            fprintf(pipe, "%f %f %f'\"", vdiagram->bpoly->points[vdiagram->bpoly->faces[ii*3+2]][0],
                                         vdiagram->bpoly->points[vdiagram->bpoly->faces[ii*3+2]][1], 
                                         vdiagram->bpoly->points[vdiagram->bpoly->faces[ii*3+2]][2]);
            fprintf(pipe, "w polygons fs transparent solid 0.01 fc 'black' notitle, ");
        }

        fprintf(pipe, "\"<echo \'");
        fprintf(pipe, "%f %f %f\\n", vdiagram->seeds[jj][0], vdiagram->seeds[jj][1], vdiagram->seeds[jj][2]);
        fprintf(pipe, "'\" pt 7 lc rgb 'black' notitle, ");
        for (int kk=0; kk<vdiagram->bpoly->Np; kk++) {
            fprintf(pipe, "\"<echo \'");
            fprintf(pipe, "%f %f %f\'\"", vdiagram->bpoly->points[kk][0], vdiagram->bpoly->points[kk][1], vdiagram->bpoly->points[kk][2]);
            fprintf(pipe, "pt 7 lc rgb 'black' notitle, ");    
        }
    }       
    fprintf(pipe, "\n");

    for (int jj=0; jj<vdiagram->N_vcells; jj++) {
        s_vcell *vcell = vdiagram->vcells[jj];
        fprintf(pipe, "set output '%s_%d.png'\n", f_name, jj);
            fprintf(pipe, "splot ");
        for (int ii=0; ii<vcell->Nf; ii++) {
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
            fprintf(pipe, "w polygons fs transparent solid 0.2 fc rgb '%s' notitle, ", colors[jj%8]);
            
            for (int kk=0; kk<vdiagram->N_vcells; kk++) {
                fprintf(pipe, "\"<echo \'");
                fprintf(pipe, "%f %f %f\\n", vdiagram->seeds[kk][0], vdiagram->seeds[kk][1], vdiagram->seeds[kk][2]);
                fprintf(pipe, "'\" pt 7 lc rgb 'black' notitle, ");
            }

            fprintf(pipe, "\"<echo \'");
            fprintf(pipe, "%f %f %f\\n", vdiagram->seeds[jj][0], vdiagram->seeds[jj][1], vdiagram->seeds[jj][2]);
            fprintf(pipe, "'\" pt 6 ps 2 lc rgb 'black' notitle, ");

            for (int kk=0; kk<vdiagram->bpoly->Np; kk++) {
                fprintf(pipe, "\"<echo \'");
                fprintf(pipe, "%f %f %f\'\"", vdiagram->bpoly->points[kk][0], vdiagram->bpoly->points[kk][1], vdiagram->bpoly->points[kk][2]);
                fprintf(pipe, "pt 7 lc rgb 'black' notitle, ");    
            }

        }

        for (int ii=0; ii<vdiagram->bpoly->Nf; ii++) {
            fprintf(pipe, "\"<echo \'");
            fprintf(pipe, "%f %f %f\\n", vdiagram->bpoly->points[vdiagram->bpoly->faces[ii*3]][0],
                                         vdiagram->bpoly->points[vdiagram->bpoly->faces[ii*3]][1], 
                                         vdiagram->bpoly->points[vdiagram->bpoly->faces[ii*3]][2]);
            fprintf(pipe, "%f %f %f\\n", vdiagram->bpoly->points[vdiagram->bpoly->faces[ii*3+1]][0],
                                         vdiagram->bpoly->points[vdiagram->bpoly->faces[ii*3+1]][1], 
                                         vdiagram->bpoly->points[vdiagram->bpoly->faces[ii*3+1]][2]);
            fprintf(pipe, "%f %f %f'\"", vdiagram->bpoly->points[vdiagram->bpoly->faces[ii*3+2]][0],
                                         vdiagram->bpoly->points[vdiagram->bpoly->faces[ii*3+2]][1], 
                                         vdiagram->bpoly->points[vdiagram->bpoly->faces[ii*3+2]][2]);
            fprintf(pipe, "w polygons fs transparent solid 0.05 fc 'black' notitle, ");
        }
        fprintf(pipe, "\n");
    }

    pclose(pipe);
}
