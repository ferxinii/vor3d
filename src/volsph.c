
#include "volsph.h"
#include "points.h"
#include "scplx.h"
#include "delaunay.h"
#include "linalg.h"
#include "hash.h"
#include <math.h>
#include <float.h>
#include <assert.h>


/* STRUCTS FOR TETRAHEDRA INFORMATION */
typedef struct tetra_lengths {
    /* Squared l*/
    double d12_2;
    double d13_2;
    double d14_2; 
    double d23_2;
    double d24_2;
    double d34_2; 
} s_tetra_lengths;

typedef struct tetra_dihedral {
    double th12;
    double th13;
    double th14;
    double th23;
    double th24;
    double th34;
} s_tetra_dihedral;

typedef struct tetra_solid_angles {
    double o1;
    double o2;
    double o3;
    double o4;
} s_tetra_solid_angles;



/* CAYLEY-MENGER AND ANGLES */
typedef enum CM_minor { CM_M12,
    CM_M13,
    CM_M14,
    CM_M23,
    CM_M24,
    CM_M34
} e_CM_minor;

static double abs_cayley_menger_minor(s_tetra_lengths i, e_CM_minor minor)
{   /* abs means that I do not multiply by (-1)^(i+j), not that it is always >0 */
    double dij, dik, dil, djk, djl, dkl;
    switch (minor) {
    case CM_M12: 
        dij=i.d12_2; dik=i.d13_2; dil=i.d14_2; djk=i.d23_2; djl=i.d24_2; dkl=i.d34_2; break;
    case CM_M13: 
        dij=i.d13_2; dik=i.d12_2; dil=i.d14_2; djk=i.d23_2; djl=i.d34_2; dkl=i.d24_2; break;
    case CM_M14: 
        dij=i.d14_2; dik=i.d12_2; dil=i.d13_2; djk=i.d24_2; djl=i.d34_2; dkl=i.d23_2; break;
    case CM_M23: 
        dij=i.d23_2; dik=i.d12_2; dil=i.d24_2; djk=i.d13_2; djl=i.d34_2; dkl=i.d14_2; break;
    case CM_M24: 
        dij=i.d24_2; dik=i.d12_2; dil=i.d23_2; djk=i.d14_2; djl=i.d34_2; dkl=i.d13_2; break;
    case CM_M34: 
        dij=i.d34_2; dik=i.d13_2; dil=i.d23_2; djk=i.d14_2; djl=i.d24_2; dkl=i.d12_2; break;
    }

    return ((dij - dik - djk)*(dij - dil - djl) - (dij + dkl - dik - djl)*(dij + dkl - djk - dil));
}

static double cayley_menger_det(s_tetra_lengths in)
{   
    double t12 = (in.d13_2 + in.d14_2 - in.d34_2);
    double t13 = (in.d12_2 + in.d14_2 - in.d24_2);
    double t14 = (in.d12_2 + in.d13_2 - in.d23_2);
    double M = 8 * in.d12_2*in.d13_2*in.d14_2 
              - 2 * (in.d12_2*t12*t12 + in.d13_2*t13*t13 + in.d14_2*t14*t14) 
              + 2 * t12*t13*t14;
    return M;
}

static s_tetra_dihedral dihedral_angles_from_lengths(s_tetra_lengths in)
{
    double M_2 = 2 * cayley_menger_det(in);
    double M12 = abs_cayley_menger_minor(in, CM_M12);
    double M13 = abs_cayley_menger_minor(in, CM_M13);
    double M14 = abs_cayley_menger_minor(in, CM_M14);
    double M23 = abs_cayley_menger_minor(in, CM_M23);
    double M24 = abs_cayley_menger_minor(in, CM_M24);
    double M34 = abs_cayley_menger_minor(in, CM_M34);

    s_tetra_dihedral out;
    out.th12 = atan2(sqrt(M_2 * in.d34_2), M34);
    out.th13 = atan2(sqrt(M_2 * in.d24_2), M24);
    out.th14 = atan2(sqrt(M_2 * in.d23_2), M23);
    out.th23 = atan2(sqrt(M_2 * in.d14_2), M14);
    out.th24 = atan2(sqrt(M_2 * in.d13_2), M13);
    out.th34 = atan2(sqrt(M_2 * in.d12_2), M12);

    out.th12 = atan2(sqrt(M_2 * in.d12_2), M34);
    out.th13 = atan2(sqrt(M_2 * in.d13_2), M24);
    out.th14 = atan2(sqrt(M_2 * in.d14_2), M23);
    out.th23 = atan2(sqrt(M_2 * in.d23_2), M14);
    out.th24 = atan2(sqrt(M_2 * in.d24_2), M13);
    out.th34 = atan2(sqrt(M_2 * in.d34_2), M12);

    return out;
}

/* In the case where all 4 vertices are known, dihedral angles are easier */
static double dihedral_angle_at_edge(s_point a, s_point b, s_point c, s_point d,
                                     double EPS_DEGEN)
{   /* Dihedral angle at edge a-b, c and d are the two opposite vertices */
    s_point edge_dir = normalize_vec(subtract_points(b, a), EPS_DEGEN);
    if (!point_is_valid(edge_dir)) return 0.0;

    /* Project c and d onto plane perpendicular to edge through a */
    s_point ac = subtract_points(c, a);
    s_point ad = subtract_points(d, a);
    s_point vc = subtract_points(ac, scale_point(edge_dir, dot_prod(ac, edge_dir)));
    s_point vd = subtract_points(ad, scale_point(edge_dir, dot_prod(ad, edge_dir)));

    return atan2(norm(cross_prod(vc, vd)), dot_prod(vc, vd));
}

static s_tetra_dihedral dihedral_angles_from_vertices(const s_point p[4], double EPS_DEGEN)
{
    return (s_tetra_dihedral){
        .th12 = dihedral_angle_at_edge(p[0], p[1], p[2], p[3], EPS_DEGEN),
        .th13 = dihedral_angle_at_edge(p[0], p[2], p[1], p[3], EPS_DEGEN),
        .th14 = dihedral_angle_at_edge(p[0], p[3], p[1], p[2], EPS_DEGEN),
        .th23 = dihedral_angle_at_edge(p[1], p[2], p[0], p[3], EPS_DEGEN),
        .th24 = dihedral_angle_at_edge(p[1], p[3], p[0], p[2], EPS_DEGEN),
        .th34 = dihedral_angle_at_edge(p[2], p[3], p[0], p[1], EPS_DEGEN),
    };
}


static s_tetra_solid_angles solid_angles_from_dihedral(s_tetra_dihedral in)
{
    s_tetra_solid_angles out;
    out.o1 = in.th12 + in.th13 + in.th14 - M_PI;
    out.o2 = in.th12 + in.th23 + in.th24 - M_PI;
    out.o3 = in.th13 + in.th23 + in.th34 - M_PI;
    out.o4 = in.th14 + in.th24 + in.th34 - M_PI;
    return out;
}



// static s_tetra_lengths tetra_lengths_from_vertices(const s_point p[4])
// {
//     s_tetra_lengths out;
//     out.d12_2 = distance_squared(p[0], p[1]);
//     out.d13_2 = distance_squared(p[0], p[2]);
//     out.d14_2 = distance_squared(p[0], p[3]);
//     out.d23_2 = distance_squared(p[1], p[2]);
//     out.d24_2 = distance_squared(p[1], p[3]);
//     out.d34_2 = distance_squared(p[2], p[3]);
//     return out;
// }



/* VOLUME INTERSECTIONS */
static double distance_radical_plane(double rA_2, double rB_2, double dAB)
{
    return 0.5 * (dAB + (rA_2 - rB_2)/dAB);
}

double volume_intersection_2_spheres(double rA, double rB, double dAB,
                                     bool check_intersection)
{
    if (check_intersection) {
        if (dAB + rB <= rA) return volume_sphere(rB);  /* B contained in A */
        if (dAB + rA <= rB) return volume_sphere(rA);  /* A contained in B */
        if (rA + rB <= dAB) return 0;
    }

    /* Volume of spherical cap of A */
    double lAB = distance_radical_plane(rA*rA, rB*rB, dAB);
    double t1 = rA - lAB;
    double VA = t1 * t1 * (2*rA + lAB);

    /* Volume of spherical cap of B */
    double lBA = dAB - lAB;
    double t2 = rB - lBA;
    double VB = t2 * t2 * (2*rB + lBA);

    return ( M_PI/3.0*(VA + VB) );
}


static double vol3_area_Dij(double ri_2, double lij_2, double th_ij)
{
    return ( (ri_2 - lij_2) * (th_ij - sin(th_ij)*cos(th_ij)) );
}


static double vol3_area_si(double ri, double lij, double lik, double th_ij, double th_ik, double psi_i)
{
    double t1 = (ri - lij) * th_ij;
    double t2 = (ri - lik) * th_ik;
    double t3 = ri * (th_ij + th_ik + psi_i - M_PI);
    return ( 2 * ri * (t1 + t2 - t3) );
}


static s_tetra_lengths tetra_lengths_from_3_spheres(double rA, double rB, double rC, double dAB, double dAC, double dBC)
{
    s_tetra_lengths out;
    out.d12_2 = dAB*dAB;
    out.d13_2 = dAC*dAC;
    out.d14_2 = rA*rA;
    out.d23_2 = dBC*dBC;
    out.d24_2 = rB*rB;
    out.d34_2 = rC*rC;
    return out;
}

static bool has_triple_intersection(double rA, double rB, double rC,
                                    double dAB, double dAC, double dBC,
                                    double EPS_DEGEN)
{   /* Checks if triple intersection happens */
    /* Place A at origin, B on x-axis */
    double Bx = dAB;
    double Cx = (dAB*dAB + dAC*dAC - dBC*dBC) / (2.0*dAB);
    double Cy2 = dAC*dAC - Cx*Cx;
    if (Cy2 < 0.0) return false;  /* degenerate: A,B,C collinear */
    double Cy = sqrt(Cy2);

    double q[3][2] = {{0, 0}, {Bx, 0}, {Cx, Cy}};
    double w[3]    = {rA*rA, rB*rB, rC*rC};

    /* 2x2 system — same as face_orthosphere_r2 but in explicit 2D coords */
    double A[2][2], b[2];
    for (int i = 1; i <= 2; i++) {
        A[i-1][0] = 2.0 * q[i][0];   /* q[0] = (0,0) so terms drop out */
        A[i-1][1] = 2.0 * q[i][1];
        b[i-1] = q[i][0]*q[i][0] + q[i][1]*q[i][1] - (w[i] - w[0]);
    }
    double c[2];
    if (solve_2x2_cramer(A, b, c, EPS_DEGEN) != 2) return false;

    return (c[0]*c[0] + c[1]*c[1] - w[0] < 0.0);
}

double volume_intersection_3_spheres(double rA, double rB, double rC, double dAB, double dAC, double dBC, 
                                     bool check_intersection, double EPS_DEGEN)
{   
    if (check_intersection) {
        /* Pairwise containment, must come first before checking intersection existence */
        if (dAB + rB <= rA) return volume_intersection_2_spheres(rB, rC, dBC, true);
        if (dAB + rA <= rB) return volume_intersection_2_spheres(rA, rC, dAC, true);
        if (dAC + rC <= rA) return volume_intersection_2_spheres(rB, rC, dBC, true);
        if (dAC + rA <= rC) return volume_intersection_2_spheres(rA, rB, dAB, true);
        if (dBC + rC <= rB) return volume_intersection_2_spheres(rA, rC, dAC, true);
        if (dBC + rB <= rC) return volume_intersection_2_spheres(rA, rB, dAB, true);
        /* Check all pairs intersect */
        if (dAB >= rA + rB) return 0.0;
        if (dAC >= rA + rC) return 0.0;
        if (dBC >= rB + rC) return 0.0;
        /* Check triple intersection */
        if (!has_triple_intersection(rA, rB, rC, dAB, dAC, dBC, EPS_DEGEN)) return 0.0;
    }
    

    /* Build tetrahedron (4th point is either of the two intersections).*/
    s_tetra_lengths d = tetra_lengths_from_3_spheres(rA, rB, rC, dAB, dAC, dBC);
    s_tetra_dihedral th = dihedral_angles_from_lengths(d);
    double rA_2 = d.d14_2, rB_2 = d.d24_2, rC_2 = d.d34_2;  /* Already computed squared radii */

    /* Distances to radical planes */
    double lAB = distance_radical_plane(rA_2, rB_2, dAB);
    double lAC = distance_radical_plane(rA_2, rC_2, dAC);
    double lBC = distance_radical_plane(rB_2, rC_2, dBC);
    double lBA = dAB - lAB;
    double lCA = dAC - lAC;
    double lCB = dBC - lBC;

    /* Contribution of A (x3) */
    double sA = vol3_area_si(rA, lAB, lAC, th.th12, th.th13, th.th14);
    double D_AB = vol3_area_Dij(rA_2, lAB*lAB, th.th12);
    double D_AC = vol3_area_Dij(rA_2, lAC*lAC, th.th13);
    double VA = rA*sA - lAB*D_AB - lAC*D_AC;

    /* Contribution of B (x3) */
    double sB = vol3_area_si(rB, lBA, lBC, th.th12, th.th23, th.th24);
    double D_BA = vol3_area_Dij(rB_2, lBA*lBA, th.th12);
    double D_BC = vol3_area_Dij(rB_2, lBC*lBC, th.th23);
    double VB = rB*sB - lBA*D_BA - lBC*D_BC;

    /* Contribution of C (x3) */
    double sC = vol3_area_si(rC, lCA, lCB, th.th13, th.th23, th.th34);
    double D_CA = vol3_area_Dij(rC_2, lCA*lCA, th.th13);
    double D_CB = vol3_area_Dij(rC_2, lCB*lCB, th.th23);
    double VC = rC*sC - lCA*D_CA - lCB*D_CB;

    return ( 1.0/3.0*(VA + VB + VC) );
}



/* VOLUME UNION */
double volume_sphere(double r)
{
    return 4.0/3.0*M_PI*r*r*r;
}

double volume_union_2_spheres(s_point pA, double rA,
                              s_point pB, double rB,
                              bool check_intersection)
{   /* Inclusion-Exclusion principle */
    double VA = volume_sphere(rA);
    double VB = volume_sphere(rB);

    double dAB = distance(pA, pB);
    double V_iAB = volume_intersection_2_spheres(rA, rB, dAB, check_intersection);
    return VA + VB - V_iAB;
}

double volume_union_3_spheres(s_point pA, double rA,
                              s_point pB, double rB,
                              s_point pC, double rC,
                              bool check_intersection, double EPS_DEGEN)
{   /* Inclusion-Exclusion principle */
    double VA = volume_sphere(rA);
    double VB = volume_sphere(rB);
    double VC = volume_sphere(rC);

    double dAB = distance(pA, pB);
    double dAC = distance(pA, pC);
    double dBC = distance(pB, pC);

    double V_iAB = volume_intersection_2_spheres(rA, rB, dAB, check_intersection);
    double V_iAC = volume_intersection_2_spheres(rA, rC, dAC, check_intersection);
    double V_iBC = volume_intersection_2_spheres(rB, rC, dBC, check_intersection);

    double V_iABC = volume_intersection_3_spheres(rA, rB, rC, dAB, dAC, dBC, 
                                                  check_intersection, EPS_DEGEN);

    return VA + VB + VC - V_iAB - V_iAC - V_iBC + V_iABC;
}


static double dihedral_from_edge_ids(const s_tetra_dihedral *dihe, int i, int j)
{   /* Note that i,j are 0-indexed*/
    if (i > j) { int tmp=i; i=j; j=tmp; }
    if (i==0 && j==1) return dihe->th12;
    if (i==0 && j==2) return dihe->th13;
    if (i==0 && j==3) return dihe->th14;
    if (i==1 && j==2) return dihe->th23;  
    if (i==1 && j==3) return dihe->th24;
    return dihe->th34;  /* if (i==2 && j==3) */
}


static double solid_angle_from_vertex_id(const s_tetra_solid_angles *sa, int v)
{   /* v is local id 0-3, solid angles are 1-indexed in the struct */
    switch(v) {
        case 0: return sa->o1;
        case 1: return sa->o2;
        case 2: return sa->o3;
        default: return sa->o4;
    }
}


static inline void sort3(int *a)
{
    if (a[0] > a[1]) { int t=a[0]; a[0]=a[1]; a[1]=t; }
    if (a[1] > a[2]) { int t=a[1]; a[1]=a[2]; a[2]=t; }
    if (a[0] > a[1]) { int t=a[0]; a[0]=a[1]; a[1]=t; }
}

static inline void sort2(int *a)
{
    if (a[0] > a[1]) { int t=a[0]; a[0]=a[1]; a[1]=t; }
}


double volume_union_spheres(const s_points *centers, const double radii[centers->N],
                            double EPS_DEGEN, double TOL_dup, s_dynarray *buff_ncellPTR)
{
    int N = centers->N;
    if (N == 1) {
        return volume_sphere(radii[0]);
    } else if (N == 2) {
        return volume_union_2_spheres(centers->p[0], radii[0], 
                                      centers->p[1], radii[1], 
                                      true);
    } else if (N == 3) {
        return volume_union_3_spheres(centers->p[0], radii[0],
                                      centers->p[1], radii[1],
                                      centers->p[2], radii[2],
                                      true, EPS_DEGEN);
    }

    double weights[N];  /* Build weights from radii */
    for (int i=0; i<N; i++) weights[i] = radii[i]*radii[i];

    /* Build Delaunay triangulation and mark alpha complex */
    s_scplx dt = construct_dt_3d(centers, weights, true, TOL_dup);

    s_hash_table ht_faces, ht_edges;
    bool *vertex_mark = malloc(dt.points.N * sizeof(bool));
    extract_alpha_complex(&dt, true, 0, EPS_DEGEN, buff_ncellPTR, 
                          &ht_faces, &ht_edges, vertex_mark);

    /* Compute volume of union */
    double vol = 0.0;
    for (s_ncell *nc = dt.head; nc; nc = nc->next) {
        s_point p[4];  extract_vertices_ncell(&dt, nc, p);
        double nc_vol = fabs(signed_volume_tetra(p));
        if (nc_vol <= EPS_DEGEN) continue;

        /* FIRST SUM: tets in K */
        if (nc->mask_alpha) { vol += nc_vol; continue; }

        /* SECOND SUM: tets NOT in K */
        s_tetra_dihedral  dihe = dihedral_angles_from_vertices(p, EPS_DEGEN);
        s_tetra_solid_angles sa = solid_angles_from_dihedral(dihe);

        /* FACES */
        for (int omit = 0; omit < 4; omit++) {
            int face[3]; extract_ids_face(nc, 2, &omit, face); sort3(face);
            bool *face_in_K = hash_get(&ht_faces, face);
            assert(face_in_K);
            if (!(*face_in_K)) continue; 

            vol += 0.5 * volume_intersection_3_spheres(
                            radii[face[0]-4], radii[face[1]-4], radii[face[2]-4],
                            distance(dt.points.p[face[0]], dt.points.p[face[1]]),
                            distance(dt.points.p[face[0]], dt.points.p[face[2]]),
                            distance(dt.points.p[face[1]], dt.points.p[face[2]]),
                            false, EPS_DEGEN);
        }

        /* EDGES */
        for (int i=0; i<3; i++) for (int j=i+1; j<4; j++) {
            int edge[2] = {nc->vertex_id[i], nc->vertex_id[j]}; sort2(edge);
            bool *edge_in_K = hash_get(&ht_edges, edge);
            assert(edge_in_K);
            if (!(*edge_in_K)) continue; 

            double theta = dihedral_from_edge_ids(&dihe, i, j);
            vol -= (theta / (2*M_PI)) 
                    * volume_intersection_2_spheres(radii[edge[0]-4], radii[edge[1]-4], 
                            distance(dt.points.p[edge[0]], dt.points.p[edge[1]]), false);
        }

        /* VERTICES: v_localid is the vertex local id (0-3) */
        for (int v=0; v<4; v++) {
            if (!vertex_mark[nc->vertex_id[v]]) continue;
            double omega = solid_angle_from_vertex_id(&sa, v); 
            vol += (omega / (4*M_PI)) * volume_sphere(radii[nc->vertex_id[v]-4]);
        }
    }

    free_complex(&dt);
    hash_free(&ht_faces);
    hash_free(&ht_edges);
    free(vertex_mark);
    return vol;
}

