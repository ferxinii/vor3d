
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

    /* Dihedral at edge ij:  tan(theta_ij) = sqrt(2 M d_ij^2) / M_kl, where kl is
     * the OPPOSITE edge ({i,j,k,l} = {1,2,3,4}).  So theta_12 uses d12 and M34,
     * etc.  (A wrong variant using d_kl instead of d_ij was deleted; verified
     * against dihedral_angles_from_vertices on skewed tetrahedra.) */
    s_tetra_dihedral out;
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

/* Per-ball spherical caps of the 2-ball lens, split by the radical plane.
 * cap_A = volume of B_A on B's side of the radical plane, cap_B symmetric.
 * cap_A + cap_B == volume_intersection_2_spheres(rA,rB,dAB, false). */
static void caps_intersection_2_spheres(double rA, double rB, double dAB,
                                        double *cap_A, double *cap_B)
{
    double lAB = distance_radical_plane(rA*rA, rB*rB, dAB);
    double t1 = rA - lAB;
    *cap_A = M_PI/3.0 * (t1 * t1 * (2*rA + lAB));

    double lBA = dAB - lAB;
    double t2 = rB - lBA;
    *cap_B = M_PI/3.0 * (t2 * t2 * (2*rB + lBA));
}

double volume_intersection_2_spheres(double rA, double rB, double dAB,
                                     bool check_intersection)
{
    if (check_intersection) {
        if (dAB + rB <= rA) return volume_sphere(rB);  /* B contained in A */
        if (dAB + rA <= rB) return volume_sphere(rA);  /* A contained in B */
        if (rA + rB <= dAB) return 0;
    }

    double cap_A, cap_B;
    caps_intersection_2_spheres(rA, rB, dAB, &cap_A, &cap_B);
    return cap_A + cap_B;
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

    /* 2x2 system -- same as face_orthosphere_r2 but in explicit 2D coords */
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

/* Per-ball caps of the triple intersection, split by the radical planes.
 * out_cap[0..2] = share of B_A,B_B,B_C (the 1/3 factor is already inside, as in
 * the scalar version); out_cap[0]+out_cap[1]+out_cap[2] ==
 * volume_intersection_3_spheres(...,false,...).  Assumes the triple actually
 * intersects (the caller does the check-path early-outs). */
static void caps_intersection_3_spheres(double rA, double rB, double rC,
                                        double dAB, double dAC, double dBC,
                                        double out_cap[3])
{
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
    out_cap[0] = 1.0/3.0 * (rA*sA - lAB*D_AB - lAC*D_AC);

    /* Contribution of B (x3) */
    double sB = vol3_area_si(rB, lBA, lBC, th.th12, th.th23, th.th24);
    double D_BA = vol3_area_Dij(rB_2, lBA*lBA, th.th12);
    double D_BC = vol3_area_Dij(rB_2, lBC*lBC, th.th23);
    out_cap[1] = 1.0/3.0 * (rB*sB - lBA*D_BA - lBC*D_BC);

    /* Contribution of C (x3) */
    double sC = vol3_area_si(rC, lCA, lCB, th.th13, th.th23, th.th34);
    double D_CA = vol3_area_Dij(rC_2, lCA*lCA, th.th13);
    double D_CB = vol3_area_Dij(rC_2, lCB*lCB, th.th23);
    out_cap[2] = 1.0/3.0 * (rC*sC - lCA*D_CA - lCB*D_CB);
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

    double cap[3];
    caps_intersection_3_spheres(rA, rB, rC, dAB, dAC, dBC, cap);
    return cap[0] + cap[1] + cap[2];
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


/* PER-BALL K-TET SPLIT (power / Laguerre partition of a tetrahedron) */

/* Virtual dihedral at edge (i,m) of the 3-ball tet {c[i],c[m],c[e],P}, where P
 * is either point of the (S_i cap S_m cap S_e) intersection -- the SAME angle
 * the 3-ball intersection volume uses (distinct from the real center-tet
 * dihedral phi).  d[.][.] are center distances. */
static double virtual_dihedral_edge(const double r[4], double d[4][4],
                                    int i, int m, int e)
{
    s_tetra_lengths L = tetra_lengths_from_3_spheres(r[i], r[m], r[e],
                            d[i][m], d[i][e], d[m][e]);
    return dihedral_angles_from_lengths(L).th12;   /* th12 == edge (i,m) */
}

/* Split a K-tetrahedron's volume among its 4 balls by the power partition:
 * out_F[i] = vol F_{i;jkl} (Mach 2011 A.16 / SI Koehl S.16).  Within ball i's
 * term for neighbour m, all angles are on edge (i,m): the two VIRTUAL dihedrals
 * theta[i][m][a], theta[i][m][b] (a,b the other two vertices) and the REAL
 * dihedral phi[i][m].  Guarantees sum_i out_F[i] == fabs(signed_volume_tetra(c))
 * (identity A.26).  Individual out_F[i] MAY be negative (orthocenter outside the
 * tet) -- do NOT clamp; only the sum is guaranteed positive.  sin(phi) is safe
 * for non-degenerate tets (dihedral in (0,pi)); the caller guards
 * nc_vol > EPS_DEGEN. */
static void tet_volume_split_power(const s_point c[4], const double r[4],
                                   double EPS_DEGEN, double out_F[4])
{
    double d[4][4] = {{0}};
    for (int a = 0; a < 4; a++)
        for (int b = a+1; b < 4; b++) {
            double dab = distance(c[a], c[b]);
            d[a][b] = dab; d[b][a] = dab;
        }

    s_tetra_dihedral phi = dihedral_angles_from_vertices(c, EPS_DEGEN);

    for (int i = 0; i < 4; i++) {
        double Fi = 0.0;
        for (int m = 0; m < 4; m++) {
            if (m == i) continue;
            int rest[2], k = 0;
            for (int t = 0; t < 4; t++) if (t != i && t != m) rest[k++] = t;

            double lam = distance_radical_plane(r[i]*r[i], r[m]*r[m], d[i][m]);
            double rdisk2 = r[i]*r[i] - lam*lam;
            double ca = cos(virtual_dihedral_edge(r, d, i, m, rest[0]));
            double cb = cos(virtual_dihedral_edge(r, d, i, m, rest[1]));
            double ph = dihedral_from_edge_ids(&phi, i, m);
            Fi += (1.0/6.0) * lam * rdisk2
                  * (2*ca*cb - (ca*ca + cb*cb)*cos(ph)) / sin(ph);
        }
        out_F[i] = Fi;
    }
}


/* Shared core for the N>=4 union volume (Edelsbrunner 8.1 / the .tex theorem).
 * Returns the scalar total.  If out_contrib != NULL (length N, indexed by
 * ORIGINAL ball id) every term is additionally accumulated into its
 * contributing ball(s), and each K-tetrahedron's volume is split among its 4
 * balls by the power partition (tet_volume_split_power); then
 * sum_i out_contrib[i] == the returned total up to A.26 rounding.  When
 * out_contrib == NULL only the scalar total is computed (K-tet uses the exact
 * per-tet volume; no split rounding).  Dropped/redundant balls keep contrib 0. */
static double union_core(const s_points *centers, const double radii[],
                         double EPS_DEGEN, double TOL_dup,
                         s_dynarray *buff_ncellPTR, double *out_contrib)
{
    int N = centers->N;
    if (out_contrib) for (int i = 0; i < N; i++) out_contrib[i] = 0.0;

    /* Build weights (w_i = r_i^2) and the seed->ball identity map.  The builder
     * drops redundant/duplicate balls, so we cannot assume scplx point i+4 ==
     * input ball i; carry a remap. */
    double *weights = malloc(N * sizeof(double));
    int *ball_of_point = malloc(N * sizeof(int));   /* seed i -> compacted vid, or -1 */
    for (int i = 0; i < N; i++) { weights[i] = radii[i]*radii[i]; ball_of_point[i] = i; }

    s_dt_builder b = dt_builder_begin(centers, weights, TOL_dup, NULL, NULL);
    s_scplx dt = dt_builder_end(&b, /*keep_big_tetra=*/true, NULL, ball_of_point, N);

    /* Inverse map: compacted scplx point id -> original ball id (-1 for the four
     * sentinels and for any dropped ball). */
    int *ball_of_vid = malloc(dt.points.N * sizeof(int));
    for (int i = 0; i < dt.points.N; i++) ball_of_vid[i] = -1;
    for (int i = 0; i < N; i++)
        if (ball_of_point[i] >= 0) ball_of_vid[ball_of_point[i]] = i;

    s_hash_table ht_faces, ht_edges;
    bool *vertex_mark = malloc(dt.points.N * sizeof(bool));
    extract_alpha_complex(&dt, true, 0, buff_ncellPTR,
                          &ht_faces, &ht_edges, vertex_mark);

    /* Compute volume of union */
    double vol = 0.0;
    for (s_ncell *nc = dt.head; nc; nc = nc->next) {
        s_point p[4];  extract_vertices_ncell(&dt, nc, p);
        double nc_vol = fabs(signed_volume_tetra(p));
        if (nc_vol <= EPS_DEGEN) continue;

        /* FIRST SUM: tets in K.  Their 4 vertices are all real balls (a tet
         * touching a sentinel has a huge orthoradius and is never in K). */
        if (nc->mask_alpha) {
            vol += nc_vol;
            if (out_contrib) {
                double rF[4];
                for (int v = 0; v < 4; v++) {
                    int bv = ball_of_vid[nc->vertex_id[v]];
                    assert(bv >= 0);
                    rF[v] = radii[bv];
                }
                double F[4];
                tet_volume_split_power(p, rF, EPS_DEGEN, F);
                for (int v = 0; v < 4; v++)
                    out_contrib[ball_of_vid[nc->vertex_id[v]]] += F[v];
            }
            continue;
        }

        /* SECOND SUM: tets NOT in K */
        s_tetra_dihedral  dihe = dihedral_angles_from_vertices(p, EPS_DEGEN);
        s_tetra_solid_angles sa = solid_angles_from_dihedral(dihe);

        /* FACES */
        for (int omit = 0; omit < 4; omit++) {
            int face[3]; extract_ids_face(nc, 2, &omit, face); sort3(face);
            bool *face_in_K = hash_get(&ht_faces, face);
            assert(face_in_K);
            if (!(*face_in_K)) continue;

            int ba = ball_of_vid[face[0]], bb = ball_of_vid[face[1]], bc = ball_of_vid[face[2]];
            assert(ba >= 0 && bb >= 0 && bc >= 0);
            double cap3[3];
            caps_intersection_3_spheres(radii[ba], radii[bb], radii[bc],
                            distance(dt.points.p[face[0]], dt.points.p[face[1]]),
                            distance(dt.points.p[face[0]], dt.points.p[face[2]]),
                            distance(dt.points.p[face[1]], dt.points.p[face[2]]), cap3);
            vol += 0.5 * (cap3[0] + cap3[1] + cap3[2]);
            if (out_contrib) {
                out_contrib[ba] += 0.5 * cap3[0];
                out_contrib[bb] += 0.5 * cap3[1];
                out_contrib[bc] += 0.5 * cap3[2];
            }
        }

        /* EDGES */
        for (int i=0; i<3; i++) for (int j=i+1; j<4; j++) {
            int edge[2] = {nc->vertex_id[i], nc->vertex_id[j]}; sort2(edge);
            bool *edge_in_K = hash_get(&ht_edges, edge);
            assert(edge_in_K);
            if (!(*edge_in_K)) continue;

            int be0 = ball_of_vid[edge[0]], be1 = ball_of_vid[edge[1]];
            assert(be0 >= 0 && be1 >= 0);
            double theta = dihedral_from_edge_ids(&dihe, i, j);
            double coef = theta / (2*M_PI);
            double cap0, cap1;
            caps_intersection_2_spheres(radii[be0], radii[be1],
                            distance(dt.points.p[edge[0]], dt.points.p[edge[1]]),
                            &cap0, &cap1);
            vol -= coef * (cap0 + cap1);
            if (out_contrib) {
                out_contrib[be0] -= coef * cap0;
                out_contrib[be1] -= coef * cap1;
            }
        }

        /* VERTICES: v_localid is the vertex local id (0-3) */
        for (int v=0; v<4; v++) {
            if (!vertex_mark[nc->vertex_id[v]]) continue;
            int bv = ball_of_vid[nc->vertex_id[v]];
            assert(bv >= 0);
            double omega = solid_angle_from_vertex_id(&sa, v);
            double cv = (omega / (4*M_PI)) * volume_sphere(radii[bv]);
            vol += cv;
            if (out_contrib) out_contrib[bv] += cv;
        }
    }

    free(weights);
    free(ball_of_point);
    free(ball_of_vid);
    free_complex(&dt);
    hash_free(&ht_faces);
    hash_free(&ht_edges);
    free(vertex_mark);
    return vol;
}


/* 2-ball caps WITH the check_intersection early-outs of
 * volume_intersection_2_spheres: a contained (redundant) ball takes the whole
 * lens, disjoint balls take nothing.  cap_A + cap_B ==
 * volume_intersection_2_spheres(rA,rB,dAB, true). */
static void caps_intersection_2_spheres_checked(double rA, double rB, double dAB,
                                                double *cap_A, double *cap_B)
{
    if (dAB + rB <= rA) { *cap_A = 0.0;              *cap_B = volume_sphere(rB); return; } /* B in A */
    if (dAB + rA <= rB) { *cap_A = volume_sphere(rA); *cap_B = 0.0;             return; } /* A in B */
    if (rA + rB <= dAB) { *cap_A = 0.0;              *cap_B = 0.0;              return; } /* disjoint */
    caps_intersection_2_spheres(rA, rB, dAB, cap_A, cap_B);
}

/* 3-ball caps WITH the check_intersection early-outs of
 * volume_intersection_3_spheres (same branch order / reductions): under
 * pairwise containment the triple reduces to the lens of the surviving two.
 * cap[0]+cap[1]+cap[2] == volume_intersection_3_spheres(...,true,...). */
static void caps_intersection_3_spheres_checked(double rA, double rB, double rC,
        double dAB, double dAC, double dBC, double EPS_DEGEN, double cap[3])
{
    cap[0] = cap[1] = cap[2] = 0.0;
    if (dAB + rB <= rA) { caps_intersection_2_spheres_checked(rB, rC, dBC, &cap[1], &cap[2]); return; }
    if (dAB + rA <= rB) { caps_intersection_2_spheres_checked(rA, rC, dAC, &cap[0], &cap[2]); return; }
    if (dAC + rC <= rA) { caps_intersection_2_spheres_checked(rB, rC, dBC, &cap[1], &cap[2]); return; }
    if (dAC + rA <= rC) { caps_intersection_2_spheres_checked(rA, rB, dAB, &cap[0], &cap[1]); return; }
    if (dBC + rC <= rB) { caps_intersection_2_spheres_checked(rA, rC, dAC, &cap[0], &cap[2]); return; }
    if (dBC + rB <= rC) { caps_intersection_2_spheres_checked(rA, rB, dAB, &cap[0], &cap[1]); return; }
    if (dAB >= rA + rB) return;
    if (dAC >= rA + rC) return;
    if (dBC >= rB + rC) return;
    if (!has_triple_intersection(rA, rB, rC, dAB, dAC, dBC, EPS_DEGEN)) return;
    caps_intersection_3_spheres(rA, rB, rC, dAB, dAC, dBC, cap);
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
    return union_core(centers, radii, EPS_DEGEN, TOL_dup, buff_ncellPTR, NULL);
}


void volume_contribution_spheres(const s_points *centers, const double radii[centers->N],
                                 double EPS_DEGEN, double TOL_dup,
                                 s_dynarray *buff_ncellPTR, double out_contrib[centers->N])
{
    int N = centers->N;
    for (int i = 0; i < N; i++) out_contrib[i] = 0.0;

    /* N=1,2,3 mirror the analytic totals of volume_union_spheres term by term,
     * splitting every inclusion-exclusion intersection into its per-ball caps.
     * Because the operations match exactly, sum_i out_contrib[i] equals the
     * analytic union total, and the per-ball values are the power (Laguerre)
     * partition (a redundant/contained ball gets 0). */
    if (N == 1) {
        out_contrib[0] = volume_sphere(radii[0]);
        return;
    }
    if (N == 2) {
        double cA, cB;
        caps_intersection_2_spheres_checked(radii[0], radii[1],
                distance(centers->p[0], centers->p[1]), &cA, &cB);
        out_contrib[0] = volume_sphere(radii[0]) - cA;
        out_contrib[1] = volume_sphere(radii[1]) - cB;
        return;
    }
    if (N == 3) {
        double rA = radii[0], rB = radii[1], rC = radii[2];
        const s_point *P = centers->p;
        double dAB = distance(P[0], P[1]);
        double dAC = distance(P[0], P[2]);
        double dBC = distance(P[1], P[2]);
        out_contrib[0] = volume_sphere(rA);
        out_contrib[1] = volume_sphere(rB);
        out_contrib[2] = volume_sphere(rC);
        double ca, cb;
        caps_intersection_2_spheres_checked(rA, rB, dAB, &ca, &cb); out_contrib[0] -= ca; out_contrib[1] -= cb;
        caps_intersection_2_spheres_checked(rA, rC, dAC, &ca, &cb); out_contrib[0] -= ca; out_contrib[2] -= cb;
        caps_intersection_2_spheres_checked(rB, rC, dBC, &ca, &cb); out_contrib[1] -= ca; out_contrib[2] -= cb;
        double c3[3];
        caps_intersection_3_spheres_checked(rA, rB, rC, dAB, dAC, dBC, EPS_DEGEN, c3);
        out_contrib[0] += c3[0]; out_contrib[1] += c3[1]; out_contrib[2] += c3[2];
        return;
    }

    /* N >= 4: general RT + alpha-complex path.
     * sum_i out_contrib[i] == volume_union_spheres total. */
    union_core(centers, radii, EPS_DEGEN, TOL_dup, buff_ncellPTR, out_contrib);
}

