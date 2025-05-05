#include "vd_3d.h"


s_vdiagram *construct_vd_from_txt(double (*f_radius_poiss)(double *), char *file_bounding_polyhedron, int max_tries);
s_vdiagram *construct_vd_from_bp(double (*f_radius_poiss)(double *), s_bound_poly *bp, int max_tries);

// void generate_file_cube_bp(const char *filename, double length);
// void generate_file_tetrahedron_bp(const char *filename, double length);
// void generate_file_sphere_bp(const char *filename, double radius, int nTheta, int nPhi);

extern void generate_file_cube_bp(const char *filename, double length);
extern void generate_file_tetrahedron_bp(const char *filename, double length);
extern void generate_file_sphere_bp(const char *filename, double radius, int nTheta, int nPhi);


// VD_3D.C
extern void clear_volumes_file(char *fname);
extern void append_volumes_to_file(s_vdiagram *vdiagram, char *fname, int id);
extern void plot_vdiagram_auto(s_vdiagram *vdiagram, char *f_name, int max_files);

