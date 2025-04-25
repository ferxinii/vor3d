#include "vd_3d.h"


s_vdiagram *construct_vd(double (*f_radius_poiss)(double *), char *file_bounding_polyhedron, int max_tries);


// void generate_file_cube_bp(const char *filename, double length);
// void generate_file_tetrahedron_bp(const char *filename, double length);
// void generate_file_sphere_bp(const char *filename, double radius, int nTheta, int nPhi);

void generate_file_cube_bp(const char *filename, double length);
void generate_file_tetrahedron_bp(const char *filename, double length);
void generate_file_sphere_bp(const char *filename, double radius, int nTheta, int nPhi);


// VD_3D.C
void clear_volumes_file(char *fname);
void append_volumes_to_file(s_vdiagram *vdiagram, char *fname);
void plot_vdiagram_auto(s_vdiagram *vdiagram, char *f_name, int max_files);

