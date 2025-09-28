#ifndef PLOT_H
#define PLOT_H

/* -----------------------------------------
 *              GNUPLOTC
 *
 *  Description: C Interface to gnuplot.
 *  Author: Fernando Muñoz.
 *  Github: ferxinii/gnuplotC
 *  License: MIT
 * --------------------------------------- */

typedef struct t_gnuplot t_gnuplot;  // Interface to gnuplot

enum gnuplot_type {
    PNG_2D,
    PDF_2D,
    EPS_2D,
    VIDEO_2D,
    PNG_3D,
    PDF_3D,
    VIDEO_3D
};

extern int GNUPLOTC_FRAMERATE;

// BASIC
t_gnuplot *gnuplot_start(enum gnuplot_type type, char *file_name, int size[2], int font_size, ...);
void gnuplot_config(t_gnuplot *interface, ...);
void gnuplot_end(t_gnuplot *interface);

// VIDEO
void activate_parallel_video_processing(int num_threads);
void next_frame(t_gnuplot *interface, ...); 
void video_to_gif(const char *file_video, const char *file_gif, int size[2], int framerate);

// GENERIC ELEMENTS
void add_datablock_from_file(t_gnuplot *interface, const char *filename, int header_lines_skip, const char *datablock_name);
void draw_datablock(t_gnuplot *interface, const char *datablock_name, const char *config);
void draw_file(t_gnuplot *interface, const char *file_name, const char *config);
void draw_command(t_gnuplot *interface, const char *command);

// 2D ELEMENTS:
void draw_point_2d(t_gnuplot *interface, double x, double y, const char *config);
void draw_2d(t_gnuplot *interface, double *x, double *y, int N, const char *config); 
void draw_array_2d(t_gnuplot *interface, double **coords, int N, const char *config);
void draw_arrows_from_array_2d(t_gnuplot *interface, double **coords, int N, int spacing, int offset, const char *config);  // Recall size command for head size (and fixed keyword), ex: size 0.2,20 fixed
void draw_errorbars_2d(t_gnuplot *interface, double *x, double *mean, double *err, int N, const char *config);
void draw_segment_2d(t_gnuplot *interface, double x0, double xf, double y0, double yf, const char *config);
void draw_function_2d(t_gnuplot *interface, double x0, double xf, int N, double (*fun)(double), const char *config);

// 3D ELEMENTS:
void draw_sphere_3d(t_gnuplot *interface, double x, double y, double z, double r, const char *type, const char *config);
void draw_solid_triangle_3d(t_gnuplot *interface, double v0[3], double v1[3], double v2[3], const char *config);





/* ---------------------------------
    PRIVATE DECLARATIONS (IGNORE)
  --------------------------------- */

#define TEMPLATES_DIR "GNUPLOTC_TEMPLATES_TMP"
#define FRAMES_DIR "GNUPLOTC_FRAMES_TMP"

#define GNUPLOTC_ARRAY_MARKER  ((char*) -1)

#define GNUPLOTC_ARRAY(arr) \
        GNUPLOTC_ARRAY_MARKER, arr

#define gnuplot_start(obj, ...) \
        gnuplot_start_impl(obj, __VA_OPT__(__VA_ARGS__,) NULL)

t_gnuplot *gnuplot_start_impl(enum gnuplot_type type, char *file_name, int size[2], int font_size, ...);

#define gnuplot_config(obj, ...) \
        gnuplot_config_impl(obj, __VA_OPT__(__VA_ARGS__,) NULL)

void gnuplot_config_impl(t_gnuplot *interface, ...);

#define next_frame(obj, ...) \
        next_frame_impl(obj, __VA_OPT__(__VA_ARGS__,) NULL)

void next_frame_impl(t_gnuplot *interface, ...);

#endif
