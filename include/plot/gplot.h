/* 
 * Header-only gnuplot interface.
 * Allows making plots and videos through an easy interface to gnuplot.
 * 
 * Copyright (c) 2026 Fernando Muñoz
 * MIT license. See bottom of file.
 *
 * REQUIREMENTS
 *  - The macros use two C23 features (__VA_OPT__ and empty {});
 *    only -Wpedantic below C23 warns about them.
 *    Tested on macOS (clang 17): default, -std=c11, -std=c17 and -std=c2x.
 *  - gnuplot (with the pngcairo and pdfcairo terminals) in PATH. Videos and GIFs also need ffmpeg.
 *
 *  Linux with a strict -std=cXX, not gnuXX (NOT TESTED, Claude's recommendation): glibc hides the
 *  POSIX functions above in strict modes, so define this before ANY #include:
 *      #define _POSIX_C_SOURCE 200809L
 *  On macOS that define hides mkdtemp instead; there, either don't define it, or also add
 *  #define _DARWIN_C_SOURCE (tested).
 */

#ifndef HLIBS_GPLOT_H
#define HLIBS_GPLOT_H
#include <stdio.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include <dirent.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <sys/types.h>

/* EXAMPLE
    double x[100], y[100], z[100];  // ...fill them

    s_gplot *g = gplot_start(GPLOT_PNG_2D, "plot.png", .size = {800, 600});
    gplot_config(g, "set grid", "set xlabel 'x'");
    gplot_draw_2d(g, 100, x, y, .title = "y");
    gplot_draw_2d(g, 100, x, z, .with = GPLOT_POINTS, .pt = 7, .color = "red");

    s_gplot_style st = { .lw = 2, .dt = 2 };  // Reuse one style: call the _impl version
    gplot_draw_2d_impl(g, 100, y, x, st);
    gplot_end(g);

    s_gplot *v = gplot_start(GPLOT_VIDEO_2D, "anim.mp4", .framerate = 30);
    gplot_config(v, "set yrange [-1:1]");  // Commands persist until changed, also across frames
    for (int f = 0; f < 60; f++) {
        if (f > 0) gplot_next_frame(v, "set title 'later frames'");  // gplot_start already began frame 1
        gplot_draw_2d(v, 100, x, y);
    }
    gplot_end(v);
    // GPLOT_GIF_2D / GPLOT_GIF_3D work the same and write an animated GIF ("anim.gif")
*/


/* ---------------------------------
    PUBLIC: TYPES
  --------------------------------- */

typedef struct gplot s_gplot;  // Interface to gnuplot

typedef enum gplot_plot {
    GPLOT_PNG_2D,
    GPLOT_PDF_2D,
    GPLOT_VIDEO_2D,       // .mp4 (H.264)
    GPLOT_GIF_2D,         // animated GIF, one palette for all frames
    GPLOT_PNG_3D,
    GPLOT_PDF_3D,
    GPLOT_VIDEO_3D,
    GPLOT_GIF_3D
} e_gplot_plot;

typedef struct gplot_conf {  // Zero fields get defaults: 800x600, font 12, 24 fps, layout 1x1
    int size[2];
    int font_size;
    int framerate;
    int video_parallel_threads;
    int layout[2];
} s_gplot_conf;

typedef enum gplot_elem {
    GPLOT_DEFAULT,  // 0: the draw function picks its natural style
    GPLOT_LINES,
    GPLOT_POINTS,
    GPLOT_LINESPOINTS,
    GPLOT_POLYGONS,
    GPLOT_PM3D
} e_gplot_elem;

typedef struct gplot_style {  // Zero / NULL fields are not emitted, so gnuplot's defaults apply
    e_gplot_elem with;
    double lw;           // line width, 0 = unset
    int dt;              // dash type, 0 = unset
    int pt;              // point type, 0 = unset (GPLOT_PT_DOT for gnuplot's pt 0).
                         // Run "test" in gnuplot to see the shapes for each number.
    double ps;           // point size, 0 = unset
    const char *color;   // "red", "#ff8800"; NULL = unset
    const char *title;   // NULL = notitle
    const char *extra;   // raw gnuplot style options, appended last
} s_gplot_style;

#define GPLOT_PT_DOT -1  // gnuplot's pt 0 (a dot), since 0 means "unset"


/* ---------------------------------
    PUBLIC: FUNCTIONS
    Call shapes. ...conf / ...style are optional fields of s_gplot_conf / s_gplot_style:
        gplot_start(GPLOT_PNG_2D, "a.png", .size = {800, 600})
        gplot_draw_2d(g, N, x, y, .with = GPLOT_POINTS, .title = "a")
    Every call taking ...style has a gplot_X_impl(..., s_gplot_style style) to pass a style variable.
    cmd, ... = any number of gnuplot commands, and GPLOT_ARRAY(cmds) for a NULL-terminated array.
    Commands persist until changed, like in plain gnuplot (also across frames and subplots).
    gnuplot exits at the first invalid command; later calls on that plot then do nothing.
    So that this cannot kill your program, gplot ignores SIGPIPE, unless you set your own handler.
  --------------------------------- */

// START / END
// s_gplot *gplot_start(e_gplot_plot type, const char *file, ...conf);
// void     gplot_end(s_gplot *g);

// CONFIG AND SUBPLOTS
// void gplot_config(s_gplot *g, const char *cmd, ...);
// void gplot_next_subplot(s_gplot *g, const char *cmd, ...);

// VIDEO AND GIF
// void gplot_next_frame(s_gplot *g, const char *cmd, ...);  // Not before the first frame
// A GIF's palette is computed from all frames, so ffmpeg holds them all in memory until gplot_end.

// DRAW: GENERIC
// void gplot_draw_command(s_gplot *g, const char *command);  // Raw plot clause, ex: "t, sin(t) w l notitle"
// void gplot_draw_file(s_gplot *g, const char *file, const char *modifiers, ...style);  // modifiers go before "w": skip, using, every, index

// DRAW: 2D
// void gplot_draw_point_2d(s_gplot *g, double x, double y, ...style);
// void gplot_draw_2d(s_gplot *g, int N, const double x[N], const double y[N], ...style);
// void gplot_draw_function_2d(s_gplot *g, double x0, double xf, int N, double (*fun)(double), ...style);
// void gplot_draw_errorbars_2d(s_gplot *g, int N, const double x[N], const double mean[N], const double err[N], ...style);  // (fixed style)
// void gplot_draw_arrows_2d(s_gplot *g, int N, const double x[N], const double y[N], int spacing, int offset, ...style);  // (fixed style) Head size via .extra, ex: "size 0.2,20 fixed"
// void gplot_draw_segment_2d(s_gplot *g, double x0, double y0, double xf, double yf, ...style);  // (fixed style)

// DRAW: 3D
// void gplot_draw_sphere_3d(s_gplot *g, double x, double y, double z, double r, ...style);
// void gplot_draw_solid_triangle_3d(s_gplot *g, const double v0[3], const double v1[3], const double v2[3], ...style);  // (fixed style)
// void gplot_draw_polytope_3d(s_gplot *g, int N, const double x[N], const double y[N], const double z[N], ...style);  // (fixed style)
// void gplot_draw_polygons_3d(s_gplot *g, int N_polys, const int offset[N_polys+1],
//                             const double x[], const double y[], const double z[],
//                             const unsigned rgb[N_polys], ...style);  // (fixed style)
//      Many polygons in one clause (fast). Polygon i has vertices offset[i] .. offset[i+1]-1.
//      rgb: one 0xRRGGBB fill color per polygon, or NULL for a single color (.color).
//      For opaque, correctly hidden faces: gplot_config(g, "set pm3d depthorder", "set style fill solid 1 border lc 'black'")

// (fixed style): the data only fits one style, so .with is ignored with a warning




/* ---------------------------------
    PRIVATE: MACROS AND PROTOTYPES BEHIND THE PUBLIC CALLS
  --------------------------------- */

#define GPLOT_ARRAY(arr) \
        GPLOT_ARRAY_MARKER, arr

#define gplot_start(type, file, ...) \
        gplot_start_impl(type, file, (s_gplot_conf){ __VA_ARGS__ })
#define gplot_config(g, ...) \
        gplot_config_impl(g, __VA_OPT__(__VA_ARGS__,) NULL)
#define gplot_next_subplot(g, ...) \
        gplot_next_subplot_impl(g, __VA_OPT__(__VA_ARGS__,) NULL)
#define gplot_next_frame(g, ...) \
        gplot_next_frame_impl(g, __VA_OPT__(__VA_ARGS__,) NULL)

#define gplot_draw_file(g, file, modifiers, ...) \
        gplot_draw_file_impl(g, file, modifiers, (s_gplot_style){ __VA_ARGS__ })
#define gplot_draw_point_2d(g, x, y, ...) \
        gplot_draw_point_2d_impl(g, x, y, (s_gplot_style){ __VA_ARGS__ })
#define gplot_draw_2d(g, N, x, y, ...) \
        gplot_draw_2d_impl(g, N, x, y, (s_gplot_style){ __VA_ARGS__ })
#define gplot_draw_function_2d(g, x0, xf, N, fun, ...) \
        gplot_draw_function_2d_impl(g, x0, xf, N, fun, (s_gplot_style){ __VA_ARGS__ })
#define gplot_draw_errorbars_2d(g, N, x, mean, err, ...) \
        gplot_draw_errorbars_2d_impl(g, N, x, mean, err, (s_gplot_style){ __VA_ARGS__ })
#define gplot_draw_arrows_2d(g, N, x, y, spacing, offset, ...) \
        gplot_draw_arrows_2d_impl(g, N, x, y, spacing, offset, (s_gplot_style){ __VA_ARGS__ })
#define gplot_draw_segment_2d(g, x0, y0, xf, yf, ...) \
        gplot_draw_segment_2d_impl(g, x0, y0, xf, yf, (s_gplot_style){ __VA_ARGS__ })
#define gplot_draw_sphere_3d(g, x, y, z, r, ...) \
        gplot_draw_sphere_3d_impl(g, x, y, z, r, (s_gplot_style){ __VA_ARGS__ })
#define gplot_draw_solid_triangle_3d(g, v0, v1, v2, ...) \
        gplot_draw_solid_triangle_3d_impl(g, v0, v1, v2, (s_gplot_style){ __VA_ARGS__ })
#define gplot_draw_polytope_3d(g, N, x, y, z, ...) \
        gplot_draw_polytope_3d_impl(g, N, x, y, z, (s_gplot_style){ __VA_ARGS__ })
#define gplot_draw_polygons_3d(g, N_polys, offset, x, y, z, rgb, ...) \
        gplot_draw_polygons_3d_impl(g, N_polys, offset, x, y, z, rgb, (s_gplot_style){ __VA_ARGS__ })

static inline s_gplot *gplot_start_impl(e_gplot_plot type, const char *file_name, s_gplot_conf conf);
static inline void gplot_end(s_gplot *interface);
static inline void gplot_config_impl(s_gplot *interface, ...);
static inline void gplot_next_subplot_impl(s_gplot *interface, ...);
static inline void gplot_next_frame_impl(s_gplot *interface, ...);

static inline void gplot_draw_command(s_gplot *interface, const char *command);
static inline void gplot_draw_file_impl(s_gplot *interface, const char *file,
                                        const char *modifiers, s_gplot_style style);
static inline void gplot_draw_point_2d_impl(s_gplot *interface, double x, double y, s_gplot_style style);
static inline void gplot_draw_2d_impl(s_gplot *interface,
                                      int N, const double x[N], const double y[N],
                                      s_gplot_style style);
static inline void gplot_draw_function_2d_impl(s_gplot *interface,
                                               double x0, double xf, int N, double (*fun)(double),
                                               s_gplot_style style);
static inline void gplot_draw_errorbars_2d_impl(s_gplot *interface,
                                                int N, const double x[N], const double mean[N], const double err[N],
                                                s_gplot_style style);
static inline void gplot_draw_arrows_2d_impl(s_gplot *interface,
                                             int N, const double x[N], const double y[N],
                                             int spacing, int offset, s_gplot_style style);
static inline void gplot_draw_segment_2d_impl(s_gplot *interface,
                                              double x0, double y0, double xf, double yf,
                                              s_gplot_style style);
static inline void gplot_draw_sphere_3d_impl(s_gplot *interface,
                                             double x, double y, double z, double r,
                                             s_gplot_style style);
static inline void gplot_draw_solid_triangle_3d_impl(s_gplot *interface,
                                                     const double v0[3], const double v1[3], const double v2[3],
                                                     s_gplot_style style);
static inline void gplot_draw_polytope_3d_impl(s_gplot *interface,
                                               int N, const double x[N], const double y[N], const double z[N],
                                               s_gplot_style style);
static inline void gplot_draw_polygons_3d_impl(s_gplot *interface, int N_polys, const int offset[],
                                               const double x[], const double y[], const double z[],
                                               const unsigned rgb[], s_gplot_style style);


/* ---------------------------------
    PRIVATE: INTERNALS (IGNORE)
  --------------------------------- */

#ifndef GPLOT_GNUPLOT_CMD
#define GPLOT_GNUPLOT_CMD "gnuplot /dev/stdin"
#endif

#define GPLOT_DBL_FMT "%.17g"

#define GPLOT_ARRAY_MARKER  ((char*) -1)  // Sentinel placed before an array by GPLOT_ARRAY


// Main interface
struct gplot {
    FILE *pipe;
    bool failed;          // gnuplot died: every later call does nothing

    e_gplot_plot type;
    int dim;
    const char *terminal;
    bool is_video;        // Also true for GIFs: frames go through ffmpeg
    bool is_gif;

    char *file_name;
    int font_size;
    int size[2];
    int layout[2];

    char *frames_dir;
    int frames_N;
    int framerate;
    int video_parallel_threads;
    FILE **procs;
    int worker;           // Index of the current pipe in procs (always 0 without workers)

    FILE *cmd_log;        // Every gnuplot command sent so far (they persist, like in plain gnuplot)
    char *cmd_log_buf;
    size_t cmd_log_len;
    size_t *cmd_sent;     // Per worker: how much of cmd_log it has received

    bool in_multiplot;    // "set multiplot" sent for the current frame
    int N_items;
    FILE *plot_cmd; 
    char *plot_cmd_buf;
    size_t plot_cmd_len;
};


/* ---------------------------------
    PRIVATE: SMALL HELPERS
  --------------------------------- */


static inline void gplot_pipe_and_newline(FILE *out, const char *cmd)
{
    fputs(cmd, out);
    fputs("\n", out);
}


static inline bool gplot_ok(s_gplot *interface)
{   // false on NULL, or once gnuplot has died. Reports the death once
    // Writes are buffered, so a death shows up a few calls late. Those calls are lost anyway
    if (!interface || interface->failed) return false;
    if (ferror(interface->pipe)) {
        fprintf(stderr, "gplot ERROR: gnuplot exited, probably after an invalid command "
                        "(see its message above). Ignoring further calls\n");
        interface->failed = true;
        return false;
    }
    return true;
}


static inline const char* gplot_elem_str(e_gplot_elem type)
{
    switch(type) {
        case GPLOT_LINES: return "lines";
        case GPLOT_POINTS: return "points";
        case GPLOT_LINESPOINTS: return "linespoints";
        case GPLOT_POLYGONS: return "polygons";
        case GPLOT_PM3D: return "pm3d";
        default: {
            fprintf(stderr, "gplot ERROR: gplot_elem_str missing case\n");
            return NULL;
        }
    }
}


static inline FILE *gplot_popen_gnuplot(void)
{   // NULL on error
    // gnuplot reads from default stdin one byte at a time (SLOW)
    // Opening it as a file buffers the input and is WAY faster
    void (*old)(int) = signal(SIGPIPE, SIG_IGN);  // A dead gnuplot makes writes fail, not kill the program
    if (old != SIG_DFL) {
        signal(SIGPIPE, old);  // Keep the user's choice
        static bool warned = false;
        if (old != SIG_IGN && !warned) {
            fprintf(stderr, "gplot WARNING: your program has its own SIGPIPE handler. gplot keeps it, "
                            "so it will run if gnuplot exits on an invalid command\n");
            warned = true;
        }
    }

    FILE *pipe = popen(GPLOT_GNUPLOT_CMD, "w");
    if (!pipe) {
        fprintf(stderr, "gplot ERROR: could not open pipe to gnuplot\n");
        return NULL;
    }
    return pipe;
}


static inline void gplot_remove_directory(const char *path) 
{
    DIR *dir = opendir(path);
    if (!dir) return;  // Directory does not exist

    struct dirent *entry;
    while ((entry = readdir(dir)) != NULL) {
        if (!strcmp(entry->d_name, ".") || !strcmp(entry->d_name, "..")) continue;  // Skip . and ..

        char fname[512];
        snprintf(fname, 512, "%s/%s", path, entry->d_name);

        struct stat statbuf;
        if (stat(fname, &statbuf) == 0) {
            // if (S_ISDIR(statbuf.st_mode)) { r2 = remove_directory(buf); } // recursive call
            if (S_ISREG(statbuf.st_mode)) {  // Test if regular file
                if(unlink(fname) != 0) { // remove file
                    fprintf(stderr, "gplot ERROR: could not remove file %s\n", fname);
                }
            } else {
                fprintf(stderr, "gplot ERROR: %s is not a regular file, not removed\n", fname);
            }
        }
    }
    closedir(dir);

    if (rmdir(path) != 0) {  // Directory must be empty before removing it!
        fprintf(stderr, "gplot ERROR: could not remove directory %s\n", path);
    }
}


static inline bool gplot_is2d(s_gplot *interface)
{
    if (interface->dim == 2) return true;
    fprintf(stderr, "gplot ERROR: Cannot use 2D method.\n");
    return false;
}


static inline bool gplot_is3d(s_gplot *interface)
{
    if (interface->dim == 3) return true;
    fprintf(stderr, "gplot ERROR: Cannot use 3D method.\n");
    return false;
}


/* ---------------------------------
    PRIVATE: PLOT COMMAND BUFFER (clauses, datablocks, styles)
  --------------------------------- */


static inline void gplot_pipe_style(FILE *out, s_gplot_style s, const char *default_with)
{   // Writes " w <with>" plus only the style fields that are set, and "notitle" unless a title is given
    fprintf(out, " w %s", s.with ? gplot_elem_str(s.with) : default_with);
    if (s.lw > 0)                  fprintf(out, " lw %g", s.lw);
    if (s.dt > 0)                  fprintf(out, " dt %d", s.dt);
    if (s.pt > 0)                  fprintf(out, " pt %d", s.pt);
    else if (s.pt == GPLOT_PT_DOT) fputs(" pt 0", out);
    if (s.ps > 0)                  fprintf(out, " ps %g", s.ps);
    if (s.color)                   fprintf(out, " lc rgb \"%s\"", s.color);
    if (s.extra)                   fprintf(out, " %s", s.extra);
    if (s.title)                   fprintf(out, " title \"%s\"", s.title);
    else                           fputs(" notitle", out);
}


static inline void gplot_ignore_with(s_gplot_style *style, const char *func_name)
{   // For draw functions whose data only fits one style (errorbars, vectors, polygons)
    if (style->with == GPLOT_DEFAULT) return;
    fprintf(stderr, "gplot WARNING: %s ignores .with\n", func_name);
    style->with = GPLOT_DEFAULT;
}


static inline int gplot_begin_clause(s_gplot *interface)
{   // Writes ", \" + newline before every clause except the first. Returns this element's index.
    // This is needed because gnuplot reading from a file (/dev/stding) cuts lines longer than 2048 characters
    if (interface->N_items > 0) fputs(", \\\n", interface->plot_cmd);
    return interface->N_items++;  // Returns the value BEFORE incrementing
}


static inline void gplot_end_clause(s_gplot *interface, s_gplot_style style, const char *default_with)
{   // Finishes the clause opened by gplot_begin_clause: " w <with> <style fields> [notitle]"
    gplot_pipe_style(interface->plot_cmd, style, default_with);
}


static inline void gplot_begin_datablock(s_gplot *interface)
{
    fprintf(interface->pipe, "$gplot_%d << EOD\n", interface->N_items);
}


static inline void gplot_end_datablock(s_gplot *interface, s_gplot_style style, const char *default_with)
{
    gplot_pipe_and_newline(interface->pipe, "EOD");

    int N = gplot_begin_clause(interface);
    fprintf(interface->plot_cmd, "$gplot_%d", N); 

    gplot_end_clause(interface, style, default_with);
}


static inline bool gplot_sync_cmds(s_gplot *interface)
{   // Sends the commands the current gnuplot has not received yet. A parallel worker
    // also catches up on the ones sent while other workers had the frames
    if (fflush(interface->cmd_log) != 0) return false;  // Updates cmd_log_buf and cmd_log_len
    size_t *sent = &interface->cmd_sent[interface->worker];
    fwrite(interface->cmd_log_buf + *sent, 1, interface->cmd_log_len - *sent, interface->pipe);
    *sent = interface->cmd_log_len;
    return true;
}


static inline void gplot_open_multiplot(s_gplot *interface)
{   // The frame's multiplot opens at its first plot, after the pending commands. gnuplot undoes
    // some settings (margins) made inside a multiplot at "unset multiplot", so commands given
    // before a frame's first plot must reach gnuplot outside it, on every worker alike
    if (interface->in_multiplot) return;
    fprintf(interface->pipe, "set multiplot layout %d,%d\n", interface->layout[0], interface->layout[1]);
    interface->in_multiplot = true;
}


static inline bool gplot_flush_plot(s_gplot *interface)
{   // Sends "plot <clauses>" to gnuplot, after any pending commands, and empties the clause buffer
    if (!gplot_sync_cmds(interface)) return false;
    gplot_open_multiplot(interface);
    const char *keyword = (interface->dim == 2) ? "plot " : "splot ";

    if (interface->N_items == 0) {  // Empty frame: one transparent point, so axes and title still exist
        fprintf(interface->pipe, "%s%s lc rgb '#FF000000' notitle\n",
                keyword, interface->dim == 2 ? "0, 0" : "0, 0, 0");
        return true;
    }

    if (fflush(interface->plot_cmd) != 0) return false;  // Updates plot_cmd_buf and plot_cmd_len
    fputs(keyword, interface->pipe);
    fwrite(interface->plot_cmd_buf, 1, interface->plot_cmd_len, interface->pipe);
    fputc('\n', interface->pipe);

    rewind(interface->plot_cmd);  // The next plot overwrites the buffer from the start
    interface->N_items = 0;
    return true;
}


/* ---------------------------------
    PRIVATE: FRAMES
  --------------------------------- */


static inline bool gplot_set_output(s_gplot *interface, const char *output)
{
    char template[] = "set terminal %s enhanced font 'Arial,%d' size %d,%d\nset output '%s'\n";
    int o = fprintf(interface->pipe,
                    template, interface->terminal, interface->font_size, 
                    interface->size[0], interface->size[1], output);

    return (o >= 0);
}


static inline bool gplot_set_frame_output(s_gplot *interface)
{   // Parallel video: points gnuplot at frames_dir/frame_<frames_N>.png
    char path[1024];
    int o = snprintf(path, sizeof(path), "%s/frame_%05d.png", interface->frames_dir, interface->frames_N);
    if (o < 0 || (size_t)o >= sizeof(path)) {
        fprintf(stderr, "gplot ERROR: frame path too long\n");
        return false;
    }
    return gplot_set_output(interface, path);
}


static inline void gplot_frame_end(s_gplot *interface)
{
    if (!gplot_flush_plot(interface)) {
        fprintf(stderr, "gplot ERROR: gplot_frame_end\n");
    }
    gplot_pipe_and_newline(interface->pipe, "unset multiplot");
    interface->in_multiplot = false;
}


/* ---------------------------------
    IMPLEMENTATION: LIFECYCLE (start, config, subplots, frames, end)
  --------------------------------- */


static inline s_gplot_conf gplot_process_conf(s_gplot_conf conf)
{   // 0 gets overwritten by the default value
    s_gplot_conf out = conf;
    if (conf.size[0] <= 0)   { out.size[0] = 800; }
    if (conf.size[1] <= 0)   { out.size[1] = 600; }
    if (conf.font_size <= 0)  { out.font_size = 12; }
    if (conf.framerate <= 0) { out.framerate = 24; } 
    if (conf.video_parallel_threads < 0) { out.video_parallel_threads = 0; }
    if (conf.layout[0] <= 0) { out.layout[0] = 1; }
    if (conf.layout[1] <= 0) { out.layout[1] = 1; }
    return out;
}


static inline void gplot_free(s_gplot *interface)
{   // Frees memory only. The caller closes gnuplot (and workers) and removes the frames directory
    if (!interface) return;
    if (interface->plot_cmd) fclose(interface->plot_cmd);
    free(interface->plot_cmd_buf);
    if (interface->cmd_log) fclose(interface->cmd_log);
    free(interface->cmd_log_buf);
    free(interface->cmd_sent);
    free(interface->file_name);
    free(interface->frames_dir);
    free(interface->procs);
    free(interface);
}


static inline bool gplot_set_type(s_gplot *interface, e_gplot_plot type)
{   // Terminal, dimension and video flag of each plot type. false for an unknown type
    interface->type = type;
    switch (type) {
        case GPLOT_PNG_2D:   { interface->terminal = "pngcairo"; interface->dim = 2; interface->is_video = false; break; }
        case GPLOT_PDF_2D:   { interface->terminal = "pdfcairo"; interface->dim = 2; interface->is_video = false; break; }
        case GPLOT_PNG_3D:   { interface->terminal = "pngcairo"; interface->dim = 3; interface->is_video = false; break; }
        case GPLOT_PDF_3D:   { interface->terminal = "pdfcairo"; interface->dim = 3; interface->is_video = false; break; }
        case GPLOT_VIDEO_2D: { interface->terminal = "pngcairo"; interface->dim = 2; interface->is_video = true;  break; }
        case GPLOT_VIDEO_3D: { interface->terminal = "pngcairo"; interface->dim = 3; interface->is_video = true;  break; }
        case GPLOT_GIF_2D:   { interface->terminal = "pngcairo"; interface->dim = 2; interface->is_video = true;  interface->is_gif = true; break; }
        case GPLOT_GIF_3D:   { interface->terminal = "pngcairo"; interface->dim = 3; interface->is_video = true;  interface->is_gif = true; break; }
        default: { fprintf(stderr, "gplot ERROR: UNIMPLEMENTED TYPE\n"); return false; }
    }
    return true;
}


static inline bool gplot_make_frames_dir(s_gplot *interface)
{   // Parallel video: a temporary directory for the PNG frames, and the worker slots (worker 0 is the current pipe)
    const char *tmp = getenv("TMPDIR");
    if (!tmp || !*tmp) tmp = ".";

    char path[1024];
    int o = snprintf(path, sizeof(path), "%s/gplot_XXXXXX", tmp);
    if (o < 0 || (size_t)o >= sizeof(path)) return false;
    if (!mkdtemp(path)) return false;
    interface->frames_dir = strdup(path);
    if (!interface->frames_dir) { rmdir(path); return false; }

    interface->procs = calloc(interface->video_parallel_threads, sizeof(FILE *));
    if (!interface->procs) return false;
    interface->procs[0] = interface->pipe;
    return true;
}


static inline const char *gplot_ffmpeg_encoder(const s_gplot *interface)
{   // GIF: one palette computed from all frames; video: H.264
    return interface->is_gif
        ? "-filter_complex \"split[a][b];[a]palettegen[p];[b][p]paletteuse\" -loop 0"
        : "-pix_fmt yuv420p -c:v libx264 -crf 18";
}


static inline bool gplot_set_ffmpeg_output(s_gplot *interface)
{   // Serial video: gnuplot writes every PNG into one ffmpeg process
    char cmd[1024];
    int o = snprintf(cmd, sizeof(cmd),
                     "| ffmpeg -y -loglevel error "
                     "-f png_pipe -r %d -i pipe: %s \"%s\"",
                     interface->framerate, gplot_ffmpeg_encoder(interface), interface->file_name);
    if (o < 0 || (size_t)o >= sizeof(cmd)) {
        fprintf(stderr, "gplot ERROR: ffmpeg command too long\n");
        return false;
    }
    return gplot_set_output(interface, cmd);
}


static inline s_gplot *gplot_start_impl(e_gplot_plot type, const char *file_name, s_gplot_conf conf)
{
    conf = gplot_process_conf(conf);

    s_gplot *out = calloc(1, sizeof(s_gplot));
    if (!out) goto error;
    if (!gplot_set_type(out, type)) goto error;
    out->file_name = strdup(file_name);
    if (!out->file_name) goto error;
    out->size[0] = conf.size[0];
    out->size[1] = conf.size[1];
    out->font_size = conf.font_size;
    out->frames_N = 1;
    out->framerate = conf.framerate;
    out->video_parallel_threads = conf.video_parallel_threads;
    out->layout[0] = conf.layout[0];
    out->layout[1] = conf.layout[1];
    out->plot_cmd = open_memstream(&out->plot_cmd_buf, &out->plot_cmd_len);
    if (!out->plot_cmd) goto error;
    out->cmd_log = open_memstream(&out->cmd_log_buf, &out->cmd_log_len);
    if (!out->cmd_log) goto error;
    int N_workers = (out->is_video && out->video_parallel_threads > 0) ? out->video_parallel_threads : 1;
    out->cmd_sent = calloc(N_workers, sizeof(size_t));
    if (!out->cmd_sent) goto error;
    gplot_pipe_and_newline(out->cmd_log, "set parametric");  // Default commands go here
    out->pipe = gplot_popen_gnuplot();
    if (!out->pipe) goto error;

    // Where the first frame goes: the file, a PNG in the frames directory, or ffmpeg
    bool ok;
    if (!out->is_video)                       ok = gplot_set_output(out, out->file_name);
    else if (out->video_parallel_threads > 0) ok = gplot_make_frames_dir(out) && gplot_set_frame_output(out);
    else                                      ok = gplot_set_ffmpeg_output(out);
    if (!ok) goto error;

    return out;

error:
    fprintf(stderr, "gplot ERROR: Error in gplot_start_impl\n");
    if (out && out->pipe) { pclose(out->pipe); }
    if (out && out->frames_dir) { rmdir(out->frames_dir); }
    gplot_free(out);
    return NULL;
}


static inline void gplot_pipe_config_vargs(FILE *out, va_list ap)
{   // Pipes configuration variable arguments which may be strings or a NULL-terminated
    // array of strings (preceded by GPLOT_ARRAY_MARKER, which GPLOT_ARRAY adds).
    // Any order and as many as needed; the last argument must be NULL (the macros add it).
    char *cmd = va_arg(ap, char *);
    while(cmd != NULL) {
        if (cmd == GPLOT_ARRAY_MARKER) {
            char **cmd_array = va_arg(ap, char **);
            for (int ii=0; cmd_array[ii] != NULL; ii++) {
                gplot_pipe_and_newline(out, cmd_array[ii]);
            }
        } else {
            gplot_pipe_and_newline(out, cmd);
        }
        cmd = va_arg(ap, char *);
    }
}


static inline void gplot_config_impl(s_gplot *interface, ...)
{
    if (!gplot_ok(interface)) return;
    va_list ap;
    va_start(ap, interface);
    gplot_pipe_config_vargs(interface->cmd_log, ap);  // Sent before the next plot
    va_end(ap);
}


static inline void gplot_next_subplot_impl(s_gplot *interface, ...)
{
    if (!gplot_ok(interface)) return;
    if (!gplot_flush_plot(interface)) {
        fprintf(stderr, "gplot ERROR: gplot_next_subplot_impl\n");
    }

    va_list ap;
    va_start(ap, interface);
    gplot_pipe_config_vargs(interface->cmd_log, ap);
    va_end(ap);
}


static inline bool gplot_switch_worker(s_gplot *interface)
{   // Parallel video: hands the finished frame to its worker and moves pipe to the worker of frame frames_N
    fflush(interface->pipe);  // Hand the finished frame to its worker now
    int slot = (interface->frames_N - 1) % interface->video_parallel_threads;
    if (!interface->procs[slot]) interface->procs[slot] = gplot_popen_gnuplot();
    if (!interface->procs[slot]) return false;
    interface->worker = slot;
    interface->pipe = interface->procs[slot];  // A busy worker blocks the write until it catches up
    return gplot_set_frame_output(interface);
}


static inline void gplot_next_frame_impl(s_gplot *interface, ...)
{
    if (!gplot_ok(interface)) return;
    gplot_frame_end(interface);
    interface->frames_N++;

    if (interface->procs && !gplot_switch_worker(interface)) {  // procs only exists for parallel video
        fprintf(stderr, "gplot ERROR: could not hand frame %d to a worker\n", interface->frames_N);  
    }

    va_list ap;
    va_start(ap, interface);
    gplot_pipe_config_vargs(interface->cmd_log, ap);
    va_end(ap);
}


static inline void gplot_frames_to_video(s_gplot *interface)
{
    char buff[1024];
    snprintf(buff, 1024,
             "ffmpeg -y -loglevel error -framerate %d -i \"%s/frame_%%05d.png\" %s \"%s\"",
             interface->framerate, interface->frames_dir, gplot_ffmpeg_encoder(interface), interface->file_name);
    system(buff);
}


static inline void gplot_end(s_gplot *interface)
{
    if (!interface) return;

    if (gplot_ok(interface)) {
        gplot_frame_end(interface);
        fputs("set output\n", interface->pipe);  // Closes the output file explicitly (waits for ffmpeg to finish)
        fflush(interface->pipe);
    }
    bool alive = gplot_ok(interface);  // Also catches a death during the last frame

    if (interface->procs) {  // Parallel video (procs only exists then). pipe is one of procs
        for (int i=0; i<interface->video_parallel_threads; i++) {
            if (interface->procs[i]) pclose(interface->procs[i]);  // Waits for that worker's last PNG
        }
        if (alive) gplot_frames_to_video(interface);
        gplot_remove_directory(interface->frames_dir);
    } else {
        pclose(interface->pipe);
    }

    gplot_free(interface);
}


/* ---------------------------------
    IMPLEMENTATION: DRAW FUNCTIONS
  --------------------------------- */


// GENERIC ELEMENTS


static inline void gplot_draw_command(s_gplot *interface, const char *cmd)
{
    if (!gplot_ok(interface) || !cmd) return;
    gplot_begin_clause(interface);
    fputs(cmd, interface->plot_cmd);
}


static inline void gplot_draw_file_impl(s_gplot *interface, const char *file,
                                        const char *modifiers, s_gplot_style style)
{
    if (!gplot_ok(interface) || !file) return;
    gplot_begin_clause(interface);
    fprintf(interface->plot_cmd, "'%s'", file);
    if (modifiers && modifiers[0] != '\0') fprintf(interface->plot_cmd, " %s", modifiers);
    gplot_end_clause(interface, style, "lines");
}


// 2D ELEMENTS


static inline void gplot_draw_point_2d_impl(s_gplot *interface, double x, double y, s_gplot_style style)
{
    if (!gplot_ok(interface)) return;
    if (!gplot_is2d(interface)) return;

    gplot_begin_datablock(interface);
    fprintf(interface->pipe,
            GPLOT_DBL_FMT " " GPLOT_DBL_FMT "\n",
            x, y);
    gplot_end_datablock(interface, style, "points");
}


static inline void gplot_draw_2d_impl(s_gplot *interface,
                                      int N, const double x[N], const double y[N],
                                      s_gplot_style style)
{
    if (!gplot_ok(interface)) return;
    if (!gplot_is2d(interface)) return;
    if (N <= 0) return;

    gplot_begin_datablock(interface);
    for (int i=0; i<N; i++) {
        fprintf(interface->pipe,
                GPLOT_DBL_FMT " " GPLOT_DBL_FMT "\n",
                x[i], y[i]);
    }
    gplot_end_datablock(interface, style, "lines");
}


static inline void gplot_draw_errorbars_2d_impl(s_gplot *interface,
                                                int N, const double x[N], const double mean[N], const double err[N],
                                                s_gplot_style style)
{
    if (!gplot_ok(interface)) return;
    if (!gplot_is2d(interface)) return;
    if (N <= 0) return;
    gplot_ignore_with(&style, "gplot_draw_errorbars_2d");

    gplot_begin_datablock(interface);
    for (int i=0; i<N; i++) {
        fprintf(interface->pipe,
                GPLOT_DBL_FMT " " GPLOT_DBL_FMT " " GPLOT_DBL_FMT " " GPLOT_DBL_FMT "\n",
                x[i], mean[i], mean[i] + err[i], mean[i] - err[i]);
    }
    gplot_end_datablock(interface, style, "errorbars");
}


static inline void gplot_draw_arrows_2d_impl(s_gplot *interface,
                                             int N, const double x[N], const double y[N],
                                             int spacing, int offset, s_gplot_style style)
{
    if (!gplot_ok(interface)) return;
    if (!gplot_is2d(interface)) return;
    if (N < 2 || offset < 0 || offset >= N-1) return;
    gplot_ignore_with(&style, "gplot_draw_arrows_2d");

    if (spacing <= 0) spacing = 1;

    gplot_begin_datablock(interface);
    for (int i = offset; i < N - 1; i += spacing) {
        fprintf(interface->pipe,
                GPLOT_DBL_FMT " " GPLOT_DBL_FMT " " GPLOT_DBL_FMT " " GPLOT_DBL_FMT "\n",
                x[i], y[i], x[i+1]-x[i], y[i+1]-y[i]);

    }
    gplot_end_datablock(interface, style, "vectors");
}


static inline void gplot_draw_segment_2d_impl(s_gplot *interface,
                                              double x0, double y0, double xf, double yf,
                                              s_gplot_style style)
{
    if (!gplot_ok(interface)) return;
    if (!gplot_is2d(interface)) return;
    gplot_ignore_with(&style, "gplot_draw_segment_2d");

    gplot_begin_datablock(interface);
    fprintf(interface->pipe,
            GPLOT_DBL_FMT " " GPLOT_DBL_FMT " " GPLOT_DBL_FMT " " GPLOT_DBL_FMT "\n",
            x0, y0, xf-x0, yf-y0);
    gplot_end_datablock(interface, style, "vectors nohead");
}


static inline void gplot_draw_function_2d_impl(s_gplot *interface,
                                               double x0, double xf, int N, double (*fun)(double),
                                               s_gplot_style style)
{
    if (!gplot_ok(interface)) return;
    if (!gplot_is2d(interface)) return;
    if (N < 2) return;

    gplot_begin_datablock(interface);
    for (int i=0; i<N; i++) {
        double x = x0 + i * (xf - x0) / (N - 1);
        fprintf(interface->pipe,
                GPLOT_DBL_FMT " " GPLOT_DBL_FMT "\n",
                x, fun(x));
    }
    gplot_end_datablock(interface, style, "lines");
}


// 3D ELEMENTS


static inline void gplot_draw_sphere_3d_impl(s_gplot *interface,
                                             double x, double y, double z, double r,
                                             s_gplot_style style)
{   // .with: pm3d, lines, points, ... (polygons? Should produce nothing if haven't "set pm3d", else, same as pm3d)
    if (!gplot_ok(interface)) return;
    if (!gplot_is3d(interface)) return;

    gplot_begin_clause(interface);
    fprintf(interface->plot_cmd,
            "[u=-pi/2:pi/2] [v=0:2*pi] " 
            GPLOT_DBL_FMT "+" GPLOT_DBL_FMT "*cos(u)*cos(v),"
            GPLOT_DBL_FMT "+" GPLOT_DBL_FMT "*cos(u)*sin(v),"
            GPLOT_DBL_FMT "+" GPLOT_DBL_FMT "*sin(u)",
            x, r, y, r, z, r);
    gplot_end_clause(interface, style, "pm3d");
}


static inline void
    gplot_draw_solid_triangle_3d_impl(s_gplot *interface,
                                      const double v0[3], const double v1[3], const double v2[3],
                                      s_gplot_style style)
{
    if (!gplot_ok(interface)) return;
    if (!gplot_is3d(interface)) return;
    gplot_ignore_with(&style, "gplot_draw_solid_triangle_3d");

    gplot_begin_datablock(interface);
    fprintf(interface->pipe,
            GPLOT_DBL_FMT " " GPLOT_DBL_FMT " " GPLOT_DBL_FMT "\n"
            GPLOT_DBL_FMT " " GPLOT_DBL_FMT " " GPLOT_DBL_FMT "\n"
            GPLOT_DBL_FMT " " GPLOT_DBL_FMT " " GPLOT_DBL_FMT "\n",
            v0[0], v0[1], v0[2],
            v1[0], v1[1], v1[2],
            v2[0], v2[1], v2[2]);
    gplot_end_datablock(interface, style, "polygons");
}


static inline void gplot_draw_polytope_3d_impl(s_gplot *interface,
                                               int N, const double x[N], const double y[N], const double z[N],
                                               s_gplot_style style)
{
    if (!gplot_ok(interface)) return;
    if (!gplot_is3d(interface)) return;
    if (N <= 0) return;
    gplot_ignore_with(&style, "gplot_draw_polytope_3d");

    gplot_begin_datablock(interface);
    for (int i=0; i<N; i++) {
        fprintf(interface->pipe,
                GPLOT_DBL_FMT " " GPLOT_DBL_FMT " " GPLOT_DBL_FMT "\n",
                x[i], y[i], z[i]);
    }
    gplot_end_datablock(interface, style, "polygons");
}



static inline void gplot_draw_polygons_3d_impl(s_gplot *interface, int N_polys, const int offset[],
                                               const double x[], const double y[], const double z[],
                                               const unsigned rgb[], s_gplot_style style)
{   // One datablock, polygons separated by blank lines; per-polygon color in column 4
    if (!gplot_ok(interface)) return;
    if (!gplot_is3d(interface)) return;
    if (N_polys <= 0) return;
    gplot_ignore_with(&style, "gplot_draw_polygons_3d");

    gplot_begin_datablock(interface);
    for (int i=0; i<N_polys; i++) {
        for (int k=offset[i]; k<offset[i+1]; k++) {
            fprintf(interface->pipe, GPLOT_DBL_FMT " " GPLOT_DBL_FMT " " GPLOT_DBL_FMT, x[k], y[k], z[k]);
            if (rgb) fprintf(interface->pipe, " %u", rgb[i] & 0xFFFFFFu);
            fputc('\n', interface->pipe);
        }
        fputc('\n', interface->pipe);
    }
    gplot_pipe_and_newline(interface->pipe, "EOD");

    int N = gplot_begin_clause(interface);
    fprintf(interface->plot_cmd, "$gplot_%d u 1:2:3%s", N, rgb ? ":4" : "");
    // gnuplot wants: w polygons <fill style (.extra)> fc <color> title. Polygons take a
    // fill color (fc), not a line color (lc), and fc is rejected before the fill style.
    FILE *out = interface->plot_cmd;
    fputs(" w polygons", out);
    if (style.lw > 0) fprintf(out, " lw %g", style.lw);
    if (style.extra)  fprintf(out, " %s", style.extra);
    if (rgb)                fputs(" fc rgb variable", out);
    else if (style.color)   fprintf(out, " fc rgb \"%s\"", style.color);
    if (style.title) fprintf(out, " title \"%s\"", style.title);
    else             fputs(" notitle", out);
}

#endif

/* MIT License.
 *
 * Copyright (c) 2026 Fernando Muñoz.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

